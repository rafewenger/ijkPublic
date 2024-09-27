/*!
 *  @file ijkgenMCtable.cpp
 *  @brief Generate Marching Cubes isosurface lookup table.
 *  - Version 0.5.1
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2001-2024 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "ijk.tpp"
#include "ijkbits.tpp"
#include "ijkcube.tpp"

#include "ijkgenMCpatch.h"
#include "ijkMCtable.h"
#include "ijkMCtable_xitIO.h"

#include "ijkMCtable_orient.tpp"

using namespace std;
using namespace IJKMCUBE_TABLE;

// global variables
const int DIM2(2);
const int DIM3(3);
int dimension = 0;
POLYTOPE_SHAPE mesh_polytope_type = CUBE;
bool flag_verbose = true;
int output_trigger = 5000;
const int MAX_HYPERCUBE_DIMENSION = 4;
const int MAX_PYRAMID_DIMENSION = 5;
const int MAX_SIMPLEX_DIMENSION = 20;
bool flag_interval_volume = false;
bool flag_nep = false;
bool use_cube_lex_order(true);
bool flag_check_all_pairs = false;   // Slow check version, checks all pairs.
IJKXIO::XIT_VERSION_TYPE xit_version = IJKXIO::XIT_VERSION_2_0;
char * output_filename = NULL;
const char * table_string = NULL;
const char * isosurface_table_string = "isosurface table";
const char * interval_volume_table_string = "interval volume table";

// If true, separate opposite vertices.
IJK::BOOLEAN_SET_VALUE flag_separate_opposite(false);

// If true, separate positive vertices.
IJK::BOOLEAN_SET_VALUE flag_separate_positive(false);

// If true, generate lookup table with triangulation based on edge groups.
IJK::BOOLEAN_SET_VALUE flag_edge_groups(true);
// If true, consistently orient all table entries.
bool flag_orient(true);         
// If true, orient so outward normal points to positive region.
bool flag_pos_orient(false);    

// Label output (.xit) file name.
bool flag_label_fname_with_triangulation_type(false);
bool flag_label_fname_with_separation_type(false);
bool flag_label_fname_with_orientation(false);

// routines
void generate_isosurface_table(ISOSURFACE_TABLE & isotable);
void generate_edge_groups_isosurface_table(ISOSURFACE_TABLE & isotable);
void generate_nep_isosurface_table(ISOSURFACE_TABLE & isotable);
void generate_interval_volume_table(ISOSURFACE_TABLE & isotable);
void orient_marching_cubes_table(ISOSURFACE_TABLE & isotable);
void check_dimension(const int d,
                     const POLYTOPE_SHAPE mesh_poly_type);
void check_num_vertices(const int d,
                        const POLYTOPE_SHAPE mesh_poly_type,
                        const int max_numv);

// I/O routines
void write_isotable
(const char * output_filename, const ISOSURFACE_TABLE & isotable);

// Misc routines
void parse_command_line(int argc, char **argv);
void memory_exhaustion();
void usage_error();
void help_msg();


int main(int argc, char **argv)
{
  std::set_new_handler(memory_exhaustion);
  ostringstream isotable_filename;

  try {

    parse_command_line(argc, argv);

    if (dimension == 0) {
      cout << "Enter dimension: ";
      cin >> dimension;
    };

    check_dimension(dimension, mesh_polytope_type);

    if (dimension == DIM2) {
      // Ignore flag_separate_opposite in dimension 2.
      flag_separate_opposite.Unset();
    }

    ISOSURFACE_TABLE isotable(dimension);

    if (!flag_interval_volume) {
      table_string = isosurface_table_string;
      isotable.SetSimplexDimension(dimension-1);
      if (!flag_nep) {
        isotable.SetBinaryEncoding();
      }
      else {
        isotable.SetBase3Encoding();
      };
    }
    else {
      assert(!flag_nep);

      table_string = interval_volume_table_string;
      isotable.SetBase3Encoding();
      isotable.SetSimplexDimension(dimension);
    };

    // Initialize simplex orientation to NO_ORIENT.
    isotable.SetIsoPolyOrientation(NO_ORIENT);

    check_num_vertices(dimension, mesh_polytope_type,
                       isotable.MaxNumVertices());

    switch(mesh_polytope_type) {

    case CUBE:
    default:
      if (use_cube_lex_order) {
        isotable.GenCube(dimension);
      }
      else {
        isotable.GenCubeOrderA(dimension);
      }
      break;

    case SIMPLEX:
      isotable.GenSimplex(dimension);
      break;

    case PYRAMID:
      isotable.GenPyramid(dimension);
      break;

    case SIMPLEX_PRISM:
      ISOTABLE_POLY base_polytope(dimension-1),
        prism(dimension);
      base_polytope.GenSimplex(dimension-1);
      generate_prism(base_polytope, prism);
      prism.SetShape(SIMPLEX_PRISM);
      isotable.Set(prism);
      break;
    };

    IJK::ERROR error;
    if (!isotable.Polytope().Check(error)) {
      cerr << "Problem constructing polytope." << endl;
      error.Print(cerr);
      cerr << "Exiting." << endl;
      exit(10);
    };

    unsigned long num_table_entries = 0;
    int nume = isotable.Polytope().NumEdges();
    int numv = isotable.Polytope().NumVertices();
    if (!flag_interval_volume) {
      if (!flag_nep) {
        int num_isosurface_vertices = nume;
        isotable.SetNumIsosurfaceVertices(num_isosurface_vertices);
        isotable.StorePolyEdgesAsIsoVertices(0);
        num_table_entries = calculate_num_entries(numv, 2);
      }
      else {
        int num_isosurface_vertices = nume+numv;
        isotable.SetNumIsosurfaceVertices(num_isosurface_vertices);
        isotable.StorePolyVerticesAsIsoVertices(0);
        isotable.StorePolyEdgesAsIsoVertices(numv);
        num_table_entries = calculate_num_entries(numv, 3);
      }
    }
    else {
      int num_isosurface_vertices = 2*nume+numv;
      isotable.SetNumIsosurfaceVertices(num_isosurface_vertices);
      for (int ie = 0; ie < nume; ie++) {
        isotable.SetIsoVertexFace(2*ie, ie);
        isotable.SetIsoVertexType(2*ie, ISOSURFACE_VERTEX::EDGE);
        isotable.SetIsoVertexLabel(2*ie, "0");
        isotable.SetIsoVertexFace(2*ie+1, ie);
        isotable.SetIsoVertexType(2*ie+1, ISOSURFACE_VERTEX::EDGE);
        isotable.SetIsoVertexLabel(2*ie+1, "1");
      };

      for (int iv = 0; iv < numv; iv++) {
        isotable.SetIsoVertexFace(iv+2*nume, iv);
        isotable.SetIsoVertexType(iv+2*nume, ISOSURFACE_VERTEX::VERTEX);
      }
      num_table_entries = calculate_num_entries(numv, 3);
    };

    isotable.SetNumTableEntries(num_table_entries);

    if (flag_verbose) {
      cout << isotable.Polytope().Dimension() << "D ";
      switch(mesh_polytope_type) {

      case SIMPLEX:
        cout << "simplex has ";
        break;

      case CUBE:
        cout << "hypercube has ";
        break;

      case PYRAMID:
        cout << "pyramid has ";
        break;

      case SIMPLEX_PRISM:
        cout << "prism with simplex base has ";
        break;
      };

      cout << isotable.Polytope().NumVertices()
           << " vertices and " << isotable.Polytope().NumEdges()
           << " edges." << endl;
      string table_string2 = table_string;
      if (table_string2.length() > 0) {
        table_string2[0] = toupper(table_string2[0]);
      }
      cout << table_string2 << " has "
           << isotable.NumTableEntries() << " entries." << endl;
    };

    if (flag_verbose) {
      cout << "Generating " << table_string << "." << endl;
    }

    if (flag_interval_volume) {
      generate_interval_volume_table(isotable);
    }
    else if (flag_nep) {
      generate_nep_isosurface_table(isotable);
    }
    else if (flag_edge_groups.Value() &&
             mesh_polytope_type == CUBE &&
             isotable.Polytope().Dimension() == DIM3) {
      generate_edge_groups_isosurface_table(isotable);
    }
    else {
      // Regular Marching Cubes isosurface lookup table.
      generate_isosurface_table(isotable);
    }

    if (!isotable.CheckTable(error)) {
      cerr << "Error detected in isosurface table." << endl;
      error.Print(cerr);
      cerr << "Exiting." << endl;
      exit(10);
    };

    if (flag_orient)
      { orient_marching_cubes_table(isotable); }

    if (flag_separate_positive.IsSetAndTrue()) {
      ISOSURFACE_TABLE isotableB(dimension);
      invert_mcube_isotable(isotable, isotableB);

      // Correct orientation, if necessary.
      if ((flag_pos_orient &&
           (isotableB.IsoPolyOrientation() == NEGATIVE_ORIENT)) ||
          (!flag_pos_orient &&
           (isotableB.IsoPolyOrientation() == POSITIVE_ORIENT))) {
        isotableB.FlipAllIsoPolyOrientations();
      }
      write_isotable(output_filename, isotableB);
    }
    else {
      write_isotable(output_filename, isotable);
    }
  }
  catch (IJK::ERROR & error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

  return 0;
}


void generate_isosurface_table(ISOSURFACE_TABLE & isotable)
{
  const int dimension = isotable.Dimension();
  const int numv = isotable.Polytope().NumVertices();
  const int num_cube_vert = IJK::compute_num_cube_vertices(dimension);
  int vertex_sign[numv];
  vector<ISOSURFACE_VERTEX_INDEX> edge_list;
  IJK::PROCEDURE_ERROR error("generate_isosurface_table");

  if (flag_separate_opposite.Value() &&
      isotable.Polytope().Shape() != CUBE) {
    error.AddMessage("Programming error.  Separate opposites requires cube.");
    error.AddMessage
      ("Polytope is a ", isotable.Polytope().ShapeString(), ".");
    throw error;
  }

  isotable.SetTableType(ISOSURFACE);
  isotable.SetBinaryEncoding();
  isotable.SetGridVertexLabelType(NEG_POS);
  isotable.SetSeparationType(SEPARATE_NEG);
  isotable.SetTriangulationType(CONVEX_HULL);
  if ((isotable.Polytope().Shape() == CUBE) &&
      (flag_separate_opposite.IsSet()))
    { isotable.SetSeparateOpposite(flag_separate_opposite.Value()); }

  for (int it = 0; it < isotable.NumTableEntries(); it++) {
    edge_list.clear();
    int num_simplices = 0;

    int it2 = it;
    if (flag_separate_opposite.Value() && dimension > 2) {
      int num_zeros, num_ones;
      IJK::count_bits(it, numv, num_zeros, num_ones);

      if (IJK::equals_reverse_bits(it, numv) && (num_ones == 2))
        { it2 = isotable.NumTableEntries() - it - 1; }
    }
    convert2base(it2, 2, vertex_sign, numv, error);

    for (int j = 0; j < numv; j++) {
      // convert 0 to -1
      vertex_sign[j] = 2*vertex_sign[j]-1;
      // if vertex_sign[j] < 0, vertex j is negative
      // if vertex_sign[j] > 0, vertex j is positive
    }

    gen_isopatch(isotable.Polytope(), vertex_sign, edge_list, num_simplices);

    isotable.SetNumSimplices(it, num_simplices);

    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.Dimension(); j++) {
        int ie = edge_list[is*isotable.Dimension() + j];
        isotable.SetSimplexVertex(it, is, j, ie);
      };
    };

    if (flag_verbose && it > 0 && it%output_trigger == 0) {
      cout << "  " << it << " out of " << isotable.NumTableEntries()
           << " isosurface table entries completed." << endl;
    };
  };

  if (flag_verbose && isotable.NumTableEntries() &&
      isotable.NumTableEntries() > output_trigger) {
    cout << "  All " << isotable.NumTableEntries()
         << " isosurface table entries completed." << endl;
  }

}


void generate_edge_groups_isosurface_table(ISOSURFACE_TABLE & isotable)
{
  const int DIM3(3);
  const int dimension = isotable.Dimension();
  const int numv = isotable.Polytope().NumVertices();
  const int num_cube_vert = IJK::compute_num_cube_vertices(dimension);
  int vertex_sign[numv];
  vector<ISOSURFACE_VERTEX_INDEX> edge_list;
  IJK::PROCEDURE_ERROR error("generate_isosurface_table");

  if (dimension != DIM3) {
    error.AddMessage
      ("Programming error. Edge groups only for cubes in 3D.");
    throw error;
  }

  if (flag_separate_opposite.Value() &&
      isotable.Polytope().Shape() != CUBE) {
    error.AddMessage("Programming error.  Separate opposites requires cube.");
    error.AddMessage
      ("Polytope is a ", isotable.Polytope().ShapeString(), ".");
    throw error;
  }

  isotable.SetTableType(ISOSURFACE);
  isotable.SetBinaryEncoding();
  isotable.SetGridVertexLabelType(NEG_POS);
  isotable.SetSeparationType(SEPARATE_NEG);
  isotable.SetTriangulationType(EDGE_GROUPS);
  isotable.SetSeparateOpposite(flag_separate_opposite.Value());

  ISOTABLE_CUBE_3D cube;
  cube.GenCube3D();

  for (int it = 0; it < isotable.NumTableEntries(); it++) {
    edge_list.clear();
    int num_simplices = 0;

    int it2 = it;
    if (flag_separate_opposite.Value() && dimension > 2) {
      int num_zeros, num_ones;
      IJK::count_bits(it, numv, num_zeros, num_ones);

      if (IJK::equals_reverse_bits(it, numv) && (num_ones == 2))
        { it2 = isotable.NumTableEntries() - it - 1; }
    }
    convert2base(it2, 2, vertex_sign, numv, error);

    for (int j = 0; j < numv; j++) {
      // convert 0 to -1
      vertex_sign[j] = 2*vertex_sign[j]-1;
      // if vertex_sign[j] < 0, vertex j is negative
      // if vertex_sign[j] > 0, vertex j is positive
    }

    gen_isopatch_edge_groups(cube, vertex_sign, edge_list, num_simplices);

    isotable.SetNumSimplices(it, num_simplices);
    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.Dimension(); j++) {
        int ie = edge_list[is*isotable.Dimension() + j];
        isotable.SetSimplexVertex(it, is, j, ie);
      };
    };

    if (flag_verbose && it > 0 && it%output_trigger == 0) {
      cout << "  " << it << " out of " << isotable.NumTableEntries()
           << " isosurface table entries completed." << endl;
    };
  };

  if (flag_verbose && isotable.NumTableEntries() &&
      isotable.NumTableEntries() > output_trigger) {
    cout << "  All " << isotable.NumTableEntries()
         << " isosurface table entries completed." << endl;
  }

}


void generate_nep_isosurface_table(ISOSURFACE_TABLE & isotable)
{
  IJK::PROCEDURE_ERROR error("generate_nep_isosurface_table");

  int vertex_sign[isotable.Polytope().NumVertices()];
  vector<ISOSURFACE_VERTEX_INDEX> isov_list;
  vector<ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE> isov_type;
  const int numv = isotable.Polytope().NumVertices();

  isotable.SetTableType(ISOSURFACE);
  isotable.SetBase3Encoding();
  isotable.SetGridVertexLabelType(NEG_EQUALS_POS);
  isotable.SetSeparationType(SEPARATE_NEG);
  isotable.SetTriangulationType(CONVEX_HULL);

  for (int it = 0; it < isotable.NumTableEntries(); it++) {
    isov_list.clear();
    int num_simplices = 0;

    convert2base(it, 3, vertex_sign, numv, error);
    for (int j = 0; j < numv; j++) {
      // convert 0 to -1, 1 to 0 and 2 to 1
      vertex_sign[j] = vertex_sign[j]-1;
      // if vertex_sign[j] < 0, vertex j is negative
      // if vertex_sign[j] == 0, vertex j is zero
      // if vertex_sign[j] > 0, vertex j is positive
    }

    gen_isopatch_nep(isotable.Polytope(), vertex_sign, 
                     isov_list, isov_type, num_simplices);
    assert(isov_list.size() == isov_type.size());

    isotable.SetNumSimplices(it, num_simplices);

    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.Dimension(); j++) {
        int k = is*isotable.Dimension() + j;
        int kf = isov_list[k];
        if (isov_type[k] == ISOSURFACE_VERTEX::VERTEX) {
          isotable.SetSimplexVertex(it, is, j, kf);
        }
        else {
          assert(isov_type[k] == ISOSURFACE_VERTEX::EDGE);
          isotable.SetSimplexVertex(it, is, j, kf+numv);
        }
      };
    };

    if (flag_verbose && it > 0 && it%output_trigger == 0) {
      cout << "  " << it << " out of " << isotable.NumTableEntries()
           << " isosurface table entries completed." << endl;
    };
  };

  if (flag_verbose && isotable.NumTableEntries() &&
      isotable.NumTableEntries() > output_trigger) {
    cout << "  All " << isotable.NumTableEntries()
         << " isosurface table entries completed." << endl;
  }

}


void generate_interval_volume_table(ISOSURFACE_TABLE & isotable)
{
  int index_base3[isotable.Polytope().NumVertices()];
  int vertex_prism_sign[2*isotable.Polytope().NumVertices()];
  vector<ISOSURFACE_VERTEX_INDEX> edge_list;
  const int numv = isotable.Polytope().NumVertices();
  const int nume = isotable.Polytope().NumEdges();
  IJK::PROCEDURE_ERROR error("generate_interval_volume_table");

  isotable.SetTableType(INTERVAL_VOLUME);
  isotable.SetBase3Encoding();
  isotable.SetGridVertexLabelType(NEG_STAR_POS);
  isotable.SetSeparationType(SEPARATE_NEG);
  isotable.SetTriangulationType(CONVEX_HULL);

  ISOTABLE_POLY prism(isotable.Dimension()+1);

  generate_prism(isotable.Polytope(), prism);
  const int num_prism_vertices = prism.NumVertices();

  for (int it = 0; it < isotable.NumTableEntries(); it++) {
    edge_list.clear();
    int num_simplices = 0;

    convert2base(it, 3, index_base3, numv, error);
    for (int i = 0; i < numv; i++) {
      if (index_base3[i] == 0) {
        vertex_prism_sign[i] = -1;
        vertex_prism_sign[i+numv] = -1;
      }
      else if (index_base3[i] == 2) {
        vertex_prism_sign[i] = 1;
        vertex_prism_sign[i+numv] = 1;
      }
      else {
        vertex_prism_sign[i] = 1;
        vertex_prism_sign[i+numv] = -1;
      }
    }

    gen_isopatch(prism, vertex_prism_sign, edge_list, num_simplices);

    isotable.SetNumSimplices(it, num_simplices);

    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
        int prism_ie = edge_list[is*isotable.NumVerticesPerSimplex() + j];
        int ie;

        if (prism_ie < nume) {
          ie = 2*prism_ie;
        }
        else if (prism_ie < 2*nume) {
          ie = 2*(prism_ie-nume)+1;
        }
        else if (prism_ie < 2*nume+numv) {
          ie = prism_ie;
        }
        else {
          error.AddMessage("Programming error. Illegal edge ", ie,
                           " in isosurface vertex list");
          throw error;
        }
        isotable.SetSimplexVertex(it, is, j, ie);

      };
    };

    if (flag_verbose && it > 0 && it%output_trigger == 0 &&
        isotable.NumTableEntries() > output_trigger) {
      cout << "  " << it << " out of " << isotable.NumTableEntries()
           << " interval volume table entries completed." << endl;
    };
  };

  if (flag_verbose && isotable.NumTableEntries() &&
      isotable.NumTableEntries() > output_trigger) {
    cout << "  All " << isotable.NumTableEntries()
         << " interval volume table entries completed." << endl;
  }
}


void sort_simplex_vertices_to_orient_first_table_entry
(ISOSURFACE_TABLE & isotable, const TABLE_INDEX istart,
 IJK::ERROR & error)
{
  if (isotable.NumSimplices(istart) < 1) {
    error.AddMessage
      ("Programming error. Table entry ", istart, " has no simplices.");
    error.AddMessage
      ("  Orientation must start from a table entry with at least one simplex.");
    throw error;
  }

  isotable.SortSimplexVertices(istart, 0);
}


void orient_marching_cubes_table(ISOSURFACE_TABLE & isotable)
{
  IJK::PROCEDURE_ERROR error("orient_marching_cubes_table");

  if (isotable.NumTableEntries() <= 0) {
    // Nothing to orient.
    return;
  }

  int istart = 1;

  // Initialize first simplex in table entry istart so that
  //   simplex normals point toward negative region.
  if (flag_interval_volume) {
    istart = 1;
    sort_simplex_vertices_to_orient_first_table_entry
      (isotable, istart, error);
  }
  else if (flag_nep) {
    istart = 2;
    sort_simplex_vertices_to_orient_first_table_entry
      (isotable, istart, error);
  }
  else {
    // Standard Marching Cubes isosurface lookup table
    //   with "+" and "-" vertex labels.
    istart = 1;
    sort_simplex_vertices_to_orient_first_table_entry
      (isotable, istart, error);
  }

  int num_components;
  isotable.OrientAllSimplicesInTableEntry(istart, num_components);
  if (num_components == 1) {

    if (flag_verbose && (isotable.NumTableEntries() > output_trigger))
      { cout << "Orienting simplices in isosurface table." << endl; }

    if (flag_verbose) {
      orient_mcube_table(cout, isotable, istart, flag_verbose, output_trigger);
    }
    else {
      orient_mcube_table(isotable, istart);
    }

    if (flag_pos_orient) {
      // Reverse all simplex orientations so that all normals
      //   point toward positive region.
      isotable.FlipAllIsoPolyOrientations();
      isotable.SetIsoPolyOrientation(POSITIVE_ORIENT);
    }
    else {
      isotable.SetIsoPolyOrientation(NEGATIVE_ORIENT);
    }

    if (!isotable.CheckTable(error)) {
      cerr << "Error detected in isosurface table after call to orient_MC_table()."
           << endl;
      error.Print(cerr);
      cerr << "Exiting." << endl;
      exit(10);
    }

    if (flag_verbose && isotable.NumTableEntries() > output_trigger)
      { cout << "Checking simplex orientations." << endl; }

    bool check_result;
    if (flag_verbose) {
      check_result = check_mcube_table_orientation
        (cout, isotable, flag_verbose, output_trigger,
         flag_check_all_pairs, error);
    }
    else {
      check_result =
        check_mcube_table_orientation
        (isotable, flag_check_all_pairs, error);
    }

    if (!check_result) {
      cerr << "*** Warning: Inconsistent orientations in isosurface lookup table."
           << endl;
      error.Print(cerr);
    }

  }
  else {
    cerr << "*** Warning: Unable to orient isosurface lookup table." << endl;
    cerr << "  Table entry " << istart
         << " had multiple connected components." << endl;
  }
}


std::string create_isotable_filename
(const char * output_filename, const ISOSURFACE_TABLE & isotable)
{
  if (output_filename == NULL) {
    std::ostringstream isotable_filename;


    // generate isosurface table file name
    isotable_filename.str("");
    if (!flag_interval_volume) 
      { isotable_filename << "iso."; }
    else 
      { isotable_filename << "ivol."; };

    if (flag_nep) 
      { isotable_filename << "nep."; }

    switch (mesh_polytope_type) {

    case SIMPLEX:
      isotable_filename << "simplex";
      break;

    case PYRAMID:
      isotable_filename << "pyramid";
      break;

    case SIMPLEX_PRISM:
      isotable_filename << "sprism";
      break;

    case CUBE:
    default:
      isotable_filename << "cube";
      break;
    };

    if (flag_label_fname_with_triangulation_type) {
      const std::string triangulation_type_label =
        isotable.Properties().TriangulationTypeLabel();
      if (triangulation_type_label != "") 
        { isotable_filename << "." << triangulation_type_label; }
    }

    if (flag_label_fname_with_separation_type) {
      const std::string separation_type_label =
        isotable.Properties().SeparationTypeLabel();
      if (separation_type_label != "") 
        { isotable_filename << "." << separation_type_label; }
    }

    if (flag_label_fname_with_orientation) {
      const std::string orientation_label =
        isotable.Properties().IsoPolyOrientationLabel();
      if (orientation_label != "") 
        { isotable_filename << "." << orientation_label; }
    }

    isotable_filename << "." << isotable.Dimension()
                      << "D.xit";

    return isotable_filename.str();
  }
  else {
    return std::string(output_filename);
  }
}


void write_isotable
(const char * output_filename, const ISOSURFACE_TABLE & isotable)
{
  const std::string isotable_filename =
    create_isotable_filename(output_filename, isotable);

  if (flag_verbose)
    cout << "Writing " << table_string << " to file "
         << isotable_filename << "." << endl;

  ofstream isotable_file(isotable_filename.c_str(), ios::out);
  if (!isotable_file) {
    cerr << "Unable to open file " << isotable_filename << endl;
    cerr << "Exiting." << endl;
    exit(15);
  };

  try {
    IJKXIO::write_xit(isotable_file, xit_version, isotable);
  }
  catch(...) {
    cerr << "Error writing file " << isotable_filename << "." << endl;
    throw;
  }
}


void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc) {

    if (argv[iarg][0] != '-')
      usage_error();

    string s = argv[iarg];
    if (s == "-d") {
      iarg++;
      if (iarg >= argc)
        usage_error();

      sscanf(argv[iarg], "%d", &dimension);
    }
    else if (strcmp(argv[iarg], "-V") == 0) {
      flag_verbose = true;
      iarg++;
      if (iarg >= argc)
        usage_error();

      sscanf(argv[iarg], "%d", &output_trigger);
    }
    else if (strcmp(argv[iarg], "-poly") == 0) {
      iarg++;
      if (iarg >= argc)
        usage_error();

      if (strcmp(argv[iarg], "cube") == 0)
        mesh_polytope_type = CUBE;
      else if (strcmp(argv[iarg], "simplex") == 0)
        mesh_polytope_type = SIMPLEX;
      else if (strcmp(argv[iarg], "pyramid") == 0)
        mesh_polytope_type = PYRAMID;
      else if (strcmp(argv[iarg], "sprism") == 0)
        mesh_polytope_type = SIMPLEX_PRISM;
      else
        usage_error();
    }
    else if (strcmp(argv[iarg], "-ivol") == 0) {
      flag_interval_volume = true;
    }
    else if (strcmp(argv[iarg], "-nep") == 0) {
      flag_nep = true;
    }
    else if (strcmp(argv[iarg], "-edge_groups") == 0) {
      flag_edge_groups.Set(true);
    }
    else if (strcmp(argv[iarg], "-chull") == 0) {
      flag_edge_groups.Set(false);
    }
    else if (strcmp(argv[iarg], "-o") == 0) {
      iarg++;
      if (iarg >= argc)
        usage_error();

      output_filename = argv[iarg];
    }
    else if (s == "-lexorder")
      { use_cube_lex_order = true; }
    else if (s == "-orderA")
      { use_cube_lex_order = false; }
    else if (s == "-sep_pos") 
      { flag_separate_positive.Set(true); }
    else if (s == "-sep_neg") 
      { flag_separate_positive.Set(false); }
    else if (s == "-sep_opp")
      { flag_separate_opposite.Set(true); }
    else if (s == "-no_sep_opp")
      { flag_separate_opposite.Set(false); }
    else if (s == "-skip_orient")
      { flag_orient = false; }
    else if (s == "-pos_orient") {
      flag_orient = true;
      flag_pos_orient = true;
    }
    else if (s == "-neg_orient") {
      flag_orient = true;
      flag_pos_orient = false;
    }    
    else if (s == "-check_all_pairs") {
      flag_check_all_pairs = true;
    }
    else if (s == "-label_fname_with_tri_type") {
      // Label output filename with triangulation type.
      flag_label_fname_with_triangulation_type= true;
    }
    else if (s == "-label_fname_with_sep_type") {
      flag_label_fname_with_separation_type = true;
    }
    else if (s == "-label_fname_with_orientation") {
      flag_label_fname_with_orientation = true;
    }
    else if (s == "-label_fname_with_all_flags") {
      flag_label_fname_with_triangulation_type= true;
      flag_label_fname_with_separation_type = true;
      flag_label_fname_with_orientation = true;
    }
    else if (s == "-version1") {
      xit_version = IJKXIO::XIT_VERSION_1_0;
    }
    else if (s == "-version2") {
      xit_version = IJKXIO::XIT_VERSION_2_0;
    }
    else {
      for (int j = 1; argv[iarg][j] != '\0'; j++) {
        if (argv[iarg][j] == 'c')
          mesh_polytope_type = CUBE;
        else if (argv[iarg][j] == 'h')
          help_msg();
        else if (argv[iarg][j] == 'q')
          flag_verbose = false;
        else if (argv[iarg][j] == 's')
          mesh_polytope_type = SIMPLEX;
        else if (argv[iarg][j] == 'v')
          flag_verbose = true;
        else {
          cerr << "Usage error. Illegal option: "
               << argv[iarg] << endl;
          usage_error();
        }
      };
    };
    iarg++;
  };


  if (flag_interval_volume && flag_nep) {
    cerr << "Error. Interval volume not yet implemented with -nep flag."
         << endl;
    exit(10);
  }

  if (flag_edge_groups.IsSetAndTrue()) {
    if (mesh_polytope_type != CUBE) {
      cerr << "Usage error.  Option -edge_groups can only be used with cube."
           << endl;
      exit(20);
    }

    if (flag_interval_volume) {
      cerr << "Usage error. Option -ivol can only be used with cube."
           << endl;
      exit(20);
    }

    if (flag_nep) {
      cerr << "Usage error. Option -edge_groups not available with -nep."
           << endl;
      exit(20);
    }
  }

  // Process flag_separate_opposite
  if (mesh_polytope_type == CUBE) {
    if (!flag_separate_opposite.IsSet()) {
      // For CUBE, default to separate opposite.
      flag_separate_opposite.Set(true);
    }
  }
  else if (flag_separate_opposite.IsSetAndTrue()) {
    cerr << "Usage error.  Option -sep_opp can only be used with cube."
         << endl;
    exit(20);
  }
}


void check_dimension(const int dimension,
                     const POLYTOPE_SHAPE mesh_polytope_type)
// exit if illegal dimension found
{
  if (dimension < 2) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Dimension must be at least 2." << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };

  if (mesh_polytope_type == CUBE && dimension > MAX_HYPERCUBE_DIMENSION) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Hypercube dimension is at most "
         << MAX_HYPERCUBE_DIMENSION << "."
         << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };

  if (mesh_polytope_type == SIMPLEX && dimension > MAX_SIMPLEX_DIMENSION) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Simplex dimension is at most " << MAX_SIMPLEX_DIMENSION << "."
         << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };

  if (mesh_polytope_type == PYRAMID && dimension > MAX_PYRAMID_DIMENSION) {
    cerr << "Illegal dimension: " << dimension << endl;
    cerr << "Pyramid dimension is at most " << MAX_PYRAMID_DIMENSION << "."
         << endl;
    cerr << "Exiting." << endl;
    exit(10);
  }
}

void check_num_vertices(const int d,
                        const POLYTOPE_SHAPE mesh_poly_type,
                        const int max_numv)
{
  int numv = 0;
  if (mesh_polytope_type == CUBE) {
    numv = (1L << dimension);
  }
  else {
    numv = dimension+1;
  };

  if (numv > max_numv) {
    cerr << "Polytope has too many vertices." << endl;
    cerr << d << "D ";
    if (mesh_polytope_type == CUBE) {
      cerr << "hypercube ";
    }
    else {
      cerr << "simplex ";
    };
    cerr << "has " << numv << " vertices.  Maximum number of vertices = "
         << max_numv << "."
         << endl;
    cerr << "Exiting." << endl;
    exit(10);
  };
}

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void usage_msg()
{
  cout << "Usage: ijkgenMCtable [-d dimension] [-chqsv]" << endl;
  cout << "    [-ivol | -nep | -edge_groups] [-V <n>]" << endl;
  cout << "    [-sep_neg|-sep_pos] [-sep_opp|-no_sep_opp]" << endl;
  cout << "    [-poly {cube|simplex|pyramid|sprism}]" << endl;
  cout << "    [-pos_orient|neg_orient|-skip_orient]" << endl;
  cout << "    [-label_fname_with_tri_type] [-label_fname_with_sep_type]"
       << endl;
  cout << "    [-label_fname_with_orientation] [-label_fname_with_all_flags]"
       << endl;
  cout << "    [-version1|-version2]" << endl;
  cout << "    [-o {output_filename}]" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help_msg()
{
  cout << "ijkgentable - generate isosurface or interval volume table" << endl;
  cout << "  Generate table of isosurface patches for all +/- vertex sign patterns." << endl;
  usage_msg();
  cout << "  -ivol: Generate interval volume lookup table." << endl;
  cout << "  -nep: Negative, equals, positive." << endl;
  cout << "        Differentiate scalar values equal to the isovalue" << endl;
  cout << "        from scalar values less (-) or greater (+) than the isovalue." << endl;
  cout << "  -edge_groups: Use edge groups isosurface lookup table."
       << "  (Default.)" << endl;
  cout << "        See \"Edge Groups:"
       << "  An approach to understanding the mesh quality" <<  endl
       << "        of Marching Cubes\" by Dietrich et al, IEEE TVCG, 2008."
       << endl;
  cout << "  -sep_pos: Isosurface separates positive vertices." << endl;
  cout << "  -sep_neg: Isosurface separates negative vertices. (Default.)"
       << endl;
  cout << "  -sep_opp: Separate opposite vertices. (Only for cube.)" << endl;
  cout << "        (Default for cube.)" << endl;
  cout << "  -no_sep_opp: Do not force separation of opposite vertices."
       << endl;
  cout << "  -skip_orient: Skip orientation of lookup table simplices." << endl;
  cout << "  -poly cube: Cube polytope (default.)" << endl;
  cout << "  -poly simplex: Simplex polytope." << endl;
  cout << "  -poly pyramid: Pyramid with cube (square) base." << endl;
  cout << "  -poly sprism: Prism with simplex base." << endl;
  cout << "  -pos_orient: Orient with normal pointing to positive region."
       << endl;
  cout << "  -neg_orient: Orient with normal pointing to negative region."
       << "  (Default.)" << endl;
  cout << "  -skip_orient: Skip orientation step." << endl;
  cout << "      Isosurface polytopes will not be consistently oriented."
       << endl;
  cout << "  -label_fname_with_tri_type:" << endl;
  cout << "     Add triangulation type (edgeGroups or cHull) to isosurface"
       << endl
       << "     lookup table file name." << endl;
  cout << "  -label_fname_with_sep_type:" << endl;
  cout << "     Add separation type (sepNeg or sepPos) to isosurface"
       << endl
       << "     lookup table file name." << endl;
  cout << "  -label_fname_with_orientation:" << endl;
  cout << "     Add orientation (negO or posO) to isosurface"
       << endl
       << "     lookup table file name." << endl;
  cout << "  -label_fname_with_all_flags:" << endl;
  cout << "     Add all flag labels to isosurface lookup table file name."
       << endl;
  cout << "  -version1: Write format .xit version 1." << endl;
  cout << "  -version2: Write format .xit version 2. (Default)" << endl;
  cout << "  -c: Cube polytope. (Equivalent to \"-poly cube\".)" << endl;
  cout << "  -s: Simplex polytope. (Equivalent to \"-poly simplex\".)" << endl;
  cout << "  -d <dimension> : Dimension is <dimension>." << endl;
  cout << "  -h : Print this help message (and exit.)" << endl;
  cout << "  -o {fname} : Output table to file {fname}." << endl;
  cout << "               Without this option, output file name is built " << endl;
  cout << "               from poly and dimension." << endl;
  cout << "  -q : Quiet mode." << endl;
  cout << "  -v : Verbose mode (default.)" << endl;
  cout << "  -V <n>: Verbose mode. Output after every <n> table entries."
       << endl;
  exit(0);
}
