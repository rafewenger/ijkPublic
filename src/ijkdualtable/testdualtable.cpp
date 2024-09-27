/*!
 *  @file testdualtable.cpp
 *  @brief Test ijkdualtable.
 *  - Version 0.1.2
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2012-2024 Rephael Wenger

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

#include <cstdlib>
#include <bitset>
#include <iostream>

#include "ijk.tpp"
#include "ijkbits.tpp"
#include "ijkcommand_line.tpp"
#include "ijkdualtable.tpp"
#include "ijkdualtableX.tpp"

using namespace std;
using namespace IJK;

// types
typedef typename IJKDUALTABLE::ISODUAL_TABLE_ENTRY<int,int,unsigned char> 
ISOTABLE_ENTRY_TYPE;
typedef typename IJKDUALTABLE::ISODUAL_AMBIG_TABLE_ENTRY<int,int,int,unsigned char> 
ISOTABLE_AMBIG_ENTRY_TYPE;
typedef typename IJKDUALTABLE::IVOLDUAL_TABLE_ENTRY<int,int,unsigned char,int,int>
IVOLTABLE_ENTRY_TYPE;
typedef typename 
IJKDUALTABLE::ISODUAL_CUBE_TABLE<int,int,int,ISOTABLE_ENTRY_TYPE>
ISODUAL_CUBE_TABLE;
typedef typename 
IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG<int,int,int,ISOTABLE_AMBIG_ENTRY_TYPE>
ISODUAL_CUBE_TABLE_AMBIG;
typedef typename 
IJKDUALTABLE::IVOLDUAL_CUBE_TABLE_BASE<4,int,int,int,IVOLTABLE_ENTRY_TYPE>
IVOLDUAL_CUBE_TABLE_BASE;
typedef typename 
IJKDUALTABLE::IVOLDUAL_CUBE_TABLE<4,int,int,int,IVOLTABLE_ENTRY_TYPE>
IVOLDUAL_CUBE_TABLE;
typedef typename 
IJKDUALTABLE::IVOLDUAL_CUBE_DOUBLE_TABLE<4,int,int,int,IVOLTABLE_ENTRY_TYPE>
IVOLDUAL_CUBE_DOUBLE_TABLE;
typedef typename IJKDUALTABLE::ISODUAL_CUBE_FACE_INFO<int,int,int> CUBE_TYPE;


// global constants
const int DIM2 = 2;
const int DIM3 = 3;
const int DIM4 = 4;
const int NUM_FACET_BITS16 = 16;

// global variables
bool flag_opposite_vertices(true);
bool flag_separate_neg(true);
bool flag_opposite(false);
bool flag_ivol(false);
bool flag_iso_dim4(false);
bool flag_dim2(true);
bool flag_dim3(true);
bool flag_double(false);
bool flag_only_ambig_facets(false);
bool flag_no_ambig_facets(false);
bool flag_only_non_manifold(false);
bool flag_no_below(false);
bool flag_no_above(false);
bool flag_no_I1(false);
bool flag_no_I2(false);
IJK::SET_VALUE<int> min_num_isov(0);
bool flag_at_most_one_lower_isov(false);
bool flag_at_most_one_upper_isov(false);
bool flag_at_most_one_lower_isov_opposite_sep(false);
bool flag_at_most_one_upper_isov_opposite_sep(false);
bool flag_only_LU_ambig_facets_at_most_one_isov(false);
bool flag_only_LU_ambig_facets_at_most_one_isov_all_sep(false);
bool flag_selective(false);

class ISO_VERTEX_INFO:
  public IJKDUALTABLE::VERTEX_CONNECTIVITY_INFO<NUM_FACET_BITS16>
{
public:
  int degree;
};

class IVOL_VERTEX_INFO {
public:
  typedef unsigned char CUBE_VERTEX_TYPE;
  typedef unsigned char CUBE_EDGE_TYPE;
  typedef unsigned char DIR_BITS_TYPE;

  const static CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX = 255;
  const static CUBE_EDGE_TYPE UNDEFINED_CUBE_EDGE = 255;

  int num_incident_poly;
  int num_incident_isopoly;
  CUBE_VERTEX_TYPE separation_vertex;
  CUBE_EDGE_TYPE separation_edge;
  int doubly_connected_facet;
  bool is_doubly_connected;

  /// Bit flag for connection direction.
  /// If bit i is 1, then vertex connects across facet i.
  DIR_BITS_TYPE connect_dir;

  /// Bit flag for connection direction on isosurface.
  /// If bit i is 1, then vertex connects across facet i
  ///   and edge across facet i is on the isosurface.
  DIR_BITS_TYPE iso_connect_dir;
};

typedef IJKDUALTABLE::DUAL_TABLE_VERTEX_INFO<int,ISO_VERTEX_INFO> 
DUAL_TABLE_ISO_VERTEX_INFO;

typedef IJKDUALTABLE::DUAL_TABLE_VERTEX_INFO<int,IVOL_VERTEX_INFO> 
DUAL_TABLE_IVOL_VERTEX_INFO;

// output routines
void output_isodualtable(const int dimension);
void output_isodualtable_ambig(const int dimension);
void output_isodualtable_ambig_vinfo(const int dimension);
void output_ivoldualtable(const int dimension);
void output_ivoldual_doubletable(const int dimension);
void output_opposite_zeros_or_ones(const int dimension);


// query routines
bool do_LU_ambig_facets_have_at_most_one_isov
(const IVOLDUAL_CUBE_TABLE_BASE & table, const int table_index);


// check routines
void check_isodual_table(const int dimension);
void check_isodual_ambig_table(const int dimension);
void check_contains_two_opposite(const int dimension);
void check_contains_two_opposite(const ISODUAL_CUBE_TABLE & table);
void check_isodual_vertex_info(const int dimension);
void check_isodual_vertex_info(const ISODUAL_CUBE_TABLE & table);
void check_ivoldual_table(const int dimension);
void check_ivoldual_doubletable(const int dimension);
void check_ivoldual_vertex_info(const int dimension);

void memory_exhaustion();
void usage_error(), help();
void parse_command_line(int argc, char **argv);


int main(int argc, char **argv)
{
  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);


    if (flag_ivol) {

      if (flag_double) {
        if (flag_dim2) { output_ivoldual_doubletable(DIM2); }
        if (flag_dim3) { output_ivoldual_doubletable(DIM3); }

        if (flag_dim2) { check_ivoldual_doubletable(DIM2); }
        if (flag_dim3) { check_ivoldual_doubletable(DIM3); }

        cout << "Passed ivoldual_doubletable check." << endl;
      }
      else {
        if (flag_dim2) { output_ivoldualtable(DIM2); }
        if (flag_dim3) { output_ivoldualtable(DIM3); }

        if (flag_dim2) { check_ivoldual_table(DIM2); }
        if (flag_dim3) { check_ivoldual_table(DIM3); }

        if (flag_dim2) { check_ivoldual_vertex_info(DIM2); }
        if (flag_dim3) { check_ivoldual_vertex_info(DIM3); }

        cout << "Passed ivoldual vertex info check." << endl;
        cout << "Passed ivoldual check." << endl;
      }

    }
    else {

      if (flag_opposite) {
        output_opposite_zeros_or_ones(DIM2);
        output_opposite_zeros_or_ones(DIM3);
      }
      else {
        if (flag_dim2) { output_isodualtable_ambig_vinfo(DIM2); }
        if (flag_dim3) { output_isodualtable_ambig_vinfo(DIM3); }
        if (flag_iso_dim4) { output_isodualtable_ambig_vinfo(DIM4); }
      }

      if (flag_dim2) { check_isodual_table(DIM2); }
      if (flag_dim3) { check_isodual_table(DIM3); }
      if (flag_dim2) { check_isodual_ambig_table(DIM2); }
      if (flag_dim3) { check_isodual_ambig_table(DIM3); }
      if (flag_dim2) { check_isodual_vertex_info(DIM2); }
      if (flag_dim3) { check_isodual_vertex_info(DIM3); }
      if (flag_iso_dim4) { check_isodual_vertex_info(DIM4); }
      if (flag_dim2) { check_contains_two_opposite(DIM2); }
      if (flag_dim3) { check_contains_two_opposite(DIM3); }
      if (flag_iso_dim4) { check_contains_two_opposite(DIM4); }

      cout << "Passed isodual table check." << endl;
      cout << "Passed isodual ambig table check." << endl;
      cout << "Passed isodual vertex info check." << endl;
      cout << "Passed check of contains_two_opposite_ones and contains_two_opposite_zeros." << endl;
    }
  }
  catch (ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(30);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}

// **************************************************
// I/O routines
// **************************************************

template <typename TABLE_TYPE>
void output_dualtable_header(const TABLE_TYPE & table)
{
  cout << "Dimension: " << table.Dimension() << endl;
  cout << "Num poly vertices: " << table.NumPolyVertices() << endl;
  cout << "Num poly edges: " << table.NumPolyEdges() << endl;
  cout << "Num poly facets: " << table.NumPolyFacets() << endl;
  cout << "Number of table entries: "
       << table.NumTableEntries() << endl;
}

template <typename TABLE_TYPE>
void output_active_facets
(const int it, const TABLE_TYPE & table, const CUBE_TYPE & cube_info)
{
  IJK::PROCEDURE_ERROR error("output_active_facets");

  int num_neg, num_pos;
  int num_active = cube_info.ComputeNumActiveCubeFacets(it);
  int k = 0;
  cout << "  Active facets (#-,#+):";
  for (int ifacet = 0; ifacet < cube_info.NumFacets(); ifacet++) {
    if (cube_info.IsCubeFacetActive(it, ifacet)) {
      cube_info.ComputeNumCubeFacetBits
        (it, ifacet, num_neg, num_pos);
      cout << "  " << ifacet;
      cout << " (" << num_neg << "," << num_pos << ")";
      k++;
    }
  }
  cout << endl;

  // Check number of active facets.
  if (k != num_active) {
    error.AddMessage
      ("Programming error.  Incorrect number of active facets.");
    error.AddMessage
      ("  Reported ", k, " but should be ", num_active, ".");
    throw error;
  }
 
}

template <typename TABLE_TYPE>
void output_entry_info
(const int it, const TABLE_TYPE & table, const CUBE_TYPE & cube_info)
{
  cout << "Entry " << it << ".";
  cout << "  Num isov: " << table.NumIsoVertices(it) << ".";
  cout << endl;

  if (table.NumIsoVertices(it) > 0) { 
    cout << "  Bipolar edges (dual isov):";
    for (int ie = 0; ie < table.NumPolyEdges(); ie++) {
      if (table.IsBipolar(it, ie)) {
        cout << " " << ie 
             << " (" << int(table.IncidentIsoVertex(it, ie)) << ")";
      }
    }
    cout << endl;
    output_active_facets(it, table, cube_info);
  }
}


template <typename AMBIG_TABLE_TYPE>
void output_ambig_info
(const int it, const AMBIG_TABLE_TYPE & ambig_table, 
 const CUBE_TYPE & cube_info)
{
  const int num_cube_facets = cube_info.NumFacets();

  if (ambig_table.IsAmbiguous(it))
    { cout << "  Is ambiguous."; }
  else
    { cout << "  Not ambiguous."; }
  cout << "  Num ambiguous facets = " 
       << int(ambig_table.NumAmbiguousFacets(it)) << ".";
  cout << endl;

  if (ambig_table.NumAmbiguousFacets(it) > 0) {
    cout << "  Ambiguous facets: ";
    for (int kf = 0; kf < num_cube_facets; kf++) {
      if (ambig_table.IsFacetAmbiguous(it, kf)) 
        { cout << " " << kf; }
    }
    cout << endl;
  }
}


// *** DEPRECATED/OBSOLETE ***
template <typename AMBIG_TABLE_TYPE>
void output_num_incident_edges_dual_to_facet
(const int it, const AMBIG_TABLE_TYPE & ambig_table, 
 const CUBE_TYPE & cube_info)
{
  for (int isov = 0; isov < ambig_table.NumIsoVertices(it); isov++) {
    for (int kf = 0; kf < cube_info.NumFacets(); kf++) {
      const int num_incident_edges_dual_to_facet =
	ambig_table.CountNumIncidentIsoEdgesDualToFacet
	(it, isov, kf, cube_info);

      if (num_incident_edges_dual_to_facet > 1) {
	cout << "  Num edges incident to iso vertex " << isov
	     << " and dual to facet " << kf << ": "
	     << num_incident_edges_dual_to_facet << endl;
      }
    }
  }
}


template <typename TABLE_TYPE>
void output_table_vertex_info
(const TABLE_TYPE & table,
 const DUAL_TABLE_ISO_VERTEX_INFO & vinfo, 
 const int it, const char * label)
{
  const int num_facets = table.Dimension()*2;

  cout << "  " << label << ":";
  for (int j = 0; j < vinfo.NumVertices(it); j++) {
    cout << "  " << vinfo.VertexInfo(it, j).degree;
  }
  cout << endl;

  cout << "  Intersects facets: ";
  for (int j = 0; j < vinfo.NumVertices(it); j++) {
    vinfo.VertexInfo(it, j).PrintIntersectsFacet
      (cout, num_facets, "  ", "");
  }
  cout << endl;
}


void output_table_vertex_info
(const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo, 
 const int it, const char * label)
{
  cout << "  " << label << ":";
  for (int j = 0; j < vinfo.NumVertices(it); j++) {
    cout << "  " << vinfo.VertexInfo(it, j).num_incident_poly;
    cout << " (" << vinfo.VertexInfo(it, j).num_incident_isopoly << ")";
  }
  cout << endl;
}


template <typename TABLE_TYPE>
void output_ivoldual_vertex_connection_direction
(const TABLE_TYPE & table, const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo, 
 const int it)
{
  const int NUM_FACETS(table.Dimension()*2);
  std::string bit_string;

  cout << "  Signed vertex connection directions:";
  for (int i = 0; i < vinfo.NumVertices(it); i++) {
    IJK::convert2bit_string
      (vinfo.VertexInfo(it,i).connect_dir, NUM_FACETS, bit_string);
    cout << "  " << bit_string;
  }
  cout << endl;
  cout << "  Signed vertex iso connection directions:";
  for (int i = 0; i < vinfo.NumVertices(it); i++) {
    IJK::convert2bit_string
      (vinfo.VertexInfo(it,i).iso_connect_dir, NUM_FACETS, bit_string);
    cout << "  " << bit_string;
  }
  cout << endl;
}

void output_separation_vertices
(const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo, 
 const int it)
{
  const IVOL_VERTEX_INFO::CUBE_VERTEX_TYPE 
    UNDEFINED_CUBE_VERTEX(IVOL_VERTEX_INFO::UNDEFINED_CUBE_VERTEX);

  bool flag_separation_vertex = false;

  for (int j = 0; j < vinfo.NumVertices(it); j++) {
    if (vinfo.VertexInfo(it, j).separation_vertex != UNDEFINED_CUBE_VERTEX)
      { flag_separation_vertex = true;  }
  }

  if (flag_separation_vertex) {
    cout << "  Separation vertices (ivolv, cubev):";
    for (int j = 0; j < vinfo.NumVertices(it); j++) {
      if (vinfo.VertexInfo(it, j).separation_vertex != UNDEFINED_CUBE_VERTEX) {
        cout << "  (" << j << "," 
             << int(vinfo.VertexInfo(it,j).separation_vertex) << ")";
      }
    }
    cout << endl;
  }
}


template <typename TABLE_TYPE>
void output_separation_edges
(const TABLE_TYPE & table, const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo, 
 const int it)
{
  const IVOL_VERTEX_INFO::CUBE_VERTEX_TYPE 
    UNDEFINED_CUBE_VERTEX(IVOL_VERTEX_INFO::UNDEFINED_CUBE_VERTEX);
  const IVOL_VERTEX_INFO::CUBE_EDGE_TYPE 
    UNDEFINED_CUBE_EDGE(IVOL_VERTEX_INFO::UNDEFINED_CUBE_EDGE);

  bool flag_separation_edge = false;

  for (int j = 0; j < vinfo.NumVertices(it); j++) {
    if (vinfo.VertexInfo(it, j).separation_edge != UNDEFINED_CUBE_EDGE)
      { flag_separation_edge = true;  }
  }

  if (flag_separation_edge) {
    cout << "  Separation edge (ivolv, cube edge):";
    for (int j = 0; j < vinfo.NumVertices(it); j++) {
      const int iedge = vinfo.VertexInfo(it, j).separation_edge;
      if (iedge != UNDEFINED_CUBE_EDGE) {
        const int iend0 = table.Cube().EdgeEndpoint(iedge,0);
        const int iend1 = table.Cube().EdgeEndpoint(iedge,1);
        cout << "  (" << j << ", [" << iend0 << "," << iend1 << "])";
      }
    }
    cout << endl;
  }
}


void output_doubly_connected_vertices
(const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo, 
 const int it)
{
  bool flag_doubly_connected_vertex = false;

  for (int j = 0; j < vinfo.NumVertices(it); j++) {
    if (vinfo.VertexInfo(it, j).is_doubly_connected)
      { flag_doubly_connected_vertex = true;  }
  }

  if (flag_doubly_connected_vertex) {
    cout << "  Doubly connected vertices (dc facet):";
    for (int j = 0; j < vinfo.NumVertices(it); j++) {
      if (vinfo.VertexInfo(it, j).is_doubly_connected) {
        const int dc_facet = vinfo.VertexInfo(it,j).doubly_connected_facet;
        cout << "  " << j << " (" << dc_facet << ")";
      }
    }
    cout << endl;
  }
}


template <typename TABLE_TYPE>
void output_isotable_vertex_info
(const TABLE_TYPE & table,
 const DUAL_TABLE_ISO_VERTEX_INFO & vinfo, const int it)
{
  output_table_vertex_info(table, vinfo, it, "Isov degrees");
}


template <typename TABLE_TYPE>
void output_ivoltable_vertex_info
(const TABLE_TYPE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo, const int it)
{
  output_table_vertex_info(vinfo, it, "Ivolv num incident poly (isopoly)");
  output_ivoldual_vertex_connection_direction(table, vinfo, it);
  output_separation_vertices(vinfo, it);
  if (table.Dimension() == DIM3) {
    output_separation_edges(table, vinfo, it);
  }
  if (table.Dimension() == DIM3) {
    output_doubly_connected_vertices(vinfo, it);
  }
}


void output_isodualtable(const int dimension)
{
  ISODUAL_CUBE_TABLE
    table(dimension, flag_separate_neg, flag_opposite_vertices);
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("output_isodualtable");

  output_dualtable_header(table);
  for (int it = 0; it < table.NumTableEntries(); it++) {
    output_entry_info(it, table, cube_info);
  }
  cout << endl;
}


void output_isodualtable_ambig(const int dimension)
{
  ISODUAL_CUBE_TABLE_AMBIG
    ambig_table(dimension, flag_separate_neg, flag_opposite_vertices);
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("output_isodualtable");

  output_dualtable_header(ambig_table);
  for (int it = 0; it < ambig_table.NumTableEntries(); it++) {
    output_entry_info(it, ambig_table, cube_info);

    if (ambig_table.NumIsoVertices(it) > 0) {
      output_ambig_info(it, ambig_table, cube_info);
    }
  }
  cout << endl;
}

void output_isodualtable_ambig_vinfo(const int dimension)
{
  ISODUAL_CUBE_TABLE_AMBIG
    ambig_table(dimension, flag_separate_neg, flag_opposite_vertices);
  DUAL_TABLE_ISO_VERTEX_INFO vinfo(ambig_table);
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("output_isodualtable");

  compute_dual_isotable_vertex_degrees(ambig_table, vinfo);
  compute_dual_cube_isotable_vertex_connectivity
    (ambig_table, vinfo);

  output_dualtable_header(ambig_table);
  for (int it = 0; it < ambig_table.NumTableEntries(); it++) {
    output_entry_info(it, ambig_table, cube_info);

    if (ambig_table.NumIsoVertices(it) > 0) {
      output_isotable_vertex_info(ambig_table, vinfo, it);
      output_ambig_info(it, ambig_table, cube_info);

      /* DEPRECATED/OBSOLETE
      output_num_incident_edges_dual_to_facet
	(it, ambig_table, cube_info);
      */
    }
  }
  cout << endl;
}

void output_opposite_zeros_or_ones(const int dimension)
{
  ISODUAL_CUBE_TABLE table
    (dimension, flag_separate_neg, flag_opposite_vertices);

  cout << "Dimension: " << dimension << endl;
  cout << "Number of table entries: "
       << table.NumTableEntries() << endl;
  for (int it = 0; it < table.NumTableEntries(); it++) {
    cout << "Entry " << it << ".";
    if (contains_two_opposite_ones(it, table.NumPolyVertices())) {
      cout << "  Opposite ones.";
    }

    if (contains_two_opposite_zeros(it, table.NumPolyVertices())) {
      cout << "  Opposite zeros.";
    }
    cout << endl;
  }
  cout << endl;

}

void output_ivoldualtable_info
(const IVOLDUAL_CUBE_TABLE_BASE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo,
 const int table_index,
 const CUBE_TYPE & cube_info)
{
  cout << "  Num vertices in lower lifted: " 
       << table.NumVerticesInLowerLifted(table_index);
  cout << "  Num vertices in upper lifted: " 
       << table.NumVerticesInUpperLifted(table_index);
  cout << endl;

  cout << "  Lower iso index: " << table.LowerIsosurfaceTableIndex(table_index);
  cout << "  Upper iso index: " << table.UpperIsosurfaceTableIndex(table_index);
  cout << endl;

  cout << "  Interval volume vert: ";
  for (int j = 0; j < table.NumIVolVertices(table_index); j++) {
    cout << "  " << j;
    cout << " (";
    if (table.IsInLowerLiftedCube(table_index,j)) { cout << "InLowerLC"; }
    else { cout << "InUpperLC"; };
    if (table.OnLowerIsosurface(table_index, j)) { cout << ",OnLowerS"; }
    else if (table.OnUpperIsosurface(table_index, j)) { cout << ",OnUpperS"; }
    cout << ")";
  }
  cout << endl;

  cout << "  Active edges: ";
  for (int ie = 0; ie < table.NumPolyEdges(); ie++) {
    if (table.EdgeHasDualIVolPoly(table_index, ie)) {
      cout << " " << ie 
           << " (" << int(table.LowerIncident(table_index, ie)) << ","
           << int(table.UpperIncident(table_index, ie)) << ")";
    }
  }
  cout << endl;

  cout << "  Active vertices: ";
  for (int iv = 0; iv < table.NumPolyVertices(); iv++) {
    if (table.IsInIntervalVolume(table_index, iv)) {
      cout << " " << iv
           << " (" << int(table.IncidentIVolVertex(table_index, iv)) << ")";
    }
  }
  cout << endl;

  if (table.IsNonManifold(table_index)) 
    { cout << "  Non-manifold interval volume." << endl; }

  if (table.IsAmbiguous(table_index)) {
    cout << "  Ambiguous facets:";
    for (int jf = 0; jf < cube_info.NumFacets(); jf++) {
      if (table.IsFacetAmbiguous(table_index, jf)) {
        cout << " " << jf;
      }
    }
    cout << endl;

    cout << "  Ambiguous facets in lower lifted cube:";
    for (int jf = 0; jf < cube_info.NumFacets(); jf++) {
      if (table.IsFacetInLowerLiftedAmbiguous(table_index, jf)) {
        cout << " " << jf;
      }
    }
    cout << endl;

    cout << "  Ambiguous facets in upper lifted cube:";
    for (int jf = 0; jf < cube_info.NumFacets(); jf++) {
      if (table.IsFacetInUpperLiftedAmbiguous(table_index, jf)) {
        cout << " " << jf;
      }
    }
    cout << endl;
  }

  cout << "  Relative vertex locations: ";
  for (int iv = 0; iv < table.NumPolyVertices(); iv++) {
    cout << " " << iv << ":" ;
    if (table.IsAboveIntervalVolume(table_index,iv)) 
      { cout << "A"; }
    else if (table.IsBelowIntervalVolume(table_index,iv)) 
      { cout << "B"; }
    else 
      { cout << "I" << table.VertexType(table_index,iv); }
  }
  cout << endl;

  output_ivoltable_vertex_info(table, vinfo, table_index);
}


void output_ivoldualtable_entry
(const IVOLDUAL_CUBE_TABLE_BASE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo,
 const int table_index,
 const CUBE_TYPE & cube_info)
{
  cout << "Entry " << table_index << ".";
  cout << "  Num vertices: " << table.NumIVolVertices(table_index);
  cout << endl;

  if (table.NumIVolVertices(table_index) > 0) 
    { output_ivoldualtable_info(table, vinfo, table_index, cube_info); }
}


void output_ivoldualtable_only_ambig_facets
(const IVOLDUAL_CUBE_TABLE_BASE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo)
{
  const int dimension = table.Dimension();
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("output_ivoldualtable");

  output_dualtable_header(table);
  for (int table_index = 0; table_index < table.NumTableEntries(); 
       table_index++) {
    if (table.NumAmbiguousFacets(table_index) > 0) 
      { output_ivoldualtable_entry(table, vinfo, table_index, cube_info); }
  }
  cout << endl;
}


void output_ivoldualtable_selective
(const IVOLDUAL_CUBE_TABLE_BASE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo)
{
  const int dimension = table.Dimension();
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("output_dualtable");

  output_dualtable_header(table);
  for (int table_index = 0; table_index < table.NumTableEntries(); 
       table_index++) {
    if (flag_only_ambig_facets && table.NumAmbiguousFacets(table_index) == 0)
      { continue; }

    if (flag_no_ambig_facets && table.NumAmbiguousFacets(table_index) > 0)
      { continue; }

    if (flag_only_non_manifold && !table.IsNonManifold(table_index))
      { continue; }

    if (flag_no_below && table.ContainsNegVertex(table_index))
      { continue; }

    if (flag_no_above && table.ContainsPosVertex(table_index))
      { continue; }

    if (flag_no_I1 && table.ContainsI1(table_index))
      { continue; }

    if (flag_no_I2 && table.ContainsI2(table_index))
      { continue; }

    if (min_num_isov.IsSet() &&
        table.NumIVolVertices(table_index) < min_num_isov.Value())
      { continue; }

    if (flag_at_most_one_lower_isov &&
        table.NumVerticesInLowerLifted(table_index) > 1)
      { continue; }

    if (flag_at_most_one_lower_isov &&
        table.NumVerticesInUpperLifted(table_index) > 1)
      { continue; }

    if (flag_only_LU_ambig_facets_at_most_one_isov) {
      if (table.NumAmbiguousFacets(table_index) == 0) { continue; }
      if (!do_LU_ambig_facets_have_at_most_one_isov(table, table_index))
        { continue; }
    }
    
    output_ivoldualtable_entry(table, vinfo, table_index, cube_info); 
  }
  cout << endl;
}


void output_ivoldualtable
(const IVOLDUAL_CUBE_TABLE_BASE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo)
{
  const int dimension = table.Dimension();
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("output_dualtable");

  if (flag_selective) {
    output_ivoldualtable_selective(table, vinfo);
  }
  else {
    output_dualtable_header(table);
    for (int table_index = 0; table_index < table.NumTableEntries(); 
         table_index++) 
      { output_ivoldualtable_entry(table, vinfo, table_index, cube_info); }
    cout << endl;
  }
}


void output_ivoldual_doubletable_entry
(const IVOLDUAL_CUBE_DOUBLE_TABLE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo,
 const int table_index,
 const CUBE_TYPE & cube_info)
{
  cout << "Entry " << table_index << ".";
  cout << "  Num vertices: " << table.NumIVolVertices(table_index);
  cout << "  Opposite entry: " << table.OppositeTableIndex(table_index);
  cout << endl;

  output_ivoldualtable_info(table, vinfo, table_index, cube_info);
}

void output_ivoldual_doubletable_selective
(const IVOLDUAL_CUBE_DOUBLE_TABLE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo)
{
  const int dimension = table.Dimension();
  CUBE_TYPE cube_info(dimension);

  output_dualtable_header(table);
  for (int table_index = 0; table_index < table.NumTableEntries(); 
       table_index++) {
    if (flag_only_ambig_facets && table.NumAmbiguousFacets(table_index) == 0)
      { continue; }

    if (flag_no_ambig_facets && table.NumAmbiguousFacets(table_index) > 0)
      { continue; }

    if (flag_only_non_manifold && !table.IsNonManifold(table_index))
      { continue; }

    if (flag_no_below && table.ContainsNegVertex(table_index))
      { continue; }

    if (flag_no_above && table.ContainsPosVertex(table_index))
      { continue; }

    if (flag_no_I1 && table.ContainsI1(table_index))
      { continue; }

    if (flag_no_I2 && table.ContainsI2(table_index))
      { continue; }

    if (flag_at_most_one_lower_isov &&
        table.NumVerticesInLowerLifted(table_index) > 1)
      { continue; }

    if (flag_at_most_one_lower_isov &&
        table.NumVerticesInUpperLifted(table_index) > 1)
      { continue; }

    if (min_num_isov.IsSet() &&
        table.NumIVolVertices(table_index) < min_num_isov.Value())
      { continue; }

    const int opposite_table_index = 
      table.OppositeTableIndex(table_index);
    if (flag_at_most_one_lower_isov_opposite_sep &&
        table.NumVerticesInLowerLifted(opposite_table_index) > 1)
      { continue; }

    if (flag_at_most_one_upper_isov_opposite_sep &&
        table.NumVerticesInUpperLifted(opposite_table_index) > 1)
      { continue; }

    if (flag_only_LU_ambig_facets_at_most_one_isov) {
      if (table.NumAmbiguousFacets(table_index) == 0) { continue; }
      if (!do_LU_ambig_facets_have_at_most_one_isov(table, table_index))
        { continue; }
    }

    if (flag_only_LU_ambig_facets_at_most_one_isov_all_sep) {
      if (table.NumAmbiguousFacets(table_index) == 0) { continue; }
      if (!do_LU_ambig_facets_have_at_most_one_isov(table, table_index) ||
          !do_LU_ambig_facets_have_at_most_one_isov
          (table, opposite_table_index))
        { continue; }
    }

    output_ivoldual_doubletable_entry(table, vinfo, table_index, cube_info); 
  }
  cout << endl;
}


void output_ivoldual_doubletable
(const IVOLDUAL_CUBE_DOUBLE_TABLE & table,
 const DUAL_TABLE_IVOL_VERTEX_INFO & vinfo)
{
  const int dimension = table.Dimension();
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("output_dualtable");

  if (flag_selective) {
    output_ivoldual_doubletable_selective(table, vinfo);
  }
  else {
    output_dualtable_header(table);
    for (int table_index = 0; table_index < table.NumTableEntries(); 
         table_index++) {
      output_ivoldual_doubletable_entry
        (table, vinfo, table_index, cube_info); 
    }
    cout << endl;
  }
}

void output_ivoldualtable(const int dimension)
{
  IVOLDUAL_CUBE_TABLE table(dimension, flag_separate_neg);
  DUAL_TABLE_IVOL_VERTEX_INFO vinfo(table);

  compute_ivoldual_table_num_incident_poly(table, vinfo);
  compute_ivoldual_table_num_incident_isosurface_poly(table, vinfo);
  determine_ivol_vertex_connection_directions(table, vinfo);
  determine_ivol_vertex_iso_connection_directions(table, vinfo);
  determine_separation_vertices(table, vinfo);
  determine_separation_edges(table, vinfo);

  if (dimension == DIM3)
    { determine_doubly_connected_ivol3D_vertices(table, vinfo); }

  output_ivoldualtable(table, vinfo);
}


void output_ivoldual_doubletable(const int dimension)
{
  IVOLDUAL_CUBE_DOUBLE_TABLE table(dimension, flag_separate_neg);
  DUAL_TABLE_IVOL_VERTEX_INFO vinfo(table);

  compute_ivoldual_table_num_incident_poly(table, vinfo);
  compute_ivoldual_table_num_incident_isosurface_poly(table, vinfo);
  determine_ivol_vertex_connection_directions(table, vinfo);
  determine_ivol_vertex_iso_connection_directions(table, vinfo);
  determine_separation_vertices(table, vinfo);
  determine_separation_edges(table, vinfo);

  if (dimension == DIM3)
    { determine_doubly_connected_ivol3D_vertices(table, vinfo); }

  output_ivoldual_doubletable(table, vinfo);
}

// **************************************************
// Query routines
// **************************************************

bool do_LU_ambig_facets_have_at_most_one_isov
(const IVOLDUAL_CUBE_TABLE_BASE & table,
 const int table_index)
{
  if (table.AmbiguousFacetBitsInLowerLifted(table_index) != 0 &&
      table.NumVerticesInLowerLifted(table_index) <= 1)
    return(true);

  if (table.AmbiguousFacetBitsInUpperLifted(table_index) != 0 &&
      table.NumVerticesInUpperLifted(table_index) <= 1)
    return(true);

  return(false);
}


// **************************************************
// Check routines
// **************************************************

template <typename TABLE_TYPEA, typename TABLE_TYPEB>
void check_isodual_table
(const TABLE_TYPEA & isodual_tableA, const TABLE_TYPEB & isodual_tableB,
 IJK::ERROR & error)
{
  if (isodual_tableA.Dimension() != isodual_tableB.Dimension()) {
    error.AddMessage("Error.  Tables A and B have different dimensions.");
    error.AddMessage
      ("  Table A has dimension: ", isodual_tableA.Dimension(), "");
    error.AddMessage
      ("  Table B has dimension: ", isodual_tableB.Dimension(), "");
    throw error;
  }

  if (isodual_tableA.NumTableEntries() != isodual_tableB.NumTableEntries()) {
    error.AddMessage
      ("Error.  Tables A and B have different numbers of table endtries.");
    throw error;
  }

  if (isodual_tableA.NumPolyEdges() != isodual_tableB.NumPolyEdges()) {
    error.AddMessage
      ("Error.  Tables A and B have different numbers of polytope edges.");
    throw error;
  }

  for (int it = 0; it < isodual_tableA.NumTableEntries(); it++) {
    for (int ie = 0; ie < isodual_tableA.NumPolyEdges(); ie++) {
      if (isodual_tableA.IsBipolar(it, ie) != 
          isodual_tableB.IsBipolar(it, ie)) {
        error.AddMessage
          ("Error.  Table entry ", it, " edge ", ie, 
           " is bipolar in one table but not in the other.");
        throw error;
      }

      if (isodual_tableA.IsBipolar(it, ie)) {
        if (isodual_tableA.IncidentIsoVertex(it, ie) !=
            isodual_tableB.IncidentIsoVertex(it, ie)) {
          error.AddMessage
            ("Error.  Table entry ", it, " edge ", ie, 
             " have different incident iso vertices.");
          throw error;
        }
      }
    }
  }
}


template <typename TABLE_TYPE>
void check_dimension
(const int dimension, const TABLE_TYPE & table,  IJK::ERROR & error)
{
  if (table.Dimension() != dimension) {
    error.AddMessage("Error.  Table has incorrect dimension.");
    error.AddMessage
      ("  Table has dimension: ", table.Dimension(), "");
    error.AddMessage("  Dimension should be ", dimension, "");
    throw error;
  }
}


void check_isodual_table(const int dimension)
{
  ISODUAL_CUBE_TABLE isodual_tableA(dimension);
  ISODUAL_CUBE_TABLE isodual_tableB;
  IJK::PROCEDURE_ERROR error("check_isodual_table");

  if (!isodual_tableA.Check(error)) { throw error; }
  
  isodual_tableB.Create(dimension);

  check_dimension(dimension, isodual_tableA, error);
  check_isodual_table(isodual_tableA, isodual_tableB, error);
}


void check_isodual_ambig_table(const int dimension)
{
  ISODUAL_CUBE_TABLE isodual_table(dimension);
  ISODUAL_CUBE_TABLE_AMBIG ambig_table(dimension);
  IJK::PROCEDURE_ERROR error("check_isodual_ambig_table");

  if (!isodual_table.Check(error)) { throw error; }
  
  check_dimension(dimension, ambig_table, error);
  check_isodual_table(isodual_table, ambig_table, error);
}


void check_isodual_vertex_info(const ISODUAL_CUBE_TABLE & table)
{
  DUAL_TABLE_ISO_VERTEX_INFO vinfo(table);
  IJK::PROCEDURE_ERROR error("check_isodual_vertex_info");

  compute_dual_isotable_vertex_degrees(table, vinfo);

  for (int it = 0; it < table.NumTableEntries(); it++) {

    int num_bipolar = 0;
    for (int ie = 0; ie < table.NumPolyEdges(); ie++) {
      if (table.IsBipolar(it, ie)) { num_bipolar++; }
    }

    int sum_of_degrees = 0;
    for (int j = 0; j < vinfo.NumVertices(it); j++) 
      { sum_of_degrees += vinfo.VertexInfo(it,j).degree; }

    if (sum_of_degrees != num_bipolar) {
      error.AddMessage
        ("Error.  Table entry ", it, ".  Incorrect vertex degrees.");
      error.AddMessage("  Sum of degrees: ", sum_of_degrees, "");
      error.AddMessage("  Num of bipolar edges: ", num_bipolar, "");
      error.AddMessage
        ("Sum of degrees does not equal number of bipolar edges.");
      throw error;
    }
  }
}


void check_isodual_vertex_info(const int dimension)
{
  ISODUAL_CUBE_TABLE table(dimension);

  check_isodual_vertex_info(table);
}


template <typename TABLE_TYPE>
void check_ivoldual_vertex_info(const TABLE_TYPE & table)
{
  DUAL_TABLE_IVOL_VERTEX_INFO vinfo(table);
  IJK::PROCEDURE_ERROR error("check_ivoldual_vertex_info");

  compute_ivoldual_table_num_incident_poly(table, vinfo);

  for (int it = 0; it < table.NumTableEntries(); it++) {

    int num_active_edges = 0;
    for (int ie = 0; ie < table.NumPolyEdges(); ie++) {
      if (table.EdgeHasDualIVolPoly(it, ie)) { num_active_edges++; }
    }

    int num_active_vertices = 0;
    for (int iv = 0; iv < table.NumPolyVertices(); iv++) {
      if (table.IsInIntervalVolume(it, iv)) { num_active_vertices++; }
    }

    int sum_of_num_incident = 0;
    for (int j = 0; j < vinfo.NumVertices(it); j++) 
      { sum_of_num_incident += 
          vinfo.VertexInfo(it,j).num_incident_poly; }

    if (sum_of_num_incident != 
        2*num_active_edges+num_active_vertices) {
      error.AddMessage
        ("Error.  Table entry ", it, 
         ".  Incorrect number of incident polytopes.");
      error.AddMessage("  Sum of num_incident: ", 
                       sum_of_num_incident, "");
      error.AddMessage
        ("  Sum should equal: ", 
         2*num_active_edges+num_active_vertices, "");
      throw error;
    }
  }
}

void check_ivoldual_vertex_info(const int dimension)
{
  IVOLDUAL_CUBE_TABLE table(dimension);
  IVOLDUAL_CUBE_DOUBLE_TABLE tableB(dimension);

  check_ivoldual_vertex_info(table);
  check_ivoldual_vertex_info(tableB);
}

void check_ivoldual_table(const int dimension)
{
  IVOLDUAL_CUBE_TABLE table(dimension, flag_separate_neg);
  CUBE_TYPE cube_info(dimension);
  IJK::PROCEDURE_ERROR error("check_ivoldual_table");

  if (!table.Check(error)) { throw error; }
  
  for (int it = 0; it < table.NumTableEntries(); it++) {
    for (int j = 0; j < table.NumIVolVertices(it); j++) {

      if (table.OnLowerIsosurface(it,j) &&
          table.OnUpperIsosurface(it,j)) {
        error.AddMessage
          ("Error.  Table entry ", it, " isosurface vertex ", j, "");
        error.AddMessage("  is on both lower and upper isosurface.");
        throw error;
      }
    }

    for (int k = 0; k < table.NumPolyEdges(); k++) {
      if (table.EdgeHasDualIVolPoly(it,k)) {
        if (table.LowerIncident(it,k) >= table.NumIVolVertices(it)) {
          error.AddMessage
            ("Error.  Table entry ", it, " edge ", k, "");
          error.AddMessage
            ("  has lower incident vertex ", table.LowerIncident(it,k), "");
          error.AddMessage
            ("  but table entry ", it, " has only ",
             table.NumIVolVertices(it), " interval volume vertices.");
          throw error;
        }

        if (table.UpperIncident(it,k) >= table.NumIVolVertices(it)) {
          error.AddMessage
            ("Error.  Table entry ", it, " edge ", k, "");
          error.AddMessage
            ("  has upper incident vertex ", table.UpperIncident(it,k), "");
          error.AddMessage
            ("  but table entry ", it, " has only ",
             table.NumIVolVertices(it), " interval volume vertices.");
          throw error;
        }

        if (table.LowerIncident(it,k) == table.UpperIncident(it,k)) {
          error.AddMessage
            ("Error.  Lower and upper incident vertices of edge ", k, "");
          error.AddMessage
            ("  in table ", it, " are same vertex ",
             table.LowerIncident(it,k), ".");
          throw error;
        }
      }
    }

    for (int ie = 0; ie < table.NumPolyEdges(); ie++) {
      if (table.EdgeHasDualIVolPoly(it,ie)) {
        const int isov0 = table.LowerIncident(it,ie);
        const int isov1 = table.UpperIncident(it,ie);

        if (table.OnUpperIsosurface(it,isov0))  {
          error.AddMessage
            ("Error. Lower incident vertex of edge ", ie,
             " is on upper isosurface.");
          throw error;
        }

        if (table.OnLowerIsosurface(it,isov1))  {
          error.AddMessage
            ("Error. Lower incident vertex of edge ", ie,
             " is on lower isosurface.");
          throw error;
        }

      }
    }

    for (int ie = 0; ie < table.NumPolyEdges(); ie++) {
      int ivolv;
      if (table.DoesLowerIsosurfaceIntersectEdge(it, ie) !=
          table.DoesLowerIsosurfaceIntersectEdge(it, ie, ivolv)) {
          error.AddMessage
            ("Error. Incompatible results from DoesLowerIsosurfaceIntersectEdge.");
          error.AddMessage("  Table entry ", it, ".  Cube edge ", ie, ".");
          throw error;
      }

      if (table.DoesUpperIsosurfaceIntersectEdge(it, ie) !=
          table.DoesUpperIsosurfaceIntersectEdge(it, ie, ivolv)) {
          error.AddMessage
            ("Error. Incompatible results from DoesUpperIsosurfaceIntersectEdge.");
          error.AddMessage("  Table entry ", it, ".  Cube edge ", ie, ".");
          throw error;
      }

      if (table.DoesLowerIsosurfaceIntersectEdge(it, ie, ivolv)) {
        if (!table.OnLowerIsosurface(it, ivolv)) {
          error.AddMessage
            ("Error. DoesLowerIsosurfaceIntersectEdge returns ivol vertex ",
             ivolv, ".");
          error.AddMessage
            ("  Vertex ", ivolv, " is not on lower isosurface.");
          throw error;
        }
      }

      if (table.DoesUpperIsosurfaceIntersectEdge(it, ie, ivolv)) {
        if (!table.OnUpperIsosurface(it, ivolv)) {
          error.AddMessage
            ("Error. DoesUpperIsosurfaceIntersectEdge returns ivol vertex ",
             ivolv, ".");
          error.AddMessage
            ("  Vertex ", ivolv, " is not on upper isosurface.");
          throw error;
        }
      }
    }

    // Check ambiguity information
    int num_ambiguous_facets = 0;
    for (int jf = 0; jf < cube_info.NumFacets(); jf++) {
      if (table.IsFacetAmbiguous(it, jf)) 
        { num_ambiguous_facets++; }
    }

    if (num_ambiguous_facets != table.NumAmbiguousFacets(it)) {
        error.AddMessage
          ("Error. Incorrect number of ambiguous facets.");
        throw error;
    }

    if (num_ambiguous_facets == 0) {
      if (table.IsAmbiguous(it)) {
        error.AddMessage
          ("Error. Configuration is ambiguous but no facets are ambiguous.");
        throw error;
      }
    }
    else {
      if (!table.IsAmbiguous(it)) {
        error.AddMessage
          ("Error. Configuration is not ambiguous but ", 
           num_ambiguous_facets, " facets are ambiguous.");
        throw error;
      }
    }
    
  }

}


void compare_table_and_double_table
(const IVOLDUAL_CUBE_TABLE & table, IVOLDUAL_CUBE_DOUBLE_TABLE & double_table,
 const bool flag_separate_neg, IJK::ERROR & error)
{
  int ioffset;

  if (flag_separate_neg == double_table.FlagSeparateNeg()) 
    { ioffset = 0; }
  else 
    { ioffset = double_table.NumNegativeTableEntries(); }

  for (int it = 0; it < table.NumTableEntries(); it++) {

    int it2 = it + ioffset;
    if (table.NumIVolVertices(it) != double_table.NumIVolVertices(it2)) {
      error.AddMessage
        ("Error.  Incorrect number of vertices in double table entry ",
         it2, ".");
      throw error;
    }

    for (int j = 0; j < table.NumIVolVertices(it); j++) {

      if (table.OnLowerIsosurface(it,j) != 
          double_table.OnLowerIsosurface(it2,j)) {
        error.AddMessage
          ("Error.  Double table entry ", it2, " ivol vertex ", j, "");
        error.AddMessage
          ("  has incorrect value for OnLowerIsosurface().");
        throw error;
      }

      if (table.OnUpperIsosurface(it,j) != 
          double_table.OnUpperIsosurface(it2,j)) {
        error.AddMessage
          ("Error.  Double table entry ", it2, " ivol vertex ", j, "");
        error.AddMessage
          ("  has incorrect value for OnUpperIsosurface().");
        throw error;
      }
    }

    for (int k = 0; k < table.NumPolyEdges(); k++) {

      if (table.EdgeHasDualIVolPoly(it,k) != 
          double_table.EdgeHasDualIVolPoly(it2,k)) {
        error.AddMessage
          ("Error.  Double table entry ", it2, " edge ", k, "");
        error.AddMessage
          ("  has incorrect value for EdgeHasDualIVolPoly().");
      }

      if (table.EdgeHasDualIVolPoly(it,k)) {

        if (table.LowerIncident(it,k) != double_table.LowerIncident(it2,k)) {
          error.AddMessage
            ("Error.  Double table entry ", it2, " edge ", k, "");
          error.AddMessage
            ("  has incorrect value ", double_table.LowerIncident(it2,k),
             " for lower incident vertex.");
          throw error;
        }

        if (table.UpperIncident(it,k) != double_table.UpperIncident(it2,k)) {
          error.AddMessage
            ("Error.  Double table entry ", it2, " edge ", k, "");
          error.AddMessage
            ("  has incorrect value ", double_table.UpperIncident(it2,k),
             " for upper incident vertex.");
          throw error;
        }
      }
    }
  }

}


void check_ivoldual_doubletable(const int dimension)
{
  IVOLDUAL_CUBE_TABLE tableA(dimension, flag_separate_neg);
  IVOLDUAL_CUBE_TABLE tableB(dimension, !flag_separate_neg);
  IVOLDUAL_CUBE_DOUBLE_TABLE double_table(dimension, flag_separate_neg);
  IJK::PROCEDURE_ERROR error("check_ivoldual_doubletable");

  if (tableA.NumTableEntries() + tableB.NumTableEntries() !=
      double_table.NumTableEntries()) {
    error.AddMessage
      ("Error.  Incorrect number of table entries ",
       double_table.NumTableEntries(), " in double table.");
    throw error;
  }

  if (double_table.NumNegativeTableEntries() !=
      double_table.NumPositiveTableEntries()) {
    error.AddMessage
      ("Error.  Incorrect number of negative/positive table entries in double table.");
    throw error;
  }

  if (double_table.NumNegativeTableEntries() +
      double_table.NumPositiveTableEntries() != 
      double_table.NumTableEntries()) {

    error.AddMessage
      ("Error.  Number of table entries ",
       double_table.NumTableEntries(), " in double table");
    error.AddMessage
      ("  does not equal the sum of the negative and positive numbers of table entries.");
    throw error;
  }

  for (int it0 = 0; it0 < double_table.NumTableEntries(); it0++) {
    int it1 = double_table.OppositeTableIndex(it0);
    int it2 = double_table.OppositeTableIndex(it1);

    if (it0 != it2) {
      error.AddMessage
        ("Error.  Incorrect values for IVOLDUAL_CUBE_DOUBLE_TABLE::OppositeTableIndex.");
      error.AddMessage
        ("  doulbe_table.OppositeTableIndex(", it0, ") = ", it1, ".");
      error.AddMessage
        ("  double_table.OppositeTableIndex(", it1, ") = ", it2, ".");
      throw error;
    }

    for (int iv = 0; iv < double_table.NumPolyVertices(); iv++) {
      if (double_table.VertexType(it0, iv) != 
          double_table.VertexType(it1, iv)) {
        error.AddMessage
          ("Error.  Incorrect IVOLDUAL_CUBE_DOUBLE_TABLE::OppositeTableIndex.");
        error.AddMessage
          ("  double_table.OppositeTableIndex(", it0, ") = ", it1, ".");
        error.AddMessage
          ("  double_table.VertexType(", it0, ",", iv, ") = ",
           double_table.VertexType(it0,iv), ".");
        error.AddMessage
          ("  double_table.VertexType(", it1, ",", iv, ") = ",
           double_table.VertexType(it1,iv), ".");
        throw error;
      }
    }

  }

  compare_table_and_double_table
    (tableA, double_table, flag_separate_neg, error);
  compare_table_and_double_table
    (tableB, double_table, !flag_separate_neg, error);
}


void check_contains_two_opposite(const ISODUAL_CUBE_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("check_contains_two_opposite") ;

  const int num_vert = table.NumPolyVertices();

  for (int it = 0; it < table.NumTableEntries(); it++) {

    int it_complement = table.NumTableEntries()-it-1;

    bool flag_ones = contains_two_opposite_ones(it, num_vert);
    bool flag_zeros = contains_two_opposite_zeros(it_complement, num_vert);

    if (flag_ones != flag_zeros) {
      error.AddMessage
        ("Programming error. Result of contains_two_opposite_ones does not");
      error.AddMessage
        ("  result of contains_two_opposite_zeros applied to complement.");
      error.AddMessage
        ("    Entry: ", it, ".  Complement: ", it_complement, ".");
      throw error;
    }
  }

}


void check_contains_two_opposite(const int dimension)
{
  ISODUAL_CUBE_TABLE table
    (dimension, flag_separate_neg, flag_opposite_vertices);

  check_contains_two_opposite(table);
}



// **************************************************
// Miscellaneous routines
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void usage_msg(std::ostream & out)
{
  out << "Usage: testdualtable [OPTIONS]" << endl;
  out << "OPTIONS:" << endl;
  out << "  [-ivol] [-double] [-no_ov] [-sep_pos] [-active] [-opposite] [-help]" << endl;
  out << "  [-no_below] [-no_above] [-no_I1] [-no_I2]" << endl;
  out << "  [-only_non_manifold] [-only_ambig_facets | -no_ambig_facets]" 
      << endl;
  out << "  [-min_num_isov <N>]" << endl;
}

void usage_error()
{
  usage_msg(cerr);
  exit(10);
}

void parse_command_line(int argc, char **argv)
{
  IJK::ERROR error;

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = argv[iarg];
    if (s == "-ivol") {
      flag_ivol = true;
    }
    else if (s == "-double") {
      flag_double = true;
    }
    else if (s == "-no_ov") {
      flag_opposite_vertices = false;
    }
    else if (s == "-sep_pos") {
      flag_separate_neg = false;
    }
    else if (s == "-opposite") {
      flag_opposite = true;
    }
    else if (s == "-dim2") {
      flag_dim2 = true;
      flag_dim3 = false;
      flag_iso_dim4 = false;
    }
    else if (s == "-dim3") {
      flag_dim2 = false;
      flag_dim3 = true;
      flag_iso_dim4 = false;
    }
    else if (s == "-dim4") {
      flag_dim2 = false;
      flag_dim3 = false;
      flag_iso_dim4 = true;
    }
    else if (s == "-only_ambig_facets") {
      flag_only_ambig_facets = true;
      flag_selective = true;
    }
    else if (s == "-no_ambig_facets") {
      flag_no_ambig_facets = true;
      flag_selective = true;
    }
    else if (s == "-non_manifold") {
      flag_only_non_manifold = true;
      flag_selective = true;
    }
    else if (s == "-no_I1") {
      flag_no_I1 = true;
      flag_selective = true;
    }
    else if (s == "-no_I2") {
      flag_no_I2 = true;
      flag_selective = true;
    }
    else if (s == "-min_num_isov") {
      int x = IJK::get_arg_int(iarg, argc, argv, error);
      min_num_isov.Set(x);
      iarg++;
      flag_selective = true;
    }
    else if (s == "-no_below") {
      flag_no_below = true;
      flag_selective = true;
    }
    else if (s == "-no_above") {
      flag_no_above = true;
      flag_selective = true;
    }
    else if (s == "-one_LU_isov") {
      flag_at_most_one_lower_isov = true;
      flag_at_most_one_upper_isov = true;
      flag_selective = true;
    }
    else if (s == "-one_LU_isov_all_sep") {
      flag_at_most_one_lower_isov = true;
      flag_at_most_one_upper_isov = true;
      flag_at_most_one_lower_isov_opposite_sep = true;
      flag_at_most_one_upper_isov_opposite_sep = true;
      flag_selective = true;
    }
    else if (s == "-only_LU_ambig_facets_one_isov") {
      flag_only_LU_ambig_facets_at_most_one_isov = true;
    }
    else if (s == "-only_LU_ambig_facets_one_isov_all_sep") {
      flag_only_LU_ambig_facets_at_most_one_isov_all_sep = true;
    }
    else if (s == "-help") {
      help();
    }
    else {
      usage_error();
    }

    iarg++;
  }

  if (iarg != argc) { usage_error(); }
}

void help()
{
  usage_msg(cout);
  cout << endl;
  cout << "testdualtable - Output dual table." << endl;
  cout << endl;
  cout << "OPTIONS:" << endl;

  cout << "-ivol:      Output interval volume table." << endl;
  cout << "-double:    Output double table.  Half the table entries have"
       << endl
       << "            surfaces/volumes separating negative cube vertices."
       << endl
       << "            Half the table entries have surface/volumes"
       << endl
       << "            separating positive cube vertices."
       << endl;
  cout << "-no_ov:     No opposite vertices.  Do not force split of diagonally opposite" << endl
       << "              cube (hypercube) vertices." << endl;
  cout << "-sep_pos:   Separate positive vertices." << endl;
  cout << "-active:    Output list of active facets." << endl;
  cout << "-opposite:  Output whether index contains opposite ones or zeros."
       << endl;
  cout << "-dim2:      Output only 2D lookup tables." << endl;
  cout << "-dim3:      Output only 3D lookup tables." << endl;
  cout << "-dim4:      Output only 4D isosurface lookup table." << endl;
  cout << "-no_below:     Do not output interval volume table entries whose vertices" << endl;
  cout << "              have negative labels." << endl;
  cout << "-no_above:     Do not output interval volume table entries whose vertices" << endl;
  cout << "              have negative labels." << endl;
  cout << "-no_I1:     Do not output interval volume table entries whose vertices" << endl;
  cout << "              have I1 labels." << endl;
  cout << "-no_I2:     Do not output interval volume table entries whose vertices" << endl;
  cout << "              have I2 labels." << endl;
  cout << "-only_ambig_facets: Output only table entries with ambig facets." 
       << endl;
  cout << "            (Currently implemented only for interval volumes.)"
       << endl;
  cout << "-no_ambig_facets: Output only table entries with no ambig facets." 
       << endl;
  cout << "            (Currently implemented only for interval volumes.)"
       << endl;
  cout << "-non_manifold: Output only non-manifold table entries."
       << endl;
  cout << "-min_num_isov <N>: Output only table entries with at least <N> vertices." 
       << endl;
  cout << "            (Currently implemented only for interval volumes.)"
       << endl;
  cout << "-help:      Print this help message." << endl;
  exit(20);
}
