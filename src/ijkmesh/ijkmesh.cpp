/*!
 *  @file ijkmesh.cpp
 *  @brief Process isosurface mesh.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2006-2024 Rephael Wenger

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

/*!
  \mainpage IJKMESH: IJK MESH PROCESSING

  IJKMESH is a program for processing meshes.  

  Operations include:
  - sorting mesh vertices and triangles in lexicographic order;
  - extracting mesh boundary polytopes;
  - merging mesh vertices with identical coordinates;
  - generating a cube or cube skeleton;
  - deleting isolated vertices;
  - deleting polytopes with only 1 or 2 vertices; 
  - removing vertices and polytopes outside a rectangular region;
  - replacing vertex coordinates;
  - translating vertex coordinates;
  - extracting polytopes and line segments incident on a specified vertex;
  - reversing orientation of simplices in a simplicial mesh;
  - converting mesh file formats;
  - adding tables of measurements to mesh output.
*/


#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <algorithm>
#include <map>
#include <utility>

#include "ijkcommand_line.tpp"
#include "ijkcoord.tpp"
#include "ijkIO.tpp"
#include "ijkprint.tpp"
#include "ijkstring.tpp"

#include "ijkmesh.h"

using namespace std;
using namespace IJK;
using namespace IJKMESH;

// types
typedef enum { EDIT, SORT, VMERGE, BOUNDARY, FACETS, CONVERT, 
               MERGE, REVERSE_ORIENT, MEASURE,
               GENCUBE, GENCUBE13, SUBSET, HELP, 
               MODIFY,        // DEPRECATED
               UNKNOWN_COMMAND } COMMAND;
const char * command_string[] =  
  {"edit", "sort", "vmerge", "boundary", "facets", "convert",
   "merge", "reverse_orient", "measure",
   "gencube", "gencube13", "subset", "help", "modify", "unknown"};


typedef enum { FIG, OFF, LINE, PLY, VTK, UNKNOWN_FILE_TYPE } FILE_TYPE;
typedef SUFFIX_TYPE_LIST<FILE_TYPE> FILE_TYPE_LIST;


typedef enum COLOR { RED, BLUE, GREEN, CYAN, MAGENTA, YELLOW, ORANGE,
                     WHITE, BLACK, RGB_COLOR }
  COLOR_NAME;

const COLOR_TYPE DEFAULT_RGBA[4] = { 0.0, 0.0, 1.0, 1.0};


// global variables
POLYMESH_TYPE polymesh;
std::vector<COORD_TYPE> vertex_coord;

int dimension(DIM3);
SET_VALUE<int> mesh_dimension(DIM2);
int output_mesh_dimension;

bool flag_simplex_file(false);
bool flag_simplex_output_file(false);
bool flag_cube_file(false);
bool flag_cube_output_file(false);
bool flag_polyfile(false);
bool flag_all_poly_have_same_num_vert(false);

/// Set only if all poly have same number of vertices.
int num_vert_per_poly(0);  


bool standard_input = false;
char * input_filename1 = NULL;
char * input_filename2 = NULL;
std::string output_filename;
FILE_TYPE input_file_type = UNKNOWN_FILE_TYPE;
FILE_TYPE output_file_type = UNKNOWN_FILE_TYPE;
FILE_TYPE_LIST file_type_list;

COMMAND command = UNKNOWN_COMMAND;
vector<COORD_TYPE> translation_coord;
vector<COORD_TYPE> scale_factor;
vector<COORD_TYPE> cube_coord;
vector<COORD_TYPE> v0_coord;
vector<int> vlist;
COLOR_TYPE rgba[4];
COLOR_NAME color_name = RGB_COLOR;

bool flag_delete_isolated_vertices = false;
bool flag_delete_poly_few_vert = false;
int min_num_poly_vert = 3;
bool flag_delete_poly_vert_duplicate = false;
bool flag_delete_poly_vert_duplicate_coord = false;
bool flag_translate_coord = false;
bool flag_scale_factor = false;
bool flag_polyline = false;
bool flag_color = false;
int first_sort_coord = 0;
std::vector<COORD_TYPE> edge_length;
bool flag_neg_orient = false;
std::vector<COORD_TYPE> minc;
std::vector<COORD_TYPE> maxc;
std::vector<int> contains_vertex_list;
std::vector<REPLACE_COORD> replace_coord_list;
const int NUM_VERT_PER_TRI(3);
const int NUM_VERT_PER_QUAD(4);
bool flag_silent = false;
COORD_TYPE long_edge_ratio(1.2);
COORD_TYPE very_long_edge_ratio(1.5);
COORD_TYPE long_quad_ratio(1.5);


/// Gencube output filename include cube coordinates.
bool flag_coord_label(true);

// Translate, scale, remove duplicates in mesh.
void edit_mesh(const int dimension, 
                 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & mesh);
COORD_TYPE compute_min_coord(const std::vector<COORD_TYPE> & vertex_coord);
COORD_TYPE compute_max_coord(const std::vector<COORD_TYPE> & vertex_coord);
int compute_total_num_poly_vert(const POLYMESH_TYPE & mesh);


// sort routines
void sort_mesh
(const int dimension, const int first_sort_coord,
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & polymesh);
void sort_mesh_vertices
(const int dimension, const int first_sort_coord,
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & polymesh);
void sort_polymesh(POLYMESH_TYPE & polymesh);

// Delete degenerate polytopes.
void delete_degenerate_poly
(int * vlist, int * num_poly_vert, int * first_poly_vert, int & num_poly);

// vertex merge routines
void merge_mesh_vertices
(const int dimension, const int first_sort_coord,
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & polymesh);


// remove duplicate edges
void remove_duplicate_edges(std::vector<int> & edge_endpoint);

// Extract routines.
void extract_boundary
(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & boundary_mesh);
void extract_facets
(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & facet_mesh);

// merge mesh routine
void merge_mesh
(const int dimension,
 const std::vector<COORD_TYPE> & vertex_coord1,
 const POLYMESH_TYPE & mesh1, 
 const std::vector<COORD_TYPE> & vertex_coord2,
 const POLYMESH_TYPE & mesh2,
 std::vector<COORD_TYPE> & vertex_coord3,
 POLYMESH_TYPE & mesh3);


// extract mesh edges
void extract_mesh_edges
(const int mesh_dimension, const POLYMESH_TYPE & mesh, 
 POLYMESH_TYPE & skeleton);

// reverse orientations
void reverse_orientations(POLYMESH_TYPE & mesh);

// replace vertex coordinates.
void replace_coord
(const int dimension, 
 const std::vector<REPLACE_COORD> & replace_coord_list,
 std::vector<COORD_TYPE> & vertex_coord);

// translate/scale mesh coordinates.
void translate_mesh
(const int dimension, const vector<COORD_TYPE> & translation_coord,
 std::vector<COORD_TYPE> & vertex_coord);
void scale_mesh
(const int dimension, const vector<COORD_TYPE> & scale_factor, 
 std::vector<COORD_TYPE> & vertex_coord);

// measure mesh cells.
void measure_mesh
(const int dimension, const int mesh_dimension,
 const std::vector<COORD_TYPE> & vertex_coord,
 const POLYMESH_TYPE & mesh, MESH_MEASUREMENTS & measurements);

// gencube routines
void gencube
(const int dimension, const vector<COORD_TYPE> & coord0, 
 const std::vector<COORD_TYPE> & edge_length, 
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh);
void gencube_skeleton
(const int dimension, const vector<COORD_TYPE> & coord0, 
 const std::vector<COORD_TYPE> & edge_length,
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh);
void gencube13
(const int dimension, const vector<COORD_TYPE> & coord0,
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh);
void gencube13_skeleton
(const int dimension, const vector<COORD_TYPE> & coord0, 
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh);


// get subset mesh
void get_subset(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & subset_mesh);

// check routines
void check_cube
(const int dimension, const std::vector<COORD_TYPE> & coord0, ERROR & error);
void check_replace_coord_usage
(const int dimension, const std::vector<REPLACE_COORD> & replace_coord);

// misc routines
bool equals(const int dimension, const COORD_TYPE * const vertex_coord, 
            const int i0, const int i1);
bool equals_edge(const int * endpoint, const int i0, const int i1);
bool equals_edge
(const std::vector<int> & endpoint, const int i0, const int i1);
bool is_degenerate(const int num_poly_vert, const int * poly_vert);
void memory_exhaustion();
void backup_file(const char * fname);
void construct_gencube_output_filename();
void parse_command_line(int argc, char **argv);
void usage_error(), help();
void split_string(const string & s, const char c,
                  string & prefix, string & suffix);
void write_color_ply
(std::ostream & out,
 const int dimension, COORD_TYPE * vertex_coord, const int numv, 
 const int * poly1_vert, const int numv_per_poly1, const int num_poly1,
 const int * poly2_vert, const int numv_per_poly2, const int num_poly2);
void write_poly_mesh
(const int dimension,
 const std::vector<COORD_TYPE> & vertex_coord, const POLYMESH_TYPE & mesh);
void write_poly_mesh
(const int dimension, 
 const COORD_TYPE * vertex_coord, const int numv, 
 const int * num_poly_vert, const int * poly_vert, const int * first_poly_vert,
 const int num_poly);
void write_line
(const int dimension, const std::vector<COORD_TYPE> & vertex_coord,
 const POLYMESH_TYPE & skeleton);
void write_poly_mesh
(const int dimension, const std::vector<COORD_TYPE> & vertex_coord, 
 const POLYMESH_TYPE & mesh, const MESH_MEASUREMENTS & measurements);


void read_poly_mesh
(const char * input_filename, int & dimension, 
 IJK::SET_VALUE<int> & mesh_dimension,
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & mesh);

void write_line
(const int dimension,
 const COORD_TYPE * vertex_coord, const int numv, 
 const int * edge_endpoint, const int nume);


/// Comparison function template.
/// Compare (array[i0], array[i0], ..., array[i0+tuple_size-1]) with
/// (array[i1], array[i1], ..., array[i1+tuple_size-1]).
template <class T> class ARRAY_LESS_THAN {
protected:
  const int tuple_size;
  const T * array;

public:
  ARRAY_LESS_THAN(const int s, const T * a):
    tuple_size(s), array(a) {};
  

  bool operator ()(const int i0, const int i1) const
  { 
    const T * p0 = i0*tuple_size+array;
    const T * p1 = i1*tuple_size+array;

    return(std::lexicographical_compare(p0, p0+tuple_size, p1, p1+tuple_size)); 
  };
};

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  std::set_new_handler(memory_exhaustion);

  try {

    parse_command_line(argc, argv);

    if (command == GENCUBE || command == GENCUBE13) {
      if (v0_coord.size() == 0) {
        v0_coord.resize(dimension, 0);
        if (cube_coord.size() == dimension) {
          for (int d = 0; d < dimension; d++) {
            v0_coord[d] = cube_coord[d]*edge_length[d];
          }
        }
      }

      if (output_filename == "") 
        { construct_gencube_output_filename(); }
    }
    else {
      read_poly_mesh
        (input_filename1, dimension, mesh_dimension, vertex_coord, polymesh);

      if (polymesh.NumPoly() == 0)
        { cerr << "Warning: No mesh polygons in input." << endl; }

    }

    check_replace_coord_usage(dimension, replace_coord_list);

    if (command == SORT) {
      if (first_sort_coord < 0 || first_sort_coord >= dimension) {
        cerr << "Usage error. Illegal first sort coordinate " 
             <<  first_sort_coord << "." << endl;
        exit(128);
      }
    }

  }
  catch(ERROR & error) {
    error.Print(cerr);
    exit(15);
  };

  try {

    // default
    output_mesh_dimension = mesh_dimension.Value();

    switch(command) {

    case MODIFY:
      cerr << "*** Warning: Command modify is deprecated." << endl;
      cerr << "    Use \"edit\" instead." << endl;

    case EDIT:
      edit_mesh(dimension, vertex_coord, polymesh);
      write_poly_mesh(dimension, vertex_coord, polymesh);
      break;

    case SORT:
      sort_mesh(dimension, first_sort_coord, vertex_coord, polymesh);
      write_poly_mesh(dimension, vertex_coord, polymesh);
      break;

    case VMERGE:
      merge_mesh_vertices(dimension, first_sort_coord, vertex_coord, polymesh);
      write_poly_mesh(dimension, vertex_coord, polymesh);
      break;
    
    case BOUNDARY: 
      {
        POLYMESH_TYPE boundary_mesh;
        extract_boundary(polymesh, boundary_mesh);
        write_poly_mesh(dimension, vertex_coord, boundary_mesh);
      }
      break;

    case FACETS: 
      {
        POLYMESH_TYPE facets_mesh;
        extract_facets(polymesh, facets_mesh);
        write_poly_mesh(dimension, vertex_coord, facets_mesh);
      }
      break;

    case CONVERT:
      output_mesh_dimension = mesh_dimension.Value();
      if (output_file_type == PLY || output_file_type == FIG ||
          output_file_type == VTK) {
        write_poly_mesh(dimension, vertex_coord, polymesh);
      }
      else if (output_file_type == LINE) {
        POLYMESH_TYPE skeleton;
        extract_mesh_edges(mesh_dimension.Value(), polymesh, skeleton);
        write_line(dimension, vertex_coord, skeleton);
      }
      else {
        const char * suffix_type_string =
          get_suffix_string
          (output_file_type, file_type_list, ".unknown");
        cerr << "Cannot convert to output file type " 
             << suffix_type_string << "." << endl;
        exit(40);
      };
      break;

    case MERGE:
      {
        int dimension2;
        IJK::SET_VALUE<int> mesh_dimension2;
        mesh_dimension2.Set(mesh_dimension.Value());
        POLYMESH_TYPE polymesh2, polymesh3;
        std::vector<COORD_TYPE> vertex_coord2, vertex_coord3;

        read_poly_mesh
          (input_filename2, dimension2, mesh_dimension2, vertex_coord2, 
           polymesh2);

        if (dimension2 != dimension) {
          cerr << "Input error. Meshes 1 and 2 lie in volumes of different dimension." << endl;
          exit(20);
        }

        if (mesh_dimension2.Value() != mesh_dimension.Value()) {
          cerr << "Input error. Meshes 1 and 2 have different dimensions."
               << endl;
          cerr << "  Mesh 1 dimension: " << mesh_dimension.Value() << endl;
          cerr << "  Mesh 2 dimension: " << mesh_dimension2.Value() << endl;
          exit(20);
        }

        merge_mesh
          (dimension, vertex_coord, polymesh, vertex_coord2, polymesh2,
           vertex_coord3, polymesh3);

        write_poly_mesh(dimension, vertex_coord3, polymesh3);
      }
      break;

    case REVERSE_ORIENT:
      {
        reverse_orientations(polymesh);
        write_poly_mesh(dimension, vertex_coord, polymesh);
        break;
      }

    case MEASURE:
      {
        if (dimension == DIM3 && flag_cube_file) {
          MESH_MEASUREMENTS measurements;
          measure_mesh
            (dimension, mesh_dimension.Value(), vertex_coord, polymesh,
             measurements);

          write_poly_mesh
            (dimension, vertex_coord, polymesh, measurements);
        }
        else {
          cerr << "Command measure only implemented for 3D hexahedral meshes." << endl;
          exit(40);
        }
        break;
      }

    case GENCUBE:
      if (output_file_type == LINE) {
        gencube_skeleton
          (dimension, v0_coord, edge_length, vertex_coord, polymesh);
        write_line(dimension, vertex_coord, polymesh);
      }
      else {
        gencube(dimension, v0_coord, edge_length, vertex_coord, polymesh);
        write_poly_mesh(dimension, vertex_coord, polymesh);
      }
      break;

    case GENCUBE13:
      if (output_file_type == LINE) {
        gencube13_skeleton(dimension, v0_coord, vertex_coord, polymesh);
        write_line(dimension, vertex_coord, polymesh);
      }
      else {
        gencube13(dimension, v0_coord, vertex_coord, polymesh);
        write_poly_mesh(dimension, vertex_coord, polymesh);
      }
      break;

    case SUBSET:
      {
        POLYMESH_TYPE subset_polymesh;
        get_subset(polymesh, subset_polymesh);
        write_poly_mesh(dimension, vertex_coord, subset_polymesh);
      }
      break;

    default:
    case UNKNOWN_COMMAND:
      cerr << "Programming error. Unknown command." << endl;
      exit(20);
    };

  }
  catch(ERROR error) {
    error.Print(cerr);
    exit(30);
  }

}


// **************************************************
// EDIT MESH
// **************************************************

void edit_mesh
(const int dimension, std::vector<COORD_TYPE> & vertex_coord, 
 POLYMESH_TYPE & mesh)
{
  IJK::PROCEDURE_ERROR error("edit_mesh");

  if (replace_coord_list.size() > 0) {
    replace_coord(dimension, replace_coord_list, vertex_coord);
  }

  if (flag_translate_coord) {
    while (translation_coord.size() < dimension) 
      { translation_coord.push_back(0); }
    translate_mesh(dimension, translation_coord, vertex_coord);
  }

  if (flag_scale_factor && scale_factor.size() > 0) {
    COORD_TYPE s = scale_factor[scale_factor.size()-1];
    while (scale_factor.size() < dimension) 
      { scale_factor.push_back(s); }
    scale_mesh(dimension, scale_factor, vertex_coord);
  }

  if (minc.size() == dimension || maxc.size() == dimension) {
    COORD_TYPE cmin = compute_min_coord(vertex_coord);
    COORD_TYPE cmax = compute_max_coord(vertex_coord);
    while (minc.size() < dimension) { minc.push_back(cmin); }
    while (maxc.size() < dimension) { maxc.push_back(cmax); }

    int numv = vertex_coord.size()/dimension;
    int num_poly = mesh.NumPoly();
    delete_vert_outside_region
      (dimension, IJK::vector2pointer(minc), IJK::vector2pointer(maxc),
       IJK::vector2pointerNC(vertex_coord), numv,
       IJK::vector2pointerNC(mesh.element),
       IJK::vector2pointerNC(mesh.list_length), 
       IJK::vector2pointerNC(mesh.first_element),
       num_poly);

    vertex_coord.resize(numv*dimension);
    mesh.list_length.resize(num_poly);
    mesh.first_element.resize(num_poly);
  }

  if (flag_delete_poly_vert_duplicate_coord) {
    int num_poly = mesh.NumPoly();
    IJK::delete_poly_vert_duplicate_coord
      (dimension, IJK::vector2pointer(vertex_coord), 
       IJK::vector2pointerNC(mesh.element), 
       IJK::vector2pointerNC(mesh.list_length), 
       IJK::vector2pointerNC(mesh.first_element), 
       num_poly);
    mesh.list_length.resize(num_poly);
    mesh.first_element.resize(num_poly);
  }
  else if (flag_delete_poly_vert_duplicate) {
    int num_poly = mesh.NumPoly();
    IJK::delete_poly_vert_duplicate
      (IJK::vector2pointerNC(mesh.element), 
       IJK::vector2pointerNC(mesh.list_length), 
       IJK::vector2pointerNC(mesh.first_element), num_poly);
    mesh.list_length.resize(num_poly);
    mesh.first_element.resize(num_poly);
  }

  if (flag_delete_poly_few_vert) {
    int num_poly = mesh.NumPoly();
    IJK::delete_poly_few_vert
      (min_num_poly_vert, 
       IJK::vector2pointerNC(mesh.element), 
       IJK::vector2pointerNC(mesh.list_length),
       IJK::vector2pointerNC(mesh.first_element), num_poly);
    mesh.list_length.resize(num_poly);
    mesh.first_element.resize(num_poly);
  }

  if (contains_vertex_list.size() > 0) {
    int num_poly = mesh.NumPoly();
    IJK::delete_poly_not_incident_on_vertices
      (contains_vertex_list,
       IJK::vector2pointerNC(mesh.element),
       IJK::vector2pointerNC(mesh.list_length), 
       IJK::vector2pointerNC(mesh.first_element),
       num_poly);
    mesh.list_length.resize(num_poly);
    mesh.first_element.resize(num_poly);
  }

  if (flag_delete_isolated_vertices) {
    int numv = vertex_coord.size()/dimension;
    int total_num_poly_vert = compute_total_num_poly_vert(mesh);

    IJK::delete_unreferenced_vertices
      (dimension, IJK::vector2pointerNC(vertex_coord), numv,
       IJK::vector2pointerNC(mesh.element), total_num_poly_vert);
    vertex_coord.resize(numv*dimension);
  }

  output_mesh_dimension = mesh_dimension.Value();
}


// **************************************************
// SORT
// **************************************************

void permute_coord
(const int dimension, const int first_sort_coord,
 COORD_TYPE * vertex_coord, const int numv);

void reverse_permute_coord
(const int dimension, const int first_sort_coord,
 COORD_TYPE * vertex_coord, const int numv);


void sort_mesh
(const int dimension, const int first_sort_coord,
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & polymesh)
{
  sort_mesh_vertices(dimension, first_sort_coord, vertex_coord, polymesh);
  sort_polymesh(polymesh);
  output_mesh_dimension = mesh_dimension.Value();
}


void sort_polymesh(POLYMESH_TYPE & polymesh)
{
  IJK::PROCEDURE_ERROR error("sort_polymesh");

  if (polymesh.NumPoly() == 0) { return; }

  if (flag_simplex_file) {
    const int nums = polymesh.NumPoly();
    sort_simplices
      (polymesh.NumPolyVert(0), nums, 
       IJK::vector2pointerNC(polymesh.element));
  }
  else {
    std::vector<int> sorted_poly;
    POLYMESH_TYPE polymesh2;

    polymesh.SortPolyVert();
    polymesh.GetSortedPolytopeIndices(sorted_poly);

    for (int i = 0; i < polymesh.NumPoly(); i++) {
      const int ipoly = sorted_poly[i];
      const int numpv = polymesh.NumPolyVert(ipoly);
      polymesh2.AddPolytopes(polymesh.VertexList(ipoly), numpv, numpv);
    }

    if (polymesh.element.size() < polymesh2.element.size()) 
      { polymesh.element.resize(polymesh2.element.size()); }

    std::copy(polymesh2.element.begin(), polymesh2.element.end(),
              polymesh.element.begin());
    std::copy(polymesh2.first_element.begin(), polymesh2.first_element.end(),
              polymesh.first_element.begin());
    std::copy(polymesh2.list_length.begin(), polymesh2.list_length.end(),
              polymesh.list_length.begin());
  }

}


// Sort points by lexicographic order of coordinates.
// Modifies polymesh vertices to reflect new point indices.
void sort_mesh_vertices
(const int dimension, const int first_sort_coord,
 COORD_TYPE * vertex_coord, const int numv,
 POLYMESH_TYPE & polymesh)
{
  permute_coord(dimension, first_sort_coord, vertex_coord, numv);

  polymesh.SortVertexCoordinates(dimension, vertex_coord, numv);
                                 
  reverse_permute_coord(dimension, first_sort_coord, vertex_coord, numv);
}

// Sort points by lexicographic order of coordinates.
// Modifies polymesh vertices to reflect new point indices.
void sort_mesh_vertices
(const int dimension, const int first_sort_coord,
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & polymesh)
{
  const int numv = vertex_coord.size()/dimension;

  sort_mesh_vertices(dimension, first_sort_coord, 
                     IJK::vector2pointerNC(vertex_coord), numv, polymesh);
}

/// Sort edges by lexicographic order of edge endpoints.
/// Does NOT retain edge orientation
void sort_edges(int * endpoint, const int nume)
{
  const int NUM_ENDPOINTS = 2;
  ARRAY_LESS_THAN<int> edge_vert_less_than(NUM_ENDPOINTS, endpoint);

  IJK::ARRAY<int> index_sorted(nume);
  IJK::ARRAY<int> index_sorted2(nume);
  
  for (int i = 0; i < nume; i++)
    { index_sorted[i] = i; }

  // sort vertices of each edge
  // does NOT preserve edge orientation
  for (int i = 0; i < nume; i++) {
    if (endpoint[2*i] > endpoint[2*i+1]) {
      std::swap(endpoint[2*i], endpoint[2*i+1]);
    }
  }

  // sort simplices by lexicographic order of vertices
  sort(index_sorted.Ptr(), index_sorted.Ptr()+nume, edge_vert_less_than);

  // reverse the permutation
  // Note: Probably could do this in place, but...
  for (int i = 0; i < nume; i++) {
    int k0 = index_sorted[i];
    index_sorted2[k0] = i;
  }

  // permute endpoint using index_sorted2
  for (int i0 = 0; i0 < nume; i0++) {
    int i1 = index_sorted2[i0];
    while (i1 != i0) {
      int * c0 = endpoint+i0*NUM_ENDPOINTS;
      int * c1 = endpoint+i1*NUM_ENDPOINTS;
      swap_ranges(c0, c0+NUM_ENDPOINTS, c1);
      swap(index_sorted2[i0], index_sorted2[i1]);
      i1 = index_sorted2[i0];
    }
  }

}

void sort_edges(std::vector<int> & endpoint)
{
  const int nume = endpoint.size()/2;

  sort_edges(IJK::vector2pointerNC(endpoint), nume);
}

/// Permute coordinates, shifting first_sort_coord to coord 0.
void permute_coord
(const int dimension, const int first_sort_coord,
 COORD_TYPE * vertex_coord, const int numv)
{
  PROCEDURE_ERROR error("permute_coord");

  if (first_sort_coord == 0) { return; }

  if (first_sort_coord < 0 || first_sort_coord >= dimension) {
    error.AddMessage("Programming error.  Illegal first sort coordinate ",
                     first_sort_coord, ".");
    error.AddMessage("  first_sort_coord must be between ", 0, " and ",
                     dimension-1, ".");
    throw error;
  }
  
  for (int i = 0; i < numv; i++) {
    COORD_TYPE * vertex_ptr = vertex_coord + i*dimension;

    for (int d = first_sort_coord; d > 0; d--) 
      { std::swap(*(vertex_ptr+d-1), *(vertex_ptr+d)); }
  }

}

void reverse_permute_coord
(const int dimension, const int first_sort_coord,
 COORD_TYPE * vertex_coord, const int numv)
{
  PROCEDURE_ERROR error("permute_coord");

  if (first_sort_coord == 0) { return; }

  if (first_sort_coord < 0 || first_sort_coord >= dimension) {
    error.AddMessage("Programming error.  Illegal first sort coordinate ",
                     first_sort_coord, ".");
    error.AddMessage("  first_sort_coord must be between ", 0, " and ",
                     dimension-1, ".");
    throw error;
  }
  
  for (int i = 0; i < numv; i++) {
    COORD_TYPE * vertex_ptr = vertex_coord + i*dimension;

    for (int d = 0; d < first_sort_coord; d++) 
      { std::swap(*(vertex_ptr+d), *(vertex_ptr+d+1)); }
  }

}


// **************************************************
// MERGE MESH VERTICES
// **************************************************

/// Sort and merge vertices with identical coordinates in mesh.
void merge_mesh_vertices
(const int dimension, const int first_sort_coord,
 COORD_TYPE * vertex_coord, int & numv, POLYMESH_TYPE & polymesh)
{
  if (numv == 0) { return; };

  sort_mesh_vertices
    (dimension, first_sort_coord, vertex_coord, numv, polymesh);

  // Merge duplicate points.
  IJK::ARRAY<int> loc(numv);
  loc[0] = 0;
  int i0 = 0;
  for (int i1 = 1; i1 < numv; i1++) {
    if (!equals(dimension, vertex_coord, i0, i1)) {
      i0++;
      if (i0 != i1) {
        copy(vertex_coord+i1*dimension,
             vertex_coord+(i1+1)*dimension,
             vertex_coord+i0*dimension);
      }
    }
    loc[i1] = i0;
  }

  numv = i0+1;

  polymesh.ReplaceElements(loc.PtrConst());

  int num_poly = polymesh.NumPoly();
  delete_degenerate_poly
    (IJK::vector2pointerNC(polymesh.element), 
     IJK::vector2pointerNC(polymesh.list_length),
     IJK::vector2pointerNC(polymesh.first_element), num_poly);

  polymesh.list_length.resize(num_poly);
  polymesh.first_element.resize(num_poly);

  sort_polymesh(polymesh);

  output_mesh_dimension = mesh_dimension.Value();
}


/// Sort and merge vertices with identical coordinates in mesh.
void merge_mesh_vertices
(const int dimension, const int first_sort_coord,
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & polymesh)
{
  int numv = vertex_coord.size()/dimension;

  merge_mesh_vertices
    (dimension, first_sort_coord, IJK::vector2pointerNC(vertex_coord), numv, 
     polymesh);

  vertex_coord.resize(numv*dimension);
}


// **************************************************
// REMOVE DUPLICATE EDGES
// **************************************************

/// Remove duplicate edges
void remove_duplicate_edges(std::vector<int> & edge_endpoint)
{
  const int nume = edge_endpoint.size()/2;

  if (nume == 0) { return; }

  sort_edges(edge_endpoint);

  // Remove duplicate edges
  int i0 = 0;
  for (int i = 1; i < nume; i++) {
    if (!equals_edge(edge_endpoint, i0, i)) {
      i0++;
      edge_endpoint[2*i0] = edge_endpoint[2*i];
      edge_endpoint[2*i0+1] = edge_endpoint[2*i+1];
    }
  }

  i0++;
  edge_endpoint.resize(i0*2);
}

// **************************************************
// EXTRACT BOUNDARY
// **************************************************

// @pre mesh is a set of k-simplices for a fixed value of k.
void extract_simplicial_boundary
(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & boundary_mesh)
{
  typedef  map<int, int, ARRAY_LESS_THAN<int> > FACET_MAP;

  const int num_simplices = mesh.NumPoly();

  boundary_mesh.Clear();

  if (num_simplices == 0) { return; }

  const int numv_per_simplex = mesh.NumPolyVert(0);

  if (numv_per_simplex < 2) { return; }

  const int numv_per_facet = numv_per_simplex-1;
  const int num_facets = num_simplices*numv_per_simplex;

  IJK::ARRAY<int> facet_vertices(num_facets*numv_per_facet);
  IJK::ARRAY<int> facet_vertices2(num_facets*numv_per_facet);

  int k = 0;
  for (int is = 0; is < num_simplices; is++) {
    for (int ifacet = 0; ifacet < numv_per_simplex; ifacet++) {
      int k0 = k;
      for (int jv = 0; jv < numv_per_simplex; jv++) {
        if (jv != ifacet) {
          facet_vertices[k] = mesh.Vertex(is,jv);
          facet_vertices2[k] = facet_vertices[k];
          k++;
        }
      }

      if (ifacet%2 == 1 && numv_per_simplex > 1) {
        // swap vertices to get positive facet orientation
        swap(facet_vertices[k0], facet_vertices[k0+1]);
      }

      // sort facet vertices
      sort(facet_vertices2.Ptr()+k0, facet_vertices2.Ptr()+k0+numv_per_facet);
    };
  };
  assert (k == num_facets*numv_per_facet);


  ARRAY_LESS_THAN<int> facet_less_than(numv_per_facet, facet_vertices2.Ptr());
  FACET_MAP facet_map(facet_less_than);

  for (int ifacet = 0; ifacet < num_facets; ifacet++) {
    FACET_MAP::iterator pos = facet_map.find(ifacet);

    if (pos != facet_map.end()) {
      pos->second++;
    }
    else {
      facet_map.insert(make_pair(ifacet, 1));
    }
  }

  for (FACET_MAP::iterator pos = facet_map.begin(); pos != facet_map.end();
       pos++) {
    if (pos->second == 1) {
      int ifacet = pos->first;
      int k0 = ifacet*numv_per_facet;
      boundary_mesh.AddPolytope(facet_vertices.Ptr()+k0, numv_per_facet);
    }
  }

}


void extract_boundary
(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & boundary_mesh)
{
  IJK::PROCEDURE_ERROR error("extract_boundary");

  if (flag_simplex_file) 
    { extract_simplicial_boundary(mesh, boundary_mesh); }
  else {
    error.AddMessage
      ("Function extract_boundary only implemented for meshes of simplices.");
    throw error;
  }
}



// **************************************************
// EXTRACT FACETS
// **************************************************

// @pre mesh is a set of k-dimensional cubes for a fixed value of k.
void extract_cube_facets
(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & facet_mesh)
{
  typedef  map<int, int, ARRAY_LESS_THAN<int> > FACET_MAP;

  const int num_cubes = mesh.NumPoly();

  facet_mesh.Clear();

  if (num_cubes == 0) { return; }

  const int numv_per_cube = mesh.NumPolyVert(0);

  if (numv_per_cube < 2) { return; }

  const int cube_dimension = 
    compute_cube_dimension_from_num_vertices(numv_per_cube);
  CUBE_TYPE cube(cube_dimension);
  const int num_facets_per_cube = cube.NumFacets();
  const int numv_per_facet = cube.NumFacetVertices();
  const int num_facets = num_cubes*num_facets_per_cube;

  IJK::ARRAY<int> facet_vertices(num_facets*numv_per_facet);
  IJK::ARRAY<int> facet_vertices2(num_facets*numv_per_facet);

  int k = 0;
  for (int icube = 0; icube < num_cubes; icube++) {
    for (int ifacet = 0; ifacet < num_facets_per_cube; ifacet++) {
      int k0 = k;
      for (int j = 0; j < numv_per_facet; j++) {
        const int jv = cube.FacetVertex(ifacet, j);
        facet_vertices[k] = mesh.Vertex(icube, jv);
        facet_vertices2[k] = facet_vertices[k];
        k++;
      }

      // sort facet vertices
      sort(facet_vertices2.Ptr()+k0, facet_vertices2.Ptr()+k0+numv_per_facet);
    };
  };
  assert (k == num_facets*numv_per_facet);

  ARRAY_LESS_THAN<int> facet_less_than(numv_per_facet, facet_vertices2.Ptr());
  FACET_MAP facet_map(facet_less_than);

  if (cube_dimension == 3) {
    // Facets are quadrilaterals.
    IJK::reorder_quad_vertices(facet_vertices.Ptr(), num_facets);
  }


  for (int ifacet = 0; ifacet < num_facets; ifacet++) {
    FACET_MAP::iterator pos = facet_map.find(ifacet);

    if (pos != facet_map.end()) {
      pos->second++;
    }
    else {
      facet_map.insert(make_pair(ifacet, 1));
      int k0 = ifacet*numv_per_facet;
      facet_mesh.AddPolytope(facet_vertices.Ptr()+k0, numv_per_facet);
    }
  }

}


void extract_facets
(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & facet_mesh)
{
  IJK::PROCEDURE_ERROR error("extract_facets");

  if (flag_cube_file) 
    { extract_cube_facets(mesh, facet_mesh); }
  else {
    error.AddMessage
      ("Function extract_facets only implemented for meshes of cubes.");
    throw error;
  }

  output_mesh_dimension = mesh_dimension.Value()-1;
}


// **************************************************
// DELETE DEGENERATE POLYTOPES
// **************************************************

void delete_degenerate_poly
(int * vlist, int * num_poly_vert, int * first_poly_vert, int & num_poly)
{
  int vlist_length = 0;
  int k = 0;
  for (int ipoly = 0; ipoly < num_poly; ipoly++) {
    if (!is_degenerate(num_poly_vert[ipoly], vlist + first_poly_vert[ipoly])) {
      if (k < ipoly) {
        int * vlist_ptr0 = vlist + first_poly_vert[ipoly];
        int * vlist_ptr1 = vlist + vlist_length;
        std::copy(vlist_ptr0, vlist_ptr0+num_poly_vert[ipoly], vlist_ptr1);
        num_poly_vert[k] = num_poly_vert[ipoly];
        first_poly_vert[k] = vlist_length;
      }
      vlist_length += num_poly_vert[k];
      k++;
    }
  }
  num_poly = k;
}


// **************************************************
// MERGE MESHES
// **************************************************

void merge_mesh
(const int dimension,
 const std::vector<COORD_TYPE> & vertex_coord1,
 const POLYMESH_TYPE & mesh1, 
 const std::vector<COORD_TYPE> & vertex_coord2,
 const POLYMESH_TYPE & mesh2,
 std::vector<COORD_TYPE> & vertex_coord3,
 POLYMESH_TYPE & mesh3)
{
  const int numv1 = vertex_coord1.size()/dimension;
  std::vector<int> poly_vert;

  vertex_coord3.resize(vertex_coord1.size()+vertex_coord2.size());

  // Copy mesh1 vertices
  std::copy(vertex_coord1.begin(), vertex_coord1.end(),
            vertex_coord3.begin());

  // Copy mesh2 vertices
  std::copy(vertex_coord2.begin(), vertex_coord2.end(),
            vertex_coord3.begin() + vertex_coord1.size());

  mesh3.Clear();

  // Copy mesh1 simplices
  mesh3.Copy(mesh1);

  // Copy mesh2 simplices.  Revise vertex indices to match index in mesh3.
  for (int ipoly = 0; ipoly < mesh2.NumPoly(); ipoly++) {
    poly_vert.clear();
    for (int k = 0; k < mesh2.NumPolyVert(ipoly); k++) 
      { poly_vert.push_back(mesh2.Vertex(ipoly, k) + numv1); }
    mesh3.AddPolytope(poly_vert);
  }

}


// **************************************************
// EXTRACT MESH EDGES
// **************************************************

void extract_simplicial_mesh_edges
(const int mesh_dimension, const POLYMESH_TYPE & mesh, 
 POLYMESH_TYPE & skeleton)
{
  const int nums = mesh.NumPoly();
  const int numv_per_simplex = mesh_dimension+1;
  const int nume_per_simplex = numv_per_simplex*(numv_per_simplex-1)/2;
  const int nume = nums*nume_per_simplex;
  const int NUM_ENDPOINTS_PER_EDGE = 2;
  int endpoint[NUM_ENDPOINTS_PER_EDGE];
  IJK::PROCEDURE_ERROR error("extract_simplicial_mesh_edges");

  skeleton.Clear();

  if (mesh.NumPoly() == 0) { return; }

  if (mesh.NumPolyVert(0) != numv_per_simplex) {
    error.AddMessage("Programming error.  Incorrect mesh dimension.");
    error.AddMessage
      ("  Mesh dimension is ", mesh_dimension, " but simplex 0 has ",
       mesh.NumPolyVert(0), " vertices.");
    error.AddMessage
      ("  Simplices in mesh with dimension ", mesh_dimension, 
       " have ", numv_per_simplex, " vertices.");
    throw error;
  }

  // List of edge endpoints.
  std::vector<int> edge_endpoint;

  // Extract mesh edges (allowing duplicates)
  int k = 0;
  for (int isimp = 0; isimp < nums; isimp++) {
    for (int j0 = 0; j0 < numv_per_simplex; j0++) {
      endpoint[0] = mesh.Vertex(isimp, j0);
      for (int j1 = j0+1; j1 < numv_per_simplex; j1++) {
        endpoint[1] = mesh.Vertex(isimp, j1);
        edge_endpoint.push_back(endpoint[0]);
        edge_endpoint.push_back(endpoint[1]);
      }
    }
  }

  // Remove duplicate edges.
  // Does NOT preserve edge orientation.
  remove_duplicate_edges(edge_endpoint);

  skeleton.AddPolytopes(edge_endpoint, NUM_ENDPOINTS_PER_EDGE);
}


void extract_mesh_edges
(const int mesh_dimension, const POLYMESH_TYPE & mesh, 
 POLYMESH_TYPE & skeleton)
{
  IJK::PROCEDURE_ERROR error("extract_mesh_edges");

  if (flag_simplex_file) 
    { extract_simplicial_mesh_edges(mesh_dimension, mesh, skeleton); }
  else {
    error.AddMessage
      ("Function extract_facets only implemented for meshes of simplices.");
    throw error;
  }
}


// **************************************************
// REPLACE COORDINATES
// **************************************************

void replace_coord
(const int dimension, 
 const REPLACE_COORD & replace_coordX,
 std::vector<COORD_TYPE> & vertex_coord)
{
  typedef typename std::vector<COORD_TYPE>::size_type SIZE_TYPE;

  const SIZE_TYPE numv = vertex_coord.size()/dimension;
  const int icoord = replace_coordX.icoord;
  IJK::PROCEDURE_ERROR error("replace_coord");

  
  if ((icoord < 0) || (icoord >= dimension)) {
    error.AddMessage
      ("Programming error. Illegal value for replace_coordX.icoord.");
    error.AddMessage
      ("  replace_coordX.icoord should be in range [0,", 
       dimension-1, "].");
    error.AddMessage
      ("  replace_coordX.icoord = ", icoord, "");
    throw error;
  }

  for (SIZE_TYPE iv = 0; iv < numv; iv++) {
    const SIZE_TYPE jc = iv*dimension + icoord;
    if (vertex_coord[jc] == replace_coordX.old_coord) {
      vertex_coord[jc] = replace_coordX.new_coord;
    }
  }

}


void replace_coord
(const int dimension, 
 const std::vector<REPLACE_COORD> & replace_coord_list,
 std::vector<COORD_TYPE> & vertex_coord)
{
  typedef typename std::vector<REPLACE_COORD>::size_type SIZE_TYPE;

  for (SIZE_TYPE i = 0; i < replace_coord_list.size(); i++) {
    replace_coord(dimension, replace_coord_list[i], vertex_coord);
  }
}

// **************************************************
// TRANSLATE MESH
// **************************************************

void translate_mesh
(const int dimension, const vector<COORD_TYPE> & translation_coord,
 std::vector<COORD_TYPE> & vertex_coord)
{
  const int numv = vertex_coord.size()/dimension;
  PROCEDURE_ERROR error("translate_mesh");

  if (translation_coord.size() != dimension) {
    error.AddMessage
      ("Programming error. Wrong number of translation coordinates.");
    throw error;
  }

  for (int i = 0; i < numv; i++) {
    for (int d = 0; d < dimension; d++) {
      vertex_coord[i*dimension+d] += translation_coord[d];
    }
  }
}


// **************************************************
// SCALE MESH
// **************************************************

void scale_mesh
(const int dimension, const vector<COORD_TYPE> & scale_factor, 
 std::vector<COORD_TYPE> & vertex_coord)
{
  const int numv = vertex_coord.size()/dimension;
  PROCEDURE_ERROR error("scale_mesh");

  if (scale_factor.size() != dimension) {
    error.AddMessage
      ("Programming error. Wrong number of scale factors.");
    throw error;
  }

  for (int i = 0; i < numv; i++) {
    for (int d = 0; d < dimension; d++) {
      vertex_coord[i*dimension+d] *= scale_factor[d];
    }
  }
}

// **************************************************
// REVERSE_ORIENT
// **************************************************

void reverse_simplex_orientations(POLYMESH_TYPE & polymesh)
{
  IJK::reverse_simplex_orientations(mesh_dimension.Value(), polymesh.element);
}

void reverse_polygon_orientations(POLYMESH_TYPE & polymesh)
{
  IJK::reverse_orientations_polygon_list
    (polymesh.element, polymesh.list_length, polymesh.first_element);
}

void reverse_cube_orientations(POLYMESH_TYPE & polymesh)
{
  const int num_vert_per_cube_facet = 
    compute_num_cube_facet_vertices(mesh_dimension.Value());

  IJK::reverse_orientations_cube_list
    (polymesh.element, num_vert_per_cube_facet);
}


void reverse_orientations(POLYMESH_TYPE & polymesh)
{
  const int DIM2(2);
  IJK::PROCEDURE_ERROR error("reverse_orientations");

  if (mesh_dimension.Value() == DIM2) {
    reverse_polygon_orientations(polymesh);
  }
  else if (flag_simplex_file) {
    reverse_simplex_orientations(polymesh);
  }
  else if (flag_cube_file) {
    reverse_cube_orientations(polymesh);
  }
  else {
    error.AddMessage
      ("Function reverse_orientations only implemented for meshes",
       "  of dimension 2 or meshes of simplices or of cubes");
    throw error;
  }
}


// **************************************************
// GET SUBSET
// **************************************************

void get_subset(const POLYMESH_TYPE & mesh, POLYMESH_TYPE & subset_mesh)
{
  for (int ipoly = 0; ipoly < mesh.NumPoly(); ipoly++) {

    for (int j = 0; j < vlist.size(); j++) {

      int jv = vlist[j];
      if (mesh.DoesPolyContainVertex(ipoly, jv)) {
        subset_mesh.AddPolytope
          (mesh.VertexList(ipoly), mesh.NumPolyVert(ipoly));
      }
    }
  }
}


// **************************************************
// MEASURE_MESH
// **************************************************


void measure_hex_mesh
(const int dimension, const int mesh_dimension,
 const std::vector<COORD_TYPE> & vertex_coord,
 const POLYMESH_TYPE & mesh, MESH_MEASUREMENTS & measurements)
{
  const int num_poly = mesh.NumPoly();
  const CUBE_TYPE cube(DIM3);
  const COORD_TYPE max_small_magnitude(0.0);
  COORD_TYPE min_jdet, max_jdet;
  COORD_TYPE min_shape, max_shape;
  int num_determinants, num_shape;
  int hex_orientation = 1;

  measurements.scaled_jacobian_determinant.resize(num_poly);
  measurements.jacobian_shape.resize(num_poly);

  if (flag_neg_orient)
    { hex_orientation = -1; }

  for (int ipoly = 0; ipoly < num_poly; ipoly++) {

    compute_min_max_hexahedron_normalized_Jacobian_determinant_3D
      (polymesh.VertexList(ipoly), hex_orientation, 
       vector2pointer(vertex_coord), cube, max_small_magnitude, 
       min_jdet, max_jdet, num_determinants);

    measurements.scaled_jacobian_determinant[ipoly] = min_jdet;

    compute_min_max_hexahedron_Jacobian_shape_3D
      (polymesh.VertexList(ipoly), vector2pointer(vertex_coord), 
       cube, max_small_magnitude, 
       min_shape, max_shape, num_shape);

    measurements.jacobian_shape[ipoly] = min_shape;
  }

}


void measure_mesh
(const int dimension, const int mesh_dimension,
 const std::vector<COORD_TYPE> & vertex_coord,
 const POLYMESH_TYPE & mesh, MESH_MEASUREMENTS & measurements)
{
  IJK::PROCEDURE_ERROR error("measure_mesh");

  if (mesh_dimension == DIM3 && flag_cube_file) {
    measure_hex_mesh
      (dimension, mesh_dimension, vertex_coord, mesh, measurements);
  }
  else {
    error.AddMessage
      ("Programming error.  Measure only implemented for 3D hexahedral mesh.");
    throw error;
  }

  flag_simplex_output_file = true;
  flag_cube_output_file = false;
}



// **************************************************
// GENERATE CUBE
// **************************************************

void gencube
(const int dimension, const vector<COORD_TYPE> & coord0, 
 const std::vector<COORD_TYPE> & edge_length, 
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh)
{
  HRECT<int,int,COORD_TYPE,COORD_TYPE> hrect(dimension);
  CUBE_TYPE cube_face_info(dimension);
  const int num_quad = hrect.NumFacets();
  const int num_vert = hrect.NumVertices();
  PROCEDURE_ERROR error("gencube");

  if (dimension < 3) {
    gencube_skeleton(dimension, coord0, edge_length, vertex_coord, mesh); 
    return;
  }

  check_cube(dimension, coord0, error);

  // dimension = 3;
  mesh_dimension.Set(dimension-1);
  output_mesh_dimension = mesh_dimension.Value();

  hrect.SetVertexCoord(IJK::vector2pointer(coord0), edge_length);

  vertex_coord.resize(hrect.NumVertices()*dimension);

  std::copy(hrect.VertexCoord(), hrect.VertexCoord()+num_vert*dimension,
            vertex_coord.begin());
  mesh.AddPolytopes
    (cube_face_info.FacetVertex(), num_quad*NUM_VERT_PER_QUAD, 
     NUM_VERT_PER_QUAD);

  // Reorder quad vertices to order around quad.
  IJK::reorder_quad_vertices(mesh.element);
}


void gencube_skeleton
(const int dimension, const vector<COORD_TYPE> & coord0, 
 const std::vector<COORD_TYPE> & edge_length,
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh)
{
  const int DIM1(1);
  const int NUM_ENDPOINTS_PER_EDGE(2);
  HRECT<int,int,COORD_TYPE,COORD_TYPE> hrect(dimension);
  CUBE_TYPE cube_face_info(dimension);
  const int num_edges = cube_face_info.NumEdges();
  PROCEDURE_ERROR error("gencube_skeleton");

  check_cube(dimension, coord0, error);

  mesh_dimension.Set(DIM1);
  output_mesh_dimension = mesh_dimension.Value();

  hrect.SetVertexCoord(IJK::vector2pointer(coord0), edge_length);

  vertex_coord.resize(hrect.NumVertices()*dimension);

  std::copy(hrect.VertexCoord(), 
            hrect.VertexCoord()+hrect.NumVertices()*dimension,
            vertex_coord.begin());
  mesh.AddPolytopes
    (cube_face_info.EdgeEndpoint(), num_edges*NUM_ENDPOINTS_PER_EDGE, 
     NUM_ENDPOINTS_PER_EDGE);

  flag_all_poly_have_same_num_vert = true;
  num_vert_per_poly = 2;
}


void gencube13
(const int dimension, const vector<COORD_TYPE> & coord0,
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh)
{
  const COORD_TYPE LENGTH3(3);
  vector<COORD_TYPE> coord1(dimension);
  HRECT<int,int,COORD_TYPE,COORD_TYPE> hrect(dimension);
  CUBE_TYPE cube_face_info(dimension);
  const int numv = hrect.NumVertices();
  const int num_quad = hrect.NumFacets();
  vector<COORD_TYPE> edge_lengthx3(edge_length.size());
  PROCEDURE_ERROR error("gencube13");

  if (dimension < 3) {
    gencube13_skeleton(dimension, coord0, vertex_coord, mesh); 
    return;
  }

  check_cube(dimension, coord0, error);

  // dimension = 3;
  mesh_dimension.Set(dimension-1);
  output_mesh_dimension = mesh_dimension.Value();

  for (int d = 0; d < dimension; d++) {
    coord1[d] = coord0[d]-edge_length[d]; 
    edge_lengthx3[d] = edge_length[d]*LENGTH3;
  }

  hrect.SetVertexCoord(IJK::vector2pointer(coord0), edge_length);

  vertex_coord.resize(2*numv*dimension);

  // Store first cube.
  std::copy(hrect.VertexCoord(), hrect.VertexCoord()+numv*dimension,
            vertex_coord.begin());
  mesh.AddPolytopes
    (cube_face_info.FacetVertex(), num_quad*NUM_VERT_PER_QUAD, 
     NUM_VERT_PER_QUAD);

  hrect.SetVertexCoord(IJK::vector2pointer(coord1), edge_lengthx3);

  // Store cube with edge length 3.
  std::copy(hrect.VertexCoord(), 
            hrect.VertexCoord()+numv*dimension,
            vertex_coord.begin()+numv*dimension);
  int quad_vert[NUM_VERT_PER_QUAD];
  for (int iquad = 0; iquad < num_quad; iquad++) {
    for (int j = 0; j < NUM_VERT_PER_QUAD; j++) 
      { quad_vert[j] = cube_face_info.FacetVertex(iquad, j) + numv; }
    mesh.AddPolytope(quad_vert, NUM_VERT_PER_QUAD);
  }

  // Reorder quad vertices to order around quad.
  IJK::reorder_quad_vertices(mesh.element);
}


void gencube13_skeleton
(const int dimension, const vector<COORD_TYPE> & coord0, 
 std::vector<COORD_TYPE> & vertex_coord,
 POLYMESH_TYPE & mesh)
{
  const int DIM1(1);
  const COORD_TYPE LENGTH3(3);
  const int NUM_ENDPOINTS_PER_EDGE(2);
  vector<COORD_TYPE> coord1(dimension);
  HRECT<int,int,COORD_TYPE,COORD_TYPE> hrect(dimension);
  CUBE_TYPE cube_face_info(dimension);
  const int numv = hrect.NumVertices();
  const int num_edges = cube_face_info.NumEdges();
  vector<COORD_TYPE> edge_lengthx3(edge_length.size());
  PROCEDURE_ERROR error("gencube13_skeleton");

  check_cube(dimension, coord0, error);

  mesh_dimension.Set(DIM1);
  output_mesh_dimension = DIM1;

  for (int d = 0; d < dimension; d++) {
    coord1[d] = coord0[d]-edge_length[d]; 
    edge_lengthx3[d] = edge_length[d]*LENGTH3;
  }

  hrect.SetVertexCoord(IJK::vector2pointer(coord0), edge_length);

  vertex_coord.resize(2*numv*dimension);

  // Store first cube.
  std::copy(hrect.VertexCoord(), 
            hrect.VertexCoord()+hrect.NumVertices()*dimension,
            vertex_coord.begin());
  mesh.AddPolytopes
    (cube_face_info.EdgeEndpoint(), num_edges*NUM_ENDPOINTS_PER_EDGE, 
     NUM_ENDPOINTS_PER_EDGE);

  hrect.SetVertexCoord(IJK::vector2pointer(coord1), edge_lengthx3);

  // Store cube with edge length 3.
  std::copy(hrect.VertexCoord(), 
            hrect.VertexCoord()+numv*dimension,
            vertex_coord.begin()+numv*dimension);
  int endpoint[NUM_ENDPOINTS_PER_EDGE];
  for (int ie = 0; ie < num_edges; ie++) {

    endpoint[0] = cube_face_info.EdgeEndpoint(ie, 0) + numv;
    endpoint[1] = cube_face_info.EdgeEndpoint(ie, 1) + numv;
    mesh.AddPolytope(endpoint, NUM_ENDPOINTS_PER_EDGE);
  }

  flag_all_poly_have_same_num_vert = true;
  num_vert_per_poly = 2;
}


// **************************************************
// READ/WRITE MESH
// **************************************************

void read_poly_mesh
(const char * input_filename, int & dimension, 
 IJK::SET_VALUE<int> & mesh_dimension, 
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & polymesh)
{
  const int DIM1(1);
  const int DIM2(2);
  const int DIM3(3);
  const int NUMV_PER_CUBE = 8;
  IJK::PROCEDURE_ERROR error("read_poly_mesh");

  read_poly_mesh_OFF(input_filename, dimension, vertex_coord,
                     polymesh.list_length, polymesh.element,
                     polymesh.first_element);

  const int num_poly = polymesh.NumPoly();

  if (num_poly == 0) { return; }

  int min_num_vert = 
    *(min_element(polymesh.list_length.begin(), polymesh.list_length.end()));
  int max_num_vert = 
    *(max_element(polymesh.list_length.begin(), polymesh.list_length.end()));

  if (min_num_vert == max_num_vert) {

    num_vert_per_poly = min_num_vert;
    flag_all_poly_have_same_num_vert = true;

    if (!mesh_dimension.IsSet()) {

      if (dimension == DIM3 && num_vert_per_poly == 4) {
        cerr << "Unable to determine mesh dimension." << endl;
        cerr << "Input polytopes could be tetrahedra or quadrilaterals."
             << endl;
        cerr << "Use option -mesh_dim <mdim>." << endl;
        exit(20);
      }

      if (dimension == DIM2) {
        mesh_dimension.Set(DIM1);
      }
      else if (num_vert_per_poly <= dimension+1) {
        // Assume poly are all simplices.
        mesh_dimension.Set(num_vert_per_poly-1);
      }
      else if (num_vert_per_poly == NUMV_PER_CUBE) {
        mesh_dimension = DIM3;
        flag_cube_file = true;
      }
      else {
        mesh_dimension.Set(dimension);
      }
    }

    if (num_vert_per_poly <= mesh_dimension.Value()+1) {
      flag_simplex_file = true;
    }
    else if (num_vert_per_poly == NUMV_PER_CUBE) {
      flag_cube_file = true;
      flag_polyfile = true;
    }
    else {
      flag_polyfile = true;
    }
  }
  else {
    flag_polyfile = true;

    if (!mesh_dimension.IsSet()) {
      mesh_dimension = dimension-1;

      cerr << "Warning: Unable to determine mesh dimension from input file."
           << endl;
      cerr << "  Assuming mesh dimension is " << mesh_dimension.Value() 
	   << "." << endl;
      cerr << "  Use option \"-mesh_dim {mdim}\" to specify mesh dimension."
           << endl;
      cerr << endl;
    }
  }

  flag_simplex_output_file = flag_simplex_file;
  flag_cube_output_file = flag_cube_file;
}


void write_poly_mesh
(const int dimension, const std::vector<COORD_TYPE> & vertex_coord, 
 const POLYMESH_TYPE & mesh)
{
  const int numv = vertex_coord.size()/dimension;

  write_poly_mesh
    (dimension, IJK::vector2pointer(vertex_coord), numv,
     IJK::vector2pointer(mesh.list_length), 
     IJK::vector2pointer(mesh.element), 
     IJK::vector2pointer(mesh.first_element), mesh.NumPoly());
}


void write_color_ply
(std::ostream & out,
 const int dimension, COORD_TYPE * vertex_coord, const int numv, 
 const int * poly1_vert, const int numv_per_poly1, const int num_poly1,
 const int * poly2_vert, const int numv_per_poly2, const int num_poly2)
{
  std::vector<unsigned char> poly1_rgb(3*numv_per_poly1);
  std::vector<unsigned char> poly2_rgb(3*numv_per_poly2);
  unsigned char rgb255[3];

  for (int j = 0; j < 3; j++) 
    { rgb255[j] = int(rgba[j]*255); }

  for (int i1 = 0; i1 < num_poly1; i1++) {
    poly1_rgb[3*i1] = rgb255[0];
    poly1_rgb[3*i1+1] = rgb255[1];
    poly1_rgb[3*i1+2] = rgb255[2];
  }

  for (int i2 = 0; i2 < num_poly2; i2++) {
    poly2_rgb[3*i2] = rgb255[0];
    poly2_rgb[3*i2+1] = rgb255[1];
    poly2_rgb[3*i2+2] = rgb255[2];
  }

  ijkoutColorFacesPLY
    (out, dimension, vertex_coord, numv, 
     poly1_vert, numv_per_poly1, num_poly1, vector2pointer(poly1_rgb),
     poly2_vert, numv_per_poly2, num_poly2, vector2pointer(poly2_rgb));
}


void write_color_poly_mesh_ply
(std::ostream & out,
 const int dimension, const COORD_TYPE * vertex_coord, const int numv, 
 const int * num_poly_vert, const int * poly_vert, const int * first_poly_vert,
 int num_poly)
{
  std::vector<unsigned char> poly_rgb(3*num_poly);
  unsigned char rgb255[3];

  for (int j = 0; j < 3; j++) 
    { rgb255[j] = int(rgba[j]*255); }

  for (int i = 0; i < num_poly; i++) {
    poly_rgb[3*i] = rgb255[0];
    poly_rgb[3*i+1] = rgb255[1];
    poly_rgb[3*i+2] = rgb255[2];
  }

  ijkoutColorFacesPLYII
    (out, dimension, vertex_coord, numv, 
     num_poly_vert, poly_vert, first_poly_vert, num_poly,
     IJK::vector2pointer(poly_rgb));
}


template <typename OSTREAM_TYPE>
void write_poly_mesh
(OSTREAM_TYPE & out, const int dimension, 
 const COORD_TYPE * vertex_coord, const int numv, 
 const int * num_poly_vert, const int * poly_vert, 
 const int * first_poly_vert, const int num_poly)
{
  const int DIM1(1);
  const int DIM2(2);
  const int DIM3(3);

  if (output_file_type == OFF) {

    ijkoutPolytopeOFF
      (out, dimension, vertex_coord, numv, 
       num_poly_vert, poly_vert, first_poly_vert, num_poly);

  }
  else if (output_file_type == PLY) {

    if (output_mesh_dimension != DIM2 || dimension != DIM3) {
      cerr << "Usage error.  Can only write 2D meshes embedded in 3D to .ply files."
           << endl;
      cerr << "  Dimension: " << dimension << endl;
      cerr << "  Output mesh dimension: " << output_mesh_dimension << endl;
      exit(150);
    }

    if (flag_color) {
      write_color_poly_mesh_ply
        (out, dimension, vertex_coord, numv, 
         num_poly_vert, poly_vert, first_poly_vert, num_poly);
    }
    else {
      ijkoutPolytopePLY
        (out, dimension, vertex_coord, numv, 
         num_poly_vert, poly_vert, first_poly_vert, num_poly);
    }
  }
  else if (output_file_type == VTK) {
    if (dimension != DIM3) {
      cerr << "Usage error.  Can only write meshes embedded in 3D to .vtk files." << endl;
      cerr << "  Mesh is embedded in dimension: " << dimension << endl;
      exit(150);
    }

    if (output_mesh_dimension != DIM3) {
      cerr << "Usage error.  Writing to .vtk files only implemented for 3D meshes." << endl;
      cerr << "  Output is not a 3D mesh." << endl;
      exit(150);
    }

    if (flag_cube_output_file) {
      ijkoutHexahedraVTK
        (out, "Hexahedral mesh", dimension, vertex_coord, numv,
         poly_vert, num_poly, true);
    }
    else if (flag_simplex_output_file) {
      ijkoutTetrahedraVTK
        (out, "Tetrahedral mesh", dimension, vertex_coord, numv,
         poly_vert, num_poly);
    }
    else {
      cerr << "Usage error.  Writing to .vtk files only implemented for meshes of hexahedra."
           << endl;
      cerr << "  Output is not a mesh of hexahedra." << endl;
      exit(150);
    }

  }
  else if (output_file_type == FIG) {

    if (dimension != DIM2) {
      cerr << "Coordinate dimension must be 2 to convert to .fig file."
           << endl;
      exit(70);
    }
    else if (num_poly > 0 && mesh_dimension.Value() != DIM1) {
      cerr << "Mesh dimension must be 1 to convert to .fig file."
           << endl;
      exit(70);
    }

    IJK::FIG_OUTPUT_PARAM<int> fig_param;
    const COORD_TYPE max_coord = get_max_abs_array_value
      (vertex_coord, dimension*numv);
    const int draw_region_width =
      fig_param.DefaultPaperDrawRegionWidth();
    fig_param.coord_scale_factor = 
      COORD_TYPE(draw_region_width)/max_coord;
    fig_param.flag_polyline = flag_polyline;

    ijkoutPolygonFIG
      (out, dimension, vertex_coord, numv,
       num_poly_vert, poly_vert, first_poly_vert, num_poly,
       fig_param);
  }
  else {
    cerr << "Illegal output file format." << endl;
    exit(80);
  }
}


void write_poly_mesh
(const int dimension, 
 const COORD_TYPE * vertex_coord, const int numv, 
 const int * num_poly_vert, const int * poly_vert, 
 const int * first_poly_vert, const int num_poly)
{

  if (standard_input) {

    // write to standard output
    write_poly_mesh
      (cout, dimension, vertex_coord, numv, num_poly_vert, 
       poly_vert, first_poly_vert, num_poly);
  }
  else {

    if (output_filename == "") {

      string filename;

      // set output file to input file
      if (command == SORT || command == VMERGE) {
	filename = input_filename1;
	backup_file(input_filename1);
      }
      else {
	string prefix, suffix;
	filename = input_filename1;
	split_string(filename, '.', prefix, suffix);
	if (suffix == "off") { filename = prefix;  };

	if (command == CONVERT) {
	  if (output_file_type == PLY)
	    { filename = filename + ".ply"; }
	  else if (output_file_type == VTK)
	    { filename = filename + ".vtk"; }
	  else if (output_file_type == FIG)
	    { filename = filename + ".fig"; }
	  else
	    { filename = filename + ".off"; }
	}
	else {
	  filename = filename + ".bnd.off";
	}
      }

      ofstream out(filename.c_str(), ios::out);
      if (!out.good()) {
	cerr << "Unable to open output file " << filename << "." << endl;
	exit(65);
      }

      write_poly_mesh
	(out, dimension, vertex_coord, numv, num_poly_vert, 
	 poly_vert, first_poly_vert, num_poly);
    
      out.close();

      if (!flag_silent) 
        { cout << "Created file " << filename << "." << endl; }
    }
    else {
      // output_filename != ""

      ofstream out(output_filename.c_str(), ios::out);
      if (!out.good()) {
	cerr << "Unable to open output file " << output_filename << "." << endl;
	exit(65);
      };

      write_poly_mesh
	(out, dimension, vertex_coord, numv, num_poly_vert, 
	 poly_vert, first_poly_vert, num_poly);

      out.close();
    }

    if (command == GENCUBE || command == GENCUBE13) {
      if (!flag_silent) 
        { cout << "Created file " << output_filename << "." << endl; }
    }
  }

}


void write_poly_mesh
(const int dimension, 
 const COORD_TYPE * vertex_coord, const int numv, 
 const int * num_poly_vert, const int * poly_vert, 
 const int * first_poly_vert, const int num_poly,
 const MESH_MEASUREMENTS & measurements)
{

  if (output_filename == "") {

    if (standard_input) {
      cerr << "Usage error.  Command measure cannot be used with standard input." << endl;
      exit(51);
    }
    else {
      string filename;
      string prefix, suffix;

      filename = input_filename1;
      split_string(filename, '.', prefix, suffix);
      if (suffix == "off") {
        filename = prefix;
      };
      filename = filename + ".vtk";

      ofstream out(filename.c_str(), ios::out);
      if (!out.good()) {
        cerr << "Unable to open output file " << filename << "." << endl;
        exit(65);
      };

      if (!flag_silent) 
        { cout << "Writing .vtk file: " << filename << endl; }

      ijkoutHexahedraAndHexDataVTK
        (out, "Hexahedral mesh", dimension, vertex_coord, numv,
         poly_vert, num_poly, true, 
         "ScaledJacobianDet", measurements.scaled_jacobian_determinant,
         "JacobianShape", measurements.jacobian_shape);

      out.close();
    }
  }
  else {

    if (output_file_type != VTK) {
      cerr << "Usage error.  Output file type must be vtk." << endl;
      exit(66);
    }

    if (dimension != DIM3 || output_mesh_dimension != DIM3 ||
        !flag_cube_file) {
      cerr << "Usage error.  Command measure only implemented for 3D hexahedral meshes." << endl;
      exit(150);
    }

    ofstream out(output_filename.c_str(), ios::out);
    if (!out.good()) {
      cerr << "Unable to open output file " << output_filename << "." << endl;
      exit(65);
    };

    ijkoutHexahedraAndHexDataVTK
      (out, "Hexahedral mesh", dimension, vertex_coord, numv,
       poly_vert, num_poly, true, 
       "ScaledJacobianDet", measurements.scaled_jacobian_determinant,
       "JacobianShape", measurements.jacobian_shape);

    out.close();
  }

}


void write_poly_mesh
(const int dimension, const std::vector<COORD_TYPE> & vertex_coord, 
 const POLYMESH_TYPE & mesh, const MESH_MEASUREMENTS & measurements)
{
  const int numv = vertex_coord.size()/dimension;

  write_poly_mesh
    (dimension, IJK::vector2pointer(vertex_coord), numv,
     IJK::vector2pointer(mesh.list_length), 
     IJK::vector2pointer(mesh.element), 
     IJK::vector2pointer(mesh.first_element), mesh.NumPoly(),
     measurements);
}


void write_line
(const int dimension,
 const COORD_TYPE * vertex_coord, const int numv, 
 const int * edge_endpoint, const int nume)
{
  if (output_file_type != LINE) {
    cerr << "Programming error.  Output file type is not line." << endl;
    exit(230);
  }

  if (output_filename != "") {

    ofstream out(output_filename.c_str(), ios::out);
    if (!out.good()) {
      cerr << "Unable to open output file " << output_filename << "." << endl;
      exit(235);
    };

    ijkoutColorLINE
      (out, dimension, vertex_coord, numv, edge_endpoint, nume, rgba);

    out.close();

    if (command == GENCUBE || command == GENCUBE13) {
      if (!flag_silent) 
        { cout << "Created file " << output_filename << "." << endl; }
    }
  }
  else {
    ijkoutLINE(cout, dimension, vertex_coord, numv, edge_endpoint, nume);
  }

}


void write_line
(const int dimension,
 const std::vector<COORD_TYPE> & vertex_coord,
 const POLYMESH_TYPE & skeleton)
{
  const int numv = vertex_coord.size()/dimension;

  write_line(dimension, IJK::vector2pointer(vertex_coord), numv,
             IJK::vector2pointer(skeleton.element), skeleton.NumPoly());
}


// **************************************************
// CLASS MESH MEMBER FUNCTIONS
// **************************************************

/// *** REPLACE by LIST_OF_LISTS::ComputeSumOfListLengths ***
int compute_total_num_poly_vert(const POLYMESH_TYPE & mesh)
{
  const int num_poly = mesh.NumPoly();
  int total_num_poly_vert = 
    std::accumulate(mesh.list_length.begin(), mesh.list_length.begin()+num_poly, 0);

  return(total_num_poly_vert);
}

COORD_TYPE compute_min_coord(const std::vector<COORD_TYPE> & vertex_coord)
{
  const int numv = vertex_coord.size()/dimension;

  if (numv == 0) { return(0); }

  COORD_TYPE min_coord = vertex_coord[0];
  for (int i = 1; i < numv*dimension; i++) {
    if (min_coord > vertex_coord[i]) 
      { min_coord = vertex_coord[i]; }
  }

  return(min_coord);
}

COORD_TYPE compute_max_coord(const std::vector<COORD_TYPE> & vertex_coord)
{
  const int numv = vertex_coord.size()/dimension;

  if (numv == 0) { return(0); }

  COORD_TYPE max_coord = vertex_coord[0];
  for (int i = 1; i < numv*dimension; i++) {
    if (max_coord < vertex_coord[i]) 
      { max_coord = vertex_coord[i]; }
  }

  return(max_coord);
}


// **************************************************
// CHECK ROUTINES
// **************************************************

void check_cube
(const int dimension, const vector<COORD_TYPE> & coord0, ERROR & error)
{
  if (dimension > 3) {
    error.AddMessage("Programming error.  Cannot generate cube for dimension ",
                     dimension, ".");
    throw error;
  }

  if (coord0.size() != dimension) {
    error.AddMessage
      ("Programming error.  Wrong number of cube coordinates for dimension ",
       dimension, ".");
    error.AddMessage("  coord0 has ", coord0.size(), " coordinates.");
    error.AddMessage("  coord0 should have ", dimension, " coordinates.");
    throw error;
  }
}


// Check usage of replace coord.
// Print usage error if usage is incorrect.
void check_replace_coord_usage
(const int dimension, const std::vector<REPLACE_COORD> & replace_coord)
{
  for (int i = 0; i < replace_coord_list.size(); i++) {
    const int icoord = replace_coord_list[i].icoord;
    if ((icoord < 0) || (icoord >= dimension)) {
      cerr << "Usage error. Incorrect coordinate <k> in \"-replace_coord <k> <old_value> <new_value>\"." << endl;
      cerr << "  Coordinate <k> should be in range [0," << dimension-1 << "]." << endl;
      cerr << "  \"-replace_coord " << icoord
	   << " " << replace_coord_list[i].old_coord
	   << " " << replace_coord_list[i].new_coord << "\""
	   << " has <k> = " << icoord << "." << endl;
      exit(52);
    }
  }
}


// **************************************************
// MISCELLANEOUS ROUTINES
// **************************************************

/// return true if all coordinates of vertex i0 equals all coordinates of i1
bool equals(const int dimension, const COORD_TYPE * const vertex_coord, 
            const int i0, const int i1) 
{
  const COORD_TYPE * const v0 = vertex_coord+dimension*i0;
  const COORD_TYPE * const v1 = vertex_coord+dimension*i1;

  return IJK::is_coord_equal(dimension, v0, v1);
}


/// return true if endpoints of vertex i0 equals endpoints of vertex i1
/// returns false if endpoints are in reverse order
bool equals_edge(const int * endpoint, const int i0, const int i1)
{
  if (endpoint[2*i0] == endpoint[2*i1] && 
      endpoint[2*i0+1] == endpoint[2*i1+1]) {
    return(true);
  }

  return(false);
}

/// return true if endpoints of vertex i0 equals endpoints of vertex i1
/// returns false if endpoints are in reverse order
bool equals_edge
(const std::vector<int> & endpoint, const int i0, const int i1)
{
  if (endpoint[2*i0] == endpoint[2*i1] && 
      endpoint[2*i0+1] == endpoint[2*i1+1]) {
    return(true);
  }

  return(false);
}


/// Return true if two polytope vertices are identical.
bool is_degenerate(const int num_poly_vert, const int * poly_vert)
{
  for (int i = 0; i+1 < num_poly_vert; i++) {
    for (int j = i+1; j < num_poly_vert; j++) {
      if (poly_vert[i] == poly_vert[j]) { return(true); }
    }
  }

  return(false);
}


void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


void backup_file(const char * fname)
{
  std::string command;
  command = "mv " + std::string(fname) + " " + std::string(fname) + "~";
  system(command.c_str());
}


COMMAND get_command_token(char * s)
  // convert string s into command token
{
  for (int i = 0; i < int(UNKNOWN_COMMAND); i++)
    if (strcmp(command_string[i], s) == 0)
      return(COMMAND(i));
  return(UNKNOWN_COMMAND);
}


int get_int(const int iarg, const int argc, char **argv)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error(); 
  }

  int x;
  if (!IJK::string2val(argv[iarg+1], x)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-integer character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  return(x);
}

float get_float(const int iarg, const int argc, char **argv)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error(); 
  }

  float x;
  if (!IJK::string2val(argv[iarg+1], x)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
    exit(50);
  }

  return(x);
}

float get_float_arg(int argc, char **argv, int iarg)
{
  if (iarg >= argc) { 
    cerr << "Usage error.  Missing float argument." << endl;
    usage_error(); 
  };
  float x;
  sscanf(argv[iarg], "%f", &x);
  return(x);
}

void get_rgb_arg(int argc, char ** argv, int iarg,
                   COLOR_TYPE * rgb)
{
  if (iarg + 4 >= argc) {
    cerr << "Usage error.  Missing arguments for "
         << argv[iarg] << "." << endl;
    usage_error();
  }

  rgb[0] = get_float_arg(argc, argv, iarg+1);
  rgb[1] = get_float_arg(argc, argv, iarg+2);
  rgb[2] = get_float_arg(argc, argv, iarg+3);
}

void set_rgba(const float r, const float g, const float b, const float a)
{
  rgba[0] = r;
  rgba[1] = g;
  rgba[2] = b;
  rgba[3] = a;
  flag_color = true;
}

void set_rgba(const float r, const float g, const float b, const float a,
              const COLOR_NAME color)
{
  rgba[0] = r;
  rgba[1] = g;
  rgba[2] = b;
  rgba[3] = a;
  flag_color = true;
  color_name = color;
}

/// Get string and convert to list of arguments.
/// Does not modify iarg.
template <typename ETYPE>
void get_multiple_arguments
(const int iarg, const int argc, char **argv, vector<ETYPE> & v)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error();
  }

  if (!IJK::string2vector(argv[iarg+1], v)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
    exit(50);
  }
}


void create_file_type_list(FILE_TYPE_LIST & file_type_list)
{
  file_type_list.push_back(make_pair(FIG, ".fig"));
  file_type_list.push_back(make_pair(OFF, ".off"));
  file_type_list.push_back(make_pair(LINE, ".line"));
  file_type_list.push_back(make_pair(PLY, ".ply"));
  file_type_list.push_back(make_pair(VTK, ".vtk"));
}


void construct_gencube_output_filename()
{
  output_filename = "cube";

  switch(color_name) {

  case RED: 
    output_filename += "R";
    break;

  case BLUE: 
    output_filename += "B";
    break;

  case GREEN: 
    output_filename += "G";
    break;

  case MAGENTA: 
    output_filename += "M";
    break;

  case CYAN: 
    output_filename += "C";
    break;

  case ORANGE: 
    output_filename += "O";
    break;

  case BLACK: 
    output_filename += "Bk";
    break;

  case WHITE: 
    output_filename += "W";
    break;
  }

  if (flag_coord_label) {
    for (int d = 0; d < cube_coord.size(); d++) {
      string coord_str;
      if (IJK::val2string(cube_coord[d], coord_str)) {
        output_filename += "." + coord_str;
      }
    }
  }

  switch (output_file_type) {

  case PLY:
    output_filename += ".ply";
    break;

  case VTK:
    output_filename += ".vtk";
    break;

  case OFF:
    output_filename += ".off";
    break;

  case LINE:
    output_filename += ".line";
    break;
  }
}

void parse_command_line(int argc, char **argv)
{
  const int DIM2(2);
  IJK::ERROR error;

  if (argc < 2) { 
    cerr << "Usage error. Missing command." << endl;
    usage_error(); 
  };

  // Set default to NULL/empty string.
  input_filename1 = NULL;
  input_filename2 = NULL;
  output_filename = "";

  // intialize rgba
  std::copy(DEFAULT_RGBA, DEFAULT_RGBA+4, rgba);

  // create file_type_list
  create_file_type_list(file_type_list);

  command = get_command_token(argv[1]);

  if (command == HELP || string(argv[1]) == "-help") 
    { help(); };

  if (command == UNKNOWN_COMMAND) {
    cerr << "Usage error. Unknown command: " << argv[1] << endl;
    usage_error(); 
  };

  int iarg = 2;
  while (iarg < argc && argv[iarg][0] == '-') {
    string option = argv[iarg];
    if (option == "-mesh_dim") {
      const int mdim = get_arg_int(iarg, argc, argv, error);
      mesh_dimension.Set(mdim);
      iarg++;
    }
    else if (option == "-vertex"){
      int iv = get_int(iarg, argc, argv);
      vlist.push_back(iv);
      iarg++;
    }
    else if (option == "-first") {

      if (command == SORT) {
        first_sort_coord = get_int(iarg, argc, argv);
        iarg++;
      }
      else {
        cerr << "Option -first cannot be used with command "
             << command_string[command] << "." << endl << endl;
        usage_error();
      }
    }
    else if (option == "-delv_isolated")
      { flag_delete_isolated_vertices = true; }
    else if (option == "-delp_ev") {
      flag_delete_poly_few_vert = true;
      min_num_poly_vert = 3;
    }
    else if (option == "-del_polyvert_dup")
      { flag_delete_poly_vert_duplicate = true; }
    else if (option == "-del_polyvert_dup_coord") 
      { flag_delete_poly_vert_duplicate_coord = true; }
    else if (option == "-minc") {
      IJK::get_arg_multiple_arguments(iarg, argc, argv, minc, error);
      iarg++;
    }
    else if (option == "-maxc") {
      IJK::get_arg_multiple_arguments(iarg, argc, argv, maxc, error);
      iarg++;
    }
    else if (option == "-containsv") {
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, contains_vertex_list, error);
      iarg++;
    }
    else if (option == "-translate") {
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, translation_coord, error);
      flag_translate_coord = true;
      iarg++;
    }
    else if (option == "-scale") {
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, scale_factor, error);
      flag_scale_factor = true;
      iarg++;
    }
    else if (option == "-cc") {
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, cube_coord, error);
      iarg++;
    }
    else if (option == "-v0coord") {
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, v0_coord, error);
      iarg++;
    }
    else if (option == "-replace_coord") {
      REPLACE_COORD replace_coordX;
      IJK::get_arg3_int_float_float
	(iarg, argc, argv, replace_coordX.icoord,
	 replace_coordX.old_coord, replace_coordX.new_coord);
      replace_coord_list.push_back(replace_coordX);
      iarg = iarg+3;
    }
    else if (option == "-o") {
      iarg++;
      if (iarg >= argc) usage_error();
      output_filename = argv[iarg];
    }
    else if (option == "-dim") {
      dimension = get_int(iarg, argc, argv);
      iarg++;
    }
    else if (option == "-elength") {
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, edge_length, error);
      iarg++;
    }
    else if (option == "-neg_orient") {
      flag_neg_orient = true;
    }
    else if (option == "-off") {
      output_file_type = OFF;
    }
    else if (option == "-ply") {
      output_file_type = PLY;
    }
    else if (option == "-fig") {
      output_file_type = FIG;
    }
    else if (option == "-vtk") {
      output_file_type = VTK;
    }
    else if (option == "-line") {
      output_file_type = LINE;
    }
    else if (option == "-polyline") {
      flag_polyline = true;
    }
    else if (option == "-rgb") {
      get_rgb_arg(argc, argv, iarg, rgba);
      iarg += 3;
    }
    else if (option == "-red") { set_rgba(1, 0, 0, 1, RED); }
    else if (option == "-green") { set_rgba(0, 1, 0, 1, GREEN); }
    else if (option == "-blue") { set_rgba(0, 0, 1, 1, BLUE); }
    else if (option == "-magenta") { set_rgba(1, 0, 1, 1, MAGENTA); }
    else if (option == "-yellow") { set_rgba(1, 1, 0, 1, YELLOW); }
    else if (option == "-cyan") { set_rgba(0, 1, 1, 1, CYAN); }
    else if (option == "-white") { set_rgba(1, 1, 1, 1, WHITE); }
    else if (option == "-black") { set_rgba(0, 0, 0, 1, BLACK); }
    else if (option == "-orange") { set_rgba(1, 0.5, 0, 1, ORANGE); }
    else if (option == "-s") { flag_silent = true; }
    else {
      cerr << "Illegal option: " << argv[iarg] << endl;
      usage_error();
    }
    iarg++;
  }

  if (command == MERGE) {
    if ((argc - iarg) != 3) {
      cerr << "Usage error: "
           << "Command merge requires two input files and one output file."
           << endl << endl;
      usage_error();
    }

    input_filename1 = argv[iarg];
    input_filename2 = argv[iarg+1];
    output_filename = argv[iarg+2];
  }
  else if (command == GENCUBE || command == GENCUBE13) {
    if ((argc-iarg) > 1) {
      cerr << "Usage error: "
           << "Command gencube has at most one output file."
           << endl << endl;
      usage_error();
    }
    else if ((argc-iarg) == 1) {
      output_filename = argv[iarg];
    }
  }
  else {

    if (standard_input) {
      if (argc != iarg) {
        cerr << "Usage error: Standard input has no input/output files."
             << endl;
        usage_error();
      }
    }
    else if (command == SORT || command == REVERSE_ORIENT ||
             command == MEASURE) {
      // Commands take one or two filenames.

      if (iarg+1 == argc) {
        input_filename1 = argv[iarg];
        if (command == SORT || command == REVERSE_ORIENT)
          { output_filename = argv[iarg]; }
      }
      else if (iarg+2 == argc) {
        input_filename1 = argv[iarg];
        output_filename = argv[iarg+1];
      }
      else if (iarg == argc) {
        cerr << "Usage error: "
             << "Command " << command_string[command]
             << " requires one input file." << endl; 
        usage_error();
      }
      else if (output_filename == "") {
        cerr << "Usage error: "
             << "Command " << command_string[command]
             << " requires one input file and one output file." << endl; 
        usage_error();
      }

    }
    else if (command == CONVERT &&
	     (output_file_type != UNKNOWN_FILE_TYPE)) {
      if ((argc - iarg) < 1) {
        cerr << "Usage error: "
             << "Command " << command_string[command]
             << " requires one input file."
             << endl << endl;
        usage_error();
      }
      else if ((argc - iarg) > 2) {
        cerr << "Usage error: "
             << "Command " << command_string[command]
             << " takes at most two file names."
             << endl << endl;
        usage_error();
      }

      input_filename1 = argv[iarg];

      if (iarg+1 < argc)
        { output_filename = argv[iarg+1]; }
    }
    else {

      if (output_filename == "") {
        if ((argc - iarg) != 2) {
          cerr << "Usage error: "
               << "Command " << command_string[command]
               << " requires one input file and one output file."
               << endl << endl;
          usage_error();
        }
      }
      else if ((argc - iarg) != 1) {
        cerr << "Usage error: "
             << "Command " << command_string[command]
             << " requires one input file."
             << endl << endl;
        usage_error();
      }

      input_filename1 = argv[iarg];

      if (iarg+1 < argc)
        { output_filename = argv[iarg+1]; }
    }
  }

  if (!mesh_dimension.IsSet()) {
    if (output_file_type == PLY) {
      // If output_file_type is set to PLY, mesh dimension should be 2.
      mesh_dimension.Set(DIM2);
    }
  }
  

  if (command == SUBSET) {
    if (vlist.size() == 0) {
      cerr << "Command " << command_string[command] 
           << " requires option -vertex <iv>." << endl;
      usage_error();
    }
  }

  input_file_type = get_suffix_type
    (input_filename1, file_type_list, UNKNOWN_FILE_TYPE);

  if (output_filename != "") {
    output_file_type = get_suffix_type
      (output_filename, file_type_list, UNKNOWN_FILE_TYPE);
  };

  if (input_file_type == FIG) {
    cerr << "Input file cannot be a .fig file." << endl << endl;
    usage_error();
  }

  if (command == CONVERT) {
    if (output_file_type == OFF) {
      cerr << "Output file cannot be a Geomview .off file." << endl << endl;
      usage_error();
    }
    else if (output_file_type == UNKNOWN_FILE_TYPE) {
      cerr << "Cannot convert to unknown output file type." << endl << endl;
      usage_error();
    }
  }

  if (command == MERGE) {
    if (input_file_type != output_file_type &&
        output_file_type != UNKNOWN_FILE_TYPE) {

      if (input_file_type != UNKNOWN_FILE_TYPE &&
          output_file_type != OFF) {

        cerr << "Usage error: Output file type must match input file type." 
             << endl;
        cerr << "Use command \"convert\" to convert file types." 
             << endl << endl;
        usage_error();
      }
    }
  }

  if (command == GENCUBE || GENCUBE13) {
    if (dimension > 3) {
      cerr << "Usage error:  Cannot generate cube for dimension > 3."
           << endl;
      exit(50);
    }
    
    if (dimension != cube_coord.size() && cube_coord.size() != 0) {
      cerr << "Usage error:  Incorrect number of cube coordinates."
           << endl;
      cerr << "  Number of cube coordinates should equal dimension ("
           << dimension << ")." << endl;
      exit(51);
    }

    if (dimension != v0_coord.size() && v0_coord.size() != 0) {
      cerr << "Usage error:  Incorrect number of vertex coordinates."
           << endl;
      cerr << "  Number of vertex coordinates should equal dimension ("
           << dimension << ")." << endl;
      exit(51);
    }

    if (edge_length.size() == 0)
      { edge_length.push_back(1); }

    if (edge_length.size() == 1) {

      const COORD_TYPE elength = edge_length[0];
      for (int d = 1; d < dimension; d++) 
        { edge_length.push_back(elength); }
    }

    if (edge_length.size() != dimension) {

      cerr << "Usage error:  Incorrect number of edge lengths."
           << endl;
      cerr << "  Number of cube coordinates should equal one or dimension ("
           << dimension << ")." << endl;
      exit(52);
    }

    if (output_file_type == UNKNOWN_FILE_TYPE)
      { output_file_type = PLY; }
  }

}

void usage_msg()
{
  cerr << "Usage: ijkmesh [COMMAND] [OPTIONS] {input file1} {input file2} {output file}" << endl;
  cerr << "COMMAND: convert || edit || sort || boundary || facets || vmerge ||"
       << endl
       << "         reverse_orient || merge || subset || measure ||"
       << endl
       << "         gencube || gencube13 || help" 
       << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "   [ -mesh_dim <mdim> ]" << endl;
  cerr << "   [ -delv_isolated ] | [ -delp_ev ]" << endl;
  cerr << "   [ -del_polyvert_dup | -del_polyvert_dup_coord ]" << endl;
  cerr << "   [ -minc \"c1 c2 ... cd\" ] | [ -maxc \"c1 c2 ... cd\" ]" << endl;
  cerr << "   [ -containsv v0 | -containsv \"v0 v1 v2 v3...\" ]" << endl;
  cerr << "   [ -replace_coord <k> <old_value> <new_value> ]" << endl;
  cerr << "   [ -translate \"c1 c2 ... cd\" ]" << endl;
  cerr << "   [ -scale \"s1 s2 ... sd\" ]" << endl;
  cerr << "   [ -cc \"x_1 x_2 ... x_d\" | -v0coord \"x_1 x_2 ... x_d\" ]"
       << endl;
  cerr << "   [ -first <coord index> | -vertex <iv> ]" << endl;
  cerr << "   [ -dim <D> | -elength <L> | -elength \"L_1 L_2 ... L_d\" ]" 
       << endl;
  cerr << "   [ -neg_orient ]" << endl;
  cerr << "   [ -rgb <R> <G> <B> ]" << endl;
  cerr << "   [ -red | -green | -blue | -magenta | -cyan | -yellow ]" << endl;
  cerr << "   [ -white | -black | -orange ]" << endl;
  cerr << "   [ -off | -ply | -fig | -vtk | -line ] [-s]" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help()
{
  cerr << "Usage: ijkmesh [COMMAND] [OPTIONS] {[input file1]} {[input file2]} [output file]" << endl;

  cerr << endl;
  cerr << "COMMAND:" << endl;
  cerr << "  convert     - Convert .off to .ply file. (Dimension must be 3.)"
       << endl;
  cerr << "                Convert .off to .vtk file. (Dimension must be 2.)" 
       << endl;
  cerr << "                Convert .off to .fig file. (Dimension must be 2.)" 
       << endl;
  cerr << "                Convert .off to .line file of mesh edges." << endl;
  cerr << "  edit        - Edit the mesh as specified by options." << endl;
  cerr << "  sort        - Sort mesh vertices in lexicographic order." << endl;
  cerr << "                Sort mesh triangles in lexicographic order." << endl;
  cerr << "  boundary    - Extract mesh boundary." << endl;
  cerr << "  facets      - Extract facets of mesh polytopes." << endl;
  cerr << "  vmerge      - Merge mesh vertices with identical coordinates." 
       << endl;
  cerr << "  reverse_orient - Reverse orientations of all simplices." << endl;
  cerr << "  merge       - Merge two geomview .off files." 
       << endl;
  cerr << "  subset      - Extract a subset of the mesh." << endl;
  cerr << "  measure     -  Add measurements to mesh output file." << endl;
  cerr << "  gencube     -  Generate a cube or cube skeleton." << endl;
  cerr << "  gencube13   -  Generate a cube or cube skeleton of edge length 1" 
       << endl
       << "                   surrounded by a cube of edge length 3." << endl;
  cerr << "  help        - Print this help message." << endl;
  cerr << endl;

  cerr << "OPTIONS:" << endl;
  cerr << "  -mesh_dim <mdim>:  Set mesh dimension to <mdim>." << endl;
  cerr << "  -delv_isolated: Delete isolated vertices.  Isolated vertices"
       << endl
       << "        are not contained in any facet." << endl;  
  cerr << "  -delp_ev: Delete polytopes with only 1 or 2 vertices." << endl;
  cerr << "  -del_polyvert_dup: Delete duplicate polytope vertices." << endl;
  cerr << "  -del_polyvert_dup_coord:" << endl;
  cerr << "        Delete polytope vertices with duplicate coordinates." 
       << endl;
  cerr << "        Note: Any polytope vertex deleted by -del_polyvert_dup"
       << endl
       << "        will also be deleted by -del_polyvert_dup_coord." << endl;
  cerr << "  -minc \"c1 c2 ... cd\": Region minimum coordinates." << endl;
  cerr << "        Remove vertices whose j'th coordinate is less than cj."
       << endl;
  cerr << "        Remove incident polytopes." << endl;
  cerr << "  -maxc \"c1 c2 ... cd\": Region maximum coordinates." << endl;
  cerr << "        Remove vertices whose j'th coordinate is greater than cj."
       << endl;
  cerr << "        Remove incident polytopes." << endl;
  cerr << "  -containsv v0:" << endl;
  cerr << "        Remove any polytope which does not contain vertex v0." 
       << endl;
  cerr << "  -containsv \"v1 v2 v3 v4...\":" << endl;
  cerr << "        Remove any polytope which does not contain any vertex in the list." << endl;
  cerr << "  -replace_coord <k> <old_value> <new_value>:" << endl;
  cerr << "        Replace <k>'th coordinate of any vertex." << endl;
  cerr << "        Replace <k>'th coordinate whose value equals <old_value>" << endl;
  cerr << "        with <new_value>." << endl;
  cerr << "        All -replace_coord commands are executed before any other commands" << endl;
  cerr << "        such as -translate or -scale." << endl;
  cerr << "  -translate \"c1 c2 ... cd\":" << endl;
  cerr << "        Translate vertex coordinates by (c1, c2, ..., cd)." << endl;
  cerr << "  -scale \"s1 s2 ... sd\":" << endl;
  cerr << "        Scale vertex coordinates by (s1, s2, ..., sd)." << endl;
  cerr << "  -first <coord index>:  Index of first coordinate used in sorting."
       << endl;
  cerr << "  -vertex <iv>:  Select polygons and line segments incident on iv."
       << endl;
  cerr << "  -dim <D> : Dimension of cube." << endl;
  cerr << "  -elength <L> :  Cube edge length." << endl;
  cerr << "  -elength \"L_1 L_2 ... L_d\":  Cube edge lengths." << endl;
  cerr << "  -cc \"x_1 x_2 ... x_d\": Grid cube coordinates." << endl;
  cerr << "  -v0coord \"x_1 x_2 ... x_d\":" << endl
       << "        Geometric coordinates of leftmost/lowest cube vertex." 
       << endl;
  cerr << "  -uniform:  Uniform triangulation." << endl;
  cerr << "  -neg_orient: Indicates that mesh elements have negative orientation."
       << endl
       << "     Useful when computing scaled Jacobian values." << endl;
  cerr << "  -rgb <R> <G> <B>: Set line color." << endl;
  cerr << "  -red | -green | -blue:  Set to given line color." << endl;
  cerr << "  -magenta | -cyan | -yellow:  Set to given line color." 
       << endl;
  cerr << "  -white | -black | -orange:  Set to given line color." << endl;
  cerr << "  -off: Output file type is off (Geomview off file format)." << endl;
  cerr << "  -ply: Output file type is ply (Stanford Polygon file format)." << endl;
  cerr << "  -fig: Output file type is fig (xfig file format)." << endl;
  exit(0);
}

void split_string(const string & s, const char c,
                  string & prefix, string & suffix)
  // split string at last occurrence of character c into prefix and suffix
{
  string::size_type i = s.rfind(c);
  if (i == string::npos) {
    prefix = s;
    suffix = "";
  }
  else {
    if (i > 0) { prefix = s.substr(0,i); }
    else { prefix = ""; };

    if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
    else { suffix = ""; };
  }
}
