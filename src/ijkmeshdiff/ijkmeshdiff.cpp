/*!
 *  @file ijkmeshdiff.cpp
 *  @brief Report differences between two meshes.
 *  - Version 0.4.0
 */

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2008-2024 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <vector>

#include "ijk.tpp"
#include "ijkIO.tpp"
#include "ijkmesh_datastruct.tpp"
#include "ijkprint.tpp"
#include "ijkstring.tpp"

using namespace std;
using namespace IJK;

typedef float COORD_TYPE;
typedef int VERTEX_INDEX_TYPE;
typedef int POLY_INDEX_TYPE;
typedef BOX<COORD_TYPE> BOUNDING_BOX;
typedef typename IJK::POLYMESH<VERTEX_INDEX_TYPE,int>
POLYMESH_TYPE;

// global variables
char * input_filename0 = NULL;
char * input_filename1 = NULL;
COORD_TYPE min_diff(0.0);
bool flag_sort(false);
bool flag_terse(false);
int max_num_poly_out(10);
int max_num_coord_out(10);

// read routine
void read_poly_mesh
(const char * input_filename, int & dimension, 
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & mesh);

// sort routines
void sort_mesh_vert
(const int dimension, std::vector<COORD_TYPE> & coord,
 POLYMESH_TYPE & mesh);
void sort_mesh_poly(const POLYMESH_TYPE & mesh,
                    std::vector<POLY_INDEX_TYPE> & sorted_poly);

// output routines
void out_diff_dimensions(const int dim0, const int dim1, const char * label);
void out_diff_coord
(const int dimension,
 const std::vector<COORD_TYPE> & vertex_coord0,
 const std::vector<COORD_TYPE> & vertex_coord1,
 const COORD_TYPE min_diff, bool & flag_diff, int & num_diff);
bool out_polymesh_diff
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1,
 int & num_poly_diff);
bool out_sorted_polymesh_diff
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1,
 const std::vector<POLY_INDEX_TYPE> & sorted_poly0,
 const std::vector<POLY_INDEX_TYPE> & sorted_poly1,
 int & num_poly_diff);

// misc routines
void memory_exhaustion();
void parse_command_line(int argc, char **argv);
void usage_error(), help_msg();

// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  int dimension[2];
  POLYMESH_TYPE polymesh0, polymesh1;
  std::vector<COORD_TYPE> vertex_coord[2];
  int num_poly_diff, num_coord_diff;

  std::set_new_handler(memory_exhaustion);

  parse_command_line(argc, argv);

  try {

    read_poly_mesh
      (input_filename0, dimension[0], vertex_coord[0], polymesh0);
    read_poly_mesh
      (input_filename1, dimension[1], vertex_coord[1], polymesh1); 

    if (flag_sort) {
      sort_mesh_vert(dimension[0], vertex_coord[0], polymesh0);
      sort_mesh_vert(dimension[1], vertex_coord[1], polymesh1);
    }
    
    bool flag_diff(false);
    if (dimension[0] != dimension[1]) {
      out_diff_dimensions(dimension[0], dimension[1], "dimension");

      return 1;
    }
    else {
      if (flag_sort) {
        std::vector<POLY_INDEX_TYPE> sorted_poly0, sorted_poly1;
        sort_mesh_poly(polymesh0, sorted_poly0);
        sort_mesh_poly(polymesh1, sorted_poly1);
        if (out_sorted_polymesh_diff
            (polymesh0, polymesh1, sorted_poly0, sorted_poly1,
             num_poly_diff))
          { return 1; }
      }
      else if (out_polymesh_diff(polymesh0, polymesh1, num_poly_diff))
        { return 1; }

      if (num_poly_diff == 0) {
        out_diff_coord(dimension[0], vertex_coord[0], vertex_coord[1],
                       min_diff, flag_diff, num_coord_diff);
      }
    }

  }
  catch (ERROR error) {
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
  }

  if (num_poly_diff > 0 || num_coord_diff > 0)
    { return 1; }
    
  return 0;
}


// **************************************************
// READ ROUTINE
// **************************************************

void read_poly_mesh
(const char * input_filename, int & dimension, 
 std::vector<COORD_TYPE> & vertex_coord, POLYMESH_TYPE & polymesh)
{
  read_poly_mesh_OFF(input_filename, dimension, vertex_coord,
                     polymesh.list_length, polymesh.element,
                     polymesh.first_element);
}


// **************************************************
// SORT ROUTINE
// **************************************************

void sort_mesh_vert
(const int dimension, std::vector<COORD_TYPE> & coord,
 POLYMESH_TYPE & mesh)
{
  mesh.SortVertexCoordinates(dimension, coord);
  mesh.SortPolyVert();
}

void sort_mesh_poly(const POLYMESH_TYPE & mesh,
                    std::vector<POLY_INDEX_TYPE> & sorted_poly)
{
  mesh.GetSortedPolytopeIndices(sorted_poly);
}


// **************************************************
// OUTPUT ROUTINES
// **************************************************

void out_diff_dimensions(const int dim0, const int dim1, const char * label)
{
  cout << "Input files have different " << label << "s." << endl;
  cout << "File " << input_filename0 << " has " << label << " "
       << dim0 << "." << endl;
  cout << "File " << input_filename1 << " has " << label << " "
       << dim1 << "." << endl;
}


bool are_polytopes_different
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1,
 const POLY_INDEX_TYPE ipoly0, const POLY_INDEX_TYPE ipoly1)
{
  const int num_poly0_vert = polymesh0.NumPolyVert(ipoly0);
  const int num_poly1_vert = polymesh1.NumPolyVert(ipoly1);

  if (num_poly0_vert != num_poly1_vert)
    { return true; }

  for (int j = 0; j < num_poly0_vert; j++) {
    if (polymesh0.Vertex(ipoly0, j) != polymesh1.Vertex(ipoly1, j))
      { return true; }
  }

  return false;
}


bool out_num_poly_diff
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1)
{
  const int num_poly0 = polymesh0.NumPoly();
  const int num_poly1 = polymesh1.NumPoly();

  if (num_poly0 != num_poly1) {
    cout << "Input files have different numbers of polytopes." << endl;

    if (!flag_terse) {
      cout << "File " << input_filename0 << " has " << num_poly0
           << " polytopes." << endl;
      cout << "File " << input_filename1 << " has " << num_poly1
           << " polytopes." << endl;
    }

    return true;
  }

  return false;
}

void out_polytope_index
(const POLY_INDEX_TYPE sort_index, const bool flag_sort)
{
  if (flag_sort) { cout << "Sort index "; }
  else { cout << "Polytope "; }  
  cout << sort_index << ".";
}


// Output differences between polytope ipoly0 in polymesh0 and
//   polytope ipoly1 in polymesh1.
// Return true if differences found.
bool out_polytope_diff
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1,
 const POLY_INDEX_TYPE ipoly0, const POLY_INDEX_TYPE ipoly1,
 const POLY_INDEX_TYPE sort_index, const bool flag_sort)
{
  const int num_poly0_vert = polymesh0.NumPolyVert(ipoly0);
  const int num_poly1_vert = polymesh1.NumPolyVert(ipoly1);

  if (are_polytopes_different
      (polymesh0, polymesh1, ipoly0, ipoly1)) {

    if (num_poly0_vert != num_poly1_vert){
      out_polytope_index(sort_index, flag_sort);
      cout << " File 1 has " << num_poly0_vert << " vertices.";
      cout << " File 2 has " << num_poly1_vert << " vertices." << endl;
    }
    else {      
      out_polytope_index(sort_index, flag_sort);
      IJK::print_list
        (cout, " File 1: ", polymesh0.VertexList(ipoly0), num_poly0_vert, "");
      IJK::print_list
        (cout, " File 2: ", polymesh1.VertexList(ipoly1), num_poly1_vert, "\n");
    }

    return true;
  }
  else {
    return false;
  }
    
}


// Version where polytopes are not sorted.
bool out_polytope_diff
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1,
 const POLY_INDEX_TYPE ipoly)
{
  return out_polytope_diff
    (polymesh0, polymesh1, ipoly, ipoly, ipoly, false);
}


// Return true if differences found.
bool out_polymesh_diff
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1,
 int & num_poly_diff)
{
  const int num_poly0 = polymesh0.NumPoly();
  const int num_poly1 = polymesh1.NumPoly();

  num_poly_diff = 0;

  // Output message if two meshes have different number of polytopes.
  if (out_num_poly_diff(polymesh0, polymesh1))
    { return true; }

  for (int ipoly = 0; ipoly < num_poly0; ipoly++) {
    if (are_polytopes_different
        (polymesh0, polymesh1, ipoly, ipoly)) {
      
      if (!flag_terse && (num_poly_diff < max_num_poly_out)) {
        if (num_poly_diff == 0)
          { cout << "Polytope differences: " << endl; }
        out_polytope_diff(polymesh0, polymesh1, ipoly);
      }

      num_poly_diff++;
    }
  }

  if (num_poly_diff == 0) {
    
    cout << "All mesh polytopes are the same." << endl;

    return false;
  }
  else {
    if (!flag_terse && (num_poly_diff > max_num_poly_out)) {
      cout << "  ..." << endl;
      cout << "  (Additional " << num_poly_diff-max_num_poly_out
           << " polytope differences not listed.)" << endl;
    }

    cout << "Found " << num_poly_diff
         << " different polytopes." << endl;

    return true;
  }
}


// Return true if differences found.
bool out_sorted_polymesh_diff
(const POLYMESH_TYPE & polymesh0, const POLYMESH_TYPE & polymesh1,
 const std::vector<POLY_INDEX_TYPE> & sorted_poly0,
 const std::vector<POLY_INDEX_TYPE> & sorted_poly1,
 int & num_poly_diff)
{
  const int num_poly0 = polymesh0.NumPoly();
  const int num_poly1 = polymesh1.NumPoly();

  num_poly_diff = 0;

  // Output message if two meshes have different number of polytopes.
  if (out_num_poly_diff(polymesh0, polymesh1))
    { return true; }  


  for (int i = 0; i < num_poly0; i++) {
    const POLY_INDEX_TYPE ipoly0 = sorted_poly0[i];
    const POLY_INDEX_TYPE ipoly1 = sorted_poly1[i];
    if (are_polytopes_different
        (polymesh0, polymesh1, ipoly0, ipoly1)) {
      
      if (!flag_terse && (num_poly_diff < max_num_poly_out)) {
        if (num_poly_diff == 0)
          { cout << "Polytope differences: " << endl; }
        out_polytope_diff(polymesh0, polymesh1,
                          ipoly0, ipoly1, i, true);
      }

      num_poly_diff++;
    }
  }

  if (num_poly_diff == 0) {
    
    cout << "All mesh polytopes are the same." << endl;

    return false;
  }
  else {
    if (!flag_terse && (num_poly_diff > max_num_poly_out)) {
      cout << "  ..." << endl;
      cout << "  (Additional " << num_poly_diff-max_num_poly_out
           << " polytope differences not listed.)" << endl;
    }

    cout << "Found " << num_poly_diff
         << " different polytopes." << endl;

    return true;
  }
}


// @paramt CTYPE_MAX Type of max_diff. Should be double.
// @param min_diff Differences below or equal to min_diff are ignored.
// @param[out] max_diff Maximum difference found between two coordinates.
template <typename CTYPE_MAX>
bool diff_coord(const int dimension,
                const std::vector<COORD_TYPE> & vertex_coord0,
                const std::vector<COORD_TYPE> & vertex_coord1,
                const int i, 
                const COORD_TYPE min_diff,
                CTYPE_MAX & max_diff)
{
  bool flag_found_diff = false;
  
  if (min_diff == 0.0) {
    for (int d = 0; d < dimension; d++) {
      int j = i*dimension+d;
      const COORD_TYPE v0_j = vertex_coord0[j];
      const COORD_TYPE v1_j = vertex_coord1[j];
      if (v0_j != v1_j) {
        flag_found_diff = true;
        const CTYPE_MAX diff = std::abs(v0_j-v1_j);
        max_diff = std::max(max_diff, diff);
      }
    }
  }
  else {
    for (int d = 0; d < dimension; d++) {
      int j = i*dimension+d;
      const COORD_TYPE v0_j = vertex_coord0[j];
      const COORD_TYPE v1_j = vertex_coord1[j];      
      if (std::abs(v0_j-v1_j) > min_diff) {
        flag_found_diff = true;
        const CTYPE_MAX diff = std::abs(v0_j-v1_j);
        max_diff = std::max(max_diff, diff);        
      }
    }
  }

  return flag_found_diff;
}


void out_diff_coord
(const int dimension,
 const std::vector<COORD_TYPE> & vertex_coord0,
 const std::vector<COORD_TYPE> & vertex_coord1,
 const COORD_TYPE min_diff, bool & flag_diff, int & num_diff)
{
  const int numv0 = vertex_coord0.size()/dimension;
  const int numv1 = vertex_coord1.size()/dimension;
  double max_diff = 0.0;
  flag_diff = false;
  num_diff = 0;

  if (numv0 != numv1) {
    cout << "Input files have different number of vertices." << endl;

    if (!flag_terse) {
      cout << "File " << input_filename0 << " has " << numv0
           << " vertices." << endl;
      cout << "File " << input_filename1 << " has " << numv1
           << " vertices." << endl;
    }
    flag_diff = true;
  }
  else {
    for (int i = 0; i < numv0; i++) {

      if (diff_coord
          (dimension, vertex_coord0, vertex_coord1, i,
           min_diff, max_diff)) {
        flag_diff = true;

        if (!flag_terse && (num_diff < max_num_coord_out)) {

          if (num_diff == 0)
            { cout << "Vertex coordinate differences:" << endl; }

          int j = i*dimension;
          cout << "Vertex " << i << ".  File 1: ";
          print_list(cout, &(vertex_coord0[j]), dimension);
          cout << ".  File 2: ";
          print_list(cout, &(vertex_coord1[j]), dimension);
          cout << "." << endl;
        }

        num_diff++;
      }
    }
  }


  if (flag_diff) {

    if (num_diff > 0) {

      if (!flag_terse && (num_diff > max_num_coord_out)) {
        cout << "  ..." << endl;
        cout << "  (Additional " << num_diff-max_num_coord_out
           << " coordinate differences not listed.)" << endl;
      }
      
      cout << "Found " << num_diff << " coordinate differences."
           << "  (Max coord diff: " << max_diff << ")"
           << endl;
    }
  }
  else {
    if (min_diff == 0.0) {
      cout << "All coordinates are the same." << endl;
    }
    else {
      cout << "All coordinates are (approximately) the same." << endl;
    }
  }

}

// **************************************************
// MISCELLANEOUS ROUTINES
// **************************************************

void memory_exhaustion()
{
   cerr << "Error: Out of memory.  Terminating program." << endl;
   exit(10);
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

void parse_command_line(int argc, char **argv)
{
  if (argc == 1) { usage_error(); };

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = argv[iarg];

    if (s == "-help") { help_msg(); }
    else if (s == "-min_diff") {
      min_diff = get_float(iarg, argc, argv);
      iarg++;
    }
    else if (s == "-sort") {
      flag_sort = true;
    }
    else if (s == "-terse") {
      flag_terse = true;
    }
    else { usage_error(); }

    iarg++;
  };

  if (iarg+1 >= argc) {
    cerr << "Error.  Missing input file name." << endl;
    usage_error();
  };

  input_filename0 = argv[iarg]; 
  iarg++;
  input_filename1 = argv[iarg];
  iarg++;

  if (iarg != argc) { usage_error(); }
}

void usage_msg()
{
  cerr << "Usage: ijkmeshdiff [OPTIONS] {off file1} {off file2}" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  [-min_diff {D}] [-sort]" << endl;
  cerr << "  [-terse] [-help]" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help_msg()
{
  cerr << "Usage: ijkmeshdiff [OPTIONS] {file1} {file2}" << endl;
  cerr << "  -min_diff <D>: Coordinates which differ by <D> are different."
       << endl;
  cerr << "  -sort:         Sort polytope vertices in increasing order."
       << endl;
  cerr << "                 Ignores any polytope structure and/or orientation."
       << endl;
  cerr << "  -help:         Print this help message." << endl;
  exit(20);
}

