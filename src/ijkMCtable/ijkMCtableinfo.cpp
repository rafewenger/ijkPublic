/*!
 *  @file ijkMCtableinfo.cpp
 *  @brief Print Marching Cubes isosurface table information.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2001-2024 Rephael Wenger

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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <cstdio> 
#include <cstdlib>
#include <cstring>


#include "ijkMCtable.h"
#include "ijkMCtable_xitIO.h"

#include "ijk.tpp"
#include "ijkbits.tpp"
#include "ijkprint.tpp"

using namespace IJK;
using namespace IJKMCUBE_TABLE;
using namespace std;

typedef float COORD_TYPE;

// global variables
const char * isotable_filename;
bool out_stat_flag = false;
bool out_poly_flag = false;
bool out_isovert_flag = false;
bool edge_table_flag = false;
vector<int> entry;

// routines
void tableinfo();
void edgetableinfo();
void read_isotable(ISOSURFACE_TABLE & isotable);
void out_poly(const ISOTABLE_POLY_BASE & poly);
void out_isosurface_vertices(const ISOSURFACE_TABLE & isotable);
void out_entry(const int it, const ISOSURFACE_TABLE & isotable);
void out_edge_entry
(const int it, const ISOSURFACE_EDGE_TABLE & iso_edge_table);
void out_stat(const ISOSURFACE_TABLE & isotable);
void out_etable_stat(const ISOSURFACE_EDGE_TABLE & iso_edge_table);
void parse_command_line(int argc, char **argv);
void help_msg(), usage_error();

int main(int argc, char **argv)
{

  try {

    parse_command_line(argc, argv);

    if (edge_table_flag) { edgetableinfo(); }
    else { tableinfo(); }

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
  };

}

// output information on isosurface table
void tableinfo()
{
  ISOSURFACE_TABLE isotable(3);

  read_isotable(isotable);

  if (out_poly_flag) { out_poly(isotable.Polytope()); };
  if (out_isovert_flag) { out_isosurface_vertices(isotable); };
  for (int i = 0; i < entry.size(); i++) 
    { out_entry(entry[i], isotable); };
  if (out_stat_flag) { out_stat(isotable); };
}


// output information on isosurface edge table
void edgetableinfo()
{
  ISOSURFACE_EDGE_TABLE iso_edge_table(3);

  read_isotable(iso_edge_table);
  iso_edge_table.GenEdgeLists();

  if (out_poly_flag) { out_poly(iso_edge_table.Polytope()); };
  if (out_isovert_flag) { out_isosurface_vertices(iso_edge_table); };
  for (int i = 0; i < entry.size(); i++) 
    { out_edge_entry(entry[i], iso_edge_table); };
  if (out_stat_flag) { out_etable_stat(iso_edge_table); };
}

//**************************************************
// READ ISOSURFACE TABLE
//**************************************************

// read isosurface lookup table
void read_isotable(ISOSURFACE_TABLE & isotable)
{
  ifstream isotable_file(isotable_filename, ios::in);

  if (!isotable_file) {
    string error_msg = 
      "Unable to open isosurface table file " + string(isotable_filename) 
      + ".";
    throw(ERROR(error_msg));
  };

  try {
    IJKXIO::read_xit(isotable_file, isotable);
  }
  catch(...) {
    cerr << "Error reading file " << isotable_filename << "." << endl;
    throw;
  };

  isotable_file.close();

  ERROR error;
  if (!isotable.Check(error)) {
    cerr << "Warning: Data structure inconsistency in isosurface table."
	 << endl;
    error.Print(cerr);
    cerr << "  Attempting to continue..." << endl << endl;
  };
}

//**************************************************
// OUTPUT ISOSURFACE TABLE POLYTOPE
//**************************************************

// output polytope
void out_poly(const ISOTABLE_POLY_BASE & poly)
{
  if (!poly.IsShapeUndefined()) {
    cout << "Polytope shape: " << poly.ShapeString() << endl;
  }
  int numv = poly.NumVertices();
  cout << "Number of polytope vertices = " << numv << "." << endl;
  cout << "Vertices:" << endl;
  for (int iv = 0; iv < numv; iv++) {
    cout << "  " << setw(3) << iv << ": ";
    poly.PrintVertexCoord(cout, iv);
    cout << endl;
  };
  cout << endl;

  int nume = poly.NumEdges();
  cout << "Number of polytope edges = " << nume << "." << endl;
  cout << "Edges: " << endl;
  for (int ie = 0; ie < nume; ie++) {
    cout << "  " << setw(3) << ie << ": ";
    poly.PrintEdgeEndpoints(cout, ie);
    cout << endl;
  };
  cout << endl;

  int numf = poly.NumFacets();
  cout << "Number of polytope facets = " << numf << "." << endl;
  cout << "Facets: " << endl;
  for (int jf = 0; jf < numf; jf++) {
    cout << "  " << setw(3) << jf << ": ";
    poly.PrintFacetVertices(cout, jf);
    cout << endl;
  }
  cout << endl;
}

//**************************************************
// OUTPUT ISOSURFACE VERTICES
//**************************************************

// output isosurface vertices
void out_isosurface_vertices(const ISOSURFACE_TABLE & isotable)
{
  cout << "Number of isosurface vertices = "
       << isotable.NumIsosurfaceVertices() << "." << endl;
  cout << "Vertices:" << endl;
  for (int i = 0; i < isotable.NumIsosurfaceVertices(); i++) {
    cout << "  " << setw(3) << i << ":  ";
    switch(isotable.IsosurfaceVertex(i).Type()) {
    case ISOSURFACE_VERTEX::VERTEX:
      cout << "Poly vertex " << isotable.IsosurfaceVertex(i).Face();
      break;

    case ISOSURFACE_VERTEX::EDGE:
      cout << "Poly edge " << isotable.IsosurfaceVertex(i).Face();
      break;

    case ISOSURFACE_VERTEX::FACET:
      cout << "Poly facet " << isotable.IsosurfaceVertex(i).Face();
      break;

    case ISOSURFACE_VERTEX::POINT:
      cout << "Point";
      if (isotable.IsosurfaceVertex(i).NumCoord() > 0) {
        cout << ".  Coordinates:";
        for (int d = 0; d < isotable.Dimension(); d++) {
          cout << " " << isotable.IsosurfaceVertex(i).Coord(d);
        }
      }

      break;
    };
    cout << ".";

    if (isotable.IsosurfaceVertex(i).IsLabelSet()) {
      cout << "  Label = \"" << isotable.IsosurfaceVertex(i).Label()
           << "\".";
    };
    cout << endl;
  };
  cout << endl;
}

//**************************************************
// OUTPUT ISOSURFACE TABLE ENTRY
//**************************************************

// return true if "it" is in range
// print error message and return false if "it" is out of range
bool CheckEntryRange(const int it, const ISOSURFACE_TABLE & isotable)
{
  if (it < 0 || it >= isotable.NumTableEntries()) {
    cerr << "Error.  Entry " << it
	 << " is not in range [0.." << isotable.NumTableEntries()-1
	 << "]." << endl;
    cerr << endl;
    return(false);
  }
  else return(true);
}


// output simplices in isosurface edge table entry
void out_simplices(const int it,
		   const ISOSURFACE_TABLE & isotable)
{
  const int num_poly_vertices = isotable.Polytope().NumVertices();
  std::vector<int> digit(num_poly_vertices);
  PROCEDURE_ERROR error("out_entry_simplices");

  if (it < 0 || it >= isotable.NumTableEntries()) {
    error.AddMessage("Table index ", it, " is out of bounds.");
    throw error;
  };

  IJK::convert2base
    (it, isotable.Base(), digit.data(), num_poly_vertices, error);
  
  cout << "  Number of simplices = "
       << isotable.NumSimplices(it) << "." << endl;
  IJK::print_list(cout, "  Polytope vertex labels: ",
                  digit.data(), num_poly_vertices, "\n");
  cout << "  Simplex vertices:" << endl;
  for (int is = 0; is < isotable.NumSimplices(it); is++) {
    cout << "  ";
    for (int iv = 0; iv < isotable.NumVerticesPerSimplex(); iv++) {
      cout << "  " << int(isotable.SimplexVertex(it, is, iv));
    };
    cout << endl;
  };
}


// output edges in isosurface edge table entry
void out_edges(const int it, 
               const ISOSURFACE_EDGE_TABLE & iso_edge_table)
{
  PROCEDURE_ERROR error("out_entry_simplices");

  if (it < 0 || it >= iso_edge_table.NumTableEntries()) {
    error.AddMessage("Table index ", it, " is out of bounds.");
    throw error;
  };

  cout << "  Number of edges = "
       << iso_edge_table.NumEdges(it) << "." << endl;
  cout << "  Edge endpoints:" << endl;
  for (int ie = 0; ie < iso_edge_table.NumEdges(it); ie++) {
    cout << "  ";
    for (int iend = 0; iend < 2; iend++) {
      cout << "  " << int(iso_edge_table.EdgeEndpoint(it, ie, iend));
    };
    cout << endl;
  };
}


// output simplices in isosurface table entry
void out_entry(const int it, const ISOSURFACE_TABLE & isotable)
{
  if (!CheckEntryRange(it, isotable)) { return; };

  cout << "Table Entry: " << it << endl;
  out_simplices(it, isotable);
  cout << endl;
}


// output isosurface edge table entry
void out_edge_entry
(const int it, const ISOSURFACE_EDGE_TABLE & iso_edge_table)
{
  if (!CheckEntryRange(it, iso_edge_table)) { return; };

  cout << "Table Entry: " << it << endl;
  out_simplices(it, iso_edge_table);
  out_edges(it, iso_edge_table);
  cout << endl;
}

//**************************************************
// OUTPUT ISOSURFACE TABLE STATISTICS
//**************************************************

// output polytope statistics
void out_poly_stat(const ISOSURFACE_TABLE & isotable)
{
  cout << "Polytope:" << endl;
  if (!isotable.Polytope().IsShapeUndefined()) {
    cout << "  Shape: " << isotable.Polytope().ShapeString() << endl;
  }
  cout << "  # Vertices = " << isotable.Polytope().NumVertices() << endl;
  cout << "  # Edges = " << isotable.Polytope().NumEdges() << endl;
  cout << "  # Facets = " << isotable.Polytope().NumFacets() << endl;
}

// output isosurface vertices statistics
void out_isovert_stat(const ISOSURFACE_TABLE & isotable)
{
  int num_isov_on_vert = 0;
  int num_isov_on_edges = 0;
  int num_isov_on_facets = 0;
  int num_isov_points = 0;
  for (int i = 0; i < isotable.NumIsosurfaceVertices(); i++) {
    switch(isotable.IsosurfaceVertex(i).Type()) {
    case ISOSURFACE_VERTEX::VERTEX:
      num_isov_on_vert++;
      break;

    case ISOSURFACE_VERTEX::EDGE:
      num_isov_on_edges++;
      break;

    case ISOSURFACE_VERTEX::FACET:
      num_isov_on_facets++;
      break;

    case ISOSURFACE_VERTEX::POINT:
      num_isov_points++;
      break;
    };
  }

  cout << "Isosurface Vertices:" << endl;
  cout << "  Total # vertices = " 
       << isotable.NumIsosurfaceVertices() << endl;
  if (num_isov_on_vert > 0)
    cout << "    # Vertices on poly vertices = " << num_isov_on_vert << endl;
  if (num_isov_on_edges > 0)
    cout << "    # Vertices on poly edges = " << num_isov_on_edges << endl;
  if (num_isov_on_facets > 0)
    cout << "    # Vertices on poly facets = " << num_isov_on_facets << endl;
  if (num_isov_points > 0)
    cout << "    # Free vertices = " << num_isov_points << endl;
}

// output isosurface simplices statistics
void out_simplices_stat(const ISOSURFACE_TABLE & isotable)
{
  long max_s = 0;
  long total_s = 0;
  for (long i = 0; i < isotable.NumTableEntries(); i++) {
    long n = isotable.NumSimplices(i);
    if (max_s < n) max_s = n;
    total_s += n;
  };
  double avg_s = ((double) total_s)/((double) isotable.NumTableEntries());

  cout << "  Max # simplices per entry = " << max_s << endl;
  cout << "  Avg # simplices per entry = " << avg_s << endl;
}

// output isosurface edges statistics
void out_edges_stat(const ISOSURFACE_EDGE_TABLE & iso_edge_table)
{
  long max_e = 0;
  long total_e = 0;
  for (long i = 0; i < iso_edge_table.NumTableEntries(); i++) {
    long n = iso_edge_table.NumEdges(i);
    if (max_e < n) max_e = n;
    total_e += n;
  };
  double avg_e = 
    ((double) total_e)/((double) iso_edge_table.NumTableEntries());

  cout << "  Max # edges per entry = " << max_e << endl;
  cout << "  Avg # edges per entry = " << avg_e << endl;
}

// output isosurface table statistics
void out_stat(const ISOSURFACE_TABLE & isotable)
{
  cout << "Dimension = " << isotable.Dimension() << endl;
  cout << "Simplex Dimension = " << isotable.SimplexDimension() << endl;
  isotable.Properties().Print(cout, "");
  out_poly_stat(isotable);
  out_isovert_stat(isotable);
  cout << "Table:" << endl;
  cout << "  # Entries = " << isotable.NumTableEntries() << endl;
  out_simplices_stat(isotable);
}

// output isosurface edge table statistics
void out_etable_stat(const ISOSURFACE_EDGE_TABLE & iso_edge_table)
{
  cout << "Dimension = " << iso_edge_table.Dimension() << endl;
  iso_edge_table.Properties().Print(cout, "");
  out_poly_stat(iso_edge_table);
  out_isovert_stat(iso_edge_table);
  cout << "Table:" << endl;
  cout << "  # Entries = " << iso_edge_table.NumTableEntries() << endl;
  out_simplices_stat(iso_edge_table);
  out_edges_stat(iso_edge_table);
}

//**************************************************
// MISC. ROUTINES
//**************************************************

void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc) {
    if (argv[iarg][0] != '-')
      break;

    if (strcmp(argv[iarg], "-poly") == 0) {
      out_poly_flag = true;
    }
    else if (strcmp(argv[iarg], "-isovert") == 0) {
      out_isovert_flag = true;
    }
    else if (strcmp(argv[iarg], "-stat") == 0) {
      out_stat_flag = true;
    }
    else if (strcmp(argv[iarg], "-isoedge") == 0) {
      edge_table_flag = true;
    }
    else if (strcmp(argv[iarg], "-entry") == 0) {
      iarg++;
      if (iarg >= argc)
        usage_error();

      int i;
      sscanf(argv[iarg], "%d", &i);
      entry.push_back(i);
    }
    else if (strcmp(argv[iarg], "-h") == 0 ||
	     strcmp(argv[iarg], "-help") == 0) {
      help_msg();
    }
    else { usage_error(); };
    iarg++;
  };

  if (iarg + 1 != argc) usage_error();

  if (entry.size() == 0 && !out_poly_flag && !out_isovert_flag) {
    // default: print statistics
    out_stat_flag = true; 
  };

  isotable_filename = argv[iarg];
}

void usage_msg()
{
  cout << "Usage: ijkMCtableinfo [-poly] [-stat] [-isovert] [-entry {num}]"
       << endl
       << "                      [-isoedge] [-h] {isotable file}" 
       << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help_msg()
{
  cout << "ijktableinfo - print ijk isosurface lookup table information."
       << endl;
  usage_msg();
  cout << "-poly:    Print out a complete description of the isosurface polytope" << endl;
  cout << "          including a list of vertex coordinates, a list of edges,"
       << endl;
  cout << "          and a list of facets." << endl;
  cout << "-stat:    Print out isosurface polytope and table statistics including" << endl;
  cout << "          number of polytope vertices, edges and facets, number of"
       << endl;
  cout << "          table entries, max number of simplices per entry and average"
       << endl;
  cout << "          number of simplices per entry.  This is the default." << endl;
  cout << "-isovert: Print out list of isosurface vertices." << endl;
  cout << "-entry <n>: Print entry n from the isosuface table. Print number of"
       << endl;
  cout << "          simplices in the entry and a list ofsimplex vertices.  Each"
       << endl;
  cout << "          number in the list is the index of an isosurface vertex."
       << endl;
  cout << "-isoedge: Print out maximum and average number of isosurface edges"
       << endl;
  cout << "          per entry." << endl;
  cout << "-h:       Print this help message (and exit)." << endl;

  exit(0);
}
