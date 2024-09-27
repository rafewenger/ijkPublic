/*!
 *  @file ijkdiffMCtable.cxx
 *  @brief Return the difference of two Marching Cubes isosurface tables.


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2023 Rephael Wenger

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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include <cassert>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cctype>

#include "ijkMCtable.h"
#include "ijkMCtable_xitIO.h"

using namespace std;
using namespace IJK;
using namespace IJKMCUBE_TABLE;

// global variables
char * fname[2];

// routines
void diff_table_polyhedra(const ISOSURFACE_TABLE isotable[2]);
void parse_command_line(int argc, char **argv);
bool equals(const int * vlist0, const int * vlist1, const int numv);
void sort_simplex_vertices(const ISOSURFACE_TABLE & isotable, 
			   const int it, const int is, int * vlist);

void usage_error();
void help_msg();

// create version of ISOSURFACE_TABLE with no constructor parameters
class ISOSURFACE_TABLE2:public ISOSURFACE_TABLE {
  public:
    ISOSURFACE_TABLE2():ISOSURFACE_TABLE(3) {};
};

int main(int argc, char **argv)
{
  char * isotable_filename = NULL;
  int * vlist[2] = {NULL, NULL};

  try {
    ISOSURFACE_TABLE2 isotable[2];
    ifstream in[2];

    parse_command_line(argc, argv);

    for (int i = 0; i < 2; i++) {
      isotable_filename = fname[i];
      in[i].open(isotable_filename);
      if (!in[i]) {
	cerr << "Unable to open isosurface table file "
	     << fname[i] << "." << endl;
	exit(30);
      };

      try {
	IJKXIO::read_xit(in[i], isotable[i]);
      }
      catch(...) {
	cerr << "Error reading file " << fname[i] << "." << endl;
	throw;
      };

      ERROR error;
      if (!isotable[i].Check(error)) {
	cerr << "Warning: Data structure inconsistency in isosurface table "
	     << fname[i] << "." << endl;
	error.Print(cerr);
	cerr << "  Attempting to continue..." << endl << endl;
      };
    };

    if (isotable[0].Dimension() != isotable[1].Dimension()) {
      cerr << "Table dimension mismatch:" << endl;
      for (int i = 0; i < 2; i++) 
	cerr << "  " << fname[i] << " has dimension " 
	     << isotable[i].Dimension() << "." << endl;
      exit(20);
    };

    diff_table_polyhedra(isotable);

    if (isotable[0].NumTableEntries() != 
	isotable[1].NumTableEntries()) {
      cerr << "Table size mismatch:" << endl;
      for (int i = 0; i < 2; i++) 
	cerr << "  " << fname[i] << " has " 
	     << isotable[i].NumTableEntries() << " entries." << endl;
      exit(20);
    };

    if (isotable[0].NumVerticesPerSimplex() != 
	isotable[1].NumVerticesPerSimplex()) {
      cerr << "Table simplices mismatch:" << endl;
      for (int i = 0; i < 2; i++) 
	cerr << "  " << fname[i] << " has " 
	     << isotable[i].NumVerticesPerSimplex() 
	     << " number of vertices per simplex." << endl;
      exit(25);
    }

    int numv_per_simplex = isotable[0].NumVerticesPerSimplex();
    for (int i = 0; i < 2; i++) {
      vlist[i] = new int[numv_per_simplex];
    };

    for (int it = 0; 
	 it < isotable[0].NumTableEntries(); it++) {
      if (isotable[0].NumSimplices(it) != isotable[1].NumSimplices(it)) {
	cerr << "Table entry " << it << " mismatch:" << endl;
	for (int i = 0; i < 2; i++) 
	  cerr << "  " << fname[i] << " has " 
	       << isotable[i].NumSimplices(it) << " simplices." << endl;
	continue;
      };

      for (int is0 = 0; is0 < isotable[0].NumSimplices(it); is0++) {
	bool found_matching_simplex = false;
	sort_simplex_vertices(isotable[0], it, is0, vlist[0]);
	for (int is1 = 0; is1 < isotable[1].NumSimplices(it); is1++) {
	  sort_simplex_vertices(isotable[1], it, is1, vlist[1]);
	  if (equals(vlist[0], vlist[1], numv_per_simplex)) {
	    found_matching_simplex = true;
	  };
	};
	if (!found_matching_simplex) {
	  cerr << "Table entry " << it << " mismatch:" << endl;
	  cerr << "  " << fname[0] << ", simplex " << is0
	       << " has no match in " << fname[1] << "." << endl;
	  continue;
	};
      };
    };

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


  for (int i = 0; i < 2; i++) {
    delete vlist[i];
    vlist[i] = NULL;
  };

  exit(0);
}

void diff_table_polyhedra(const ISOSURFACE_TABLE isotable[2])
// output differences in isosurface table polyhedra
// halt if vertex/edge differences are found
// do not halt if facet differences are found
{
  int dimension = isotable[0].Polytope().Dimension();

  assert(isotable[0].Polytope().Dimension() ==
	 isotable[1].Polytope().Dimension());

  if (isotable[0].Polytope().NumVertices() != 
      isotable[1].Polytope().NumVertices()) {
    cerr << "Polytope number of vertices mismatch:" << endl;
    for (int i = 0; i < 2; i++) 
      cerr << "  " << fname[i] << " has " 
	   << isotable[i].Polytope().NumVertices() 
	   << " vertices." << endl;
    exit(20);
  };

  if (isotable[0].Polytope().NumEdges() != 
      isotable[1].Polytope().NumEdges()) {
    cerr << "Polytope number of edges mismatch:" << endl;
    for (int i = 0; i < 2; i++) 
      cerr << "  " << fname[i] << " has " 
	   << isotable[i].Polytope().NumEdges() << " edges." << endl;
    exit(20);
  };

  for (int iv = 0; iv < isotable[0].Polytope().NumVertices(); iv++) {
    for (int ic = 0; ic < dimension; ic++) {
      int coord[2];
      coord[0] = isotable[0].Polytope().VertexCoord(iv, ic);
      coord[1] = isotable[1].Polytope().VertexCoord(iv, ic);
      if (coord[0] != coord[1]) {
	cerr << "Polytope vertex coordinate mismatch:" << endl;
	for (int i = 0; i < 2; i++) 
	  cerr << "  " << fname[i] << ", vertex " << iv
	       << ", coordinate " << ic << " = " 
	       << coord[i] << "." << endl;
	exit(20);
      };
    };
  };

  for (int ie = 0; ie < isotable[0].Polytope().NumEdges(); ie++) {
    int iend[2][2];

    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) {
	iend[i][j] = isotable[i].Polytope().EdgeEndpoint(ie, j);
      };

    if ((iend[0][0] != iend[1][0] && iend[0][0] != iend[1][1]) ||
	(iend[0][1] != iend[1][0] && iend[0][1] != iend[1][1])) {
      cerr << "Polytope edge mismatch:" << endl;
      for (int i = 0; i < 2; i++) 
	cerr << "  " << fname[i] << ", edge " << ie
	     << " = (" << iend[i][0] << ", " << iend[i][1] << ")." << endl;
      exit(20);
    };
  };

  // output facet differences
  // do not halt if facet differences are found
  if (isotable[0].Polytope().NumFacets() != 
      isotable[1].Polytope().NumFacets()) {
    cerr << "Polytope number of facets mismatch:" << endl;
    for (int i = 0; i < 2; i++) 
      cerr << "  " << fname[i] << " has " 
	   << isotable[i].Polytope().NumFacets() << " facets." << endl;
    return;
  };

  for (int jf = 0; jf < isotable[0].Polytope().NumFacets(); jf++) {
    if (isotable[0].Polytope().NumFacetVertices(jf) != 
	isotable[1].Polytope().NumFacetVertices(jf)) {
      cerr << "Polytope number of facet vertices mismatch:" << endl;
      for (int i = 0; i < 2; i++) 
	cerr << "  " << fname[i] << ", facet " << jf
	     << " has " << isotable[i].Polytope().NumFacetVertices(jf)
	     << " vertices." << endl;
      continue;
    };

    if (isotable[0].Polytope().Facet(jf) != 
	isotable[1].Polytope().Facet(jf)) {
      cerr << "Polytope facet vertices mismatch:" << endl;
      for (int iv = 0; iv < isotable[0].Polytope().NumVertices(); 
	   iv++) {
	if (isotable[0].Polytope().IsVertexInFacet(jf, iv) !=
	    isotable[1].Polytope().IsVertexInFacet(jf, iv)) {
	  for (int i = 0; i < 2; i++) 
	    if (isotable[i].Polytope().IsVertexInFacet(jf, iv)) 
	      cerr << "  " << fname[i] << ", facet " << jf
		   << " has vertex " << iv << "." << endl;
	    else
	      cerr << "  " << fname[i] << ", facet " << jf
		   << " does not have vertex " << iv << "." << endl;
	};
      };
    };
  };

}

void sort_simplex_vertices(const ISOSURFACE_TABLE & isotable, 
			   const int it, const int is, int * vlist)
  // store simplex vertices in vlist and sort
  // Precondition: vlist is preallocated to size at least
  //               isotable.NumVerticesPerSimplex
{
  int numv_per_simplex = isotable.NumVerticesPerSimplex();
  for (int iv = 0; iv < numv_per_simplex; iv++)
    vlist[iv] = isotable.SimplexVertex(it, is, iv);
  sort(vlist, vlist+numv_per_simplex);
}

bool equals(const int * vlist0, const int * vlist1, const int numv)
  // return true if two lists are equal
{
  for (int iv = 0; iv < numv; iv++) {
    if (vlist0[iv] != vlist1[iv]) {
      return(false);
    };
  };
  return(true);
}

void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc) {
    if (strcmp(argv[1], "-h") == 0) 
      { help_msg(); }
    else
      { break; }
  };

  if (argc != 3) { usage_error(); };

  fname[0] = argv[1];
  fname[1] = argv[2];
}

void usage_error()
{
  cerr << "Usage: ijkdifftable [-h] {file1} {file2}" << endl;
  exit(10);
}

void help_msg()
{
  cout << "ijkdifftable - output the difference of two isosurface tables" << endl;
  cout << "  Output differences in polyhedra and in table entries." << endl;
  cout << "  Ignore simplex order in table entries." << endl;
  cout << "  Ignore vertex order in comparing isosurface simplices." << endl;
  cout << endl;
  cout << "Usage: diffisotable [-h] {file1} {file2}" << endl;
  cout << "  -h : Print this help message (and exit.)" << endl;
  exit(0);
}
