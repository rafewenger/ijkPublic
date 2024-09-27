/// Test ijktable_ambig.

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2012 Rephael Wenger

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

#include <cstdlib>
#include <iostream>

#include "ijkcube.tpp"
#include "ijkstring.tpp"


#include "ijkMCtable_ambig.h"

using namespace std;
using namespace IJK;

// global variables
int dimension(3);
bool use_cube_lex_order(true);

// output routines
void output_ambig_table(const int dimension);

void memory_exhaustion();
void usage_error();
void parse_command_line(int argc, char **argv);


int main(int argc, char **argv)
{
  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    output_ambig_table(dimension);
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

void output_ambig_table(const int dimension)
{
  IJKMCUBE_TABLE::ISOSURFACE_TABLE_AMBIG_INFO table;
  IJK::CUBE_INFO<int,int> cube(dimension);

  table.SetCubeAmbiguityTable(dimension, use_cube_lex_order);

  cout << "Dimension: " << dimension << endl;
  cout << "Number of table entries: "
       << table.NumTableEntries() << endl;
  for (int it = 0; it < table.NumTableEntries(); it++) {
    cout << "Entry " << it << ":";
    if (table.IsAmbiguous(it)) 
      { cout << "  Ambiguous."; }
    else
      { cout << "  Not ambiguous."; }
    cout << "  Num ambiguous facets: " << int(table.NumAmbiguousFacets(it))
         << ".";
    if (table.NumAmbiguousFacets(it) > 0) {
      cout << "  Ambiguous facets:";
      for (int k = 0; k < cube.NumFacets(); k++) {
        if (table.IsFacetAmbiguous(it, k)) {
          cout << "  " << k;
        }
      }
    }
    cout << endl;
  }
  cout << endl;
}


// **************************************************
// Miscellaneous routines
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

void usage_msg()
{
  cerr << "Usage: outambigable [-orderA | -lexorder] <dimension>" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void parse_command_line(int argc, char **argv)
{

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = argv[iarg];
    if (s == "-orderA") {
      use_cube_lex_order = false;
    }
    else if (s == "-lexorder") {
      use_cube_lex_order = true;
    }
    else {
      usage_error();
    }

    iarg++;
  }

  if (argc != iarg+1) { usage_error(); }

  if (!IJK::string2val(argv[iarg], dimension)) {
    cerr << "Error in argument 1: " << argv[iarg] << endl;
    cerr << "Non-integer character in string: " << argv[iarg] << endl;
    exit(50);
  }
  iarg++;

  
  if (iarg != argc) { usage_error(); }
}
