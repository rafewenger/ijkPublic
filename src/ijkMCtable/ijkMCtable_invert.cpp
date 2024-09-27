/*!
 *  @file ijkMCtable_invert.cpp
 *  @brief Swap table entry it0 and it1 where it1 is index it0
 *    with all '+' mapped to '-' and all '-' mapped to '+'.
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2023-2024 Rephael Wenger

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

using namespace IJK;
using namespace IJKMCUBE_TABLE;
using namespace std;

typedef float COORD_TYPE;

// global variables
const char * input_isotable_filename(NULL);
const char * output_isotable_filename(NULL);
bool flag_replace_file(false);

// routines
void read_isotable(const char * input_filename,
                   ISOSURFACE_TABLE & isotable);
void write_isotable(const char * output_filename,
                    const ISOSURFACE_TABLE & isotable);
void parse_command_line(int argc, char **argv);
void help_msg(), usage_error();

int main(int argc, char **argv)
{
  ISOSURFACE_TABLE isotableA, isotableB;
  
  try {
    parse_command_line(argc, argv);

    read_isotable(input_isotable_filename, isotableA);
    invert_mcube_isotable(isotableA, isotableB);
    write_isotable(output_isotable_filename, isotableB);
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


//**************************************************
// READ ISOSURFACE TABLE
//**************************************************

// Read isosurface lookup table.
void read_isotable
(const char * input_filename,
 ISOSURFACE_TABLE & isotable)
{
  IJK::PROCEDURE_ERROR error("read_isotable");
  
  ifstream isotable_file(input_filename, ios::in);

  if (!isotable_file) {
    error.AddMessage
      ("Unable to open isosurface table file ",
       string(input_filename), ".");
    throw error;
  };

  try {
    IJKXIO::read_xit(isotable_file, isotable);
  }
  catch(...) {
    cerr << "Error reading file " << input_filename << "." << endl;
    throw;
  };

  isotable_file.close();

  if (!isotable.Check(error)) {
    cerr << "Warning: Data structure inconsistency in isosurface table."
	 << endl;
    error.Print(cerr);
    cerr << "  Attempting to continue..." << endl << endl;
  };
}


//**************************************************
// WRITE ISOSURFACE TABLE
//**************************************************

// Write isosurface lookup table.
void write_isotable
(const char * output_filename,
 const ISOSURFACE_TABLE & isotable)
{
  IJK::PROCEDURE_ERROR error("write_isotable");
  
  ofstream isotable_file(output_filename, ios::out);

  if (!isotable_file) {
    error.AddMessage
      ("Unable to open isosurface table file ",
       output_filename, " for writing.");
    error.AddMessage("Exiting...");
    throw error;
  };

  cout << "Writing isosurface table to file " 
       << output_filename << "." << endl;

  try {
    // *** NEED TO REVISE: SHOULD USE SAME VERSION AS INPUTE FILE ***
    IJKXIO::write_xit
      (isotable_file, IJKXIO::XIT_VERSION_2_0, isotable);
  }
  catch(...) {
    cerr << "Error writing file " << output_filename
         << "." << endl;
    isotable_file.close();
    throw;
  }
  
  isotable_file.close();
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

    const std::string s = argv[iarg];
    if (s == "-replace")
      { flag_replace_file = true; }
    else {
      cerr << "Illegal option: " << s << endl;
      usage_error();
    }

    iarg++;
  };

  if (iarg >= argc)
    { usage_error(); }

  input_isotable_filename = argv[iarg];
  iarg++;

  if (iarg+1 != argc)
    { usage_error(); }

  output_isotable_filename = argv[iarg];

  if (flag_replace_file) {
    if (output_isotable_filename == NULL) {
      output_isotable_filename = input_isotable_filename;
    }
    else {
      cerr << "Warning: Flag -replace but output isotable filename provided."
           << endl;
      cerr << "  Ignoring flag -replace." << endl;
    }
  }
  else if (iarg >= argc) {
    cerr << "Usage error. No output filename provided." << endl;
    cerr << "  Include output filename or use flag -replace." << endl;
    exit(10);
  }
}


void usage_msg()
{
  cout << "Usage: ijkmcube_table_invert [-h] [-replace] {input isotable file} [{output isotable file}]"
       << endl;
}

 
void usage_error()
{
  usage_msg();
  exit(10);
}

 
void help_msg()
{
  cout << "ijkmcube_table_invert - Invert isosurface lookup table."
       << endl;
  usage_msg();
  cout << "-replace: Replace input file with new isotable." << endl;
  cout << "-h:       Print this help message (and exit)." << endl;

  exit(0);
}
