/*!
 *  @file ijkMCtable.cpp
 *  @brief Class containing a table of Marching Cubes isosurface patches
 *    in a given polytope.
 *  - Version 0.5.0
 */


/*
  IJK: Isosurface Jeneration Kode
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

#include <ctype.h>
#include <limits>
#include <limits.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <set>
#include <vector>

#include "ijkMCtable.h"

#include "ijkbits.tpp"
#include "ijksimplex.tpp"

using namespace IJK;
using namespace IJKMCUBE_TABLE;
using namespace std;

#ifndef LONG_BIT

#define LONG_BIT (CHAR_BIT * sizeof(long))

#endif


// ***************************************************************
// TABLE MANIPULATION ROUTINES
// ***************************************************************

void IJKMCUBE_TABLE::invert_mcube_isotable
(const ISOSURFACE_TABLE & isotableA, ISOSURFACE_TABLE & isotableB)
{
  const TABLE_INDEX num_table_entries =
    isotableA.NumTableEntries();
  isotableB.CopyPolytope(isotableA);
  isotableB.CopyIsosurfaceVertices(isotableA);
  isotableB.SetNumTableEntries(num_table_entries);
  isotableB.SetSimplexDimension(isotableA.SimplexDimension());
  isotableB.CopyProperties(isotableA);

  // Inverted table has opposite separation type and orientation.
  isotableB.SetSeparationType(isotableA.OppositeSeparationType());
  isotableB.SetIsoPolyOrientation(isotableA.OppositeIsoPolyOrientation());

  const int num_vertices_per_simplex =
    isotableA.NumVerticesPerSimplex();
  
  for (TABLE_INDEX itB = 0; itB < num_table_entries; itB++) {
    const TABLE_INDEX itA = num_table_entries - itB - 1;
    const int num_simplices = isotableA.NumSimplices(itA);

    isotableB.SetSimplexVertices
      (itB, isotableA.SimplexVertices(itA),
       isotableA.NumSimplices(itA));
  }
}



// *******************************************************************
// ISOSURFACE_VERTEX member functions
// *******************************************************************

// Constructor.
ISOSURFACE_VERTEX::ISOSURFACE_VERTEX()
{ 
  coord = NULL; 
  num_coord = 0; 
  is_label_set = false; 
}

ISOSURFACE_VERTEX::~ISOSURFACE_VERTEX()
  // destructor
{
  if (coord != NULL) { 
    delete [] coord;
    coord = NULL;
  };

  num_coord = 0;
  is_label_set = false;
}

void ISOSURFACE_VERTEX::SetNumCoord(const int numc)
{
  if (coord != NULL)
    delete [] coord;
  coord = NULL;

  num_coord = numc;
  coord = new COORD_TYPE[numc];
}


// *******************************************************************
// ISOSURFACE_TABLE_PROPERTIES member functions
// *******************************************************************

void ISOSURFACE_TABLE_PROPERTIES::Init()
{
  lookup_table_type = table_type_list.UndefinedValue();
  encoding = encoding_list.UndefinedValue();
  grid_vertex_label_type = grid_vertex_label_type_list.UndefinedValue();
  isosurface_separation_type = 
    isosurface_separation_type_list.UndefinedValue();
  isosurface_triangulation_type =
    isosurface_triangulation_type_list.UndefinedValue();
  iso_poly_orientation =
    iso_poly_orientation_list.UndefinedValue();
  separate_opposite =
    separate_opposite_type_list.UndefinedValue();
}


std::string ISOSURFACE_TABLE_PROPERTIES::SeparationTypeLabel() const
{
  if (isosurface_separation_type == SEPARATE_NEG) 
    { return std::string("sepNeg"); }
  else if (isosurface_separation_type == SEPARATE_POS) 
    { return std::string("sepPos"); }
  else 
    { return std::string(""); }
}


std::string ISOSURFACE_TABLE_PROPERTIES::IsoPolyOrientationLabel() const
{
  if (iso_poly_orientation == POSITIVE_ORIENT) 
    { return std::string("posO"); }
  else if (iso_poly_orientation == NEGATIVE_ORIENT) 
    { return std::string("negO"); }
  else 
    { return std::string(""); }
}

std::string ISOSURFACE_TABLE_PROPERTIES::TriangulationTypeLabel() const
{
  if (isosurface_triangulation_type == CONVEX_HULL)
    { return std::string("cHull"); }
  else if (isosurface_triangulation_type == EDGE_GROUPS)
    { return std::string("edgeGroups"); }
  else 
    { return std::string(""); }
}


// Return opposite separation type.
ISOSURFACE_SEPARATION_TYPE ISOSURFACE_TABLE_PROPERTIES::
OppositeSeparationType() const
{
  if (isosurface_separation_type == SEPARATE_NEG)
    { return SEPARATE_POS; }
  else if (isosurface_separation_type == SEPARATE_POS)
    { return SEPARATE_NEG; }
  else {
    // All other values are unchanged.
    return isosurface_separation_type;
  }
}


// Return opposite orientation type.
ISO_POLY_ORIENTATION ISOSURFACE_TABLE_PROPERTIES::
OppositeIsoPolyOrientation() const
{
  if (iso_poly_orientation == NEGATIVE_ORIENT)
    { return POSITIVE_ORIENT; }
  else if (iso_poly_orientation == POSITIVE_ORIENT)
    { return NEGATIVE_ORIENT; }
  {
    // All other values are unchanged.
    return iso_poly_orientation;
  }
}


// Copy properties from isotable_properties.
void ISOSURFACE_TABLE_PROPERTIES::Copy
(const ISOSURFACE_TABLE_PROPERTIES & isotable_properties)
{
  lookup_table_type = isotable_properties.TableType();
  encoding = isotable_properties.Encoding();
  grid_vertex_label_type = isotable_properties.GridVertexLabelType();
  isosurface_triangulation_type = 
    isotable_properties.TriangulationType();
  isosurface_separation_type = isotable_properties.SeparationType();
  separate_opposite = isotable_properties.SeparateOpposite();
  iso_poly_orientation = isotable_properties.IsoPolyOrientation();
}


// Return true if isotable has properties specified in properties.
bool ISOSURFACE_TABLE_PROPERTIES::Check
(const ISOSURFACE_TABLE_PROPERTIES & properties,
 IJK::ERROR & error) const
{
  const bool flag0 =
    CheckTableType(properties.TableType(), error);
  const bool flag1 =
    CheckEncoding(properties.Encoding(), error);
  const bool flag2 = 
    CheckTriangulationType(properties.TriangulationType(), error);
  const bool flag3 =
    CheckSeparationType(properties.SeparationType(), error);
  const bool flag4 =
    CheckOrientation(properties.IsoPolyOrientation(), error);
  const bool flag5 =
    CheckSeparateOpposite(properties.SeparateOpposite(), error);

  return (flag0 && flag1 && flag2 && flag3 && flag4 && flag5);
}


// Return true if lookup_table_type matches table_type.
bool ISOSURFACE_TABLE_PROPERTIES::CheckTableType
(const LOOKUP_TABLE_TYPE table_type, IJK::ERROR & error) const
{
  if (table_type_list.IsUndefined(table_type)) {
    // Nothing to check.
    return true;
  }

  if (table_type != TableType()) {
    error.AddMessage("Incorrect isosurface lookup table type.");
    error.AddMessage
      ("  Isotable type: ", TableTypeString());
    error.AddMessage
      ("  Expected table type: ", TableTypeString(table_type));

    return false;
  }

  return true;
}


// Return true if this->encoding matches encoding.
bool ISOSURFACE_TABLE_PROPERTIES::CheckEncoding
(const ENCODING encoding, IJK::ERROR & error) const
{
  if (encoding_list.IsUndefined(encoding)) {
    // Nothing to check.
    return true;
  }

  if (encoding != Encoding()) {
    error.AddMessage("Incorrect isosurface lookup table encoding.");
    error.AddMessage
      ("  Isotable encoding: ", EncodingString());
    error.AddMessage("  Expected encoding: ", EncodingString(encoding));

    return false;
  }

  return true;
}


// Return true if isosurface_triangulation_type matches tri_type.
bool ISOSURFACE_TABLE_PROPERTIES::CheckTriangulationType
(const ISOSURFACE_TRIANGULATION_TYPE tri_type,
 IJK::ERROR & error) const
{
  if (isosurface_triangulation_type_list.IsUndefined(tri_type)) {
    // Nothing to check.
    return true;
  }

  if (tri_type == UNKNOWN_ISOSURFACE_TRIANGULATION_TYPE) {
    // Nothing to check.
    return true;
  }

  if (tri_type != TriangulationType()) {
    error.AddMessage("Incorrect isosurface lookup table triangulation type.");
    error.AddMessage
      ("  Isotable triangulation type: ", TriangulationTypeString());
    error.AddMessage
      ("  Expected triangulation type: ", 
       TriangulationTypeString(tri_type));

    return false;
  }

  return true;
}


// Return true if iso poly orientation matches _orientation.
bool ISOSURFACE_TABLE_PROPERTIES::CheckOrientation
(const ISO_POLY_ORIENTATION _orientation,
 IJK::ERROR & error) const
{
  if (iso_poly_orientation_list.IsUndefined(_orientation)) {
    // Nothing to check.
    return true;
  }

  if (_orientation == NO_ORIENT) {
    // Nothing to check.
    return true;
  }

  if (_orientation != IsoPolyOrientation()) {
    error.AddMessage("Incorrect isosurface lookup table polytope orientation.");
    error.AddMessage
      ("  Isotable polytope orientation: ", IsoPolyOrientationString());
    error.AddMessage
      ("  Expected polytope orientation: ", 
       IsoPolyOrientationString(_orientation));

    return false;
  }

  return true;
}


// Return true if isosurface_separation_type matches _sep_type.
bool ISOSURFACE_TABLE_PROPERTIES::CheckSeparationType
(const ISOSURFACE_SEPARATION_TYPE _sep_type, IJK::ERROR & error) const
{
  if (isosurface_separation_type_list.IsUndefined(_sep_type)) {
    // Nothing to check.
    return true;
  }

  if (_sep_type == UNKNOWN_SEPARATION_TYPE) {
    // Nothing to check.
    return true;
  }

  if (_sep_type != SeparationType()) {
    error.AddMessage("Incorrect isosurface lookup table separation type.");
    error.AddMessage
      ("  Isotable separation type: ", SeparationTypeString());
    error.AddMessage
      ("  Expected separation type: ", SeparationTypeString(_sep_type));

    return false;
  }

  return true;
}


// Return true if this->separate_opposite matches _sep_opposite.
bool ISOSURFACE_TABLE_PROPERTIES::CheckSeparateOpposite
(const SEPARATE_OPPOSITE_TYPE _sep_opposite, IJK::ERROR & error) const
{
  if (separate_opposite_type_list.IsUndefined(_sep_opposite)) {
    // Nothing to check.
    return true;
  }

  if (_sep_opposite != SeparateOpposite()) {
    error.AddMessage("Incorrect isosurface lookup table separate opposite flag.");
    error.AddMessage
      ("  Isotable separate opposite: ", SeparateOppositeString());
    error.AddMessage
      ("  Expected separate opposit: ", SeparateOppositeString(_sep_opposite));

    return false;
  }

  return true;
}


// *******************************************************************
// ISOSURFACE_TABLE
// *******************************************************************

// Constructor.
ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::ISOSURFACE_TABLE_ENTRY()
{
  num_simplices = 0;
  simplex_vertex_list = NULL;
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::~ISOSURFACE_TABLE_ENTRY()
  // destructor
{
  FreeAll();
}

bool ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::Check
(ERROR & error_msg) const
{
  if (num_simplices < 0) {
    error_msg.ClearAll();
    error_msg.AddMessage
      ("Number of simplices in isosurface table entry must be non-negative.");
    return(false);
  }

  if (num_simplices > 0 && simplex_vertex_list == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for simplex vertex list not allocated.");
    return(false);
  }

  return(true);
}

void ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::FreeAll()
  // free all memory
{
  delete [] simplex_vertex_list;
  simplex_vertex_list = NULL;
  num_simplices = 0;
}


// Default constructor. dimension = 3.
ISOSURFACE_TABLE::ISOSURFACE_TABLE(): polytope(3)
{
  Init(3, 2);
}


// Constructor.
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
ISOSURFACE_TABLE::ISOSURFACE_TABLE(const int d): polytope(d)
{
  Init(d, d-1);
}


// Constructor.
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
ISOSURFACE_TABLE::ISOSURFACE_TABLE
(const int dimension, const int simplex_dimension): polytope(dimension)
{
  Init(dimension, simplex_dimension);
}


// Constructor.
// @param d Dimension of space containing isosurface.  Should be 2, 3 or 4.
void ISOSURFACE_TABLE::Init(const int dimension, const int simplex_dimension)
{
  const char * procname = "ISOSURFACE_TABLE::Init";

  this->simplex_dimension = simplex_dimension;

  max_num_vertices = 20;
  // Note: Even tables for polytope of this size are probably impossible 
  //   to compute/store

  num_isosurface_vertices = 0;
  isosurface_vertex = NULL;

  num_table_entries = 0;
  entry = NULL;
  is_table_allocated = false;
  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal polytope dimension.");
}


// Destructor.
ISOSURFACE_TABLE::~ISOSURFACE_TABLE()
{
  FreeAll();
}


void ISOSURFACE_TABLE::SetEncoding(const ENCODING encoding)
{
  this->table_properties.encoding = encoding;
}


void ISOSURFACE_TABLE::SetTableType
(const LOOKUP_TABLE_TYPE lookup_table_type)
{
  this->table_properties.lookup_table_type = lookup_table_type;
}


void ISOSURFACE_TABLE::SetGridVertexLabelType
(const GRID_VERTEX_LABEL_TYPE grid_vertex_label_type)
{
  this->table_properties.grid_vertex_label_type = grid_vertex_label_type;
}


void ISOSURFACE_TABLE::SetSeparationType
(const ISOSURFACE_SEPARATION_TYPE separation_type)
{
  this->table_properties.isosurface_separation_type = separation_type;
}

void ISOSURFACE_TABLE::SetTriangulationType
(const ISOSURFACE_TRIANGULATION_TYPE triangulation_type)
{
  this->table_properties.isosurface_triangulation_type = 
    triangulation_type;
}

void ISOSURFACE_TABLE::SetIsoPolyOrientation
(const ISO_POLY_ORIENTATION iso_poly_orientation)
{
  this->table_properties.iso_poly_orientation = iso_poly_orientation;
}

void ISOSURFACE_TABLE::SetSeparateOpposite
(const SEPARATE_OPPOSITE_TYPE separate_opposite)
{
  this->table_properties.separate_opposite = separate_opposite;
}

void ISOSURFACE_TABLE::SetSeparateOpposite(const bool flag)
{
  this->table_properties.SetSeparateOpposite(flag);
}

void ISOSURFACE_TABLE::SetNumIsosurfaceVertices(const int num_vertices)
{
  if (isosurface_vertex != NULL)
    delete [] isosurface_vertex;
  isosurface_vertex = NULL;

  num_isosurface_vertices = num_vertices;
  isosurface_vertex = new ISOSURFACE_VERTEX[num_vertices];
}


// Check allocation of array isosurface_vertices
// procname = calling procedure name, for error messages
// vstart = first vertex
// numv = number of vertices required
void ISOSURFACE_TABLE::CheckIsoVerticesAlloc
(const char * procname, const int vstart, const int numv)
{
  if (numv == 0) return;

  if (isosurface_vertex == NULL) {
    throw PROCEDURE_ERROR
      (procname, "Set number of isosurface vertices before storing vertices.");
  }

  if (numv+vstart > NumIsosurfaceVertices()) {
    throw PROCEDURE_ERROR
      (procname, "Illegal isosurface vertex index.");
  }
}

void ISOSURFACE_TABLE::StorePolyVerticesAsIsoVertices(const int vstart)
  // store polytope vertices as isosurface vertices
  // store polytope vertices starting at isosurface vertex vstart
{
  const int num_polyv = Polytope().NumVertices();
  const char * procname = "ISOSURFACE_TABLE::StorePolyVerticesAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polyv);

  for (int iv = 0; iv < num_polyv; iv++) {
    SetIsoVertexType(iv+vstart, ISOSURFACE_VERTEX::VERTEX);
    SetIsoVertexFace(iv+vstart, iv);
  }
}


// @brief Set isosurface vertex.
void ISOSURFACE_TABLE::SetIsosurfaceVertex
(const int iv, const ISOSURFACE_VERTEX & isosurface_vertex)
{
  IJK::PROCEDURE_ERROR error
    ("ISOSURFACE_TABLE::SetIsosurfaceVertex");
  
  if (NumIsosurfaceVertices() == 0) {
    error.AddMessage
      ("Programming error. Call ISOSURFACE_TABLE::SetNumIsosurfaceVertices()");
    error.AddMessage
      ("  before calling ISOSURFACE_TABLE::SetIsosurfaceVertex().");
    throw error;
  }

  if (iv >= NumIsosurfaceVertices()) {
    error.AddMessage
      ("Programming error. Illegal isosurface vertex index ", iv, ".");
    error.AddMessage
      ("  Isosurface vertices should be in range [0..",
       NumIsosurfaceVertices()-1, "].");
    error.AddMessage
      ("  Check call to ISOSURFACE_TABLE:SetNumIsosurfaceVertices().");
    throw error;
  }
  
  SetIsoVertexType(iv, isosurface_vertex.Type());
  SetIsoVertexFace(iv, isosurface_vertex.Face());
  SetIsoVertexNumCoord(iv, isosurface_vertex.NumCoord());
  for (int ic = 0; ic < isosurface_vertex.NumCoord(); ic++)
    { SetIsoVertexCoord(iv, ic, isosurface_vertex.Coord(ic)); }
  if (isosurface_vertex.IsLabelSet())
    { SetIsoVertexLabel(iv, isosurface_vertex.Label()); }
}


// Copy isosurface vertices from isotable.
void ISOSURFACE_TABLE::CopyIsosurfaceVertices
(const ISOSURFACE_TABLE & isotable)
{
  const int num_isosurface_vertices =
    isotable.NumIsosurfaceVertices();
  
  SetNumIsosurfaceVertices(num_isosurface_vertices);
  for (int iv = 0; iv < num_isosurface_vertices; iv++) {
    SetIsosurfaceVertex(iv, isotable.IsosurfaceVertex(iv));
  }
}


void ISOSURFACE_TABLE::StorePolyEdgesAsIsoVertices(const int vstart)
  // store polytope edges as isosurface vertices
  // store polytope edges starting at isosurface vertex vstart
{
  const int num_polye = Polytope().NumEdges();
  const char * procname = "ISOSURFACE_TABLE::StorePolyEdgesAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polye);

  for (int ie = 0; ie < num_polye; ie++) {
    SetIsoVertexType(ie+vstart, ISOSURFACE_VERTEX::EDGE);
    SetIsoVertexFace(ie+vstart, ie);
  }
}


// Store polytope facets as isosurface vertices.
// store polytope facets starting at isosurface vertex vstart.
void ISOSURFACE_TABLE::StorePolyFacetsAsIsoVertices(const int vstart)
{
  const int num_polyf = Polytope().NumFacets();
  const char * procname = "ISOSURFACE_TABLE::StorePolyFacetsAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polyf);

  for (int jf = 0; jf < num_polyf; jf++) {
    SetIsoVertexType(jf+vstart, ISOSURFACE_VERTEX::FACET);
    SetIsoVertexFace(jf+vstart, jf);
  }
}


// Allocate table
void ISOSURFACE_TABLE::SetNumTableEntries(const int num_table_entries)
{
  const char * procname = "ISOSURFACE_TABLE::SetNumTableEntries";

  if (entry != NULL) delete [] entry;
  entry = NULL;
  this->num_table_entries = 0;
  is_table_allocated = false;

  entry = new ISOSURFACE_TABLE_ENTRY[num_table_entries];
  if (entry == NULL && num_table_entries > 0)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for isosurface table.");

  this->num_table_entries = num_table_entries;
  is_table_allocated = true;
}


void ISOSURFACE_TABLE::SetNumSimplices(const TABLE_INDEX it, const int nums)
  // set number of simplices in table entry it
  // it = table entry
  // nums = number of simplices
{
  const char * procname = "ISOSURFACE_TABLE::SetNumSimplices";

  if (!IsTableAllocated() || entry == NULL) {
    throw PROCEDURE_ERROR
      (procname, "Table must be allocated before entering table entries.");
  };

  if (it < 0 || it >= NumTableEntries())
    throw PROCEDURE_ERROR(procname, "Illegal table index.");
  if (nums < 0)
    throw PROCEDURE_ERROR
      (procname, "Number of simplices must be non-negative.");

  entry[it].num_simplices = 0;
  delete entry[it].simplex_vertex_list;
  entry[it].simplex_vertex_list = NULL;

  if (nums > 0)
    entry[it].simplex_vertex_list = 
      new ISOSURFACE_VERTEX_INDEX[nums*NumVerticesPerSimplex()];

  entry[it].num_simplices = nums;
}


// Set simplex vertex.
// it = index table entry.  In range [0..NumTableEntries()-1].
// is = index simplex.  
// k = k'th simplex vertex.  In range [0..NumVerticesPerSimplex()-1].
// isov = index of isosurface vertex
void ISOSURFACE_TABLE::SetSimplexVertex
(const TABLE_INDEX it, const int is, const int k, 
 const ISOSURFACE_VERTEX_INDEX isov)
{
  entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+k] = isov;
}


// Set simplex vertices.
void ISOSURFACE_TABLE::SetSimplexVertices
(const TABLE_INDEX it, const ISOSURFACE_VERTEX_INDEX simplex_vertices[],
 const int num_simplices)
{
  SetNumSimplices(it, num_simplices);

  std::copy(simplex_vertices,
            simplex_vertices+num_simplices*NumVerticesPerSimplex(),
            entry[it].simplex_vertex_list);
}


// Return true if facet vertex labels are identical
//   in table entries table_indexA and table_indexB.
bool ISOSURFACE_TABLE::AreAllFacetVertexLabelsIdentical
(const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
 const int ifacet) const
{
  const int num_poly_vertices = Polytope().NumVertices();
  std::vector<int> digitA(num_poly_vertices);
  std::vector<int> digitB(num_poly_vertices);
  IJK::PROCEDURE_ERROR error
    ("ISOSURFACE_TABLE::AreAllFacetVertexLabelsIdentical");

  IJK::convert2base
    (table_indexA, Base(), digitA.data(), num_poly_vertices, error);
  IJK::convert2base
    (table_indexB, Base(), digitB.data(), num_poly_vertices, error);
  
  for (int j = 0; j < Polytope().NumFacetVertices(ifacet); j++) {
    const int jv = Polytope().FacetVertex(ifacet, j);

    if (digitA[jv] != digitB[jv])
      { return false; }
  }

  return true;
}


// Sort simplex vertices in increasing order.
void ISOSURFACE_TABLE::SortSimplexVertices
(const TABLE_INDEX it, const int isimplex)
{
  const int numv_per_simplex = NumVerticesPerSimplex();
  ISOSURFACE_VERTEX_INDEX * simplex_vert =
    entry[it].simplex_vertex_list + isimplex*numv_per_simplex;
  
  std::sort(simplex_vert, simplex_vert+numv_per_simplex);
}


// Flip all isosurface polytope orientations from +1 to -1 or from -1 to +1.
void ISOSURFACE_TABLE::FlipIsoPolyOrientation
(const TABLE_INDEX it, const int ipoly)
{
  const int numv_per_simplex = NumVerticesPerSimplex();
  ISOSURFACE_VERTEX_INDEX * simplex_vert =
    entry[it].simplex_vertex_list + ipoly*numv_per_simplex;

  if (numv_per_simplex < 2) {
    // Nothing to switch.
    return;
  }
  
  const int ilast = numv_per_simplex-1;
  std::swap(simplex_vert[ilast], simplex_vert[ilast-1]);
}


// Flip all isosurface polytope orientations at table entry table_index.
void ISOSURFACE_TABLE::FlipAllIsoPolyOrientations
(const TABLE_INDEX table_index)
{
  for (int isimplex = 0; isimplex < NumSimplices(table_index);
       isimplex++)
    { FlipIsoPolyOrientation(table_index, isimplex); }
}


// Flip all isosurface polytope orientations at in isosurface lookup table.
void ISOSURFACE_TABLE::FlipAllIsoPolyOrientations()
{
  for (TABLE_INDEX table_index = 0;
       table_index < NumTableEntries();
       table_index++)
    { FlipAllIsoPolyOrientations(table_index); }

  // Flip property iso_poly_orientation.
  SetIsoPolyOrientation(OppositeIsoPolyOrientation());
}


// Orient simplices in table entry.
void ISOSURFACE_TABLE::OrientSimplicesInTableEntry
(const TABLE_INDEX table_index, const TABLE_INDEX istart)
{
  IJK::orient_simplices
    (entry[table_index].simplex_vertex_list,
     NumVerticesPerSimplex(), NumSimplices(table_index), istart);
}


// Orient simplices in table entry.
// - Version that uses and sets array is_orientedA[].
void ISOSURFACE_TABLE::OrientSimplicesInTableEntry
(const TABLE_INDEX table_index, const TABLE_INDEX istart,
 std::vector<bool> & is_oriented)
{
  IJK::orient_simplices
    (entry[table_index].simplex_vertex_list,
     NumVerticesPerSimplex(), NumSimplices(table_index), istart,
     is_oriented);
}


// Orient all simplices in table entry.
void ISOSURFACE_TABLE::OrientAllSimplicesInTableEntry
(const TABLE_INDEX table_index, int & num_components)
{
  const int num_vert_per_simplex = SimplexDimension()+1;
  
  IJK::orient_all_simplices
    (entry[table_index].simplex_vertex_list, num_vert_per_simplex,
     NumSimplices(table_index), num_components);
}


// Return false if simplices in table entry have
//   inconsistent orientations.
bool ISOSURFACE_TABLE::AreSimplicesConsistentlyOriented
(const TABLE_INDEX table_index,
 int & isimplexA, int & isimplexB) const
{
  const int num_vert_per_simplex = SimplexDimension()+1;
  
  return IJK::are_simplices_consistently_oriented
    (entry[table_index].simplex_vertex_list, num_vert_per_simplex,
     NumSimplices(table_index), isimplexA, isimplexB);
}


bool ISOSURFACE_TABLE::CheckDimension(const int d) const
  // check dimension
{
  if (d < 1)
    return(false);
  else
    return(true);
}


bool ISOSURFACE_TABLE::CheckTable(ERROR & error_msg) const
  // check table
{
  if (polytope.NumVertices() >= LONG_BIT) {
    error_msg.ClearAll();
    error_msg.AddMessage("Too many polytope vertices");
    return(false);
  }

  if (polytope.NumVertices() > MaxNumVertices()) {
    error_msg.ClearAll();
    error_msg.AddMessage("Too many polytope vertices");
    return(false);
  }

  if (polytope.NumVertices() < 1) {
    error_msg.ClearAll();
    error_msg.AddMessage("Polytope must have at least one vertex.");
    return(false);
  }

  if (entry == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for isosurface table not allocated.");
    return(false);
  }

  for (int it = 0; it < NumTableEntries(); it++)
    if (!entry[it].Check(error_msg)) {
      error_msg.AddMessage("Error detected at isosurface table entry ", 
                           it, ".");
      return(false);
    }

  for (int jt = 0; jt < NumTableEntries(); jt++)
    for (int is = 0; is < NumSimplices(jt); is++)
      for (int iv = 0; iv < NumVerticesPerSimplex(); iv++) {
        int iso_v = int(SimplexVertex(jt, is, iv));
        if (iso_v < 0 || iso_v >= NumIsosurfaceVertices()) {
          error_msg.ClearAll();
          error_msg.AddMessage
            ("Illegal isosurface vertex ", iso_v, " in isosurface table entry ", jt, "."); 
          return(false);
        }
      };

  return(true);
}

bool ISOSURFACE_TABLE::Check(ERROR & error_msg) const
{
  if (!Polytope().Check(error_msg)) return(false);
  if (!CheckTable(error_msg)) return(false);
  return(true);
}


// Return false and set error message if table_index
//   is not in range [0..NumTableEntries()-1].
bool ISOSURFACE_TABLE::CheckTableIndex
(const int table_index, IJK::ERROR & error) const
{
  if (table_index < 0) {
    error.AddMessage
      ("Programming error. Illegal negative table index ",
       table_index, ".");
    error.AddMessage
      ("  Table index should be non-negative.");
    return false;
  }

  if (NumTableEntries() == 0) {
    error.AddMessage("Programming error. No table entries.");
    error.AddMessage
      ("  Call SetNumTableEntries() to create table entries.");
    return false;
  }

  if (table_index >= NumTableEntries()) {
    error.AddMessage("Programming error. Table index ",
                     table_index, " out of bounds.");
    error.AddMessage
      ("  Number of table entries: ", num_table_entries);
    error.AddMessage
      ("  Table index must be less than number of table entries.");
    return false;
  }

  return true;
}


// Free all memory
void ISOSURFACE_TABLE::FreeAll()
{
  if (entry != NULL) {
    for (int i = 0; i < num_table_entries; i++)
      entry[i].FreeAll();
    delete [] entry;
    entry = NULL;
  };
  num_table_entries = 0;
  is_table_allocated = false;

  polytope.FreeAll();

  delete [] isosurface_vertex;
  isosurface_vertex = NULL;
  num_isosurface_vertices = 0;
}


// *******************************************************************
// ISOSURFACE_EDGE_TABLE member functions
// *******************************************************************

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::
ISOSURFACE_EDGE_TABLE_ENTRY()
  // constructor
{
  num_edges = 0;
  edge_endpoint_list = NULL;
}


ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::
~ISOSURFACE_EDGE_TABLE_ENTRY()
  // destructor
{
  FreeAll();
}

bool ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::Check
(ERROR & error_msg) const
{
  if (num_edges < 0) {
    error_msg.ClearAll();
    error_msg.AddMessage
      ("Number of edges in isosurface table entry must be non-negative.");
    return(false);
  }

  if (num_edges > 0 && edge_endpoint_list == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for edge endpoint list not allocated.");
    return(false);
  }

  return(true);
}

void ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::FreeAll()
  // free all memory
{
  delete [] edge_endpoint_list;
  edge_endpoint_list = NULL;
  num_edges = 0;
}

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE
(const int d):ISOSURFACE_TABLE(d)
  // constructor
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(d);
}

void ISOSURFACE_EDGE_TABLE::Init(const int d)
  // constructor
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  edge_entry = NULL;
}

ISOSURFACE_EDGE_TABLE::~ISOSURFACE_EDGE_TABLE()
  // destructor
{
  FreeAll();
}

void ISOSURFACE_EDGE_TABLE::SetNumTableEntries(const int num_table_entries)
  // allocate table
{
  const char * procname = "ISOSURFACE_EDGE_TABLE::SetNumTableEntries";

  ISOSURFACE_TABLE::SetNumTableEntries(num_table_entries);

  edge_entry = new ISOSURFACE_EDGE_TABLE_ENTRY[num_table_entries];
  if (edge_entry == NULL)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for isosurface edge table.");
}

bool ISOSURFACE_EDGE_TABLE::CheckTable(ERROR & error_msg) const
  // check table
{
  if (!ISOSURFACE_TABLE::CheckTable(error_msg)) return(false);

  if (edge_entry == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage
      ("Memory for isosurface table edge entries not allocated.");
    return(false);
  }

  for (int it = 0; it < NumTableEntries(); it++)
    if (!edge_entry[it].Check(error_msg)) {
      error_msg.AddMessage
	("Error detected at isosurface table entry ", it, ".");
      return(false);
    }

  for (int jt = 0; jt < NumTableEntries(); jt++)
    for (int ie = 0; ie < NumEdges(jt); ie++)
      for (int iend = 0; iend < 2; iend++) {
        int iso_v = int(EdgeEndpoint(jt, ie, iend));
        if (iso_v < 0 || iso_v >= NumIsosurfaceVertices()) {
          error_msg.ClearAll();
          error_msg.AddMessage
            ("Illegal isosurface vertex ", iso_v, 
	     " in isosurface table edge entry ", jt, "."); 
          return(false);
        }
      };

  return(true);
}


bool ISOSURFACE_EDGE_TABLE::Check(ERROR & error_msg) const
{
  if (!Polytope().Check(error_msg)) return(false);
  if (!CheckTable(error_msg)) return(false);
  return(true);
}


void ISOSURFACE_EDGE_TABLE::FreeAll()
  // free all memory
{
  if (edge_entry != NULL) {
    for (int i = 0; i < num_table_entries; i++) {
      edge_entry[i].FreeAll();
    }
    delete [] edge_entry;
    edge_entry = NULL;  
  };

  ISOSURFACE_TABLE::FreeAll();
}


namespace {
  class EDGE_CONTAINER { 
  public:
    EDGE_INDEX v[2]; 

    EDGE_CONTAINER(){};
    bool operator < (const EDGE_CONTAINER & e) const {
      if (v[0] < e.v[0]) { return(true); }
      else if (v[0] == e.v[0] && v[1] < e.v[1]) { return(true); }
      else {return(false); }
    }
  };

}

void ISOSURFACE_EDGE_TABLE::GenEdgeLists()
  // for each table entry, generate edge lists from simplex lists
{
  PROCEDURE_ERROR error("ISOSURFACE_EDGE_TABLE::GenEdgeLists");
  IJK::ARRAY<EDGE_INDEX> vlist(Dimension());
  typedef set<EDGE_CONTAINER> EDGE_SET;
  EDGE_SET eset;

  if (!IsTableAllocated()) {
    error.AddMessage("Programming error: Isosurface table not allocated.");
    throw error;
  };

  for (int it = 0; it < NumTableEntries(); it++) {
    eset.clear();

    for (int is = 0; is < NumSimplices(it); is++) {

      // get simplex vertices
      for (int iv = 0; iv < Dimension(); iv++) {
        vlist[iv] = SimplexVertex(it, is, iv);
      }
      sort(vlist.Ptr(), vlist.Ptr()+Dimension());

      // store simplex edges
      for (int i0 = 0; i0 < Dimension(); i0++)
        for (int i1 = i0+1; i1 < Dimension(); i1++) {
          EDGE_CONTAINER e;
          e.v[0] = vlist[i0];
          e.v[1] = vlist[i1];

          eset.insert(e);
        }
    }

    edge_entry[it].FreeAll();
    if (eset.size() > 0) {
      edge_entry[it].edge_endpoint_list = new EDGE_INDEX[2*eset.size()];
      if (edge_entry[it].edge_endpoint_list == NULL) {
        error.AddMessage("Unable to allocate memory for edge list in table entry ", it, ".");
        throw error;
      };
      edge_entry[it].num_edges = eset.size();

      int ie = 0;
      for (EDGE_SET::iterator edge_iter = eset.begin(); 
           edge_iter != eset.end(); 
           edge_iter++) {
        edge_entry[it].edge_endpoint_list[2*ie] = edge_iter->v[0];
        edge_entry[it].edge_endpoint_list[2*ie+1] = edge_iter->v[1];
        ie++;
      }
    }
  }

}


// *******************************************************************
// UTILITY FUNCTIONS
// *******************************************************************

unsigned long IJKMCUBE_TABLE::calculate_num_entries
(const int num_vert, const int num_colors)
  // calculate num table entries = (num_colors)^(num_vert)
{
  IJK::PROCEDURE_ERROR error("calculate_num_entries");

  unsigned long num_table_entries = 1;
  IJK::int_power(num_colors, num_vert, num_table_entries, error);

  return(num_table_entries);
}

