/*!
 *  \file ijkdual_datastruct.cpp
 *  @brief ijkdual data structures.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2023 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more detailrans.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <assert.h>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "ijkcoord.tpp"

#include "ijkdual_types.h"
#include "ijkdual_datastruct.h"
#include "ijkmesh.tpp"

using namespace IJK;
using namespace IJKDUAL;


// *******************************************************************
// ISOPOLY_INFO_TYPE MEMBER FUNCTIONS
// *******************************************************************

// Initialize ISOPOLY_INFO_TYPE.
void ISOPOLY_INFO_TYPE::Init()
{
  grid_edge_isov_position_coef = 0.5;
  is_diagonal_in_envelope[0] = true;
  is_diagonal_in_envelope[1] = true;
}

// *******************************************************************
// CLASS RANDOM_POS_PARAM
// *******************************************************************

// Initialize RANDOM_POS_PARAM
void RANDOM_POS_PARAM::Init()
{
  distribution = UNDEFINED_DISTRIBUTION;
  random_seed = 0;
  flag_gen_isov_apply_offset = true;
  flag_gen_isov_separated_by_cube_center = true;
  position_offset = 0.0;
  flag_gen_isov_on_cube_boundary = false;
  boundary_width = DefaultBoundaryWidth();
}


// *******************************************************************
// CLASS DUALISO_POSITION_TIME
// *******************************************************************

// constructor
IJKDUAL::DUALISO_POSITION_TIME::DUALISO_POSITION_TIME()
{
  Clear();
}

void IJKDUAL::DUALISO_POSITION_TIME::Clear()
{
  basic = 0.0;
  separate_by_cube_center = 0.0;
  separate_isov_near_cube_facets = 0.0;
  move_isov_away_from_cube_faces = 0.0;
  isov_dual_to_isopoly = 0.0;
  total = 0.0;
}


void IJKDUAL::DUALISO_POSITION_TIME::Add
(const DUALISO_POSITION_TIME & position_time)
{
  basic += position_time.basic;
  separate_by_cube_center += position_time.separate_by_cube_center;
  separate_isov_near_cube_facets += position_time.separate_isov_near_cube_facets;
  move_isov_away_from_cube_faces += position_time.move_isov_away_from_cube_faces;
  isov_dual_to_isopoly += position_time.isov_dual_to_isopoly;
  total += position_time.total;
}


// *******************************************************************
// CLASS DUALISO_TIME
// *******************************************************************

// constructor
IJKDUAL::DUALISO_TIME::DUALISO_TIME()
{
  Clear();
}

void IJKDUAL::DUALISO_TIME::Clear()
{
  preprocessing = 0.0;
  extract = 0.0;
  merge = 0.0;
  split_isov = 0.0;
  position.Clear();
  collapse = 0.0;
  envelope = 0.0;
  triangulation = 0.0;
  total = 0.0;
}

void IJKDUAL::DUALISO_TIME::Add(const DUALISO_TIME & dualiso_time)
{
  preprocessing += dualiso_time.preprocessing;
  extract += dualiso_time.extract;
  merge += dualiso_time.merge;
  split_isov += dualiso_time.split_isov;
  position.Add(dualiso_time.position);
  collapse += dualiso_time.collapse;
  envelope += dualiso_time.envelope;
  triangulation += dualiso_time.triangulation;  
  total += dualiso_time.total;
}

// **************************************************
// INFO CLASSES
// **************************************************

IJKDUAL::GRID_INFO::GRID_INFO()
{
  Clear();
}

void IJKDUAL::GRID_INFO::Clear()
{
  num_cubes = 0;
}

void IJKDUAL::SCALAR_INFO::Init(const int dimension)
{
  this->dimension = 0;
  SetDimension(dimension);

  Clear();
}

void IJKDUAL::SCALAR_INFO::FreeAll()
{
  dimension = 0;
}

void IJKDUAL::SCALAR_INFO::Clear()
{
  num_non_empty_cubes = 0;
  num_bipolar_edges = 0;
}

void IJKDUAL::SCALAR_INFO::SetDimension(const int dimension)
{
  FreeAll();

  this->dimension = dimension;

  Clear();
}

void IJKDUAL::SCALAR_INFO::Copy(const SCALAR_INFO & info)
{
  Init(info.Dimension());
  num_non_empty_cubes = info.num_non_empty_cubes;
  num_bipolar_edges = info.num_bipolar_edges;
}

/// Copy assignment.
const SCALAR_INFO &  IJKDUAL::SCALAR_INFO::operator =
(const SCALAR_INFO & right)
{
  if (&right != this) {
    FreeAll();
    Copy(right);
  }

  return *this;
}

IJKDUAL::SCALAR_INFO::~SCALAR_INFO()
{
  dimension = 0;
  Clear();
}

IJKDUAL::ISOV_INFO::ISOV_INFO()
{
  Clear();
}

void IJKDUAL::ISOV_INFO::Clear()
{
  num_times_isov_moved_away_from_grid_cube_facets = 0;
}

IJKDUAL::MULTI_ISOV_INFO::MULTI_ISOV_INFO()
{
  Clear();
}

void IJKDUAL::MULTI_ISOV_INFO::Clear()
{
  num_cubes_single_isov = 0;
  num_cubes_multi_isov = 0;
  num_non_manifold_split = 0;
  num_1_2_changed = 0;
  num_connect_changed = 0;
  num_ambig_ridge_cubes_changed = 0;
  num_non_ambig_ridge_cubes_changed = 0;
}

IJKDUAL::TRIANGULATION_INFO::TRIANGULATION_INFO()
{
  Clear();
}

void IJKDUAL::TRIANGULATION_INFO::Clear()
{
  num_iso_cubes_tri_total = 0;
  num_iso_cubes_tri_no_add = 0;
  num_iso_cubes_tri_add_interior1 = 0;
  num_iso_cubes_with_diag_outside_envelope = 0;
  num_tri_add_interior1_with_diag_outside_envelope = 0;
  num_tri_no_add_with_diag_outside_envelope = 0;
}

IJKDUAL::DUALISO_INFO::DUALISO_INFO()
{
  Clear();
}

IJKDUAL::DUALISO_INFO::DUALISO_INFO(const int dimension):scalar(dimension)
{
  Clear();
}

void IJKDUAL::DUALISO_INFO::Clear()
{
  grid.Clear();
  scalar.Clear();
  time.Clear();
  isov.Clear();
  multi_isov.Clear();
  triangulation.Clear();
}


// *******************************************************************
// MERGE DATA
// *******************************************************************

void IJKDUAL::MERGE_DATA::Init
(const int dimension, const AXIS_SIZE_TYPE * axis_size,
 const MERGE_INDEX num_obj_per_vertex, const MERGE_INDEX num_obj_per_edge)
{
  this->num_obj_per_vertex = num_obj_per_vertex;
  this->num_obj_per_edge = num_obj_per_edge;
  this->num_obj_per_grid_vertex = 
    dimension*num_obj_per_edge + num_obj_per_vertex;
  compute_num_grid_vertices(dimension, axis_size, num_vertices);
  num_edges = dimension*num_vertices;
  vertex_id0 = num_obj_per_edge*num_edges;
  MERGE_INDEX num_obj = 
    num_obj_per_vertex*num_vertices + num_obj_per_edge*num_edges;
  INTEGER_LIST<MERGE_INDEX,MERGE_INDEX>::Init(num_obj);
}

bool IJKDUAL::MERGE_DATA::Check(ERROR & error) const
{
  if (MaxNumInt() < 
      NumObjPerVertex()*NumVertices() + NumObjPerEdge()*NumEdges()) {
    error.AddMessage("Not enough allocated memory.");
    return(false);
  };

  return(true);
}

