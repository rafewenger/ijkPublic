/*!
 *  @file ijkmcube.cxx
 *  @brief Marching cubes/hypercubes isosurface generation.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2023 Rephael Wenger

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


#include <assert.h>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "ijkmcube.h"
#include "ijkmcube_util.h"
#include "ijkmcube_extract.h"
#include "ijkoctree.h"
#include "ijkMCtable.h"

#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijkpoly.tpp"
#include "ijkmerge.tpp"

using namespace IJK;
using namespace IJKMCUBE;
using namespace IJKMCUBE_TABLE;



// **************************************************
// MARCHING CUBES (HYPERCUBES)
// **************************************************

// Marching Cubes Algorithm.
//  - Calls NEP version, multires version, octree/minmax region versions, etc.
//    based on mc_data flags.
void IJKMCUBE::marching_cubes
(const MC_DATA & mc_data, const SCALAR_TYPE isovalue, 
 MC_ISOSURFACE & mc_isosurface, MCUBE_INFO & mcube_info)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = mc_data.ScalarGrid().AxisSize();
  float merge_time = 0.0;
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!mc_data.Check(error)) { throw error; };

  mc_isosurface.Clear();
  mcube_info.time.Clear();

  ISOSURFACE_TOPOLOGY isosurface_topology = mc_data.IsosurfaceTopology();
  INTERPOLATION_TYPE interpolation_type = mc_data.InterpolationType();

  if (mc_data.UseNEP()) {
    int nep_num_dup = mc_data.NEPNumDup();
    NEP_ISO_MERGE_DATA merge_data(dimension, axis_size);
    clock_t t1 = clock();
    merge_time = clock2seconds(t1-t_start);

    if (mc_data.UseOctree()) {
      marching_cubes_nep
        (mc_data.ScalarGrid(), mc_data.isotable.cube_nep,
         isovalue, nep_num_dup,
         mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
         merge_data, *(mc_data.Octree()), mcube_info);
    }
    else if (mc_data.UseMinmaxRegions()) {
      marching_cubes_nep
        (mc_data.ScalarGrid(), mc_data.isotable.cube_nep,
         isovalue, nep_num_dup,
         mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
         merge_data, *(mc_data.MinmaxRegions()), mcube_info);
    }
    else {
      marching_cubes_nep
        (mc_data.ScalarGrid(), mc_data.isotable.cube_nep,
         isovalue, nep_num_dup,
         mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
         merge_data, mcube_info);
    }

    mcube_info.time.extract += merge_time;
  }
  else if (mc_data.UseMultires()) {
    multires_cubes
      (mc_data.ScalarGrid(), *(mc_data.MultiresGrid()), mc_data.isotable,
       isovalue, mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
       mc_data.MergeEdgesParameters(), mcube_info);
  }
  else if (mc_data.CubeContainingSimplexFlag()) {

    ISO_MERGE_DATA merge_data(dimension, axis_size);
    clock_t t1 = clock();
    merge_time = clock2seconds(t1-t_start);

    marching_cubes
      (mc_data.ScalarGrid(), mc_data.isotable.cube,
       isovalue, mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
       mc_isosurface.cube_containing_simplex, merge_data, mcube_info);
  }
  else if (isosurface_topology == ISOTABLE_TOPOLOGY) {
    if (mc_data.UseOctree() || mc_data.UseMinmaxRegions() ||
        mc_data.EdgeRepresentation() == EDGE_ID) {

      ISO_MERGE_DATA merge_data(dimension, axis_size);
      clock_t t1 = clock();
      merge_time = clock2seconds(t1-t_start);

      if (mc_data.UseOctree()) {
        marching_cubes
          (mc_data.ScalarGrid(), mc_data.isotable.cube,
           isovalue, mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
           merge_data, *(mc_data.Octree()), mcube_info);
      }
      else if (mc_data.UseMinmaxRegions()) {
        marching_cubes
          (mc_data.ScalarGrid(), mc_data.isotable.cube,
           isovalue, mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
           merge_data, *(mc_data.MinmaxRegions()), mcube_info);
      }
      else {
        // Edge representation = EDGE_ID
        marching_cubes
          (mc_data.ScalarGrid(), mc_data.isotable.cube,
           isovalue, mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
           merge_data, mcube_info);
      }

      mcube_info.time.extract += merge_time;
    }
    else {
      // Edge representation = EDGE_ENDPOINT_PAIR
      marching_cubes
        (mc_data.ScalarGrid(), mc_data.isotable.cube,
         isovalue, mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
         mc_data.MergeEdgesParameters(), mcube_info);
    }
  }
  else if (isosurface_topology == CUBE_DECIDER_TOPOLOGY &&
           mc_data.EdgeRepresentation() == EDGE_ID) {

    ISO_MERGE_DATA merge_data(dimension, axis_size);
    clock_t t1 = clock();
    merge_time = clock2seconds(t1-t_start);

    marching_cubes
      (mc_data.ScalarGrid(), mc_data.isotable.cube,
       mc_data.isotable.ambig_cube, isovalue, isosurface_topology,
       mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
       merge_data, mcube_info);

    mcube_info.time.extract += merge_time;
  }
  else {

    marching_cubes
      (mc_data.ScalarGrid(), mc_data.isotable,
       isovalue, isosurface_topology, interpolation_type,
       mc_isosurface.simplex_vert, mc_isosurface.vertex_coord,
       mc_data.MergeEdgesParameters(), mcube_info);
  }

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = clock2seconds(t_end-t_start);
}

// **************************************************
// MARCHING CUBES INTERVAL VOLUME
// **************************************************

// Marching Cubes Interval Volume Algorithm.
void IJKMCUBE::MCVol
(const MC_DATA & mc_data, 
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 MC_ISOSURFACE & mc_ivol, MCUBE_INFO & mcube_info)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  const AXIS_SIZE_TYPE * axis_size = mc_data.ScalarGrid().AxisSize();
  PROCEDURE_ERROR error("MCVol");

  clock_t t_start = clock();

  if (!mc_data.Check(error)) { throw error; };

  mc_ivol.Clear();
  mcube_info.time.Clear();

  IVOL_MERGE_DATA merge_data(dimension, axis_size);
  MCVol
    (mc_data.ScalarGrid(), mc_data.isotable.cube,
     isovalue0, isovalue1, mc_ivol.simplex_vert, mc_ivol.vertex_coord,
     merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = clock2seconds(t_end-t_start);
}


// **************************************************
// MARCHING CUBES (HYPERCUBES)
// **************************************************

// Extract isosurface using Marching Cubes algorithm.
// - Returns list of isosurface simplex vertices.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);
  

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


// Extract isosurface using Marching Cubes algorithm.
// - Returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 std::vector<VERTEX_INDEX> & cube_list,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  cube_list.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_and_cube_info
    (scalar_grid, isotable, isovalue, iso_simplices, cube_list, mcube_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);
  

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}



// Extract isosurface using Marching Cubes algorithm.
// - Slower, but less memory intensive version.
// - Represents edges containing isosurface vertices as pairs of grid vertices.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & iso_simplices,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_edges_parameters, 
 MCUBE_INFO & mcube_info)
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  iso_simplices.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  // (endpoint[2*i], endpoint[2*i+1]) = endpoints of edge containing i'th isosurface vertex.
  std::vector<VERTEX_INDEX> endpoint;
  extract_iso_simplices
    (scalar_grid, isotable, isovalue, iso_simplices, endpoint,
     mcube_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, endpoint, BINARY, iso_simplices,
     vertex_coord, merge_edges_parameters, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


// Extract isosurface using Marching Cubes algorithm.
// - Version that does not return MCUBE_INFO.
// - Version that creates internal data structure to handle edge merging.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  MCUBE_INFO mcube_info;
  ISO_MERGE_DATA merge_data(dimension, axis_size);
  marching_cubes(scalar_grid, isotable, isovalue, 
                 simplex_vert, vertex_coord, merge_data, mcube_info);
}


// Marching Cubes Algorithm with specified isosurface topology and
//   specified interpolation to compute new mesh vertex scalar values.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const POLY_ISOTABLE & poly_isotable,
 const SCALAR_TYPE isovalue,
 const ISOSURFACE_TOPOLOGY isosurface_topology,
 const INTERPOLATION_TYPE interpolation_type,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)
{
  const int dimension = poly_isotable.cube.Dimension();
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding
      (poly_isotable.cube, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  // (endpoint[2*i], endpoint[2*i+1]) = endpoints of edge containing i'th isosurface vertex
  std::vector<VERTEX_INDEX> endpoint;
  MC_MESH_VERTEX_LIST new_mesh_vertices(dimension);

  extract_iso_simplices
    (scalar_grid, poly_isotable, isovalue, isosurface_topology,
     simplex_vert, endpoint, new_mesh_vertices, mcube_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, new_mesh_vertices, endpoint, BINARY, 
     interpolation_type, simplex_vert, vertex_coord, 
     merge_parameters, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

// Marching Cubes Algorithm.
// @param isosurface_topology = Isosurface topology.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const POLY_ISOTABLE & poly_isotable,
 const SCALAR_TYPE isovalue,
 const ISOSURFACE_TOPOLOGY isosurface_topology,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)
{
  marching_cubes(scalar_grid, poly_isotable, isovalue,
                 isosurface_topology, LINEAR_INTERPOLATION,
                 simplex_vert, vertex_coord, merge_parameters, mcube_info);
}


// Marching Cubes Algorithm.
// @param isosurface_topology = Isosurface topology.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const ISOSURFACE_TABLE & cube_isotable,
 const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info,
 const SCALAR_TYPE isovalue,
 const ISOSURFACE_TOPOLOGY isosurface_topology,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)
{
  const int dimension = cube_isotable.Dimension();
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(cube_isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  if (!check_isotable_ambig_info(cube_isotable, ambig_info, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  // (endpoint[2*i], endpoint[2*i+1]) = endpoints of edge containing i'th isosurface vertex
  std::vector<VERTEX_INDEX> endpoint;

  switch(isosurface_topology) {

  case ISOTABLE_TOPOLOGY:
    extract_iso_simplices
      (scalar_grid, cube_isotable, isovalue, simplex_vert, endpoint, 
       mcube_info);
    break;

  case CUBE_DECIDER_TOPOLOGY:
    extract_iso_simplices_cube_decider
      (scalar_grid, cube_isotable, ambig_info, isovalue, 
       simplex_vert, endpoint, mcube_info);
    break;

  case ASYMPTOTIC_DECIDER_TOPOLOGY:
  case LINEAR_TOPOLOGY:
    if (dimension <= 2) {
      extract_iso_simplices_cube_decider
        (scalar_grid, cube_isotable, ambig_info, isovalue, 
         simplex_vert, endpoint, mcube_info);
    }
    else {
      if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY) {
        error.AddMessage("Programming error.  Asymptotic decider requires pyramid isosurface lookup table.");
        throw error;
      }
      else {
        error.AddMessage("Programming error.  Linear topology requires pyramid isosurface lookup table.");
        throw error;
      }
    }
    break;

  default:
    error.AddMessage("Programming error. Unknown isosurface topology.");
    throw error;
  };

  merge_and_position_vertices
    (scalar_grid, isovalue, endpoint, BINARY, 
     simplex_vert, vertex_coord, merge_parameters, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


// Marching Cubes Algorithm.
// @param isosurface_topology = Isosurface topology.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, 
 const ISOSURFACE_TABLE & cube_isotable,
 const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info,
 const SCALAR_TYPE isovalue,
 const ISOSURFACE_TOPOLOGY isosurface_topology,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
{
  const int dimension = cube_isotable.Dimension();
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(cube_isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  switch(isosurface_topology) {

  case ISOTABLE_TOPOLOGY:
    extract_iso_simplices
      (scalar_grid, cube_isotable, isovalue, iso_simplices, mcube_info);
    break;

  case CUBE_DECIDER_TOPOLOGY:
    extract_iso_simplices_cube_decider
      (scalar_grid, cube_isotable, ambig_info, isovalue, 
       iso_simplices, mcube_info);
    break;

  case ASYMPTOTIC_DECIDER_TOPOLOGY:
  case LINEAR_TOPOLOGY:
    if (dimension <= 2) {
      extract_iso_simplices_cube_decider
        (scalar_grid, cube_isotable, ambig_info, isovalue, 
         iso_simplices, mcube_info);
    }
    else {
      if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY) {
        error.AddMessage("Programming error.  Asymptotic decider requires pyramid isosurface lookup table.");
        throw error;
      }
      else {
        error.AddMessage("Programming error.  Linear topology requires pyramid isosurface lookup table.");
        throw error;
      }
    }
    break;

  default:
    error.AddMessage("Programming error. Unknown isosurface topology.");
    throw error;
  };

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

// Marching Cubes Algorithm.
// - Version using minmax to extract isosurface simplices.
// @param isosurface_topology = Isosurface topology.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
 MCUBE_INFO & mcube_info)
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;
  extract_iso_simplices_from_minmax
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, minmax);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


// Marching Cubes Algorithm.
// - Version using octree to extract isosurface simplices.
// @param isosurface_topology = Isosurface topology.
void IJKMCUBE::marching_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree, 
 MCUBE_INFO & mcube_info)
{
  PROCEDURE_ERROR error("marching_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_from_octree
    (scalar_grid, isotable, isovalue, iso_simplices, mcube_info, octree);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, BINARY,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


// **************************************************
// NEP MARCHING CUBES (NEP: NEGATIVE-EQUALS-POSITIVE)
// **************************************************


// Extract isosurface using Marching Cubes algorithm.
// - Use NEP (negative-equals-positive) isosurface lookup table.
// - Returns list of isosurface simplex vertices
//   and list of isosurface vertex coordinates.
// nep_num_dup = nep control for number of duplicate isosurface simplices
//    - 0 : Do not create duplicate isosurface simplices
//    - 1 : Merge duplicate isosurface simplices into one simplex
//    - 2 : Allow duplicate isosurface simplices
void IJKMCUBE::marching_cubes_nep
(const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
{
  PROCEDURE_ERROR error("marching_cubes_nep");

  if (!check_nep_num_dup(nep_num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BASE3, error))
    { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, mcube_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, NEP,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


void IJKMCUBE::marching_cubes_nep
(const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const MINMAX_REGIONS & minmax, 
 MCUBE_INFO & mcube_info)
  // same as previous function but uses minmax to extract isosurface simplices
{
  PROCEDURE_ERROR error("marching_cubes_nep");

  if (!check_nep_num_dup(nep_num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BASE3, error))
    { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_from_minmax_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, mcube_info,
     minmax);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, NEP,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}


void IJKMCUBE::marching_cubes_nep
(const MC_SCALAR_GRID_BASE & scalar_grid, const NEP_ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue, const int nep_num_dup,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, const IJKOCTREE::OCTREE & octree,
 MCUBE_INFO & mcube_info)
  // same as previous function but uses minmax to extract isosurface simplices
{
  PROCEDURE_ERROR error("marching_cubes_nep");

  if (!check_nep_num_dup(nep_num_dup, error)) { throw error; }
  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BASE3, error))
    { throw error; };

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  extract_iso_simplices_from_octree_nep
    (scalar_grid, isotable, isovalue, nep_num_dup, iso_simplices, mcube_info,
     octree);

  merge_and_position_vertices
    (scalar_grid, isovalue, iso_simplices, NEP,
     simplex_vert, vertex_coord, merge_data, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = float(t_end-t_start)/CLOCKS_PER_SEC;
}

// **************************************************
// MARCHING CUBES INTERVAL VOLUME
// **************************************************

void IJKMCUBE::MCVol
(const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 std::vector<VERTEX_INDEX> & simplex_vert, 
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, 
 MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();

  clock_t t_start = clock();

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  clock_t t0 = clock();
  std::vector<ISO_VERTEX_INDEX> iso_simplices;

  VERTEX_INDEX num_mixed_cubes;
  extract_ivol_simplices
    (scalar_grid, isotable, isovalue0, isovalue1, 
     iso_simplices, num_mixed_cubes);
  mcube_info.scalar.num_non_empty_cubes = num_mixed_cubes;
  assert(iso_simplices.size()%isotable.NumVerticesPerSimplex() == 0);
  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> ivol_vlist;
  merge_identical
    (iso_simplices, ivol_vlist, simplex_vert, merge_data);
  clock_t t2 = clock();

  // num_ivolv = number of interval volume vertices
  int num_ivolv = ivol_vlist.size();

  int numc = num_ivolv*dimension;
  vertex_coord.resize(numc);
  position_ivol_vertices_linear
    (scalar_grid, isovalue0, isovalue1, ivol_vlist, &(vertex_coord[0]));
  clock_t t3 = clock();

  // store times
  mcube_info.time.extract = float(t1-t0)/CLOCKS_PER_SEC;
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
  mcube_info.time.total = float(t3-t_start)/CLOCKS_PER_SEC;
}


// **************************************************
// MULTIRES CUBES
// **************************************************

void IJKMCUBE::multires_cubes
(const MC_SCALAR_GRID_BASE & scalar_grid,  const MULTIRES_GRID & multires_grid, 
 const POLY_ISOTABLE & poly_isotable,
 const SCALAR_TYPE isovalue, 
 std::vector<VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)
{
  PROCEDURE_ERROR error("multires_cubes");

  clock_t t_start = clock();

  if (!check_isotable_encoding(poly_isotable, IJKMCUBE_TABLE::BINARY, error))
    { throw error; };

  if (!multires_grid.IsProcessed()) {
    error.AddMessage("Call multires_grid.Process() before calling multires_cubes.");
    throw error;
  }

  simplex_vert.clear();
  vertex_coord.clear();
  mcube_info.time.Clear();

  // (endpoint[2*i], endpoint[2*i+1]) = endpoints of edge containing i'th isosurface vertex
  std::vector<VERTEX_INDEX> endpoint;
  multires_extract
    (scalar_grid, multires_grid, poly_isotable, isovalue, 
     simplex_vert, endpoint, mcube_info);

  merge_and_position_vertices
    (scalar_grid, isovalue, endpoint, BINARY, simplex_vert,
     vertex_coord, merge_parameters, mcube_info);

  // store times
  clock_t t_end = clock();
  mcube_info.time.total = clock2seconds(t_end-t_start);
}


// **************************************************
// ERROR CHECKING
// **************************************************

bool IJKMCUBE::check_nep(const ISOSURFACE_TABLE & isotable, 
                         const MERGE_DATA & merge_data,
                         const int num_dup, ERROR & error)
{
  if (!check_nep_num_dup(num_dup, error)) { return false; }
  if (!check_isotable_encoding(isotable, IJKMCUBE_TABLE::BASE3, error))
    { return false; };
  if (!merge_data.Check(NEP, error)) { return false; };

  return true;
}
