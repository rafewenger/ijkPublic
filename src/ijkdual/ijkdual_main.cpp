/*!
 *  @file ijkdual_main.cpp
 *  @brief Generate isosurface using dual contouring algorithm.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2024 Rephael Wenger

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


#include <iostream>
#include <new>

#include "ijkdualIO.h"

#include "ijkdual.tpp"
#include "ijkdual_triangulate.tpp"

using namespace IJK;
using namespace IJKDUAL;

using namespace std;

// local subroutines
void memory_exhaustion();
void construct_dual_isosurface
(const IO_INFO & io_info, const DUALISO_DATA & dualiso_data,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);

  DUALISO_TIME dualiso_time;
  IO_TIME io_time = {0.0, 0.0};
  IO_INFO io_info;
  IJK::ERROR error;

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv, io_info);

    DUALISO_SCALAR_GRID full_scalar_grid;
    NRRD_HEADER nrrd_header;

    ELAPSED_TIME wall_time;
    read_nrrd_file(io_info.input_filename, full_scalar_grid,  nrrd_header);
    io_time.read_nrrd_time = wall_time.getElapsed();

    check_and_reset_input(full_scalar_grid, io_info);

    // *** Spacing already read in read_nrrd_file?!?
    nrrd_header.GetSpacing(io_info.grid_spacing);

    // set DUAL datastructures and flags
    DUALISO_DATA dualiso_data;

    // Note: dualiso_data.SetScalarGrid must be called before set_mesh_data.
    dualiso_data.SetScalarGrid
      (full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
       io_info.flag_supersample, io_info.supersample_resolution);
    dualiso_data.Set(io_info);
    warn_non_manifold(io_info);
    report_num_cubes(full_scalar_grid, io_info, dualiso_data);

    construct_dual_isosurface
      (io_info, dualiso_data, dualiso_time, io_time);

    if (io_info.flag_report_time) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cout << endl;
      report_time(io_info, io_time, dualiso_time, total_elapsed_time);
    };

  } 
  catch (ERROR & error) {
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

// forward declaration
template <typename DUALISO_DATA_TYPE, typename DUAL_ISOSURFACE_TYPE,
          typename DUAL_ISOVERT_TYPE, typename ISOPOLY_INFO_TYPE>
void rescale_and_triangulate
(const DUALISO_DATA_TYPE & dualiso_data,
 const SCALAR_TYPE isovalue, DUAL_ISOSURFACE_TYPE & dual_isosurface,
 const std::vector<DUAL_ISOVERT_TYPE> & isov_list,
 std::vector<ISOPOLY_INFO_TYPE> & isopoly_info, 
 DUALISO_INFO & dualiso_info);


/*!
 *  @brief Construct dual triangulated isosurface.
 */
void construct_dual_triangulated_isosurface
(const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue, 
 OUTPUT_INFO & output_info,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time)
{
  const int dimension = dualiso_data.ScalarGrid().Dimension();
  const int num_facet_vertices = dualiso_data.ScalarGrid().NumFacetVertices();
  const int num_cubes = dualiso_data.ScalarGrid().ComputeNumCubes();
  const bool allow_multiple_isov =
    dualiso_data.AllowMultipleIsoVertices();  

  DUALISO_INFO dualiso_info(dimension);
  dualiso_info.grid.num_cubes = num_cubes;

  // Dual contouring.  

  std::clock_t t_start = std::clock();

  DUAL_TRIANGULATED_ISOSURFACE dual_isosurface(dimension, num_facet_vertices);
  
  std::vector<ISOPOLY_INFO_TYPE> isopoly_info;

  if (allow_multiple_isov) {
    // isov_list[kw] = Information on isosurface vertex kw.
    // Includes cube index, lookup table index,
    //   and isosurface patch index.
    std::vector<IJKDUAL::DUAL_ISOVERT> isov_list;

    dual_contouring_multi_isov
      (dualiso_data, isovalue, dual_isosurface,
       isov_list, isopoly_info, dualiso_info);

    rescale_and_triangulate
      (dualiso_data, isovalue, dual_isosurface,
       isov_list, isopoly_info, dualiso_info);
  }
  else {
    // Only one isosurface vertex per grid cube.
    // active_cube_list[kw] =
    //   Index of cube containing isosurface vertex kw.
    std::vector<CUBE_INDEX> active_cube_list;

    dual_contouring_single_isov
      (dualiso_data, isovalue, dual_isosurface,
       active_cube_list, isopoly_info, dualiso_info);

    rescale_and_triangulate
      (dualiso_data, isovalue, dual_isosurface,
       active_cube_list, isopoly_info, dualiso_info);    
  }

  
  // store times
  std::clock_t t_end = std::clock();
  IJK::clock2seconds(t_end-t_start, dualiso_info.time.total);

  dualiso_time.Add(dualiso_info.time);

  // *** SHOULD REALLY JUST USE max x coordinate. ***
  const COORD_TYPE max_coord =
    get_max_abs_array_value(dual_isosurface.vertex_coord);
  output_info.fig_param.coord_scale_factor =
    output_info.fig_param.DefaultPaperDrawRegionWidth()/max_coord;

  if (dualiso_data.mesh_type == MIXED_MESH) {
    if (dimension == DIM3) {
      // Remove any quadrilaterals that are collapsed
      //   to a single vertex representing deletion of the quad.
      delete_collapsed_quads(dual_isosurface.isopoly_vert);
    }
  }

  output_dual_isosurface
    (output_info, dualiso_data, dual_isosurface, dualiso_info, io_time);
}


/*!
 *  @brief Construct dual isosurface.
 *  - Outputs cubes (line segments, quads, hexahedra, ...)
 */
void construct_dual_cube_complex_isosurface
(const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue, 
 OUTPUT_INFO & output_info,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time)
{
  const int dimension = dualiso_data.ScalarGrid().Dimension();
  const int num_facet_vertices = dualiso_data.ScalarGrid().NumFacetVertices();
  const int num_cubes = dualiso_data.ScalarGrid().ComputeNumCubes();

  DUALISO_INFO dualiso_info(dimension);
  dualiso_info.grid.num_cubes = num_cubes;

  // Dual contouring.  

  DUAL_ISOSURFACE dual_isosurface(dimension, num_facet_vertices);

  dual_contouring(dualiso_data, isovalue, dual_isosurface, dualiso_info);

  dualiso_time.Add(dualiso_info.time);

  rescale_vertex_coord
    (dimension, dualiso_data.ScalarGrid().SpacingPtrConst(),
     dual_isosurface.vertex_coord);

  // *** SHOULD REALLY JUST USE max x coordinate. ***
  const COORD_TYPE max_coord =
    get_max_abs_array_value(dual_isosurface.vertex_coord);
  output_info.fig_param.coord_scale_factor =
    output_info.fig_param.DefaultPaperDrawRegionWidth()/max_coord;

  output_dual_cube_complex_isosurface
    (output_info, dualiso_data, dual_isosurface, dualiso_info, io_time);
}


void construct_dual_isosurface
(const IO_INFO & io_info, const DUALISO_DATA & dualiso_data,
 DUALISO_TIME & dualiso_time, IO_TIME & io_time)
{
  const int dimension = dualiso_data.ScalarGrid().Dimension();

  io_time.write_time = 0;
  for (unsigned int i = 0; i < io_info.isovalue.size(); i++) {

    const SCALAR_TYPE isovalue = io_info.isovalue[i];

    OUTPUT_INFO output_info;
    set_output_info(io_info, dimension, i, output_info);

    if (dimension == DIM3 && 
        (dualiso_data.mesh_type == SIMPLICIAL_COMPLEX ||
         dualiso_data.mesh_type == MIXED_MESH)) {
      construct_dual_triangulated_isosurface
        (dualiso_data, isovalue, output_info, dualiso_time, io_time);
    }
    else {
      construct_dual_cube_complex_isosurface
        (dualiso_data, isovalue, output_info, dualiso_time, io_time);
    }
  }
}


template <typename DUALISO_DATA_TYPE, typename DUAL_ISOSURFACE_TYPE,
          typename DUAL_ISOVERT_TYPE, typename ISOPOLY_INFO_TYPE>
void rescale_and_triangulate
(const DUALISO_DATA_TYPE & dualiso_data,
 const SCALAR_TYPE isovalue, DUAL_ISOSURFACE_TYPE & dual_isosurface,
 const std::vector<DUAL_ISOVERT_TYPE> & isov_list,
 std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
 DUALISO_INFO & dualiso_info)
{
  const int dimension = dualiso_data.ScalarGrid().Dimension();

  if (dualiso_data.flag_check_envelope) {
    dual_isosurface.SetQuadVertexOrder(CIRCULAR_VERTEX_ORDER);
    compute_isoquad_diagonals_deep_in_envelope
      (dualiso_data.ScalarGrid(), dual_isosurface, isov_list,
       dualiso_data, isopoly_info, dualiso_info);
  }

  
  if (dualiso_data.flag_tri4_quad ||
      (dualiso_data.flag_check_envelope &&
       dualiso_data.flag_allow_tri4_envelope)) {

    std::clock_t t0 = std::clock();

    bool flag_add_all_tri4_interior_vertices = false;

    if (dualiso_data.flag_force_add_all_interior_vertices)
      { flag_add_all_tri4_interior_vertices = true; }
    else if (dualiso_data.flag_check_envelope &&
             dualiso_data.flag_allow_tri2_envelope &&
             dualiso_data.flag_allow_tri4_envelope &&
             (dualiso_data.envelope_quad_tri_method ==
              ENVELOPE_MAX_MIN_ANGLE))
      { flag_add_all_tri4_interior_vertices = true; }

    else if (dualiso_data.mesh_type == SIMPLICIAL_COMPLEX) {
      if ((dualiso_data.quad_tri_method == TRI4_MAX_MIN_ANGLE) ||
          (dualiso_data.quad_tri_method == TRI4_ALL_QUADS))
        { flag_add_all_tri4_interior_vertices = true; }
    }

    if (flag_add_all_tri4_interior_vertices) {
      add_tri4_interior_vertices
        (dualiso_data, isovalue, isopoly_info, dual_isosurface);
    }

    std::clock_t t1 = std::clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.position.isov_dual_to_isopoly);
    dualiso_info.time.position.total +=
      dualiso_info.time.position.isov_dual_to_isopoly;
  }

  std::clock_t t0 = std::clock();

  // Rescale before triangulation.
  rescale_vertex_coord
    (dimension, dualiso_data.ScalarGrid().SpacingPtrConst(),
     dual_isosurface.vertex_coord);

  if (dimension == DIM3) {

    if (dualiso_data.flag_check_envelope) {
      int num_tri2_split = 0;
      int num_tri4_split = 0;
      if (dualiso_data.flag_allow_tri2_envelope) {
        if (dualiso_data.flag_allow_tri4_envelope) {

          if (dualiso_data.envelope_quad_tri_method ==
              ENVELOPE_PREFER_TRI2) {
            tri2_or_tri4_quads_with_diagonals_outside_envelope_prefer_tri2
              (dualiso_data, isovalue, isopoly_info, dual_isosurface,
               num_tri2_split, num_tri4_split);
          }
          else {
            tri2_or_tri4_quads_with_diagonals_outside_envelope
              (dualiso_data, isopoly_info, dual_isosurface,
               num_tri2_split, num_tri4_split);
          }

          dualiso_info.triangulation.num_tri_no_add_with_diag_outside_envelope =
            num_tri2_split;
          dualiso_info.triangulation.num_tri_add_interior1_with_diag_outside_envelope =
            num_tri4_split;
        }
        else {
          tri2_dual_quad_list_with_exactly_one_diagonal_outside_envelope
            (dualiso_data, isopoly_info, dual_isosurface, num_tri2_split);

          dualiso_info.triangulation.num_tri_no_add_with_diag_outside_envelope =
            num_tri2_split;
        }
      }
      else if (dualiso_data.flag_allow_tri4_envelope) {
        tri4_quads_with_diagonals_outside_envelope
          (dualiso_data, isovalue, isopoly_info,
           dual_isosurface, num_tri4_split);

        dualiso_info.triangulation.num_tri_add_interior1_with_diag_outside_envelope =
            num_tri4_split;
      }

      dualiso_info.triangulation.num_iso_cubes_tri_no_add +=
        num_tri2_split;
      dualiso_info.triangulation.num_iso_cubes_tri_add_interior1 +=
        num_tri4_split;
      dualiso_info.triangulation.num_iso_cubes_tri_total +=
        num_tri2_split + num_tri4_split;
    }

    if (dualiso_data.mesh_type == SIMPLICIAL_COMPLEX) {
      triangulate_isosurface_quadrilaterals
        (dualiso_data, isovalue, isopoly_info, dual_isosurface,
         dualiso_info.triangulation);
    }

    // Reset to COORDINATE_VERTEX_ORDER, if necessary.
    dual_isosurface.SetQuadVertexOrder(COORDINATE_VERTEX_ORDER);
  }
  
  std::clock_t t1 = std::clock();
  IJK::clock2seconds(t1-t0, dualiso_info.time.triangulation);
}


void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

