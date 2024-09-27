/*!
 *  @file ijkmcube.cxx
 *  @brief Marching Cubes algorithm to generate isosurface from scalar field.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2001-2023 Rephael Wenger

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

#include "ijkmcubeIO.h"


using namespace IJK;
using namespace IJKMCUBE;

using namespace std;

// local subroutines
void memory_exhaustion();
void construct_isosurface
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time);
void construct_interval_volume
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  time_t start_time;
  time(&start_time);
  string isotable_filename;

  MCUBE_TIME mcube_time;
  IO_TIME io_time = {0.0, 0.0, 0.0};
  IO_INFO io_info;
  IJK::ERROR error;

  try {

    std::set_new_handler(memory_exhaustion);

    get_isotable_directory(io_info.isotable_directory);
    parse_command_line(argc, argv, io_info);

    MC_SCALAR_GRID full_scalar_grid;
    NRRD_HEADER nrrd_header;

    ELAPSED_TIME wall_time;
    read_nrrd_file
      (io_info.input_filename, full_scalar_grid,  nrrd_header);
    io_time.read_nrrd_time = wall_time.getElapsed();

    check_input(io_info, full_scalar_grid);

    // *** Spacing already read in read_nrrd_file?!?
    nrrd_header.GetSpacing(io_info.grid_spacing);

    // Set mc datastructures and flags.
    MC_DATA mc_data;

    // Note: mc_data.SetScalarGrid must be called before set_mc_data.
    mc_data.SetScalarGrid
      (full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
       io_info.flag_supersample, io_info.supersample_resolution);
    set_mc_data(io_info, mc_data, mcube_time);

    // Note: All flags in mc_data should be set before calling 
    //       read isosurface lookup tables.
    // read isosurface lookup tables
    read_poly_isotable(io_info, mc_data, io_time);

    // Report number of cubes.
    report_num_cubes(full_scalar_grid, io_info, mc_data);

    // construct isosurface or interval volume
    if (mc_data.FlagIntervalVolume()) {
      construct_interval_volume(io_info, mc_data, mcube_time, io_time);
    }
    else {
      construct_isosurface(io_info, mc_data, mcube_time, io_time);
    }

    if (io_info.flag_report_time) {

      time_t end_time;
      time(&end_time);
      double total_elapsed_time = difftime(end_time, start_time);

      cout << endl;
      report_time(io_info, io_time, mcube_time, total_elapsed_time);
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

void construct_isosurface
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
  const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

  io_time.write_time = 0;
  for (unsigned int i = 0; i < io_info.isovalue.size(); i++) {

    const SCALAR_TYPE isovalue = io_info.isovalue[i];

    MC_ISOSURFACE mc_isosurface;
    MCUBE_INFO mcube_info(dimension);
    mcube_info.grid.num_cubes = num_cubes;
    SNAP_INFO snap_info(dimension);
    snap_info.grid.num_cubes = num_cubes;

    if (mc_data.Snap()) {

      if (io_info.UseCubeList()) {
        std::vector<VERTEX_INDEX> cube_list;

        float preprocessing_time;
        IJKSNAPMC::get_nonempty_snap_cubes
          (mc_data.ScalarGrid(), mc_data.isotable.cube_nep, 
           isovalue, cube_list, preprocessing_time);

        mcube_time.preprocessing += preprocessing_time;
        mcube_time.total += preprocessing_time;
      
        snapMC(mc_data, isovalue, cube_list, mc_isosurface, snap_info);
      }
      else {
        snapMC(mc_data, isovalue, mc_isosurface, snap_info);
      }

      mcube_time.Add(snap_info.time);
    }
    else {
      marching_cubes(mc_data, isovalue, mc_isosurface, mcube_info);
      mcube_time.Add(mcube_info.time);
    }


    OUTPUT_INFO output_info;
    if (mc_data.UseNEP()) {
      set_output_info(mc_data.isotable.cube_nep, io_info, i, output_info);
    }
    else {
      set_output_info(mc_data.isotable.cube, io_info, i, output_info);
    }

    VERTEX_INDEX nums = 
      mc_isosurface.simplex_vert.size()/numv_per_simplex;

    int subsample_resolution = 1;
    int supersample_resolution = 1;
    if (io_info.flag_subsample) 
      { subsample_resolution = io_info.subsample_resolution; }
    if (io_info.flag_supersample) 
      { supersample_resolution = io_info.supersample_resolution; }

    rescale_vertex_coord(subsample_resolution, supersample_resolution,
                         io_info.grid_spacing, mc_isosurface.vertex_coord);

    // *** SHOULD REALLY JUST USE max x coordinate. ***
    const COORD_TYPE max_coord =
      get_max_abs_array_value(mc_isosurface.vertex_coord);
    output_info.fig_param.coord_scale_factor =
      output_info.fig_param.DefaultPaperDrawRegionWidth()/max_coord;

    ELAPSED_TIME wall_time;

    if (mc_data.Snap()) {
      output_snapmc_isosurface
        (output_info, mc_data, mc_isosurface, snap_info);
    }
    else if (mc_data.UseNEP()) {
      output_nep_isosurface
        (output_info, mc_data, mc_isosurface, mcube_info);
    }
    else if (io_info.flag_color_alternating &&
             mc_isosurface.cube_containing_simplex.size() == nums) {
      output_isosurface_color_alternating
        (output_info, mc_data, mc_isosurface, mcube_info);
    }
    else {
      output_isosurface
        (output_info, mc_data, mc_isosurface, mcube_info);
    }

    io_time.write_time += wall_time.getElapsed();
  }
}

void construct_interval_volume
(const IO_INFO & io_info, const MC_DATA & mc_data,
 MCUBE_TIME & mcube_time, IO_TIME & io_time)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
  const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

  io_time.write_time = 0;
  for (unsigned int i = 0; i+1 < io_info.isovalue.size(); i++) {

    MC_ISOSURFACE mc_ivol;
    MCUBE_INFO mcube_info(dimension);
    mcube_info.grid.num_cubes = num_cubes;
    SNAP_INFO snap_info(dimension);
    snap_info.grid.num_cubes = num_cubes;

    MCVol(mc_data, io_info.isovalue[i], io_info.isovalue[i+1], 
          mc_ivol, mcube_info);

    mcube_time.Add(mcube_info.time);

    OUTPUT_INFO output_info;
    set_output_info(mc_data.isotable.cube, io_info, i, output_info);

    VERTEX_INDEX nums = 
      mc_ivol.simplex_vert.size()/numv_per_simplex;

    int subsample_resolution = 1;
    int supersample_resolution = 1;
    if (io_info.flag_subsample) 
      { subsample_resolution = io_info.subsample_resolution; }
    if (io_info.flag_supersample) 
      { supersample_resolution = io_info.supersample_resolution; }

    rescale_vertex_coord
      (subsample_resolution, supersample_resolution, io_info.grid_spacing,
       mc_ivol.vertex_coord);

    ELAPSED_TIME wall_time;

    output_isosurface(output_info, mc_data, mc_ivol, mcube_info);

    io_time.write_time += wall_time.getElapsed();
  }

}

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}

