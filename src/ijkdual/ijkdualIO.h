/*!
 *  @file ijkdualIO.h
 *  @brief IO classes and routines for ijkdual.
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
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJKDUALIO_H_
#define _IJKDUALIO_H__

#include <ctime>
#include <string>

#include "ijk.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkIO.tpp"
#include "ijkisoIO.tpp"
#include "ijkstring.tpp"

#include "ijkdualIO.tpp"

#include "ijkdual_types.h"
#include "ijkdual_datastruct.h"


namespace IJKDUAL {

  // **************************************************
  //! @name TYPE DEFINITIONS
  // **************************************************

  //@{

  typedef float COLOR_TYPE;           //!< Color type.

  //! Nrrd header.
  typedef typename IJK::NRRD_DATA<int, IJKDUAL::AXIS_SIZE_TYPE> NRRD_HEADER;

  typedef typename IJK::OUTPUT_DEFINITIONS::OUTPUT_FORMAT OUTPUT_FORMAT;

  /// Command line options.
  typedef enum
     {SUBSAMPLE_OPT, SUPERSAMPLE_OPT, POSITION_OPT, CUBE_CENTER_OPT,
      MANIFOLD_OPT, MULTI_ISOV_OPT, SINGLE_ISOV_OPT,
      TRIMESH_OPT, TRIMESH_UNIFORM_OPT,
      TRIMESH_MAX_MIN_ANGLE_OPT,
      TRIMESH_TRI4_MAX_MIN_ANGLE_OPT,
      ENVELOPE_TRI2_OR_TRI4_OPT, ENVELOPE_ONLY_TRI2_OPT,
      ENVELOPE_OFFSET_OPT,
      ENVELOPE_USE_RELATIVE_POSITIONS_OPT,
      SEPARATE_ISOV_OPT,

      // More options
      SEP_NEG_OPT, SEP_POS_OPT,
      SELECT_SPLIT_OPT, CONNECT_AMBIG_OPT,
      TRIMESH_SPLIT_MAX_OPT,
      TRIMESH_TRI4_ALL_QUADS_OPT,
      TRI4_IN_ENVELOPE_OPT, TRI4_ON_GRID_EDGE_OPT, TRI4_CENTROID_OPT, 
      MIN_DISTANCE_TO_CUBE_FACE_OPT,
      POSITION_OFFSET_OPT,
      NO_POSITION_OFFSET_OPT,
      SEP_ISOV_NEAR_CUBE_FACETS_OPT,
      MOVE_ALL_ISOV_AWAY_FROM_CUBE_FACETS_OPT,
      
      SEPARATE_ISOV_BY_CUBE_CENTER_OPT,
      ENVELOPE_TRI2_OR_TRI4_PREFER_TRI2_OPT,
      ENVELOPE_ONLY_TRI4_OPT,
      QEI_INTERPOLATE_SCALAR_OPT, QEI_BILINEAR_OPT,
      QEI_AVERAGE_OPT,
      INTERIOR_VERTEX_RATIO_OPT,
      NO_COMMENTS_OPT,

      // Testing options
      POSITION_TESTING_OPT, RANDOM_SEED_OPT,
      RANDOM_POS_OPT, RANDOM_POS_UNIFORM_OPT,
      RANDOM_POS_DISTRIBUTION_OPT,      
      RANDOM_SAMPLE_ENTIRE_CUBE_OPT,
      RANDOM_SAMPLE_CUBE_BOUNDARY_OPT,
      RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT,
      SPLIT_NON_MANIFOLD_AMBIG_OPT,
      FORCE_ADD_ALL_INTERIOR_VERTICES_OPT,

      // Usage or I/O options
      USAGE_OPT, MORE_OPTIONS_OPT, TESTING_OPTIONS_OPT,
            ALL_OPTIONS_OPT,
      HELP_OPT, HELP_MORE_OPT, HELP_TESTING_OPT, HELP_ALL_OPT,
      OFF_OPT, PLY_OPT, FIG_OPT,
      OUTPUT_FILENAME_OPT, OUTPUT_FILENAME_PREFIX_OPT, STDOUT_OPT,
      LABEL_WITH_ISOVALUE_OPT,
      NO_WRITE_OPT, SILENT_OPT, TERSE_OPT, NO_WARN_OPT,
      INFO_OPT, TIME_OPT, DETAILED_TIME_OPT, OUT_ISOV_OPT, UNKNOWN_OPT}
    COMMAND_LINE_OPTION_TYPE;

   /// Command line option group types.
   typedef enum {
     REGULAR_OPTG, EXTENDED_OPTG, TESTING_OPTG
   } COMMAND_LINE_OPTION_GROUP;

  //@}


  // **************************************************
  //! @name IO INFORMATION
  // **************************************************

  //@{

  /// IO information.
  class IO_INFO:
    public IJK::IO_PARAM_BASE
    <SCALAR_TYPE,COORD_TYPE,DUALISO_DATA_PARAM> {

  public:
    IO_INFO() {};
    ~IO_INFO() {};

  };


  /// Output information.
  class OUTPUT_INFO:public IJK::OUTPUT_PARAM_BASE<IO_INFO> {

  public:

    /// @brief Set \a dimension and \a num_vertices_per_isopoly.
    /// - \a num_vertices_per_isopoly depends upon \a dimension and \a mesh_type.
    void SetDimension(const int d);
  };

  //@}


  // **************************************************
  //! @name TIMING FUNCTIONS/CLASSES
  // **************************************************

  ///@{

  // *** SHOULD BE MOVED TO COMMON .h FILE ***
  /// Elapsed CPU time.
  class ELAPSED_CPU_TIME {

  protected:
    clock_t t;

  public:
    ELAPSED_CPU_TIME() { t = clock(); };

    clock_t getElapsed() {
      clock_t old_t = t;
      t = clock();
      return(t - old_t);
    };
  };

  /// Elapsed wall time.
  class ELAPSED_TIME {

  protected:
    time_t t;

  public:
    ELAPSED_TIME() { time(&t);  };

    double getElapsed() {
      time_t old_t = t;
      time(&t);
      return(difftime(t,old_t));
    };
  };

  /// IO time.
  struct IO_TIME {
    double read_nrrd_time;  ///< Wall time to read nrrd file.
    double write_time;      ///< Wall time to write output.
  };

  ///@}


  // **************************************************
  //! @name PARSE COMMAND LINE
  // **************************************************

  ///@{

  /// Parse the command line.
  void parse_command_line(int argc, char **argv, IO_INFO & io_info);

  /*!
   *  @brief Check input information in io_info and reset, if necessary.
   *  - Exit if usage error found.
   *  - Print warnings when appropriate.
   *  - Reset input parameters, if necessary.
   */
  void check_and_reset_input
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   IO_INFO & io_info);

  ///@}


  // **************************************************
  //! @name RESCALE ROUTINES
  // **************************************************

  ///@{

  /// Rescale vertex coordinates by grid_spacing.
  void rescale_vertex_coord
    (const int dimension, const COORD_TYPE * grid_spacing,
     std::vector<COORD_TYPE> & vertex_coord);

  /// @brief Rescale vertex coordinates by grid_spacing.
  /// @pre grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord(const std::vector<COORD_TYPE> & grid_spacing,
                            std::vector<COORD_TYPE> & vertex_coord);

  /// Rescale subsampled/supersampled vertex coordinates.
  void rescale_vertex_coord
  (const int subsample_resolution, const int supersample_resolution, COORD_ARRAY & vertex_coord);

  /// @brief Rescale vertex coordinates by subsample_resolution and supersample_resolution
  ///        and by grid_spacing.
  /// @pre grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord
    (const int subsample_resolution, const int supersample_resolution,
     const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord);

  ///@}


  // **************************************************
  //! @name WRITE_DUAL_MESH
  // **************************************************

  /*!
   *  @brief Write dual mesh to output file(s).
   *  - Output file(s) and format(s) depend on flags in \a output_info.
   *  - May write to more than one file.
   */
  void write_dual_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist);

  /*!
   *  @brief Write dual mesh and record output time.
   *  - Output file(s) and format(s) depend on flags in \a output_info.
   */
  void write_dual_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     IO_TIME & io_time);


  /*!
   *  @brief Write dual isosurface mesh of quad and triangles.
   *  @param output_info Output information.
   *  @param output_format Output format.
   *  @param vertex_coord List of vertex coordinates.
   */
  void write_dual_quad_tri_mesh
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert);

  /*!
   *  @brief Write dual isosurface mesh of quad and triangles.
   *  @param output_info Output information.
   *  @param vertex_coord List of vertex coordinates.
   */
  void write_dual_quad_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert);

  /// Write dual isosurface mesh of quad and triangles and record output time.
  void write_dual_quad_tri_mesh
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   IO_TIME & io_time);

  /*!
   *  @brief Write dual isosurface triangular mesh and color vertices
   *    with output format output_format.
   *  @param output_info Output information.
   *  @param output_format Output format.
   *  @param vertex_coord List of vertex coordinates.
   *  @param tri_vert[] List of triangle vertices.
   *         tri_vert[3*i+k] is k'th vertex of triangle i.
   */
  void write_dual_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /*!
   *  @brief Write dual isosurface triangular mesh, color vertices.
   *  @param output_info Output information.
   *  @param vertex_coord List of vertex coordinates.
   *  @param tri_vert[] List of triangle vertices.
   *         tri_vert[3*i+k] is k'th vertex of triangle i.
   */
  void write_dual_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// Write dual isosurface triangular mesh, color vertices and record output time.
  void write_dual_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
   IO_TIME & io_time);

  /*!
   *  @brief Write dual isosurface mesh of quad and triangles and color vertices with output format output_format.
   *  @param output_info Output information.
   *  @param output_format Output format.
   *  @param vertex_coord List of vertex coordinates.
   */
  void write_dual_quad_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /*!
   *  @brief Write dual isosurface mesh of quad and triangles and color vertices.
   *  @param output_info Output information.
   *  @param vertex_coord List of vertex coordinates.
   */
  void write_dual_quad_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  /// @brief Write dual isosurface mesh of quad and triangles, color vertices, and report output time.
  void write_dual_quad_tri_mesh_color_vertices
  (const OUTPUT_INFO & output_info,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
   IO_TIME & io_time);

  ///@}


  // **************************************************
  //! @name SET ROUTINES
  // **************************************************

  ///@{

  /*!
   *  @brief Set dualiso_data based on io_info.
   *  @pre Scalar field in dualiso_data must be set before
   *    this routines is called.
   */
  void set_dualiso_data
  (const IO_INFO & io_info, DUALISO_DATA & dualiso_data, 
   DUALISO_TIME & dualiso_time);

  /*!
   *  @brief Set output information, copying information from io_info.
   *  @param dimension Dimension of volume containing the isosurface.
   *  @param i Index of isovalue being currently processed.
   */
  void set_output_info
  (const IO_INFO & io_info, const int dimension, const int i,
   OUTPUT_INFO & output_info);

  ///@}


  // **************************************************
  //! @name CREATE STRING FROM enum TYPE
  // **************************************************

  ///@{

  /// @brief Convert enum type VERTEX_POSITION_METHOD to string.
  void vpos2string(const VERTEX_POSITION_METHOD vpos,
                   std::string & vpos_string);

  /// @brief Convert enum type QUAD_TRI_METHOD to string.
  void tri_method2string
    (const QUAD_TRI_METHOD quad_tri_method,
     std::string & method_string);

  /// @brief Convert enum type INTERIOR_VERTEX_POSITION_METHOD to string.
  void splitv_pos2string
  (const IJK::INTERIOR_VERTEX_POSITION_METHOD position_method,
   std::string & splitv_pos_string);

  /// @brief Convert enum type POLY_EDGE_INTERSECTION_METHOD to string.
  void qei2string
  (const IJK::POLY_EDGE_INTERSECTION_METHOD qei_method,
   std::string & qei_string);
  
  ///@}


  // **************************************************
  //! @name CREATE MESH FILE COMMENTS
  // **************************************************

  ///@{

  /// @brief Add isosurface output comments to comment.
  void add_meshfile_comments
  (const OUTPUT_INFO & output_info,
   std::vector<std::string> & comment);

  ///@}


  // **************************************************
  //! @name REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // **************************************************

  ///@{

  /// Report number of cubes.
  void report_num_cubes
    (const DUALISO_GRID & full_grid, const IO_INFO & io_info, 
     const DUALISO_DATA & dualiso_data);

  void report_num_cubes
    (const DUALISO_GRID & full_grid, const IO_INFO & io_info, 
     const DUALISO_GRID & dualiso_data_grid);

  /// Print warning that isosurface mesh may not be a manifold.
  void warn_non_manifold(const IO_INFO & io_info);

  ///@}


  // **************************************************
  //! @name REPORT TIMING INFORMATION
  // **************************************************

  ///@{

  void report_dualiso_time
    (const IO_INFO & io_info, const DUALISO_TIME & dualiso_time, 
     const char * mesh_type_string);

  void report_time
    (const IO_INFO & io_info, const IO_TIME & io_time, 
     const DUALISO_TIME & dualiso_time, const double total_elapsed_time);

  ///@}


  // **************************************************
  //! @name USAGE/HELP MESSAGES
  // **************************************************

  ///@{

  /// Print usage message and exit with non-zero exit code.
  void usage_error();

  /// Print usage message and exit.
  void usage(std::ostream & out, const int return_code);

  /// Print option and exit.
  void print_options
    (std::ostream & out, const COMMAND_LINE_OPTION_GROUP group,
     const int return_code);

  /// Print usage message with all options and exit.
  void usage_all(std::ostream & out, const int return_code);

  /// Print help message and exit.
  void help();

  /// Print help messages for options in group.
  void print_options_help(const COMMAND_LINE_OPTION_GROUP group);

  /// Print help message for all options and exit.
  void help_all();

  ///@}

}

#endif
