/*!
 *  @file ijkmcubeIO.h
 *  @brief IO classes and routines for ijkmcube.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009-2023 Rephael Wenger

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

#ifndef _IJKMCUBEIO_
#define _IJKMCUBEIO_

#include <ctime>
#include <string>

#include "ijk.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkIO.tpp"
#include "ijkisoIO.tpp"
#include "ijkstring.tpp"

#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"
#include "ijksnapmc.h"
#include "ijkNrrd.h"


namespace IJKMCUBE {

// **************************************************
// TYPE DEFINITIONS
// **************************************************

  //! Nrrd header.
  typedef typename IJK::NRRD_DATA<int, IJKMCUBE::AXIS_SIZE_TYPE> NRRD_HEADER;

  typedef float COLOR_TYPE;           /// Color type.

  using IJKSNAPMC::SNAP_INFO;

  /// Commmand line options.
  typedef enum {
    REGION_OPT, OCTREE_OPT, USE_CUBE_LIST_OPT, NEP_OPT, NEP_DUP_OPT,
    IVOL_OPT, SNAP_OPT, HIGHRES_OPT,
    TOPOLOGY_OPT, INTERPOLATE_OPT, NUMRES_LEVELS_OPT,
    SUBSAMPLE_OPT, SUPERSAMPLE_OPT,
    LOG_INTERVAL_OPT,
    COLOR_ALTERNATING_OPT,

    // Options specifying isosurface lookup table properties.
    EDGE_GROUPS_OPT, CHULL_OPT, 
    SEP_POS_OPT, SEP_NEG_OPT, SEP_OPP_OPT, NO_SEP_OPP_OPT,
    POS_ORIENT_OPT, NEG_ORIENT_OPT,

    USAGE_OPT, MORE_OPTIONS_OPT, ALL_OPTIONS_OPT,
    HELP_OPT, HELP_MORE_OPT, HELP_ALL_OPT,
    OFF_OPT, PLY_OPT, FIG_OPT,
    ISOTABLE_DIR_OPT, ISOTABLE_PATH_OPT,
    OUTPUT_FILENAME_OPT, OUTPUT_FILENAME_PREFIX_OPT, STDOUT_OPT,
    LABEL_WITH_ISOVALUE_OPT,
    NO_WRITE_OPT, SILENT_OPT, VERBOSE_OPT, TIME_OPT, UNKNOWN_OPT,

    // More (extended) options
    EDGE1_OPT, EDGE2_OPT,
    NO_COMMENTS_OPT
  }
  COMMAND_LINE_OPTION_TYPE;

  typedef enum {
    REGULAR_OPTG, EXTENDED_OPTG, TESTING_OPTG
  } COMMAND_LINE_OPTION_GROUP;


// **************************************************
// IO INFORMATION
// **************************************************

  ///@{

  /// IO Information
  class IO_INFO:
    public IJK::IO_PARAM_BASE
      <SCALAR_TYPE,COORD_TYPE,MC_DATA_FLAGS> {

    protected:

      /// Initialize data structure.
      void Init();


    public:

    std::string isotable_directory;
    std::string isotable_path;

    bool flag_color_alternating;

    // Isosurface lookup table properties.
    IJKMCUBE_TABLE::ISOSURFACE_TABLE_PROPERTIES isotable_properties;

    public:
      IO_INFO() { Init(); };
      ~IO_INFO() {};
  };


  /// Output information.
  class OUTPUT_INFO:public IJK::OUTPUT_PARAM_BASE<IO_INFO> {

  protected:
    void Init() {};

  public:

    /// Return number of vertices in each simplex.
    int NumVerticesPerSimplex() const
    { return num_vertices_per_isopoly; }

    OUTPUT_INFO() { Init(); };
    ~OUTPUT_INFO() {};
  };


// **************************************************
// TIMING FUNCTIONS/CLASSES
// **************************************************

  // *** SHOULD BE MOVED TO COMMON .h FILE ***/
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
    double read_table_time; ///< Wall time to read isosurface lookup table.
    double read_nrrd_time;  ///< Wall time to read nrrd file.
    double write_time;      ///< Wall time to write output.
  };

// **************************************************
//! @name PARSE COMMAND LINE
// **************************************************

  //@{

  /// Parse the command line.
  void parse_command_line(int argc, char **argv, IO_INFO & io_info);

  /*!
   *  @brief Check input information in io_info.
   *  - Exit if usage error found.
   *  - Print warnings when appropriate.
   */
  void check_input
    (const IO_INFO & io_info, 
     const MC_SCALAR_GRID_BASE & scalar_grid);

  //@}


// **************************************************
//! @name READ ISOSURFACE LOOKUP TABLE(S)
// **************************************************

  //@{

  /// @brief Read cube isosurface lookup table.
  void read_cube_isotable
    (const int dimension, const IO_INFO & io_info, 
     ISOSURFACE_TABLE & cube_isotable, IO_TIME & io_time);

  /// @brief Read cube, pyramid and simplex isosurface lookup tables.
  /// Only reads pyramid and simplex tables if needed.
  void read_poly_isotable
    (const IO_INFO & io_info, MC_DATA & mc_data, IO_TIME & io_time);

  /// @brief Determine isosurface table type 
  ///   from ISOSURFACE_TABLE_PROPERTIES.
  ISOTABLE_TYPE get_isotable_type
    (const IJKMCUBE_TABLE::ISOSURFACE_TABLE_PROPERTIES & isotable_properties);

  /// Get default isotable directory and get directory from environment.
  void get_isotable_directory(std::string & isotable_directory);

  /// @brief Read isosurface lookup table.
  /// - Construct lookup table file name from dimension, poly_name 
  ///   and isotable_type.
  void read_isosurface_table
    (const int dimension, const char * poly_name,
     const IO_INFO & io_info,
     ISOSURFACE_TABLE & isotable, IO_TIME & io_time);

  /// @overload
  /// @brief Read isosurface lookup table. (No io_time.)
  /// - Version that does not return io_time.
  void read_isosurface_table
    (const int dimension, const char * poly_name,
     const IO_INFO & io_info, ISOSURFACE_TABLE & isotable);

  /*!
   *  @brief Open isosurface lookup table.
   *  - Constructs filename from dimension, poly_name and isotable_type.
   *  @param[out] isotable_filename File name of isosurface lookup table.
   *  @param[out] isotable_pathname Path name of isosurface lookup table.
   *  @param[out] isotable_file Isosurface table file stream.
   *  @param[out] path_list List of paths searched for lookup table.
   */
  void open_isotable_file
    (const int dimension, const char * poly_name,
     const IO_INFO & io_info, const bool flag_use_local_directory,
     std::string & isotable_filename,
     std::string & isotable_pathname,
     std::ifstream & isotable_file,
     std::vector<std::string> & path_list);

  //@}


// **************************************************
//! @name OUTPUT ISOSURFACE
// **************************************************

  //@{

  void output_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const MCUBE_INFO & mcube_info);

  void output_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const MCUBE_INFO & mcube_info);

  void output_nep_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const MCUBE_INFO & mcube_info);

  void output_nep_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const MCUBE_INFO & mcube_info);

  void output_snapmc_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const SNAP_INFO & snap_info);

  void output_snapmc_isosurface
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const SNAP_INFO & snap_info);

  void output_isosurface_color
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const MCUBE_INFO & mcube_info);

  void output_isosurface_color
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
     const MCUBE_INFO & mcube_info);

  void output_isosurface_color_alternating
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const MC_ISOSURFACE & mc_isosurface,
     const MCUBE_INFO & mcube_info);

  //@}


// **************************************************
//! @name RESCALE ROUTINES
// **************************************************

  //@{

  /// Rescale subsampled/supersampled vertex coordinates.
  /// Also rescale to reflect grid spacing.
  void rescale_vertex_coord
    (const OUTPUT_INFO & output_info, std::vector<COORD_TYPE> & vertex_coord);

  /// Rescale vertex coordinates by grow and shrink factor and by grid_spacing.
  /// Precondition: grid_spacing.size() equals vertex dimension.
  void rescale_vertex_coord
    (const int subsample_resolution, const int supersample_resolution,
     const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord);

  //@}


// **************************************************
//! @name WRITE_MESH
// **************************************************

  //@{

  void write_mcube_mesh
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist);

  void write_mcube_mesh_color
    (const OUTPUT_INFO & output_info,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist,
     const COLOR_TYPE * front_color, const COLOR_TYPE * back_color);

  //@}


// **************************************************
//! @name SET ROUTINES
// **************************************************

  //@{

  /*!
   *  @brief Set mc_data based on io_info.
   *  @pre Scalar field in mc_data must be set before
   *    this routines is called.
   */
  void set_mc_data
    (const IO_INFO & io_info, MC_DATA & mc_data, MCUBE_TIME & mcube_time);

  /// Set output_info based on isotable, io_info and isovalue index i.
  void set_output_info
    (const ISOSURFACE_TABLE & isotable, const IO_INFO & io_info, 
     const int i, OUTPUT_INFO & output_info);

  /// Set simplices in alternating cubes to have different colors.
  void set_color_alternating
    (const MC_GRID & grid, const std::vector<VERTEX_INDEX> & cube_list, 
     COLOR_TYPE * color);

  //@}


  // **************************************************
  //! @name CREATE STRING FROM enum TYPE
  // **************************************************

  //@{

  /// Convert enum type ISOSURFACE_TOPOLOGY to string.
  void topology2string
    (const ISOSURFACE_TOPOLOGY & isosurface_topology,
     std::string & topology_string);

  //@}


  // **************************************************
  //! @name CREATE MESH FILE COMMENTS
  // **************************************************

  //@{

  void add_meshfile_comments
    (const OUTPUT_INFO & output_info,
     std::vector<std::string> & comment);

  //@}


  // **************************************************
  //! @name REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // **************************************************

  //@{

  void report_num_cubes
    (const MC_GRID & full_grid, const IO_INFO & io_info, 
     const MC_DATA & mc_data);

  void report_iso_info
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const MCUBE_INFO & mcube_info);

  void report_nep_info
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const MCUBE_INFO & mcube_info);

  void report_snap_info
    (const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
     const std::vector<COORD_TYPE> & vertex_coord, 
     const std::vector<VERTEX_INDEX> & slist, 
     const SNAP_INFO & snap_info);

  //@}


// **************************************************
//! @name REPORT TIMING INFORMATION
// **************************************************

  //@{

  void report_mcube_time
    (const IO_INFO & io_info, const MCUBE_TIME & mcube_time, 
     const char * mesh_type_string);

  void report_time
    (const IO_INFO & io_info, const IO_TIME & io_time, 
     const MCUBE_TIME & mcube_time, const double total_elapsed_time);

  //@}


// **************************************************
//! @name USAGE/HELP MESSAGES
// **************************************************

  //@{

  /// Print usage message and exit.
  void usage(std::ostream & out, const int return_code);

  /// Print usage message and exit with non-zero exit code.
  void usage_error();

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

  //@}

}

#endif
