/*!
 *  @file ijkisoIO.tpp
 *  @brief Templates for isosurface IO routines.
 *  - Class IO_PARAM.
 *  - Output formats: Geomview .off, OpenInventor .iv (3D), Fig .fig (2D)
 *    Stanford .ply, and Visualization Toolkit, .vtk.
 *  - Requires compilation with c++17.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2023 Rephael Wenger

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

#ifndef _IJKIO_PARAM_
#define _IJKIO_PARAM_

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkIO.tpp"
#include "ijkstring.tpp"


namespace IJK {

  // *****************************************************************
  // Read Nrrd file into SCALAR_GRID data structure
  // *****************************************************************

  /*!
   *  @brief Read Nrrd file into scalar_grid.
   *  @tparam SCALAR_GRID_TYPE Scalar grid type.
   *    - Must include SetSpacing() member function.
   */
  template <typename SCALAR_GRID_TYPE, typename DATA_TYPE,
            typename AXIS_SIZE_TYPE>
  void read_nrrd_file
  (const char * input_filename, SCALAR_GRID_TYPE & scalar_grid,
   IJK::NRRD_DATA<DATA_TYPE,AXIS_SIZE_TYPE> & nrrd_header)
  {
    typedef typename SCALAR_GRID_TYPE::SPACING_TYPE SPACING_TYPE;

    IJK::GRID_NRRD_IN<DATA_TYPE,AXIS_SIZE_TYPE> nrrd_in;
    IJK::PROCEDURE_ERROR error("read_nrrd_file");

    nrrd_in.ReadScalarGrid(input_filename, scalar_grid, nrrd_header, error);
    if (nrrd_in.ReadFailed()) { throw error; }

    std::vector<SPACING_TYPE> grid_spacing;
    nrrd_header.GetSpacing(grid_spacing);
    scalar_grid.SetSpacing(grid_spacing);
  }

  
  /*!
   *  @brief Read Nrrd file into scalar_grid.
   *  - Version using std::string for input_filename.
   */
  template <typename SCALAR_GRID_TYPE, typename DATA_TYPE,
            typename AXIS_SIZE_TYPE>
  void read_nrrd_file
  (const std::string & input_filename, SCALAR_GRID_TYPE & scalar_grid,
   IJK::NRRD_DATA<DATA_TYPE,AXIS_SIZE_TYPE> & nrrd_header)
  {
    read_nrrd_file(input_filename.c_str(), scalar_grid, nrrd_header);
  }


  // *****************************************************************
  // Class IO_PARAM
  // *****************************************************************
  
  /// IO Parameters
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename ISO_DATA_PARAM>
  class IO_PARAM_BASE:public ISO_DATA_PARAM {

  public:

    typedef _SCALAR_TYPE SCALAR_TYPE;
    typedef _COORD_TYPE COORD_TYPE;
    typedef typename OUTPUT_DEFINITIONS::OUTPUT_FORMAT OUTPUT_FORMAT;


  protected:

    OUTPUT_DEFINITIONS output_definitions;

    /// Initialize data structure.
    void Init();


  public:

    // Input parameters.
    std::vector<SCALAR_TYPE> isovalue;     ///< List of isovalues.
    std::vector<std::string> isovalue_string;
    std::vector<COORD_TYPE> grid_spacing;
    int subsample_resolution;
    int supersample_resolution;
    int region_length;

    /// @brief List of high resolution arguments.
    /// - Currently ignored with ijkdual.
    std::vector<std::string> high_resolution_string;

    // File/directory names.
    std::string input_filename;
    std::string output_filename;
    std::string output_filename_prefix;
    std::string output_off_filename;
    std::string output_ply_filename;
    std::string output_fig_filename;
    std::string output_iv_filename;
    std::string report_isov_filename;

    // Flags.
    bool flag_output_off;    ///< Output Geomview .off file.
    bool flag_output_ply;    ///< Output PLY file.
    bool flag_output_fig;    ///< Output FIG file.
    bool flag_report_time;
    bool flag_report_time_detailed;
    bool flag_report_info;
    bool flag_use_stdout;
    bool flag_nowrite;
    bool flag_silent;
    bool flag_terse;
    bool flag_verbose;
    bool flag_no_warn;
    bool flag_subsample;
    bool flag_report_all_isov;
    bool flag_supersample;
    bool flag_no_comments;

    /// If true, output filename includes isovalue.
    bool flag_label_with_isovalue;

    bool are_output_filenames_set;

    /// If true, triangulate mesh.
    bool flag_trimesh;

    /*!
     *  @brief If true, color isosurface boundary vertices. 
     *  - Currently not used (always false).
     */
    bool flag_color_boundary_vert;

    // Flags for options set.

    /// Is quadrilateral edge intersection method set?
    bool is_qei_method_set;

    bool is_tri4_position_method_set;
    bool is_file_format_set;
    bool is_random_seed_set;

    /// Return true if output_filename is not "".
    bool IsOutputFilenameSet() const
    { return (output_filename != ""); }

    /// Return number of selected output formats.
    int NumOutputFormats() const;

    /*!
     *  @brief Return format of filename suffix.
     *  - Return default_format is filename suffix does not correspond
     *    to any file types in output_format_suffix.
     */
    OUTPUT_FORMAT GetFileType
    (const std::string & filename, const OUTPUT_FORMAT & default_format)
    { return output_definitions.GetFileType(filename, default_format); }

    /// @brief Return suffix for given output format.
    const char * const GetOutputFormatSuffix
      (const OUTPUT_FORMAT output_format) const
    { return output_definitions.GetOutputFormatSuffix(output_format); }


    /// Copy IO_PARAM_BASE
    void Set(const IO_PARAM_BASE & io_param_base);

    /// @brief Set output format to output_format.
    /// - Note: More than one output format may be set.
    void SetOutputFormat(const OUTPUT_FORMAT output_format);

    /// @brief Set output filename for output format that is flagged true.
    /// @pre At most one output format is flagged true.
    void SetOutputFilename(const char * output_filename);

    /// @brief Set output filename for output format that is flagged true.
    /// @pre At most one output format is flagged true.
    void SetOutputFilename(const std::string & output_filename)
    { SetOutputFilename(output_filename.c_str()); }

    /// Set output filename for a given output format.
    void SetOutputFilename
    (const OUTPUT_FORMAT output_format, const char * output_filename);

    /// Set output filename for a given output format.
    void SetOutputFilename
    (const OUTPUT_FORMAT output_format, const std::string & output_filename)
    { SetOutputFilename(output_format, output_filename.c_str()); }

    /// @brief Construct output filenames.
    /// @param i Construct filenames for isovalue i.
    void ConstructOutputFilenames
      (const int i, const bool flag_interval_volume = false);


  public:
    IO_PARAM_BASE() { Init(); };
    ~IO_PARAM_BASE() { Init(); };


    // Print routines (mainly for debugging.)

    /// Print IO_PARAM_BASE information. (May not be complete.)
    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

  };


  // *****************************************************************
  // Class OUTPUT_PARAM_BASE
  // *****************************************************************

  /*!
   *  @brief Output parameters.
   *  @tparam _IO_PARAM IO parameters class.
   *    - Usually derived from IO_PARAM_BASE.
   *    - Must include type definition of SCALAR_TYPE.
   */
  template <typename _IO_PARAM>
  class OUTPUT_PARAM_BASE:public _IO_PARAM {

  public:
    typedef typename _IO_PARAM::SCALAR_TYPE SCALAR_TYPE;
    typedef typename _IO_PARAM::COORD_TYPE COORD_TYPE;
    typedef int FIG_UNITS_TYPE;


  protected:
    void Init();

  public:
    /// Dimension of scalar grid.
    int dimension;

    /// Dimension of mesh.
    int mesh_dimension;

    /// Number of vertices per output isosurface polytope.
    int num_vertices_per_isopoly;

    /// Isovalues in this output.
    std::vector<SCALAR_TYPE> output_isovalue;

    /// Fig (xfig) output parameters.
    IJK::FIG_OUTPUT_PARAM<FIG_UNITS_TYPE> fig_param;

    /// Constructor.
    OUTPUT_PARAM_BASE() { Init(); };

    /// Destructor.
    ~OUTPUT_PARAM_BASE() { Init(); };

    /// Return dimension of scalar grid.
    int Dimension() const
    { return dimension; }

    // Add comment routines.

    /// Add isovalue comment to array comment[].
    void AddIsovalueComment(std::vector<std::string> & comment) const;

    /// Add high_resolution_region comment(s) to array comment[].
    void AddHighResolutionComments(std::vector<std::string> & comment) const;

    /// @brief Add subsample resolution comment to array comment[].
    /// - Only adds comment if flag_subsample is true.
    void AddSubsampleResolutionComment
      (std::vector<std::string> & comment) const;

    /// @brief Add supersample resolution comment to array comment[].
    /// - Only adds comment if flag_supersample is true.
    void AddSupersampleResolutionComment
      (std::vector<std::string> & comment) const;

    // Print routines (mainly for debugging.)

    /// Print IO_PARAM_BASE information. (May not be complete.)
    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;
  };


  // *****************************************************************
  // Class IO_PARAM_BASE member functions
  // *****************************************************************

  // Initialize IO_PARAM_BASE.
  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::Init()
  {
    // Initialize input parameters.
    isovalue.clear();
    isovalue_string.clear();
    input_filename .clear();
    subsample_resolution = 2;
    supersample_resolution = 2;
    region_length = 1;

    // Initialize IO flags.
    flag_output_off = false;
    flag_output_ply = false;
    flag_output_fig = false;
    flag_report_time = false;
    flag_report_time_detailed = false;
    flag_report_info = false;
    flag_use_stdout = false;
    flag_nowrite = false;
    flag_silent = false;
    flag_terse = false;
    flag_verbose = false;
    flag_no_warn = false;
    flag_no_comments = false;
    flag_supersample = false;
    flag_subsample = false;
    flag_report_all_isov = false;
    flag_label_with_isovalue = false;
    are_output_filenames_set = false;
    flag_trimesh = false;

    // color isosurface boundary vertices
    flag_color_boundary_vert = false;

    is_qei_method_set = false;
    is_tri4_position_method_set = false;
    is_file_format_set = false;
    is_random_seed_set = false;
  }


  // Return number of output formats.
  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  int IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  NumOutputFormats() const
  {
    int num_output_formats = 0;

    if (flag_output_off) { num_output_formats++; }
    if (flag_output_ply) { num_output_formats++; }
    if (flag_output_fig) { num_output_formats++; }

    return(num_output_formats);
  }

  
  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  Set(const IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM> & io_param)
  {
    *this = io_param;
  }


  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  SetOutputFormat(const OUTPUT_FORMAT output_format)
  {
    IJK::PROCEDURE_ERROR error("IO_PARAM_BASE::SetOutputFormat");

    switch(output_format) {

    case OUTPUT_DEFINITIONS::OFF:
      flag_output_off = true;
      break;

    case OUTPUT_DEFINITIONS::PLY:
      flag_output_ply = true;
      break;

    case OUTPUT_DEFINITIONS::FIG:
      flag_output_fig = true;
      break;

    default:
      error.AddMessage
        ("Programming error. Unable to set output format to ",
         GetOutputFormatSuffix(output_format), ".");
      throw error;
    }
  }


  // Set output filename for output format flagged true.
  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  SetOutputFilename(const char * output_filename)
  {
    IJK::PROCEDURE_ERROR error("IO_PARAM_BASE::SetOutputFilename");

    int num_output_formats = 0;

    if (flag_output_off) {
      output_off_filename = output_filename;
      num_output_formats++;
      are_output_filenames_set = true;
    }

    if (flag_output_ply) {
      output_ply_filename = output_filename;
      num_output_formats++;
      are_output_filenames_set = true;
    }

    if (flag_output_fig) {
      output_fig_filename = output_filename;
      num_output_formats++;
      are_output_filenames_set = true;
    }

    if (num_output_formats > 1) {
      error.AddMessage
        ("Programming error.  More than one output format is set.");
      error.AddMessage
        ("  Cannot set same output filename for more than one output format.");
      throw error;
    }
  }


  // Set output filename for a given output format.
  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  SetOutputFilename
  (const OUTPUT_FORMAT output_format, const char * output_filename)
  {
    IJK::PROCEDURE_ERROR error("IO_PARAM_BASE::SetOutputFilename");

    switch(output_format) {

    case OUTPUT_DEFINITIONS::OFF:
      output_off_filename = output_filename;
      break;

    case OUTPUT_DEFINITIONS::PLY:
      output_ply_filename = output_filename;
      break;

    case OUTPUT_DEFINITIONS::FIG:
      output_fig_filename = output_filename;
      break;

    default:
      error.AddMessage
        ("Programming error.  Unknown file type ",
         GetOutputFormatSuffix(output_format), ".");
      throw error;
    }

    are_output_filenames_set = true;
  }


  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  ConstructOutputFilenames(const int i, const bool flag_interval_volume)
  {
    typedef typename std::vector<std::string>::size_type SIZE_TYPE;

    using std::string;

    string prefix, suffix;
    string ofilename;

    // create output filename
    if (output_filename_prefix == "") {

      const string fname =
        std::filesystem::path(input_filename).filename().string();

      // construct output filename
      IJK::split_string(fname, '.', prefix, suffix);
      if (suffix == "nrrd" || suffix == "nhdr") { ofilename = prefix; }
      else { ofilename = string(input_filename); }
    }
    else {
      ofilename = output_filename_prefix;
    }

    if (flag_interval_volume) {
      if (flag_label_with_isovalue || isovalue_string.size() > 2) {
        if (SIZE_TYPE(i+1) < isovalue_string.size()) {
          ofilename += string(".") + string("ivol=") + isovalue_string[i] +
                       string("_") + isovalue_string[i+1];
        }
      }
    }
    else {
      if (flag_label_with_isovalue || isovalue_string.size() > 1) {
        ofilename += string(".") + string("isov=") + isovalue_string[i];
      }
    }

    output_off_filename = ofilename + ".off";
    output_ply_filename = ofilename + ".ply";
    output_fig_filename = ofilename + ".fig";
    output_iv_filename = ofilename + ".iv";
  }


  /// Print IO_INFO information. (May not be complete.)
  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix << "flag_output_off: "
        << IJK::bool2string(flag_output_off) << "\n";
    out << line_prefix << "flag_ouput_ply: "
        << IJK::bool2string(flag_output_ply) << "\n";

    // SHOULD PRINT LOTS MORE FIELDS.
  }


  // *****************************************************************
  // Class OUTPUT_PARAM member functions
  // *****************************************************************

  template <typename _IO_PARAM>
  void OUTPUT_PARAM_BASE<_IO_PARAM>::Init()
  {
    const int DIM2(2);
    const int DIM3(3);
    const int NUM_VERT_PER_TRIANGLE(3);

    dimension = DIM3;
    mesh_dimension = DIM2;
    num_vertices_per_isopoly = NUM_VERT_PER_TRIANGLE;
  }


  template <typename _IO_PARAM>
  void OUTPUT_PARAM_BASE<_IO_PARAM>::
  AddIsovalueComment(std::vector<std::string> & comment) const
  {
    if (this->output_isovalue.size() == 1) {
      std::string isovalue_string;
      IJK::create_string
        ("Isovalue: ", this->output_isovalue[0], isovalue_string);
      comment.push_back(isovalue_string);
    }
    else if (this->output_isovalue.size() > 1) {
      std::string isovalue0_string, isovalue1_string;
      std::string isovalue_string;
      IJK::val2string(this->output_isovalue[0], isovalue0_string);
      IJK::val2string(this->output_isovalue[1], isovalue1_string);
      isovalue_string = "Isovalues: " + isovalue0_string + " " + isovalue1_string;
      comment.push_back(isovalue_string);
    }
  }


  // Add high_resolution_region comment(s) to array comment[].
  template <typename _IO_PARAM>
  void OUTPUT_PARAM_BASE<_IO_PARAM>::
  AddHighResolutionComments(std::vector<std::string> & comment) const
  {
    for (unsigned int i = 0; i < this->high_resolution_string.size(); i++) {
      std::string str_i;
      IJK::val2string(i, str_i);
      std::string highres_string_i =
        std::string("HighResRegion") + str_i + ": " +
        this->high_resolution_string[i];
      comment.push_back(highres_string_i);
    }
  }


  // Add subsample resolution to array comment[].
  template <typename _IO_PARAM>
  void OUTPUT_PARAM_BASE<_IO_PARAM>::
  AddSubsampleResolutionComment
    (std::vector<std::string> & comment) const
  {
    if (this->flag_subsample) {
      std::string str;
      IJK::create_string("SubsampleResolution: ",
                         this->subsample_resolution, str);
      comment.push_back(str);
    }
  }


  // Add supersample resolution to array comment[].
  template <typename _IO_PARAM>
  void OUTPUT_PARAM_BASE<_IO_PARAM>::
  AddSupersampleResolutionComment
    (std::vector<std::string> & comment) const
  {
    if (this->flag_supersample) {
      std::string str;
      IJK::create_string("SupersampleResolution: ",
                         this->supersample_resolution, str);
      comment.push_back(str);
    }
  }


  template <typename _IO_PARAM>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void OUTPUT_PARAM_BASE<_IO_PARAM>::
  Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix << "dimension: " << dimension << "\n";
    out << line_prefix << "mesh_dimension: "
        << mesh_dimension << "\n";
    out << line_prefix << "num_vertices_per_isopoly: "
        << num_vertices_per_isopoly << "\n";
    _IO_PARAM::Print(out, line_prefix);
  }

}

#endif
  
