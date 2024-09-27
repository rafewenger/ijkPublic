/*!
 *  @file ijkgrid_nrrd.tpp
 *  @brief ijk templates for reading nrrd file into a scalar/vector grid.
 *  - Version 0.5.0
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

#ifndef _IJKGRID_NRRD_
#define _IJKGRID_NRRD_

#include <string>

#include "ijk.tpp"
#include "ijkNrrd.h"

namespace IJK {

  // **************************************************
  // TEMPLATE CLASS NRRD_DATA
  // **************************************************

  /// C++ class allocating/freeing Nrrd data
  template <typename DTYPE, typename ATYPE>
  class NRRD_DATA {

  protected:
    Nrrd * data;

  public:
    NRRD_DATA();
    ~NRRD_DATA();

    // Get functions
    DTYPE Dimension() const
    { return(data->dim); }            ///< Return dimension
    ATYPE AxisSize                    ///< Return # vertices on axis d.
    (const DTYPE d) const 
    { return(data->axis[d].size); }
    Nrrd * DataPtr()                  ///< Return Nrrd data structure
    { return(data); };
    Nrrd * DataPtrConst() const       ///< Return Nrrd data structure
    { return(data); };

    /// Get grid spacing.
    template <typename STYPE>
    void GetSpacing(std::vector<STYPE> & grid_spacing) const;

    // Set functions

    /// Set volume dimension and axis size.
    template <typename DTYPE2, typename ATYPE2>
    void SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size);

    /// Set spacing.
    template <typename SPACING_TYPE>
    void SetSpacing(const SPACING_TYPE spacing[]);
    template <typename SPACING_TYPE>
    void SetSpacing(const std::vector<SPACING_TYPE> & spacing);

    void CopyHeader                   ///< Copy header.
    (const Nrrd * from_data);
    void CopyAxisInfo                 ///< Copy axis information.
    (const Nrrd * from_data);
    void CopyKeyValues                ///< Copy key values.
    (const Nrrd * from_data);
    void CopyComments                 ///< Copy comments.
    (const Nrrd * from_data);

    /// Add key and value.
    void AddKeyValue(const char * key, const char * value_string);
    void AddKeyValue(const char * key, const std::string & value_string);
    void AddKeyValue
    (const std::string & key, const std::string & value_string);

    /// Return true if dimension and axis sizes match.
    template <typename DTYPE2, typename ATYPE2>
    bool CheckSize
    (const DTYPE2 dimension, const ATYPE2 * axis_size,
     IJK::ERROR & error) const;

    /// Return true if dimension and axis sizes match.
    bool CheckAxisSize
    (const Nrrd * nrrd_data2, IJK::ERROR & error) const;

    /// Return true if scalar grid dimension and axis sizes match.
    template <class GTYPE>
    bool CheckScalarGridSize
    (const GTYPE & grid, IJK::ERROR & error) const;

    /// Return true if vector grid size matches nrrd size.
    /// if true, nrrd dimension = vector grid dimension + 1, and
    ///          nrrd axis_size[0] = grid.VectorLength(), and
    ///          nrrd_axis_size[d+1] = grid.AxisSize[d].
    template <class GTYPE>
    bool CheckVectorGridSize
    (const GTYPE & grid, IJK::ERROR & error) const;
  };

  // **************************************************
  // TEMPLATE CLASS GRID_NRRD_IN
  // **************************************************

  template <typename DTYPE, typename ATYPE>
  class GRID_NRRD_IN:protected NRRD_DATA<DTYPE, ATYPE> {

  protected:
    bool read_failed;            ///< True if read failed.

  public:
    GRID_NRRD_IN();              ///< Constructor
    ~GRID_NRRD_IN();             ///< Destructor
    
    const Nrrd * DataPtrConst() const       ///< Return Nrrd data structure
    { return(this->data); };

    // Get functions
    bool ReadFailed() const      ///< Return true if read failed.
    { return(read_failed); };

    /// Read nrrd file.
    void Read(const char * input_filename, IJK::ERROR & error);

    /// @brief Read scalar grid.
    template <typename SCALAR_GRID>
    void ReadScalarGrid(const char * input_filename, SCALAR_GRID & grid,
                        IJK::ERROR & error);

    /// @overload
    /// @brief Read scalar grid. (Version using std::string filename)
    template <typename SCALAR_GRID>
    void ReadScalarGrid(const std::string & input_filename, SCALAR_GRID & grid,
                        IJK::ERROR & error);

    /// @overload
    /// @brief Read scalar grid and return axis info.
    template <typename SCALAR_GRID, typename DTYPE2, typename ATYPE2>
    void ReadScalarGrid
    (const char * input_filename, SCALAR_GRID & grid,
     NRRD_DATA<DTYPE2,ATYPE2> & header, IJK::ERROR & error);

    /// Read vector grid.
    template <typename VECTOR_GRID>
    void ReadVectorGrid(const char * input_filename, VECTOR_GRID & grid,
                        IJK::ERROR & error);

    /// Read vector grid and return axis info.
    template <typename VECTOR_GRID, typename DTYPE2, typename ATYPE2>
    void ReadVectorGrid
    (const char * input_filename, VECTOR_GRID & grid,
     NRRD_DATA<DTYPE2,ATYPE2> & header, IJK::ERROR & error);
  };


  // **************************************************
  //! @name CHECK FUNCTIONS
  // **************************************************

  //!@{

  //! Return false and set error message if dimension is negative 
  //! or larger than NRRD_DIM_MAX.
  template <typename DTYPE>
  inline bool check_nrrd_dimension(const DTYPE dimension, IJK::ERROR & error)
  {
    if (dimension > NRRD_DIM_MAX) {
      error.AddMessage
        ("Dimension error. Dimension ", dimension, " is too large.");  
      error.AddMessage("  Max allowable dimension is ", NRRD_DIM_MAX, ".");
      return(false);
    }
    else if (dimension < 0) {
      error.AddMessage
        ("Dimension error.  Dimension ", dimension, " is negative.");
      error.AddMessage("  Dimension must be non-negative.");
      
      return(false);
    }

    return(true);
  }

  //!@}


  // **************************************************
  // FUNCTION add_nrrd_message
  // **************************************************

  /// Add nrrd error message.
  inline void add_nrrd_message
  (const char * msg_header, IJK::ERROR & error)
  {
    char * s = biffGetDone(NRRD);

    // Make a copy of string s
    std::string nrrd_error_msg = s;

    free(s);      // free error message string

    error.AddMessage(msg_header, nrrd_error_msg);
  }

  inline void add_nrrd_message(IJK::ERROR & error)
  {
    add_nrrd_message("  Nrrd error: ", error);
  }

  // **************************************************
  /// @name NRRD SET/COPY FUNCTIONS
  // **************************************************

  ///@{

  /// Set dimension and axis size of Nrrd data.
  template <typename DTYPE, typename ATYPE>
  inline void set_nrrd_size
  (const DTYPE dimension, const ATYPE * axis_size, Nrrd * nrrd_data)
  {
    size_t nrrd_axis_size[NRRD_DIM_MAX];
    IJK::PROCEDURE_ERROR error("set_nrrd_size");

    if (!check_nrrd_dimension(dimension, error)) { throw error; };

    nrrd_data->dim = dimension;

    for (DTYPE d = 0; d < dimension; d++) 
      { nrrd_axis_size[d] = axis_size[d]; }

    nrrdAxisInfoSet_nva(nrrd_data, nrrdAxisInfoSize, nrrd_axis_size);
  }


  /// Copy dimension from one Nrrd data structure to another.
  inline void copy_nrrd_dimension(const Nrrd * from_nrrd, Nrrd * to_nrrd)
  {
    to_nrrd->dim = from_nrrd->dim;
  }


  /// Copy axis information from one Nrrd data structure to another.
  inline void copy_nrrd_axis_info(const Nrrd * from_nrrd, Nrrd * to_nrrd)
  {
    IJK::PROCEDURE_ERROR error("copy_nrrd_axis_info");

    // Make sure that dimension is copied.
    copy_nrrd_dimension(from_nrrd, to_nrrd);

    bool copy_failed =
      nrrdAxisInfoCopy(to_nrrd, from_nrrd, NULL, 0);

    if (copy_failed) {
      error.AddMessage("Error copying Nrrd axis information.");
      add_nrrd_message(error);
      throw error;
    }
  }


  /// Set spacing of nrrd data.
  template <typename DTYPE, typename SPACING_TYPE>
  inline void set_nrrd_spacing
  (const DTYPE dimension, const SPACING_TYPE spacing[], Nrrd * nrrd_data)
  {
    double nrrd_spacing[NRRD_DIM_MAX];
    IJK::PROCEDURE_ERROR error("set_nrrd_spacing");

    if (!check_nrrd_dimension(dimension, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) 
      { nrrd_spacing[d] = spacing[d]; }

    nrrdAxisInfoSet_nva(nrrd_data, nrrdAxisInfoSpacing, nrrd_spacing);
  }


  /// Copy key values from one Nrrd data structure to another.
  inline void copy_nrrd_key_values(const Nrrd * from_nrrd, Nrrd * to_nrrd)
  {
    IJK::PROCEDURE_ERROR error("copy_nrrd_key_values");

    bool copy_failed = nrrdKeyValueCopy(to_nrrd, from_nrrd);

    if (copy_failed) {
      error.AddMessage("Error copying Nrrd key values.");
      add_nrrd_message(error);
      throw error;
    }
  }


  /// Copy comments from one Nrrd data structure to another.
  inline void copy_nrrd_comments(const Nrrd * from_nrrd, Nrrd * to_nrrd)
  {
    IJK::PROCEDURE_ERROR error("copy_nrrd_comments");

    bool copy_failed = nrrdCommentCopy(to_nrrd, from_nrrd);

    if (copy_failed) {
      error.AddMessage("Error copying Nrrd comments.");
      add_nrrd_message(error);
      throw error;
    }
  }

  /// Copy header information from one Nrrd data structure to another.
  /// - Copy dimension, axis size, key values and comments.
  inline void copy_nrrd_header(const Nrrd * from_nrrd, Nrrd * to_nrrd)
  {
    copy_nrrd_dimension(from_nrrd, to_nrrd);
    copy_nrrd_axis_info(from_nrrd, to_nrrd);
    copy_nrrd_key_values(from_nrrd, to_nrrd);
    copy_nrrd_comments(from_nrrd, to_nrrd);
  }

  ///@}


  // **************************************************
  /// @name WRITE FUNCTIONS
  // **************************************************

  ///@{

  /// Wrap scalar data in nrrd_data.
  template <typename DTYPE, typename ATYPE>
  void wrap_scalar_grid_data
  (Nrrd * nrrd_data, const float * scalar, 
   const DTYPE dimension, const ATYPE * axis_size)
  {
    size_t nrrd_axis_size[NRRD_DIM_MAX];
    IJK::PROCEDURE_ERROR error("wrap_scalar_grid_data");

    if (!check_nrrd_dimension(dimension, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) 
      { nrrd_axis_size[d] = axis_size[d]; }

    nrrdWrap_nva(nrrd_data, (void *)(scalar), nrrdTypeFloat,
                 dimension, nrrd_axis_size);
  }

  /// Wrap scalar data in nrrd_data.
  template <typename DTYPE, typename ATYPE>
  void wrap_scalar_grid_data
  (Nrrd * nrrd_data, const int * scalar, 
   const DTYPE dimension, const ATYPE * axis_size)
  {
    size_t nrrd_axis_size[NRRD_DIM_MAX];
    IJK::PROCEDURE_ERROR error("wrap_scalar_grid_data");

    if (!check_nrrd_dimension(dimension, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) 
      { nrrd_axis_size[d] = axis_size[d]; }

    nrrdWrap_nva(nrrd_data, (void *)(scalar), nrrdTypeInt,
                 dimension, nrrd_axis_size);
  }


  /// Wrap vector data in nrrd_data.
  template <typename DTYPE, typename ATYPE, typename LTYPE>
  void wrap_vector_grid_data
  (Nrrd * nrrd_data, const float * scalar, 
   const DTYPE dimension, const ATYPE * axis_size, 
   const LTYPE vector_length)
  {
    size_t nrrd_axis_size[NRRD_DIM_MAX];
    IJK::PROCEDURE_ERROR error("wrap_vector_grid_data");

    if (!check_nrrd_dimension(dimension+1, error)) { throw error; };

    nrrd_axis_size[0] = vector_length;
    for (DTYPE d = 0; d < dimension; d++) 
      { nrrd_axis_size[d+1] = axis_size[d]; }

    nrrdWrap_nva(nrrd_data, (void *)(scalar), nrrdTypeFloat,
                 dimension+1, nrrd_axis_size);
  }

  /// Write scalar grid in nrrd file.
  template <typename GTYPE>
  void write_scalar_grid_nrrd
  (const char * output_filename, const GTYPE & grid)
  {
    IJK::PROCEDURE_ERROR error("write_scalar_grid_nrrd");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();

    wrap_scalar_grid_data(data, grid.ScalarPtrConst(),
                          grid.Dimension(), grid.AxisSize());
    bool save_failed = nrrdSave(output_filename, data, NULL);

    nrrdNix(data);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }

  }

  /// Write scalar grid in nrrd file. C++ string version.
  template <typename GTYPE>
  void write_scalar_grid_nrrd
  (const std::string & output_filename, const GTYPE & grid)
  {
    write_scalar_grid_nrrd(output_filename.c_str(), grid);
  }

  /// Write scalar grid in nrrd file. Add header information.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_scalar_grid_nrrd
  (const char * output_filename, const GTYPE & grid, 
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    IJK::PROCEDURE_ERROR error("write_scalar_grid_nrrd");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();

    wrap_scalar_grid_data(data, grid.ScalarPtrConst(),
                          grid.Dimension(), grid.AxisSize());

    // Check that grid sizes match nrrd_header sizes.
    if (!nrrd_header.CheckAxisSize(data, error)) { 
      nrrdNix(data);
      throw error; 
    };

    copy_nrrd_header(nrrd_header.DataPtrConst(), data);

    bool save_failed = nrrdSave(output_filename, data, NULL);
    nrrdNix(data);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }
  }

  /// \brief Write scalar grid in nrrd file. Add header information. 
  /// C++ string version.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_scalar_grid_nrrd
  (const std::string & output_filename, const GTYPE & grid, 
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    write_scalar_grid_nrrd(output_filename.c_str(), grid, nrrd_header);
  }

  /// Write scalar grid in nrrd file. Compress data using gzip.
  template <typename GTYPE>
  void write_scalar_grid_nrrd_gzip
  (const char * output_filename, const GTYPE & grid)
  {
    IJK::PROCEDURE_ERROR error("write_scalar_grid_nrrd_gzip");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();
    NrrdIoState * nio = nrrdIoStateNew();

    nrrdIoStateEncodingSet(nio, nrrdEncodingGzip);

    wrap_scalar_grid_data(data, grid.ScalarPtrConst(),
                          grid.Dimension(), grid.AxisSize());
    bool save_failed = nrrdSave(output_filename, data, nio);

    nrrdNix(data);
    nrrdIoStateNix(nio);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }

  }

  /// \brief Write scalar grid in nrrd file. Compress data using gzip. 
  /// C++ string version.
  template <typename GTYPE>
  void write_scalar_grid_nrrd_gzip
  (const std::string & output_filename, const GTYPE & grid)
  {
    write_scalar_grid_nrrd_gzip(output_filename.c_str(), grid);
  }

  /// \brief Write scalar grid in nrrd file. Add header information.
  /// Compress data using gzip.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_scalar_grid_nrrd_gzip
  (const char * output_filename, const GTYPE & grid, 
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    IJK::PROCEDURE_ERROR error("write_scalar_grid_nrrd_gzip");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();
    NrrdIoState * nio = nrrdIoStateNew();

    nrrdIoStateEncodingSet(nio, nrrdEncodingGzip);

    wrap_scalar_grid_data(data, grid.ScalarPtrConst(),
                          grid.Dimension(), grid.AxisSize());

    // Check that grid sizes match nrrd_header sizes.
    if (!nrrd_header.CheckAxisSize(data, error)) { 
      nrrdNix(data);
      throw error; 
    };

    copy_nrrd_header(nrrd_header.DataPtrConst(), data);

    bool save_failed = nrrdSave(output_filename, data, nio);
    nrrdNix(data);
    nrrdIoStateNix(nio);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }
  }

  /// \brief Write scalar grid in nrrd file. Compress data using gzip.
  /// Add header information.
  /// C++ string version.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_scalar_grid_nrrd_gzip
  (const std::string & output_filename, const GTYPE & grid, 
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    write_scalar_grid_nrrd_gzip(output_filename.c_str(), grid, nrrd_header);
  }

  /// Write vector grid in nrrd file.
  template <typename GTYPE>
  void write_vector_grid_nrrd
  (const char * output_filename, const GTYPE & grid)
  {
    IJK::PROCEDURE_ERROR error("write_vector_grid_nrrd");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();
    wrap_vector_grid_data
      (data, grid.VectorPtrConst(),
       grid.Dimension(), grid.AxisSize(), grid.VectorLength());

    bool save_failed = nrrdSave(output_filename, data, NULL);
    nrrdNix(data);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }
  }

  /// \brief Write vector grid in nrrd file.
  /// C++ string version.
  template <typename GTYPE>
  void write_vector_grid_nrrd
  (const std::string & output_filename, const GTYPE & grid)
  {
    write_vector_grid_nrrd(output_filename.c_str(), grid);
  }

  /// Write vector grid in nrrd file. Add header information.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_vector_grid_nrrd
  (const char * output_filename, const GTYPE & grid,
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    IJK::PROCEDURE_ERROR error("write_vector_grid_nrrd");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();
    wrap_vector_grid_data
      (data, grid.VectorPtrConst(),
       grid.Dimension(), grid.AxisSize(), grid.VectorLength());

    // Check that grid sizes match nrrd_header sizes.
    if (!nrrd_header.CheckAxisSize(data, error)) {
      nrrdNix(data);
      throw error;
    };

    copy_nrrd_header(nrrd_header.DataPtrConst(), data);

    bool save_failed = nrrdSave(output_filename, data, NULL);
    nrrdNix(data);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }
  }

  /// \brief Write vector grid in nrrd file. Add header information.
  /// C++ string version.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_vector_grid_nrrd
  (const std::string & output_filename, const GTYPE & grid,
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    write_vector_grid_nrrd
      (output_filename.c_str(), grid, nrrd_header);
  }

  /// Write vector grid in nrrd file. Compress data using gzip.
  template <typename GTYPE>
  void write_vector_grid_nrrd_gzip
  (const char * output_filename, const GTYPE & grid)
  {
    IJK::PROCEDURE_ERROR error("write_vector_grid_nrrd_gzip");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();
    NrrdIoState * nio = nrrdIoStateNew();

    nrrdIoStateEncodingSet(nio, nrrdEncodingGzip);

    wrap_vector_grid_data
      (data, grid.VectorPtrConst(),
       grid.Dimension(), grid.AxisSize(), grid.VectorLength());

    bool save_failed = nrrdSave(output_filename, data, nio);
    nrrdNix(data);
    nrrdIoStateNix(nio);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }
  }

  /// \brief Write vector grid in nrrd file. Compress data using gzip.
  /// Add header information.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_vector_grid_nrrd_gzip
  (const char * output_filename, const GTYPE & grid,
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    IJK::PROCEDURE_ERROR error("write_vector_grid_nrrd_gzip");

    if (output_filename == NULL) {
      error.AddMessage("Programming error: Empty output filename.");
      throw error;
    }

    Nrrd * data = nrrdNew();
    NrrdIoState * nio = nrrdIoStateNew();

    nrrdIoStateEncodingSet(nio, nrrdEncodingGzip);

    wrap_vector_grid_data
      (data, grid.VectorPtrConst(),
       grid.Dimension(), grid.AxisSize(), grid.VectorLength());

    copy_nrrd_header(nrrd_header.DataPtrConst(), data);

    bool save_failed = nrrdSave(output_filename, data, nio);
    nrrdNix(data);
    nrrdIoStateNix(nio);

    if (save_failed) {
      error.AddMessage("Unable to save nrrd data to ", output_filename, ".");
      add_nrrd_message(error);
      throw error;
    }
  }

  /// \brief Write vector grid in nrrd file. Compress data using gzip.
  /// Add header information.
  /// C++ string version.
  template <typename GTYPE, typename DTYPE, typename ATYPE>
  void write_vector_grid_nrrd_gzip
  (const std::string & output_filename, const GTYPE & grid,
   const NRRD_DATA<DTYPE,ATYPE> & nrrd_header)
  {
    write_vector_grid_nrrd_gzip
      (output_filename.c_str(), grid, nrrd_header);
  }

  ///@}


  // **************************************************
  // CLASS NRRD_DATA MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename ATYPE>
  NRRD_DATA<DTYPE,ATYPE>::NRRD_DATA()
  {
    data = NULL;
    data = nrrdNew();
  }

  /// Destructor.
  template <typename DTYPE, typename ATYPE>
  NRRD_DATA<DTYPE,ATYPE>::~NRRD_DATA()
  {
    if (data != NULL) { nrrdNuke(data); }
    data = NULL;
  }

  /// Get grid spacing.
  template <typename DTYPE, typename ATYPE>
  template <typename STYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  GetSpacing(std::vector<STYPE> & grid_spacing) const
  {
    IJK::PROCEDURE_ERROR error("GRID_NNRD_IN::GetSpacing");

    grid_spacing.clear();

    if (this->DataPtrConst() == NULL) {
      error.AddMessage("Programming error.  Nrrd data structure is empty.");
      throw error;
    }

    DTYPE dimension = this->Dimension();
    if (dimension < 1) { return; }

    if (!check_nrrd_dimension(dimension, error)) { throw error; };

    double nrrd_spacing[NRRD_DIM_MAX];
    nrrdAxisInfoGet_nva
      (this->DataPtrConst(), nrrdAxisInfoSpacing, nrrd_spacing);

    grid_spacing.resize(dimension);
    for (DTYPE d = 0; d < dimension; d++) {
      if (AIR_EXISTS(nrrd_spacing[d])) 
        { grid_spacing[d] = nrrd_spacing[d]; }
      else {
        // default to 1.
        grid_spacing[d] = 1;
      }
    }
        
  }

  // Set dimensions and axis size.
  template <typename DTYPE, typename ATYPE>
  template <typename DTYPE2, typename ATYPE2>
  void NRRD_DATA<DTYPE,ATYPE>::
  SetSize(const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    set_nrrd_size(dimension, axis_size, this->data);
  }


  // Set spacing.
  template <typename DTYPE, typename ATYPE>
  template <typename SPACING_TYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  SetSpacing(const SPACING_TYPE spacing[])
  {
    set_nrrd_spacing(this->Dimension(), spacing, this->data);
  }


  // Set spacing.
  // - Version using C++ STL vector for spacing.
  template <typename DTYPE, typename ATYPE>
  template <typename SPACING_TYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  SetSpacing(const std::vector<SPACING_TYPE> & spacing)
  {
    set_nrrd_spacing
      (this->Dimension(), IJK::vector2pointer(spacing), this->data);
  }

  // Copy header.
  template <typename DTYPE, typename ATYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  CopyHeader(const Nrrd * from_data)
  {
    copy_nrrd_header(from_data, this->data);
  }

  /// Copy axis information.
  template <typename DTYPE, typename ATYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  CopyAxisInfo(const Nrrd * from_data)
  {
    copy_nrrd_axis_info(from_data, this->data);
  }

  /// Copy key values
  template <typename DTYPE, typename ATYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  CopyKeyValues(const Nrrd * from_data)
  {
    copy_nrrd_key_values(from_data, this->data);
  }

  /// Copy comments
  template <typename DTYPE, typename ATYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  CopyComments(const Nrrd * from_data)
  {
    copy_nrrd_comments(from_data, this->data);
  }

  /// Add key and value.
  template <typename DTYPE, typename ATYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  AddKeyValue(const char * key, const char * value_string)
  {
    nrrdKeyValueAdd(this->DataPtr(), key, value_string);
  }

  /// Add key and value.
  template <typename DTYPE, typename ATYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  AddKeyValue(const char * key, const std::string & value_string)
  {
    AddKeyValue(key, value_string.c_str());
  }

  /// Add key and value.
  template <typename DTYPE, typename ATYPE>
  void NRRD_DATA<DTYPE,ATYPE>::
  AddKeyValue
  (const std::string & key, const std::string & value_string)
  {
    AddKeyValue(key.c_str(), value_string.c_str());
  }

  /// Return true if dimension and axis sizes match.
  /// Otherwise, return false and set error message.
  template <typename DTYPE, typename ATYPE>
  template <typename DTYPE2, typename ATYPE2>
  bool NRRD_DATA<DTYPE,ATYPE>::
  CheckSize(const DTYPE2 dimension, const ATYPE2 * axis_size,
            IJK::ERROR & error) const
  {
    if (this->Dimension() != dimension) {
      error.AddMessage("Incorrect Nrrd dimension.");
      error.AddMessage
        ("  Nrrd dimension = ", this->Dimension(), ".");
      error.AddMessage("  Should be = ", dimension, ".");
      return(false);
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (AxisSize(d) != axis_size[d]) {
        error.AddMessage("Incorrect axis_size[", d, "].");
        error.AddMessage
          ("  axis_size[", d, "] = ", this->AxisSize(d), ".");
        error.AddMessage("  Should be = ", axis_size[d], ".");
        return(false);
      }
    }

    return(true);
  }

  /// Return true if dimension and axis sizes match.
  /// Otherwise, return false and set error message.
  template <typename DTYPE, typename ATYPE>
  bool NRRD_DATA<DTYPE,ATYPE>::
  CheckAxisSize(const Nrrd * data2, IJK::ERROR & error) const
  {
    const DTYPE dimension = Dimension();

    if (data2 == NULL) {
      error.AddMessage("Programming error. Nrrd data structure data2 has not been allocated.");
      error.AddMessage
        ("Call \"data = nrrdNew()\" before calling this procedure.");
      throw error;
    }

    if (dimension != data2->dim) {
      error.AddMessage("Programming error. Dimension ", dimension, 
                       " of this nrrd data structure");
      error.AddMessage("  does not match dimension ", data2->dim, 
                       " of data2.");
      return (false);
    }

    if (dimension < 1) { return(true); };

    if (!check_nrrd_dimension(dimension, error)) { throw error; };

    size_t axis_size[NRRD_DIM_MAX];
    nrrdAxisInfoGet_nva
      (this->DataPtrConst(), nrrdAxisInfoSize, axis_size);

    return(CheckSize(dimension, axis_size, error));
  }

  /// Return true if dimension and axis sizes match.
  /// Otherwise, return false and set error message.
  template <typename DTYPE, typename ATYPE>
  template <typename GTYPE>
  bool NRRD_DATA<DTYPE,ATYPE>::
  CheckScalarGridSize(const GTYPE & grid, IJK::ERROR & error) const
  {
    return(CheckSize(grid.Dimension(), grid.AxisSize(), error));
  }

  /// Return true if vector grid size matches nrrd size.
  template <typename DTYPE, typename ATYPE>
  template <typename GTYPE>
  bool NRRD_DATA<DTYPE,ATYPE>::
  CheckVectorGridSize(const GTYPE & grid, IJK::ERROR & error) const
  {
    const DTYPE dimension = Dimension();

    if (dimension < 1) {
      error.AddMessage("Programming error: Illegal nrrd dimension ",
                       Dimension(), ".");
      error.AddMessage("  Nrrd dimension must be at least 1 for vector data.");
      throw error;
    }

    IJK::ARRAY<size_t> axis_size(dimension);
    axis_size[0] = grid.VectorLength();
    for (DTYPE d = 1; d < dimension; d++)
      { axis_size[d] = grid.AxisSize(d-1); }

    return(CheckSize(dimension, axis_size.PtrConst(), error));
  }

  // **************************************************
  // CLASS GRID_NRRD_IN MEMBER FUNCTIONS
  // **************************************************

  // Constructor.
  template <typename DTYPE, typename ATYPE>
  GRID_NRRD_IN<DTYPE,ATYPE>::GRID_NRRD_IN()
  {
    read_failed = false;
  }

  
  // Destructor.
  template <typename DTYPE, typename ATYPE>
  GRID_NRRD_IN<DTYPE,ATYPE>::~GRID_NRRD_IN()
  {
    read_failed = false;
  }

  
  // Read nrrd file.
  template <typename DTYPE, typename ATYPE>
  void GRID_NRRD_IN<DTYPE,ATYPE>::
  Read(const char * input_filename, IJK::ERROR & read_error)
  {
    IJK::PROCEDURE_ERROR error("GRID_NRRD_IN::Read");

    if (input_filename == NULL) {
      read_failed = true;
      error.AddMessage("Programming error: Empty input filename.");
      throw error;
    }

    // *** NOTE: SHOULD PROBABLY CALL nrrdNuke BEFORE nrrdLoad ***

    read_failed = nrrdLoad(this->data, input_filename, NULL);

    if (read_failed) {
      read_error.AddMessage("Error reading: ", input_filename);
      add_nrrd_message(read_error);
      return;
    }
  }

  
  // Read scalar grid.
  template <typename DTYPE, typename ATYPE>
  template <typename SCALAR_GRID>
  void GRID_NRRD_IN<DTYPE,ATYPE>::
  ReadScalarGrid(const char * input_filename, SCALAR_GRID & grid,
                 IJK::ERROR & read_error)
  {
    Read(input_filename, read_error);
    if (ReadFailed()) { return; }

    const DTYPE dimension = this->Dimension();

    if (dimension < 1) {
      grid.SetSize(0, (ATYPE *)(NULL));
      return;
    }

    size_t size[NRRD_DIM_MAX];
    nrrdAxisInfoGet_nva(this->data, nrrdAxisInfoSize, size);
    grid.SetSize(dimension, size);
    nrrd2scalar(this->data, grid.ScalarPtr());
  };


  // Read scalar grid. (Version using std::string filename.)
  template <typename DTYPE, typename ATYPE>
  template <typename SCALAR_GRID>
  void GRID_NRRD_IN<DTYPE,ATYPE>::
  ReadScalarGrid
  (const std::string & input_filename, SCALAR_GRID & grid,
   IJK::ERROR & error)
  {
    ReadScalarGrid(input_filename.c_str(), grid, error);
  }

  
  // Read scalar grid and return axis info.
  template <typename DTYPE, typename ATYPE>
  template <typename SCALAR_GRID, typename DTYPE2, typename ATYPE2>
  void GRID_NRRD_IN<DTYPE,ATYPE>::
  ReadScalarGrid(const char * input_filename, SCALAR_GRID & grid,
                 NRRD_DATA<DTYPE2,ATYPE2> & header, IJK::ERROR & error)
  {
    ReadScalarGrid(input_filename, grid, error);

    if (ReadFailed()) { return; };

    header.CopyHeader(this->DataPtrConst());
  }

  
  // Read vector grid.
  template <typename DTYPE, typename ATYPE>
  template <typename VECTOR_GRID>
  void GRID_NRRD_IN<DTYPE,ATYPE>::
  ReadVectorGrid(const char * input_filename, VECTOR_GRID & grid,
                 IJK::ERROR & read_error)
  {
    IJK::PROCEDURE_ERROR error("GRID_NRRD_IN::ReadVectorGrid");

    Read(input_filename, read_error);
    if (ReadFailed()) { return; }

    if (this->Dimension() < 1) {
      error.AddMessage("Illegal nrrd dimension for vector nrrd file.");
      error.AddMessage("  Dimension in nrrd file must be at least 1.");
      throw error;
    }

    const DTYPE dimension = this->Dimension()-1;

    size_t size[NRRD_DIM_MAX];
    nrrdAxisInfoGet_nva(this->data, nrrdAxisInfoSize, size);

    if (dimension == 0) {
      grid.SetSize(dimension, (ATYPE *) NULL, size[0]);
      return;
    }

    grid.SetSize(dimension, size+1, size[0]);
    nrrd2scalar(this->data, grid.VectorPtr());
  }

  /// Read vector grid and return axis info.
  template <typename DTYPE, typename ATYPE>
  template <typename VECTOR_GRID, typename DTYPE2, typename ATYPE2>
  void GRID_NRRD_IN<DTYPE,ATYPE>::
  ReadVectorGrid(const char * input_filename, VECTOR_GRID & grid,
                 NRRD_DATA<DTYPE2,ATYPE2> & header, IJK::ERROR & error)
  {
    ReadVectorGrid(input_filename, grid, error);

    if (ReadFailed()) { return; };

    header.CopyHeader(this->DataPtrConst());
  }

}

#endif

