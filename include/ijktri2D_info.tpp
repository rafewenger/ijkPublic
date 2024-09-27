/*!
 * @file ijktri2D_info.tpp
 * @brief ijk templates for storing and setting triangulations of 2D polygons.
 * - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2019-2023 Rephael Wenger

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

#ifndef _IJKTRI2D_INFO_
#define _IJKTRI2D_INFO_

#define _USE_MATH_DEFINES

#include <bitset>
#include <cmath>

#include "ijk.tpp"
#include "ijktri2D_encoding.tpp"

#include "ijk_macros.h"

namespace IJK {

  // ***************************************************************
  // CLASSES FOR POLYGON TRIANGULATION RESULTS
  // ***************************************************************

  /// Data structure for results of triangulation of polygons.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  class POLYGON_TRIANGULATION_RESULT {

  public:
    typedef COS_TYPE_X COS_TYPE;
    typedef NTYPE NUMBER_TYPE;

    static const int BitSetSize() { return(BIT_SET_SIZE); };


  public:
    /// Cosine of min angle in triangulations.
    COS_TYPE_X cos_min_triangulation_angle;

    /// If true, all triangulations have a zero length edge.
    bool flag_zero;

    /// Encoding of triangulation as 1's and 0's.
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> triangulation_encoding;

    /// Number of triangulation triangles.
    NTYPE num_triangles;

    /*!
     *  @brief Index of vertex used in triangulation.
     *  - Usually all (most) triangles are incident
     *    on polygon vertex tri_vertex_index.
     *  - Note: This index has values 0,1,...,(num_poly_vert-1).
     *  - Index is with respect to the polygon.  It does not represent
     *    the vertex index with respect to all the vertices in the mesh.
     */
    NTYPE tri_vertex_index;

    /// Number of interior vertices used in the triangulation.
    NTYPE num_interior_tri_vertices;

    /// Return number of interior triangulation vertices.
    NTYPE NumInteriorTriVertices() const
    { return(num_interior_tri_vertices); }

    /*!
     *  @brief Set cos_min_triangulation_angle, num_triangles, tri_vertex_index.
     *  - Set num_interior_tri_vertices to 0.
     *  - Set flag_zero to false.
     */
    void SetNoInterior(const COS_TYPE_X cos_min, const NTYPE numt,
                       const NTYPE triv_index)
    {
      cos_min_triangulation_angle = cos_min;
      num_triangles = numt;
      tri_vertex_index = triv_index;
      num_interior_tri_vertices = 0;
      flag_zero = false;
    }


    /*!
     *  @brief Set cos_min_triangulation_angle, num_triangles, and num_interior_tri_vertices.
     *  - Set flag_zero to false.
     *  - Does not modify tri_vertex_index.
     */
    void SetInterior(const COS_TYPE_X cos_min, const NTYPE numt,
                     const NTYPE num_interior_triv)
    {
      cos_min_triangulation_angle = cos_min;
      num_triangles = numt;
      num_interior_tri_vertices = num_interior_triv;
      flag_zero = false;
    }
      

    /// Clear.
    void Clear() {
      cos_min_triangulation_angle = 1.0;

      flag_zero = false;
      num_triangles = 0;
      tri_vertex_index = 0;
      num_interior_tri_vertices = 0;
    };


    /// Constructor
    POLYGON_TRIANGULATION_RESULT() 
    { Clear(); }

    /// Copy.
    template <typename TRI_RESULT_TYPE>
    const POLYGON_TRIANGULATION_RESULT & Copy
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };


    /*!
     *  @brief Copy POLYGON_TRIANGULATION_RESULT.
     *  - Specify copying of this specific class.
     *  - Useful if copying is overwritten in inherited classes.
     */
    template <typename TRI_RESULT_TYPE>
    const POLYGON_TRIANGULATION_RESULT & CopyPolyTriangulationResult
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };


    // Print routines (mainly for debugging.)

    /// Print triangulation information.
    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

    ///  Print cosine of the min angle in the triangulation.
    template <typename OSTREAM_TYPE, typename STYPE0>
    void PrintCosMinTriangulationAngle
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

    ///  Print min angle in the triangulation in degrees.
    template <typename OSTREAM_TYPE, typename STYPE0>
    void PrintMinTriangulationAngle
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

  };


  /*!
   *  @brief Data structure for results of triangulation of polygons.
   *  - NEW VERSION using ear bits to specify triangulation.
   *  @tparam MAX_NUM_POLYGON_VERTICES Maximum number of polygon vertices.
   */
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  class POLYGON_TRIANGULATION_RESULT_E:
    public POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE, COS_TYPE_X,NTYPE> {

  protected:
    static const int MAX_NUM_POLYGON_VERTICES = BIT_SET_SIZE/2;

  public:
    typedef COS_TYPE_X COS_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:

    /*!
     *  @brief Bits to flag which vertices are in triangulation ears.
     *  - If bit i is 0, then i'th polygon vertex is in a triangulation ear.
     *  - Index is with respect to the polygon.  It does not represent
     *    the vertex index with respect to all the vertices in the mesh.
     */
    std::bitset<MAX_NUM_POLYGON_VERTICES> ear_flag;

    /*!
     *  @brief Return true if vertex ivert and two adjacent polygon vertices form a triangulation ear.
     *  @param jv_loc Location of vertex in polygon.
     *    - jv_loc is in range [0,NumPolygonVertices(0-1)].
     */
    template <typename ITYPE2>
    bool IsEar(const ITYPE2 jv_loc) const
    { return(bool(ear_flag[jv_loc])); }

    /*!
     *  @brief Return true if triangulation has some ear, i.e. if ear_flag is not zero.
     */
    bool HasEar() const
    { return((ear_flag != 0)); }

    /*!
     *  @brief Return location of first ear in polygon.
     *  - Return false if no ear.
     */
    template <typename ITYPE1>
    bool GetFirstEar(ITYPE1 & iloc) const;

    /*!
     *  @brief Set ear_flag[jv_loc] to 0 or 1.
     *  @param jv_loc Location of vertex in polygon.
     *    - jv_loc is in range [0,NumPolygonVertices(0-1)].
     */
    template <typename ITYPE2>
    void SetEarFlag(const ITYPE2 jv_loc, const bool flag)
    { ear_flag[jv_loc] = flag; }

    /*!
     *  @brief Rotate ear_flag bits k locations to left.
     *  - Leftmost bits are appended on right.
     */
    template <typename ITYPE2, typename NTYPE2>
    void RotateLeftEarFlags(const ITYPE2 k, const NTYPE2 num_polyv);

    /*!
     *  @brief Rotate ear_flag bits k locations to right.
     *  - Rightmost bits are appended on left.
     */
    template <typename ITYPE2, typename NTYPE2>
    void RotateRightEarFlags(const ITYPE2 k, const NTYPE2 num_polyv);

    /// Clear ear_flag.  Set all bits to 0.
    void ClearEarFlag()
    { ear_flag.reset(); }
      

    /// Clear.
    void Clear() {
      POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::Clear();
      ClearEarFlag();
    };


    /// Constructor
    POLYGON_TRIANGULATION_RESULT_E() 
    { Clear(); }

    /// Copy.
    template <typename TRI_RESULT_TYPE>
    const POLYGON_TRIANGULATION_RESULT_E & Copy
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };

    /*!
     *  @brief Copy POLYGON_TRIANGULATION_RESULT_E.
     *  - Specify copying of this specific class.
     *  - Useful if copying is overwritten in inherited classes.
     */
    template <typename TRI_RESULT_TYPE>
    const POLYGON_TRIANGULATION_RESULT_E & CopyPolyTriangulationResultE
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };

    // Print routines (mainly for debugging.)

    /// Print first ear, if it exists.
    template <typename OSTREAM_TYPE, typename STYPE0, typename STYPE1>
    void PrintEar(OSTREAM_TYPE & out, 
                  const STYPE0 & prefix, const STYPE1 & suffix) const;

  };


  /*!
   *  @brief Data structure for results of triangulation of polygons with bit set size 16.
   *  @tparam MAX_NUM_VERTICES Maximum number of polygon vertices.
   */
  template <typename COS_TYPE_X, typename NTYPE>
  class POLYGON_TRIANGULATION_RESULT_E16:
    public POLYGON_TRIANGULATION_RESULT_E<16, COS_TYPE_X, NTYPE> {
  };


  // ***************************************************************
  // CLASSES FOR POLYGON TRIANGULATION SETTINGS
  // ***************************************************************

  /*!
   *  @brief Settings for polygon triangulation routines.
   *  - Note: Because all fields are boolean, this IS NOT a template class.
   *  - Member functions are defined as inline, to avoid compiler multi-definition erros.
   */
  class POLYGON_TRIANGULATION_SETTINGS {

  public:

    bool flag_triangulate_from_split_vertices;
    bool flag_triangulate_from_original_vertices;
    bool flag_triangulate_from_interior_vertex;
    bool flag_split_ear_triangulate_from_interior_vertex;
    bool flag_split_ear_triangulate_from_original_vertices;

    /// @brief Try all possible triangulations from polygon vertices.
    /// - Does not include triangulations from interior vertex.
    bool flag_all_triangulations_from_polygon_vertices;

    /// #brief Try all possible triangulations from hexagon vertices.
    /// - Ignored if flag_all_triangulations_from_polygon_vertices is true.
    bool flag_all_triangulations_from_hexagon_vertices;


  public:

    /// Return true if all triangulation flags are zero.
    inline bool AreAllTriangulationFlagsZero() const;

    /// Set all flags to flag.
    inline void SetAll(const bool flag);

    /*!
     *  @brief Set triangulation flags for split vertices, original vertices or interior vertex to flag.
     *  - Does not modify/initialize other flags.
     */
    inline void SetVorI(const bool flag);

    /// @brief Set triangulation flags for split vertices and original vertices.
    /// - Does not modify/initialize other flags.
    inline void SetV(const bool flag);

    /// #brief Set triangulation flags for splitting ears.
    /// - Does not modify/initialize other flags.
    inline void SetSplitEar(const bool flag);

    /// Or current settings with settings2.
    inline void OrSettings(const POLYGON_TRIANGULATION_SETTINGS & settings2);
    
    /// Set all flags to false.
    inline void Clear() { SetAll(false); };


    /// Constructor.
    POLYGON_TRIANGULATION_SETTINGS() { Clear(); }


    // Print routines (mainly for debugging.)

    /// Print triangulation settings.
    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;
  };


  // ***************************************************************
  // Member functions for POLYGON_TRIANGULATION_RESULT
  // ***************************************************************

  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
  PrintCosMinTriangulationAngle
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix << "cos_min_triangulation_angle: "
        << cos_min_triangulation_angle << "\n";
  }

  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
  PrintMinTriangulationAngle
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    // Only compute min_triangulation_angle for output.
    // You should never need/use angles for mesh processing.
    const COS_TYPE_X min_triangulation_angle =
      acos(cos_min_triangulation_angle) * (180.0/M_PI);

    out << line_prefix << "min triangulation angle (degrees): "
        << min_triangulation_angle << "\n";
  }


  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::Print
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    PrintCosMinTriangulationAngle(out, line_prefix);
    PrintMinTriangulationAngle(out, line_prefix);
    out << line_prefix << "flag_zero: " << int(flag_zero) << "\n";
    out << line_prefix << "num_triangles: " << num_triangles << "\n";
    out << line_prefix << "tri_vertex_index: " << tri_vertex_index << "\n";
    out << line_prefix << "num_interior_tri_vertices: "
        << num_interior_tri_vertices << "\n";

    std::string s0 = line_prefix;
    s0 = s0 + "triangulation representation: ";
    triangulation_encoding.PrintTriangulationRepresentation
      (out, s0, "\n");
  }


  // ***************************************************************
  // Member functions for POLYGON_TRIANGULATION_RESULT_E
  // ***************************************************************


  // Return location of first ear in polygon.
  // - Undefined if triangulation has no ear.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename ITYPE1>
  bool POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
  GetFirstEar(ITYPE1 & iloc) const
  {
    // Initialize iloc.
    iloc = 0;

    if (!HasEar()) { return(false); }

    while (iloc < MAX_NUM_POLYGON_VERTICES) {

      if (IsEar(iloc)) { return(true); }

      iloc++;
    }

    iloc = 0;
    return(false);
  }


  // Rotate ear_flag bits k locations to left.
  // - Leftmost bits are appended on right.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename ITYPE2, typename NTYPE2>
  void POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
  RotateLeftEarFlags(const ITYPE2 k, const NTYPE2 num_polyv)
  {
    const std::bitset<MAX_NUM_POLYGON_VERTICES> old_ear_flag(ear_flag);

    ear_flag.reset();
    for (NTYPE i = 0; i < num_polyv; i++) {
      if (old_ear_flag[i]) {
        const ITYPE2 k2 = (k+i) % num_polyv;
        SetEarFlag(k2, true);
      }
    }
  }


  // Rotate ear_flag bits k locations to right.
  // - Rightmost bits are appended on left.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename ITYPE2, typename NTYPE2>
  void POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
  RotateRightEarFlags(const ITYPE2 k, const NTYPE2 num_polyv)
  {
    const std::bitset<MAX_NUM_POLYGON_VERTICES> old_ear_flag(ear_flag);

    ear_flag.reset();
    for (NTYPE i = 0; i < num_polyv; i++) {
      if (old_ear_flag[i]) {
        const ITYPE2 k2 = ((i+num_polyv)-k) % num_polyv;
        SetEarFlag(k2, true);
      }
    }
  }


  // Print first ear, if it exists.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0, typename STYPE1>
  void POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
  PrintEar(OSTREAM_TYPE & out, 
           const STYPE0 & prefix, const STYPE1 & suffix) const
  {
    NTYPE iloc_ear;

    if (GetFirstEar(iloc_ear)) {
      out << prefix << iloc_ear << suffix;
    }
  }


  // ***************************************************************
  // MEMBER FUNCTIONS FOR POLYGON TRIANGULATION SETTINGS
  // ***************************************************************


  // Return true if all flags are zero.
  inline bool POLYGON_TRIANGULATION_SETTINGS::AreAllTriangulationFlagsZero() const
  {
    const bool flag = (flag_triangulate_from_split_vertices ||
                       flag_triangulate_from_original_vertices ||
                       flag_triangulate_from_interior_vertex ||
                       flag_split_ear_triangulate_from_interior_vertex ||
                       flag_split_ear_triangulate_from_original_vertices ||
                       flag_all_triangulations_from_polygon_vertices ||
                       flag_all_triangulations_from_hexagon_vertices);

    return (!flag);
  }


  inline void POLYGON_TRIANGULATION_SETTINGS::SetAll(const bool flag)
  {
    SetVorI(flag);
    SetSplitEar(flag);
    flag_all_triangulations_from_polygon_vertices = flag;
    flag_all_triangulations_from_hexagon_vertices = flag;
  }


  // Set triangulation flags for split vertices, original vertices 
  //    or interior vertex to flag.
  // - Does not modify/initialize other flags.
  inline void POLYGON_TRIANGULATION_SETTINGS::SetVorI(const bool flag)
  {
    flag_triangulate_from_split_vertices = flag;
    flag_triangulate_from_original_vertices = flag;
    flag_triangulate_from_interior_vertex = flag;
  }


  // Set triangulation flags for split vertices or original vertices.
  // - Does not modify/initialize other flags.
  inline void POLYGON_TRIANGULATION_SETTINGS::SetV(const bool flag)
  {
    flag_triangulate_from_split_vertices = flag;
    flag_triangulate_from_original_vertices = flag;
  }


  // Set triangulation flags for splitting ears.
  // - Does not modify/initialize other flags.
  inline void POLYGON_TRIANGULATION_SETTINGS::SetSplitEar(const bool flag)
  {
    flag_split_ear_triangulate_from_interior_vertex = flag;
    flag_split_ear_triangulate_from_original_vertices = flag;
  }


  // Or current settings with settings2.
  inline void POLYGON_TRIANGULATION_SETTINGS::OrSettings
  (const POLYGON_TRIANGULATION_SETTINGS & settings2)
  {
    IJK_OR_VARIABLE(flag_triangulate_from_split_vertices, settings2);
    IJK_OR_VARIABLE(flag_triangulate_from_original_vertices, settings2);
    IJK_OR_VARIABLE(flag_triangulate_from_interior_vertex, settings2);
    IJK_OR_VARIABLE(flag_split_ear_triangulate_from_interior_vertex, settings2);
    IJK_OR_VARIABLE(flag_split_ear_triangulate_from_original_vertices, settings2);
    IJK_OR_VARIABLE(flag_all_triangulations_from_polygon_vertices, settings2);
    IJK_OR_VARIABLE(flag_all_triangulations_from_hexagon_vertices, settings2);
  }


  template <typename OSTREAM_TYPE, typename STYPE0>
  void POLYGON_TRIANGULATION_SETTINGS::Print
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix << "flag_triangulate_from_split_vertices: " 
        << int(flag_triangulate_from_split_vertices) << "\n";
    out << line_prefix << "flag_triangulate_from_original_vertices: "
        << int(flag_triangulate_from_original_vertices) << "\n";
    out << line_prefix << "flag_triangulate_from_interior_vertex: "
        << int(flag_triangulate_from_interior_vertex) << "\n";
    out << line_prefix << "flag_split_ear_triangulate_from_interior_vertex: "
        << int(flag_split_ear_triangulate_from_interior_vertex) << "\n";
    out << line_prefix << "flag_split_ear_triangulate_from_original_vertices: "
        << int(flag_split_ear_triangulate_from_original_vertices) << "\n";
    out << line_prefix << "flag_all_triangulations_from_polygon_vertices: "
        << int(flag_all_triangulations_from_polygon_vertices) << "\n";
    out << line_prefix << "flag_all_triangulations_from_hexagon_vertices: "
        << int(flag_all_triangulations_from_hexagon_vertices) << "\n";
  }

}

#endif
