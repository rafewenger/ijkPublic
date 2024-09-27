/// \file ijktriangulate.tpp
/// @brief ijk templates for triangulating polytopes (DEPRECATED).
/// - Version 0.4.0
/// - Replaced by ijktri2D_info.tpp.

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2019-2021 Rephael Wenger

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

#ifndef _IJKTRIANGULATE_
#define _IJKTRIANGULATE_

#include <bitset>

#include "ijk.tpp"
#include "ijktri2D_encoding.tpp"


namespace IJK {

  // ***************************************************************
  // CLASSES FOR POLY TRIANGULATION RESULTS
  // ***************************************************************

  /// Data structure for results of triangulation of polygons.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  class POLY_TRIANGULATION_RESULT {

  public:
    typedef COS_TYPE_X COS_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:
    /// Cosine of min angle in triangulations.
    COS_TYPE_X cos_min_triangulation_angle;

    /// If true, all triangulations have a zero length edge.
    bool flag_zero;

    /// Encoding of triangulation as 1's and 0's.
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> triangulation_encoding;

    /// Number of triangulation triangles.
    NTYPE num_triangles;

    /// Index of vertex used in triangulation.
    /// - Usually all (most) triangles are incident 
    ///   on polygon vertex tri_vertex_index.
    /// - Note: This index has values 0,1,...,(num_poly_vert-1).
    /// - Index is with respect to the polygon.  It does not represent
    ///   the vertex index with respect to all the vertices in the mesh.
    NTYPE tri_vertex_index;

    /// Number of interior vertices used in the triangulation.
    NTYPE num_interior_tri_vertices;

    /// Return number of interior triangulation vertices.
    NTYPE NumInteriorTriVertices() const
    { return(num_interior_tri_vertices); }

    /// Set cos_min_triangulation_angle, num_triangles, tri_vertex_index.
    /// - Set num_interior_tri_vertices to 0.
    /// - Set flag_zero to false.
    void SetNoInterior(const COS_TYPE_X cos_min, const NTYPE numt,
                       const NTYPE triv_index)
    {
      cos_min_triangulation_angle = cos_min;
      num_triangles = numt;
      tri_vertex_index = triv_index;
      num_interior_tri_vertices = 0;
      flag_zero = false;
    }


    /// Set cos_min_triangulation_angle, num_triangles,
    ///   and num_interior_tri_vertices.
    /// - Set flag_zero to false.
    /// - Does not modify tri_vertex_index.
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
    POLY_TRIANGULATION_RESULT() 
    { Clear(); }

    /// Copy.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULT & Copy
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };

    /// Copy POLY_TRIANGULATION_RESULT.
    /// - Specify copying of this specific class.
    /// - Useful if copying is overwritten in inherited classes.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULT & CopyPolyTriangulationResult
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };


    // Print routines (mainly for debugging.)

    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

    template <typename OSTREAM_TYPE, typename STYPE0>
    void PrintCosMinTriangulationAngle
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

  };


  /// Data structure for results of triangulation of polygons.
  /// - NEW VERSION using ear bits to specify triangulation.
  /// @tparam MAX_NUM_POLYGON_VERTICES Maximum number of polygon vertices.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  class POLY_TRIANGULATION_RESULT_E:
    public POLY_TRIANGULATION_RESULT<BIT_SET_SIZE, COS_TYPE_X,NTYPE> {

  protected:
    static const int MAX_NUM_POLYGON_VERTICES = BIT_SET_SIZE/2;

  public:
    typedef COS_TYPE_X COS_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:
    /// Bits to flag which vertices are in triangulation ears.
    /// - If bit i is 0, then i'th polygon vertex is in a triangulation ear.
    /// - Index is with respect to the polygon.  It does not represent
    ///   the vertex index with respect to all the vertices in the mesh.
    std::bitset<MAX_NUM_POLYGON_VERTICES> ear_flag;

    /// Return true if vertex ivert and two adjacent polygon vertices
    ///   form a triangulation ear.
    /// @param jv_loc Location of vertex in polygon.
    ///   - jv_loc is in range [0,NumPolygonVertices(0-1)].
    template <typename ITYPE2>
    bool IsEar(const ITYPE2 jv_loc) const
    { return(bool(ear_flag[jv_loc])); }

    /// Return true if triangulation has some ear, 
    ///   i.e. if ear_flag is not zero.
    bool HasEar() const
    { return((ear_flag != 0)); }

    /// Return location of first ear in polygon.
    /// - Return false if no ear.
    template <typename ITYPE1>
    bool GetFirstEar(ITYPE1 & iloc) const;

    /// Set ear_flag[jv_loc] to 0 or 1.
    /// @param jv_loc Location of vertex in polygon.
    ///   - jv_loc is in range [0,NumPolygonVertices(0-1)].
    template <typename ITYPE2>
    void SetEarFlag(const ITYPE2 jv_loc, const bool flag)
    { ear_flag[jv_loc] = flag; }

    /// Rotate ear_flag bits k locations to left.
    /// - Leftmost bits are appended on right.
    template <typename ITYPE2, typename NTYPE2>
    void RotateLeftEarFlags(const ITYPE2 k, const NTYPE2 num_polyv);

    /// Rotate ear_flag bits k locations to right.
    /// - Rightmost bits are appended on left.
    template <typename ITYPE2, typename NTYPE2>
    void RotateRightEarFlags(const ITYPE2 k, const NTYPE2 num_polyv);

    /// Clear ear_flag.  Set all bits to 0.
    void ClearEarFlag()
    { ear_flag.reset(); }
      

    /// Clear.
    void Clear() {
      POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::Clear();
      ClearEarFlag();
    };


    /// Constructor
    POLY_TRIANGULATION_RESULT_E() 
    { Clear(); }

    /// Copy.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULT_E & Copy
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };

    /// Copy POLY_TRIANGULATION_RESULT_E.
    /// - Specify copying of this specific class.
    /// - Useful if copying is overwritten in inherited classes.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULT_E & CopyPolyTriangulationResultE
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };

  };


  /// Data structure for results of triangulation of polygons 
  ///   with bit set size 16.
  /// @tparam MAX_NUM_VERTICES Maximum number of polygon vertices.
  template <typename COS_TYPE_X, typename NTYPE>
  class POLY_TRIANGULATION_RESULT_E16:
    public POLY_TRIANGULATION_RESULT_E<16, COS_TYPE_X, NTYPE> {
  };


  // ***************************************************************
  // Member functions for POLY_TRIANGULATION_RESULT
  // ***************************************************************

  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
  PrintCosMinTriangulationAngle
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix << "cos_min_triangulation_angle: "
        << cos_min_triangulation_angle << "\n";
  }

  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::Print
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    PrintCosMinTriangulationAngle(out, line_prefix);
    out << line_prefix << "flag_zero: " << int(flag_zero) << "\n";
    out << line_prefix << "num_triangles: " << num_triangles << "\n";
    out << line_prefix << "tri_vertex_index: " << tri_vertex_index << "\n";
    out << line_prefix << "num_interior_tri_vertices: "
        << num_interior_tri_vertices << "\n";
  }


  // ***************************************************************
  // Member functions for POLY_TRIANGULATION_RESULT_E
  // ***************************************************************


  // Return location of first ear in polygon.
  // - Undefined if triangulation has no ear.
  template <int BIT_SET_SIZE, typename COS_TYPE_X, typename NTYPE>
  template <typename ITYPE1>
  bool POLY_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
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
  void POLY_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
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
  void POLY_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE_X,NTYPE>::
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


  // ***************************************************************
  // DEPRECATED CLASSES
  // ***************************************************************

  /// Extended data structure for results of triangulation of polygons.
  /// - *** DEPRECATED ***
  template <int BIT_SET_SIZE, typename COS_TYPE, typename NTYPE>
  class POLY_TRIANGULATION_RESULTX:
    public POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE, NTYPE> {

  public:

    /// Number of ears in the triangulation.
    NTYPE num_tri_ears;

    /// List of ears.
    /// - Note: This index has values 0,1,...,(num_poly_vert-1).
    /// - Index is with respect to the polygon.  It does not represent
    ///   the vertex index with respect to all the vertices in the mesh.
    NTYPE ear_list[BIT_SET_SIZE];

    /// Number of split polygon edges.
    NTYPE num_split_edges;

    /// Constructor
    POLY_TRIANGULATION_RESULTX() {
      num_tri_ears = 0;
      num_split_edges = 0;
    };

    /// Copy.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULTX & Copy
    (const TRI_RESULT_TYPE & right) {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
      return *this;
    };

    /// Copy.
    template <typename COS_TYPE2, typename NTYPE2>
    const POLY_TRIANGULATION_RESULTX & Copy
    (const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE2,NTYPE2> & right) {
      if (&right != this) {         // avoid self-assignment
        this->CopyPolyTriangulationResult(right);
      }
      return *this;
    };

    /// Copy and set number of split edges.
    template <const int NUM_SPLIT_EDGES, typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULTX & CopySetNumSplitEdges
    (const TRI_RESULT_TYPE & right) {
      Copy(right);
      num_split_edges = NUM_SPLIT_EDGES;
      return *this;
    };

    /// Copy and set number of split edges to 0.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULTX & Copy0Split
    (const TRI_RESULT_TYPE & right)
    {
      return(CopySetNumSplitEdges<0>(right));
    };

    /// Copy and set number of split edges to 1.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULTX & Copy1Split
    (const TRI_RESULT_TYPE & right)
    {
      return(CopySetNumSplitEdges<1>(right));
    };

    /// Copy and set number of split edges to 2.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULTX & Copy2Split
    (const TRI_RESULT_TYPE & right)
    {
      return(CopySetNumSplitEdges<2>(right));
    };

    /// Copy and set number of split edges to 3.
    template <typename TRI_RESULT_TYPE>
    const POLY_TRIANGULATION_RESULTX & Copy3Split
    (const TRI_RESULT_TYPE & right)
    {
      return(CopySetNumSplitEdges<3>(right));
    };

  };


  /// *** DEPRECATED ***
  /// Set equal to POLY_TRIANGULATION_RESULT. 
  template <int BIT_SET_SIZE, typename COS_TYPE, typename NTYPE>
  class POLY_TRIANGULATION_RESULTX2:
    public POLY_TRIANGULATION_RESULTX<BIT_SET_SIZE, COS_TYPE, NTYPE> {
  };

}

#endif
