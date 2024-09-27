/*!
 *  @file ijkdistance.tpp
 *  @brief ijk templates to compute distances and intersections between points, line segments and triangles.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2022 Rephael Wenger

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

#ifndef _IJKDISTANCE_
#define _IJKDISTANCE_

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"

// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // *****************************************************************
  //! @name Triangle properties.
  // *****************************************************************

  ///@{

  template <typename CTYPE0, typename CTYPE1>
  class TRIANGLE3D_PROPERTIES {
    
  protected:
    static const int NUM_VERTICES = 3;
    static const int DIM3 = 3;

    /// 3 pointers to vertex coordinates.
    CTYPE0 const * vertex_coord[NUM_VERTICES];

    /// @brief Triangle normal.
    /// - Vector orthogonal to triangle edges.
    /// - Always a unit vector.
    CTYPE1 normal[DIM3];

    /// @brief Array storing edge direction coordinates.
    CTYPE1 edge_dir_coord[DIM3*NUM_VERTICES];

    /// @brief 3 pointers to edge directions
    /// - edge_dir[i] is direction (vertex_coord[i+1]-vertex_cood[i]).
    /// - Could be (0,0,...,0) if two vertex coordinates are
    ///   almost identical.
    CTYPE1 * edge_dir[NUM_VERTICES];

    /// Edge lengths.
    CTYPE1 edge_length[NUM_VERTICES];

    /// @brief Array storing edge orthogonal direction coordinates.
    CTYPE1 edge_orth_dir_coord[DIM3*NUM_VERTICES];

    /// @brief 3 pointers to directions orthogonal to edge_dir[] and normal.
    /// - Could be (0,0,...,0) if edge_dir[] or normal[] is (0,0,...,0).
    /// - edge_orth_dir[j] points away from vertex (j+2)%DIM3.
    CTYPE1 * edge_orth_dir[NUM_VERTICES];

    /// Initialize.
    void Init();


  protected:

    /**
       @brief Print coordinates.
       - Mainly for debugging.
       - Extended version.  Define left/right delimiters and separator.
       @param c0 Left delimiter.
       @param c1 Separator.
       @param c2 Right delimiter.
    */
    template <typename OSTREAM_TYPE, typename CTYPEX>
    void _PrintCoordX
    (OSTREAM_TYPE & out, const CTYPEX * coord,
     const char c0, const char c1, const char c2) const;

    /**
       @brief Print coordinates.
       - Mainly for debugging.
       - Version using delimiters '(' and ')' and separator ','.
    */
    template <typename OSTREAM_TYPE, typename CTYPEX>
    void _PrintCoord
    (OSTREAM_TYPE & out, const CTYPEX * coord) const
    { _PrintCoordX(out, coord, '(', ',', ')'); }


    /**
       @brief Print coordinates.
       - Mainly for debugging.
       - Version adding prefix and suffix strings.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE, typename CTYPEX,
              typename STYPE0, typename STYPE1>
    void _PrintCoord
    (OSTREAM_TYPE & out, const STYPE0 & s0, 
     const CTYPEX * coord, const STYPE1 & s1) const;


  public:
    
    typedef CTYPE0 COORD_TYPE;
    typedef CTYPE1 NORMAL_TYPE;
    typedef CTYPE1 DIR_TYPE;


  public:
    TRIANGLE3D_PROPERTIES(){ Init(); }

    // Get functions.

    /// Return pointer to array vertex_coord[].
    const CTYPE0 * const * VertexCoord() const
    { return(vertex_coord); }

    /// Return pointer to coordinates of i'th vertex.
    template <typename ITYPE>
    const CTYPE0 * VertexCoord(const ITYPE i) const
    { return(vertex_coord[i]); }

    /// Return pointer to normal vector.
    const CTYPE1 * Normal() const
    { return(normal); }

    /// @brief Return index of edge (i0,i1).
    /// - If i1 = (i0+1)%DIM3, then return i0.
    /// - Else return i1.
    template <typename ITYPE>
    ITYPE EdgeIndex(const ITYPE i0, const ITYPE i1) const
    {
      if (((i0+1)%NUM_VERTICES) == i1) { return(i0); }
      else { return(i1); }
    }

    /// Return pointer to array edge_dir[].
    const CTYPE1 * const * EdgeDir() const
    { return(edge_dir); }

    /// Return pointer to i'th edge direction.
    template <typename ITYPE>
    const CTYPE1 * EdgeDir(const ITYPE i) const
    { return(edge_dir[i]); }

    /// @brief Return pointer to direction of edge between i0 and i1.
    /// - Note: Edge may point from i0 to i1 or from i1 to i0.
    template <typename ITYPE>
    const CTYPE1 * EdgeDir(const ITYPE i0, const ITYPE i1) const
    { return(EdgeDir(EdgeIndex(i0,i1))); }

    /// @brief Return pointer to array edge_orth_dir[].
    CTYPE1 * const * EdgeOrthDir() const
    { return(edge_orth_dir); }

    /// @brief Return pointer to i'th direction orthogonal to normal and edge_dir[i].
    template <typename ITYPE>
    const CTYPE1 * EdgeOrthDir(const ITYPE i) const
    { return(edge_orth_dir[i]); }

    /// @brief Return pointer to array edge_length[].
    const CTYPE1 * EdgeLength() const
    { return(edge_length); }

    /// @brief Return length of i'th edge.
    template <typename ITYPE>
    CTYPE1 EdgeLength(const ITYPE i) const
    { return(edge_length[i]); }


    // Set functions

    /// Set pointers to i'th vertex coord.
    template <typename ITYPE>
    void SetVertexCoord(const ITYPE i, const CTYPE0 * coord)
    { vertex_coord[i] = coord; }

    /**
     *  @brief Set pointers to i'th vertex coord from vertex list.
     *  - Polygon has at most three vertices (typically exactly three.)
     *  @param vertex_coord[] Array of vertex coordinates.
     *  @param vertex_list[] Array of vertex indices.
     *  @param num_poly_vert Number of polygon vertices.
     *    - @pre num_poly_vert <= 3.
     */
    template <typename ITYPE, typename IVTYPE, typename NTYPE0>
    void SetVertexCoordFromVertexList
    (const ITYPE i, const CTYPE0 vertex_coord_list[],
     const IVTYPE * vertex_list, const NTYPE0 num_poly_vert);

    /// Copy coord into vertex normal.
    template <typename _CTYPE>
    void CopyNormal(const _CTYPE * coord)
    { std::copy(coord, coord+DIM3, normal); }

    /// Copy coord into i'th edge direction.
    template <typename ITYPE, typename _CTYPE>
    void CopyEdgeDir(const ITYPE i, const _CTYPE * coord)
    { std::copy(coord, coord+DIM3, edge_dir[i]); }


    // Compute functions

    /// Compute edge directions.
    /// @pre All values of vertex_coord[] are set.
    void ComputeEdgeDir();

    /// Compute normal.
    /// - Also computes edge_orth_dir[] for each edge.
    /// - Sets edge_orth_dir[j] to point away from vertex (j+2)%DIM3.
    /// @pre All values of vertex_coord[] are set.
    /// @pre All edge directions are set.
    void ComputeNormal();

    /**
     *  Compute signed distance of three triangle vertices to plane.
     *  - Return 0 if triangle intersects plane.
     *  - Return signed distance to closest point if triangle 
     *    does not intersect plan.
     *  @param[out] signed_distance[i] Signed distance of vertex_coord[i] 
     *     to plane.
     *  @param[out] jindex If triangle intersects plane, index (0, 1 or 2)  
     *    of point separated from two others.
     *    - If triangle does not intersect plane, closest point to plane.
     */
    template <typename CTYPEQ, typename CTYPEN, typename DIST_TYPE,
              typename ITYPE>
    DIST_TYPE ComputeSignedDistanceToPlane
    (const CTYPEQ q0, const CTYPEN qnormal,
     DIST_TYPE signed_distance[3], ITYPE & jindex) const;


    // Print routines (mainly for debugging.)

    /**
       @brief Print vertex coordinates.
       - Mainly for debugging.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void PrintVertexCoord
    (OSTREAM_TYPE & out, const STYPE0 & s0,
     const ITYPE2 i, const STYPE1 & s1) const
    { _PrintCoord(out, s0, vertex_coord[i], s1); }

    /**
       @brief Print normal direction.
       - Mainly for debugging.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE,
              typename STYPE0, typename STYPE1>
    void PrintNormal
    (OSTREAM_TYPE & out, const STYPE0 & s0, const STYPE1 & s1) const
    { _PrintCoord(out, s0, normal, s1); }

    /**
       @brief Print i'th edge direction.
       - Mainly for debugging.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void PrintEdgeDir
    (OSTREAM_TYPE & out, const STYPE0 & s0,
     const ITYPE2 i, const STYPE1 & s1) const
    { _PrintCoord(out, s0, edge_dir[i], s1); }

    /**
       @brief Print i'th direction orthogonal to normal and EdgeDir(i).
       - Mainly for debugging.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void PrintEdgeOrthDir
    (OSTREAM_TYPE & out, const STYPE0 & s0,
     const ITYPE2 i, const STYPE1 & s1) const
    { _PrintCoord(out, s0, edge_orth_dir[i], s1); }

    /**
       @brief Print all three edge lengths.
    */
    template <typename OSTREAM_TYPE, typename STYPE0, typename STYPE1>
    void PrintEdgeLengths(OSTREAM_TYPE & out, 
                          const STYPE0 & s0, const STYPE1 & s1) const;

    /**
       @brief Print vertex, edge_dir, edge_length, edge_orth_dir[] coordinates.
       @param line_prefix Print line_prefix at beginning of each line.
    */
    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;
    

    // Check functions.

    /// Check that all vertex coordinates are set, i.e., that each vertex_coord[i] is not NULL.
    bool CheckVertexCoord(IJK::ERROR & error) const;

  };
  
  //@}


  // *****************************************************************
  //! @name Distance type.
  // *****************************************************************

  ///@{

  /**
   *  @brief Type of distance.
   *  - V2E_DIST Vertex of first poly to edge interior of second poly.
   *  - E2V_DIST Edge interior of first poly to vertex of second poly.
   *  - V2T_DIST Vertex of first poly to triangle interior of second poly.
   *  - T2V_DIST Triangle interior of first poly to vertex of second poly.
   *  - E2T_DIST Edge interior of first poly to triangle interior 
   *    of second poly.
   *  - T2E_DIST Triangle interior of first poly to edge interior 
   *    of second poly.
   *  - If distance type is E2T_DIST or T2E_DIST, then edge should
   *   intersect triangle and distance should be 0.
   */
  typedef enum { V2V_DIST, V2E_DIST, E2V_DIST, E2E_DIST, 
                 V2T_DIST, T2V_DIST, E2T_DIST, T2E_DIST,
                 DIST_PAIR_TYPE_NOT_SET }
    DISTANCE_PAIR_TYPE;

  template <typename DIST_TYPE, typename ITYPEF>
  class DISTANCE_DATA {

  protected:

    /// Distance to iclosest_poly.
    IJK::SET_VALUE<DIST_TYPE> distance;

    /// Type of distance pair.
    DISTANCE_PAIR_TYPE distance_pair_type;

    /**
     *  @brief Index of closest poly faces.
     *  - Can represent vertex index or edge index or be undefined.
     */
    ITYPEF iclosest_face_index[2];

    /// Initialize.
    void Init();

    /// Set function.
    template <typename DIST_TYPE2, typename ITYPE2>
    void _SetDistance
    (const DISTANCE_PAIR_TYPE distance_pair_type,
     const DIST_TYPE2 dist, const ITYPE2 iface0, const ITYPE2 iface1);


  public:
    DISTANCE_DATA() { Init(); }

    // get functions

    /// Return distance to closest poly.
    DIST_TYPE Distance() const
    { return(distance.Value()); }

    /// Return true if distance is set.
    bool IsDistanceSet() const
    { return(distance.IsSet()); }

    /// Return distance pair type.
    DISTANCE_PAIR_TYPE DistancePairType() const
    { return(distance_pair_type); }

    /**
     *  @brief Return face index.
     *  - Note: This index is with respect to the list of poly
     *    vertices and edges.
     *  - If defined, the index is between 
     *    0 and (number of poly vertices minus 1) or between
     *    (number of poly edges minus 1.)
     */
    template <typename ITYPE2>
    ITYPEF IndexOfFace(const ITYPE2 j) const
    { return(iclosest_face_index[j]); }

    // Set functions


    /// Set distance from first poly to second poly.
    /// - Does not set distance type or face information.
    template <typename DIST_TYPE2>
    void SetDistance(const DIST_TYPE2 dist)
    { _SetDistance(DIST_PAIR_TYPE_NOT_SET, dist, 0, 0); }

    /// Set distance to closest of type vertex to vertex.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceV2V
    (const DIST_TYPE2 dist, const ITYPE2 iv0, const ITYPE2 iv1)
    { _SetDistance(V2V_DIST, dist, iv0, iv1); }

    /// Set distance to closest of type vertex to edge.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceV2E
    (const DIST_TYPE2 dist, const ITYPE2 iv0, const ITYPE2 ie1)
    { _SetDistance(V2E_DIST, dist, iv0, ie1); }

    /// Set distance to closest of type edge to vertex.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceE2V
    (const DIST_TYPE2 dist, const ITYPE2 ie0, const ITYPE2 iv1)
    { _SetDistance(E2V_DIST, dist, ie0, iv1); }

    /// Set distance to closest of type edge to edge.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceE2E
    (const DIST_TYPE2 dist, const ITYPE2 ie0, const ITYPE2 ie1)
    { _SetDistance(E2E_DIST, dist, ie0, ie1); }

    /// Set distance to closest of type vertex to triangle.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceV2T
    (const DIST_TYPE2 dist, const ITYPE2 iv0)
    { _SetDistance(V2T_DIST, dist, iv0, 0); }

    /// Set distance to closest of type triangle to vertex.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceT2V
    (const DIST_TYPE2 dist, const ITYPE2 iv1)
    { _SetDistance(V2T_DIST, dist, 0, iv1); }

    /// Set distance to closest of type edge to triangle.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceE2T
    (const DIST_TYPE2 dist, const ITYPE2 ie0)
    { _SetDistance(E2T_DIST, dist, ie0, 0); }


    /// Set distance to closest of type triangle to edge.
    template <typename DIST_TYPE2, typename ITYPE2>
    void SetDistanceT2E
    (const DIST_TYPE2 dist, const ITYPE2 ie1)
    { _SetDistance(T2E_DIST, dist, 0, ie1); }

    /**
     *  @brief Swap first and second poly data.
     *  - Swap iclosest_facet_index[0] and iclosest_facet_index[1].
     *  - If distance_pair_type is V2E or V2T or E2T, 
     *    change to E2V or T2V or T2E.
     *  - If distance_pair_type is E2V or T2V or T2E, 
     *    change to E2V or T2V or E2T.
     */
    void SwapPolyData();

    /// Copy right to this.
    template <typename DIST_DATA_TYPE2>
    void Copy(const DIST_DATA_TYPE2 & right);


    // Print routines (mainly for debugging.)

    /**
       @brief Print distance type and facet info.
       - Mainly for debugging.
    */
    template <typename OSTREAM_TYPE>
    void PrintDistancePairType(OSTREAM_TYPE & out) const;

    /**
       @brief Print distance type and facet info.
       - Mainly for debugging.
    */
    template <typename OSTREAM_TYPE, typename STYPE0, typename STYPE1>
    void PrintDistancePairType
    (OSTREAM_TYPE & out, const STYPE0 & s0, const STYPE1 & s1) const;

    /**
       @brief Print distance, distance pair type and facet info.
       - Mainly for debugging.
       @param line_prefix Print line_prefix at beginning of each line.
    */
    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

  };


  ///@}


  // *****************************************************************
  //! @name Compute orthogonal vector.
  // *****************************************************************

  ///@{

  /**
   *  @brief Compute unit vector orthogonal to two unit vectors.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *     @pre dimension >= 3.
   *  @param dir0[] Direction 0. 
   *    @pre dir0[] is a unit vector or a zero vector.
   *  @param dir1[] Direction 1. 
   *    - dir1[] could be parallel to dir0[].
   *    @pre dir1[] is a unit vector or a zero vector.
   *  @param[out] v_orth[] Component of \a v orthogonal 
   *    to \a dir0[] and \a dir1[].
   *    - v_orth[] is guaranteed to be non-zero.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPEV>
  void compute_orthogonalII_vector
  (const DTYPE dimension, 
   const CTYPE0 dir0[], const CTYPE1 dir1[],
   CTYPEV v_orth[])
  {
    const DTYPE THREE(3);

    // Vector orthogonal to dir0[].
    IJK::ARRAY<CTYPEV> orth_dir0(dimension);
    IJK::ARRAY<CTYPEV> v_init(dimension);
    CTYPEV mag;

    compute_orthogonal_vector_component
      (dimension, dir1, dir0, orth_dir0.Ptr());
    normalize_vector_robust
      (dimension, orth_dir0.Ptr(), 0, orth_dir0.Ptr(), mag);

    const DTYPE d0 = get_index_of_max_abs_coord(dimension, dir0);
    const DTYPE d1 = 
      get_index_of_max_abs_coord(dimension, orth_dir0.Ptr());

    set_coord(dimension, 1, v_init.Ptr());
    v_init[d0] = 0;
    v_init[d1] = 0;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << endl << "-----" << endl;
    IJK::print_coord3D(cerr, "  dir0: ", dir0, "\n");
    cerr << "  max component: " << d0 << endl;
    IJK::print_coord3D(cerr, "  dir1: ", dir1, "\n");
    IJK::print_coord3D(cerr, "  initial non-planar vector: ", v_orth, "\n");
    IJK::print_coord3D(cerr, "  orth_dir0: ", orth_dir0.Ptr(), "\n");
    cerr << "  max component: " << d1 << endl;
    */

    // v_orth[] is far from the plane spanned by dir0 and dir1.
    compute_orthogonalII_vector_component
      (dimension, v_init.Ptr(), dir0, orth_dir0.Ptr(), v_orth);

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "  v_orth: ", v_orth, "\n");
    */

    bool flag_zero;
    normalize_vector(dimension, v_orth, 0, v_orth, mag, flag_zero);

    if (flag_zero) {
      IJK::PROCEDURE_ERROR error("compute_orthogonalII_vector");
      // This should never happen, but...
      if (dimension < THREE) {
        error.AddMessage
          ("Programming error. Illegal dimension ", dimension, ".");
        error.AddMessage
          ("  Dimension must be at least ", THREE, ".");
        throw error;
      }
      else {
        error.AddMessage
          ("Programming error. Unable to compute vector normal to dir0[] and dir1[].");
        throw error;
      }
    }

  }

  ///@}

  // *****************************************************************
  //! @name Compute triangle normal.
  // *****************************************************************

  ///@{

  /**
   *  @brief Compute unit vector orthogonal to triangle edge directions.
   *  - If dimension is 3, orient normal in same direction as u01 x u12.
   *  - If dimension is greater than 3, oriented normal has direction
   *      matching orientation of u01' x u12' where u01' and u12' are
   *      first three coordinates of u01 and u12.
   *  - If u01, u12 and u20 are zero vectors,
   *    returns a unit vector in an arbitrary direction.
   *  - If three coordinates are (almost) collinear, 
   *    returns an arbitrary unit vector that is perpendicular
   *    to the collinear points.
   *  @pre Dimension is at least 3.
   *  @param u01[] First edge direction.
   *  @param u12[] Secondt edge direction.
   *  @param u20[] Third edge direction.
   *  @param edge_length[i] Length of edge i.
   *  @param[out] normal[] Unit vector normal to the triangle.
   *    - Guaranteed to be non-zero.
   */
  template <typename DTYPE, typename CTYPE, typename LTYPE, 
            typename CTYPEN>
  void compute_triangle_normal
  (const DTYPE dimension,
   const CTYPE u01[], const CTYPE u12[], const CTYPE u20[],
   const LTYPE edge_length[3],
   CTYPEN normal[])
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    const int DIM3(3);
    const CTYPE * u[NUM_VERT_PER_TRIANGLE] = { u01, u12, u20 };
    IJK::PROCEDURE_ERROR error("compute_triangle_normal");

    if (dimension < DIM3) {
      error.AddMessage
        ("Programming error. Dimension must be at least ", DIM3, ".");
      error.AddMessage("  dimension = ", dimension, "");
      throw error;
    }

    // Choose longest edges to compute normal.
    DTYPE imin_mag = 0;
    if (edge_length[imin_mag] > edge_length[1]) { imin_mag = 1; }
    if (edge_length[imin_mag] > edge_length[2]) { imin_mag = 2; }

    const DTYPE d0 = (imin_mag+1)%NUM_VERT_PER_TRIANGLE;
    const DTYPE d1 = (d0+1)%NUM_VERT_PER_TRIANGLE;

    compute_orthogonalII_vector(dimension, u[d0], u[d1], normal);

    // If dimension > 3, this function will use compute the determinant
    //   of the 3x3 matrix consisting of the first 3 coordinates
    //   of each vector.
    CTYPEN det;
    determinant_3x3(u[d0], u[d1], normal, det);

    if (det < 0) {
      // Reverse normal orientation.
      multiply_coord(dimension, -1, normal, normal);
    }

  }

  ///@}


  // *****************************************************************
  //! @name Orient vector
  // *****************************************************************

  /// Orient vector u1 away from vector u0[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void orient_vector_away_from
  (const DTYPE dimension, const CTYPE0 u0[], CTYPE1 u1[])
  {
    CTYPE1 product;
    IJK::compute_inner_product_3D(u0, u1, product);
    if (product > 0)  {
      // Reverse u1 orientation.
      multiply_coord(dimension, -1, u1, u1);
    }
  }


  // *****************************************************************
  //! @name Compute distance between lines.
  // *****************************************************************

  ///@{

  /**
   *  @brief Compute distance between two lines.
   *  - Lines are p0+t*u0 and q0+s*u1
   *  @param[out] v_orth[] Unit vector orthogonal to (p0,p1) and (q0,q1).
   *    - v_orth[] points in same direction as (q0-p0).
   */
  template <typename DTYPE, 
            typename CTYPEP0, typename CTYPEP1,
            typename CTYPEU0, typename CTYPEU1,
            typename CTYPEV, typename DIST_TYPE>
  void compute_line_to_line_distance
  (const DTYPE dimension, 
   const CTYPEP0 p0[], const CTYPEU0 u0[],
   const CTYPEP1 p1[], const CTYPEU1 u1[],
   CTYPEV v_orth[], DIST_TYPE & dist)
  {
    IJK::ARRAY<CTYPEV> u(dimension);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  In " << __func__ << endl;
    IJK::print_coord3D(cerr, "    p0: ", p0, "\n");
    IJK::print_coord3D(cerr, "    u0: ", u0, "\n");
    IJK::print_coord3D(cerr, "    p1: ", p1, "\n");
    IJK::print_coord3D(cerr, "    u1: ", u1, "\n");
    */

    compute_orthogonalII_vector
      (dimension, u0, u1, v_orth);

    IJK::subtract_coord(dimension, p1, p0, u.Ptr());

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "    (p1-p0): ", u.Ptr(), "\n");
    IJK::print_coord3D(cerr, "    v_orth:  ", v_orth, "\n");
    */


    IJK::compute_inner_product
      (dimension, u.Ptr(), v_orth, dist);

    if (dist < 0) {
      IJK::multiply_coord(dimension, -1, v_orth, v_orth);
      dist = -dist;
    }
  }

  ///@}


  // *****************************************************************
  //! @name Distance of plane to triangle or line segment.
  // *****************************************************************

  ///@{

  /**
   *  @brief Return signed distance of triangle to hyperplane.
   *  - Return 0 if triangle intersects plane.
   *  @param coord[] Pointers to 3 triangle coordinates.
   *  @param q0_coord Coord of point on plane.
   *  @param normal[] Unit normal.
   *  @param[out] signed_distance[i] Signed distance of coord[i] to plane.
   *  @param[out] jindex If triangle intersects plane, index (0, 1 or 2)  
   *    of point separated from two others.
   *    - If triangle does not intersect plane, closest point to plane.
   */
  template <typename DTYPE,
            typename CTYPEQ, typename CTYPEN, typename CTYPE_PTR,
            typename DIST_TYPE, typename ITYPE>
  DIST_TYPE compute_signed_distance_triangle_to_hyperplane
  (const DTYPE dimension, const CTYPE_PTR coord[], 
   const CTYPEQ q0_coord[], const CTYPEN normal[],
   DIST_TYPE signed_distance[],
   ITYPE & jindex)
  {
    const int NUM_VERT_PER_TRIANGLE(3);

    IJK::compute_signed_distance_to_hyperplane
      (dimension, coord[0], q0_coord, normal, signed_distance[0]);
    IJK::compute_signed_distance_to_hyperplane
      (dimension, coord[1], q0_coord, normal, signed_distance[1]);
    IJK::compute_signed_distance_to_hyperplane
      (dimension, coord[2], q0_coord, normal, signed_distance[2]);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  In " << __func__ << endl;
    */

    // Check if triangle intersects plane.
    for (int d0 = 0; d0 < NUM_VERT_PER_TRIANGLE; d0++) {
      const int d1 = (d0+1)%NUM_VERT_PER_TRIANGLE;
      const int d2 = (d1+1)%NUM_VERT_PER_TRIANGLE;
      if (signed_distance[d0] > 0) {
        if ((signed_distance[d1] <= 0) &&
            (signed_distance[d2] <= 0)) {

          jindex = d0;
          return(0);
        }
      }
      else if (signed_distance[d0] < 0) {
        if ((signed_distance[d1] >= 0) &&
            (signed_distance[d2] >= 0)) {

          jindex = d0;
          return(0);
        }
      }
    }

     // Either all signed_distances are positive, or all are negative
    //   or all are 0.
    jindex = get_index_of_min_abs_coord(NUM_VERT_PER_TRIANGLE, signed_distance);
    return(signed_distance[jindex]);
  }


  /**
   *  @brief Return signed distance of line segment to hyperplane.
   *  - Return 0 if line segment intersects hyperplane.
   *  @param coordp0[] Coordinates of line segment endpoint 0.
   *  @param coordp1[] Coordinates of line segment endpoint 1.
   *  @param q0_coord Coord of point on yperplane.
   *  @param normal[] Unit normal.
   *  @param[out] signed_distance[i] Signed distance of endpoint i to plane.
   *  @param[out] jindex If line segment does not intersect plane
   *      closest point to plane.
   *    - Otherwise, 0.
   */
  template <typename DTYPE,
            typename CTYPEQ, typename CTYPEN, 
            typename CTYPEP0, typename CTYPEP1,
            typename DIST_TYPE, typename ITYPE>
  DIST_TYPE compute_signed_distance_line_segment_to_hyperplane
  (const DTYPE dimension,
   const CTYPEP0 coordp0[], const CTYPEP1 coordp1[],
   const CTYPEQ q0_coord[], const CTYPEN normal[],
   DIST_TYPE signed_distance[2], ITYPE & jindex)
  {
    const int TWO(2);

    // Initialize jindex.
    jindex = 0;

    IJK::compute_signed_distance_to_hyperplane
      (dimension, coordp0, q0_coord, normal, signed_distance[0]);
    IJK::compute_signed_distance_to_hyperplane
      (dimension, coordp1, q0_coord, normal, signed_distance[1]);

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "        coordp0: ", coordp0, "\n");
    IJK::print_coord3D(cerr, "        coordp1: ", coordp1, "\n");
    IJK::print_coord3D(cerr, "        q0_coord: ", q0_coord, "\n");
    cerr << "        signed_distance[0]: " << signed_distance[0] << endl;
    cerr << "        signed_distance[1]: " << signed_distance[1] << endl;
    */

    if (signed_distance[0]*signed_distance[1] <= 0)
      { return(0); }
    else {
      jindex = get_index_of_min_abs_coord(TWO, signed_distance);
      return(signed_distance[jindex]);
    }
  }


  /**
   *  @brief Return signed distance of line segment to hyperplane.
   *  - Version that does not return signed_distance[] at each endoint
   *    and does not return jindex.
   *  - Return 0 if line segment intersects plane.
   *  @tparam DIST_TYPE Distance type. Must be given in function call.
   *  @param coordp0[] Coordinates of first line segment endpoint.
   *  @param coordp1[] Coordinates of second line segment endpoint.
   *  @param q0_coord Coord of point on plane.
   *  @param normal[] Unit normal.
   */
  template <typename DIST_TYPE, typename DTYPE, 
            typename CTYPEQ, typename CTYPEN, 
            typename CTYPEP0, typename CTYPEP1>
  DIST_TYPE compute_signed_distance_line_segment_to_hyperplane
  (const DTYPE dimension, 
   const CTYPEP0 coordp0[], const CTYPEP1 coordp1[],
   const CTYPEQ q0_coord[], const CTYPEN normal[])
  {
    const int TWO(2);
    DIST_TYPE signed_distance[TWO];
    int jindex;
    const DIST_TYPE dist = 
      compute_signed_distance_line_segment_to_hyperplane
      (dimension, coordp0, coordp1, q0_coord, normal, 
       signed_distance, jindex);

    return(dist);
  }

  ///@}


  // *****************************************************************
  //! @name Intersect line segment and hyperplane.
  // *****************************************************************

  /// Return intersection of line segment and hyperplane.
  /// - Returns closest endpoint of (coordp0[], coordp1[]) to hyperplane,
  ///   if dist2plane0 and dist2plane1 are both non-negative
  ///   or both non-positive.
  /// @param dist2plane0 Signed distance of coord0[] to hyperplane.
  /// @param dist2plane1 Signed distance of coord1[] to hyperplane.
  template <typename DTYPE, typename CTYPEP,
            typename DIST_TYPE, typename CTYPEQ, typename CTYPEN,
            typename CTYPEI>
  void intersect_line_segment_and_hyperplane
  (const DTYPE dimension, const CTYPEP coordp0[], const CTYPEP coordp1[],
   const DIST_TYPE dist2plane0, const DIST_TYPE dist2plane1,
   const CTYPEQ coordq[], const CTYPEN normal[], 
   CTYPEI intersection_coord[])
  {
    const DTYPE TWO(2);
    const CTYPEP * const coordp[TWO] = { coordp0, coordp1 };
    const DIST_TYPE dist2plane[TWO] = { dist2plane0, dist2plane1 };
    DIST_TYPE interpolation_coef[TWO];

    if ((dist2plane0 >= 0 && dist2plane1 >= 0) || 
        (dist2plane0 <= 0 && dist2plane1 <= 0)) {
      const DTYPE imin = IJK::get_index_of_min_abs_coord(TWO, dist2plane);
      IJK::copy_coord(dimension, coordp[imin], intersection_coord);
      return;
    }

    const DTYPE imax = IJK::get_index_of_max_abs_coord(TWO, dist2plane);
    const DTYPE imin = (imax+1)%TWO;
    const DIST_TYPE max_dist = std::abs(dist2plane[imax]);
    const DIST_TYPE min_dist = std::abs(dist2plane[imin]);

    if (max_dist == 0) {
      // Should never happen but just in case.
      IJK::copy_coord(dimension, coordp[imax], intersection_coord);
    }

    interpolation_coef[imax] = min_dist/max_dist;
    interpolation_coef[imin] = 1;

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "      coordq: ", coordq, "\n");
    IJK::print_coord3D(cerr, "      normal: ", normal, "\n");
    cerr << "     imax: " << imax << "  max_dist: " << max_dist << endl;
    cerr << "     imin: " << imin << "  min_dist: " << min_dist << endl;
    */

    // Normalize interpolation_coord[] to sum to 1.
    const DIST_TYPE sum_weights = interpolation_coef[0] + interpolation_coef[1];
    interpolation_coef[0] = interpolation_coef[0]/sum_weights;
    interpolation_coef[1] = interpolation_coef[1]/sum_weights;

    IJK::linear_interpolate_coord
      (dimension, interpolation_coef[0], coordp0, coordp1, 
       intersection_coord);

    return;
  }


  // *****************************************************************
  //! @name Do objects intersect or cross.
  // *****************************************************************

  ///@{

  /**
   *  @brief Return true if hyperplane separates points, i.e. if line segment (coordp0,coordp1) intersects plane.
   *  @tparam DIST_TYPE Distance type. Must be given in function call.
   *  @param q0_coord Coord of point on plane.
   *  @param normal[] Unit normal.
   *  @param coordp0[] Coordinates of first line segment endpoint.
   *  @param coordp1[] Coordinates of second line segment endpoint.
   */
  template <typename DIST_TYPE,
            typename DTYPE,
            typename CTYPEQ, typename CTYPEN, 
            typename CTYPEP0, typename CTYPEP1>
  bool does_hyperplane_separate_points
  (const DTYPE dimension,
   const CTYPEQ q0_coord[], const CTYPEN normal[],
   const CTYPEP0 coordp0[], const CTYPEP1 coordp1[])
  {
    // Note: Point coordinates, coordp0, coordp1, precede
    //   hyperplane coordinates q_0_coord, normal in function call.
    if (compute_signed_distance_line_segment_to_hyperplane<DIST_TYPE>
        (dimension, coordp0, coordp1, q0_coord, normal) == 0)

      { return(true); }
    else
      { return(false); }
  }


  /**
   *  @brief Return true if line segment projections cross.
   *  @tparam DIST_TYPE Distance type. Must be given in function call.
   *  @param coordp0[] Coordinates of line segment p endpoint 0.
   *  @param coordp1[] Coordinates of line segment p endpoint 1.
   *  @param dirp[] Unit vector in direction (coordp1[]-coordp0[])
   *              or zero vector.
   *  @param pdist Distance from coordp0[] to coordp1[].
   *  @param coordq0[] Coordinates of line segment q endpoint 0.
   *  @param coordq1[] Coordinates of line segment q endpoint 1.
   *  @param u_q[] Unit vector in direction (coordq1[]-coordq0[]).
   *              or zero vector.
   *  @param qdist Distance from coordq0[] to coordq1[].
   *  @param projection_dir[] Projection direction.
   *    @pre projection_dir is a unit vector.
   */
  template <typename DIST_TYPE, typename DTYPE,
            typename CTYPEP0, typename CTYPEP1, typename DIRP_TYPE,
            typename CTYPEQ0, typename CTYPEQ1,typename DIRQ_TYPE,
            typename CDIST_TYPE, typename DIR_TYPE>
  bool do_projected_line_segments_cross
  (const DTYPE dimension,
   const CTYPEP0 coordp0[], const CTYPEP1 coordp1[],
   const DIRP_TYPE dirp[], const CDIST_TYPE pdist,
   const CTYPEQ0 coordq0[], const CTYPEQ1 coordq1[],
   const DIRQ_TYPE dirq[], const CDIST_TYPE qdist,
   const DIR_TYPE projection_dir[])
  {
    IJK::ARRAY<DIST_TYPE> orthp(dimension);
    IJK::ARRAY<DIST_TYPE> orthq(dimension);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "      In " << __func__ << endl;
    IJK::print_coord3D(cerr, "      coordp0: ", coordp0, "\n");
    IJK::print_coord3D(cerr, "      coordp1: ", coordp1, "\n");
    IJK::print_coord3D(cerr, "      dirp:    ", dirp, "\n");
    IJK::print_coord3D(cerr, "      coordq0: ", coordq0, "\n");
    IJK::print_coord3D(cerr, "      coordq1: ", coordq1, "\n");
    IJK::print_coord3D(cerr, "      dirq:    ", dirq, "\n");
    */

    if (pdist == 0 || qdist == 0)
      { return(false); }

    // dirp[] and dirq[] are unit vectors orthogonal to projection_dir.

    compute_orthogonalII_vector
      (dimension, projection_dir, dirp, orthp.Ptr());
    compute_orthogonalII_vector
      (dimension, projection_dir, dirq, orthq.Ptr());

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "      orthp:    ", orthp.Ptr(), "\n");
    IJK::print_coord3D(cerr, "      orthq:    ", orthq.Ptr(), "\n");
    */

    if (does_hyperplane_separate_points<DIST_TYPE>
        (dimension, coordp0, orthp.Ptr(), coordq0, coordq1) &&
        does_hyperplane_separate_points<DIST_TYPE>
        (dimension, coordq0, orthq.Ptr(), coordp0, coordp1))
      { return(true); }
    else
      { return(false); }
  }

  /** 
   *  @brief Compute abs((u0*projection_dir)/(u1*projection_dir)).
   *  - Return false if (u0*projection_dir) and rsiqn*(u1*projection_dir)
   *    have different signs.
   *  - Return false if (abs(u0*projection_dir) >= abs(u1*projection_dir)).
   *  @param[out] ratio = abs((u0*projection_dir)/(u1*projection_dir)).
   *  - Set to 0 if function returns false.
   */
  template <typename DTYPE, typename CTYPEU0, typename CTYPEU1,
            typename CTYPEP, typename SIGN_TYPE, typename RTYPE>
  bool compute_projection_length_ratio
  (const DTYPE dimension, const CTYPEU0 u0[], const CTYPEU1 u1[],
   const CTYPEP projection_dir[], const SIGN_TYPE rsign, RTYPE & ratio)
  {
    RTYPE ax, ay;

    // Initialize.
    ratio = 0;

    IJK::compute_inner_product(dimension, u0, projection_dir, ax);
    IJK::compute_inner_product(dimension, u1, projection_dir, ay);
    ay = rsign*ay;

    if ((ax < 0) && (ay > 0)) { return(false); }
    if ((ax > 0) && (ay < 0)) { return(false); }
    ax = std::abs(ax);
    ay = std::abs(ay);

    if (ay == 0 || (ax >= ay)) { return(false); }
    ratio = ax/ay;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "      ax: " << ax << "  ay: " << ay << "  ratio: " << ratio
         << endl;
    */

    return(true);
  }


  /** 
   *  @brief Return true if edge ieA segment intersects triangle interior.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE, typename ITYPEF>
  bool does_edge_intersect_triangle_interior
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   const DIST_TYPE dist2planeB[3],
   const ITYPEF ieA)
  {
    typedef typename TRI_TYPEB::COORD_TYPE COORDB_TYPE;

    const ITYPEF DIM3(3);
    const ITYPEF NUM_VERT_PER_TRIANGLE(3);

    const ITYPEF ij0 = ieA;
    const ITYPEF ij1 = (ieA + 1)%DIM3;
    const ITYPEF ij2 = (ieA + 2)%DIM3;
    DIST_TYPE qcoord[DIM3];
    DIST_TYPE vq[DIM3], vB20[DIM3], vB21[DIM3];
    const COORDB_TYPE * const * coordB = triB_properties.VertexCoord();
    
    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  In " << __func__ << endl;
    cerr << "    triA_properties:" << endl;
    triA_properties.Print(cerr, "      ");
    cerr << "    triB_properties:" << endl;
    triB_properties.Print(cerr, "      ");
    cerr << "    coordA[" << ij0 << "]:";
    IJK::print_coord3D
      (cerr, " ", triA_properties.VertexCoord(ij0), "\n");
    cerr << "    coordA[" << ij1 << "]:";
    IJK::print_coord3D
      (cerr, " ", triA_properties.VertexCoord(ij1), "\n");
    */

    intersect_line_segment_and_hyperplane
      (DIM3, triA_properties.VertexCoord(ij0), 
       triA_properties.VertexCoord(ij1),
       dist2planeB[ij0], dist2planeB[ij1], 
       triB_properties.VertexCoord(2), triB_properties.Normal(), qcoord);

    // Find b0, b1 such that:
    //   - b0*coordB[0] + b1*coordB[1] + (1-b0-b1)*coordB[2] = qcoord.
    // Equivalently:
    //   - b0*(coordB[0]-coordB[2]) + b1*(coordB[1]-coordB[2]) = 
    //       (qcoord - coordB[2]).

    IJK::subtract_coord_3D(qcoord, coordB[2], vq);
    IJK::subtract_coord_3D(coordB[0], coordB[2], vB20);
    IJK::subtract_coord_3D(coordB[1], coordB[2], vB21);

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "    qcoord: ", qcoord, "\n");
    IJK::print_coord3D(cerr, "    vq: ", vq, "\n");
    */
    

    DIST_TYPE b0, b1;

    // b1 = vq*edge_orth_dir[2]/vB21*edge_orth_dir[2],
    if (!compute_projection_length_ratio
        (DIM3, vq, vB21, triB_properties.EdgeOrthDir(2), 1, b1)) 
      { return(false); }

    // b0 = vq*edge_orth_dir[1]/vB20*edge_orth_dir[1],
    if (!compute_projection_length_ratio
        (DIM3, vq, triB_properties.EdgeDir(2), vB21, 1, b0))
      { return(false); }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "      b0: " << b0 << "  b1: " << b1 << "  b2: " << b2 << endl;
    */

    if ((b0 <= 0) || (b0 >= 1)) { return(false); }
    if ((b1 <= 0) || (b1 >= 1)) { return(false); }

    const DIST_TYPE b2 = 1-b0-b1;
    if ((b2 <= 0) || (b2 >= 1)) { return(false); }

    return(true);
  }

  ///@}


  // *****************************************************************
  //! @name Compute distance between triangle faces.
  // *****************************************************************

  ///@{

  /**
   *  @brief Compute minimum distance between vertices of two triangles.
   */
  template <typename DTYPE, typename CTYPEA, typename CTYPEB,
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_vertices_of_two_triangles
  (const DTYPE dimension,
   const CTYPEA * const coordA[3], const CTYPEB * const coordB[3],
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distance_data)
  {
    const DTYPE NUM_VERT_PER_TRIANGLE(3);
    
    ITYPEF ivA, ivB;
    DIST_TYPE min_dist;

    // Initialize.
    IJK::compute_distance
      (dimension, coordA[0], coordB[0], min_dist);
    ivA = 0;
    ivB = 0;

    for (DTYPE i = 0; i < NUM_VERT_PER_TRIANGLE; i++) {
      for (DTYPE j = 0; j  < NUM_VERT_PER_TRIANGLE; j++) {
        DIST_TYPE dist;
        IJK::compute_distance
          (dimension, coordA[i], coordB[j], dist);
        if (dist < min_dist) { 
          ivA = i;
          ivB = j;
          min_dist = dist; 
        }
      }
    }

    distance_data.SetDistanceV2V(min_dist, ivA, ivB);
  }



  /**
   *  @brief Compute minimum distance between vertices of two triangles.
   */
  template <typename DTYPE, typename CTYPEA, typename CTYPEB,
            typename DIST_TYPE>
  void compute_min_distance_between_vertices_of_two_triangles
  (const DTYPE dimension,
   const CTYPEA * const coordA[3], const CTYPEB * const coordB[3],
   DIST_TYPE & min_dist)
  {
    DISTANCE_DATA<DIST_TYPE,int> distance_data;

    compute_min_distance_between_vertices_of_two_triangles
      (coordA, coordB, distance_data);

    min_dist = distance_data.Distance();
  }


  /**
   *  @brief Compute minimum distance between vertices of two triangles.
   *  - Version using TRIANGLE_PROPERTIES
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_vertices_of_two_triangles
  (const TRI_TYPEA & triA_properties, const TRI_TYPEB & triB_properties,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distance_data)
  {
    const int DIM3(3);

    compute_min_distance_between_vertices_of_two_triangles
      (DIM3, triA_properties.VertexCoord(), triB_properties.VertexCoord(),
       distance_data);
  }


  /**
   *  @brief Compute minimum distance between vertices of two triangles.
   *  - Version using TRIANGLE_PROPERTIES
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE>
  void compute_min_distance_between_vertices_of_two_triangles
  (const TRI_TYPEA & triA_properties, const TRI_TYPEB & triB_properties,
   DIST_TYPE & min_dist)
  {
    const int DIM3(3);

    compute_min_distance_between_vertices_of_two_triangles
      (DIM3, triA_properties.VertexCoord(), triB_properties.VertexCoord(),
       min_dist);
  }


  /**
   *  @brief Compute minimum distance between vertices of one triangles and edges of the other.
   *  @param min_distv
   *    - Minimum distance between triangle vertices.
   *  @param[out] Minimum distance between vertices of one triangle
   *    and edgs of the other.
   *  @param[out] flag_min_dist_vert2edge True if minimum distance
   *    is between a vertex and an edge interior.
   */
  template <typename DTYPE, typename CTYPEA, 
            typename CTYPEB, typename CTYPEB2,
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_triangle_vert_and_triangle_edges
  (const DTYPE dimension,
   const CTYPEA * const coordA[3], const CTYPEB * const coordB[3],
   const CTYPEB2 * const edgeB_dir[3],
   const DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceV2V_data,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceV2E_data)
  {
    const DTYPE NUM_VERT_PER_TRIANGLE(3);

    IJK::ARRAY<DIST_TYPE> diff(dimension);
    IJK::ARRAY<DIST_TYPE> v_orth(dimension);
    DIST_TYPE mag[NUM_VERT_PER_TRIANGLE];

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */
    
    distanceV2E_data.Copy(distanceV2V_data);

    for (DTYPE i = 0; i < NUM_VERT_PER_TRIANGLE; i++) {
      for (DTYPE j0 = 0; j0 < NUM_VERT_PER_TRIANGLE; j0++) {
        const DTYPE j1 = (j0+1)%NUM_VERT_PER_TRIANGLE;

        if (does_hyperplane_separate_points<DIST_TYPE>
            (dimension, coordA[i], edgeB_dir[j0], 
             coordB[j0], coordB[j1])) {
          IJK::subtract_coord
            (dimension, coordA[i], coordB[j0], diff.Ptr());
          IJK::compute_orthogonal_vector_component
            (dimension, diff.Ptr(), edgeB_dir[j0], v_orth.Ptr());
          DIST_TYPE dist;
          IJK::compute_magnitude(dimension, v_orth.Ptr(), dist);
          if (dist < distanceV2E_data.Distance()) {
            distanceV2E_data.SetDistanceV2E(dist, i, j0);
          }
        }
      }
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "Leaving " << __func__ << endl;
    */
  }


  /**
   *  @brief Compute minimum distance between vertices of one triangles and edges of the other.
   *  - Version using TRIANGLE_PROPERTIES.
   *  @param min_distv
   *    - Minimum distance between triangle vertices.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB, 
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_triangle_vert_and_triangle_edges
  (const TRI_TYPEA & triA_properties, const TRI_TYPEB & triB_properties,
   const DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceV2V_data,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceV2E_data)
  {
    const int DIM3(3);

    compute_min_distance_between_triangle_vert_and_triangle_edges
      (DIM3, triA_properties.VertexCoord(), triB_properties.VertexCoord(),
       triB_properties.EdgeDir(), distanceV2V_data, distanceV2E_data);
  }


  /**
   *  @brief Compute minimum distance between edge interiors of two triangles.
   *  - Assumes that min distance between vertices of one triangle and
   *    edges of the other have already been computed.
   *  @param distanceVE_data
   *    - Information on minimum distance between vertices of one triangle
   *      and edges of the other. Could be V2E or E2V.
   */
  template <typename DTYPE, typename CTYPEA, typename CTYPEB,
            typename CTYPEA2, typename CTYPEB2,
            typename LTYPEA, typename LTYPEB,
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_edge_interiors_of_two_triangles
  (const DTYPE dimension,
   const CTYPEA * const coordA[3], const CTYPEA2 * const edgeA_dir[3],
   const LTYPEA edgeA_length[3],
   const CTYPEB * const coordB[3], const CTYPEB2 * const edgeB_dir[3],
   const LTYPEB edgeB_length[3],
   const DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceVE_data,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceE2E_data)
  {
    const int NUM_VERT_PER_TRIANGLE(3);

    IJK::ARRAY<DIST_TYPE> v_orth(dimension);
    DIST_TYPE magA[NUM_VERT_PER_TRIANGLE];
    DIST_TYPE magB[NUM_VERT_PER_TRIANGLE];


    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  In " << __func__ << endl;
    */

    
    distanceE2E_data.Copy(distanceVE_data);

    for (DTYPE jA0 = 0; jA0 < NUM_VERT_PER_TRIANGLE; jA0++) {
      const DTYPE jA1 = (jA0+1)%NUM_VERT_PER_TRIANGLE;

      for (DTYPE jB0 = 0; jB0 < NUM_VERT_PER_TRIANGLE; jB0++) {
        const DTYPE jB1 = (jB0+1)%NUM_VERT_PER_TRIANGLE;

        DIST_TYPE dist;
        compute_line_to_line_distance
          (dimension, coordA[jA0], edgeA_dir[jA0],
           coordB[jB0], edgeB_dir[jB0], v_orth.Ptr(), dist);

        if (do_projected_line_segments_cross<DIST_TYPE>
            (dimension, coordA[jA0], coordA[jA1],
             edgeA_dir[jA0], edgeA_length[jA0], 
             coordB[jB0], coordB[jB1],
             edgeB_dir[jB0], edgeB_length[jB0], v_orth.Ptr())) {
          if (dist < distanceE2E_data.Distance()) {
            distanceE2E_data.SetDistanceE2E(dist, jA0, jB0); 

            // *** DEBUG ***
            /*
            using namespace std;
            cerr << "    e2e min dist: " << dist << endl;
            IJK::print_coord3D
              (cerr, "    edgeA_dir: ", edgeA_dir[jA0], "\n");
            IJK::print_coord3D
              (cerr, "    edgeB_dir: ", edgeB_dir[jB0], "\n");
            IJK::print_coord3D
              (cerr, "    v_orth: ", v_orth.Ptr(), "\n");
            */
          }
        }
      }
    }
  }


  /**
   *  @brief Compute minimum distance between edge interiors of two triangles.
   *  - Version using TRIANGLE_PROPERTIES.
   *  @param min_dist_vert2edge
   *    - Minimum distance between vertices of one triangle
   *      and edges of the other.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB, 
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_edge_interiors_of_two_triangles
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   const DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceVE_data,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distanceE2E_data)
  {
    const int DIM3(3);

    compute_min_distance_between_edge_interiors_of_two_triangles
      (DIM3, triA_properties.VertexCoord(), triA_properties.EdgeDir(),
       triA_properties.EdgeLength(),
       triB_properties.VertexCoord(), triB_properties.EdgeDir(),
       triB_properties.EdgeLength(),
       distanceVE_data, distanceE2E_data);
  }


  /**
   *  @brief Return true if point projects onto triangle along normal direction.
   *  @tparam DIST_TYPE Distance type. Must be given in function call.
   *  @param edgeB_orth_dir[i] Direction orthogonal to edgeB_dir[i] and
   *     normal direction.
   */
  template <typename DIST_TYPE, typename DTYPE, 
            typename CTYPEA, typename CTYPEB,
            typename CTYPEB2, typename CTYPEBN>
  bool does_point_project_onto_triangle
  (const DTYPE dimension,
   const CTYPEA coordA[], 
   const CTYPEB * const coordB[3], 
   const CTYPEB2 * const edgeB_orth_dir[3],
   const CTYPEBN normalB[])
  {
    const int NUM_VERT_PER_TRIANGLE(3);

    typedef CTYPEB2 CTYPE;

    IJK::ARRAY<CTYPE> v_orth(dimension);
    IJK::ARRAY<CTYPE> proj_coordA(dimension);

    IJK::subtract_coord
      (dimension, coordA, coordB[0], v_orth.Ptr());

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    IJK::print_coord3D(cerr, "  Vector B to A: ", v_orth.Ptr(), "\n");
    */

    IJK::compute_orthogonalII_vector
      (dimension, v_orth.Ptr(), normalB, v_orth.Ptr());

    IJK::subtract_coord
      (dimension, coordA, v_orth.Ptr(), proj_coordA.Ptr());
    IJK::add_coord
      (dimension, coordB[0], proj_coordA.Ptr(), proj_coordA.Ptr());
    
    for (int j0 = 0; j0 < NUM_VERT_PER_TRIANGLE; j0++) {
      const int j2 = (j0+2)%NUM_VERT_PER_TRIANGLE;

      if (does_hyperplane_separate_points<DIST_TYPE>
          (dimension, coordB[j0], edgeB_orth_dir[j0],
           coordB[j2], proj_coordA.Ptr()))
        { return(false); }
    }

    return(true);
  }


  /**
   *  @brief Compute minimum distance between triangle vertices and triangle.
   *  @param min_dist_vert2edge
   *    - Minimum distance between vertices of one triangle
   *      and edges of the other.
   */
  template <typename DTYPE, typename CTYPEA, 
            typename CTYPEB, typename CTYPEB2, typename CTYPEB3, 
            typename DIST_TYPE>
  void compute_min_distance_between_triangle_vertices_and_triangle
  (const DTYPE dimension,
   const CTYPEA * const coordA[3],
   const CTYPEB * const coordB[3], 
   const CTYPEB2 * const edgeB_orth_dir[3],
   const CTYPEB3 normalB[],
   const DIST_TYPE min_dist_vert2edge,
   DIST_TYPE & min_dist_vert2triangle)
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    const int DIM3(3);

    typedef CTYPEB2 CTYPE;

    IJK::ARRAY<CTYPE> proj_coordA(dimension);
    IJK::ARRAY<CTYPE> v_orth(dimension);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */
    
    min_dist_vert2triangle = min_dist_vert2edge;

    for (int i = 0; i < NUM_VERT_PER_TRIANGLE; i++) {

      DIST_TYPE dist;
      compute_distance_to_hyperplane
        (dimension, coordA[i], coordB[0], normalB, dist);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << endl;
      IJK::print_coord3D
        (cerr, "  Distance from ", coordA[i], " to plane: ");
      cerr << dist << endl;
      */

      if (dist >= min_dist_vert2triangle) { 
        // Point coordA[i] is farther than min_dist_vert2triangle.
        // Skip point.
        continue; 
      }

      if (does_point_project_onto_triangle<DIST_TYPE>
          (dimension, coordA[i], coordB, edgeB_orth_dir, normalB)) {

        // Point projects onto triangle.
        // Update min_dist_vert2triangle.
        min_dist_vert2triangle = dist;
      }
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  Leaving " << __func__ << endl;
    */
  }


  /**
   *  @brief Compute minimum distance between triangle vertices and triangle.
   *  - Version using TRIANGLE_PROPERTIES.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB, typename DIST_TYPE>
  void compute_min_distance_between_triangle_vertices_and_triangle
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   const DIST_TYPE min_dist_vert2edge,
   DIST_TYPE & min_dist_vert2triangle)
  {
    const int DIM3(3);

    compute_min_distance_between_triangle_vertices_and_triangle
      (DIM3, triA_properties.VertexCoord(), triB_properties.VertexCoord(),
       triB_properties.EdgeOrthDir(), triB_properties.Normal(),
       min_dist_vert2edge, min_dist_vert2triangle);
  }

  ///@}


  // *****************************************************************
  //! @name Triangle distance.
  // *****************************************************************

  ///@{

  /**
   *  @brief Return min distance between edges of two triangles.
   *  - Returns 0 if triangle edges intersect.
   *  - Computes min distance between two sets of triangle vertices and
   *    between vertices of one triangle and edges of the other.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_triangle_edges
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distance_data)
  {
    DISTANCE_DATA<DIST_TYPE,ITYPEF> distanceV2V;
    DISTANCE_DATA<DIST_TYPE,ITYPEF> distanceV2E;
    DISTANCE_DATA<DIST_TYPE,ITYPEF> distanceE2V;

    compute_min_distance_between_vertices_of_two_triangles
      (triA_properties, triB_properties, distanceV2V);

    compute_min_distance_between_triangle_vert_and_triangle_edges
      (triA_properties, triB_properties, distanceV2V, distanceV2E);

    compute_min_distance_between_triangle_vert_and_triangle_edges
      (triB_properties, triA_properties, distanceV2V, distanceE2V);
    if (distanceE2V.DistancePairType() == V2E_DIST) 
      { distanceE2V.SwapPolyData(); }

    if (distanceV2E.Distance() <= distanceE2V.Distance()) {
      compute_min_distance_between_edge_interiors_of_two_triangles
        (triA_properties, triB_properties, distanceV2E, distance_data);
    }
    else {
      compute_min_distance_between_edge_interiors_of_two_triangles
        (triA_properties, triB_properties, distanceE2V, distance_data);
    }

  }


  /**
   *  @brief Return min distance between edges of two triangles.
   *  - Returns 0 if triangle edges intersect.
   *  - Computes min distance between two sets of triangle vertices and
   *    between vertices of one triangle and edges of the other.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE>
  void compute_min_distance_between_triangle_edges
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   DIST_TYPE & min_dist)
  {
    DISTANCE_DATA<DIST_TYPE,int> distance_data;

    compute_min_distance_between_triangle_edges
      (triA_properties, triB_properties, distance_data);

    min_dist = distance_data.Distance();
  }


  /**
   *  @brief Return min distance between two non-intersecting triangles.
   *  - Returns 0 if triangles intersect.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_non_intersecting_triangles
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distance_data)
  {
    DISTANCE_DATA<DIST_TYPE,ITYPEF> distanceE2E;

    DIST_TYPE min_dist_edge2edge;
    DIST_TYPE min_dist_vertA2triangleB, min_dist_vertB2triangleA;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  In " << __func__ << endl;
    */

    compute_min_distance_between_triangle_edges
      (triA_properties, triB_properties, distanceE2E);
    min_dist_edge2edge = distanceE2E.Distance();

    // *** DEBUG ***
    /*
    using namespace std;
    distanceE2E.Print(cerr, "  ");
    */
    
    compute_min_distance_between_triangle_vertices_and_triangle
      (triA_properties, triB_properties, min_dist_edge2edge,
       min_dist_vertA2triangleB);
    compute_min_distance_between_triangle_vertices_and_triangle
      (triA_properties, triB_properties, min_dist_edge2edge,
       min_dist_vertB2triangleA);

    const DIST_TYPE min_distVT = 
      std::min(min_dist_vertA2triangleB, min_dist_vertB2triangleA);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  min_distVT: " << min_distVT << endl;
    */

    if (distanceE2E.Distance() <= min_distVT)
      { distance_data.Copy(distanceE2E); }
    else
      { distance_data.SetDistance(min_distVT); }

    return;
  }


  /**
   *  @brief Return min distance between two triangles.
   *  - Returns 0 if triangles intersect.
   *  - Intersection test based on "Fast Triangle-Triangle Intersection Tests"
   *    by Devillers and Guigue, 2006.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE, typename ITYPEF>
  void compute_min_distance_between_triangles
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   DISTANCE_DATA<DIST_TYPE,ITYPEF> & distance_data)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_TRIANGLE(3);

    DIST_TYPE dist2planeA[NUM_VERT_PER_TRIANGLE];
    DIST_TYPE dist2planeB[NUM_VERT_PER_TRIANGLE];
    ITYPEF jA0, jB0, jA1, jB1, jA2, jB2;
    DIST_TYPE min_dist;

    const DIST_TYPE distA = 
      triA_properties.ComputeSignedDistanceToPlane
      (triB_properties.VertexCoord(0), triB_properties.Normal(), 
       dist2planeB, jA0);
    const DIST_TYPE distB = 
      triB_properties.ComputeSignedDistanceToPlane
      (triA_properties.VertexCoord(0), triA_properties.Normal(), 
       dist2planeA, jB0);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "    distA: " << distA
         << "  distB: " << distB << endl;
    */

    if ((distA == 0) && (distB == 0)) {
      jA1 = (jA0+1)%DIM3;
      jA2 = (jA1+1)%DIM3;
      jB1 = (jB0+1)%DIM3;
      jB2 = (jB1+1)%DIM3;

      // *** DEBUG ***
      /*
      cerr << "    dist2planeB[" << jA0 << "]: "
           << dist2planeB[jA0] << endl;
      cerr << "    dist2planeA[" << jB0 << "]: "
           << dist2planeA[jB0] << endl;
      */

      if (dist2planeB[jA0] < 0)
        { std::swap(jA1, jA2); }

      if (dist2planeA[jB0] < 0)
        { std::swap(jB1, jB2); }

      const ITYPEF ieA = triA_properties.EdgeIndex(jA0,jA1);
      const ITYPEF ieB = triB_properties.EdgeIndex(jB0,jB1);
      if (does_edge_intersect_triangle_interior
          (triA_properties, triB_properties, dist2planeB, ieA)) {
        distance_data.SetDistanceE2T(0, ieA);
        return;
      }
      else if (does_edge_intersect_triangle_interior
               (triB_properties, triA_properties, dist2planeA, ieB)) {
        distance_data.SetDistanceT2E(0, ieB);
        return;
      }
      else {
        // Triangles do not intersect.
        // Each triangle intersects the plane containing the other triangle.
        compute_min_distance_between_triangle_edges
          (triA_properties, triB_properties, distance_data);
        return;
      }
    }

    // Triangles do not intersect.
    compute_min_distance_between_non_intersecting_triangles
      (triA_properties, triB_properties, distance_data);

    return;
  }

  /**
   *  @brief Return min distance between two triangles.
   *  - Returns 0 if triangles intersect.
   *  - Intersection test based on "Fast Triangle-Triangle Intersection Tests"
   *    by Devillers and Guigue, 2006.
   */
  template <typename TRI_TYPEA, typename TRI_TYPEB,
            typename DIST_TYPE>
  void compute_min_distance_between_triangles
  (const TRI_TYPEA & triA_properties,
   const TRI_TYPEB & triB_properties,
   DIST_TYPE & min_dist)
  {
    DISTANCE_DATA<DIST_TYPE,int> distance_data;

    compute_min_distance_between_triangles
      (triA_properties, triB_properties, distance_data);

    min_dist = distance_data.Distance();
  }

  ///@}


  // *****************************************************************
  //! @name Class TRIANGLE3D_PROPERTIES member functions.
  // *****************************************************************

  ///@{

  // Initialize TRIANGLE3D_PROPERTIES.
  template <typename CTYPE0, typename CTYPE1>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::Init()
  {
    for (int i = 0; i < NUM_VERTICES; i++) {
      vertex_coord[i] = 0;
      edge_dir[i] = edge_dir_coord + i*DIM3;
      edge_orth_dir[i] = edge_orth_dir_coord + i*DIM3;
    }
  }


  // Set pointers to i'th vertex coord from vertex list.
  // - Polygon has at most three vertices (typically exactly three.)
  template <typename CTYPE0, typename CTYPE1>
  template <typename ITYPE, typename IVTYPE, typename NTYPE0>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::
  SetVertexCoordFromVertexList
  (const ITYPE i, const CTYPE0 vertex_coord_list[],
   const IVTYPE * vertex_list, const NTYPE0 num_poly_vert)
  {
    IJK::PROCEDURE_ERROR error
      ("TRIANGLE3D_PROPERTIES::SetVertexCoordFromVertexList");

    if (num_poly_vert > NUM_VERTICES) {
      error.AddMessage
        ("Error. num_poly_vert can be at most ", NUM_VERTICES, ".");
      error.AddMessage("  num_poly_vert = ", num_poly_vert, "");
      throw error;
    }
    else if (num_poly_vert < 1) {
      error.AddMessage
        ("Error. num_poly_vert must be a positive integer.");
      error.AddMessage("  num_poly_vert = ", num_poly_vert, "");
      throw error;
    }

    if (vertex_coord_list == NULL) {
      error.AddMessage("Programming error. Array vertex_coord_list[] is NULL.");
      throw error;
    }

    if (num_poly_vert == NUM_VERTICES) {
      for (int j = 0; j < NUM_VERTICES; j++) {
        const IVTYPE jv = vertex_list[j];
        SetVertexCoord(j, vertex_coord_list+jv*DIM3);
      }
    }
    else {
      // Degenerate case where polygon has one or two vertices.

      const COORD_TYPE * prev_coord = NULL;
      for (int j = 0; j < NUM_VERTICES; j++) {
        if (j < num_poly_vert) {
          const IVTYPE jv = vertex_list[j];
          SetVertexCoord(j, vertex_coord_list+jv*DIM3);
          prev_coord = vertex_coord_list + jv*DIM3;
        }
        else {
          // Copy coordinates of previous vertex.
          SetVertexCoord(j, prev_coord);
        }
      }
    }
  }


  // Compute edge directions.
  // @pre All values of vertex_coord[] are set.
  template <typename CTYPE0, typename CTYPE1>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::ComputeEdgeDir()
  {
    IJK::PROCEDURE_ERROR error("TRIANGLE3D_PROPERTIES::ComputeEdgeDir");

    if (!CheckVertexCoord(error)) { throw error; }

    compute_triangle_edge_vectors_oriented
      (DIM3, vertex_coord[0], vertex_coord[1], vertex_coord[2], 
       edge_dir[0], edge_dir[1], edge_dir[2]);

    for (int i = 0; i < NUM_VERTICES; i++) {
      normalize_vector_robust
        (DIM3, edge_dir[i], 0, edge_dir[i], edge_length[i]); 
    }

  }


  // Compute normal.
  // - Also computes edge_orth_dir[] for each edge.
  // @pre All values of vertex_coord[] are set.
  // @pre All edge directions are set.
  template <typename CTYPE0, typename CTYPE1>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::ComputeNormal()
  {
    IJK::PROCEDURE_ERROR error("TRIANGLE3D_PROPERTIES::ComputeNormal");

    if (!CheckVertexCoord(error)) { throw error; }

    compute_triangle_normal
      (DIM3, edge_dir[0], edge_dir[1], edge_dir[2], edge_length, normal);

    for (int j = 0; j < NUM_VERTICES; j++) {
      compute_orthogonalII_vector
        (DIM3, normal, edge_dir[j], edge_orth_dir[j]);
    }

    const int j0 = get_index_of_max_coord(NUM_VERTICES, edge_length);
    const int j1 = (j0+1)%DIM3;
    const int j2 = (j0+2)%DIM3;

    // Orient edge_orth_dir[j0] away from vertex (j+2)%DIM3.
    orient_vector_away_from(DIM3, edge_dir[j1], edge_orth_dir[j0]);

    // - Because edge j0 has longest length, edge_orth_dir[j1] and
    //   edge_orth_dir[j2] should point away from edge_orth_dir[j0].
    orient_vector_away_from
      (DIM3, edge_orth_dir[j0], edge_orth_dir[j1]);
    orient_vector_away_from
      (DIM3, edge_orth_dir[j0], edge_orth_dir[j2]);
  }


  // Compute signed distance of three triangle vertices to plane.
  // - Return 0 if triangle intersects plane.
  // - Return signed distance to closest point if triangle 
  //   does not intersect plan.
  template <typename CTYPE0, typename CTYPE1>
  template <typename CTYPEQ, typename CTYPEN, typename DIST_TYPE,
            typename ITYPE>
  DIST_TYPE TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::
  ComputeSignedDistanceToPlane
    (const CTYPEQ q0, const CTYPEN qnormal,
     DIST_TYPE signed_distance[3], ITYPE & jindex) const
  {
    const DIST_TYPE dist = 
      compute_signed_distance_triangle_to_hyperplane
      (DIM3, VertexCoord(), q0, qnormal, signed_distance, jindex);

    return(dist);
  }


  // Check that all vertex coordinates are set, i.e., that each vertex_coord[i] is not NULL.
  template <typename CTYPE0, typename CTYPE1>
  bool TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::
  CheckVertexCoord(IJK::ERROR & error) const
  {
    for (int i = 0; i < NUM_VERTICES; i++) {
      if (vertex_coord[i] == NULL) {
        error.AddMessage
          ("Programming error. NULL pointer to vertex coordinates.");
        error.AddMessage
          ("  vertex_coord[", i, "] = NULL.");
        return(false);
      }
    }

    return(true);
  }


  // Print coordinates.
  // - Extended version.  Define left/right delimiters and separator.
  template <typename CTYPE0, typename CTYPE1>
  template <typename OSTREAM_TYPE, typename CTYPEX>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::_PrintCoordX
  (OSTREAM_TYPE & out, const CTYPEX * coord,
   const char c0, const char c1, const char c2) const
  {
    out << c0 << coord[0] << c1 << coord[1] << c1 << coord[2] << c2;
  }


  // Print coordinates.
  template <typename CTYPE0, typename CTYPE1>
  template <typename OSTREAM_TYPE, typename CTYPEX,
            typename STYPE0, typename STYPE1>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::_PrintCoord
  (OSTREAM_TYPE & out, const STYPE0 & s0, 
   const CTYPEX * coord, const STYPE1 & s1) const
  {
    out << s0;
    _PrintCoord(out, coord);
    out << s1;
  }


  // Print all three edge lengths.
  template <typename CTYPE0, typename CTYPE1>
  template <typename OSTREAM_TYPE, typename STYPE0, typename STYPE1>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::PrintEdgeLengths
  (OSTREAM_TYPE & out, const STYPE0 & s0, const STYPE1 & s1) const
  {
    out << s0;
    for (int i = 0; i < NUM_VERTICES; i++) {
      out << EdgeLength(i);
      out << "  ";
    }
    out << s1;
  }


  // Print vertex, edge_dir, edge_length, edge_orth_dir[] coordinates.
  template <typename CTYPE0, typename CTYPE1>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void TRIANGLE3D_PROPERTIES<CTYPE0,CTYPE1>::Print
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix;
    out << "Vertex coordinates: ";
    for (int i = 0; i < NUM_VERTICES; i++) 
      { PrintVertexCoord(out, " ", i, ""); }
    out << "\n";

    out << line_prefix;
    PrintNormal(out, "Normal direction: ", "\n");

    out << line_prefix;
    out << "Edge directions: ";
    for (int i = 0; i < NUM_VERTICES; i++) 
      { PrintEdgeDir(out, " ", i, ""); }
    out << "\n";

    out << line_prefix;
    PrintEdgeLengths(out, "Edge lengths: ", "\n");

    out << line_prefix;
    out << "Edge orthogonal direction: ";
    for (int i = 0; i < NUM_VERTICES; i++) 
      { PrintEdgeOrthDir(out, " ", i, ""); }
    out << "\n";
  }

  ///@}


  // *****************************************************************
  //! @name = Class DISTANCE_DATA member functions.
  // *****************************************************************
  
  ///@{

  template <typename DIST_TYPE, typename ITYPEF>
  void DISTANCE_DATA<DIST_TYPE,ITYPEF>::Init()
  {
    distance = std::numeric_limits<DIST_TYPE>::max();
    distance_pair_type = DIST_PAIR_TYPE_NOT_SET;
    iclosest_face_index[0] = 0;
    iclosest_face_index[1] = 0;
  }


  // Set distance to closest.
  template <typename DIST_TYPE, typename ITYPEF>
  template <typename DIST_TYPE2, typename ITYPE2>
  void DISTANCE_DATA<DIST_TYPE,ITYPEF>::_SetDistance
  (const DISTANCE_PAIR_TYPE distance_pair_type,
   const DIST_TYPE2 dist, const ITYPE2 iface0, const ITYPE2 iface1)
  {
    this->distance_pair_type = distance_pair_type;
    distance.Set(dist);
    iclosest_face_index[0] = iface0;
    iclosest_face_index[1] = iface1;
  }


  // Swap first and second poly data.
  // - Swap iclosest_facet_index[0] and iclosest_facet_index[1].
  // - Change distance_pair_type.
  template <typename DIST_TYPE, typename ITYPEF>
  void DISTANCE_DATA<DIST_TYPE,ITYPEF>::SwapPolyData()
  {
    std::swap(iclosest_face_index[0], iclosest_face_index[1]);

    if (DistancePairType() == V2E_DIST) 
      { distance_pair_type = E2V_DIST; }
    else if (DistancePairType() == E2V_DIST) 
      { distance_pair_type = V2E_DIST; }
    else if (DistancePairType() == V2T_DIST) 
      { distance_pair_type = T2V_DIST; }
    else if (DistancePairType() == T2V_DIST) 
      { distance_pair_type = V2T_DIST; }
    else if (DistancePairType() == E2T_DIST) 
      { distance_pair_type = T2E_DIST; }
    else if (DistancePairType() == T2E_DIST) 
      { distance_pair_type = E2T_DIST; }
  }


  // Copy distance_data2 to this.
  template <typename DIST_TYPE, typename ITYPEF>
  template <typename DIST_DATA_TYPE2>
  void DISTANCE_DATA<DIST_TYPE,ITYPEF>::Copy
  (const DIST_DATA_TYPE2 & right)
  {
    if (&right != this) {
      *this = right;
    }
    
  }


  // Print distance type.
  template <typename DIST_TYPE, typename ITYPEF>
  template <typename OSTREAM_TYPE>
  void DISTANCE_DATA<DIST_TYPE,ITYPEF>::
  PrintDistancePairType(OSTREAM_TYPE & out) const
  {
    if (DistancePairType() == V2V_DIST) {
      out << "Distance between vertices "
          << IndexOfFace(0) << " and "
          << IndexOfFace(1) << ".";
    }
    else if (DistancePairType() == V2E_DIST) {
      out << "Distance between vertex "
          << IndexOfFace(0) << " and edge "
          << IndexOfFace(1) << ".";
    }
    else if (DistancePairType() == E2V_DIST) {
      out << "Distance between edge "
          << IndexOfFace(0) << " and vertex "
          << IndexOfFace(1) << ".";
    }
    else if (DistancePairType() == E2E_DIST) {
      out << "Distance between interiors of edges "
          << IndexOfFace(0) << " and "
          << IndexOfFace(1) << ".";
    }
    else if (DistancePairType() == V2T_DIST) {
      out << "Distance between vertex "
          << IndexOfFace(0) << " and triangle.";
    }
    else if (DistancePairType() == T2V_DIST) {
      out << "Distance between triangle and vertex "
          << IndexOfFace(1) << ".";
    }
    else if (DistancePairType() == E2T_DIST) {
      out << "Distance between edge "
          << IndexOfFace(0) << " and triangle.";
    }
    else if (DistancePairType() == T2E_DIST) {
      out << "Distance between triangle and edge "
          << IndexOfFace(1) << ".";
    }
    else {
      out << "Distance pair type not set.";
    }
  }


  // Print distance type.
  template <typename DIST_TYPE, typename ITYPEF>
  template <typename OSTREAM_TYPE, typename STYPE0, typename STYPE1>
  void DISTANCE_DATA<DIST_TYPE,ITYPEF>::PrintDistancePairType
  (OSTREAM_TYPE & out, const STYPE0 & s0, const STYPE1 & s1) const
  {
    out << s0;
    PrintDistancePairType(out);
    out << s1;
  }


  // Print distance, distance type and facet info.
  template <typename DIST_TYPE, typename ITYPEF>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void DISTANCE_DATA<DIST_TYPE,ITYPEF>::Print
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix << "Distance: " << Distance() << "\n";
    PrintDistancePairType(out, line_prefix, "\n");
  }

  ///@}

}


#endif
