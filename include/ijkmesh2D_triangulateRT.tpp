/*!
 *  @file ijkmesh2D_triangulateRT.tpp
 *  @brief ijk template functions for triangulating by replacing triangles.
 *  - Replace triangle with triangle and edge.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2021-2022 Rephael Wenger

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

#ifndef _IJKMESH2D_TRIANGULATE_RT_
#define _IJKMESH2D_TRIANGULATE_RT_

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkmesh2D_triangulate.tpp"
#include "ijkmesh2D_angle_triangle.tpp"

#include <algorithm>
#include <numeric>
#include <vector>

namespace IJK {

  // *****************************************************************
  /// @name Add split vertex and replace triangle
  // *****************************************************************

  ///@{

  /**
     @brief Add vertex and replace triangle with split edge by triangle edge.
  */
  template <typename MESH_TYPE, typename CTYPE2, typename CTYPE3,
            typename CTYPE4, typename CTYPE5,
            typename ITYPEA, typename ITYPEB,
            typename ITYPEV,
            typename COS_TYPE2, typename NTYPE2,
            typename IPOLY_TYPE, int BIT_SET_SIZE>
  void add_vertex_replace_triangle_with_splitE1_by_triangle_edge
  (MESH_TYPE & mesh, std::vector<CTYPE2> & vcoord,
   const std::vector<CTYPE3> & new_triangle_vcoord,
   const std::vector<CTYPE4> & interior_vcoord0,
   const std::vector<CTYPE5> & interior_vcoord2,
   const ITYPEA ibase_half_edge, const ITYPEB ilong_half_edge,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE2,NTYPE2> 
   poly_tri_info[3],
   ITYPEV & iv_new,
   IPOLY_TYPE ipoly_new[3])
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    IPOLY_TYPE ipolyY0, ipolyY1;
    ITYPEV iv_split;

    const DIMENSION_TYPE dimension = mesh.Dimension();

    // If true, ilong_half_edge is next half edge of ibase_half_edge.
    const bool flag_is_next_half_edge =
      mesh.IsNextHalfEdge(ibase_half_edge, ilong_half_edge);
    const HALF_EDGE_INDEX_TYPE ilong_half_edgeX =
      mesh.IndexOfNextHalfEdgeAroundEdge(ilong_half_edge);

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "*** Replacing triangle ", ilong_half_edge, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "*** Splitting half edge ", ilong_half_edge, "\n");
    mesh.PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "  Adjacent polygon ", ilong_half_edgeX, "\n");
    */

    // First, split long edge.
    HALF_EDGE_INDEX_TYPE ilong_half_edgeY;
    HALF_EDGE_INDEX_TYPE ibase_half_edgeY;
    mesh.MESH_TYPE::MESH2D_SPLIT_II_BASE_TYPE::AddVertexSplitPolygonEdgeII
      (dimension, new_triangle_vcoord, vcoord, ilong_half_edge,
       iv_split, ipolyY1, ipolyY0, ilong_half_edgeY);
    mesh.SetVertexType(iv_split, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);
    if (flag_is_next_half_edge) {
      ibase_half_edgeY = mesh.IndexOfPrevHalfEdgeInPolygon(ilong_half_edgeY);
    }
    else {
      ibase_half_edgeY = mesh.IndexOfKthNextHalfEdgeInPolygon
        (ilong_half_edgeY, 2);
    }

    mesh.SetNumSplitEdgePolygonVertices(ipolyY0);

    // Next, merge two triangle split edges.
    HALF_EDGE_INDEX_TYPE ibase_half_edge_new;
    mesh.MESH_TYPE::MESH2D_SPLIT_II_BASE_TYPE::
      AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
      (dimension, new_triangle_vcoord, vcoord,
       ibase_half_edgeY, iv_new, ibase_half_edge_new);

    const HALF_EDGE_INDEX_TYPE isplit_half_edge_new = 
      mesh.GetIndexOfAdjacentHalfEdgeInPolygon
      (ibase_half_edge_new, !flag_is_next_half_edge);

    ipoly_new[0] = ipolyY0;
    ipoly_new[1] = mesh.IndexOfPolygonContainingHalfEdge(ibase_half_edge_new);
    ipoly_new[2] = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(isplit_half_edge_new);

    mesh.SetVertexType(iv_new, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);

    mesh.SetNumSplitEdgePolygonVertices(ipoly_new[0]);
    mesh.SetNumSplitEdgePolygonVertices(ipoly_new[1]);
    mesh.SetNumSplitEdgePolygonVertices(ipoly_new[2]);

    // *** DEBUG **
    /*
    using namespace std;
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  New polygon: ", ipolyY0, "\n");
    poly_tri_info[0].Print(cerr, "    ", mesh.NumPolygonVertices(ipolyY0));
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  New polygon: ", ipoly_new[2], "\n");
    poly_tri_info[2].Print(cerr, "    ", mesh.NumPolygonVertices(ipoly_new[2]));
    */

    mesh.CopyTriInfoAndInteriorCoord
      (poly_tri_info[0], interior_vcoord0, ipoly_new[0]);
    mesh.CopyTriInfo(poly_tri_info[1], ipoly_new[1]);
    mesh.CopyTriInfoAndInteriorCoord
      (poly_tri_info[2], interior_vcoord2, ipoly_new[2]);
  }
   

  /**
     @brief Add vertex and replace triangle with split edge by triangle edge.
     - Version that does not return new vertex or polygon indices.
  */
  template <typename MESH_TYPE, typename CTYPE2, typename CTYPE3,
            typename CTYPE4, typename CTYPE5,
            typename ITYPEA, typename ITYPEB,
            typename COS_TYPE2, typename NTYPE2, int BIT_SET_SIZE>
  void add_vertex_replace_triangle_with_splitE1_by_triangle_edge
  (MESH_TYPE & mesh, std::vector<CTYPE2> & vcoord,
   const std::vector<CTYPE3> & new_triangle_vcoord, 
   const std::vector<CTYPE4> & interior_vcoord0,
   const std::vector<CTYPE5> & interior_vcoord2,
   const ITYPEA ibase_half_edge,
   const ITYPEB ilong_half_edge,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE2,NTYPE2> 
   poly_tri_info[3])
  {
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE THREE_POLYGONS(3);
    POLYGON_INDEX_TYPE ipoly_new[THREE_POLYGONS];
    VERTEX_INDEX_TYPE iv_new;

    add_vertex_replace_triangle_with_splitE1_by_triangle_edge
      (mesh, vcoord, new_triangle_vcoord, interior_vcoord0, interior_vcoord2,
       ibase_half_edge, ilong_half_edge, poly_tri_info, iv_new, ipoly_new);
  }

  ///@}


  // *****************************************************************
  /// @name Replace triangles
  // *****************************************************************

  ///@{

  /**
     @brief Replace one triangle to max min triangulation angle.
     @pre All checks have been done to ensure that itriangle
          is a candidate for replacement.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param ihalf_edge Triangle base.
       - Replace vertex opposite triangle base.
     @param flag_allow_centroidx2  If true, try two different centroids
         as interior triangulation vertex.
  */
  template <typename MESH_TYPE, typename CTYPE, 
            typename ITYPE, typename MTYPE>
  void replace_triangleI_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edge0,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_max_min_all)
  {
    mesh.ReplaceTriangleWithTriangleEdgeToMaxMinTriangulationAngle
      (ihalf_edge0, max_small_magnitude, 
       triangulation_settings, flag_allow_centroidx2,
       flag_max_min_all, vcoord);
  }


  /**
     @brief Replace one triangle to max min triangulation angle.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param ihalf_edge Triangle base.
       - Replace vertex opposite triangle base.
     @param flag_allow_centroidx2  If true, try two different centroids
         as interior triangulation vertex.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename PTYPE, typename COS_TYPE, typename MTYPE>
  void replace_triangleI_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold,
   const PTYPE itriangle,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE THREE(3);
    const DIMENSION_TYPE dimension = mesh.Dimension();

    COS_TYPE cos_min_angle;
    NUMBER_TYPE iangle;
    bool flag_zero;

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord, mesh.PolygonVertexList(itriangle), 
       max_small_magnitude, cos_min_angle, iangle, flag_zero);

    VERTEX_INDEX_TYPE iv_apex;
    HALF_EDGE_INDEX_TYPE ihalf_edge0;
    mesh.GetTriangleVertexAndOpposingEdge
      (itriangle, iangle, iv_apex, ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge1 =
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2 =
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge1);

    const POLYGON_INDEX_TYPE ipolyA =
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge1);
    const POLYGON_INDEX_TYPE ipolyC =
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge2);

    if (flag_triangle_has_min_angle) {
      if (cos_min_angle < 
          mesh.PolygonCosMinTriangulationAngle(ipolyA)) { return; }

      if (cos_min_angle <
          mesh.PolygonCosMinTriangulationAngle(ipolyC)) { return; }

      if (cos_min_angle <= cos_angle_threshold)
        { return; } 
    }
    else {
      COS_TYPE cos_min_angle_poly3 = 
        std::max(mesh.PolygonCosMinTriangulationAngle(ipolyA),
                   cos_min_angle);
      cos_min_angle_poly3 = 
        std::max(cos_min_angle_poly3,
                 mesh.PolygonCosMinTriangulationAngle(ipolyC));

      if (cos_min_angle_poly3 <= cos_angle_threshold) 
        { return; }
    }

    replace_triangleI_to_max_min_triangulation_angle
      (mesh, vcoord, ihalf_edge0, max_small_magnitude, 
       triangulation_settings,
       flag_allow_centroidx2, flag_max_min_all);
  }


  /**
     @brief Replace triangles to max min triangulation angle.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param flag_allow_centroidx2  If true, try two different centroids
         as interior triangulation vertex.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void replace_triangles_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      if (!mesh.IsTriangle(ipoly)) { continue; }

      replace_triangleI_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold, ipoly, max_small_magnitude, 
         triangulation_settings, flag_allow_centroidx2,
         flag_triangle_has_min_angle, flag_max_min_all);
    }
  }


  /**
     @brief Replace triangles to max min triangulation angle.
     - Version with array cos_angle_threshold[].
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param flag_allow_centroid  If true, allow adding triangulation
         vertices at centroids of adjacent polygons.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void replace_triangles_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      replace_triangles_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude, 
         triangulation_settings,
         flag_allow_centroidx2,
         flag_triangle_has_min_angle, flag_max_min_all);

    }
  }

  ///@}


  // *****************************************************************
  /// @name Replace long triangles
  // *****************************************************************

  ///@{

  /**
     @brief Replace "long" triangle to max min triangulation angle.
     - Replace only if each long triangle edge is longest edge
       in the adjacent polygon.
     - Replace triangle if it or its long edge adjacent polygons
       have triangulation angle less than angle_threshold.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @pre itriangle is the polygon index of a triangle.
     @pre itriangle is not a deleted polygon.
     @param flag_allow_centroid  If true, allow adding triangulation
         vertices at centroids of adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename PTYPE, typename COS_TYPE, typename MTYPE>
  void replace_long_triangle_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold,
   const PTYPE itriangle,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE THREE(3);
    const DIMENSION_TYPE dimension = mesh.Dimension();

    COS_TYPE cos_min_angle;
    NUMBER_TYPE iangle;
    bool flag_zero;

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord, mesh.PolygonVertexList(itriangle), 
       max_small_magnitude, cos_min_angle, iangle, flag_zero);

    VERTEX_INDEX_TYPE iv_apex;
    HALF_EDGE_INDEX_TYPE ihalf_edge0;
    mesh.GetTriangleVertexAndOpposingEdge
      (itriangle, iangle, iv_apex, ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge1 =
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2 =
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge1X =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2X =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintPolygonVertices(cerr, "  Triangle: ", itriangle, "\n");
    */

    if (!is_longest_polygon_half_edge
        (dimension, vcoord, mesh, ihalf_edge1X)) { return; }
    if (!is_longest_polygon_half_edge
        (dimension, vcoord, mesh, ihalf_edge1X)) { return; }

    const POLYGON_INDEX_TYPE ipolyA =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge1X);
    const POLYGON_INDEX_TYPE ipolyC =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge2X);

    COS_TYPE cos_min_angle_poly3 = 
      std::max(mesh.PolygonCosMinTriangulationAngle(ipolyA),
               cos_min_angle);
    cos_min_angle_poly3 = 
      std::max(cos_min_angle_poly3,
               mesh.PolygonCosMinTriangulationAngle(ipolyC));

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  cos_min_angle_poly3: " << cos_min_angle_poly3
         << "  cos_angle_threshold: " << cos_angle_threshold << endl;
    */


    if (cos_min_angle_poly3 <= cos_angle_threshold) 
      { return; }

    replace_triangleI_to_max_min_triangulation_angle
      (mesh, vcoord, ihalf_edge0, max_small_magnitude, triangulation_settings,
       flag_allow_centroidx2, flag_max_min_all);
  }


  /**
     @brief Replace long triangles to max min triangulation angle.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param flag_allow_centroid  If true, allow adding triangulation
         vertices at centroids of adjacent polygons.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void replace_long_triangles_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      if (!mesh.IsTriangle(ipoly)) { continue; }

      replace_long_triangle_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold, ipoly, max_small_magnitude, 
         triangulation_settings, flag_allow_centroidx2, 
         flag_max_min_all);
    }
  }


  /**
     @brief Replace long triangles to improve max min triangulation angle.
     - Version with array cos_angle_threshold[].
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param flag_allow_centroid  If true, allow adding triangulation
         vertices at centroids of adjacent polygons.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void replace_long_triangles_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      replace_long_triangles_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude, 
         triangulation_settings, flag_allow_centroidx2, 
         flag_max_min_all);
    }
  }

  ///@}


  // *****************************************************************
  /// @name Replace triangles with one split edge
  // *****************************************************************

  ///@{

  /**
     @brief Replace triangle with one split edge to max min triangulation angle.
     @pre All checks have been done to ensure that itriangle
          is a candidate for replacement.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @pre itriangle is the polygon index of a triangle.
     @pre itriangle is not a deleted polygon.
  */
  template <typename MESH_TYPE, typename CTYPE, 
            typename ITYPEA, typename ITYPEB, typename MTYPE>
  bool replace_triangleI_with_one_split_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPEA ibase_half_edge,
   const ITYPEB ilongest_half_edge,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::COS_TYPE COS_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();

    // If true, ilong_half_edge is next half edge of ibase_half_edge.
    const bool flag_is_next_half_edge =
      mesh.IsNextHalfEdge(ibase_half_edge, ilongest_half_edge);
    const HALF_EDGE_INDEX_TYPE isplit_half_edge = 
      mesh.GetIndexOfAdjacentHalfEdgeInPolygon
      (ibase_half_edge, !flag_is_next_half_edge);
    const POLYGON_INDEX_TYPE ipoly0 =
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ilongest_half_edge);
    const POLYGON_INDEX_TYPE itriangle =
      mesh.IndexOfPolygonContainingHalfEdge(ilongest_half_edge);
    const POLYGON_INDEX_TYPE ipoly2 =
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(isplit_half_edge);
    const VERTEX_INDEX_TYPE iv_apex =
      mesh.GetIndexOfHalfEdgeEndpoint
      (ilongest_half_edge, flag_is_next_half_edge);
    const VERTEX_INDEX_TYPE iv_split =
      mesh.GetIndexOfHalfEdgeEndpoint
      (isplit_half_edge, !flag_is_next_half_edge);
    POLYGON_TRIANGULATION_RESULT_E16<COS_TYPE, NUMBER_TYPE> poly_tri_result[3];
    const CTYPE * split_vcoord = &(vcoord[0]) + iv_split*dimension;
    COS_TYPE cos_min_angle;
    bool flag_zero;

    std::vector<CTYPE> ilongest_midpoint_coord(dimension);
    std::vector<CTYPE> midpoint_coord(dimension);
    std::vector<CTYPE> centroid0_coord(dimension);
    std::vector<CTYPE> centroid2_coord(dimension);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintPolygonVertices(cerr, "  Triangle: ", itriangle, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Longest half edge: ", ilongest_half_edge, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Split half edge: ", isplit_half_edge, "\n");
    cerr <<  "  Apex vertex: " << iv_apex << endl;
    */

    if (mesh.IsBoundaryEdge(ilongest_half_edge)) { return(false); }
    if (mesh.IsBoundaryEdge(isplit_half_edge)) { return(false); }

    if (mesh.AreThreeHalfEdgesAroundVertex(iv_apex)) {
      // Polygons ipolyA and ipolyC already share an edge.
      // Replacing the triangle would cause them to share two edges.
      return(false);
    }

    mesh.ComputeEdgeMidpoint
      (vcoord, ilongest_half_edge, ilongest_midpoint_coord);
    IJK::compute_midpoint
      (dimension, IJK::vector2pointer(ilongest_midpoint_coord), split_vcoord, 
       IJK::vector2pointerNC(midpoint_coord));

    mesh.ComputePolygonCentroid(vcoord, ipoly0, centroid0_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly2, centroid2_coord);

    // Replace triangle itriangle with triangle and edge.
    compute_cos_max_min_replace_splitE1T_triangulation_angle
      (dimension, vcoord, midpoint_coord, 
       centroid0_coord, centroid2_coord, 
       mesh, ibase_half_edge, ilongest_half_edge,
       max_small_magnitude, triangulation_settings, 
       poly_tri_result, cos_min_angle, flag_zero);

    const COS_TYPE cos_min_angle_poly3 = 
      mesh.PolygonIIICosMinTriangulationAngle(ipoly0, itriangle, ipoly2);

    if (!flag_zero) {
      if (cos_min_angle < mesh.PolygonCosMinTriangulationAngle(itriangle) ||
          (flag_max_min_all &&
           cos_min_angle < cos_min_angle_poly3)) {

        add_vertex_replace_triangle_with_splitE1_by_triangle_edge
          (mesh, vcoord, midpoint_coord, centroid0_coord, centroid2_coord,
           ibase_half_edge, ilongest_half_edge, poly_tri_result);
        
        return(true);
      }
    }

    return(false);
  }


  /**
     @brief Replace triangle with one split edge to max min triangulation angle.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @pre itriangle is the polygon index of a triangle.
     @pre itriangle is not a deleted polygon.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename PTYPE, typename COS_TYPE, typename MTYPE>
  void replace_triangleI_with_one_split_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold,
   const PTYPE itriangle,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE THREE(3);
    const DIMENSION_TYPE dimension = mesh.Dimension();
    const COS_TYPE cos_min_triangle_angle =
      mesh.PolygonCosMinTriangulationAngle(itriangle);

    const HALF_EDGE_INDEX_TYPE ilongest_half_edge = 
      compute_longest_polygon_half_edge
      (mesh.Dimension(), vcoord, mesh, itriangle);

    const bool flag_next_edge_is_split =
      mesh.IsNextEdgeSplit(ilongest_half_edge);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << " (C)" << endl;
    */

    if (!flag_next_edge_is_split &&
        !mesh.IsPrevEdgeSplit(ilongest_half_edge)) {

      // No split edge.  (Check not necessary, but just in case.)
      return;
    }

    const HALF_EDGE_INDEX_TYPE ibase_half_edge = 
      mesh.GetIndexOfAdjacentHalfEdgeInPolygon
      (ilongest_half_edge, !flag_next_edge_is_split);
    const HALF_EDGE_INDEX_TYPE isplit_half_edge = 
      mesh.GetIndexOfAdjacentHalfEdgeInPolygon
      (ibase_half_edge, !flag_next_edge_is_split);

    if (flag_next_edge_is_split) {
      if (!mesh.AreTwoHalfEdgesAroundFromVertex(isplit_half_edge)) 
        { return; }
    }
    else {
      const HALF_EDGE_INDEX_TYPE isplit_half_edgeX = 
        mesh.IndexOfNextHalfEdgeAroundEdge(isplit_half_edge);

      if (!mesh.AreTwoHalfEdgesAroundFromVertex(isplit_half_edgeX)) 
        { return; }
    }

    const POLYGON_INDEX_TYPE ipolyA =
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ilongest_half_edge);
    const POLYGON_INDEX_TYPE ipolyC =
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(isplit_half_edge);

    if (flag_triangle_has_min_angle) {
      if (cos_min_triangle_angle < 
          mesh.PolygonCosMinTriangulationAngle(ipolyA)) { return; }

      if (cos_min_triangle_angle <
          mesh.PolygonCosMinTriangulationAngle(ipolyC)) { return; }

      if (cos_min_triangle_angle <= cos_angle_threshold)
        { return; } 
    }
    else {
      COS_TYPE cos_min_angle_poly3 = 
        std::max(mesh.PolygonCosMinTriangulationAngle(ipolyA),
                   cos_min_triangle_angle);
      cos_min_angle_poly3 = 
        std::max(cos_min_angle_poly3,
                 mesh.PolygonCosMinTriangulationAngle(ipolyC));

      if (cos_min_angle_poly3 <= cos_angle_threshold) 
        { return; }
    }

    replace_triangleI_with_one_split_edge_to_max_min_triangulation_angle
      (mesh, vcoord, ibase_half_edge, ilongest_half_edge, max_small_magnitude,
       triangulation_settings, flag_max_min_all);
  }


  /**
     @brief Replace triangles with one split edge to max min triangulation angle.
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void replace_triangles_with_one_split_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << " (B)" << endl;
    */

    for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintPolygonIndexAndVertices(cerr, "  Polygon: ", ipoly, "\n");
      cerr << "    Num split vertices: "
           << mesh.NumSplitEdgeVertices(ipoly) << endl;
      */

      if (!mesh.IsTriangleWithOneSplitEdge(ipoly)) { continue; }

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "*** Calling replace triangleT ***" << endl;
      */

      replace_triangleI_with_one_split_edge_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold, ipoly, max_small_magnitude, 
         triangulation_settings,
         flag_triangle_has_min_angle, flag_max_min_all);
    }
  }


  /**
     @brief Replace triangles with one split edge to max min triangulation angle.
     - Version with array cos_angle_threshold[].
     @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
        has been called on all polygons to determine triangulation.
     @param flag_triangle_has_min_angle  If true, replace only if
         min triangle angle is less than min angle in adjacent polygons.
     @param flag_max_min_all If true, replace triangle if it
         improves the min triangulation of the triangle and
         its two neighbors.
       - Otherwise, replace triangle only if it increase the min
         triangulation angle of the triangle.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void replace_triangles_with_one_split_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      replace_triangles_with_one_split_edge_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude, 
         triangulation_settings,
         flag_triangle_has_min_angle, flag_max_min_all);

    }
  }

  ///@}

}


#endif

