/*!
 *  @file ijkmesh2D_triangulateE.tpp
 *  @brief ijk template functions for triangulating by splitting edges.
 *   - Version 0.4.0
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

#ifndef _IJKMESH2D_TRIANGULATE_E_
#define _IJKMESH2D_TRIANGULATE_E_

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkmesh2D_triangulate.tpp"

#include <algorithm>
#include <numeric>
#include <vector>

// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // *****************************************************************
  /// @name Add split vertices, split polygon edges.
  // *****************************************************************

  /// @{

  /// Add split vertex and split polygon edge.
  /// - Allow polygon triangulations using interior vertices.
  template <int BIT_SET_SIZE,
            typename MESH_TYPE, typename CTYPEV, typename CTYPES,
            typename CTYPEI0, typename CTYPEI1,
            typename ITYPEH, 
            typename COS_TYPE2, typename NTYPE2,
            typename ITYPEV, typename ITYPEP>
  void add_vertex_split_polygon_edgeII_allow_IV
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPES> & split_vcoord, 
   const std::vector<CTYPEI0> & interior_coord0,
   const std::vector<CTYPEI1> & interior_coord1, 
   const ITYPEH ihalf_edgeA0, 
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> & polyA0_tri_info,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> & polyA1_tri_info,
   ITYPEV & iv_split, ITYPEP & ipolyB0, ITYPEP & ipolyB1)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;

    static_assert(BIT_SET_SIZE == MESH_TYPE::BitSetSize(), 
                  "Programming error.");

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const HALF_EDGE_INDEX_TYPE ihalf_edgeA1 =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA0);

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "*** Splitting edge ", ihalf_edgeA0, "");
    IJK::print_coord3D(cerr, " at location ", split_vcoord, "\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Polygon 0: ", 
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA0), "\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Polygon 1: ", 
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA1), "\n");
    */

    mesh.MESH_TYPE::MESH2D_SPLIT_II_BASE_TYPE::AddVertexSplitPolygonEdgeII
      (dimension, split_vcoord, vcoord, ihalf_edgeA0, 
       iv_split, ipolyB0, ipolyB1);
    mesh.SetVertexType(iv_split, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);

    mesh.SetNumSplitEdgePolygonVertices(ipolyB0);
    mesh.SetNumSplitEdgePolygonVertices(ipolyB1);

    mesh.CopyTriInfoAndInteriorCoord
      (polyA0_tri_info, interior_coord0, ipolyB0);
    if (ipolyB0 != ipolyB1) {
      mesh.CopyTriInfoAndInteriorCoord
        (polyA1_tri_info, interior_coord1, ipolyB1);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Polygon ", ipolyB0, "\n");
    polyA0_tri_info.Print(cerr, "    ");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Polygon ", ipolyB1, "\n");
    polyA1_tri_info.Print(cerr, "    ");
    */
  }


  /// Add split vertex and split polygon edge.
  /// - Allow polygon triangulations using interior vertices.
  /// - Version that does not return iv_split, ipolyB0, or ipolyB1.
  template <int BIT_SET_SIZE, 
            typename MESH_TYPE, typename CTYPEV, typename CTYPES,
            typename CTYPEI0, typename CTYPEI1, 
            typename ITYPEH, typename COS_TYPE2, typename NTYPE2>
  void add_vertex_split_polygon_edgeII_allow_IV
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPES> & split_vcoord, 
   const std::vector<CTYPEI0> & interior_coord0,
   const std::vector<CTYPEI1> & interior_coord1, 
   const ITYPEH ihalf_edgeA0, 
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> & polyA0_tri_info,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> & polyA1_tri_info)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    VERTEX_INDEX_TYPE iv_split;
    POLYGON_INDEX_TYPE ipolyB[2];

    add_vertex_split_polygon_edgeII_allow_IV
      (mesh, vcoord, split_vcoord, interior_coord0, interior_coord1, 
       ihalf_edgeA0, polyA0_tri_info, polyA1_tri_info, 
       iv_split, ipolyB[0], ipolyB[1]);
  }


  // *****************************************************************
  // Add two split vertices and split two edges of a polygon.
  // *****************************************************************

  /// Add two split vertices and split two polygon edges.
  template <int BIT_SET_SIZE, typename MESH_TYPE, 
            typename CTYPEV, typename CTYPEA, typename CTYPEB,
            typename CTYPEI0, typename CTYPEI1, typename CTYPEI2, 
            typename ITYPEH,
            typename COS_TYPE2, typename NTYPE2,
            typename ITYPEV, typename ITYPEH2>
  void add_vertices_split_two_edges
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord, 
   const std::vector<CTYPEI0> & interior_vcoord0,
   const std::vector<CTYPEI1> & interior_vcoord1, 
   const std::vector<CTYPEI2> & interior_vcoord2, 
   const ITYPEH ihalf_edgeA1,
   const ITYPEH ihalf_edgeB1,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> poly_tri_info[3],
   ITYPEV iv_splitA, ITYPEV & iv_splitB,
   ITYPEH2 & ihalf_edgeA_new,
   ITYPEH2 & ihalf_edgeB_new)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    IJK::PROCEDURE_ERROR error("add_vertices_split_two_edges");

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const NUMBER_TYPE THREE(3);

    if (!mesh.CheckAreAdjacentPolygonsDifferent
        (ihalf_edgeA1, ihalf_edgeB1, error)) 
      { throw error; }

    mesh.AddTwoVerticesSplitTwoEdgesOfPolygon
      (dimension, splitA_vcoord, splitB_vcoord, vcoord,
       ihalf_edgeA1, ihalf_edgeB1, iv_splitA, iv_splitB,
       ihalf_edgeA_new, ihalf_edgeB_new);
    mesh.SetVertexType(iv_splitA, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);
    mesh.SetVertexType(iv_splitB, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);

    POLYGON_INDEX_TYPE ipoly_new[THREE];
    ipoly_new[1] = mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA_new);
    ipoly_new[0] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edgeA_new);
    ipoly_new[2] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edgeB_new);

    mesh.SetNumSplitEdgePolygonVertices(ipoly_new[1]);

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "*** Splitting two edges, ", ihalf_edgeA1, "");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, " and ", ihalf_edgeB1, "\n");
    */

    if (!mesh.IsBoundaryEdge(ihalf_edgeA_new)) {
      mesh.SetNumSplitEdgePolygonVertices(ipoly_new[0]);
      mesh.CopyTriInfoAndInteriorCoord
        (poly_tri_info[0], interior_vcoord0, ipoly_new[0]);
    }
    mesh.CopyTriInfoAndInteriorCoord
      (poly_tri_info[1], interior_vcoord1, ipoly_new[1]);
    if (!mesh.IsBoundaryEdge(ihalf_edgeB_new)) {
      mesh.SetNumSplitEdgePolygonVertices(ipoly_new[2]);
      mesh.CopyTriInfoAndInteriorCoord
        (poly_tri_info[2], interior_vcoord2, ipoly_new[2]);
    }

  }


  /// @brief Add two split vertices and split two polygon edges.
  /// - Version that does not return split vertices or
  ///   new polygon indices.
  template <int BIT_SET_SIZE, typename MESH_TYPE, 
            typename CTYPEV, typename CTYPEA, typename CTYPEB,
            typename CTYPEI0, typename CTYPEI1, typename CTYPEI2, 
            typename ITYPEH,
            typename COS_TYPE2, typename NTYPE2>
  void add_vertices_split_two_edges
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord, 
   const std::vector<CTYPEI0> & interior_vcoord0,
   const std::vector<CTYPEI1> & interior_vcoord1, 
   const std::vector<CTYPEI2> & interior_vcoord2, 
   const ITYPEH ihalf_edgeA1,
   const ITYPEH ihalf_edgeB1,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> poly_tri_info[3])
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;

    VERTEX_INDEX_TYPE iv_splitA, iv_splitB;
    HALF_EDGE_INDEX_TYPE ihalf_edgeA_new, ihalf_edgeB_new;

    add_vertices_split_two_edges
      (mesh, vcoord, splitA_vcoord, splitB_vcoord,
       interior_vcoord0, interior_vcoord1, interior_vcoord2, 
       ihalf_edgeA1, ihalf_edgeB1, poly_tri_info, 
       iv_splitA, iv_splitB, ihalf_edgeA_new, ihalf_edgeB_new);
  }




  // *****************************************************************
  // Add three split vertices and split three edges of a polygon.
  // *****************************************************************

  /**
     @brief Add three split vertices and split three polygon edges.
     @param interior_vcoord0 Coordinates of potential vertex in interior
       of polygon incident on all three split edges.
     @param interior_vcoord1 Coordinates of potential vertex in interior
       of polygon incident on ihalf_edgeA.
     @param interior_vcoord2 Coordinates of potential vertex in interior
       of polygon incident on ihalf_edgeB.
     @param interior_vcoord3 Coordinates of potential vertex in interior
       of polygon incident on ihalf_edgeC.
     @param ihalf_edgeA First split half edge.
     @param ihalf_edgeB Second split half edge.
       @pre ihalf_edgeB must be in same polygon as ihalf_edgeA.
     @param ihalf_edgeC Third split half edge.
       @pre ihalf_edgeC must be in same polygon as ihalf_edgeA.
  */
  template <int BIT_SET_SIZE, typename MESH_TYPE, 
            typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEI0, typename CTYPEI1, 
            typename CTYPEI2, typename CTYPEI3,
            typename ITYPEH,
            typename COS_TYPE2, typename NTYPE2,
            typename ITYPEV, typename ITYPEH2>
  void add_vertices_split_three_edges
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord, 
   const std::vector<CTYPEC> & splitC_vcoord, 
   const std::vector<CTYPEI0> & interior_vcoord0,
   const std::vector<CTYPEI1> & interior_vcoord1, 
   const std::vector<CTYPEI2> & interior_vcoord2, 
   const std::vector<CTYPEI3> & interior_vcoord3, 
   const ITYPEH ihalf_edgeA,
   const ITYPEH ihalf_edgeB,
   const ITYPEH ihalf_edgeC,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> poly_tri_info[4],
   ITYPEV iv_split[3], ITYPEH2 ihalf_edge_new[3])
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    IJK::PROCEDURE_ERROR error("add_vertices_split_three_edges");

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const NUMBER_TYPE THREE(3);

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const POLYGON_INDEX_TYPE iadj_polyA = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeA);
    const POLYGON_INDEX_TYPE iadj_polyB = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeB);
    const POLYGON_INDEX_TYPE iadj_polyC = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeC);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */

    if (!mesh.CheckAreAdjacentPolygonsDifferent
        (ihalf_edgeA, ihalf_edgeB, ihalf_edgeC, error)) 
      { throw error; }

    mesh.AddThreeVerticesSplitThreeEdgesOfPolygon
      (dimension, splitA_vcoord, splitB_vcoord, splitC_vcoord,
       vcoord, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       iv_split, ihalf_edge_new);
    mesh.SetVertexType(iv_split[0], MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);
    mesh.SetVertexType(iv_split[1], MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);
    mesh.SetVertexType(iv_split[2], MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);

    POLYGON_INDEX_TYPE ipoly_new, iadj_poly_new[THREE];
    ipoly_new = mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge_new[0]);
    iadj_poly_new[0] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edge_new[0]);
    iadj_poly_new[1] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edge_new[1]);
    iadj_poly_new[2] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edge_new[2]);


    mesh.SetNumSplitEdgePolygonVertices(ipoly_new);
    mesh.CopyTriInfoAndInteriorCoord
      (poly_tri_info[0], interior_vcoord0, ipoly_new);

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "*** Splitting three edges, ", ihalf_edgeA, "");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, " and ", ihalf_edgeB, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, " and ", ihalf_edgeC, "\n");
    IJK::print_coord3D
      (cerr, "  interior_vcoord1: ", interior_vcoord1, "\n");
    IJK::print_coord3D
      (cerr, "  interior_vcoord2: ", interior_vcoord2, "\n");
    IJK::print_coord3D
      (cerr, "  interior_vcoord3: ", interior_vcoord3, "\n");
    */

    if (!mesh.IsBoundaryEdge(ihalf_edge_new[0])) {
        mesh.SetNumSplitEdgePolygonVertices(iadj_poly_new[0]);
        mesh.CopyTriInfoAndInteriorCoord
          (poly_tri_info[1], interior_vcoord1, iadj_poly_new[0]);
    }

    if (!mesh.IsBoundaryEdge(ihalf_edge_new[1])) {
        mesh.SetNumSplitEdgePolygonVertices(iadj_poly_new[1]);
        mesh.CopyTriInfoAndInteriorCoord
          (poly_tri_info[2], interior_vcoord2, iadj_poly_new[1]);
    }

    if (!mesh.IsBoundaryEdge(ihalf_edge_new[2])) {
        mesh.SetNumSplitEdgePolygonVertices(iadj_poly_new[2]);
        mesh.CopyTriInfoAndInteriorCoord
          (poly_tri_info[3], interior_vcoord3, iadj_poly_new[2]);
    }

  }


  /// @brief Add three split vertices and split three polygon edges.
  /// - Version that does not return split vertices or
  ///   new polygon indices.
  template <int BIT_SET_SIZE, typename MESH_TYPE, 
            typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEI0, typename CTYPEI1, 
            typename CTYPEI2, typename CTYPEI3,
            typename ITYPEH, typename COS_TYPE2, typename NTYPE2>
  void add_vertices_split_three_edges
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord, 
   const std::vector<CTYPEC> & splitC_vcoord, 
   const std::vector<CTYPEI0> & interior_vcoord0,
   const std::vector<CTYPEI1> & interior_vcoord1, 
   const std::vector<CTYPEI2> & interior_vcoord2, 
   const std::vector<CTYPEI3> & interior_vcoord3, 
   const ITYPEH ihalf_edgeA,
   const ITYPEH ihalf_edgeB,
   const ITYPEH ihalf_edgeC,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> poly_tri_info[4])
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;

    const NUMBER_TYPE THREE(3);
    VERTEX_INDEX_TYPE iv_split[THREE];
    HALF_EDGE_INDEX_TYPE ihalf_edge_new[THREE];

    add_vertices_split_three_edges
      (mesh, vcoord, splitA_vcoord, splitB_vcoord, splitC_vcoord,
       interior_vcoord0, interior_vcoord1, 
       interior_vcoord2, interior_vcoord3,
       ihalf_edgeA, ihalf_edgeB, ihalf_edgeC, poly_tri_info, 
       iv_split, ihalf_edge_new);
  }





  // *****************************************************************
  // Add three split vertices and split three polygon edges
  //   that are incident on the same vertex.
  // *****************************************************************

  /// Add three split vertices and split three polygon edges
  ///   that are incident on the same vertex.
  /// - Allow polygon triangulations using interior vertices.
  template <int BIT_SET_SIZE, 
            typename MESH_TYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename ITYPE2, typename ITYPE3, typename ITYPE4, 
            typename ITYPE5, typename ITYPE6,
            typename COS_TYPE2, typename NTYPE2>
  void add_vertices_split_three_incident_edges_allow_IV
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord, 
   const std::vector<CTYPEC> & splitC_vcoord, 
   const std::vector<CTYPE0> & interior_vcoord0,
   const std::vector<CTYPE1> & interior_vcoord1, 
   const std::vector<CTYPE2> & interior_vcoord2, 
   const std::vector<CTYPE3> & interior_vcoord3, 
   const ITYPE2 ihalf_edgeB1,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> poly_tri_info[4],
   ITYPE3 & iv_splitA, ITYPE4 & iv_splitB, ITYPE5 & iv_splitC,
   ITYPE6 ipoly[4])
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC1 =
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edgeB1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC0 =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeC1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB2 =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeA2 =
      mesh.IndexOfPrevHalfEdgeInPolygon(ihalf_edgeB2);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeA3 =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA2);

    // Note: ihalf_edgeA2 and ihalf_edgeC1 are unlinked from mesh
    //   when mesh.AddVerticesSplitThreeIncidentEdgesIV() is called.
    // Must call mesh.IsBoundaryEdge() here, before calling
    //   mesh.AddVerticesSplitThreeIncidentEdgesIV().
    const bool flag_boundary_edgeA = mesh.IsBoundaryEdge(ihalf_edgeA2);
    const bool flag_boundary_edgeC = mesh.IsBoundaryEdge(ihalf_edgeC1);

    // *** DEBUG ***
    /*
    using namespace std;
    int iloc;
    const POLYGON_INDEX_TYPE ipoly1 = 
    mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB1);
    const POLYGON_INDEX_TYPE ipoly2 = 
    mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB2);
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "*** Splitting three incident edges, ", ihalf_edgeC1, "");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, " and ", ihalf_edgeB2, "");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, " and ", ihalf_edgeA3, "\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Splitting two adjacent edges of polygon ", ipoly1, "\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Splitting two adjacent edges of polygon ", ipoly2, "\n");
    poly_tri_info[0].PrintEar(cerr, "  ear 0: ", "\n");
    poly_tri_info[1].PrintEar(cerr, "  ear 1: ", "\n");
    poly_tri_info[2].PrintEar(cerr, "  ear 2: ", "\n");
    poly_tri_info[3].PrintEar(cerr, "  ear 3: ", "\n");
    */
    
    mesh.AddVerticesSplitThreeIncidentEdgesIV
      (dimension, splitA_vcoord, splitB_vcoord, splitC_vcoord,
       vcoord, ihalf_edgeB1, iv_splitA, iv_splitB, iv_splitC, ipoly);
    mesh.SetVertexType(iv_splitA, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);
    mesh.SetVertexType(iv_splitB, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);
    mesh.SetVertexType(iv_splitC, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);

    if (!flag_boundary_edgeC)
      { mesh.SetNumSplitEdgePolygonVertices(ipoly[0]); }
    mesh.SetNumSplitEdgePolygonVertices(ipoly[1]);
    mesh.SetNumSplitEdgePolygonVertices(ipoly[2]);
    if (!flag_boundary_edgeA)
      { mesh.SetNumSplitEdgePolygonVertices(ipoly[3]); }

    if (!flag_boundary_edgeC) {
      mesh.CopyTriInfoAndInteriorCoord
        (poly_tri_info[0], interior_vcoord0, ipoly[0]);
    }
    mesh.CopyTriInfoAndInteriorCoord
      (poly_tri_info[1], interior_vcoord1, ipoly[1]);
    mesh.CopyTriInfoAndInteriorCoord
      (poly_tri_info[2], interior_vcoord2, ipoly[2]);


    if (!flag_boundary_edgeA) {
      mesh.CopyTriInfoAndInteriorCoord
        (poly_tri_info[3], interior_vcoord3, ipoly[3]);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintPolygonIndexAndVertices(cerr, "  New poly0: ", ipoly[0], "\n");
    mesh.PrintPolygonIndexAndVertices(cerr, "  New poly1: ", ipoly[1], "\n");
    mesh.PrintPolygonIndexAndVertices(cerr, "  New poly2: ", ipoly[2], "\n");
    mesh.PrintPolygonIndexAndVertices(cerr, "  New poly3: ", ipoly[3], "\n");
    */
  }


  /// Add three split vertices and split three polygon edges
  ///   that are incident on the same vertex.
  /// - Allow polygon triangulations using interior vertices.
  /// - Version that does not return iv_split, ipoly[].
  template <int BIT_SET_SIZE, 
            typename MESH_TYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename ITYPE2,
            typename COS_TYPE2, typename NTYPE2>
  void add_vertices_split_three_incident_edges_allow_IV
  (MESH_TYPE & mesh, std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord, 
   const std::vector<CTYPEC> & splitC_vcoord, 
   const std::vector<CTYPE0> & interior_vcoord0,
   const std::vector<CTYPE1> & interior_vcoord1, 
   const std::vector<CTYPE2> & interior_vcoord2, 
   const std::vector<CTYPE3> & interior_vcoord3, 
   const ITYPE2 ihalf_edgeB1,
   POLYGON_TRIANGULATION_RESULT_E
   <BIT_SET_SIZE, COS_TYPE2, NTYPE2> poly_tri_info[4])
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE FOUR_POLYGONS(4);
    VERTEX_INDEX_TYPE iv_splitA, iv_splitB, iv_splitC;
    POLYGON_INDEX_TYPE ipoly[FOUR_POLYGONS];

    add_vertices_split_three_incident_edges_allow_IV
      (mesh, vcoord, splitA_vcoord, splitB_vcoord, splitC_vcoord,
       interior_vcoord0, interior_vcoord1, 
       interior_vcoord2, interior_vcoord3,
       ihalf_edgeB1, poly_tri_info, 
       iv_splitA, iv_splitB, iv_splitC, ipoly);
  }

  ///@}


  // *****************************************************************
  /// @name Split boundary edges to max min triangulation.
  // *****************************************************************

  /// @{

  /*!
   *  @brief Split boundary edge to improve max min triangulation.
   *  - Return true if edge is split.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE>
  bool split_boundary_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edgeA, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const CTYPE epsilon_cos = 0.001)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::COS_TYPE COS_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const POLYGON_INDEX_TYPE ipoly0 =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    std::vector<CTYPE>  midpoint_coord(dimension);
    POLYGON_TRIANGULATION_RESULT_E16<COS_TYPE, ITYPE> poly_tri_result;

    VERTEX_INDEX_TYPE iv_split;
    COS_TYPE cos_min_angle;
    bool flag_zero;
    std::vector<CTYPE> centroid0_coord(dimension);

    /* OBSOLETE. NO REASON TO CALL THIS ON A BOUNDARY EDGE.
    if (mesh.DoPolygonsShareThreeVertices(ihalf_edgeA)) { return(false); }
    */

    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeA, midpoint_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly0, centroid0_coord);

    compute_split_half_edge_triangulation_to_max_min_angleP
      (dimension, vcoord, midpoint_coord, centroid0_coord,
       mesh, ihalf_edgeA, max_small_magnitude, 
       triangulation0_settings, poly_tri_result);

    if (!flag_zero) {
      if (poly_tri_result.cos_min_triangulation_angle < 
          mesh.PolygonCosMinTriangulationAngle(ipoly0)) {

        const VERTEX_INDEX_TYPE iv_split =
          mesh.AddSplitVertex(midpoint_coord, vcoord);

        POLYGON_INDEX_TYPE ipoly_new0;
        POLYGON_INDEX_TYPE ipoly_new1;
        mesh.SplitPolygonEdgeII
          (ihalf_edgeA, iv_split, ipoly_new0, ipoly_new1);

        mesh.CopyTriInfoAndInteriorCoord
          (poly_tri_result, centroid0_coord, ipoly_new0);

        return(true);
      }
    }

    return(false);
  }

  ///@}


  // *****************************************************************
  /// @name Split interior edges to max min triangulation.
  // *****************************************************************

  ///@{

  /*!
   * Split interior edge to improve max min triangulation.
   * - Return true if edge is split.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  @pre ihalf_edge0 is an interior edge.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if min angle with splitting
   *     is larger than the min angle of polygon ipoly0.
   *  @param flag_must_improve_smallest_polygon_angles If true, only split edge
   *     if splitting improves smallest polygon angle.
   *  - If both flag_must_improve_both_polygon_angles and 
   *     flag_must_improve_smallest_polygon_angles are false,
   *     then will split if improves min angle in polygon ipoly0 and does
   *     not (substantially) decrease the min angle in polygon ipoly1, 
   *     even if min angle in polygon ipoly1 is less than min angle in polygon ipoly0.
   *  @param epsilon_cos If new cos of min angle ipoly1 is at least
   *     epsilon_cos - PolygonCosMinTriangulationAngle(ipoly1),
   *     then the new triangulation does not substantially decrease 
   *     the min angle in polygon ipoly1.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE>
  bool split_interior_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edge0, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const bool flag_must_improve_both_polygon_angles,
   const bool flag_must_improve_smaller_polygon_angle = true,
   const CTYPE epsilon_cos = 0.001)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::COS_TYPE COS_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const POLYGON_INDEX_TYPE ipoly0 =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge1 =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge0);
    const POLYGON_INDEX_TYPE ipoly1 =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge1);
    std::vector<CTYPE>  midpoint_coord(dimension);
    POLYGON_TRIANGULATION_RESULT_E16<COS_TYPE, ITYPE> poly_tri_result[2];
    ITYPE tri_vertex0_index;
    COS_TYPE cosA, cosB;

    VERTEX_INDEX_TYPE iv_split;
    COS_TYPE cos_min_angle;
    bool flag_zero;
    std::vector<CTYPE> centroid0_coord(dimension);
    std::vector<CTYPE> centroid1_coord(dimension);

    if (mesh.DoPolygonsShareThreeVertices(ihalf_edge0)) { return(false); }

    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edge0, midpoint_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly0, centroid0_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly1, centroid1_coord);

    compute_split_edge_triangulation_to_max_min_angle
      (dimension, vcoord, midpoint_coord, 
       centroid0_coord, centroid1_coord, mesh,
       ihalf_edge0, 
       max_small_magnitude,
       triangulation0_settings, triangulation1_settings,
       poly_tri_result[0], poly_tri_result[1],
       cos_min_angle, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge: ", ihalf_edge0, "\n");
    mesh.PrintPolygonIndexAndVertices(cerr, "  ipoly0: ", ipoly0, "\n");
    mesh.PolygonTriangulationInfo(ipoly0).Print(cerr, "    ");
    mesh.PrintPolygonIndexAndVertices(cerr, "  ipoly1: ", ipoly1, "\n");
    mesh.PolygonTriangulationInfo(ipoly1).Print(cerr, "    ");
    cerr << "  flag_zero: " << int(flag_zero) << endl;
    cerr << "  cos_min_angle: " << cos_min_angle
         << "  angle: " << std::acos(cos_min_angle)*180.0/M_PI << endl;    
    cerr << "  poly_tri_result[0]: " << endl;;
    poly_tri_result[0].Print(cerr, "    ");
    cerr << "  poly_tri_result[1]: " << endl;
    poly_tri_result[1].Print(cerr, "    ");
    cerr << "  epsilon_cos: " << epsilon_cos << endl;
    */

    if (!flag_zero) {

      if (flag_must_improve_smaller_polygon_angle ||
          flag_must_improve_both_polygon_angles) {
        if (cos_min_angle < mesh.PolygonCosMinTriangulationAngle(ipoly0) - epsilon_cos) {
          if (!flag_must_improve_both_polygon_angles ||
              cos_min_angle < mesh.PolygonCosMinTriangulationAngle(ipoly1) - epsilon_cos) {

            // *** DEBUG ***
	    /*
            using namespace std;
            cerr << "In " << __func__ << endl;
            mesh.PrintHalfEdgeIndexAndEndpoints
              (cerr, "  *** Splitting half edge (A): ", ihalf_edge0, "\n");
	    */

            add_vertex_split_polygon_edgeII_allow_IV
              (mesh, vcoord, midpoint_coord, 
               centroid0_coord, centroid1_coord,
               ihalf_edge0, poly_tri_result[0], poly_tri_result[1]);

            return(true);
          }
        }
      }
      else {
        if ((poly_tri_result[0].cos_min_triangulation_angle <
             mesh.PolygonCosMinTriangulationAngle(ipoly0) - epsilon_cos) &&
            (cos_min_angle <
             mesh.PolygonCosMinTriangulationAngle(ipoly1) - epsilon_cos)) {

          // *** DEBUG ***
          /*
          using namespace std;
            cerr << "In " << __func__ << endl;
          mesh.PrintHalfEdgeIndexAndEndpoints
            (cerr, "  *** Splitting half edge (B): ", ihalf_edge0, "\n");
          */
        
          add_vertex_split_polygon_edgeII_allow_IV
            (mesh, vcoord, midpoint_coord, 
             centroid0_coord, centroid1_coord,
             ihalf_edge0, poly_tri_result[0], poly_tri_result[1]);

          return(true);
        }
      }
    }
    return(false);
  }

///@}


  // *****************************************************************
  /// @name Split edges to max min triangulation.
  // *****************************************************************

  ///@{

  /*!
   *  Split edge to improve max min triangulation.
   *  - Return true if edge is split.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
  */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE>
  bool split_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edge,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const bool flag_must_improve_both_polygon_angles,
   const bool flag_must_improve_smaller_polygon_angle = true,
   const CTYPE epsilon_cos = 0.001)
  {
    bool flag;

    if (mesh.IsBoundaryEdge(ihalf_edge)) { 
      flag = 
        split_boundary_edge_to_max_min_triangulation_angle
        (mesh, vcoord, ihalf_edge, max_small_magnitude, 
         triangulation0_settings, epsilon_cos);
    }
    else {
      flag =
        split_interior_edge_to_max_min_triangulation_angle
        (mesh, vcoord, ihalf_edge, max_small_magnitude, 
         triangulation0_settings, triangulation1_settings,
         flag_must_improve_both_polygon_angles,
         flag_must_improve_smaller_polygon_angle,
         epsilon_cos);
    }
        
    return(flag);
  }


  /*!
   *  Split longest polygon edge to improve max min triangulation.
   *  - Return true if edge is split.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear and 
   *    has an interior triangulation vertex.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE>
  bool split_longest_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ipoly, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const bool flag_must_improve_both_polygon_angles)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;

    bool flag;

    const HALF_EDGE_INDEX_TYPE ilongest_half_edge = 
      compute_longest_polygon_half_edge
      (mesh.Dimension(), vcoord, mesh, ipoly);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Splitting longest edge of polygon ", ipoly, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "    Longest half edge: ", ilongest_half_edge, "\n");
    */

    flag = split_edge_to_max_min_triangulation_angle
        (mesh, vcoord, ilongest_half_edge, max_small_magnitude, 
         triangulation0_settings, triangulation1_settings,
         flag_must_improve_both_polygon_angles);

    return(flag);
  }


  /*!
   *  Split longest polygon edges to maximize minimum triangulation angle.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear and 
   *    has an interior triangulation vertex.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split longest polygon edge
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void split_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const bool flag_must_improve_both_polygon_angles)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (NUMBER_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      if (mesh.PolygonCosMinTriangulationAngle(ipoly) > cos_angle_threshold) {
        split_longest_edge_to_max_min_triangulation_angle
          (mesh, vcoord, ipoly, max_small_magnitude, 
           triangulation0_settings, triangulation1_settings,
           flag_must_improve_both_polygon_angles);
      }
    }

  }


  /*!
   *  Split longest polygon edges to maximize minimum triangulation angle.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  - Version with array cos_angle_threshold[].
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split longest polygon edge
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void split_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const bool flag_must_improve_both_polygon_angles)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      split_longest_polygon_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude,
         triangulation0_settings, triangulation1_settings,
         flag_must_improve_both_polygon_angles);
    }
  }


  /*!
   *  Split longest or second longest polygon edge to maximize
   *    minimum triangulation angle.
   *  - Return true if some edge is split.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE, typename RTYPE>
  bool split_longest_or_second_longest_edge_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ipoly, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const RTYPE R,
   const bool flag_must_improve_both_polygon_angles)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;

    const RTYPE Rsquared = R*R;
    HALF_EDGE_INDEX_TYPE ilongest_half_edge, isecond_longest_half_edge;
    CTYPE longest_length_squared;
    CTYPE second_longest_length_squared;
    CTYPE shortest_length_squared, L;

    bool flag;

    compute_two_longest_polygon_half_edges
      (mesh.Dimension(), vcoord, mesh, ipoly, 
       ilongest_half_edge, isecond_longest_half_edge,
       longest_length_squared, second_longest_length_squared);

    mesh.ComputeShortestHalfEdge(vcoord, ipoly, shortest_length_squared);
    L = shortest_length_squared*Rsquared;

    flag = split_edge_to_max_min_triangulation_angle
      (mesh, vcoord, ilongest_half_edge, max_small_magnitude, 
       triangulation0_settings, triangulation1_settings, 
       flag_must_improve_both_polygon_angles);

    if (flag) { 
      // ilongest_half_edge was split.
      return(true);
    }

    if (second_longest_length_squared < L) {
      // Second longest edge is not much longer than shortest.  
      //    Don't try to split.
      return(false);
    }
    

    flag = split_edge_to_max_min_triangulation_angle
      (mesh, vcoord, isecond_longest_half_edge, max_small_magnitude, 
       triangulation0_settings, triangulation1_settings,
       flag_must_improve_both_polygon_angles);
        
    return(flag);
  }


  /*!
   *  Split longest or second longest polygon edges to maximize 
   *    minimum triangulation angle.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split longest polygon edge
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_longest_or_second_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const RTYPE R, const bool flag_must_improve_both_polygon_angles)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (NUMBER_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      if (mesh.PolygonCosMinTriangulationAngle(ipoly) > cos_angle_threshold) {
        split_longest_or_second_longest_edge_to_max_min_triangulation_angle
          (mesh, vcoord, ipoly, max_small_magnitude, 
           triangulation0_settings, triangulation1_settings, 
           R, flag_must_improve_both_polygon_angles);
      }
    }

  }


  /*!
   *  Split longest or second longest polygon edges to maximize minimum triangulation angle.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  - Version with array cos_angle_threshold[].
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split longest polygon edge
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_longest_or_second_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const RTYPE R,
   const bool flag_must_improve_both_polygon_angles)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      split_longest_or_second_longest_polygon_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude, 
         triangulation0_settings, triangulation1_settings,
         R, flag_must_improve_both_polygon_angles);
    }
  }


  /*!
   *  Split long polygon edges to improve max min triangulation.
   *  - Return true if some edge is split.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   *  @param flag_must_improve_both_polygon_angles If true, only split edge
   *     if splitting improves min angles in both polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE, typename RTYPE>
  bool split_long_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ipoly, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const RTYPE R,
   const bool flag_must_improve_both_polygon_angles,
   const bool flag_must_improve_smaller_polygon_angle = true,
   const CTYPE epsilon_cos = 0.001)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;

    const int THREE(3);
    const RTYPE Rsquared = R*R;
    HALF_EDGE_INDEX_TYPE ihalf_edge[THREE];
    CTYPE length_squared[THREE];
    CTYPE shortest_length_squared;

    // Stores longest edges in decreasing order by length.
    compute_three_longest_polygon_half_edges
      (mesh.Dimension(), vcoord, mesh, ipoly, ihalf_edge, length_squared);

    mesh.ComputeShortestHalfEdge(vcoord, ipoly, shortest_length_squared);
    const CTYPE L = shortest_length_squared*Rsquared;

    for (int j = 0; j < THREE; j++) {
      if (length_squared[j] >= L) {

        // *** DEBUG ***
        /*
        using namespace std;
        mesh.PrintHalfEdgeIndexAndEndpoints
          (cerr, "Attempting to split half edge ", ihalf_edge[j], "\n");
        */

        if (split_edge_to_max_min_triangulation_angle
            (mesh, vcoord, ihalf_edge[j], max_small_magnitude, 
             triangulation0_settings, triangulation1_settings, 
             flag_must_improve_both_polygon_angles,
             flag_must_improve_smaller_polygon_angle, epsilon_cos)) {

          // Edge was split.
          return true;
        }
      }
      else {
        // All remaining edges are shorter than sqrt(L).
        return false;
      }
    }

    // No edge split.
    return(false);
  }

  ///@}


  // *****************************************************************
  /// @name Split two longest edges to max min triangulation.
  // *****************************************************************

  ///@{

  /*!
   *  Split two polygon edges to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param triangulation_settings[i] Indicates which types of triangulations 
   *     to consider for i'th polygon.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE, typename ITYPE,
            typename MTYPE>
  bool split_two_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edgeA1, const ITYPE ihalf_edgeB1,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::COS_TYPE COS_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const NUMBER_TYPE THREE_POLYGONS(3);
    POLYGON_INDEX_TYPE ipoly[THREE_POLYGONS];
    POLYGON_TRIANGULATION_RESULT_E16<COS_TYPE, ITYPE> 
      poly_tri_result[THREE_POLYGONS];
    bool flag_split_edge[2];
    COS_TYPE cos_min_angle;
    bool flag_zero;
    IJK::PROCEDURE_ERROR error
      ("split_two_edges_to_max_min_triangulation_angle");

    if (!mesh.CheckAreHalfEdgesInSamePolygon
        (ihalf_edgeA1, ihalf_edgeB1, error)) 
      { throw error; }

    ipoly[0] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edgeA1);
    ipoly[1] = mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA1);
    ipoly[2] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edgeB1);

    if (!mesh.CheckAreAdjacentPolygonsDifferent
        (ihalf_edgeA1, ihalf_edgeB1, error)) {
      // ipoly[1] shares ihalf_edgeA1 and ihalf_edgeB1 with ipoly[0].
      // Don't split.

      return(false);
    }

    std::vector<CTYPE> midpointA_coord(dimension);
    std::vector<CTYPE> midpointB_coord(dimension);
    std::vector<CTYPE> centroid0_coord(dimension);
    std::vector<CTYPE> centroid1_coord(dimension);
    std::vector<CTYPE> centroid2_coord(dimension);

    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeA1, midpointA_coord);
    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeB1, midpointB_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[0], centroid0_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[1], centroid1_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[2], centroid2_coord);

    compute_cos_max_min_split_1or2_edges_polyIII_triangulation_angle
      (dimension, vcoord, midpointA_coord, midpointB_coord,
       centroid0_coord, centroid1_coord, centroid2_coord,
       mesh, ihalf_edgeA1, ihalf_edgeB1, max_small_magnitude, 
       triangulation_settings, poly_tri_result, 
       flag_split_edge, cos_min_angle, flag_zero);

    if (flag_zero) { return(false); }

    if (!flag_split_edge[0] || !flag_split_edge[1]) { 
      // Only split if two edge split maximizes min triangulation angle.
      return(false); 
    }

    const COS_TYPE cos_min_angle_poly3 = 
      mesh.PolygonIIICosMinTriangulationAngle(ipoly[0], ipoly[1], ipoly[2]);


    if (cos_min_angle < mesh.PolygonCosMinTriangulationAngle(ipoly[1]) ||
        (flag_max_min_all &&
         cos_min_angle < cos_min_angle_poly3)) {

      add_vertices_split_two_edges
        (mesh, vcoord, 
         midpointA_coord, midpointB_coord,
         centroid0_coord, centroid1_coord, centroid2_coord,
         ihalf_edgeA1, ihalf_edgeB1, poly_tri_result);

      return(true);
    }

    return(false);
  }


  /*!
   *  Split two longest polygon edges to maximize minimum triangulation angle.
   *  - Return true if some edge is split.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation_settings[i] Indicates which types of triangulations 
   *     to consider for i'th polygon.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE, typename COS_TYPE,
            typename ITYPE, typename MTYPE, typename RTYPE>
  bool split_two_longest_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const ITYPE ipoly1, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const RTYPE R,
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const RTYPE Rsquared = R*R;
    HALF_EDGE_INDEX_TYPE ilongest_half_edge, isecond_longest_half_edge;
    CTYPE longest_length_squared;
    CTYPE second_longest_length_squared;
    CTYPE shortest_length_squared, L;

    compute_two_longest_polygon_half_edges
      (mesh.Dimension(), vcoord, mesh, ipoly1, 
       ilongest_half_edge, isecond_longest_half_edge,
       longest_length_squared, second_longest_length_squared);

    mesh.ComputeShortestHalfEdge(vcoord, ipoly1, shortest_length_squared);
    L = shortest_length_squared*Rsquared;

    if (second_longest_length_squared < L) {
      // Second longest edge is not much longer than shortest.  
      //    Don't try to split.
      return(false);
    }

    const POLYGON_INDEX_TYPE ipoly0 = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ilongest_half_edge);
    const POLYGON_INDEX_TYPE ipoly2 = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (isecond_longest_half_edge);

    if (mesh.PolygonCosMinTriangulationAngle(ipoly1) < cos_angle_threshold) {

      if (flag_poly_angle_below_threshold) {
        // Min angle in ipoly1 is above threshold. Don't attempt split.
        return(false);
      }
      else {
        // Attempt to split two edges if adjacent polygon angles 
        //   are above threshold.

        if (mesh.PolygonIICosMinTriangulationAngle(ipoly0, ipoly2) <
            cos_angle_threshold) {
          // Min angles in all three polygons are above threshold. 
          // Don't attempt split.
          return(false);
        }
      }
    }

    const bool flag = split_two_edges_to_max_min_triangulation_angle
      (mesh, vcoord, ilongest_half_edge, isecond_longest_half_edge,
       max_small_magnitude, triangulation_settings, flag_max_min_all);

    return(flag);
  }


  /*!
   *  Split two longest polygon edges to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation_settings[i] Indicates which types of triangulations 
   *     to consider for i'th polygon.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_two_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const RTYPE R, 
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (NUMBER_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintPolygonIndexAndVertices
        (cerr, "  Polygon: ", ipoly, "");
      cerr << "  Num mesh vertices: " << mesh.NumVertices() << endl;
      */

      split_two_longest_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold, ipoly, max_small_magnitude, 
         triangulation_settings, R, 
         flag_poly_angle_below_threshold, flag_max_min_all);
    }

  }


  /*!
   *  Split two longest polygon edges to maximize minimum triangulation angle.
   *  - Version with array cos_angle_threshold[].
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation_settings[i] Indicates which types of triangulations 
   *     to consider for i'th polygon.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_two_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const RTYPE R,
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      split_two_longest_polygon_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude, 
         triangulation_settings, R, 
         flag_poly_angle_below_threshold, flag_max_min_all);
    }
  }

  ///@}


  // *****************************************************************
  /// @name Split three longest edges to max min triangulation.
  // *****************************************************************

  ///@{

  /*!
   *  Split three polygon edges to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param inner_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygon with 3 split edges.
   *  @param outer_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygons with 1 split edge.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE, typename ITYPE,
            typename MTYPE>
  bool split_three_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edgeA, const ITYPE ihalf_edgeB,
   const ITYPE ihalf_edgeC,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::COS_TYPE COS_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const NUMBER_TYPE THREE_SPLIT_EDGES(3);
    const NUMBER_TYPE FOUR_POLYGONS(4);
    POLYGON_INDEX_TYPE ipoly[FOUR_POLYGONS];
    POLYGON_TRIANGULATION_RESULT_E16<COS_TYPE, ITYPE> 
      poly_tri_result[FOUR_POLYGONS];
    bool flag_split_edge[THREE_SPLIT_EDGES];
    COS_TYPE cos_min_angle;
    bool flag_zero;
    IJK::PROCEDURE_ERROR error
      ("split_three_edges_to_max_min_triangulation_angle");

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */

    if (!mesh.CheckAreHalfEdgesInSamePolygon
        (ihalf_edgeA, ihalf_edgeB, ihalf_edgeC, error)) 
      { throw error; }
    
    ipoly[0] = mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    ipoly[1] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edgeA);
    ipoly[2] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edgeB);
    ipoly[3] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (ihalf_edgeC);

    if (!mesh.CheckAreAdjacentPolygonsDifferent
        (ihalf_edgeA, ihalf_edgeB, ihalf_edgeC, error)) {
      // Some polygon shares two edges with ipoly[0].
      // Don't split.

      return(false);
    }

    std::vector<CTYPE> midpointA_coord(dimension);
    std::vector<CTYPE> midpointB_coord(dimension);
    std::vector<CTYPE> midpointC_coord(dimension);
    std::vector<CTYPE> centroid0_coord(dimension);
    std::vector<CTYPE> centroid1_coord(dimension);
    std::vector<CTYPE> centroid2_coord(dimension);
    std::vector<CTYPE> centroid3_coord(dimension);

    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeA, midpointA_coord);
    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeB, midpointB_coord);
    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeC, midpointC_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[0], centroid0_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[1], centroid1_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[2], centroid2_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[3], centroid3_coord);

    /* OBSOLETE
    compute_cos_max_min_split_exactly_three_edges_polyIV_triangulation_angle
    */
    compute_split_three_edges_triangulation_to_max_min_angle
      (dimension, vcoord, 
       midpointA_coord, midpointB_coord, midpointC_coord,
       centroid0_coord, centroid1_coord, 
       centroid2_coord, centroid3_coord,
       mesh, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       max_small_magnitude, 
       inner_poly_triangulation_settings,
       outer_poly_triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);

    if (flag_zero) { return(false); }

    const COS_TYPE cos_min_angle_poly4 = 
      mesh.PolygonIVCosMinTriangulationAngle
      (ipoly[0], ipoly[1], ipoly[2], ipoly[3]);


    if (cos_min_angle < mesh.PolygonCosMinTriangulationAngle(ipoly[0]) ||
        (flag_max_min_all &&
         cos_min_angle < cos_min_angle_poly4)) {

      add_vertices_split_three_edges
        (mesh, vcoord, 
         midpointA_coord, midpointB_coord, midpointC_coord,
         centroid0_coord, centroid1_coord, 
         centroid2_coord, centroid3_coord,
         ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
         poly_tri_result);

      return(true);
    }

    return(false);
  }


  /*!
   *  Split three longest polygon edges to maximize minimum triangulation angle.
   *  - Return true if some edge is split.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param inner_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygon with 3 split edges.
   *  @param outer_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygons with 1 split edge.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE, typename COS_TYPE,
            typename ITYPE, typename MTYPE, typename RTYPE>
  bool split_three_longest_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const ITYPE ipoly, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const RTYPE R,
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE THREE_SPLIT_EDGES(3);
    const RTYPE Rsquared = R*R;
    HALF_EDGE_INDEX_TYPE ihalf_edge_long[THREE_SPLIT_EDGES];
    CTYPE length_squared[THREE_SPLIT_EDGES];
    POLYGON_INDEX_TYPE iadj_poly[THREE_SPLIT_EDGES];
    CTYPE shortest_length_squared, L;

    compute_three_longest_polygon_half_edges
      (mesh.Dimension(), vcoord, mesh, ipoly, 
       ihalf_edge_long, length_squared);

    mesh.ComputeShortestHalfEdge(vcoord, ipoly, shortest_length_squared);
    L = shortest_length_squared*Rsquared;

    if (length_squared[THREE_SPLIT_EDGES-1] < L) {
      // Third longest edge is not much longer than shortest.  
      //    Don't try to split.
      return(false);
    }

    for (NUMBER_TYPE j = 0; j < THREE_SPLIT_EDGES; j++) {
      iadj_poly[j] = 
        mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge
        (ihalf_edge_long[j]);
    }

    if (mesh.PolygonCosMinTriangulationAngle(ipoly) < cos_angle_threshold) {

      if (flag_poly_angle_below_threshold) {
        // Min angle in ipoly1 is above threshold. Don't attempt split.
        return(false);
      }
      else {
        // Attempt to split three edges if adjacent polygon angles 
        //   are above threshold.

        if (mesh.PolygonIIICosMinTriangulationAngle(iadj_poly)
            < cos_angle_threshold) {
          // Min angles in all three polygons are above threshold. 
          // Don't attempt split.
          return(false);
        }
      }
    }

    const bool flag = split_three_edges_to_max_min_triangulation_angle
      (mesh, vcoord, ihalf_edge_long[0], ihalf_edge_long[1],
       ihalf_edge_long[2], max_small_magnitude, 
       inner_poly_triangulation_settings, 
       outer_poly_triangulation_settings, 
       flag_max_min_all);

    return(flag);
  }


  /*!
   *  Split three longest polygon edges to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param inner_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygon with 3 split edges.
   *  @param outer_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygons with 1 split edge.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_three_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const RTYPE R, 
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (NUMBER_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintPolygonIndexAndVertices
        (cerr, "  Polygon: ", ipoly, "");
      cerr << "  Num mesh vertices: " << mesh.NumVertices() << endl;
      */

      split_three_longest_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold, ipoly, max_small_magnitude, 
         inner_poly_triangulation_settings, 
         outer_poly_triangulation_settings, 
         R, flag_poly_angle_below_threshold, flag_max_min_all);
    }

  }


  /*!
   *  Split three longest polygon edges to maximize minimum triangulation angle.
   *  - Version with array cos_angle_threshold[].
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param inner_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygon with 3 split edges.
   *  @param outer_poly_triangulation_settings Indicates types 
   *     of triangulations to consider for polygons with 1 split edge.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *  @param flag_max_min_all If true, only split edges
   *     if splitting improves min angles in all polygons.
   *     - Otherwise, split edge if splitting improves min angle
   *     in polygon containing half edge ihalf_edge.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_three_longest_polygon_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const RTYPE R,
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      split_three_longest_polygon_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude, 
         inner_poly_triangulation_settings, 
         outer_poly_triangulation_settings, 
         R, flag_poly_angle_below_threshold, flag_max_min_all);
    }
  }

  ///@}


  // *****************************************************************
  /// @name Split three triangle edges to max min triangulation.
  // *****************************************************************

  ///@{


  /*!
   *  Split three triangle edges to maximize minimum triangulation angle.
   *  - Return true if some edge is split.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   */
  template <typename MESH_TYPE, typename CTYPE, typename COS_TYPE,
            typename ITYPE, typename MTYPE>
  bool split_three_triangle_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const ITYPE ipoly, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE THREE(3);
    HALF_EDGE_INDEX_TYPE ihalf_edge[THREE];
    POLYGON_INDEX_TYPE iadj_poly[THREE];
    POLYGON_TRIANGULATION_SETTINGS triangle_triangulation_settings;

    if (!mesh.IsTriangle(ipoly)) {
      // Polygon ipoly is not a triangle.
      return(false);
    }

    triangle_triangulation_settings.flag_all_triangulations_from_polygon_vertices = true;
    triangle_triangulation_settings.flag_all_triangulations_from_hexagon_vertices = true;
    triangle_triangulation_settings.SetV(true);

    for (NUMBER_TYPE j = 0; j < THREE; j++) {
      ihalf_edge[j] = mesh.HalfEdgeIndex(ipoly, j);
      iadj_poly[j] = 
        mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge[j]);
    }

    // Attempt to split three edges if adjacent polygon angles 
    //   are below threshold.

    if (mesh.PolygonIIICosMinTriangulationAngle(iadj_poly)
        < cos_angle_threshold) {
      // Min angles in all three polygons are above threshold. 
      // Don't attempt split.
      return(false);
    }

    const bool flag = split_three_edges_to_max_min_triangulation_angle
      (mesh, vcoord, ihalf_edge[0], ihalf_edge[1],
       ihalf_edge[2], max_small_magnitude, 
       triangle_triangulation_settings, 
       outer_poly_triangulation_settings,
       true);

    return(flag);
  }



  /*!
   *  Split three triangle edges to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void split_three_triangle_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (NUMBER_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintPolygonIndexAndVertices
        (cerr, "  Polygon: ", ipoly, "");
      cerr << "  Num mesh vertices: " << mesh.NumVertices() << endl;
      */

      split_three_triangle_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold, ipoly, max_small_magnitude, 
         outer_poly_triangulation_settings);
    }

  }


  /*!
   *  Split three triangle edges to maximize minimum triangulation angle.
   *  - Version with array cos_angle_threshold[].
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split edges
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param outer_poly_triangulation_settings 
   *     Indicates which types of triangulations to consider 
   *     for adjacent polygons.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void split_three_triangle_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      split_three_triangle_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude, 
         outer_poly_triangulation_settings);
    }
  }

  ///@}


  // *****************************************************************
  /// @name Split long polygon edges.
  // *****************************************************************

  ///@{

  /*!
   *  Split edge.
   *  - Always split, whether or not splitting maximizes
   *    the min triangulation angle.
   *  - Not actually a good idea. Included for testing purposes only.
   *  - Consider all possible triangulations.
   */
  template <typename MESH_TYPE, 
            typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4,
            typename ITYPE, typename MTYPE>
  void split_edge
  (MESH_TYPE & mesh, std::vector<CTYPE1> & vcoord,
   const std::vector<CTYPE2> & split_vcoord, 
   const std::vector<CTYPE3> & interior_coord0,
   const std::vector<CTYPE4> & interior_coord1, 
   const ITYPE ihalf_edge,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const MTYPE max_small_magnitude)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "*** Splitting long edge: ", ihalf_edge, "\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Old polygon: ", 
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge), "\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Old polygon: ", 
       mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge), "\n");
    */

    const DIMENSION_TYPE dimension = mesh.Dimension();
    VERTEX_INDEX_TYPE iv_split;
    POLYGON_INDEX_TYPE inew_poly0, inew_poly1;
    mesh.MESH_TYPE::MESH2D_SPLIT_II_BASE_TYPE::AddVertexSplitPolygonEdgeII
      (dimension, split_vcoord, vcoord, ihalf_edge,
       iv_split, inew_poly0, inew_poly1);
    mesh.SetVertexType(iv_split, MESH_TYPE::VERTEX_TYPE::TRIV_SPLITE);

    mesh.SetNumSplitEdgePolygonVertices(inew_poly0);
    mesh.ComputeCosMaxMinAngle
      (vcoord, inew_poly0, max_small_magnitude, triangulation0_settings);

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  New poly0: ", inew_poly0, "\n");
    mesh.PolygonTriangulationInfo(inew_poly0).Print
      (cerr, "    ", mesh.NumPolygonVertices(inew_poly0));
    */

    if (inew_poly0 != inew_poly1) {
      mesh.SetNumSplitEdgePolygonVertices(inew_poly1);
      mesh.ComputeCosMaxMinAngle
        (vcoord, inew_poly1, max_small_magnitude, triangulation1_settings);

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintPolygonIndexAndVertices
        (cerr, "  New poly1: ", inew_poly1, "\n");
      mesh.PolygonTriangulationInfo(inew_poly1).Print
        (cerr, "    ", mesh.NumPolygonVertices(inew_poly1));
      */

    }

  }


  /*!
   *  Split long polygon edges.
   *  - Return true if some edge is split.
   *  - Always split, whether or not splitting maximizes
   *    the min triangulation angle.
   *  - Not actually a good idea. Included for testing purposes only.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param R Split edge if ratio of edge length to shortest polygon edge
   *     in each polygon is at least R.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE, typename RTYPE>
  bool split_long_polygon_edges
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ipoly, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const RTYPE R)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const int THREE(3);
    const DIMENSION_TYPE dimension = mesh.Dimension();
    const RTYPE Rsquared = R*R;
    HALF_EDGE_INDEX_TYPE ihalf_edge[THREE];
    CTYPE length_squared[THREE];
    POLYGON_INDEX_TYPE ishortest, ishortestX;
    CTYPE shortest_length_squared, shortest_length_squaredX;

    // Stores longest edges in decreasing order by length.
    compute_three_longest_polygon_half_edges
      (mesh.Dimension(), vcoord, mesh, ipoly, ihalf_edge, length_squared);

    if (!mesh.ComputeShortestNonSplitHalfEdge
        (vcoord, ipoly, ishortest, shortest_length_squared))
      { return false; }

    const CTYPE L = shortest_length_squared*Rsquared;

    std::vector<CTYPE> centroid_coord(dimension);
    std::vector<CTYPE> centroidX_coord(dimension);
    std::vector<CTYPE> midpoint_coord(dimension);

    mesh.ComputePolygonCentroid(vcoord, ipoly, centroid_coord);


    for (int j = 0; j < THREE; j++) {
      if (length_squared[j] >= L) {

        // *** DEBUG ***
        /*
        using namespace std;
        mesh.PrintHalfEdgeIndexAndEndpoints
          (cerr, "Attempting to split half edge ", ihalf_edge[j], "\n");
        */

        mesh.ComputeEdgeMidpoint(vcoord, ihalf_edge[j], midpoint_coord);
        
        if (mesh.IsBoundaryEdge(ihalf_edge[j])) {
          split_edge(mesh, vcoord, midpoint_coord, 
                     centroid_coord, centroid_coord, ihalf_edge[j],
                     triangulation0_settings, triangulation1_settings,
                     max_small_magnitude);

          return true;
        }
        else {

          // *** DEBUG ***
          /*
          using namespace std;
          cerr << "  Half edge not boundary." << endl;
          */

          const HALF_EDGE_INDEX_TYPE ihalf_edgeX =
            mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge[j]);
          const POLYGON_INDEX_TYPE ipolyX =
            mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeX);
          if (!mesh.ComputeShortestNonSplitHalfEdge
              (vcoord, ipolyX, ishortestX, shortest_length_squaredX))
            { return false; }
          const CTYPE LX = shortest_length_squaredX*Rsquared;

          // *** DEBUG ***
          /*
          using namespace std;
          cerr << "  length_squared[j]: " << length_squared[j] << endl;
          cerr << "  shortest_length_squaredX: " 
               << shortest_length_squaredX << endl;
          cerr << "  LX: " << LX << endl;
          */

          if (length_squared[j] >= LX) {

            mesh.ComputePolygonCentroid(vcoord, ipolyX, centroidX_coord);
            
            split_edge(mesh, vcoord, midpoint_coord, 
                       centroid_coord, centroidX_coord, ihalf_edge[j],
                       triangulation0_settings, triangulation1_settings,
                       max_small_magnitude);

            return true;
          }
        }
      }
      else {
        // All remaining edges are shorter than sqrt(L).
        return false;
      }
    }

    // No edge split.
    return false;
  }


  /*!
   *  Split long polygon edges.
   *  - Always split, whether or not splitting maximizes
   *    the min triangulation angle.
   *  - Not actually a good idea. Included for testing purposes only.
   */
  template <typename MESH_TYPE, typename CTYPE, typename MTYPE,
            typename RTYPE>
  void split_long_polygon_edges
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings,
   const RTYPE R)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (NUMBER_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      split_long_polygon_edges
        (mesh, vcoord, ipoly, max_small_magnitude, 
         triangulation0_settings, triangulation1_settings, R);
    }

  }

  ///@}


  // *****************************************************************
  /// @name Split polygon edges to reduce edge ratio.
  // *****************************************************************

  ///@{

  /*!
   *  Split longest polygon edge to improve ratio of shortest to longest edge.
   *  - Return true if some edge is split.
   *  - Always split, whether or not splitting maximizes
   *    the min triangulation angle.
   *  - Don't split triangle edges.
   *  - Allow triangulation with interior triangulation vertex
   *    located at polygon centroid.
   *  - Allow triangulations that split off one ear.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on ipoly to determine triangulation.
   *  @param triangulation0_settings Indicates types of triangulations 
   *    to consider for polygon ipoly0.
   *  @param triangulation1_settings Indicates types of triangulations 
   *    to consider for polygon ipoly1.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename ITYPE, typename MTYPE>
  bool split_longest_polygon_edge_to_improve_edge_ratio
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ipoly, const MTYPE max_small_magnitude, 
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const int TWO(3);
    const DIMENSION_TYPE dimension = mesh.Dimension();
    HALF_EDGE_INDEX_TYPE ilongest_half_edge;
    HALF_EDGE_INDEX_TYPE isecond_longest_half_edge;
    CTYPE longest_length_squared, second_longest_length_squared;
    POLYGON_INDEX_TYPE ishortest, ishortestX;
    CTYPE shortest_length_squared, shortest_length_squaredX;

    compute_two_longest_polygon_half_edges
      (mesh.Dimension(), vcoord, mesh, ipoly, 
       ilongest_half_edge, isecond_longest_half_edge,
       longest_length_squared, second_longest_length_squared);

    if (!mesh.ComputeShortestNonSplitHalfEdge
        (vcoord, ipoly, ishortest, shortest_length_squared))
      { return false; }

    std::vector<CTYPE> centroid_coord(dimension);
    std::vector<CTYPE> centroidX_coord(dimension);
    std::vector<CTYPE> midpoint_coord(dimension);

    mesh.ComputePolygonCentroid(vcoord, ipoly, centroid_coord);

    const CTYPE edge_ratio_squared = 
      shortest_length_squared/longest_length_squared;

    // Ratio of (longest_length_squared/4) to second longest length squared.
    const CTYPE new_edge_ratio_squared =
      (longest_length_squared/4.0)/second_longest_length_squared;

    if (edge_ratio_squared < new_edge_ratio_squared) {

      // *** DEBUG ***
      /*
        using namespace std;
        mesh.PrintHalfEdgeIndexAndEndpoints
        (cerr, "Attempting to split half edge ", ilongest_half_edge, "\n");
      */

      mesh.ComputeEdgeMidpoint(vcoord, ilongest_half_edge, midpoint_coord);
        
      if (mesh.IsBoundaryEdge(ilongest_half_edge)) {
        split_edge(mesh, vcoord, midpoint_coord, 
                   centroid_coord, centroid_coord, ilongest_half_edge,
                   triangulation0_settings, triangulation1_settings,
                   max_small_magnitude);

        return true;
      }
      else {
        const HALF_EDGE_INDEX_TYPE ihalf_edgeX =
          mesh.IndexOfNextHalfEdgeAroundEdge(ilongest_half_edge);
        const POLYGON_INDEX_TYPE ipolyX =
          mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeX);

        HALF_EDGE_INDEX_TYPE ilongest_half_edgeX;
        HALF_EDGE_INDEX_TYPE isecond_longest_half_edgeX;
        CTYPE longest_length_squaredX, second_longest_length_squaredX;

        if (mesh.IsTriangle(ipolyX)) { return false; }

        compute_two_longest_polygon_half_edges
          (mesh.Dimension(), vcoord, mesh, ipolyX, 
           ilongest_half_edgeX, isecond_longest_half_edgeX,
           longest_length_squaredX, second_longest_length_squaredX);

        if (!mesh.ComputeShortestNonSplitHalfEdge
            (vcoord, ipolyX, ishortestX, shortest_length_squaredX))
          { return false; }

        const CTYPE edge_ratio_squaredX = 
          shortest_length_squaredX/longest_length_squaredX;

        // Ratio of (longest_length_squaredX/4) 
        //   to second longest length squared in polyX.
        const CTYPE new_edge_ratio_squaredX =
          (longest_length_squaredX/4.0)/second_longest_length_squaredX;

        if ((ihalf_edgeX == ilongest_half_edgeX) &&
            (edge_ratio_squaredX < new_edge_ratio_squaredX)) {

          mesh.ComputePolygonCentroid(vcoord, ipolyX, centroidX_coord);

          // *** DEBUG ***
          /*
          using namespace std;
          cerr << "In " << __func__ << endl;
          mesh.PrintHalfEdgeIndexAndEndpoints
            (cerr, "  Splitting half edge ", ilongest_half_edge, "\n");
          */
            
          split_edge(mesh, vcoord, midpoint_coord, 
                     centroid_coord, centroidX_coord, 
                     ilongest_half_edge,
                     triangulation0_settings, triangulation1_settings,
                     max_small_magnitude);

          return true;
        }
      }
    }

    // No edge split.
    return false;
  }

  
  /*!
   *  Split longest polygon edge to improve ratio of shortest to longest edge.
   *  - Always split, whether or not splitting maximizes
   *    the min triangulation angle.
   *  - Don't split triangle edges.
   */
  template <typename MESH_TYPE, typename CTYPE, typename MTYPE>
  void split_longest_polygon_edge_to_improve_edge_ratio
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation0_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation1_settings)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const NUMBER_TYPE num_poly = mesh.NumPolygons();

    for (NUMBER_TYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (mesh.IsPolygonDeleted(ipoly)) { continue; }
      if (mesh.IsTriangle(ipoly)) { continue; }

      split_longest_polygon_edge_to_improve_edge_ratio
        (mesh, vcoord, ipoly, max_small_magnitude, 
         triangulation0_settings, triangulation1_settings);
    }
  }

  ///@}


  // *****************************************************************
  /// @name Split three edges incident on one vertex to max min triangulation angle.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Split three polygon edges incident on one vertex to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param triangulation_settings Indicates which types of triangulations to consider.
   *  @param flag_max_min_all If true, split edges if the split improves
   *      the min triangulation of the two inner polygons
   *      and their two neighbors.
   *    - Otherwise, split edges only if the split improves
   *      the min triangulation of the two inner polygons.
   */
  template <typename MESH_TYPE, typename CTYPE, typename ITYPE,
            typename MTYPE>
  bool split_three_edges_incident_on_toV_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edgeB1, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::COS_TYPE COS_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const NUMBER_TYPE FOUR_POLYGONS(4);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB2 =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC1 =
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edgeB1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeA2 =
      mesh.IndexOfPrevHalfEdgeInPolygon(ihalf_edgeB2);

    POLYGON_INDEX_TYPE ipoly[FOUR_POLYGONS];
    POLYGON_TRIANGULATION_RESULT_E16<COS_TYPE, ITYPE>
      poly_tri_result[FOUR_POLYGONS];
    COS_TYPE cos_min_angle;
    bool flag_zero;

    if (mesh.AreTwoHalfEdgesAroundFromVertex(ihalf_edgeB2)) {
      // Two polygons share two edges. Do not split.
      return false;
    }

    if (mesh.AreThreeHalfEdgesAroundFromVertex(ihalf_edgeB2)) {
      // *** SHOULD BE ABLE TO HANDLE THIS CASE, BUT... ***
      // FromVertex(ihalf_edgeB2) is incident on exactly three polygons.
      return(false);
    }

    ipoly[0] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeC1);
    ipoly[1] = mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB1);
    ipoly[2] = mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB2);
    ipoly[3] = mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeA2);

    std::vector<CTYPE> midpointA_coord(dimension);
    std::vector<CTYPE> midpointB_coord(dimension);
    std::vector<CTYPE> midpointC_coord(dimension);
    std::vector<CTYPE> centroid0_coord(dimension);
    std::vector<CTYPE> centroid1_coord(dimension);
    std::vector<CTYPE> centroid2_coord(dimension);
    std::vector<CTYPE> centroid3_coord(dimension);

    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeA2, midpointA_coord);
    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeB1, midpointB_coord);
    mesh.ComputeEdgeMidpoint(vcoord, ihalf_edgeC1, midpointC_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[0], centroid0_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[1], centroid1_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[2], centroid2_coord);
    mesh.ComputePolygonCentroid(vcoord, ipoly[3], centroid3_coord);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints(cerr, "  Half edge: ", ihalf_edgeB1, "\n");
    */

    compute_split_three_incident_edges_triangulation_to_max_min_angle
      (dimension, vcoord, midpointA_coord, midpointB_coord, midpointC_coord,
       centroid0_coord, centroid1_coord, centroid2_coord, centroid3_coord,
       mesh, ihalf_edgeB1, max_small_magnitude,
       inner_poly_triangulation_settings, outer_poly_triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Polygon 1: ", mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB1), "\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Polygon 2: ", mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA2), "\n");
    cerr << "  cos_min_angle: " << cos_min_angle << endl;
    cerr << "  poly_tri_result[1].IsEar(0): " << int(poly_tri_result[1].IsEar(0)) << endl;
    cerr << "  poly_tri_result[1].IsEar(1): " << int(poly_tri_result[1].IsEar(1)) << endl;
    cerr << "  poly_tri_result[1].IsEar(2): " << int(poly_tri_result[1].IsEar(2)) << endl;
    cerr << "  poly_tri_result[1].IsEar(3): " << int(poly_tri_result[1].IsEar(3)) << endl;
    cerr << "  poly_tri_result[2].IsEar(0): " << int(poly_tri_result[2].IsEar(0)) << endl;
    cerr << "  poly_tri_result[2].IsEar(1): " << int(poly_tri_result[2].IsEar(1)) << endl;
    cerr << "  poly_tri_result[2].IsEar(2): " << int(poly_tri_result[2].IsEar(2)) << endl;
    cerr << "  poly_tri_result[2].IsEar(3): " << int(poly_tri_result[2].IsEar(3)) << endl;
    cerr << "  poly_tri_result[3]:" << endl;
    poly_tri_result[3].Print(cerr, "    ");
    */


    if (flag_zero) { return(false); }

    const COS_TYPE cos_min_angle_poly4 =
      mesh.PolygonIVCosMinTriangulationAngle
      (ipoly[0], ipoly[1], ipoly[2], ipoly[3]);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  cos_min_angle_poly4: " << cos_min_angle_poly4 << endl;
    */

    if (cos_min_angle < mesh.PolygonCosMinTriangulationAngle(ipoly[1]) ||
        cos_min_angle < mesh.PolygonCosMinTriangulationAngle(ipoly[2]) ||
        (flag_max_min_all &&
         cos_min_angle < cos_min_angle_poly4)) {

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "In " << __func__ << endl;
      for (int i = 0; i < 4; i++) {
        mesh.PrintPolygonIndexAndVertices
          (cerr, "  *** Splitting polygon ", ipoly[i], "\n");
      }
      */

      add_vertices_split_three_incident_edges_allow_IV
        (mesh, vcoord,
         midpointA_coord, midpointB_coord, midpointC_coord,
         centroid0_coord, centroid1_coord,
         centroid2_coord, centroid3_coord,
         ihalf_edgeB1, poly_tri_result);

      return(true);
    }

    return(false);
  }


  /*!
   *  @brief Split three polygon edges incident on one vertex to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split longest polygon edge
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation_settings Indicates which types of triangulations to consider.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *     - Otherwise, split if adjacent polygon triangulation angles
   *     are below threshold.
   *  @param flag_max_min_all If true, split edges if the split improves
   *      the min triangulation of the two inner polygons 
   *      and their two neighbors.
   *    - Otherwise, split edges only if the split improves
   *      the min triangulation of the two inner polygons.
   */
  template <typename MESH_TYPE, typename CTYPE, typename ITYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  bool split_three_edges_incident_on_toV_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const ITYPE ihalf_edgeB1,
   const COS_TYPE cos_angle_threshold, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const RTYPE Rsquared,
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const DIMENSION_TYPE dimension = mesh.Dimension();
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB2 = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB1);
    const POLYGON_INDEX_TYPE ipoly1 =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB1);
    const POLYGON_INDEX_TYPE ipoly2 =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB2);
    CTYPE shortest_length1_squared, shortest_length2_squared;

    if (mesh.IsBoundaryEdge(ihalf_edgeB1)) { return(false); }

    mesh.ComputeShortestHalfEdge(vcoord, ipoly1, shortest_length1_squared);
    mesh.ComputeShortestHalfEdge(vcoord, ipoly2, shortest_length2_squared);
    const CTYPE L1 = shortest_length1_squared*Rsquared;
    const CTYPE L2 = shortest_length2_squared*Rsquared;
    const CTYPE lengthB_squared = 
      compute_edge_length_squared(dimension, vcoord, mesh, ihalf_edgeB1);

    if (lengthB_squared < L1 || lengthB_squared < L2) {
      // Edge is not much longer than shortest edge.  Don't try to split.
      return(false);
    }

    const HALF_EDGE_INDEX_TYPE ihalf_edgeC1 =
                          mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edgeB1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeA2 =
                          mesh.IndexOfPrevHalfEdgeInPolygon(ihalf_edgeB2);
    const CTYPE lengthC_squared = 
      compute_edge_length_squared(dimension, vcoord, mesh, ihalf_edgeC1);
    const CTYPE lengthA_squared = 
      compute_edge_length_squared(dimension, vcoord, mesh, ihalf_edgeA2);

    if (lengthC_squared < L1 || lengthA_squared < L2) {
      // Incident edges are not much longer than shortest edge.  Don't try to split.
      return(false);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints(cerr, "  Half edge: ", ihalf_edgeB1, "\n");
    mesh.PrintPolygonIndexAndVertices(cerr, "  ipoly1: ", ipoly1, "\n");
    mesh.PrintPolygonIndexAndVertices(cerr, "  ipoly2: ", ipoly2, "\n");
    cerr << "  cos_angle_threshold: " << cos_angle_threshold << endl;
    cerr << "  PolygonCosMinTriangulationAngle(" << ipoly1 << ")" 
         << mesh.PolygonCosMinTriangulationAngle(ipoly1) << endl;
    cerr << "  PolygonCosMinTriangulationAngle(" << ipoly2 << ")" 
         << mesh.PolygonCosMinTriangulationAngle(ipoly2) << endl;
    */

    if (mesh.PolygonCosMinTriangulationAngle(ipoly1) > cos_angle_threshold ||
        mesh.PolygonCosMinTriangulationAngle(ipoly2) > cos_angle_threshold) {

      const bool flag = split_three_edges_incident_on_toV_to_max_min_triangulation_angle
        (mesh, vcoord, ihalf_edgeB1, max_small_magnitude,
         inner_poly_triangulation_settings, outer_poly_triangulation_settings, 
         flag_max_min_all);
      return(flag);
    }
    else if (!flag_poly_angle_below_threshold) {
      // Attempt to split three incident edges if adjacent polygon angles 
      //   are below threshold.

      const POLYGON_INDEX_TYPE ipoly0 = 
        mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeC1);
      const POLYGON_INDEX_TYPE ipoly3 = 
        mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeA2);

      if ((mesh.PolygonCosMinTriangulationAngle(ipoly0) > cos_angle_threshold) ||
          (mesh.PolygonCosMinTriangulationAngle(ipoly3) > cos_angle_threshold)) {

        const bool flag = split_three_edges_incident_on_toV_to_max_min_triangulation_angle
          (mesh, vcoord, ihalf_edgeB1, max_small_magnitude,
           inner_poly_triangulation_settings, outer_poly_triangulation_settings,
           flag_max_min_all);
        return(flag);
      }
    }


    return(false);
  }

  /*!
   *  @brief Split three polygon edges incident on one vertex to maximize minimum triangulation angle.
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split longest polygon edge
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation_settings Indicates which types of triangulations to consider.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *     - Otherwise, split if adjacent polygon triangulation angles
   *     are below threshold.
   *  @param flag_max_min_all If true, split edges if the split improves
   *      the min triangulation of the two inner polygons 
   *      and their two neighbors.
   *    - Otherwise, split edges only if the split improves
   *      the min triangulation of the two inner polygons.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_three_incident_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const COS_TYPE cos_angle_threshold, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const RTYPE Rsquared,
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const NUMBER_TYPE num_half_edges = mesh.NumHalfEdges();

    for (HALF_EDGE_INDEX_TYPE ihalf_edge0 = 0; ihalf_edge0 < num_half_edges; 
         ihalf_edge0++) {
 
      const POLYGON_INDEX_TYPE ipoly0 =
        mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0);

      if (mesh.IsPolygonDeleted(ipoly0)) { continue; }

      split_three_edges_incident_on_toV_to_max_min_triangulation_angle
        (mesh, vcoord, ihalf_edge0, cos_angle_threshold, max_small_magnitude,
         inner_poly_triangulation_settings, outer_poly_triangulation_settings, 
         Rsquared, 
         flag_poly_angle_below_threshold, flag_max_min_all); 
    }

  }


  /*!
   *  @brief Split three polygon edges incident on one vertex to maximize minimum triangulation angle.
   *  - Version with array cos_angle_threshold[].
   *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
   *     has been called on all polygons to determine triangulation.
   *  @param cos_angle_threshold Attempt to split longest polygon edge
   *     if cos_min_triangulation_angle is above cos_angle_threshold.
   *  @param triangulation_settings Indicates which types of triangulations to consider.
   *  @param flag_poly_angle_below_threshold If true, only split edges
   *     if focus polygon triangulation angle is below threshold.
   *     - Otherwise, split if adjacent polygon triangulation angles
   *     are below threshold.
   *  @param flag_max_min_all If true, split edges if the split improves
   *      the min triangulation of the two inner polygons 
   *      and their two neighbors.
   *    - Otherwise, split edges only if the split improves
   *      the min triangulation of the two inner polygons.
   */
  template <typename MESH_TYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE, typename RTYPE>
  void split_three_incident_edges_to_max_min_triangulation_angle
  (MESH_TYPE & mesh, std::vector<CTYPE> & vcoord,
   const std::vector<COS_TYPE> & cos_angle_threshold, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   const RTYPE R,
   const bool flag_poly_angle_below_threshold,
   const bool flag_max_min_all)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    const RTYPE Rsquared = R*R;

    for (NUMBER_TYPE i = 0; i < cos_angle_threshold.size(); i++) {
      split_three_incident_edges_to_max_min_triangulation_angle
        (mesh, vcoord, cos_angle_threshold[i], max_small_magnitude,
         inner_poly_triangulation_settings, outer_poly_triangulation_settings, 
         Rsquared, flag_poly_angle_below_threshold, flag_max_min_all); 
    }
  }

  ///@}

}


#endif

