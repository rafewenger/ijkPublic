/*!
 *  @file ijkmesh2D_angle.tpp
 *  @brief ijk template functions for computing triangulations and triangulation angles.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2019-2022 Rephael Wenger

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

#ifndef _IJKMESH2D_ANGLE_
#define _IJKMESH2D_ANGLE_

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkmesh2D_geom.tpp"
#include "ijkmesh2D_split.tpp"
#include "ijktri2D_info.tpp"
#include "ijktri2D_angle.tpp"

#include <algorithm>
#include <numeric>
#include <vector>

// *** DEBUG ***
#include "ijkprint.tpp"

namespace IJK {

  // *****************************************************************
  //! @name Set triangulation encoding
  // *****************************************************************

  ///@{

  /*!
   *  @brief Set fan triangulation encoding in split edge polygon.
   *  @param ihalf_edgeA Split half edge.
   *  @param tri_vertex_index Index (0..num_polyv-1) of triangulation vertex.
   */
  template <typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE1,
            typename ENCODING_TYPE>
  void set_split_edge_fan_triangulation_encoding
  (const MESH_TYPE & mesh, const ITYPE0 ihalf_edgeA,
   const ITYPE1 tri_vertex_index,
   ENCODING_TYPE & tri_encoding)
  {
    typedef typename MESH_TYPE::INDEX_TYPE INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const INDEX_TYPE iloc_new = 
      mesh.LocationInSplitEdgePolygon(ihalf_edgeA, tri_vertex_index);

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NUMBER_TYPE numv = mesh.NumPolygonVertices(ipoly);

    // Set triangulation_encoding of triangulation in polygon
    //   with split half edge.
    // Note: Polygon with split half edge has (numv+1) vertices.
    tri_encoding.SetFan(iloc_new, numv+1);
  }

  ///@}


  // *****************************************************************
  //! @name Create list of pointers to polygon vertices.
  // *****************************************************************

  /*!
     @brief Create list of pointers to vertices of polygon.
     @param vcoord[] List of mesh vertices.
     @param ipoly Polygon index.
  */
  template <typename DTYPE, typename CTYPEV,
            typename MESH_TYPE, typename ITYPEP, typename CTYPE_POLYV>
  void create_poly_vertex_list
  (const DTYPE dimension, const CTYPEV vcoord[], 
   const MESH_TYPE & mesh, const ITYPEP ipoly,
   std::vector<CTYPE_POLYV> & poly_vcoord_ptr)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);

    poly_vcoord_ptr.resize(num_poly_vert);

    for (NUMBER_TYPE i = 0; i < num_poly_vert; i++) {
      const VERTEX_INDEX_TYPE iv = mesh.PolygonVertex(ipoly, i);
      poly_vcoord_ptr[i] = (vcoord + iv*dimension);
    }

  }


  /*!
   *  @brief Create list of pointers to vertices of polygon with split edge.
   *  @param vcoord[] List of mesh vertices.
   *  @param split_coord[] Location of vertex splitting edge.
   *  @param ihalf_edgeA Index of split edge.
   */
  template <typename DTYPE, typename CTYPEV,
            typename MESH_TYPE, typename ITYPEH, typename CTYPE_POLYV>
  void create_split_edge_poly_vertex_list
  (const DTYPE dimension,
   const CTYPEV vcoord[], const CTYPEV split_coord[],
   const MESH_TYPE & mesh, const ITYPEH ihalf_edgeA,
   std::vector<CTYPE_POLYV> & poly_vcoord_ptr)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const int NUM_SPLIT_EDGES(1);
    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
    const NUMBER_TYPE num_split_poly_vert = num_poly_vert+NUM_SPLIT_EDGES;
    IJK::PROCEDURE_ERROR error("create_split_edge_poly_vertex_list");

    poly_vcoord_ptr.resize(num_split_poly_vert);

    NUMBER_TYPE i1 = 0;
    for (NUMBER_TYPE i0 = 0; i0 < num_poly_vert; i0++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge = 
        mesh.HalfEdgeIndex(ipoly, i0);
      const VERTEX_INDEX_TYPE iv = mesh.FromVertexIndex(ihalf_edge);
      poly_vcoord_ptr[i1] = (vcoord + iv*dimension);
      i1++;

      if (ihalf_edge == ihalf_edgeA) {
        poly_vcoord_ptr[i1] = split_coord;
        i1++;
      }
    }

    if (i1 != num_split_poly_vert) {
      error.AddMessage("Programming error.");
      error.AddMessage("  Incorrect number of pointers to vertex coordinates added to poly_vcoord_ptr[].");
      error.AddMessage("  Expected ", num_split_poly_vert, " pointers. ",
                       " Added ", i1, " pointers.");
      throw error;
    }

  }
   

  /*!
   *  @brief Create list of pointers to vertices of polygon with split edge.
   *  - For each vertex i, set flag_not_ear[i] to true if vertex i
   *    is an internal vertex and is incident on exactly two polygons.
   *  @param vcoord[] List of mesh vertices.
   *  @param split_coord[] Location of vertex splitting edge.
   *  @param ihalf_edgeA Index of split edge.
   *  @param[out] flag_not_ear[i] True if i'th vertex is an interior
   *    vertex incident on exactly two polygons.
   */
  template <typename DTYPE, typename CTYPEV,
            typename MESH_TYPE, typename ITYPEH,
            typename CTYPE_POLYV, typename BIT_SET_TYPE>
  void create_split_edge_poly_vertex_list_and_flag_not_ear
  (const DTYPE dimension,
   const CTYPEV vcoord[], const CTYPEV split_coord[],
   const MESH_TYPE & mesh, const ITYPEH ihalf_edgeA,
   std::vector<CTYPE_POLYV> & poly_vcoord_ptr,
   BIT_SET_TYPE & flag_not_ear)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const int NUM_SPLIT_EDGES(1);
    const POLYGON_INDEX_TYPE ipoly =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
    const NUMBER_TYPE num_split_poly_vert = num_poly_vert+NUM_SPLIT_EDGES;
    IJK::PROCEDURE_ERROR error("create_split_edge_poly_vertex_list");

    poly_vcoord_ptr.resize(num_split_poly_vert);
    flag_not_ear.reset();

    NUMBER_TYPE i1 = 0;
    for (NUMBER_TYPE i0 = 0; i0 < num_poly_vert; i0++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge =
        mesh.HalfEdgeIndex(ipoly, i0);
      const VERTEX_INDEX_TYPE iv = mesh.FromVertexIndex(ihalf_edge);
      poly_vcoord_ptr[i1] = (vcoord + iv*dimension);

      if (mesh.AreTwoHalfEdgesAroundFromVertex(ihalf_edge) &&
          !mesh.IsBoundaryVertex(iv))
        { flag_not_ear[i1] = true; }

      i1++;

      if (ihalf_edge == ihalf_edgeA) {
        poly_vcoord_ptr[i1] = split_coord;
        flag_not_ear[i1] = true;
        i1++;
      }
    }

    if (i1 != num_split_poly_vert) {
      error.AddMessage("Programming error.");
      error.AddMessage("  Incorrect number of pointers to vertex coordinates added to poly_vcoord_ptr[].");
      error.AddMessage("  Expected ", num_split_poly_vert, " pointers. ",
                       " Added ", i1, " pointers.");
      throw error;
    }

  }


  /*!
   *  @brief Create list of pointers to vertices of polygon with two split edges.
   *  @param vcoord[] List of mesh vertices.
   *  @param split_coordA[] Location of vertex splitting ihalf_edgeA.
   *  @param split_coordB[] Location of vertex splitting ihalf_edgeB.
   *  @param ihalf_edgeA Index of first split edge.
   *  @param ihalf_edgeB Index of second split edge.
   */
  template <typename DTYPE, typename CTYPEV,
            typename CTYPEA, typename CTYPEB,
            typename MESH_TYPE, typename ITYPEH, typename CTYPE_POLYV>
  void create_two_split_edges_poly_vertex_list
  (const DTYPE dimension, const CTYPEV vcoord[],
   const CTYPEA splitA_coord[], const CTYPEB splitB_coord[],
   const MESH_TYPE & mesh,
   const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
   std::vector<CTYPE_POLYV> & poly_vcoord_ptr)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const int NUM_SPLIT_EDGES(2);
    const POLYGON_INDEX_TYPE ipoly =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
    const NUMBER_TYPE num_split_poly_vert = num_poly_vert+NUM_SPLIT_EDGES;
    IJK::PROCEDURE_ERROR error("create_two_split_edges_poly_vertex_list");

    poly_vcoord_ptr.resize(num_split_poly_vert);

    NUMBER_TYPE i1 = 0;
    for (NUMBER_TYPE i0 = 0; i0 < num_poly_vert; i0++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge =
        mesh.HalfEdgeIndex(ipoly, i0);
      const VERTEX_INDEX_TYPE iv = mesh.FromVertexIndex(ihalf_edge);
      poly_vcoord_ptr[i1] = (vcoord + iv*dimension);
      i1++;

      if (ihalf_edge == ihalf_edgeA) {
        poly_vcoord_ptr[i1] = splitA_coord;
        i1++;
      }
      else if (ihalf_edge == ihalf_edgeB) {
        poly_vcoord_ptr[i1] = splitB_coord;
        i1++;
      }
    }

    if (i1 != num_split_poly_vert) {
      error.AddMessage("Programming error.");
      error.AddMessage("  Incorrect number of pointers to vertex coordinates added to poly_vcoord_ptr[].");
      error.AddMessage("  Expected ", num_split_poly_vert, " pointers. ",
                       " Added ", i1, " pointers.");
      throw error;
    }

  }


  /*!
    *  @brief Create list of pointers to vertices of polygon with two split edges.
    *  - For each vertex i, set flag_not_ear[i] to true if vertex i
    *    is an internal vertex and is incident on exactly two polygons.
    *  @param vcoord[] List of mesh vertices.
    *  @param splitA_coord[] Location of vertex splitting ihalf_edgeA.
    *  @param splitB_coord[] Location of vertex splitting ihalf_edgeB.
    *  @param ihalf_edgeA Index of first split edge.
    *  @param ihalf_edgeB Index of second split edge.
    *  @param[out] flag_not_ear[i] True if i'th vertex is an interior
    *    vertex incident on exactly two polygons.
    */
   template <typename DTYPE, typename CTYPEV,
             typename CTYPEA, typename CTYPEB,
             typename MESH_TYPE, typename ITYPEH,
             typename CTYPE_POLYV, typename BIT_SET_TYPE>
   void create_two_split_edges_poly_vertex_list_and_flag_not_ear
   (const DTYPE dimension, const CTYPEV vcoord[],
    const CTYPEA splitA_coord[], const CTYPEB splitB_coord[],
    const MESH_TYPE & mesh,
    const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
    std::vector<CTYPE_POLYV> & poly_vcoord_ptr,
    BIT_SET_TYPE & flag_not_ear)
   {
     typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
     typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
     typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
     typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

     const int NUM_SPLIT_EDGES(2);
     const POLYGON_INDEX_TYPE ipoly =
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
     const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
     const NUMBER_TYPE num_split_poly_vert = num_poly_vert+NUM_SPLIT_EDGES;
     IJK::PROCEDURE_ERROR error
       ("create_two_split_edges_poly_vertex_list_flag_not_ear");

     poly_vcoord_ptr.resize(num_split_poly_vert);
     flag_not_ear.reset();

     NUMBER_TYPE i1 = 0;
     for (NUMBER_TYPE i0 = 0; i0 < num_poly_vert; i0++) {
       const HALF_EDGE_INDEX_TYPE ihalf_edge =
         mesh.HalfEdgeIndex(ipoly, i0);
       const VERTEX_INDEX_TYPE iv = mesh.FromVertexIndex(ihalf_edge);
       poly_vcoord_ptr[i1] = (vcoord + iv*dimension);

       if (mesh.AreTwoHalfEdgesAroundFromVertex(ihalf_edge) &&
            !mesh.IsBoundaryVertex(iv))
          { flag_not_ear[i1] = true; }

       i1++;

       if (ihalf_edge == ihalf_edgeA) {
         poly_vcoord_ptr[i1] = splitA_coord;
         flag_not_ear[i1] = true;
         i1++;
       }
       else if (ihalf_edge == ihalf_edgeB) {
         poly_vcoord_ptr[i1] = splitB_coord;
         flag_not_ear[i1] = true;
         i1++;
       }
     }

     if (i1 != num_split_poly_vert) {
       error.AddMessage("Programming error.");
       error.AddMessage("  Incorrect number of pointers to vertex coordinates added to poly_vcoord_ptr[].");
       error.AddMessage("  Expected ", num_split_poly_vert, " pointers. ",
                        " Added ", i1, " pointers.");
       throw error;
     }

   }


   /*!
    *  @brief Create list of pointers to vertices of polygon with three split edges.
    *  @param vcoord[] List of mesh vertices.
    *  @param splitA_coord[] Location of vertex splitting ihalf_edgeA.
    *  @param splitB_coord[] Location of vertex splitting ihalf_edgeB.
    *  @param splitC_coord[] Location of vertex splitting ihalf_edgeC.
    *  @param ihalf_edgeA Index of first split edge.
    *  @param ihalf_edgeB Index of second split edge.
    *  @param ihalf_edgeC Index of third split edge.
    */
   template <typename DTYPE, typename CTYPEV,
             typename CTYPEA, typename CTYPEB, typename CTYPEC,
             typename MESH_TYPE, typename ITYPEH, typename CTYPE_POLYV>
   void create_three_split_edges_poly_vertex_list
   (const DTYPE dimension, const CTYPEV vcoord[],
    const CTYPEA splitA_coord[], const CTYPEB splitB_coord[],
    const CTYPEC splitC_coord[],
    const MESH_TYPE & mesh,
    const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
    const ITYPEH ihalf_edgeC,
    std::vector<CTYPE_POLYV> & poly_vcoord_ptr)
   {
     typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
     typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
     typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
     typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

     const int NUM_SPLIT_EDGES(3);
     const POLYGON_INDEX_TYPE ipoly =
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
     const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
     const NUMBER_TYPE num_split_poly_vert = num_poly_vert+NUM_SPLIT_EDGES;
     IJK::PROCEDURE_ERROR error("create_three_split_edges_poly_vertex_list");

     poly_vcoord_ptr.resize(num_split_poly_vert);

     NUMBER_TYPE i1 = 0;
     for (NUMBER_TYPE i0 = 0; i0 < num_poly_vert; i0++) {
       const HALF_EDGE_INDEX_TYPE ihalf_edge =
         mesh.HalfEdgeIndex(ipoly, i0);
       const VERTEX_INDEX_TYPE iv = mesh.FromVertexIndex(ihalf_edge);
       poly_vcoord_ptr[i1] = (vcoord + iv*dimension);
       i1++;

       if (ihalf_edge == ihalf_edgeA) {
         poly_vcoord_ptr[i1] = splitA_coord;
         i1++;
       }
       else if (ihalf_edge == ihalf_edgeB) {
         poly_vcoord_ptr[i1] = splitB_coord;
         i1++;
       }
       else if (ihalf_edge == ihalf_edgeC) {
         poly_vcoord_ptr[i1] = splitC_coord;
         i1++;
       }
     }

     if (i1 != num_split_poly_vert) {
       error.AddMessage("Programming error.");
       error.AddMessage("  Incorrect number of pointers to vertex coordinates added to poly_vcoord_ptr[].");
       error.AddMessage("  Expected ", num_split_poly_vert, " pointers. ",
                        " Added ", i1, " pointers.");
       throw error;
     }

   }


   /*!
     *  @brief Create list of pointers to vertices of polygon with three split edges.
     *  @param vcoord[] List of mesh vertices.
     *  @param splitA_coord[] Location of vertex splitting ihalf_edgeA.
     *  @param splitB_coord[] Location of vertex splitting ihalf_edgeB.
     *  @param splitC_coord[] Location of vertex splitting ihalf_edgeC.
     *  @param ihalf_edgeA Index of first split edge.
     *  @param ihalf_edgeB Index of second split edge.
     *  @param ihalf_edgeC Index of third split edge.
     *  @param[out] flag_not_ear[i] True if i'th vertex is an interior
     *    vertex incident on exactly two polygons.
     */
   template <typename DTYPE, typename CTYPEV,
             typename CTYPEA, typename CTYPEB, typename CTYPEC,
             typename MESH_TYPE, typename ITYPEH,
             typename CTYPE_POLYV, typename BIT_SET_TYPE>
   void create_three_split_edges_poly_vertex_list_and_flag_not_ear
   (const DTYPE dimension, const CTYPEV vcoord[],
    const CTYPEA splitA_coord[], const CTYPEB splitB_coord[],
    const CTYPEC splitC_coord[],
    const MESH_TYPE & mesh,
    const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
    const ITYPEH ihalf_edgeC,
    std::vector<CTYPE_POLYV> & poly_vcoord_ptr,
    BIT_SET_TYPE & flag_not_ear)
   {
     typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
     typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
     typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
     typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

     const int NUM_SPLIT_EDGES(3);
     const POLYGON_INDEX_TYPE ipoly =
         mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
     const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
     const NUMBER_TYPE num_split_poly_vert = num_poly_vert+NUM_SPLIT_EDGES;
     IJK::PROCEDURE_ERROR error
       ("create_three_split_edges_poly_vertex_list_and_flag_not_ear");

     poly_vcoord_ptr.resize(num_split_poly_vert);

     NUMBER_TYPE i1 = 0;
     for (NUMBER_TYPE i0 = 0; i0 < num_poly_vert; i0++) {
       const HALF_EDGE_INDEX_TYPE ihalf_edge =
         mesh.HalfEdgeIndex(ipoly, i0);
       const VERTEX_INDEX_TYPE iv = mesh.FromVertexIndex(ihalf_edge);
       poly_vcoord_ptr[i1] = (vcoord + iv*dimension);

       if (mesh.AreTwoHalfEdgesAroundFromVertex(ihalf_edge) &&
            !mesh.IsBoundaryVertex(iv))
          { flag_not_ear[i1] = true; }

       i1++;

       if (ihalf_edge == ihalf_edgeA) {
         poly_vcoord_ptr[i1] = splitA_coord;
         flag_not_ear[i1] = true;
         i1++;
       }
       else if (ihalf_edge == ihalf_edgeB) {
         poly_vcoord_ptr[i1] = splitB_coord;
         flag_not_ear[i1] = true;
         i1++;
       }
       else if (ihalf_edge == ihalf_edgeC) {
         poly_vcoord_ptr[i1] = splitC_coord;
         flag_not_ear[i1] = true;
         i1++;
       }
     }

     if (i1 != num_split_poly_vert) {
       error.AddMessage("Programming error.");
       error.AddMessage("  Incorrect number of pointers to vertex coordinates added to poly_vcoord_ptr[].");
       error.AddMessage("  Expected ", num_split_poly_vert, " pointers. ",
                        " Added ", i1, " pointers.");
       throw error;
     }

   }


  // *****************************************************************
  /// @name Compute cos of min angle in a given polygon triangulation.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute cosine of the min angle in a polygon fan triangulation.
   *  @param tri_vertex_index Fan triangulation from vertex tri_vertex_index.
   */
  template <typename DTYPE, typename CTYPE, typename MESH_TYPE,
            typename ITYPE0, typename ITYPE1, typename MTYPE,
            typename COS_TYPE>
  void compute_cos_min_polygon_fan_triangulation_angle
  (const DTYPE dimension, const CTYPE vcoord[], const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   const ITYPE1 tri_vertex_index, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const NUMBER_TYPE numv = mesh.NumPolygonVertices(ipoly);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << " poly " << ipoly << endl;
    */

    compute_cos_min_fan_triangulation_angle
      (dimension, vcoord, numv, mesh.PolygonVertexList(ipoly),
       max_small_magnitude, tri_vertex_index, cos_min_angle, flag_zero);
  }


  /// @brief Compute cosine of the min angle in the polygon triangulation from interior vertex.
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE0, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_polygon_interior_vertex_triangulation_angle
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const COORD_TYPE1 interior_vcoord[], const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const NUMBER_TYPE numv = mesh.NumPolygonVertices(ipoly);

    // Initialize
    cos_min_angle = -1;
    flag_zero = false;

    for (ITYPE0 i0 = 0; i0 < numv; i0++) {
      const ITYPE0 i1 = (i0+1)%numv;

      const VERTEX_INDEX_TYPE iv0 = mesh.IndexOfVertexInPolygon(ipoly, i0);
      const VERTEX_INDEX_TYPE iv1 = mesh.IndexOfVertexInPolygon(ipoly, i1);

      const COORD_TYPE0 * vcoord0 = vcoord + iv0*dimension;
      const COORD_TYPE0 * vcoord1 = vcoord + iv1*dimension;

      COS_TYPE cosA;
      bool flagA_zero;
      IJK::compute_cos_min_triangle_angle
        (dimension, interior_vcoord, vcoord0, vcoord1, max_small_magnitude, 
         cosA, flagA_zero);
      if (flagA_zero) { flag_zero = true; }
      if (cosA > cos_min_angle) { cos_min_angle = cosA; }
    }

  }

  ///@}


  // *****************************************************************
  /// @name Compute cos of min angle in triangulation of a modified polygon.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute cosine of the min angle in interior vertex triangulation.
   *  - Triangulation formed by vcoordX[] and every edge of a polygon
   *    except ihalf_edge0.
   *  @param ihalf_edge0 Index of the half edge not forming a triangle.
   *    - Polygon is IndexOfPolygonContainingHalfEdge(ihalf_edge0).
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE0, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_polygon_minus_edge_triangulation_angle
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const COORD_TYPE1 vcoordX[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0);

    const NUMBER_TYPE nume = mesh.NumPolygonEdges(ipoly);

    // Initialize
    cos_min_angle = -1;
    flag_zero = false;

    HALF_EDGE_INDEX_TYPE ihalf_edge = ihalf_edge0;
    
    for (NUMBER_TYPE i = 1; i < nume; i++) {

      ihalf_edge = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge);

      const VERTEX_INDEX_TYPE iv0 = mesh.FromVertexIndex(ihalf_edge);
      const VERTEX_INDEX_TYPE iv1 = mesh.ToVertexIndex(ihalf_edge);

      const COORD_TYPE0 * vcoord0 = vcoord + iv0*dimension;
      const COORD_TYPE0 * vcoord1 = vcoord + iv1*dimension;

      COS_TYPE cosX;
      bool flagX_zero;
      IJK::compute_cos_min_triangle_angle
        (dimension, vcoordX, vcoord0, vcoord1, max_small_magnitude, 
         cosX, flagX_zero);

      if (flagX_zero) { flag_zero = true; }
      if (cosX > cos_min_angle) { cos_min_angle = cosX; }
    }

  }


  /*!
   *  @brief Compute cosine of the min angle in interior vertex triangulation.
   *  - Triangulation formed by vcoordX[] and every edge of a polygon
   *    except ihalf_edge0.
   *  - Version using C++ STL vector vcoord[].
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE0, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_polygon_minus_edge_triangulation_angle
  (const DTYPE dimension, 
   const std::vector<COORD_TYPE0> & vcoord,
   const COORD_TYPE1 vcoordX[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_polygon_minus_edge_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), vcoordX,
       mesh, ihalf_edge0, max_small_magnitude,
       cos_min_angle, flag_zero);
  }


  /*!
    *  @brief Compute cosine of the min angle in fan triangulation.
    *  - Triangulation formed by one vertex and every edge of a polygon
    *    except ihalf_edge0 and NextHalfEdgeInPolygon(ihalf_edge0).
    *  @param ihalf_edge0 Index of the first half edge not forming a triangle.
    *    - Polygon is IndexOfPolygonContainingHalfEdge(ihalf_edge0).
    */
  template <typename DTYPE, typename COORD_TYPE0,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1, 
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_polygon_minus_two_edges_triangulationV_angle
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, const ITYPE1 tri_vertex_index, 
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const VERTEX_INDEX_TYPE ivX =
      mesh.PolygonVertex(ipoly, tri_vertex_index);
    const COORD_TYPE0 * vcoordX = vcoord + ivX*dimension;

    const NUMBER_TYPE nume = mesh.NumPolygonEdges(ipoly);

    // Initialize
    cos_min_angle = -1;
    flag_zero = false;

    HALF_EDGE_INDEX_TYPE ihalf_edge = 
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    
    for (NUMBER_TYPE i = 2; i < nume; i++) {

      ihalf_edge = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge);

      const VERTEX_INDEX_TYPE iv0 = mesh.FromVertexIndex(ihalf_edge);
      const VERTEX_INDEX_TYPE iv1 = mesh.ToVertexIndex(ihalf_edge);

      if (ivX == iv0 || ivX == iv1) { 
        // Skip if ivX equals iv0 or iv1.
        continue; 
      }

      const COORD_TYPE0 * vcoord0 = vcoord + iv0*dimension;
      const COORD_TYPE0 * vcoord1 = vcoord + iv1*dimension;

      COS_TYPE cosX;
      bool flagX_zero;
      IJK::compute_cos_min_triangle_angle
        (dimension, vcoordX, vcoord0, vcoord1, max_small_magnitude, 
         cosX, flagX_zero);

      if (flagX_zero) { flag_zero = true; }
      if (cosX > cos_min_angle) { cos_min_angle = cosX; }
    }

  }


  /*!
   *  @brief Compute cosine of the min angle in interior vertex triangulation.
   *  - Triangulation formed by vcoordX[] and every edge of a polygon
   *    except ihalf_edge0 and NextHalfEdgeInPolygon(ihalf_edge0).
   *  @param ihalf_edge0 Index of the first half edge not forming a triangle.
   *    - Polygon is IndexOfPolygonContainingHalfEdge(ihalf_edge0).
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE0, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_polygon_minus_two_edges_triangulationX_angle
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const COORD_TYPE1 vcoordX[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0);

    const NUMBER_TYPE nume = mesh.NumPolygonEdges(ipoly);

    // Initialize
    cos_min_angle = -1;
    flag_zero = false;

    HALF_EDGE_INDEX_TYPE ihalf_edge = 
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    
    for (NUMBER_TYPE i = 2; i < nume; i++) {

      ihalf_edge = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge);

      const VERTEX_INDEX_TYPE iv0 = mesh.FromVertexIndex(ihalf_edge);
      const VERTEX_INDEX_TYPE iv1 = mesh.ToVertexIndex(ihalf_edge);

      const COORD_TYPE0 * vcoord0 = vcoord + iv0*dimension;
      const COORD_TYPE0 * vcoord1 = vcoord + iv1*dimension;

      COS_TYPE cosX;
      bool flagX_zero;
      IJK::compute_cos_min_triangle_angle
        (dimension, vcoordX, vcoord0, vcoord1, max_small_magnitude, 
         cosX, flagX_zero);

      if (flagX_zero) { flag_zero = true; }
      if (cosX > cos_min_angle) { cos_min_angle = cosX; }
    }

  }

  /*!
   *  @brief Compute cosine of the min angle in interior vertex triangulation.
   *  - Triangulation formed by vcoordX[] and every edge of a polygon
   *    except two half edges inciden on the jloc vertex of the polygon.
   *  @param ipoly Index of polygon to be triangulated.
   *  @param jloc Location in PolygonVertexList() of vertex.
   *    - Edges incident on vertex at jloc are not include in triangulation.
   */
  /// Compute cosine of the min angle in the set of triangles
  ///   formed by vcoordX[] and every edge of a polygon except
  ///   the two edges incident on the j'th vertex of the polygon.
  /// @param jloc Location of vertex to be removed.
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_polygon_minus_vertex_triangulation_angle
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const COORD_TYPE1 vcoordX[], const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const ITYPE1 jloc,
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const HALF_EDGE_INDEX_TYPE ihalf_edge0 = 
      mesh.HalfEdgeIndex(ipoly, jloc);
    const NUMBER_TYPE nume = mesh.NumPolygonEdges(ipoly);

    // Initialize
    cos_min_angle = -1;
    flag_zero = false;

    HALF_EDGE_INDEX_TYPE ihalf_edge = ihalf_edge0;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    IJK::print_coord3D(cerr, "  vcoordX: ", vcoordX, "\n");
    cerr << "  Polygon " << ipoly << ": ";
    IJK::print_list(cerr, mesh.PolygonVertexList(ipoly),
                    mesh.NumPolygonVertices(ipoly));
    cerr << endl;
    cerr << "  jloc: " << jloc 
         << "  vertex: " << mesh.PolygonVertex(ipoly,jloc)
         << endl;
    */

    for (NUMBER_TYPE i = 2; i < nume; i++) {

      ihalf_edge = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge);

      const VERTEX_INDEX_TYPE iv0 = mesh.FromVertexIndex(ihalf_edge);
      const VERTEX_INDEX_TYPE iv1 = mesh.ToVertexIndex(ihalf_edge);

      const COORD_TYPE0 * vcoord0 = vcoord + iv0*dimension;
      const COORD_TYPE0 * vcoord1 = vcoord + iv1*dimension;

      COS_TYPE cosX;
      bool flagX_zero;
      IJK::compute_cos_min_triangle_angle
        (dimension, vcoordX, vcoord0, vcoord1, max_small_magnitude, 
         cosX, flagX_zero);

      if (flagX_zero) { flag_zero = true; }
      if (cosX > cos_min_angle) { cos_min_angle = cosX; }
    }

  }

  ///@}


  // *****************************************************************
  /// @name Compute polygon triangulation that maximizes the min angle.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute fan triangulation vertex which maximizes min triangulation angle.
   *  - Compute the cosine of the max min triangulation angle.
   *  - Version that has input parameter mesh.
   */
  template <typename DTYPE, typename CTYPE, typename MESH_TYPE,
            typename ITYPE0, typename ITYPE1, typename MTYPE,
            typename COS_TYPE>
  void compute_fan_triangulation_to_max_min_angle_M
  (const DTYPE dimension, const CTYPE vcoord[], const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   ITYPE1 & tri_vertex_index, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    tri_vertex_index = 0;
    compute_cos_min_polygon_fan_triangulation_angle
      (dimension, vcoord, mesh, ipoly, max_small_magnitude, 
       0, cos_min_angle, flag_zero);
       
    for (ITYPE1 i = 1; i < mesh.NumPolygonVertices(ipoly); i++) {
      COS_TYPE cosA;
      bool flagA_zero;

      compute_cos_min_polygon_fan_triangulation_angle
        (dimension, vcoord, mesh, ipoly, max_small_magnitude, 
         i, cosA, flagA_zero);

      if (!flagA_zero) {
        if (flag_zero || cos_min_angle > cosA) {
          cos_min_angle = cosA;
          flag_zero = false;
          tri_vertex_index = i;
        }
      }
    }
  }


  /*!
   *  @brief Compute fan triangulation vertex which maximizes min triangulation angle.
   *  - Compute the cosine of the max min triangulation angle.
   *  - Version that has input parameter mesh.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT_E.
   */
  template <typename DTYPE, typename CTYPE, typename MESH_TYPE,
            typename ITYPE0, typename ITYPE1, typename MTYPE,
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_fan_triangulation_to_max_min_angle_M
  (const DTYPE dimension, const CTYPE vcoord[], const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,ITYPE1> & result)
  {
    compute_fan_triangulation_to_max_min_angle_M
      (dimension, vcoord, mesh, ipoly, max_small_magnitude,
       result.tri_vertex_index,
       result.cos_min_triangulation_angle,
       result.flag_zero);

    result.num_interior_tri_vertices = 0;
    result.triangulation_encoding.SetFan
      (result.tri_vertex_index, mesh.NumPolygonVertices(ipoly));
  }


  /*!
   *  @brief Compute the triangulation that maximizes the min polygon triangulation angle.
   *  - Consider fan triangulation (triangulation from a single polygon vertex)
   *      or from an interior vertex.
   *  - Version that has input parameter mesh.
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE,
            typename ITYPE0, typename ITYPE1, typename NTYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_fan_or_interior_vertex_triangulation_to_max_min_angle_M
  (const DTYPE dimension, 
   const COORD_TYPE0 vcoord[], const COORD_TYPE1 interior_vcoord[],
   const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, NTYPE & num_interior, 
   ITYPE1 & tri_vertex_index, bool & flag_zero)
  {
    COS_TYPE cosA, cosB;
    bool flagA_zero, flagB_zero;
    
    compute_fan_triangulation_to_max_min_angle_M
      (dimension, vcoord, mesh, ipoly, max_small_magnitude,
       tri_vertex_index, cosA, flagA_zero);

    compute_cos_min_polygon_interior_vertex_triangulation_angle
      (dimension, vcoord, interior_vcoord, mesh, ipoly, max_small_magnitude,
       cosB, flagB_zero);

    ITYPE0 index_selected;
    IJK::select_min(cosA, flagA_zero,
                    cosB, flagB_zero,
                    cos_min_angle, index_selected,
                    flag_zero);

    if (index_selected == 0) 
      { num_interior = 0; }
    else 
      { num_interior = 1; }
  }


  /*!
   *  @brief Compute the triangulation that maximizes the min polygon triangulation angle.
   *  - Version using C++ STL vector vertex_coord[].
   *  - Version that has input parameter mesh.
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE,
            typename ITYPE0, typename ITYPE1, typename NTYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_fan_or_interior_vertex_triangulation_to_max_min_angle_M
  (const DTYPE dimension, 
   const std::vector<COORD_TYPE0> & vcoord, const COORD_TYPE1 interior_vcoord[],
   const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, NTYPE & num_interior, 
   ITYPE1 & tri_vertex_index, bool & flag_zero)
  {
    compute_fan_or_interior_vertex_triangulation_to_max_min_angle
      (dimension, IJK::vector2pointer(vcoord), interior_vcoord,
       mesh, ipoly, max_small_magnitude, cos_min_angle, num_interior,
       tri_vertex_index, flag_zero);
  }


  /*!
   *  @brief Compute the triangulation that maximizes the min polygon triangulation angle.
   *  - Consider triangulations from a single polygon vertex
   *      or from an interior vertex.
   *  - Version that has input parameter mesh.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT_E.
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE,
            typename ITYPE0, typename NTYPE,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_fan_or_interior_vertex_triangulation_to_max_min_angle_M
  (const DTYPE dimension, 
   const COORD_TYPE0 vcoord[], const COORD_TYPE1 interior_vcoord[],
   const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, NTYPE> & result)
  {
    compute_fan_or_interior_vertex_triangulation_to_max_min_angle_M
      (dimension, vcoord, interior_vcoord, mesh, ipoly, max_small_magnitude,
       result.cos_min_triangulation_angle,
       result.num_interior_tri_vertices,
       result.tri_vertex_index,
       result.flag_zero);

    if (result.num_interior_tri_vertices == 0) {
      result.triangulation_encoding.SetFan
        (result.tri_vertex_index, mesh.NumPolygonVertices(ipoly));
    }
  }

  ///@}


  // *****************************************************************
  /// @name Compute polygon triangulations when some vertex is replaced.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute cosine of the min angle in the polygon triangulation
   *    that maximizes the min angle where j'th vertex is replaced
   *    by new_vertex_coord[].
   *  - Triangulation is from vertex new_vertex_coord[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_polygon_triangulationA_angle_replace_vertex
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 new_vertex_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 ipoly, const ITYPE1 jloc,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    // Consider only triangulation from new_vertex_coord[].
    compute_cos_min_polygon_minus_vertex_triangulation_angle
      (dimension, vcoord, new_vertex_coord, mesh, ipoly, jloc,
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /*!
   *  Compute cosine of the min angle in the polygon triangulation
   *    from an interior vertex that maximizes the min angle
   *  - j'th vertex is replaced by new_vertex_coord[].
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_polygon_interior_vertex_triangulation_angle_replaceV
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 new_vertex_coord[], const CTYPE2 interior_vcoord[], 
   const MESH_TYPE & mesh, const ITYPE0 ipoly, const ITYPE1 jloc_replace,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const NUMBER_TYPE nume = mesh.NumPolygonEdges(ipoly);
    const ITYPE1 j1 = (jloc_replace+1)%nume;
    const ITYPE1 j2 = (jloc_replace+(nume-1))%nume;
    const VERTEX_INDEX_TYPE iv1 = mesh.PolygonVertex(ipoly, j1);
    const VERTEX_INDEX_TYPE iv2 = mesh.PolygonVertex(ipoly, j2);
    const CTYPE0 * vcoord1 = vcoord + iv1*dimension;
    const CTYPE0 * vcoord2 = vcoord + iv2*dimension;
    

    COS_TYPE cosA, cosB;
    bool flagA_zero, flagB_zero;

    compute_cos_min_polygon_minus_vertex_triangulation_angle
      (dimension, vcoord, interior_vcoord, mesh, ipoly, jloc_replace,
       max_small_magnitude, cos_min_angle, flag_zero);

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord1, new_vertex_coord, interior_vcoord, 
       max_small_magnitude, cosA, flagA_zero);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord2, new_vertex_coord, interior_vcoord, 
       max_small_magnitude, cosB, flagB_zero);

    flag_zero = (flag_zero || flagA_zero || flagB_zero);
    cos_min_angle = std::max(cos_min_angle, cosA);
    cos_min_angle = std::max(cos_min_angle, cosB);
  }

  ///@}


  // *****************************************************************
  /// @name Compute triangulations for split edge polygons, basic routines.
  // *****************************************************************

  ///@{

  /// @brief Compute cosine of the min angle in the polygon triangulation from vertex splitting an edge.
  /// @param ihalf_edge Fan triangulation from vertex splitting half edge ihalf_edge.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_polygon_triangulation_angle_from_splitE_vertex
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 split_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    // Consider only triangulation from split vertex.
    compute_cos_min_polygon_minus_edge_triangulation_angle
      (dimension, vcoord, split_coord, mesh, ihalf_edge,
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute cosine of the min angle in the polygon triangulation from vertex splitting an edge.
   *  - Version using C++ STL vector vcoord[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_polygon_triangulation_angle_from_splitE_vertex
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 split_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_polygon_triangulation_angle_from_splitE_vertex
      (dimension, IJK::vector2pointer(vcoord), split_coord,
       mesh, ihalf_edge, max_small_magnitude, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute cosine of the min angle in the polygon triangulation from vertex splitting an edge.
   *  - Version using C++ STL vector vcoord[].
   *  - Version using C++ STL vector split_coord[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_polygon_triangulation_angle_from_splitE_vertex
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & split_coord, const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_polygon_triangulation_angle_from_splitE_vertex
      (dimension, vcoord, IJK::vector2pointer(split_coord),
       mesh, ihalf_edge, max_small_magnitude, cos_min_angle, flag_zero);
  }


  /// @brief Compute cosine of the min angle in the polygon triangulation from vertex splitting an edge.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_split_edge_two_triangles_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 split_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge, const ITYPE1 iv0,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const VERTEX_INDEX_TYPE iv1 = mesh.FromVertexIndex(ihalf_edge);
    const VERTEX_INDEX_TYPE iv2 = mesh.ToVertexIndex(ihalf_edge);
    const CTYPE0 * vcoord0 = vcoord + iv0*dimension;
    const CTYPE0 * vcoord1 = vcoord + iv1*dimension;
    const CTYPE0 * vcoord2 = vcoord + iv2*dimension;
    COS_TYPE cosA, cosB;
    bool flagA_zero, flagB_zero;

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, split_coord, max_small_magnitude,
       cosA, flagA_zero);
    compute_cos_min_triangle_angle
      (dimension, vcoord0, split_coord, vcoord2, max_small_magnitude,
       cosB, flagB_zero);

    cos_min_angle = std::max(cosA, cosB);
    flag_zero = (flagA_zero || flagB_zero);
  }


  /*!
   *  Compute cosine of the min angle in triangulation of split edge polygon from interior vertex.
   *  @param ihalf_edge Split half edge.
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename COORD_TYPE2,
            typename MESH_TYPE, typename ITYPE0, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_splitE_polygon_interior_vertex_triangulation_angle
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const COORD_TYPE1 split_coord[], const COORD_TYPE2 interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge);
    const VERTEX_INDEX_TYPE jv1 = mesh.ToVertexIndex(ihalf_edge);
    const COORD_TYPE0 * vcoord0 = vcoord + jv0*dimension;
    const COORD_TYPE0 * vcoord1 = vcoord + jv1*dimension;

    COS_TYPE cosA, cosB;
    bool flagA_zero, flagB_zero;

    // Consider only triangulation from split vertex.
    compute_cos_min_polygon_minus_edge_triangulation_angle
      (dimension, vcoord, interior_vcoord, mesh, ihalf_edge,
       max_small_magnitude, cos_min_angle, flag_zero);

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, split_coord, interior_vcoord, max_small_magnitude,
       cosA, flagA_zero);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord1, split_coord, interior_vcoord, max_small_magnitude,
       cosB, flagB_zero);

    flag_zero = (flag_zero || flagA_zero || flagB_zero);
    cos_min_angle = std::max(cos_min_angle, cosA);
    cos_min_angle = std::max(cos_min_angle, cosB);
  }


  /// Compute cosine of the min angle in the triangulation
  ///   of a polygon with a split edge from an interior vertex.
  /// - Triangulation has ear separating vertex FromVertexIndex(ihalf_edge1).
  /// - Other triangles are all incident on interior_vcoord[].
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename COORD_TYPE2,
            typename MESH_TYPE, typename ITYPE0, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_splitE_polygon_triangulation_angle_cut_earI0
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const COORD_TYPE1 split_coord[], const COORD_TYPE2 interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge1, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const ITYPE0 ihalf_edge0 = mesh.IndexOfPrevHalfEdgeInPolygon(ihalf_edge1);
    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE jv1 = mesh.FromVertexIndex(ihalf_edge1);
    const VERTEX_INDEX_TYPE jv2 = mesh.ToVertexIndex(ihalf_edge1);
    const COORD_TYPE0 * vcoord0 = vcoord + jv0*dimension;
    const COORD_TYPE0 * vcoord1 = vcoord + jv1*dimension;
    const COORD_TYPE0 * vcoord2 = vcoord + jv2*dimension;

    COS_TYPE cosA;
    bool flagA_zero;

    // Consider triangles not containing edges ihalf_edge0 or ihalf_edge1.
    compute_cos_min_polygon_minus_two_edges_triangulationX_angle
      (dimension, vcoord, interior_vcoord, mesh, ihalf_edge0,
       max_small_magnitude, cos_min_angle, flag_zero);

    // Compute cos min angle of three other triangles.
    // Three triangles are all incident on split_coord.
    IJK::compute_cos_min_pentagon_triangulation_angle
      (dimension, split_coord, vcoord2, interior_vcoord, vcoord0, vcoord1,
       max_small_magnitude, cosA, flagA_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    IJK::print_coord3D(cerr, "  Interior coord: ", interior_vcoord, "");
    cerr << "  cosA: " << cosA << endl;
    */

    flag_zero = (flag_zero || flagA_zero);
    cos_min_angle = std::max(cos_min_angle, cosA);
  }


  /// Compute cosine of the min angle in the triangulation
  ///   of a polygon with a split edge from an interior vertex.
  /// - Triangulation has ear separating vertex ToVertexIndex(ihalf_edge0).
  /// - Other triangles are all incident on interior_vcoord[].
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename COORD_TYPE2,
            typename MESH_TYPE, typename ITYPE0, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_splitE_polygon_triangulation_angle_cut_earI1
  (const DTYPE dimension, const COORD_TYPE0 vcoord[], 
   const COORD_TYPE1 split_coord[], const COORD_TYPE2 interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const ITYPE0 ihalf_edge1 = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE jv1 = mesh.FromVertexIndex(ihalf_edge1);
    const VERTEX_INDEX_TYPE jv2 = mesh.ToVertexIndex(ihalf_edge1);
    const COORD_TYPE0 * vcoord0 = vcoord + jv0*dimension;
    const COORD_TYPE0 * vcoord1 = vcoord + jv1*dimension;
    const COORD_TYPE0 * vcoord2 = vcoord + jv2*dimension;

    COS_TYPE cosA, cosB, cosC;
    bool flagA_zero, flagB_zero, flagC_zero;

    // Consider triangles not containing edges ihalf_edge0 or ihalf_edge1.
    compute_cos_min_polygon_minus_two_edges_triangulationX_angle
      (dimension, vcoord, interior_vcoord, mesh, ihalf_edge0,
       max_small_magnitude, cos_min_angle, flag_zero);

    // Compute cos min angle of three other triangles.
    // Three triangles are all incident on split_coord.
    IJK::compute_cos_min_pentagon_triangulation_angle
      (dimension, split_coord, vcoord1, vcoord2, interior_vcoord, vcoord0,
       max_small_magnitude, cosA, flagA_zero);

    flag_zero = (flag_zero || flagA_zero);
    cos_min_angle = std::max(cos_min_angle, cosA);
  }


  /// Compute cosine of the min angle in the triangulation
  ///   of a polygon with a split edge from an interior vertex.
  /// - Triangulation has ear separating vertex FromVertexIndex(ihalf_edgeB).
  /// - Other triangles are all incident on an original polygon vertex.
  /// @param ihalf_edgeB Half edge containing split vertex.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1, 
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_splitE_polygon_triangulation_angle_cut_earV0
  (const DTYPE dimension, 
   const COORD_TYPE0 vcoord[], const COORD_TYPE1 split_coord[], 
   const MESH_TYPE & mesh, const ITYPE0 ihalf_edgeB, 
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, ITYPE1 & tri_vertex_index, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB);
    const ITYPE1 ilocB = mesh.LocationOfHalfEdgeInPolygon(ihalf_edgeB);
    const NUMBER_TYPE numv = mesh.NumPolygonVertices(ipoly);

    // ilocB is the location of the cut ear.
    // Start at ilocB+1.
    ITYPE1 iloc = (ilocB+1) % numv;

    // Initialize
    cos_min_angle = 1;
    flag_zero = false;
    tri_vertex_index = iloc;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintPolygonIndexAndVertices(cerr, "  ipoly: ", ipoly, ".\n");
    cerr << "  iloc: " << iloc << endl;
    */

    for (NUMBER_TYPE i = 1; i < numv; i++) {

      COS_TYPE cosX;
      bool flagX_zero;

      compute_cos_min_polygon_triangulation_angle_cut_earV0
        (dimension, vcoord, split_coord, mesh, ihalf_edgeB, 
         iloc, max_small_magnitude, cosX, flagX_zero);

      if (!flagX_zero) {
        if (flag_zero || (cosX  < cos_min_angle)) { 
          flag_zero = false;
          cos_min_angle = cosX;
          tri_vertex_index = iloc;
        }
      }

      iloc++;
      if (iloc == numv) { iloc = 0; }
    }

  }


  /// Compute cosine of the min angle in the polygon triangulation 
  ///   that maximizes the min angle.
  /// - Triangulation has one ear.  All other triangles are incident
  ///   on an interior vertex.
  /// - Ear cuts FromVertex(ihalf_edge0) or ToVertex(ihalf_edge0).
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_cos_min_splitE_polygon_triangulation_angle_cut_earI
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 split_coord[], const CTYPE2 interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, 
   const MTYPE max_small_magnitude,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & 
   poly_tri_result)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const NUMBER_TYPE numv = mesh.NumPolygonVertices(ipoly);

    ITYPE1 triB_vertex_index;
    COS_TYPE cosD0, cosD1;
    bool flagD0_zero, flagD1_zero;

    poly_tri_result.Clear();

    compute_cos_min_splitE_polygon_triangulation_angle_cut_earI0
      (dimension, vcoord, split_coord, interior_vcoord, mesh, ihalf_edge0,
       max_small_magnitude, cosD0, flagD0_zero);

    if (!flagD0_zero) {
      if (cosD0 < poly_tri_result.cos_min_triangulation_angle) {
        poly_tri_result.Clear();
        poly_tri_result.SetInterior(cosD0, numv+1, 1);
        const NUMBER_TYPE ifrom = 
          mesh.LocationOfFromVertexInPolygon(ihalf_edge0);

        poly_tri_result.SetEarFlag(ifrom, true);
      }
    }

    compute_cos_min_splitE_polygon_triangulation_angle_cut_earI1
      (dimension, vcoord, split_coord, interior_vcoord, mesh, ihalf_edge0,
       max_small_magnitude, cosD1, flagD1_zero);

    if (!flagD1_zero) {
      if (cosD1 < poly_tri_result.cos_min_triangulation_angle) {
        poly_tri_result.Clear();
        poly_tri_result.SetInterior(cosD1, numv+1, 1);
        const NUMBER_TYPE ito = 
          mesh.LocationOfToVertexInPolygon(ihalf_edge0);

        poly_tri_result.SetEarFlag(ito, true);
      }
    }

  }

  ///@}


  // *****************************************************************
  /// @name Compute triangulations for split edge polygons based on triangulation settings.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute triangulation to max min angle of polygon with split half edge.
   *  - Version using list of pointers to vertex coordinates.
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_half_edge_triangulation_to_max_min_angleP
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 split_coord[], const CTYPE2 interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,ITYPE1> & 
   poly_tri_result)
  {
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef const CTYPE0 * CTYPE_PTR;

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge);
    std::bitset<BIT_SET_SIZE> flag_not_ear;
    std::vector<CTYPE_PTR> poly_vcoord_ptr;

    poly_tri_result.Clear();

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintPolygonIndexAndVertices(cerr, "  Poly ", ipoly, "\n");
    */

    create_split_edge_poly_vertex_list_and_flag_not_ear
      (dimension, vcoord, split_coord, mesh, ihalf_edge,
       poly_vcoord_ptr, flag_not_ear);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "    vertex order in poly_vcoord_ptr[]: ";
    for (int i = 0; i < mesh.NumPolygonVertices(ipoly)+1; i++) {
      cerr << " " << (poly_vcoord_ptr[i] - vcoord)/dimension;
    }
    cerr << endl;
    cerr << "    flag_not_ear: ";
    for (int i = 0; i < mesh.NumPolygonVertices(ipoly)+1; i++) {
      cerr << " " << int(flag_not_ear[i]);
    }
    cerr << endl;
    */


    if (mesh.IsTriangle(ipoly)) {
      // Don't split ears or use interior vertex on triangles.

      POLYGON_TRIANGULATION_SETTINGS triangulation_settings2 =
        triangulation_settings;

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintPolygonIndexAndVertices
        (cerr, "Checking edge split for triangle: ", ipoly, "\n");
      mesh.PrintHalfEdgeIndexAndEndpoints
        (cerr, "  Splitting half edge: ", ihalf_edge, "\n");
      IJK::print_coord3D(cerr, "  split_coord:", split_coord, "\n");
      IJK::print_coord3D(cerr, "  interior_vcoord:", interior_vcoord, "\n");
      */
      
      triangulation_settings2.SetSplitEar(false);
      triangulation_settings2.flag_triangulate_from_interior_vertex = false;
      compute_triangulation_to_max_min_angleP
        (dimension, poly_vcoord_ptr, interior_vcoord, flag_not_ear, max_small_magnitude,
         triangulation_settings2, poly_tri_result);
    }
    else {
      compute_triangulation_to_max_min_angleP
        (dimension, poly_vcoord_ptr, interior_vcoord, flag_not_ear, max_small_magnitude,
         triangulation_settings, poly_tri_result);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "  Triangulation result:" << endl;
      poly_tri_result.Print(cerr, "    ");
      */
    }
  }


  /*!
   *  @brief Compute triangulation to max min angle of polygon with split half edge.
   *  - Version using list of pointers to vertex coordinates.
   *  - Version using C++ STL vector vcoord[].
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_half_edge_triangulation_to_max_min_angleP
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord,
   const CTYPE1 split_coord[], const CTYPE2 interior_vcoord[],
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeA, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,ITYPE1> & 
   poly_tri_result)
  {
    compute_split_half_edge_triangulation_to_max_min_angleP
      (dimension, IJK::vector2pointer(vcoord), 
       split_coord, interior_vcoord,
       mesh, ihalf_edgeA, max_small_magnitude,
       triangulation_settings, poly_tri_result);
  }


  /*!
   *  @brief Compute triangulation to max min angle of polygon with split half edge.
   *  - Version using list of pointers to vertex coordinates.
   *  - Version using C++ STL vector vcoord[].
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_half_edge_triangulation_to_max_min_angleP
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord,
   const std::vector<CTYPE1> & split_coord, 
   const std::vector<CTYPE2> & interior_vcoord, 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeA, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,ITYPE1> & 
   poly_tri_result)
  {
    compute_split_half_edge_triangulation_to_max_min_angleP
      (dimension, vcoord, IJK::vector2pointer(split_coord),
       IJK::vector2pointer(interior_vcoord),
       mesh, ihalf_edgeA, max_small_magnitude,
       triangulation_settings, poly_tri_result);
  }


  /// Compute cosine of the min angle in the polygon triangulation 
  ///   that maximizes the min angle.
  /// @param triangulation_settings 
  ///        Indicates types of triangulations to consider.
  /// @param[out] selected_interior_vcoord[] Selected internal
  ///             vertex coordiates
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename CTYPE3, typename CTYPE4,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_split_half_edge_polygon_triangulation_angleP_IVx2
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 split_coord[], 
   const CTYPE2 interior_vcoordI[], 
   const CTYPE3 interior_vcoordII[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & 
   poly_tri_result,
   std::vector<CTYPE4> & selected_interior_vcoord)
  {
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    typedef const CTYPE0 * CTYPE_PTR;

    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge);
    std::vector<CTYPE_PTR> poly_vcoord_ptr;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << "." << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge: ", ihalf_edgeA, "\n");
    */

    poly_tri_result.Clear();

    create_split_edge_poly_vertex_list
      (dimension, vcoord, split_coord, mesh, ihalf_edge,
       poly_vcoord_ptr);

    if (mesh.IsTriangle(ipoly)) {
      // Don't split ears or use interior vertex on triangles.

      POLYGON_TRIANGULATION_SETTINGS triangulation_settings2 =
        triangulation_settings;

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintPolygonIndexAndVertices
        (cerr, "Checking edge split for triangle (IVx2): ", ipoly, "\n");
      */

      triangulation_settings2.SetSplitEar(false);
      triangulation_settings2.flag_triangulate_from_interior_vertex = false;
      compute_cos_max_min_polygon_triangulation_angleP_IVx2
        (dimension, poly_vcoord_ptr.size(), 
         IJK::vector2pointer(poly_vcoord_ptr), 
         interior_vcoordI, interior_vcoordII, max_small_magnitude, 
         triangulation_settings2, 
         poly_tri_result, selected_interior_vcoord);
    }
    else {
      compute_cos_max_min_polygon_triangulation_angleP_IVx2
        (dimension, poly_vcoord_ptr.size(), 
         IJK::vector2pointer(poly_vcoord_ptr), 
         interior_vcoordI, interior_vcoordII, max_small_magnitude, 
         triangulation_settings, 
         poly_tri_result, selected_interior_vcoord);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintPolygonIndexAndVertices
      (cerr, "Split half edge of polygon: ", 
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA), "\n");
    poly_tri_result.Print(cerr, "  ", numv1);
    */
  }



  /*!
   *  @brief Compute triangulation to max min angle of two polygons with a split edge.
   *  @param triangulationA_settings Indicates types of triangulations
   *    to consider for polygon A.
   *  @param triangulationB_settings Indicates types of triangulations
   *    to consider for polygon B.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_edge_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 split_coord[],
   const CTYPE2 interior_coordA[],
   const CTYPE3 interior_coordB[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulationA_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulationB_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & resultA,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & resultB,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const HALF_EDGE_INDEX_TYPE ihalf_edge1 = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge0);
    bool flagA_zero, flagB_zero;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__;
    mesh.PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "  polygon ", ihalf_edge0, "\n");
      */

    // Compute min angle in triangulation of polygon
    //   containing half_edge0.
    compute_split_half_edge_triangulation_to_max_min_angleP
      (dimension, vcoord, split_coord, interior_coordA, mesh, ihalf_edge0, 
       max_small_magnitude, triangulationA_settings, resultA);

    cos_min_angle = resultA.cos_min_triangulation_angle;
    flag_zero = resultA.flag_zero;

    if (ihalf_edge1 != ihalf_edge0) {
      // Compute min angle in triangulation of polygon
      //   containing half_edge1.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, split_coord, interior_coordB, mesh, ihalf_edge1, 
         max_small_magnitude, triangulationB_settings, resultB);

      cos_min_angle = 
        std::max(cos_min_angle, resultB.cos_min_triangulation_angle);
      flag_zero = (flag_zero || resultB.flag_zero);
    }
  }


  /*!
   *  @brief Compute triangulation to max min angle of two polygons with split edge.
   *  - Version using C++ STL vector vcoord[].
   *  @param triangulationA_settings Indicates types of triangulations
   *    to consider for polygon A.
   *  @param triangulationB_settings Indicates types of triangulations
   *    to consider for polygon B.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_edge_triangulation_to_max_min_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 split_coord[],
   const CTYPE2 interior_coordA[],
   const CTYPE3 interior_coordB[], 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulationA_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulationB_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & resultA,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & resultB,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_split_edge_triangulation_to_max_min_angle
      (dimension, IJK::vector2pointer(vcoord), split_coord,
       interior_coordA, interior_coordB, mesh, ihalf_edge0, 
       max_small_magnitude, 
       triangulationA_settings, triangulationB_settings,
       resultA, resultB, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute triangulation to max min angle of two polygons with split edge.
   *  - Version using C++ STL vector vcoord[].
   *  - Version using C++ STL vector split_coord[].
   *  @param triangulationA_settings Indicates types of triangulations
   *    to consider for polygon A.
   *  @param triangulationB_settings Indicates types of triangulations
   *    to consider for polygon B.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_edge_triangulation_to_max_min_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & split_coord, 
   const std::vector<CTYPE2> & interior_coordA,
   const std::vector<CTYPE3> & interior_coordB, 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edge0, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulationA_settings,
   const POLYGON_TRIANGULATION_SETTINGS & triangulationB_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & resultA,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> & resultB,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_split_edge_triangulation_to_max_min_angle
      (dimension, vcoord, IJK::vector2pointer(split_coord),
       IJK::vector2pointer(interior_coordA), 
       IJK::vector2pointer(interior_coordB),
       mesh, ihalf_edge0, max_small_magnitude, 
       triangulationA_settings, triangulationB_settings,
       resultA, resultB, cos_min_angle, flag_zero);
  }

  ///@}


  // *****************************************************************
  /// @name Compute triangulations for polygons with two split edges.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute triangulation to max min angle of polygon with two split half edges.
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, 
            typename CTYPEV, typename CTYPEA, 
            typename CTYPEB, typename CTYPEI,
            typename MESH_TYPE, typename ITYPEH,
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_split_two_half_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPEV vcoord[], 
   const CTYPEA splitA_coord[], const CTYPEB splitB_coord[], 
   const CTYPEI interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> & 
   poly_tri_result)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef const CTYPEV * CTYPE_PTR;

    const NUMBER_TYPE TWO_SPLIT_EDGES(2);
    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
    const NUMBER_TYPE numv2 = num_poly_vert+TWO_SPLIT_EDGES;
    std::vector<CTYPE_PTR> vcoord_ptr(numv2);
    std::bitset<BIT_SET_SIZE> flag_not_ear;
    IJK::PROCEDURE_ERROR error
      ("compute_split_two_half_edges_triangulation_to_max_min_angle");

    if (!mesh.CheckAreHalfEdgesInSamePolygon
        (ihalf_edgeA, ihalf_edgeB, error))
      { throw error; }

    if (!mesh.CheckAreHalfEdgesDifferent
        (ihalf_edgeA, ihalf_edgeB, error))
      { throw error; }

    poly_tri_result.Clear();

    create_two_split_edges_poly_vertex_list_and_flag_not_ear
      (dimension, vcoord, splitA_coord, splitB_coord, mesh,
         ihalf_edgeA, ihalf_edgeB, vcoord_ptr, flag_not_ear);

    compute_triangulation_to_max_min_angleP
      (dimension, numv2, IJK::vector2pointer(vcoord_ptr),
       interior_vcoord, flag_not_ear, max_small_magnitude,
       triangulation_settings, poly_tri_result);
  }

  /*!
   *  @brief Compute triangulation to max min angle of three polygons with two split edges.
   *  - Split exactly two polygon edges in polygon 1.
   */
  template <typename DTYPE, 
            typename CTYPEV, typename CTYPEA, typename CTYPEB,
            typename CTYPEI0, typename CTYPEI1, typename CTYPEI2,
            typename MESH_TYPE, typename ITYPEH, 
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_split_two_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPEV vcoord[], 
   const CTYPEA splitA_coord[], const CTYPEB splitB_coord[],
   const CTYPEI0 interior_vcoord0[],
   const CTYPEI1 interior_vcoord1[],
   const CTYPEI2 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA1, 
   const ITYPEH ihalf_edgeB1, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> 
   poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const HALF_EDGE_INDEX_TYPE ihalf_edgeA0 = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB2 = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB1);
    bool flagA_zero, flagB_zero;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    triangulation_settings[1].Print(cerr, "  ");
    */

    // Compute min angle in triangulation of polygon
    //   containing ihalf_edgeA1 and ihalf_edgeB1.
    compute_split_two_half_edges_triangulation_to_max_min_angle
      (dimension, vcoord, splitA_coord, splitB_coord, interior_vcoord1, 
       mesh, ihalf_edgeA1, ihalf_edgeB1, max_small_magnitude, 
       triangulation_settings[1], poly_tri_result[1]);

    cos_min_angle = poly_tri_result[1].cos_min_triangulation_angle;
    flag_zero = poly_tri_result[1].flag_zero;

    if (ihalf_edgeA0 != ihalf_edgeA1) {
      // Compute min angle in triangulation of polygon
      //   containing half_edgeA0.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitA_coord, interior_vcoord0, mesh, 
         ihalf_edgeA0, max_small_magnitude, triangulation_settings[0], 
         poly_tri_result[0]);

      cos_min_angle = 
        std::max(cos_min_angle, poly_tri_result[0].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[0].flag_zero);
    }

    if (ihalf_edgeB2 != ihalf_edgeB1) {
      // Compute min angle in triangulation of polygon
      //   containing half_edgeB2.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitB_coord, interior_vcoord2, mesh, 
         ihalf_edgeB2, max_small_magnitude, triangulation_settings[2], 
         poly_tri_result[2]);

      cos_min_angle = 
        std::max(cos_min_angle, poly_tri_result[2].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[2].flag_zero);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  First half edge: ", ihalf_edgeA1, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Second half edge: ", ihalf_edgeB1, "\n");
    cerr << "  Poly 1 cos min angle: " 
         << poly_tri_result[1].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge A0:  ", ihalf_edgeA0, "\n");
    cerr << "  Poly 0 cos min angle: " 
         << poly_tri_result[0].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge B2:  ", ihalf_edgeB2, "\n");
    cerr << "  Poly 2 cos min angle: " 
         << poly_tri_result[2].cos_min_triangulation_angle << endl;
    */
  }


  /*!
   *  Compute cosine of the min angle in the polygon triangulation
   *    that maximizes the min angle.
   *  - Split 1 or 2 polygon edges in polygon 1.
   */
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB,
            typename CTYPEI0, typename CTYPEI1, typename CTYPEI2,
            typename MESH_TYPE, typename ITYPEH,
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_cos_max_min_split_1or2_edges_polyIII_triangulation_angle
  (const DTYPE dimension, const CTYPEV vcoord[], 
   const CTYPEA splitA_coord[], const CTYPEB splitB_coord[],
   const CTYPEI0 interior_vcoord0[],
   const CTYPEI1 interior_vcoord1[],
   const CTYPEI2 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA1, 
   const ITYPEH ihalf_edgeB1, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> 
   poly_tri_result[3],
   bool flag_split_edge[2],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const HALF_EDGE_INDEX_TYPE ihalf_edgeA0 = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA1);
    POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, NTYPER> 
      poly_tri012_result[3];
    POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, NTYPER> 
      poly_tri01_result[2];
    POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, NTYPER> 
      poly_tri12_result[2];
    COS_TYPE cos012, cos01, cos12;
    bool flag012_zero, flag01_zero, flag12_zero;
    
    // Initialize.
    flag_split_edge[0] = false;
    flag_split_edge[1] = false;

    compute_split_two_edges_triangulation_to_max_min_angle
      (dimension, vcoord, splitA_coord, splitB_coord, 
       interior_vcoord0, interior_vcoord1, interior_vcoord2,
       mesh, ihalf_edgeA1, ihalf_edgeB1,
       max_small_magnitude, triangulation_settings,
       poly_tri012_result, cos012, flag012_zero);

    compute_split_edge_triangulation_to_max_min_angle
      (dimension, vcoord, splitA_coord, interior_vcoord0, interior_vcoord1,
       mesh, ihalf_edgeA0, max_small_magnitude,
       triangulation_settings[0], triangulation_settings[1],
       poly_tri01_result[0], poly_tri01_result[1], cos01, flag01_zero);

    compute_split_edge_triangulation_to_max_min_angle
      (dimension, vcoord, splitB_coord, interior_vcoord1, interior_vcoord2,
       mesh, ihalf_edgeB1, max_small_magnitude,
       triangulation_settings[1], triangulation_settings[2],
       poly_tri12_result[0], poly_tri12_result[1], cos12, flag12_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;    
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  First half edge: ", ihalf_edgeA1, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Second half edge: ", ihalf_edgeB1, "\n");
    cerr << "cos012:" << cos01 << endl;
    cerr << "poly_tri012_result[1]:" << endl;
    poly_tri012_result[1].Print(cerr, "  ");
    poly_tri012_result[1].PrintEar(cerr, "  ear: ", "\n");
    cerr << "cos01:" << cos01 << endl;
    poly_tri01_result[0].Print(cerr, "  ");
    poly_tri01_result[1].Print(cerr, "  ");
    cerr << "cos12:" << cos12 << endl;
    poly_tri12_result[0].Print(cerr, "  ");
    poly_tri12_result[1].Print(cerr, "  ");
    */

    NUMBER_TYPE index_selected;

    if (!mesh.IsBoundaryEdge(ihalf_edgeA1) &&
        (poly_tri012_result[0].cos_min_triangulation_angle >=
         poly_tri012_result[1].cos_min_triangulation_angle) &&
        (poly_tri01_result[0].cos_min_triangulation_angle >=
         poly_tri01_result[1].cos_min_triangulation_angle)) {
      // Triangulation of poly 0 always has angle larger than
      //   triangulation of poly 1.
      // Split only ihalf_edgeA1.
      index_selected = 0;

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintHalfEdgeIndexAndEndpoints
        (cerr, "  Split only half edge (A) ", ihalf_edgeA1, ".\n");
      */
    }
    else if (!mesh.IsBoundaryEdge(ihalf_edgeB1) &&
             (poly_tri012_result[2].cos_min_triangulation_angle >=
              poly_tri012_result[1].cos_min_triangulation_angle) &&
             (poly_tri12_result[1].cos_min_triangulation_angle >=
              poly_tri12_result[0].cos_min_triangulation_angle)) {
      // Triangulation of poly 2 always has angle larger than
      //   triangulation of poly 1.
      // Split only half_edgeB1.
      index_selected = 1;

      // *** DEBUG ***
      /*
      using namespace std;
      mesh.PrintHalfEdgeIndexAndEndpoints
        (cerr, "  Split only half edge (B) ", ihalf_edgeB1, ".\n");
      */
    }
    else {
      IJK::select_minIII(cos01, flag01_zero,
                         cos12, flag12_zero,
                         cos012, flag012_zero,
                         cos_min_angle, index_selected, flag_zero);
    }



    /* DEBUG
    IJK::select_minIII(cos01, flag01_zero,
                       cos12, flag12_zero,
                       cos012, flag012_zero,
                       cos_min_angle, index_selected, flag_zero);
    */
    
    switch(index_selected) {

    case 0:
      poly_tri_result[0].Copy(poly_tri01_result[0]);
      poly_tri_result[1].Copy(poly_tri01_result[1]);
      poly_tri_result[2].Clear();
      flag_split_edge[0] = true;
      flag_split_edge[1] = false;
      cos_min_angle = cos01;
      flag_zero = flag01_zero;
      break;

    case 1:
      poly_tri_result[0].Clear();
      poly_tri_result[1].Copy(poly_tri12_result[0]);
      poly_tri_result[2].Copy(poly_tri12_result[1]);
      flag_split_edge[0] = false;
      flag_split_edge[1] = true;
      cos_min_angle = cos12;
      flag_zero = flag12_zero;
      break;

    default:
    case 2:
      // poly_tri_result already set.
      poly_tri_result[0].Copy(poly_tri012_result[0]);
      poly_tri_result[1].Copy(poly_tri012_result[1]);
      poly_tri_result[2].Copy(poly_tri012_result[2]);
      flag_split_edge[0] = true;
      flag_split_edge[1] = true;
      cos_min_angle = cos012;
      flag_zero = flag012_zero;
      break;
    }

  }


  /// Compute cosine of the min angle in the polygon triangulation 
  ///   that maximizes the min angle.
  /// - Split 1 or 2 polygon edges in polygon 1.
  /// - Version using C++ STL vector vcoord[].
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB,
            typename CTYPEI0, typename CTYPEI1, typename CTYPEI2,
            typename MESH_TYPE, typename ITYPEH,
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_cos_max_min_split_1or2_edges_polyIII_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPEV> & vcoord,
   const CTYPEA splitA_coord[], const CTYPEB splitB_coord[],
   const CTYPEI0 interior_vcoord0[],
   const CTYPEI1 interior_vcoord1[],
   const CTYPEI2 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA1, 
   const ITYPEH ihalf_edgeB1, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> 
   poly_tri_result[3],
   bool flag_split_edge[2],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_split_1or2_edges_polyIII_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord),
       splitA_coord, splitB_coord,
       interior_vcoord0, interior_vcoord1, interior_vcoord2,
       mesh, ihalf_edgeA1, ihalf_edgeB1, max_small_magnitude,
       triangulation_settings, poly_tri_result,
       flag_split_edge, cos_min_angle, flag_zero);
  }


  /// Compute cosine of the min angle in the polygon triangulation 
  ///   that maximizes the min angle.
  /// - Split 1 or 2 polygon edges in polygon 1.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector split_coord*[], interior_vcoord*[].
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB,
            typename CTYPEI0, typename CTYPEI1, typename CTYPEI2,
            typename MESH_TYPE, typename ITYPEH,
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_cos_max_min_split_1or2_edges_polyIII_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_coord,
   const std::vector<CTYPEB> & splitB_coord,
   const std::vector<CTYPEI0> & interior_vcoord0,
   const std::vector<CTYPEI1> & interior_vcoord1,
   const std::vector<CTYPEI2> & interior_vcoord2,
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA1, 
   const ITYPEH ihalf_edgeB1, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> 
   poly_tri_result[3],
   bool flag_split_edge[2],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_split_1or2_edges_polyIII_triangulation_angle
      (dimension, vcoord,
       IJK::vector2pointer(splitA_coord), 
       IJK::vector2pointer(splitB_coord),
       IJK::vector2pointer(interior_vcoord0), 
       IJK::vector2pointer(interior_vcoord1), 
       IJK::vector2pointer(interior_vcoord2),
       mesh, ihalf_edgeA1, ihalf_edgeB1, max_small_magnitude,
       triangulation_settings, poly_tri_result,
       flag_split_edge, cos_min_angle, flag_zero);
  }

  ///@}


  // *****************************************************************
  /// @name Compute triangulations for polygons with three split edges.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute triangulation to max min angle of polygon with three split half edges.
   */
  template <typename DTYPE, 
            typename CTYPEV, typename CTYPEA, 
            typename CTYPEB, typename CTYPEI,
            typename MESH_TYPE, typename ITYPEH,
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_split_three_half_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPEV vcoord[], 
   const CTYPEA splitA_coord[], 
   const CTYPEB splitB_coord[], 
   const CTYPEB splitC_coord[], 
   const CTYPEI interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA, 
   const ITYPEH ihalf_edgeB,
   const ITYPEH ihalf_edgeC,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> & 
   poly_tri_result)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef const CTYPEV * CTYPE_PTR;

    const NUMBER_TYPE THREE_SPLIT_EDGES(3);
    const POLYGON_INDEX_TYPE ipoly = 
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NUMBER_TYPE num_poly_vert = mesh.NumPolygonVertices(ipoly);
    const NUMBER_TYPE numv2 = num_poly_vert+THREE_SPLIT_EDGES;
    std::vector<CTYPE_PTR> vcoord_ptr(numv2);
    std::bitset<BIT_SET_SIZE> flag_not_ear;
    IJK::PROCEDURE_ERROR error
      ("compute_split_three_half_edges_triangulation_to_max_min_angle");

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */

    if (!mesh.CheckAreHalfEdgesInSamePolygon
        (ihalf_edgeA, ihalf_edgeB, ihalf_edgeC, error))
      { throw error; }

    if (!mesh.CheckAreHalfEdgesDifferent
        (ihalf_edgeA, ihalf_edgeB, ihalf_edgeC, error))
      { throw error; }

    poly_tri_result.Clear();

    create_three_split_edges_poly_vertex_list_and_flag_not_ear
       (dimension, vcoord, splitA_coord, splitB_coord, splitC_coord,
        mesh, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
        vcoord_ptr, flag_not_ear);
    
    compute_triangulation_to_max_min_angleP
      (dimension, numv2, IJK::vector2pointer(vcoord_ptr), 
       interior_vcoord, flag_not_ear, max_small_magnitude,
       triangulation_settings, poly_tri_result);

  }


  /*!
   *  @brief Compute triangulation to max min angle of four polygons with three split edges.
   *  - Split exactly three polygon edges in polygon 0.
   */
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEI0, typename CTYPEI1, 
            typename CTYPEI2, typename CTYPEI3,
            typename MESH_TYPE, typename ITYPEH, 
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_split_three_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPEV vcoord[], 
   const CTYPEA splitA_coord[], 
   const CTYPEB splitB_coord[],
   const CTYPEC splitC_coord[],
   const CTYPEI0 interior_vcoord0[],
   const CTYPEI1 interior_vcoord1[],
   const CTYPEI2 interior_vcoord2[],
   const CTYPEI3 interior_vcoord3[],
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA,
   const ITYPEH ihalf_edgeB, 
   const ITYPEH ihalf_edgeC,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> 
   poly_tri_result[4],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const HALF_EDGE_INDEX_TYPE ihalf_edgeAX = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeBX = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeCX = 
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeC);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */

    // Compute min angle in triangulation of polygon
    //   containing ihalf_edgeA and ihalf_edgeB and ihalf_edgeC.
    compute_split_three_half_edges_triangulation_to_max_min_angle
      (dimension, vcoord, splitA_coord, splitB_coord, splitC_coord,
       interior_vcoord0, mesh, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       max_small_magnitude, inner_poly_triangulation_settings,
       poly_tri_result[0]);



    cos_min_angle = poly_tri_result[0].cos_min_triangulation_angle;
    flag_zero = poly_tri_result[0].flag_zero;

    if (!mesh.IsBoundaryEdge(ihalf_edgeA)) {
      // Compute min angle in triangulation of polygon
      //   containing half_edgeAX.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitA_coord, interior_vcoord1, mesh, 
         ihalf_edgeAX, max_small_magnitude,
         outer_poly_triangulation_settings,
         poly_tri_result[1]);

      cos_min_angle = 
        std::max(cos_min_angle, poly_tri_result[1].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[1].flag_zero);
    }

    if (!mesh.IsBoundaryEdge(ihalf_edgeB)) {
      // Compute min angle in triangulation of polygon
      //   containing half_edgeB.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitB_coord, interior_vcoord2, mesh, 
         ihalf_edgeBX, max_small_magnitude, 
         outer_poly_triangulation_settings,
         poly_tri_result[2]);

      cos_min_angle = 
        std::max(cos_min_angle, poly_tri_result[2].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[2].flag_zero);
    }

    if (!mesh.IsBoundaryEdge(ihalf_edgeC)) {
      // Compute min angle in triangulation of polygon
      //   containing half_edgeB.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitC_coord, interior_vcoord3, mesh, 
         ihalf_edgeCX, max_small_magnitude, 
         outer_poly_triangulation_settings,
         poly_tri_result[3]);

      cos_min_angle = 
        std::max(cos_min_angle, poly_tri_result[3].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[3].flag_zero);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  First half edge: ", ihalf_edgeA1, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Second half edge: ", ihalf_edgeB1, "\n");
    cerr << "  Poly 1 cos min angle: " 
         << poly_tri_result[1].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge A0:  ", ihalf_edgeA0, "\n");
    cerr << "  Poly 0 cos min angle: " 
         << poly_tri_result[0].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge B2:  ", ihalf_edgeB2, "\n");
    cerr << "  Poly 2 cos min angle: " 
         << poly_tri_result[2].cos_min_triangulation_angle << endl;
    */
  }


  /*!
   *  @brief Compute triangulation to max min angle of four polygons with three split edges.
   *  - Split exactly three polygon edges in polygon 0.
   *  - Version using C++ STL vector vcoord[].
   */
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEI0, typename CTYPEI1, 
            typename CTYPEI2, typename CTYPEI3,
            typename MESH_TYPE, typename ITYPEH, 
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_split_three_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const std::vector<CTYPEV> & vcoord,
   const CTYPEA splitA_coord[], 
   const CTYPEB splitB_coord[],
   const CTYPEC splitC_coord[],
   const CTYPEI0 interior_vcoord0[],
   const CTYPEI1 interior_vcoord1[],
   const CTYPEI2 interior_vcoord2[],
   const CTYPEI3 interior_vcoord3[],
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA,
   const ITYPEH ihalf_edgeB, 
   const ITYPEH ihalf_edgeC,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> 
   poly_tri_result[4],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_split_three_edges_triangulation_to_max_min_angle
      (dimension, IJK::vector2pointer(vcoord),
       splitA_coord, splitB_coord, splitC_coord,
       interior_vcoord0, interior_vcoord1, 
       interior_vcoord2, interior_vcoord3,
       mesh, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       max_small_magnitude,
       inner_poly_triangulation_settings,
       outer_poly_triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute triangulation to max min angle of four polygons with three split edges.
   *  - Split exactly three polygon edges in polygon 0.
   *  - Version using C++ STL vector vcoord[].
   *  - Version using C++ STL vector split_coord*[], interior_vcoord*[].
   */
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEI0, typename CTYPEI1, 
            typename CTYPEI2, typename CTYPEI3,
            typename MESH_TYPE, typename ITYPEH, 
            typename MTYPE, typename COS_TYPE, typename NTYPER,
            int BIT_SET_SIZE>
  void compute_split_three_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const std::vector<CTYPEV> & vcoord,
   const std::vector<CTYPEA> & splitA_coord, 
   const std::vector<CTYPEB> & splitB_coord, 
   const std::vector<CTYPEC> & splitC_coord, 
   const std::vector<CTYPEI0> & interior_vcoord0,
   const std::vector<CTYPEI1> & interior_vcoord1,
   const std::vector<CTYPEI2> & interior_vcoord2,
   const std::vector<CTYPEI3> & interior_vcoord3,
   const MESH_TYPE & mesh, 
   const ITYPEH ihalf_edgeA,
   const ITYPEH ihalf_edgeB, 
   const ITYPEH ihalf_edgeC,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPER> 
   poly_tri_result[4],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_split_three_edges_triangulation_to_max_min_angle
      (dimension, vcoord,
       IJK::vector2pointer(splitA_coord), 
       IJK::vector2pointer(splitB_coord), 
       IJK::vector2pointer(splitC_coord), 
       IJK::vector2pointer(interior_vcoord0),
       IJK::vector2pointer(interior_vcoord1),
       IJK::vector2pointer(interior_vcoord2), 
       IJK::vector2pointer(interior_vcoord3), 
       mesh, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       max_small_magnitude,
       inner_poly_triangulation_settings,
       outer_poly_triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);
  }


  ///@}


  // *****************************************************************
  /// @name Compute triangulations for polygons with three incident split edges.
  // *****************************************************************

  ///@{

  /*!
    *  @brief Compute triangulation to max min angle of four polygons with three split edges incident on the same vertex.
    *  @param inner_poly_triangulation_settings Indicates types
    *    of triangulations to consider for polygons 1 and 2.
    *  @param outer_poly_triangulation_settings Indicates types
    *    of triangulations to consider for polygons 0 and 3.
    */
  template <typename DTYPE, typename CTYPEV,
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_three_incident_edges_triangulation_to_max_min_angle_OBSOLETE
  (const DTYPE dimension, const CTYPEV vcoord[], 
   const CTYPEA splitA_coord[],
   const CTYPEB splitB_coord[],
   const CTYPEC splitC_coord[],
   const CTYPE0 interior_vcoord0[],
   const CTYPE1 interior_vcoord1[],
   const CTYPE2 interior_vcoord2[],
   const CTYPE3 interior_vcoord3[],
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeB1,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> 
   poly_tri_result[4],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE
      NUMBER_TYPE;

    const NUMBER_TYPE FOUR(4);
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

    // Compute min angle in triangulation of polygon 1.
    compute_split_two_half_edges_triangulation_to_max_min_angle
      (dimension, vcoord, splitB_coord, splitC_coord,
       interior_vcoord1, mesh, ihalf_edgeB1, ihalf_edgeC1, 
       max_small_magnitude, inner_poly_triangulation_settings,
       poly_tri_result[1]);

    cos_min_angle = poly_tri_result[1].cos_min_triangulation_angle;
    flag_zero = poly_tri_result[1].flag_zero;

    // Compute min angle in triangulation of polygon 2.
    compute_split_two_half_edges_triangulation_to_max_min_angle
      (dimension, vcoord, splitA_coord, splitB_coord,
       interior_vcoord2, mesh, ihalf_edgeA2, ihalf_edgeB2, 
       max_small_magnitude, inner_poly_triangulation_settings,
       poly_tri_result[2]);

    cos_min_angle = 
      std::max(cos_min_angle, poly_tri_result[2].cos_min_triangulation_angle);
    flag_zero = (flag_zero || poly_tri_result[2].flag_zero);

    if (!mesh.IsBoundaryEdge(ihalf_edgeC1)) {
      // Compute min angle in triangulation of polygon 0.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitC_coord, interior_vcoord0, mesh, 
         ihalf_edgeC0, max_small_magnitude, 
         outer_poly_triangulation_settings, poly_tri_result[0]);

      cos_min_angle = 
        std::max(cos_min_angle, poly_tri_result[0].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[0].flag_zero);
    }

    if (!mesh.IsBoundaryEdge(ihalf_edgeA2)) {
      // Compute min angle in triangulation of polygon 3.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitA_coord, interior_vcoord3, mesh, 
         ihalf_edgeA3, max_small_magnitude, 
         outer_poly_triangulation_settings, poly_tri_result[3]);

      cos_min_angle = 
        std::max(cos_min_angle, poly_tri_result[3].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[3].flag_zero);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  First half edge: ", ihalf_edgeA1, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Second half edge: ", ihalf_edgeB1, "\n");
    cerr << "  Poly 1 cos min angle: " 
         << poly_tri_result[1].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge A0:  ", ihalf_edgeA0, "\n");
    cerr << "  Poly 0 cos min angle: " 
         << poly_tri_result[0].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge B2:  ", ihalf_edgeB2, "\n");
    cerr << "  Poly 2 cos min angle: " 
         << poly_tri_result[2].cos_min_triangulation_angle << endl;
    */

  }


  // *** NEW VERSION ***
  /*!
    *  @brief Compute triangulation to max min angle of four polygons with three split edges incident on the same vertex.
    *  @param inner_poly_triangulation_settings Indicates types
    *    of triangulations to consider for polygons 1 and 2.
    *  @param outer_poly_triangulation_settings Indicates types
    *    of triangulations to consider for polygons 0 and 3.
    */
  template <typename DTYPE, typename CTYPEV,
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_three_incident_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPEV vcoord[],
   const CTYPEA splitA_coord[],
   const CTYPEB splitB_coord[],
   const CTYPEC splitC_coord[],
   const CTYPE0 interior_vcoord0[],
   const CTYPE1 interior_vcoord1[],
   const CTYPE2 interior_vcoord2[],
   const CTYPE3 interior_vcoord3[],
   const MESH_TYPE & mesh,
   const ITYPE0 ihalf_edgeB1,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1>
   poly_tri_result[4],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE
      NUMBER_TYPE;

    const NUMBER_TYPE FOUR(4);
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

    // Compute min angle in triangulation of polygon 1.
    compute_split_two_half_edges_triangulation_to_max_min_angle
      (dimension, vcoord, splitB_coord, splitC_coord,
       interior_vcoord1, mesh, ihalf_edgeB1, ihalf_edgeC1,
       max_small_magnitude, inner_poly_triangulation_settings,
       poly_tri_result[1]);

    cos_min_angle = poly_tri_result[1].cos_min_triangulation_angle;
    flag_zero = poly_tri_result[1].flag_zero;

    // Compute min angle in triangulation of polygon 2.
    compute_split_two_half_edges_triangulation_to_max_min_angle
      (dimension, vcoord, splitA_coord, splitB_coord,
       interior_vcoord2, mesh, ihalf_edgeA2, ihalf_edgeB2,
       max_small_magnitude, inner_poly_triangulation_settings,
       poly_tri_result[2]);

    cos_min_angle =
      std::max(cos_min_angle, poly_tri_result[2].cos_min_triangulation_angle);
    flag_zero = (flag_zero || poly_tri_result[2].flag_zero);

    if (!mesh.IsBoundaryEdge(ihalf_edgeC1)) {
      // Compute min angle in triangulation of polygon 0.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitC_coord, interior_vcoord0, mesh,
         ihalf_edgeC0, max_small_magnitude,
         outer_poly_triangulation_settings, poly_tri_result[0]);

      cos_min_angle =
        std::max(cos_min_angle, poly_tri_result[0].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[0].flag_zero);
    }

    if (!mesh.IsBoundaryEdge(ihalf_edgeA2)) {
      // Compute min angle in triangulation of polygon 3.
      compute_split_half_edge_triangulation_to_max_min_angleP
        (dimension, vcoord, splitA_coord, interior_vcoord3, mesh,
         ihalf_edgeA3, max_small_magnitude,
         outer_poly_triangulation_settings, poly_tri_result[3]);

      cos_min_angle =
        std::max(cos_min_angle, poly_tri_result[3].cos_min_triangulation_angle);
      flag_zero = (flag_zero || poly_tri_result[3].flag_zero);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  First half edge: ", ihalf_edgeA1, "\n");
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Second half edge: ", ihalf_edgeB1, "\n");
    cerr << "  Poly 1 cos min angle: "
         << poly_tri_result[1].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge A0:  ", ihalf_edgeA0, "\n");
    cerr << "  Poly 0 cos min angle: "
         << poly_tri_result[0].cos_min_triangulation_angle << endl;
    mesh.PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Half edge B2:  ", ihalf_edgeB2, "\n");
    cerr << "  Poly 2 cos min angle: "
         << poly_tri_result[2].cos_min_triangulation_angle << endl;
    */

  }


  /*!
   *  @brief Compute triangulation to max min angle of four polygons with three split edges incident on the same vertex.
   *  - Version using C++ STL vector vcoord[].
   */
  template <typename DTYPE, typename CTYPEV,
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_three_incident_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const std::vector<CTYPEV> & vcoord, 
   const CTYPEA splitA_coord[],
   const CTYPEB splitB_coord[],
   const CTYPEC splitC_coord[],
   const CTYPE0 interior_vcoord0[],
   const CTYPE1 interior_vcoord1[],
   const CTYPE2 interior_vcoord2[],
   const CTYPE3 interior_vcoord3[],
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeB1, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> poly_tri_result[4],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_split_three_incident_edges_triangulation_to_max_min_angle
      (dimension, IJK::vector2pointer(vcoord), 
       splitA_coord, splitB_coord, splitC_coord,
       interior_vcoord0, interior_vcoord1, 
       interior_vcoord2, interior_vcoord3,
       mesh, ihalf_edgeB1, max_small_magnitude, 
       inner_poly_triangulation_settings, outer_poly_triangulation_settings,
       poly_tri_result,
       cos_min_angle, flag_zero);
  }

  /*!
   *  @brief Compute triangulation to max min angle of four polygons with three split edges incident on the same vertex.
   *  - Version using C++ STL vector vcoord[].
   *  - Version using C++ STL vector split*_coord[], interior_coord*[].
   */
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_split_three_incident_edges_triangulation_to_max_min_angle
  (const DTYPE dimension, const std::vector<CTYPEV> & vcoord, 
   const std::vector<CTYPEA> & splitA_coord,
   const std::vector<CTYPEB> & splitB_coord,
   const std::vector<CTYPEC> & splitC_coord,
   const std::vector<CTYPE0> & interior_vcoord0, 
   const std::vector<CTYPE1> & interior_vcoord1, 
   const std::vector<CTYPE2> & interior_vcoord2, 
   const std::vector<CTYPE3> & interior_vcoord3, 
   const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeB1, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & inner_poly_triangulation_settings,
   const POLYGON_TRIANGULATION_SETTINGS & outer_poly_triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE, ITYPE1> poly_tri_result[4],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_split_three_incident_edges_triangulation_to_max_min_angle
      (dimension, vcoord, 
       IJK::vector2pointer(splitA_coord),
       IJK::vector2pointer(splitB_coord),
       IJK::vector2pointer(splitC_coord),
       IJK::vector2pointer(interior_vcoord0),
       IJK::vector2pointer(interior_vcoord1),
       IJK::vector2pointer(interior_vcoord2),
       IJK::vector2pointer(interior_vcoord3),
       mesh, ihalf_edgeB1, max_small_magnitude, 
       inner_poly_triangulation_settings, outer_poly_triangulation_settings,
       poly_tri_result,
       cos_min_angle, flag_zero);
  }

  ///@}


  // *****************************************************************
  /// @name Compute triangulations for polygons with replaced vertex.
  // *****************************************************************


  /// Compute cosine of the min angle in the polygon triangulation 
  ///   that maximizes the min angle.
  /// - Version using list of pointers to vertex coordinates.
  /// @param triangulation_settings 
  ///   Indicates types of triangulations to consider.
  template <typename DTYPE, typename CTYPEV, typename MESH_TYPE, 
            typename ITYPEP, typename ITYPER,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE,
            typename NTYPE>
  void compute_cos_max_min_replace_vertex_polygon_triangulation_angleP
  (const DTYPE dimension, 
   const CTYPEV vcoord[], 
   const CTYPEV new_triangle_vcoord[],
   const CTYPEV interior_vcoord[], 
   const MESH_TYPE & mesh, 
   const ITYPEP ipoly, 
   const ITYPER iloc_replace, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   poly_tri_result)
  {
    typedef const CTYPEV * CTYPE_PTR;

    std::vector<CTYPE_PTR> poly_vcoord_ptr;
    IJK::PROCEDURE_ERROR error
      ("compute_cos_max_min_replace_vertex_polygon_triangulation_angleP");

    if (!mesh.CheckPolygonLocation(ipoly, iloc_replace, error))
      { throw error; }

    poly_tri_result.Clear();

    create_poly_vertex_list
      (dimension, vcoord, mesh, ipoly, poly_vcoord_ptr);

    poly_vcoord_ptr[iloc_replace] = new_triangle_vcoord;

    // *** DEBUG ***
    /*
    using namespace std;
    for (int i = 0; i < poly_vcoord_ptr.size(); i++) {
      cerr << "poly_vcoord[" << i << "]:";
      IJK::print_coord3D(cerr, "  ", poly_vcoord_ptr[i], "\n");
    }
    cerr << "  iloc_replace: " << iloc_replace << endl;
    */

    compute_triangulation_to_max_min_angleP
      (dimension, poly_vcoord_ptr, interior_vcoord, max_small_magnitude,
       triangulation_settings, poly_tri_result);
  }


  // *****************************************************************
  /// @name Compute triangulations for merged polygons
  // *****************************************************************

  ///@{

  /// @brief Compute cosine of the min angle in the triangulation of two merged polygons.
  /// - Merge polygons by deleting a shared edge.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_merge_two_polygons_triangulation_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 interior_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeA, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    COS_TYPE cosA, cosB;
    bool flagA_zero, flagB_zero;

    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);

    compute_cos_min_polygon_minus_edge_triangulation_angle
      (dimension, vcoord, interior_coord, mesh,
       ihalf_edgeA, max_small_magnitude, cosA, flagA_zero);
    compute_cos_min_polygon_minus_edge_triangulation_angle
      (dimension, vcoord, interior_coord, mesh,
       ihalf_edgeB, max_small_magnitude, cosB, flagB_zero);

    flag_zero = (flagA_zero || flagB_zero);
    cos_min_angle = std::max(cosA, cosB);
  }


  /// @brief Compute cosine of the min angle in the triangulation of two merged polygons.
  /// - Version using C++ STL vector vcoord[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_merge_two_polygons_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 interior_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeA, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_merge_two_polygons_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), interior_coord,
       mesh, ihalf_edgeA, max_small_magnitude, cos_min_angle,
       flag_zero);
  }


  /// @brief Compute cosine of the min angle in the triangulation of two merged polygons.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector interior_coord[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_merge_two_polygons_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & interior_coord, const MESH_TYPE & mesh, 
   const ITYPE0 ihalf_edgeA, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_merge_two_polygons_triangulation_angle
      (dimension, vcoord, IJK::vector2pointer(interior_coord),
       mesh, ihalf_edgeA, max_small_magnitude, cos_min_angle,
       flag_zero);
  }


  /// @brief Compute cosine of the min angle in the triangulation of a triangle merged with a polygon.
  /// - Triangle is jv1 and endpoints of ihalf_edge.
  /// - Polygon is polygon containing ihalf_edgeB.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_merge_triangle_polygon_triangulation_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 interior_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 jv1, const ITYPE1 ihalf_edge, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge);
    const VERTEX_INDEX_TYPE jv2 = mesh.ToVertexIndex(ihalf_edge);
    const CTYPE0 * vcoord0 = vcoord + jv0*dimension;
    const CTYPE0 * vcoord1 = vcoord + jv1*dimension;
    const CTYPE0 * vcoord2 = vcoord + jv2*dimension;
    COS_TYPE cosA, cosB, cosC;
    bool flagA_zero, flagB_zero, flagC_zero;

    compute_cos_min_polygon_minus_edge_triangulation_angle
      (dimension, vcoord, interior_coord, mesh,
       ihalf_edge, max_small_magnitude, cosA, flagA_zero);

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, interior_coord, max_small_magnitude,
       cosB, flagB_zero);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord2, interior_coord, max_small_magnitude,
       cosC, flagC_zero);

    flag_zero = (flagA_zero || flagB_zero || flagC_zero);
    cos_min_angle = std::max(cosA, cosB);
    cos_min_angle = std::max(cos_min_angle, cosC);
  }

  /// @brief Compute cosine of the min angle in the triangulation of a triangle merged with a polygon.
  /// - Version using C++ STL vector vcoord[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_merge_triangle_polygon_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 interior_coord[], const MESH_TYPE & mesh, 
   const ITYPE0 jv0, const ITYPE1 ihalf_edge, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_merge_triangle_polygon_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), interior_coord,
       mesh, jv0, ihalf_edge, max_small_magnitude, cos_min_angle,
       flag_zero);
  }


  /// @brief Compute cosine of the min angle in the triangulation of a triangle merged with a polygon.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector interior_coord[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MESH_TYPE, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_merge_triangle_polygon_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & interior_coord, const MESH_TYPE & mesh, 
   const ITYPE0 jv0, const ITYPE1 ihalf_edge, const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_merge_triangle_polygon_triangulation_angle
      (dimension, vcoord, IJK::vector2pointer(interior_coord),
       mesh, jv0, ihalf_edge, max_small_magnitude, cos_min_angle,
       flag_zero);
  }

  ///@}

}

#endif

