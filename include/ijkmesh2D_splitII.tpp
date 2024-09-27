/*!
 *  @file ijkmesh2D_splitII.tpp
 *  @brief ijk template classes for splitting edges and polygons 
 *    in 2D polygonal mesh data structures.
 *  - Uses half-edges.  Totally distinct from POLYMESH.
 *  - Version setting coordinate vectors.
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

#ifndef _IJKMESH2D_SPLIT_II__
#define _IJKMESH2D_SPLIT_II_

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijklist.tpp"
#include "ijkmesh2D_split.tpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // *****************************************************************
  // Class MESH2D_SPLIT_II_BASE
  // *****************************************************************

  /// Mesh of 2D polygons supporting polygon splits/replacement.
  /// - Version setting coordinate vectors.
  /// \tparam VTYPE Vertex type.  Should be derived from VERTEX_BASE.
  /// \tparam HALFE_TYPE Half edge type.  
  ///         Should be derived from HALF_EDGE_BASE.
  /// \tparam PTYPE POLYGON type.  
  ///         Should be derived from POLYGON_BASE.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  class MESH2D_SPLIT_II_BASE:
    public MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE> {

  public:

    typedef MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>
    MESH2D_SPLIT_BASE_TYPE;

    typedef typename MESH2D_SPLIT_BASE_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH2D_SPLIT_BASE_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH2D_SPLIT_BASE_TYPE::HALF_EDGE_INDEX_TYPE 
    HALF_EDGE_INDEX_TYPE;
    typedef typename MESH2D_SPLIT_BASE_TYPE::NUMBER_TYPE NUMBER_TYPE;


  protected:

    void Init();  /// Initialize

  public:
    /// constructor
    MESH2D_SPLIT_II_BASE() { Init(); };

    /// Copy.
    template <typename MTYPE>
    void Copy(const MTYPE & meshA);

    // Split functions.

    /// Add vertex and split edge contained in two polygon.
    /// - Set vertex coordinates.
    /// @param[out] ihalf_edgeA_new Index of first new half edge that
    ///   replaces ihalf_edgeA.
    ///   - IndexOfNextHalfEdgeInPolygon(ihalf_edgeA_new) is index 
    ///     of the second half edge that replaces ihalf_edgeA.
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5, typename ITYPE6,
              typename CTYPE2, typename CTYPE3>
    void AddVertexSplitPolygonEdgeII
    (const DTYPE dimension, const CTYPE2 split_vcoord[], 
     std::vector<CTYPE3> & vertex_coord,
     const ITYPE2 ihalf_edgeA, 
     ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1,
     ITYPE6 & ihalf_edgeA_new);

    /// Add vertex and split edge contained in two polygon.
    /// - Set vertex coordinates.
    /// - Version using C++ STL vector split_vcoord[].
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5, typename ITYPE6,
              typename CTYPE2, typename CTYPE3>
    void AddVertexSplitPolygonEdgeII
    (const DTYPE dimension, const std::vector<CTYPE2> & split_vcoord,
     std::vector<CTYPE3> & vertex_coord,
     const ITYPE2 ihalf_edgeA, 
     ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1,
     ITYPE6 & ihalf_edgeA_new);

    /// Add vertex and split edge contained in two polygon.
    /// - Set vertex coordinates.
    /// - Version that does not return ihalf_edgeA_new.
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5,
              typename CTYPE2, typename CTYPE3>
    void AddVertexSplitPolygonEdgeII
    (const DTYPE dimension, const CTYPE2 split_vcoord[], 
     std::vector<CTYPE3> & vertex_coord,
     const ITYPE2 ihalf_edgeA,
     ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1);

    /// Add vertex and split edge contained in two polygon.
    /// - Set vertex coordinates.
    /// - Version that does not return ihalf_edgeA_new.
    /// - Version using C++ STL vector split_vcoord[].
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5,
              typename CTYPE2, typename CTYPE3>
    void AddVertexSplitPolygonEdgeII
    (const DTYPE dimension, const std::vector<CTYPE2> & split_vcoord, 
     std::vector<CTYPE3> & vertex_coord,
     const ITYPE2 ihalf_edgeA,
     ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1);

    /// Add two vertices and split two polygon edges.
    /// - Splits half edges in 3 polygons.
    /// - Set vertex coordinates.
    template <typename DTYPE, typename CTYPEA, typename CTYPEB,
              typename CTYPEV, typename ITYPEH, typename ITYPEV,
              typename ITYPEH2>
    void AddTwoVerticesSplitTwoEdgesOfPolygon
    (const DTYPE dimension, 
     const CTYPEA splitA_vcoord[], const CTYPEB splitB_vcoord[],
     std::vector<CTYPEV> & vertex_coord,
     const ITYPEH ihalf_edgeA, 
     const ITYPEH ihalf_edgeB, 
     ITYPEV & iv_splitA, ITYPEV & iv_splitB,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new);

    /// Add two vertices and split two polygon edges.
    /// - Splits half edges in 3 polygons.
    /// - Version using C++ STL vectors splitA_vcoord[] 
    ///   and splitB_vcoord[].
    template <typename DTYPE, typename CTYPEA, typename CTYPEB,
              typename CTYPEV, typename ITYPEH, typename ITYPEV,
              typename ITYPEH2>
    void AddTwoVerticesSplitTwoEdgesOfPolygon
    (const DTYPE dimension, 
     const std::vector<CTYPEA> & splitA_vcoord, 
     const std::vector<CTYPEB> & splitB_vcoord,
     std::vector<CTYPEV> & vertex_coord,
     const ITYPEH ihalf_edgeA, 
     const ITYPEH ihalf_edgeB, 
     ITYPEV & iv_splitA, ITYPEV & iv_splitB,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new);

    /// Add three vertices and split three polygon edges.
    /// - Splits half edges in 4 polygons.
    /// - Set vertex coordinates.
    template <typename DTYPE, 
              typename CTYPEA, typename CTYPEB, typename CTYPEC,
              typename CTYPEV, typename ITYPEH, typename ITYPEV,
              typename ITYPEH2>
    void AddThreeVerticesSplitThreeEdgesOfPolygon
    (const DTYPE dimension, 
     const CTYPEA splitA_vcoord[], const CTYPEB splitB_vcoord[],
     const CTYPEC splitC_vcoord[],
     std::vector<CTYPEV> & vertex_coord,
     const ITYPEH ihalf_edgeA, 
     const ITYPEH ihalf_edgeB, 
     const ITYPEH ihalf_edgeC,
     ITYPEV & iv_splitA, ITYPEV & iv_splitB, ITYPEV & iv_splitC,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
     ITYPEH2 & ihalf_edgeC_new);

    /// Add three vertices and split three polygon edges.
    /// - Splits half edges in 4 polygons.
    /// - Version using C++ STL vectors splitA_vcoord[] 
    ///   and splitB_vcoord[].
    template <typename DTYPE, 
              typename CTYPEA, typename CTYPEB, typename CTYPEC,
              typename CTYPEV, typename ITYPEH, typename ITYPEV,
              typename ITYPEH2>
    void AddThreeVerticesSplitThreeEdgesOfPolygon
    (const DTYPE dimension, 
     const std::vector<CTYPEA> & splitA_vcoord, 
     const std::vector<CTYPEB> & splitB_vcoord,
     const std::vector<CTYPEC> & splitC_vcoord,
     std::vector<CTYPEV> & vertex_coord,
     const ITYPEH ihalf_edgeA, 
     const ITYPEH ihalf_edgeB, 
     const ITYPEH ihalf_edgeC,
     ITYPEV & iv_splitA, ITYPEV & iv_splitB, ITYPEV & iv_splitC,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
     ITYPEH2 & ihalf_edgeC_new);

    /// Add three vertices and split three polygon edges.
    /// - Splits half edges in 4 polygons.
    /// - Version using C++ STL vectors splitA_vcoord[] 
    ///   and splitB_vcoord[].
    /// - Version using arrays iv_split[] and ihalf_edge_new[].
    template <typename DTYPE, 
              typename CTYPEA, typename CTYPEB, typename CTYPEC,
              typename CTYPEV, typename ITYPEH, typename ITYPEV,
              typename ITYPEH2>
    void AddThreeVerticesSplitThreeEdgesOfPolygon
    (const DTYPE dimension, 
     const std::vector<CTYPEA> & splitA_vcoord, 
     const std::vector<CTYPEB> & splitB_vcoord,
     const std::vector<CTYPEC> & splitC_vcoord,
     std::vector<CTYPEV> & vertex_coord,
     const ITYPEH ihalf_edgeA, 
     const ITYPEH ihalf_edgeB, 
     const ITYPEH ihalf_edgeC,
     ITYPEV iv_split[3], ITYPEH2 ihalf_edge_new[3]);

    /// Add vertices and split three edges that are incident
    ///   on the same vertex.
    /// - Set vertex coordinates.
    /// - Modifies four polygons.
    /// @pre ihalf_edgeB1 is not a boundary edge.
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5, typename ITYPE6,
              typename CTYPEA, typename CTYPEB, typename CTYPEC,
              typename CTYPEV>
    void AddVerticesSplitThreeIncidentEdgesIV
    (const DTYPE dimension, 
     const CTYPEA splitA_vcoord[], 
     const CTYPEB splitB_vcoord[],
     const CTYPEC splitC_vcoord[],
     std::vector<CTYPEV> & vertex_coord,
     const ITYPE2 ihalf_edgeA1,
     ITYPE3 & iv_splitA, ITYPE4 & iv_splitB, ITYPE5 & iv_splitC,
     ITYPE6 ipoly[4]);

    /// Add vertices and split three edges that are incident
    ///   on the same vertex.
    /// - Set vertex coordinates.
    /// - Modifies four polygons.
    /// - Version using C++ STL vector splitA_vcoord[] and splitB_vcoord[].
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5, typename ITYPE6,
              typename CTYPEA, typename CTYPEB, typename CTYPEC,
              typename CTYPEV>
    void AddVerticesSplitThreeIncidentEdgesIV
    (const DTYPE dimension, 
     const std::vector<CTYPEA> & splitA_vcoord, 
     const std::vector<CTYPEB> & splitB_vcoord, 
     const std::vector<CTYPEC> & splitC_vcoord, 
     std::vector<CTYPEV> & vertex_coord,
     const ITYPE2 ihalf_edgeA1,
     ITYPE3 & iv_splitA, ITYPE4 & iv_splitB, ITYPE5 & iv_splitC,
     ITYPE6 ipoly[4]);

    /// Add vertices and split triangle and its edges.
    /// - Set vertex coordinates.
    /// - Modifies four polygons (one triangle and three adjacent polygons.)
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5,
              typename CTYPE2, typename CTYPE3>
    void AddVerticesSplitTriangleAndEdgesIV
    (const DTYPE dimension, 
     const CTYPE2 splitA_vcoord[], 
     const CTYPE2 splitB_vcoord[], 
     const CTYPE2 splitC_vcoord[], 
     std::vector<CTYPE3> & vertex_coord,
     const ITYPE2 ihalf_edgeA,
     ITYPE3 iv_split[3], ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3]);

    /// Add vertices and split triangle and its edges.
    /// - Set vertex coordinates.
    /// - Version using C++ STL vector split*_vcoord[] 
    ///   and interior_coord*[].
    template <typename DTYPE, typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5,
              typename CTYPE2, typename CTYPE3>
    void AddVerticesSplitTriangleAndEdgesIV
    (const DTYPE dimension, 
     const std::vector<CTYPE2> & splitA_vcoord, 
     const std::vector<CTYPE2> & splitB_vcoord, 
     const std::vector<CTYPE2> & splitC_vcoord, 
     std::vector<CTYPE3> & vertex_coord,
     const ITYPE2 ihalf_edgeA,
     ITYPE3 iv_split[3], ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3]);

    /// Add new triangle vertex and replace triangle with triangle edge.
    /// - Set vertex coordinates.
    /// @param ihalf_edge0 Base edge of triangle.
    ///   - Replace vertex opposite ihalf_edge0 with iv_new.
    /// @param[out] Replace ihalf_edge0 with ihalf_edge_new0.
    template <typename DTYPE, typename CTYPE2, typename CTYPE3, 
              typename ITYPE2, typename ITYPE3>
    void AddVertexReplaceTriangleWithTriangleEdge
    (const DTYPE dimension, const CTYPE2 new_triangle_vcoord[], 
     std::vector<CTYPE3> & vertex_coord, 
     const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
     ITYPE3 & ihalf_edge_new0);

    /// Add new triangle vertex and replace triangle with triangle edge.
    /// - Set vertex coordinates.
    /// - Version using C++ STL vector new_triangle_vcoord[].
    template <typename DTYPE, typename CTYPE2, typename CTYPE3, 
              typename ITYPE2, typename ITYPE3>
    void AddVertexReplaceTriangleWithTriangleEdge
    (const DTYPE dimension, const std::vector<CTYPE2> & new_triangle_vcoord,
     std::vector<CTYPE3> & vertex_coord,
     const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
     ITYPE3 & ihalf_edge_new0);



    /// Add vertex and replace split edge triangle with triangle edge.
    /// - Triangle has two split edges.
    /// - Set vertex coordinates.
    template <typename DTYPE, typename CTYPE2, typename CTYPE3,
              typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
    (const DTYPE dimension, const CTYPE2 new_triangle_vcoord[], 
     std::vector<CTYPE3> & vertex_coord, 
     const ITYPE2 ihalf_edge0, ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0);

    /// Add vertex and replace split edge triangle with triangle edge.
    /// - Triangle has two split edges.
    /// - Set vertex coordinates.
    /// - Version using C++ STL vector new_triangle_vcoord[].
    template <typename DTYPE, typename CTYPE2, typename CTYPE3,
              typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
    (const DTYPE dimension, const std::vector<CTYPE2> & new_triangle_vcoord,
     std::vector<CTYPE3> & vertex_coord, 
     const ITYPE2 ihalf_edge0, ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0);

  };


  // *****************************************************************
  // Class MESH2D_A
  // *****************************************************************

  /// Mesh of 2D polygons.
  /// - Use VERTEX_BASE, HALF_EDGE_BASE, POLYGON_BASE
  ///   to represent vertices, half edges and polygons.
  template <typename ITYPE, typename NTYPE>
  class MESH2D_SPLIT_IIA:
    public MESH2D_SPLIT_II_BASE<VERTEX_BASE<ITYPE,ITYPE,NTYPE>, 
                                HALF_EDGE_BASE<ITYPE,ITYPE,ITYPE>,
                                POLYGON_BASE,ITYPE,NTYPE>
  {
  public:
    MESH2D_SPLIT_IIA() {};
  };


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE member functions
  // *****************************************************************

  // constructor
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  void MESH2D_SPLIT_II_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  Init()
  {
    // Nothing to do, so far.
  }

  
  // Copy.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename MTYPE>
  void MESH2D_SPLIT_II_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  Copy(const MTYPE & meshA)
  {
    MESH2D_SPLIT_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::Copy(meshA);
  }



  // *****************************************************************
  // Class MESH2D_SPLIT_II_BASE Split member functions
  // *****************************************************************

  // Add vertex and split edge contained in two polygon.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5, typename ITYPE6,
            typename CTYPE2, typename CTYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexSplitPolygonEdgeII
  (const DTYPE dimension, const CTYPE2 split_vcoord[], 
   std::vector<CTYPE3> & vertex_coord,
   const ITYPE2 ihalf_edgeA, 
   ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1,
   ITYPE6 & ihalf_edgeA_new)
  {
    IJK::insert_coord(dimension, split_vcoord, vertex_coord);
    MESH2D_SPLIT_BASE_TYPE::AddVertexSplitPolygonEdgeII
      (ihalf_edgeA, iv_split, ipoly0, ipoly1, ihalf_edgeA_new);
  }


  // Add vertex and split edge contained in two polygon.
  // - Set vertex coordinates.
  // - Version using C++ STL vector split_vcoord[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5, typename ITYPE6,
            typename CTYPE2, typename CTYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexSplitPolygonEdgeII
  (const DTYPE dimension, const std::vector<CTYPE2> & split_vcoord,
   std::vector<CTYPE3> & vertex_coord,
   const ITYPE2 ihalf_edgeA, 
   ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1,
   ITYPE6 & ihalf_edgeA_new)
  {
    AddVertexSplitPolygonEdgeII
      (dimension, IJK::vector2pointer(split_vcoord),
       vertex_coord, ihalf_edgeA,
       iv_split, ipoly0, ipoly1, ihalf_edgeA_new);
  }


  // Add vertex and split edge contained in two polygon.
  // - Set vertex coordinates.
  // - Version that does not return ihalf_edgeA_new.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5,
            typename CTYPE2, typename CTYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexSplitPolygonEdgeII
  (const DTYPE dimension, const CTYPE2 split_vcoord[], 
   std::vector<CTYPE3> & vertex_coord,
   const ITYPE2 ihalf_edgeA,
   ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1)
  {
    HALF_EDGE_INDEX_TYPE ihalf_edgeA_new;
    AddVertexSplitPolygonEdgeII
      (dimension, split_vcoord, vertex_coord, ihalf_edgeA, 
       iv_split, ipoly0, ipoly1, ihalf_edgeA_new);
  }


  // Add vertex and split edge contained in two polygon.
  // - Set vertex coordinates.
  // - Version that does not return ihalf_edgeA_new.
  // - Version using C++ STL vector split_vcoord[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5,
            typename CTYPE2, typename CTYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexSplitPolygonEdgeII
  (const DTYPE dimension, const std::vector<CTYPE2> & split_vcoord, 
   std::vector<CTYPE3> & vertex_coord,
   const ITYPE2 ihalf_edgeA, 
   ITYPE3 & iv_split, ITYPE4 & ipoly0, ITYPE5 & ipoly1)
  {
    AddVertexSplitPolygonEdgeII
      (dimension, IJK::vector2pointer(split_vcoord), vertex_coord,
       ihalf_edgeA, iv_split, ipoly0, ipoly1);
  }


  // Add two vertices and split two polygon edges.
  // - Splits half edges in 3 polygons.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename CTYPEA, typename CTYPEB,
            typename CTYPEV, typename ITYPEH, typename ITYPEV,
            typename ITYPEH2>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddTwoVerticesSplitTwoEdgesOfPolygon
  (const DTYPE dimension, 
   const CTYPEA splitA_vcoord[], const CTYPEB splitB_vcoord[],
   std::vector<CTYPEV> & vertex_coord,
   const ITYPEH ihalf_edgeA, 
   const ITYPEH ihalf_edgeB, 
   ITYPEV & iv_splitA, ITYPEV & iv_splitB,
   ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new)
  {
    IJK::insert_coord(dimension, splitA_vcoord, vertex_coord);
    IJK::insert_coord(dimension, splitB_vcoord, vertex_coord);

    MESH2D_SPLIT_BASE_TYPE::AddTwoVerticesSplitTwoEdgesOfPolygon
      (ihalf_edgeA, ihalf_edgeB, iv_splitA, iv_splitB, 
       ihalf_edgeA_new, ihalf_edgeB_new);
  }


  // Add two vertices and split two edges of one polygon.
  // - Version using C++ STL vectors splitA_vcoord[] 
  //   and splitB_vcoord[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename CTYPEA, typename CTYPEB,
            typename CTYPEV, typename ITYPEH, typename ITYPEV,
            typename ITYPEH2>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddTwoVerticesSplitTwoEdgesOfPolygon
  (const DTYPE dimension, 
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord,
   std::vector<CTYPEV> & vertex_coord,
   const ITYPEH ihalf_edgeA, 
   const ITYPEH ihalf_edgeB, 
   ITYPEV & iv_splitA, ITYPEV & iv_splitB,
   ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new)
  {
    AddTwoVerticesSplitTwoEdgesOfPolygon
      (dimension, IJK::vector2pointer(splitA_vcoord), 
       IJK::vector2pointer(splitB_vcoord),
       vertex_coord, ihalf_edgeA, ihalf_edgeB, 
       iv_splitA, iv_splitB, ihalf_edgeA_new, ihalf_edgeB_new);
  }


  // Add three vertices and split three polygon edges.
  // - Splits half edges in 4 polygons.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE,
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEV, typename ITYPEH, typename ITYPEV,
            typename ITYPEH2>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddThreeVerticesSplitThreeEdgesOfPolygon
  (const DTYPE dimension, 
   const CTYPEA splitA_vcoord[], const CTYPEB splitB_vcoord[],
   const CTYPEC splitC_vcoord[],
   std::vector<CTYPEV> & vertex_coord,
   const ITYPEH ihalf_edgeA, 
   const ITYPEH ihalf_edgeB, 
   const ITYPEH ihalf_edgeC, 
   ITYPEV & iv_splitA, ITYPEV & iv_splitB, ITYPEV & iv_splitC,
   ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
   ITYPEH2 & ihalf_edgeC_new)
  {
    IJK::insert_coord(dimension, splitA_vcoord, vertex_coord);
    IJK::insert_coord(dimension, splitB_vcoord, vertex_coord);
    IJK::insert_coord(dimension, splitC_vcoord, vertex_coord);

    MESH2D_SPLIT_BASE_TYPE::AddThreeVerticesSplitThreeEdgesOfPolygon
      (ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       iv_splitA, iv_splitB, iv_splitC,
       ihalf_edgeA_new, ihalf_edgeB_new, ihalf_edgeC_new);
  }


  // Add three vertices and split three edges of one polygon.
  // - Version using C++ STL vectors splitA_vcoord[] 
  //   and splitB_vcoord[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEV, typename ITYPEH, typename ITYPEV,
            typename ITYPEH2>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddThreeVerticesSplitThreeEdgesOfPolygon
  (const DTYPE dimension, 
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord,
   const std::vector<CTYPEC> & splitC_vcoord,
   std::vector<CTYPEV> & vertex_coord,
   const ITYPEH ihalf_edgeA, 
   const ITYPEH ihalf_edgeB, 
   const ITYPEH ihalf_edgeC, 
   ITYPEV & iv_splitA, ITYPEV & iv_splitB, ITYPEV & iv_splitC,
   ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
   ITYPEH2 & ihalf_edgeC_new)
  {
    AddThreeVerticesSplitThreeEdgesOfPolygon
      (dimension, 
       IJK::vector2pointer(splitA_vcoord), 
       IJK::vector2pointer(splitB_vcoord),
       IJK::vector2pointer(splitC_vcoord),
       vertex_coord, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       iv_splitA, iv_splitB, iv_splitC,
       ihalf_edgeA_new, ihalf_edgeB_new, ihalf_edgeC_new);
  }


  // Add three vertices and split three edges of one polygon.
  // - Version using C++ STL vectors splitA_vcoord[] 
  //   and splitB_vcoord[].
  // - Version using arrays iv_split[] and ihalf_edge_new[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEV, typename ITYPEH, typename ITYPEV,
            typename ITYPEH2>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddThreeVerticesSplitThreeEdgesOfPolygon
  (const DTYPE dimension, 
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord,
   const std::vector<CTYPEC> & splitC_vcoord,
   std::vector<CTYPEV> & vertex_coord,
   const ITYPEH ihalf_edgeA, 
   const ITYPEH ihalf_edgeB, 
   const ITYPEH ihalf_edgeC, 
   ITYPEV iv_split[3], ITYPEH2 ihalf_edge_new[3])
  {
    AddThreeVerticesSplitThreeEdgesOfPolygon
      (dimension, splitA_vcoord, splitB_vcoord, splitC_vcoord,
       vertex_coord, ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       iv_split[0], iv_split[1], iv_split[2],
       ihalf_edge_new[0], ihalf_edge_new[1], ihalf_edge_new[2]);
  }


  // Add vertices and split three edges that are incident
  //   on the same vertex.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5, typename ITYPE6,
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEV>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVerticesSplitThreeIncidentEdgesIV
  (const DTYPE dimension, 
   const CTYPEA splitA_vcoord[], 
   const CTYPEB splitB_vcoord[],
   const CTYPEC splitC_vcoord[],
   std::vector<CTYPEV> & vertex_coord,
   const ITYPE2 ihalf_edgeA1,
   ITYPE3 & iv_splitA, ITYPE4 & iv_splitB, ITYPE5 & iv_splitC,
   ITYPE6 ipoly[4])
  {
    IJK::insert_coord(dimension, splitA_vcoord, vertex_coord);
    IJK::insert_coord(dimension, splitB_vcoord, vertex_coord);
    IJK::insert_coord(dimension, splitC_vcoord, vertex_coord);
    MESH2D_SPLIT_BASE_TYPE::AddVerticesSplitThreeIncidentEdgesIV
      (ihalf_edgeA1, iv_splitA, iv_splitB, iv_splitC, ipoly);
  }


  // Add vertices and split three edges that are incident
  //   on the same vertex.
  // - Set vertex coordinates.
  // - Version using C++ STL vector splitA_vcoord[] and splitB_vcoord[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5, typename ITYPE6,
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename CTYPEV>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVerticesSplitThreeIncidentEdgesIV
  (const DTYPE dimension, 
   const std::vector<CTYPEA> & splitA_vcoord, 
   const std::vector<CTYPEB> & splitB_vcoord, 
   const std::vector<CTYPEC> & splitC_vcoord, 
   std::vector<CTYPEV> & vertex_coord,
   const ITYPE2 ihalf_edgeA1,
   ITYPE3 & iv_splitA, ITYPE4 & iv_splitB, ITYPE5 & iv_splitC,
   ITYPE6 ipoly[4])
  {
    AddVerticesSplitThreeIncidentEdgesIV
      (dimension, 
       IJK::vector2pointer(splitA_vcoord), 
       IJK::vector2pointer(splitB_vcoord), 
       IJK::vector2pointer(splitC_vcoord), 
       vertex_coord,
       ihalf_edgeA1, iv_splitA, iv_splitB, iv_splitC, ipoly);
  }


  // Add vertices and split triangle and its edges.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5,
            typename CTYPE2, typename CTYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVerticesSplitTriangleAndEdgesIV
  (const DTYPE dimension, 
   const CTYPE2 splitA_vcoord[], 
   const CTYPE2 splitB_vcoord[], 
   const CTYPE2 splitC_vcoord[], 
   std::vector<CTYPE3> & vertex_coord,
   const ITYPE2 ihalf_edgeA,
   ITYPE3 iv_split[3], ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3])
  {
    IJK::insert_coord(dimension, splitA_vcoord, vertex_coord);
    IJK::insert_coord(dimension, splitB_vcoord, vertex_coord);
    IJK::insert_coord(dimension, splitC_vcoord, vertex_coord);

    MESH2D_SPLIT_BASE_TYPE::AddVerticesSplitTriangleAndEdgesIV
      (ihalf_edgeA, iv_split, itriangle_new, ipoly_new);
  }


  // Add vertices and split triangle and its edges.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5,
            typename CTYPE2, typename CTYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVerticesSplitTriangleAndEdgesIV
  (const DTYPE dimension, 
   const std::vector<CTYPE2> & splitA_vcoord, 
   const std::vector<CTYPE2> & splitB_vcoord, 
   const std::vector<CTYPE2> & splitC_vcoord, 
   std::vector<CTYPE3> & vertex_coord,
   const ITYPE2 ihalf_edgeA,
   ITYPE3 iv_split[3], ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3])
  {
    AddVerticesSplitTriangleAndEdgesIV
      (dimension, IJK::vector2pointer(splitA_vcoord),
       IJK::vector2pointer(splitB_vcoord),
       IJK::vector2pointer(splitC_vcoord),
       vertex_coord, ihalf_edgeA, 
       iv_split, itriangle_new, ipoly_new);
  }


  // Add new triangle vertex and replace triangle with triangle edge.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename CTYPE2, typename CTYPE3, 
            typename ITYPE2, typename ITYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexReplaceTriangleWithTriangleEdge
  (const DTYPE dimension, const CTYPE2 new_triangle_vcoord[], 
   std::vector<CTYPE3> & vertex_coord, 
   const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
   ITYPE3 & ihalf_edge_new0)
  {
    IJK::insert_coord(dimension, new_triangle_vcoord, vertex_coord);
    MESH2D_SPLIT_BASE_TYPE::AddVertexReplaceTriangleWithTriangleEdge
      (ihalf_edge0, iv_new, ihalf_edge_new0);
  }


  // Add new triangle vertex and replace triangle with triangle edge.
  // - Set vertex coordinates.
  // - Version using C++ STL vector new_triangle_vcoord[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename CTYPE2, typename CTYPE3, 
            typename ITYPE2, typename ITYPE3>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexReplaceTriangleWithTriangleEdge
  (const DTYPE dimension, const std::vector<CTYPE2> & new_triangle_vcoord, 
   std::vector<CTYPE3> & vertex_coord, 
   const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
   ITYPE3 & ihalf_edge_new0)
  {
    AddVertexReplaceTriangleWithTriangleEdge
      (dimension, IJK::vector2pointer(new_triangle_vcoord), vertex_coord,
       ihalf_edge0, iv_new, ihalf_edge_new0);
  }


  // Add vertex and replace split edge triangle with triangle edge.
  // - Triangle has two split edges.
  // - Set vertex coordinates.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename CTYPE2, typename CTYPE3,
            typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
  (const DTYPE dimension, const CTYPE2 new_triangle_vcoord[], 
   std::vector<CTYPE3> & vertex_coord, 
   const ITYPE2 ihalf_edge0, ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0)
  {
    IJK::insert_coord(dimension, new_triangle_vcoord, vertex_coord);
    MESH2D_SPLIT_BASE_TYPE::AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
      (ihalf_edge0, iv_new, ihalf_edge_new0);

  }


  // Add vertex and replace split edge triangle with triangle edge.
  // - Triangle has two split edges.
  // - Set vertex coordinates.
  /// - Version using C++ STL vector new_triangle_vcoord[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename DTYPE, typename CTYPE2, typename CTYPE3,
            typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
  (const DTYPE dimension, const std::vector<CTYPE2> & new_triangle_vcoord,
   std::vector<CTYPE3> & vertex_coord, 
   const ITYPE2 ihalf_edge0, ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0)
  {
    AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
      (dimension, IJK::vector2pointer(new_triangle_vcoord), vertex_coord,
       ihalf_edge0, iv_new, ihalf_edge_new0);
  }


}


#endif
