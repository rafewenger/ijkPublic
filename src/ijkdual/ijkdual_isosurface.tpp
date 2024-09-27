/*!
 *  \file ijkdual_isosurface.tpp
 *  @brief Data structure storing isosurface.
 *  - Version 0.4.0
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

#ifndef _IJKDUAL_ISOSURFACE_TPP_
#define _IJKDUAL_ISOSURFACE_TPP_

#include <string>
#include <vector>

#include "ijk.tpp"

#include "ijkdual_types.h"


namespace IJKDUAL {

  // ***************************************************************
  // DUAL CONTOURING ISOSURFACE CLASS
  // ***************************************************************

  /*!
   *  @brief Dual contouring isosurface.
   *  - Representation of isosurface returned by Dual Contouring algorithm.
   *  @tparam ISOPOLY_INFO_TYPE Class storing isosurface 
   *    polytope information.
   *    - Should be class derived from IJK::GRID_EDGE or 
   *      IJK::ISODUAL_POLY_INFO.
   */
  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>
  class DUAL_ISOSURFACE_BASE {

  protected:
    DTYPE dimension;
    NTYPE numv_per_isopoly;

    void Init(const int dimension, 
              const ISO_VERTEX_INDEX numv_per_isopoly);

    /// @brief Vertex order of vertices of each quad/cube/hypercube
    ///    in isopoly_vert[].
    CUBE_VERTEX_ORDER cube_vertex_order;

    
  public:

    /// List of isosurface polytope vertices.
    std::vector<ISO_VERTEX_INDEX_TYPE> isopoly_vert;

    /// List of vertex coordinates.
    std::vector<CTYPE> vertex_coord;

    
  public:

    /*!
     *  @brief Constructor.
     *  - Set dimension of space containing the isosurface and
     *    set number of vertices in each isosurface polytope.
     */
    DUAL_ISOSURFACE_BASE
      (const int dimension, const int numv_per_isopoly)
      { Init(dimension, numv_per_isopoly); }

    /// Return dimension of space containing the isosurface.
    int Dimension() const
    { return(dimension); }

    /// Return number of vertices in each isosurface polytope.
    int NumVerticesPerIsoPoly() const
      { return(numv_per_isopoly); };

    /// Return number of isosurface polytopes.
    ISO_VERTEX_INDEX_TYPE NumIsoPoly() const
      { return(isopoly_vert.size()/NumVerticesPerIsoPoly()); };

    /// Return number of isosurface vertices.
    ISO_VERTEX_INDEX_TYPE NumIsoVert() const
    { return(vertex_coord.size()/dimension); }

    /// Return cube_vertex_order.
    CUBE_VERTEX_ORDER CubeVertexOrder()const
    { return cube_vertex_order; }

    /*!
     *  @brief Reorder isosurface quadrilateral vertices in isopoly_vert().
     *  - Returns the new order of the vertices.
     *  - If the vertices of each quad are in COORDINATE_ORDER,
     *    reorder the vertices of each quad to be in CIRCULAR_ORDER.
     *  - If the vertices of each quad are in CIRCULAR_ORDER,
     *    reorder the vertices of each quad to be in COORDINATE_ORDER.
     *  @pre NumVerticesPerIsoPoly() = 4, (i.e., isosurface polytopes are quadrilaterals.)
     */
    CUBE_VERTEX_ORDER ReorderQuadVertices();

    /*!
     *  @brief Set quad vertex order to CUBE_VERTEX_ORDER.
     *  - Reorder vertices in each quad, if necessary.
     *  @pre NumVerticesPerIsoPoly() = 4, (i.e., isosurface polytopes are quadrilaterals.)
     */
    void SetQuadVertexOrder(const CUBE_VERTEX_ORDER vertex_order);

    void Clear();

    /*!
     *  @brief Return false if NumVerticesPerIsoPoly() is not numv.
     *  - Return false and set error message if NumVerticesPerIsoPoly()
     *    is not numv.
     */
    bool CheckNumVerticesPerIsoPoly
      (const int numv, IJK::ERROR & error) const;

    /*!
     *  @brief Return false if isosurface polytopes are not quadrilaterals,
     *     i.e. do not have 4 vertices.
     *  - Return false and set error message if isosurface polytopes
     *    do not have 4 vertices.
     */
    bool CheckIsoPolyAreQuads(IJK::ERROR & error) const;

    /*!
     *  @brief Return false if isopoly_info.size() is incorrect.
     *  - Return false and set error message in error if isopoly_info.size() 
     *   does not match number of poly.
     */
    template <typename _ISOPOLY_INFO_TYPE>
    bool CheckIsoPolyInfo
    (const std::vector<_ISOPOLY_INFO_TYPE> & isopoly_info,
    IJK::ERROR & error) const;
  };


  /*!
   *  @brief Dual contouring isosurface.
   * - Representation of isosurface returned by Dual Contouring algorithm.
   * @tparam ISOPOLY_INFO_TYPE Class storing isosurface polytope information.
   *   - Should be class derived from IJK::GRID_EDGE or 
   *     IJK::ISODUAL_POLY_INFO.
   */
  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>  
  class DUAL_TRIANGULATED_ISOSURFACE_BASE:
    public DUAL_ISOSURFACE_BASE<DTYPE,ISO_VERTEX_INDEX_TYPE,
                                ISOPOLY_INFO_TYPE,CTYPE,NTYPE> {

  protected:
    NTYPE numv_per_simplex;

    /// Initialize.
    void Init(const int dimension);

    /// True if isosurface vertices have been added dual to isosurface polytopes.
    bool flag_isov_dual_to_isopoly;

    /*!
     *  @brief Index of first isosurface vertex dual to an isosurface polytope.
     *  - Used for triangulating isosurface quadrilaterals into four triangles.
     *  - i'th isosurface vertex after first_isov_dual_to_iso_poly
     *    is dual to i'th isosurface polytope.
     */
    ISO_VERTEX_INDEX_TYPE first_isov_dual_to_isopoly;


  public:

    /// @brief List of vertices of each simplex in simplicial mesh.
    /// - Currently only used on 2D surfaces embedded in 3D.
    VERTEX_INDEX_ARRAY simplex_vert;

  public:
    DUAL_TRIANGULATED_ISOSURFACE_BASE
    (const int dimension, const int numv_per_isopoly):
      DUAL_ISOSURFACE_BASE<DTYPE,ISO_VERTEX_INDEX_TYPE,
                           ISOPOLY_INFO_TYPE,CTYPE,NTYPE>
      (dimension, numv_per_isopoly)
      { Init(dimension); };

    /// Return number of vertices per simplex.
    int NumVerticesPerSimplex() const
      { return(numv_per_simplex); };

    /// Return first_isov_dual_to_isopoly.
    ISO_VERTEX_INDEX_TYPE FirstIsovDualToIsopoly() const
    { return first_isov_dual_to_isopoly; }

    /// Return flag_isov_dual_to_isopoly.
    bool FlagIsovDualToIsopoly() const
    { return flag_isov_dual_to_isopoly; }

    /*!
     *  @brief Allocate isov dual to isopoly.
     *  - Sets first_isov_dual_to_isopoly and 
     *    sets flag_isov_dual_to_isopoly to true.
     */
    void AllocateIsovDualToIsopoly();
    
    /// Clear the data structure.
    void Clear();

  };


  // ***************************************************************
  // DUAL_ISOSURFACE_BASE MEMBER FUNCTIONS
  // ***************************************************************

  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>  
  void DUAL_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::
  Init(const int dimension, const ISO_VERTEX_INDEX numv_per_isopoly)
  {
    this->dimension = dimension;
    this->numv_per_isopoly = numv_per_isopoly;

    // Initial vertex order is by coordinates, increasing x,
    //   then increasing y, then ....
    //   - 2D mesh; (0,0), (1,0), (0,1), (1,1).
    // Note: This is NOT the typical circular order used to represent quads (polygons).
    cube_vertex_order = COORDINATE_VERTEX_ORDER;
  }


  // Reorder vertices in each isosurface quadrilateral.
  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>  
  CUBE_VERTEX_ORDER DUAL_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::    
  ReorderQuadVertices()
  {
    IJK::PROCEDURE_ERROR
      error("DUAL_ISOSURFACE_BASE::ReorderQuadVertices()");

    if (!CheckIsoPolyAreQuads(error))
      { throw error; }

    IJK::reorder_quad_vertices(isopoly_vert);

    if (cube_vertex_order == COORDINATE_VERTEX_ORDER)
      { cube_vertex_order = CIRCULAR_VERTEX_ORDER; }
    else {
      // cube_vertex_order == CIRCULAR_VERTEX_ORDER.
      cube_vertex_order = COORDINATE_VERTEX_ORDER;
    }

    return cube_vertex_order;
  }


  // Set quad vertex order to vertex_order.
  // - Reorder quad vertices, if necessary.
  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>    
  void DUAL_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::
  SetQuadVertexOrder(const CUBE_VERTEX_ORDER vertex_order)
  {
    IJK::PROCEDURE_ERROR
      error("DUAL_ISOSURFACE_BASE::SetQuadVertexOrder");

    if (!CheckIsoPolyAreQuads(error))
      { throw error; }

    if (this->cube_vertex_order != vertex_order)
      { ReorderQuadVertices(); }
  }

  
  // Return false if NumVerticesPerIsoPoly() is not numv.
  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>    
  bool DUAL_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::    
  CheckNumVerticesPerIsoPoly(const int numv, IJK::ERROR & error) const
  {
    if (NumVerticesPerIsoPoly() != numv) {
      error.AddMessage
        ("Programming error. Incorrect number of vertices per isosurface polytope.");
      error.AddMessage("Number of vertices per isosurface polytope is ",
                       NumVerticesPerIsoPoly(), ".");
      error.AddMessage("Expected number is ", numv, ".");
      return(false);
    }

    return(true);
  }


  // Return false if isosurface polytopes are not quads,
  //   i.e., do not have 4 vertices.
  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>    
  bool DUAL_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::        
  CheckIsoPolyAreQuads(IJK::ERROR & error) const
  {
    const int NUM_VERTICES_PER_QUAD(4);

    if (NumVerticesPerIsoPoly() != NUM_VERTICES_PER_QUAD) {
      error.AddMessage("Programming error. Isosurface polytopes are not quads.");
      error.AddMessage("Number of vertices per isosurface polytope is ",
                       NumVerticesPerIsoPoly(), ".");
      error.AddMessage("Expected number is 4.");
      return(false);
    }

    return(true);
  }


  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>    
  void DUAL_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::
  Clear()
  {
    isopoly_vert.clear();
    vertex_coord.clear();
  }


  // Return false if isopoly_info.size() is incorrect.
  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>
  template <typename _ISOPOLY_INFO_TYPE>
  bool DUAL_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::    
  CheckIsoPolyInfo
  (const std::vector<_ISOPOLY_INFO_TYPE> & isopoly_info,
   IJK::ERROR & error) const
  {
    typedef typename std::vector<ISOPOLY_INFO_TYPE>::size_type SIZE_TYPE;
    const SIZE_TYPE num_iso_poly = NumIsoPoly();

    if (isopoly_info.size() != num_iso_poly) {
      if (isopoly_info.size() == 0) {
        error.AddMessage
          ("Programming error. Array DUAL_ISOSURFACE_BASE::isopoly_info not allocated.");
        return(false);
      }
      else {
        error.AddMessage
          ("Programming error. Size of array DUAL_ISOSURFACE_BASE::isopoly_info");
        error.AddMessage
          ("  does not match number of isosurface polytopes stored in isopoly_vert.");
        error.AddMessage("  Number of isosurface polytopes: ", num_iso_poly, "");
        error.AddMessage
          ("  Size of array DUAL_ISOSURFACE_BASE::isopoly_info: ", 
           isopoly_info.size(), "");
        return(false);
      }
    }

    return(true);
  }  


  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>    
  void DUAL_TRIANGULATED_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::    
  Init(const int dimension)
  {
    numv_per_simplex = dimension+1;
    flag_isov_dual_to_isopoly = false;
  }


  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>    
  void DUAL_TRIANGULATED_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::    
  AllocateIsovDualToIsopoly()
  {
    typedef typename COORD_ARRAY::size_type SIZE_TYPE;
    
    IJK::PROCEDURE_ERROR
      error("DUAL_TRIANGULATED_ISOSURFACE_BASE::AllocateIsovDualToIsopoly");
    
    if (flag_isov_dual_to_isopoly) {
      error.AddMessage
        ("Programming error. Vertex coordinates already allocated.");
      error.AddMessage
        ("  Coordinates of isosurface vertices dual to isosurface polytopes");
      error.AddMessage
        ("  are already allocated.");
      throw error;
    }

    const SIZE_TYPE new_size =
      this->vertex_coord.size() + this->Dimension()*this->NumIsoPoly();

    first_isov_dual_to_isopoly = this->NumIsoVert();
    this->vertex_coord.resize(new_size);
    flag_isov_dual_to_isopoly = true;
  }


  template <typename DTYPE, typename ISO_VERTEX_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE,
            typename NTYPE>    
  void DUAL_TRIANGULATED_ISOSURFACE_BASE
  <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::    
  Clear()
  {
    DUAL_ISOSURFACE_BASE
      <DTYPE,ISO_VERTEX_INDEX_TYPE,ISOPOLY_INFO_TYPE,CTYPE,NTYPE>::
      Clear();
    
    simplex_vert.clear();
    first_isov_dual_to_isopoly = 0;
    flag_isov_dual_to_isopoly = false;
  }

}

#endif
