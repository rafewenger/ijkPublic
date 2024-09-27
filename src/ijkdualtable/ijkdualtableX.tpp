/*!
 *  \file ijkdualtableX.tpp
 *  @brief Extended dual lookup tables of isosurface and interval volume vertices.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2023 Rephael Wenger

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

#ifndef _IJKDUALTABLE_X_TPP_
#define _IJKDUALTABLE_X_TPP_

#include "ijkdualtable.tpp"
#include "ijkdualtable_vert.tpp"

#include <vector>


namespace IJKDUALTABLE {

  // **************************************************
  // DUAL TABLE VERTEX INFORMATION
  // **************************************************

  /* OBSOLETE/MOVED
  /// Information about each of the isosurface/interval volume vertices
  ///   for each table entry.
  template <typename NTYPE, typename VINFO_TYPE>
  class DUAL_TABLE_VERTEX_INFO:
    public IJK::LIST_OF_LISTS<VINFO_TYPE,NTYPE> {

  public:
    typedef VINFO_TYPE VERTEX_INFO_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:
    DUAL_TABLE_VERTEX_INFO() {};
    template <typename DUAL_TABLE_TYPE>
    DUAL_TABLE_VERTEX_INFO(const DUAL_TABLE_TYPE & table);

    template <typename TI_TYPE>
    NTYPE NumVertices(const TI_TYPE ientry) const
    { return(this->ListLength(ientry)); };

    template <typename DUAL_TABLE_TYPE>
    void Set(const DUAL_TABLE_TYPE & table);

    /// Return reference to vertex info.
    template <typename TI_TYPE>
    const VINFO_TYPE & VertexInfo
    (const TI_TYPE ientry, const NTYPE isov) const
    { return(this->ElementRefConst(ientry, isov)); }

    /// Return non-constant reference to vertex info.
    template <typename TI_TYPE>
    VINFO_TYPE & VertexInfoNC
    (const TI_TYPE ientry, const NTYPE isov)
    { return(this->ElementRef(ientry, isov)); }
  };
  */


  // **************************************************
  // CLASS DUAL_TABLE_VERTEX_INFO MEMBER FUNCTIONS
  // **************************************************

  /* OBSOLETE/MOVED
  template <typename NTYPE, typename VINFO_TYPE>
  template <typename DUAL_TABLE_TYPE>
  DUAL_TABLE_VERTEX_INFO<NTYPE,VINFO_TYPE>::
  DUAL_TABLE_VERTEX_INFO(const DUAL_TABLE_TYPE & table)
  {
    Set(table);
  }

  template <typename NTYPE, typename VINFO_TYPE>
  template <typename DUAL_TABLE_TYPE>
  void DUAL_TABLE_VERTEX_INFO<NTYPE,VINFO_TYPE>::
  Set(const DUAL_TABLE_TYPE & table)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    std::vector<NTYPE> num_isov(table.NumTableEntries());

    for (TABLE_INDEX i = 0; i < table.NumTableEntries(); i++) 
      { num_isov[i] = table.Entry(i).NumVertices(); }
      
    this->CreateLists(num_isov);
  }
  */


  // **************************************************
  // COMPUTE VERTEX DEGREE
  // **************************************************


  /* OBSOLETE/MOVED
  /// Compute number of isosurface polytopes incident on each vertex.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_dual_isotable_vertex_degrees
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & dual_table_vertex_info)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { dual_table_vertex_info.VertexInfoNC(it, j).degree = 0; }
        

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.IsBipolar(it, ie)) {
          const NTYPE isov = table.IncidentIsoVertex(it, ie);
          dual_table_vertex_info.VertexInfoNC(it, isov).degree++;
        }
      }
    }

  }


  /// Compute number of interval volume polytopes incident on each vertex.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_ivoldual_table_num_incident_poly
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & dual_table_vertex_info)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { dual_table_vertex_info.VertexInfoNC(it, j).num_incident_poly = 0; }
        

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.EdgeHasDualIVolPoly(it, ie)) {
          const NTYPE isov0 = table.LowerIncident(it, ie);
          const NTYPE isov1 = table.UpperIncident(it, ie);
          dual_table_vertex_info.VertexInfoNC(it, isov0).num_incident_poly++;
          dual_table_vertex_info.VertexInfoNC(it, isov1).num_incident_poly++;
        }
      }

      for (NTYPE iv = 0; iv < table.NumPolyVertices(); iv++) {
        if (table.IsInIntervalVolume(it, iv)) {
          const NTYPE isov = table.IncidentIVolVertex(it, iv);
          dual_table_vertex_info.VertexInfoNC(it, isov).num_incident_poly++;
        }
      }
    }

  }


  /// Compute number of isosurface polytopes incident on each vertex.
  /// - Note: An isosurface polytope is a facet of an interval volume polytope.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_ivoldual_table_num_incident_isosurface_poly
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::IVOL_VERTEX_INDEX_TYPE IVOLV_TYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    IVOLV_TYPE ivolv;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (IVOLV_TYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { vinfo.VertexInfoNC(it, j).num_incident_isopoly = 0; }

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {

        if (table.DoesLowerIsosurfaceIntersectEdge(it, ie, ivolv))
          { vinfo.VertexInfoNC(it, ivolv).num_incident_isopoly++; }

        if (table.DoesUpperIsosurfaceIntersectEdge(it, ie, ivolv))
          { vinfo.VertexInfoNC(it, ivolv).num_incident_isopoly++; }
      }
    }
  }
  */


  // **************************************************
  // DETERMINE IVOL VERTEX CONNECTION DIRECTIONS
  // **************************************************

  /// Determine interval volume vertex connection directions.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_ivol_vertex_connection_directions
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE
      VINFO_TYPE;
    typedef typename VINFO_TYPE::CUBE_VERTEX_TYPE CUBE_VERTEX_TYPE;
    typedef typename VINFO_TYPE::CUBE_EDGE_TYPE CUBE_EDGE_TYPE;
    typedef typename VINFO_TYPE::DIR_BITS_TYPE DIR_BITS_TYPE;

    const NTYPE NUM_FACETS(table.Cube().NumFacets());
    const NTYPE NUM_EDGES(table.Cube().NumEdges());
    const NTYPE NUM_VERTICES(table.Cube().NumVertices());
    NTYPE ivolv[2], iend[2];
    DIR_BITS_TYPE connect_dir[2];
    std::vector<DIR_BITS_TYPE> bit_mask(NUM_FACETS);


    // Set bit_mask[]
    for (NTYPE i = 0; i < NUM_FACETS; i++) 
      { bit_mask[i] = (DIR_BITS_TYPE(1) << i); }

    // Initialize connect_dir.
    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {
      for (NTYPE ivolv = 0; ivolv < table.NumIVolVertices(it); ivolv++) {
        vinfo.VertexInfoNC(it,ivolv).connect_dir = 0;
      }
    }

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      for (CUBE_EDGE_TYPE ie = 0; ie < NUM_EDGES; ie++) {

        if (table.EdgeHasDualIVolPoly(it, ie)) {
          ivolv[0] = table.LowerIncident(it, ie);
          ivolv[1] = table.UpperIncident(it, ie);
          iend[0] = table.Cube().EdgeEndpoint(ie, 0);
          iend[1] = table.Cube().EdgeEndpoint(ie, 1);
          connect_dir[0] = vinfo.VertexInfo(it, ivolv[0]).connect_dir;
          connect_dir[1] = vinfo.VertexInfo(it, ivolv[1]).connect_dir;
          DTYPE edge_dir = table.Cube().EdgeDir(ie);

          for (DTYPE d = 1; d < table.Dimension(); d++) {
            const DTYPE dir2 = (edge_dir+d)%table.Dimension();

            if ((iend[0] & bit_mask[dir2]) == 0) {
              connect_dir[0] = (connect_dir[0] | bit_mask[dir2]);
              connect_dir[1] = (connect_dir[1] | bit_mask[dir2]);
            }
            else {
              const DTYPE jfacet = dir2 + table.Dimension();
              connect_dir[0] = (connect_dir[0] | bit_mask[jfacet]);
              connect_dir[1] = (connect_dir[1] | bit_mask[jfacet]);
            }
          }

          vinfo.VertexInfoNC(it, ivolv[0]).connect_dir = connect_dir[0];
          vinfo.VertexInfoNC(it, ivolv[1]).connect_dir = connect_dir[1];
        }
      }


      for (CUBE_VERTEX_TYPE iv = 0; iv < NUM_VERTICES; iv++) {

        if (table.IsInIntervalVolume(it, iv)) {
          ivolv[0] = table.IncidentIVolVertex(it, iv);
          connect_dir[0] = vinfo.VertexInfo(it, ivolv[0]).connect_dir;

          for (DTYPE d = 0; d < table.Dimension(); d++) {
            if ((iv & bit_mask[d]) == 0) {
              connect_dir[0] = (connect_dir[0] | bit_mask[d]);
            }
            else {
              const DTYPE jfacet = d + table.Dimension();
              connect_dir[0] = (connect_dir[0] | bit_mask[jfacet]);
            }
          }

          vinfo.VertexInfoNC(it, ivolv[0]).connect_dir = connect_dir[0];
        }
      }

    }
  }


  /// Determine interval volume vertex isosurface connection directions,
  /// i.e., connection directions of edges which lie on the isosurface.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_ivol_vertex_iso_connection_directions
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_TYPE::IVOL_VERTEX_INDEX_TYPE IVOLV_TYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE
      VINFO_TYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename VINFO_TYPE::CUBE_VERTEX_TYPE CUBE_VERTEX_TYPE;
    typedef typename VINFO_TYPE::CUBE_EDGE_TYPE CUBE_EDGE_TYPE;
    typedef typename VINFO_TYPE::DIR_BITS_TYPE DIR_BITS_TYPE;

    const NTYPE NUM_CUBE_FACETS(table.Cube().NumFacets());
    CUBE_VERTEX_TYPE iend[2];
    IVOLV_TYPE ivolv;
    NTYPE jfacet;
    DIR_BITS_TYPE iso_connect_dir;
    std::vector<DIR_BITS_TYPE> bit_mask(NUM_CUBE_FACETS);


    // Set bit_mask[]
    for (NTYPE i = 0; i < NUM_CUBE_FACETS; i++) 
      { bit_mask[i] = (DIR_BITS_TYPE(1) << i); }

    // Initialize connect_dir.
    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {
      for (NTYPE ivolv = 0; ivolv < table.NumIVolVertices(it); ivolv++) {
        vinfo.VertexInfoNC(it,ivolv).iso_connect_dir = 0;
      }
    }

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      for (CUBE_EDGE_TYPE ie = 0; ie < table.NumPolyEdges(); ie++) {

        iend[0] = table.Cube().EdgeEndpoint(ie, 0);
        iend[1] = table.Cube().EdgeEndpoint(ie, 1);
        const DTYPE edge_dir = table.Cube().EdgeDir(ie);

        if (table.DoesLowerIsosurfaceIntersectEdge(it, ie, ivolv)) {

          if (table.IsBelowIntervalVolume(it, iend[1])) 
            { std::swap(iend[0], iend[1]); }

          for (DTYPE d = 1; d < table.Dimension(); d++) {
            const DTYPE dir2 = (edge_dir+d)%table.Dimension();

            if ((iend[0] & bit_mask[dir2]) == 0) { jfacet = dir2; }
            else { jfacet = dir2 + table.Dimension(); }

            iso_connect_dir = vinfo.VertexInfoNC(it, ivolv).iso_connect_dir;
            vinfo.VertexInfoNC(it, ivolv).iso_connect_dir = 
              (iso_connect_dir | bit_mask[jfacet]);
          }
        }

        if (table.DoesUpperIsosurfaceIntersectEdge(it, ie, ivolv)) {

          if (table.IsAboveIntervalVolume(it, iend[0])) 
            { std::swap(iend[0], iend[1]); }

          for (DTYPE d = 1; d < table.Dimension(); d++) {

            const DTYPE dir2 = (edge_dir+d)%table.Dimension();

            if ((iend[1] & bit_mask[dir2]) == 0) { jfacet = dir2; }
            else { jfacet = dir2 + table.Dimension(); }

            iso_connect_dir = vinfo.VertexInfoNC(it, ivolv).iso_connect_dir;
            vinfo.VertexInfoNC(it, ivolv).iso_connect_dir = 
              (iso_connect_dir | bit_mask[jfacet]);
          }
        }
      }
    }
  }


  // **************************************************
  // DETERMINE DOUBLY CONNECTED VERTICES
  // **************************************************

  /// Determine interval volume vertices which are doubly connected
  ///   to some adjacent cube.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_doubly_connected_ivol3D_vertices
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::CUBE_VERTEX_TYPE
      CUBE_VERTEX_TYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {
      const NTYPE num_ivolv = table.Entry(it).NumVertices();

      for (NTYPE j = 0; j < num_ivolv; j++) {
        NTYPE dc_facet;
        if (is_ivol3D_vertex_doubly_connected(table, it, j, dc_facet)) {
          vinfo.VertexInfoNC(it,j).is_doubly_connected = true; 
          vinfo.VertexInfoNC(it,j).doubly_connected_facet = dc_facet;
        }
        else {
          vinfo.VertexInfoNC(it,j).is_doubly_connected = false; 
          vinfo.VertexInfoNC(it,j).doubly_connected_facet = 0;
        }
      }
    }
  }


  // **************************************************
  // DETERMINE SEPARATION VERTICES
  // **************************************************

  /// Return true if vertex with isosurface connection 
  ///   directions iso_connect_dir separates a cube vertex.
  /// @param[out] separation_vertex Separated cube vertex.
  template <typename DTYPE, typename DIR_BITS_TYPE, typename VTYPE>
  bool separates_cube_vertex
  (const DTYPE dimension, const DIR_BITS_TYPE iso_connect_dir, 
   VTYPE & separation_vertex)
  {
    // Initialize
    separation_vertex = 0;

    VTYPE num_iso_connect_dir = 0;

    for (DTYPE d = 0; d < dimension; d++) {
      const DIR_BITS_TYPE mask = (DIR_BITS_TYPE(1) << d);
      const DIR_BITS_TYPE mask2 = (DIR_BITS_TYPE(1) << (d+dimension));
      if ((mask & iso_connect_dir) == 0) {
        if ((mask2 & iso_connect_dir) != 0) {
          separation_vertex = (separation_vertex | mask);
          num_iso_connect_dir++;
        }
      }
      else if ((mask2 & iso_connect_dir) == 0) {
        num_iso_connect_dir++;
      }
      else {
        // Does not separate a cube vertex.
        separation_vertex = 0;
        return(false);
      }
    }

    if (num_iso_connect_dir == dimension) 
      { return(true); }
    else {
      // Does not separate a cube vertex.
      separation_vertex = 0;
      return(false);
    }
  }


  /// Determine separation cube vertex for ivol vertex incident
  ///   incident on dimension isosurface polytopes.
  /// @pre Call compute_ivoldual_table_num_incident_poly
  ///   and compute_ivoldual_table_num_incident_isosurface_poly
  ///   before calling this routine.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_separation_vertices
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE
      VINFO_TYPE;
    typedef typename VINFO_TYPE::CUBE_VERTEX_TYPE CUBE_VERTEX_TYPE;
    typedef typename VINFO_TYPE::DIR_BITS_TYPE DIR_BITS_TYPE;

    const DTYPE dimension = table.Dimension();
    const CUBE_VERTEX_TYPE UNDEFINED_CUBE_VERTEX =
      DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::UNDEFINED_CUBE_VERTEX;

    // Initialize
    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {
      const NTYPE num_ivolv = table.Entry(it).NumVertices();

      for (NTYPE j = 0; j < num_ivolv; j++) {
        const DIR_BITS_TYPE iso_connect_dir =
          vinfo.VertexInfo(it, j).iso_connect_dir;
        CUBE_VERTEX_TYPE separation_vertex;

        if (separates_cube_vertex
            (dimension, iso_connect_dir, separation_vertex)) {
          vinfo.VertexInfoNC(it, j).separation_vertex = separation_vertex;
        }
        else {
          vinfo.VertexInfoNC(it, j).separation_vertex = UNDEFINED_CUBE_VERTEX;
        }
      }
    }
  }


  // **************************************************
  // DETERMINE SEPARATION EDGES
  // **************************************************

  /// Return true if vertex with isosurface connection 
  ///   directions iso_connect_dir separates a cube edge.
  /// @param[out] separation_vertex Separated cube vertex.
  template <typename DTYPE, typename DIR_BITS_TYPE, 
            typename CUBE_TYPE, typename ETYPE>
  bool separates_cube_edge
  (const DTYPE dimension, const DIR_BITS_TYPE iso_connect_dir, 
   const CUBE_TYPE & cube, ETYPE & separation_edge)
  {
    // Initialize
    separation_edge = 0;

    ETYPE num_iso_connect_dir = 0;
    ETYPE num_iso_connect_in_both_dir = 0;
    ETYPE iend0 = 0;
    DTYPE edge_dir = 0;

    for (DTYPE d = 0; d < dimension; d++) {
      const DIR_BITS_TYPE mask = (DIR_BITS_TYPE(1) << d);
      const DIR_BITS_TYPE mask2 = (DIR_BITS_TYPE(1) << (d+dimension));

      if ((mask & iso_connect_dir) == 0) {
        if ((mask2 & iso_connect_dir) != 0) {
          iend0 = (iend0 | mask);
          num_iso_connect_dir++;
        }
      }
      else if ((mask2 & iso_connect_dir) == 0) {
        num_iso_connect_dir++;
      }
      else {
        // Connects in both positive and negative directions.
        num_iso_connect_in_both_dir++;
        num_iso_connect_dir += 2;
        edge_dir = d;
      }
    }

    if (num_iso_connect_in_both_dir != 1) { return(false); }

    if (num_iso_connect_dir != dimension+1) { return(false); }

    separation_edge = cube.IncidentEdge(iend0, edge_dir);

    return(true);
  }


  /// Determine separation cube edges for ivol vertex whose incident
  ///   isosurface polytopes separate a cube edge.
  /// @pre Call compute_ivoldual_table_num_incident_poly
  ///   and compute_ivoldual_table_num_incident_isosurface_poly
  ///   before calling this routine.
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void determine_separation_edges
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE
      VINFO_TYPE;
    typedef typename VINFO_TYPE::CUBE_EDGE_TYPE CUBE_EDGE_TYPE;
    typedef typename VINFO_TYPE::DIR_BITS_TYPE DIR_BITS_TYPE;

    const DTYPE dimension = table.Dimension();
    const DTYPE DIM2(2);
    const CUBE_EDGE_TYPE UNDEFINED_CUBE_EDGE =
      DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE::UNDEFINED_CUBE_EDGE;

    // Initialize
    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {
      const NTYPE num_ivolv = table.Entry(it).NumVertices();

      for (NTYPE j = 0; j < num_ivolv; j++) {
        const DIR_BITS_TYPE iso_connect_dir =
          vinfo.VertexInfo(it, j).iso_connect_dir;
        CUBE_EDGE_TYPE separation_edge;

        if (dimension > DIM2 &&
            separates_cube_edge
            (dimension, iso_connect_dir, table.Cube(), separation_edge)) {
          vinfo.VertexInfoNC(it, j).separation_edge = separation_edge;
        }
        else {
          vinfo.VertexInfoNC(it, j).separation_edge = UNDEFINED_CUBE_EDGE;
        }

      }
    }
  }

}

#endif


