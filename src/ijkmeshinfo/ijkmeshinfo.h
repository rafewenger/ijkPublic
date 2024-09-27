/*!
 *  @file ijkmeshinfo.h
 *  @brief Compute mesh information.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2010-2024 Rephael Wenger

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

/*!
    @mainpage IJKMESHINFO: IJK MESH INFORMATION

    IJKMESHINFO is a program for computing and retrieving information
    about a polytope mesh.
    The mesh may be in embedded in any dimension and
    the mesh dimension may be arbitrary.
    The program supports arbitrary polygonal meshes (triangles, quads, etc.)
    embedded in any dimension, simplices (triangles, tetrahedra, etc.)
    embedded in any dimension, and boxes (quads, hexahedra, etc.)
    embedded in any dimension.
 */

#ifndef _IJKMESHINFO_
#define _IJKMESHINFO_

#include <iostream>
#include <string>
#include  <utility>

#include "ijk.tpp"
#include "ijkcube.tpp"
#include "ijkdatatable.tpp"
#include "ijkhash.tpp"
#include "ijkmesh_datastruct.tpp"
#include "ijkmesh_faces.tpp"


namespace IJKMESHINFO {

  // **************************************************
  // Global constants
  // **************************************************

  const int DIM2(2);
  const int DIM3(3);
  const int NUM_VERT_PER_TETRAHEDRON(4);
  const int NUM_VERT_PER_HEXAHEDRON(8);


  // **************************************************
  // Types
  // **************************************************

  typedef int VERTEX_INDEX;
  typedef int POLY_INDEX;
  typedef double COORD_TYPE;    // Use double to detect self-intersections.
  typedef COORD_TYPE * COORD_TYPE_PTR;
  typedef float ANGLE_TYPE;
  typedef IJK::BOX<COORD_TYPE> BOUNDING_BOX;
  typedef int NUM_TYPE;

  typedef IJK::POLYMESH<VERTEX_INDEX,int> POLYMESH_TYPE;
  typedef typename IJK::VERTEX_POLY_INCIDENCE_WITH_VLOC<int,int,int> 
  VERTEX_POLY_INCIDENCE_TYPE;
  typedef typename IJK::VERTEX_ADJACENCY_LIST<int,int>
  VERTEX_ADJACENCY_LIST_TYPE;
  typedef typename IJK::UNORIENTED_EDGE_HASH_TABLE<int,int>
  EDGE_HASH_TABLE_TYPE;
  typedef typename IJK::CUBE_FACE_INFO<int,int,int> CUBE_TYPE;

  typedef typename IJK::FACE_INFO_BASE<int,int> FACET_INFO;
  typedef std::vector<FACET_INFO> FACET_INFO_ARRAY;
  typedef std::pair<FACET_INFO,FACET_INFO> FACET_INFO_PAIR;

  typedef std::vector<FACET_INFO_PAIR> FACET_INFO_PAIRS_ARRAY;

  typedef std::pair<POLY_INDEX,POLY_INDEX> INTERSECTING_POLY_PAIR;
  typedef std::vector<INTERSECTING_POLY_PAIR> INTERSECTING_POLY_ARRAY;
  
  
  // **************************************************
  // Class POLY_DATA
  // **************************************************

  class POLY_DATA {

  public:
    /// True if polytope is degenerate (e.g. has two identical vertices.)
    bool is_degenerate;

    /// True if polytope is duplicate of some other polytope.
    bool is_duplicate;

    /// True if polytope contains a non-manifold edge.
    bool contains_non_manifold_edge;

    /// True if polytope contains a non-manifold facet.
    bool contains_non_manifold_facet;

    /// True if polytope contains a boundary facet.
    bool contains_boundary_facet;

    /// True if polytope shares a facet with another polytope
    ///   with different orientation.
    bool orientation_conflict;

  protected:
    void Init();

  public:

    POLY_DATA() { Init(); }

    // Get functions

    /// Return true if polytope is a duplicate of some other polytope.
    bool IsDuplicate() const
    { return(is_duplicate); }

    /// Return true if polytope is degenerate.
    bool IsDegenerate() const
    { return(is_degenerate); }

    /// Return true if polytope contains a non-manifold edge.
    bool ContainsNonManifoldEdge() const
    { return(contains_non_manifold_edge); }

    /// Return true if polytope contains a non-manifold facet.
    bool ContainsNonManifoldFacet() const
    { return(contains_non_manifold_facet); }

    /// Return true if polytope contains a boundary facet.
    bool ContainsBoundaryFacet() const
    { return(contains_boundary_facet); }

    /// Return true if polytope shares a facet with another polytope
    ///   with different orientation.
    bool OrientationConflict() const
    { return(orientation_conflict); }
  };


  // **************************************************
  // Class VERTEX_DATA
  // **************************************************

  class VERTEX_DATA {

  protected:
    void Init();

  public:

    /// True if vertex is on surface/volume boundary.
    bool on_boundary;

  public:
    VERTEX_DATA() { Init(); }

    /// Return true if vertex is on surface/volume boundary.
    bool OnBoundary() const
    { return(on_boundary); }

    /// Return true if vertex is internal to surface/volume.
    bool IsInternal() const
    { return(!OnBoundary()); }

  };


  // **************************************************
  // Class EDGE_DATA
  // **************************************************

  class EDGE_DATA {

  protected:
    void Init();

  public:

    /// Endpoints
    int endpoint[2];

    /// True if edge is on surface/volume boundary.
    bool on_boundary;

  public:
    EDGE_DATA() { Init(); }

    /// Return true if edge is on surface boundary or
    ///   incident on a volume boundary facet.
    /// - Note: Both vertices of an edge could be on a volume boundary,
    ///   but the edge could still be internal.
    bool OnBoundary() const
    { return(on_boundary); }

    /// Return true if edge is on surface boundary or
    ///   incident on a volume boundary facet.
    bool IsInternal() const
    { return(!OnBoundary()); }
  };


  // **************************************************
  // Class MESH_INFO
  // **************************************************

  /// Mesh information.
  class MESH_INFO {

  protected:
    void Init();

  public:
    MESH_INFO(){ Init(); };

    int num_poly_with_duplicate_vertices;
    int num_duplicate_poly;
    int num_non_manifold_facets;
    int num_non_manifold_edges;
    int num_non_manifold_vertices;
    int num_deep_non_manifold_vertices;
    int num_poly_with_orientation_conflicts;
    int num_hex_facet_pairs_sharing_exactly_two_edges;

    /// Return true if all non-manifold numbers are zero.
    bool AreAllNonManifoldZero() const;

    /// Return true if all numbers are zero.
    bool AreAllZero() const;
  };


  // **************************************************
  // Class MESH_DATA
  // **************************************************

  class MESH_DATA {

  protected:
    COORD_TYPE max_small_magnitude;

  protected:
    void Init();

  public:
    typedef int NUMBER_TYPE;
    typedef int DIMENSION_TYPE;
    typedef POLYMESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    int dimension;
    int mesh_dimension;
    int num_vertices;

    /// Mesh orientation. +1 or -1.
    int orientation;

    /// True if boundary facets have been identified.
    bool are_boundary_facets_identified;

    /// True if nonmanifold facets have been identified.
    bool are_non_manifold_facets_identified;

    /// True if hex sharing exactly two facets edges have been identified.
    bool are_hex_sharing_exactly_two_facet_edges_identified;

    /// Polymesh.
    POLYMESH_TYPE polymesh;

    /// Vertex poly incidence.
    VERTEX_POLY_INCIDENCE_TYPE vertex_poly_incidence;

    /// Poly data.
    std::vector<POLY_DATA> poly_data;

    /// Vertex data.
    std::vector<VERTEX_DATA> vertex_data;

    /// Vertex adjacency list;
    VERTEX_ADJACENCY_LIST_TYPE vertex_adjacency_list;

    /// Edge data.
    std::vector<EDGE_DATA> edge_data;

    /// Edge hash table
    EDGE_HASH_TABLE_TYPE edge_hash_table;


  public:
    MESH_DATA() { Init(); }

    /// Return number of polytopes.
    int NumPoly() const
    { return(polymesh.NumPoly()); }

    /// Return number of vertices.
    int NumVertices() const
    { return(num_vertices); }

    /// Return number of vertices of polytope i.
    int NumPolyVert(const int ipoly) const
    { return(polymesh.NumPolyVert(ipoly)); }

    /// Return pointer to first vertex of polytope ipoly.
    const VERTEX_INDEX * VertexList(const int ipoly) const
    { return(polymesh.VertexList(ipoly)); }

    /// Return jth vertex of polygon ipoly.
    const VERTEX_INDEX Vertex(const int ipoly, const int j) const
    { return(polymesh.Vertex(ipoly, j)); }

    /// Return the value used for max_small_magnitude.
    COORD_TYPE MaxSmallMagnitude() const
    { return max_small_magnitude; }

    /// Return true if polytope ipoly contains vertex iv.
    bool DoesPolyContainVertex(const int ipoly, const VERTEX_INDEX iv) const
    { return(polymesh.DoesPolyContainVertex(ipoly, iv)); }

    /// @brief Return true if polytope ipoly contains vertex iv.
    /// - Version which returns location of vertex iv in list.
    bool DoesPolyContainVertex
    (const int ipoly, const VERTEX_INDEX iv, int & iloc) const
    { return(polymesh.DoesPolyContainVertex(ipoly, iv, iloc)); }

    /// @brief Set max_small_magnitude.
    /// - max_small_magnitude is used in normalizing vectors
    ///   and calculating cosine of angles.
    /// - Vector lengths or distances less than or equal 
    ///   to max_small_magnitude are considered insignificant
    ///   and are set to 0.
    /// @param magnitude Set max_small_magnitude to magnitude.
    /// @pre magnitude >= 0.
    void SetMaxSmallMagnitude(const COORD_TYPE magnitude);

    /// @brief Add edge (iv0,iv1) to edge_data and edge_hash_table.
    /// - Return edge index.
    /// @pre iv0 < iv1.
    int AddEdge(const VERTEX_INDEX iv0, const VERTEX_INDEX iv1);
  };


  // **************************************************
  // Class GRID_OF_BINS_3D
  // *************************************************

  class GRID_OF_BINS_3D {

  protected:
    std::vector<int> ** bin;
    int num_bins_along_axis[DIM3];
    int num_bins;
    COORD_TYPE minC[DIM3];
    COORD_TYPE maxC[DIM3];
    bool bins_are_empty;

    void Init(const int k);

  public:
    GRID_OF_BINS_3D(const int k) { Init(k); };
    ~GRID_OF_BINS_3D();

    int Dimension() const { return(DIM3); };
    int NumBinsAlongAxis(const int d) const 
    { return(num_bins_along_axis[d]); };
    int NumBins() const { return(num_bins); };
    const COORD_TYPE * MinCoord() const
    { return(minC); }
    const COORD_TYPE * MaxCoord() const
    { return(maxC); }

    int LocateBinCoord(const int d, const COORD_TYPE c) const;
    int ComputeBinIndex(const int binCoord[DIM3]) const;

    /// Compute min/max bin coordinates containing triangle (v0,v2,v2).
    void ComputeMinMaxBinCoord
      (const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
       const COORD_TYPE v2[DIM3],
       int min_grid_coord[DIM3], int max_grid_coord[DIM3]) const;

    const std::vector<int> * Bin(const int i)
    { return(bin[i]); };

    void SetMinCoord(const COORD_TYPE minC[DIM3]);
    void SetMaxCoord(const COORD_TYPE maxC[DIM3]);

    void InsertTri(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
                   const COORD_TYPE v2[DIM3], const int triangle_index);

  };


  // **************************************************
  // Class ANGLE_TABLE
  // **************************************************

  template <class NTYPE, class DTYPE>
    class DATA_COLUMN: public IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE> {

    protected:
      bool is_included;
      bool is_hidden;
      DTYPE sum;

      void Init() { is_included = false; is_hidden = false; sum = 0; };

    public:
    DATA_COLUMN<NTYPE,DTYPE>() :
      IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>() { Init(); };
      DATA_COLUMN<NTYPE,DTYPE>
      (const std::string & label, const NTYPE num_rows) :
        IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>(label, num_rows) { Init(); };

    // set functions
    void Include();
    void Hide() { is_hidden = true; };
    void Show() { is_hidden = false; };
    void ComputeSum();

    // get functions
    bool IsIncluded() const { return(is_included); };
    bool IsHidden() const { return(is_hidden); };
    DTYPE Sum() const { return(sum); };

    // write functions
    void WriteLabel(std::ostream & out, const std::string & separator) const;
    void WriteData(std::ostream & out, const std::string & separator, 
		   const NTYPE width, const NTYPE irow) const;
    void WriteNormalizedData
      (std::ostream & out, const std::string & separator, const NTYPE width,
       const double normalization_factor, const NTYPE irow) const;
  };


  class ANGLE_TABLE:public IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE> {

  public:
    static const NUM_TYPE NUM_ROWS = 180;

    DATA_COLUMN<NUM_TYPE, NUM_TYPE> angle;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> min_polygon_angle_freq;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> max_polygon_angle_freq;

  public:
    ANGLE_TABLE():
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(NUM_ROWS),
      angle("angle", NUM_ROWS), 
      min_polygon_angle_freq("min-poly-angle", NUM_ROWS),
      max_polygon_angle_freq("max-poly-angle", NUM_ROWS)
        {};

    // set routines
    void HideAllExceptAngleColumn();
    void ComputeSum();

    // write routines
    void WriteColumnLabels
      (std::ostream & out, const std::string & separator) const;
    void WriteColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width) const;
    void WriteNormalizedColumnData
      (std::ostream & out, const std::string & separator, const NUM_TYPE width,
       const double normalization_factor) const;

  };


  class EDGE_LENGTH_TABLE:
    public IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE> {

  public:
    static const NUM_TYPE NUM_ROWS = 30;

    DATA_COLUMN<NUM_TYPE, COORD_TYPE> edge_length;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> edge_length_freq;

  public:
    EDGE_LENGTH_TABLE():
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(NUM_ROWS),
      edge_length("edge-length", NUM_ROWS), 
      edge_length_freq("edge-length", NUM_ROWS)
        {};

    // set routines
    void HideAllExceptEdgeLengthColumn();
    void ComputeSum();

    // write routines
    void WriteColumnLabels
      (std::ostream & out, const std::string & separator) const;
    void WriteColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width) const;
    void WriteNormalizedColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width,
       const double normalization_factor) const;
  };


  class JACOBIAN_TABLE:
    public IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE> {

  public:
    static const NUM_TYPE NUM_ROWS = 30;

    DATA_COLUMN<NUM_TYPE, COORD_TYPE> jacobian;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> min_jacobian_freq;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> max_jacobian_freq;

  public:
    JACOBIAN_TABLE():
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(NUM_ROWS),
      jacobian("jacobian", NUM_ROWS), 
      min_jacobian_freq("min-jacobian", NUM_ROWS),
      max_jacobian_freq("max-jacobian", NUM_ROWS)
        {};

    JACOBIAN_TABLE(const NUM_TYPE num_rows):
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(num_rows),
      jacobian("jacobian", num_rows), 
      min_jacobian_freq("min-jacobian", num_rows),
      max_jacobian_freq("max-jacobian", num_rows)
        {};

    // set routines
    void HideAllExceptJacobianColumn();
    void ComputeSum();

    // write routines
    void WriteColumnLabels
      (std::ostream & out, const std::string & separator) const;
    void WriteColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width) const;
    void WriteNormalizedColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width,
       const double normalization_factor) const;

  };


  class JACOBIAN_SHAPE_TABLE:
    public IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE> {

  public:
    static const NUM_TYPE NUM_ROWS = 21;

    DATA_COLUMN<NUM_TYPE, COORD_TYPE> jshape;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> min_jshape_freq;
    DATA_COLUMN<NUM_TYPE, NUM_TYPE> max_jshape_freq;

  public:
    JACOBIAN_SHAPE_TABLE():
      IJKDATATABLE::DATA_TABLE_BASE<NUM_TYPE>(NUM_ROWS),
      jshape("jshape", NUM_ROWS), 
      min_jshape_freq("min-jshape", NUM_ROWS),
      max_jshape_freq("max-jshape", NUM_ROWS)
        {};

    // set routines
    void HideAllExceptShapeColumn();
    void ComputeSum();

    // write routines
    void WriteColumnLabels
      (std::ostream & out, const std::string & separator) const;
    void WriteColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width) const;
    void WriteNormalizedColumnData
      (std::ostream & out, const std::string & separator, 
       const NUM_TYPE width,
       const double normalization_factor) const;
  };


  // **************************************************
  // Check functions
  // **************************************************

  /// Return false and set error message if dimension != _DIM_.
  template <const int _DIM_>
  bool check_dimension(const int dimension, IJK::ERROR & error)
  {
    if (dimension == _DIM_)
      { return(true); }
    else {
      error.AddMessage
        ("Programming error.  Dimension should be ", _DIM_, ".");
      error.AddMessage("  Dimension = ", dimension, ".");
      return(false);
    }
  }


  /// Return false and set error message if mesh_data.dimension != _DIM_.
  template <const int _DIM_, typename MESH_DATA_TYPE>
  bool check_dimension
  (const MESH_DATA_TYPE & mesh_data, IJK::ERROR & error)
  {
    const int dimension = mesh_data.dimension;
    return(check_dimension<_DIM_>(mesh_data.dimension, error));
  }


  /// Return false and set error message if mesh_dim != _MESH_DIM_.
  template <const int _MESH_DIM_>
  bool check_mesh_dimension(const int mesh_dim, IJK::ERROR & error)
  {
    if (mesh_dim == _MESH_DIM_)
      { return(true); }
    else {
      error.AddMessage
        ("Programming error.  Mesh dimension should be ", _MESH_DIM_, ".");
      error.AddMessage
        ("  Mesh dimension = ", mesh_dim, ".");
      return(false);
    }
  }


  /// Return false and set error message 
  ///   if mesh_data.mesh_dimension != _MESH_DIM_.
  template <const int _MESH_DIM_, typename MESH_DATA_TYPE>
  bool check_mesh_dimension
  (const MESH_DATA_TYPE & mesh_data, IJK::ERROR & error)
  {
    const int mesh_dimension = mesh_data.mesh_dimension;
    return(check_mesh_dimension<_MESH_DIM_>(mesh_data.mesh_dimension, error));
  }


  /// Return false and set error message if boundary facets are not identified.
  template <typename MESH_DATA_TYPE>
  bool check_boundary_facets
  (const MESH_DATA_TYPE & mesh_data, IJK::ERROR & error)
  {
    if (mesh_data.are_boundary_facets_identified)
      { return(true); }
    else {
      error.AddMessage("Programming error.  Need to compute boundary facets.");
      return(false);
    }
  }


  // **************************************************
  // DATA_COLUMN member functions
  // **************************************************

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::Include()
    {
      is_included = true;

      if (this->NumRows() > 0 && this->data == NULL) 
        { this->data = new DTYPE[this->NumRows()]; }
    }

  template <class NTYPE, class DTYPE>
    void DATA_COLUMN<NTYPE,DTYPE>::ComputeSum()
    {
      if (IsIncluded()) 
        { sum = IJKDATATABLE::DATA_COLUMN<NTYPE,DTYPE>::Sum(); }
      else
        { sum = 0; }
    }

  template <class NTYPE, class DTYPE>
  void DATA_COLUMN<NTYPE,DTYPE>::WriteLabel
  (std::ostream & out, const std::string & separator) const
  {
    if (IsIncluded() && !IsHidden()) 
      { out << separator << this->Label(); }
  }

  template <class NTYPE, class DTYPE>
  void DATA_COLUMN<NTYPE,DTYPE>::WriteData
  (std::ostream & out, const std::string & separator, 
   const NTYPE width, const NTYPE irow) const
  {
    if (IsIncluded() && !IsHidden()) {
      out << separator;
      out.width(width);
      out << this->data[irow]; 
    }
  }

  template <class NTYPE, class DTYPE>
  void DATA_COLUMN<NTYPE,DTYPE>::WriteNormalizedData
  (std::ostream & out, const std::string & separator, const NTYPE width,
   const double normalization_factor, const NTYPE irow) const
  {
    if (IsIncluded() && !IsHidden()) {
      double v = normalization_factor*double(this->data[irow])/Sum();
      out << separator;
      out.width(width);
      out << v; 
    }
  }

};

#endif
