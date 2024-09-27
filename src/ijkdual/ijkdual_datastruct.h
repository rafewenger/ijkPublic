/*!
 *  \file ijkdual_datastruct.h
 *  @brief Data structure definitions for ijkdual.
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

#ifndef _IJKDUAL_DATASTRUCT_H_
#define _IJKDUAL_DATASTRUCT_H_

#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkisopoly.tpp"
#include "ijkmerge.tpp"
#include "ijkmesh.tpp"
#include "ijkscalar_grid.tpp"
#include "ijktri_interior_vertex.tpp"


#include "ijkdual_isosurface.tpp"
#include "ijkdualtable.tpp"
#include "ijkdualtable_vert.tpp"

#include "ijkdual_constants.h"
#include "ijkdual_types.h"


namespace IJKDUAL {

  class MERGE_DATA;
  class MULTIRES_GRID;


  // ***************************************************************
  // GRID DATA STRUCTURES
  // ***************************************************************

  typedef IJK::GRID_PLUS<int, AXIS_SIZE_TYPE, VERTEX_INDEX, VERTEX_INDEX> 
    GRID_PLUS;                     ///< Regular grid.
  typedef IJK::GRID_SPACING<COORD_TYPE, GRID_PLUS>
    DUALISO_GRID;                  ///< Grid and spacing information.
  typedef IJK::SCALAR_GRID_BASE<DUALISO_GRID, SCALAR_TYPE> 
    DUALISO_SCALAR_GRID_BASE;      ///< Marching Cubes base scalar grid.
  typedef IJK::SCALAR_GRID_WRAPPER<DUALISO_GRID, SCALAR_TYPE>
    DUALISO_SCALAR_GRID_WRAPPER;   ///< Marching Cubes scalar grid wrapper.
  typedef IJK::SCALAR_GRID<DUALISO_GRID, SCALAR_TYPE> 
    DUALISO_SCALAR_GRID;           ///< Marching Cubes scalar grid.


  // ***************************************************************
  // ISOSURFACE LOOKUP TABLES
  // ***************************************************************

  typedef IJKDUALTABLE::ISODUAL_TABLE_ENTRY<int,int,unsigned char> ISODUAL_TABLE_ENTRY;

  typedef IJKDUALTABLE::ISODUAL_AMBIG_TABLE_ENTRY<int,int,unsigned char,FACET_BITS_TYPE> 
  ISODUAL_AMBIG_TABLE_ENTRY;

  typedef IJKDUALTABLE::ISODUAL_CUBE_TABLE
  <int,int,TABLE_INDEX,ISODUAL_TABLE_ENTRY>
  ISODUAL_CUBE_TABLE;

  typedef IJKDUALTABLE::ISODUAL_CUBE_TABLE_AMBIG
  <int,int,TABLE_INDEX,ISODUAL_AMBIG_TABLE_ENTRY>
  ISODUAL_CUBE_TABLE_AMBIG;


  // ***************************************************************
  // ISOSURFACE LOOKUP TABLE VERTEX CONNECTIVITY INFORMATION
  // ***************************************************************

  class DUAL_TABLE_VERTEX_CONNECTIVITY_ENTRY:
    public IJKDUALTABLE::
    VERTEX_CONNECTIVITY_INFO_CUBE<INTERSECTS_FACET_BIT_SET_SIZE>
  {
  };

  typedef typename IJKDUALTABLE::DUAL_TABLE_VERTEX_INFO
    <int, DUAL_TABLE_VERTEX_CONNECTIVITY_ENTRY>
    DUAL_TABLE_VERTEX_CONNECTIVITY;


  // ***************************************************************
  // CUBE DATA STRUCTURES
  // ***************************************************************

  typedef IJK::CUBE_FACE_INFO<int, int, int>
    DUALISO_CUBE_FACE_INFO;


  // ***************************************************************
  // DUAL ISOSURFACE VERTICES
  // ***************************************************************

  typedef IJK::DUAL_ISOVERT
  <ISO_VERTEX_INDEX, FACET_VERTEX_INDEX, TABLE_INDEX,
   ISO_VERTEX_INDEX>
  DUAL_ISOVERT;

  typedef IJK::DUAL_ISOVERT_BFLAG
  <ISO_VERTEX_INDEX, FACET_VERTEX_INDEX, TABLE_INDEX,
   ISO_VERTEX_INDEX>
  DUAL_ISOVERT_BFLAG;


  // ***************************************************************
  // GRID CUBE ISOVERT
  // ***************************************************************

  typedef IJK::GRID_CUBE_ISOVERT
  <ISO_VERTEX_INDEX, TABLE_INDEX, FACET_VERTEX_INDEX,
   ISO_VERTEX_INDEX> GRID_CUBE_ISOVERT;


  // ***************************************************************
  // GRID EDGE
  // ***************************************************************

  typedef IJK::GRID_EDGE<VERTEX_INDEX, DIRECTION_TYPE> GRID_EDGE_TYPE;

  /// Array of grid edges.
  typedef std::vector<GRID_EDGE_TYPE> GRID_EDGE_ARRAY;


  // ***************************************************************
  // ISOSURFACE POLYTOPE INFORMATION
  // ***************************************************************

  class ISOPOLY_INFO_TYPE:
    public IJK::ISODUAL_POLY_INFO<GRID_EDGE_TYPE>
  {

  public:

    /*!
     *  @brief Coefficient indicating position of potential isosurface vertex on grid edge.
     *  - Location = (coef)*(endpoint0 coord) + (1-coef)*(endpoint1 coord).
     */
    COORD_TYPE grid_edge_isov_position_coef;


    /*!
     *  @brief Indicates whether diagonal is inside envelope.
     *  - is_diagonal_in_envelope[0]: If true, then 
     *    quadrilateral diagonal (0,2) is in envelope.
     *  - is_diagonal_in_envelope[1]: If true, 
     *    then quadrilateral diagonal (1,3) is in envelope.
     *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
     *    for envelope definition and references.
     */
    bool is_diagonal_in_envelope[2];


  protected:

    /// Initialize.
    void Init();


  public:

    /// @brief Constructor.
    ISOPOLY_INFO_TYPE() { Init(); }

    /// @brief Constructor.
    template <typename VTYPE, typename DTYPE>
    ISOPOLY_INFO_TYPE
    (const VTYPE iend0, const VTYPE iend1, const DTYPE dir):
      IJK::ISODUAL_POLY_INFO<GRID_EDGE_TYPE>(iend0,iend1,dir)
    { Init(); }

    /// Return coefficient indicating position of potential isosurface vertex on grid edge.
    COORD_TYPE GridEdgeIsovPositionCoef() const
    { return(grid_edge_isov_position_coef); }

    /*!
     *  @brief Return true if quadrilateral diagonal (i,i+2) is in envelope.
     *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
     *    for envelope definition and references.
     *  @param i Indicates diagonal (0,2) or (1,3).
     *  @pre i is 0 or 1.
     */
    bool IsDiagonalInEnvelope(const int i) const
    { return(is_diagonal_in_envelope[i]); }

  };


  // ***************************************************************
  // GRID CUBE DATA
  // ***************************************************************

  typedef IJK::GRID_CUBE_ISOVERT<ISO_VERTEX_INDEX, TABLE_INDEX, 
                                 FACET_VERTEX_INDEX, ISO_VERTEX_INDEX>
  GRID_CUBE_DATA;

  typedef IJK::GRID_CUBE_ISOVERT_AND_CUBE_COORD
  <3, ISO_VERTEX_INDEX, TABLE_INDEX, FACET_VERTEX_INDEX, ISO_VERTEX_INDEX,
   GRID_COORD_TYPE>
  GRID_CUBE_DATA_PLUS;


  // ***************************************************************
  // CLASS INTERIOR_VERTEX_PARAM
  // ***************************************************************

  typedef IJK::INTERIOR_VERTEX_PARAM_T<COORD_TYPE>
  INTERIOR_VERTEX_PARAM;
  
  
  // ***************************************************************
  // DUAL CONTOURING ISOSURFACE CLASS
  // ***************************************************************

  typedef DUAL_ISOSURFACE_BASE
  <int,ISO_VERTEX_INDEX,GRID_EDGE_TYPE,COORD_TYPE,int>
  DUAL_ISOSURFACE;

  typedef DUAL_TRIANGULATED_ISOSURFACE_BASE
  <int,ISO_VERTEX_INDEX,GRID_EDGE_TYPE,COORD_TYPE,int>  
  DUAL_TRIANGULATED_ISOSURFACE_E;

  typedef DUAL_TRIANGULATED_ISOSURFACE_BASE
  <int,ISO_VERTEX_INDEX,ISOPOLY_INFO_TYPE,COORD_TYPE,int>  
  DUAL_TRIANGULATED_ISOSURFACE;
  

  // ***************************************************************
  // RANDOM POS PARAMETERS
  // ***************************************************************

  /*!
   *  @brief Parameters for random positioning of isosurface vertices.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  class RANDOM_POS_PARAM {

  protected:

    ///  Initialize.
    void Init();

  public:
    /// Constuctor.
    RANDOM_POS_PARAM() { Init(); }

    /// Random distribution type.
    RANDOM_DISTRIBUTION distribution;

    /// Random seed.
    RANDOM_SEED_TYPE random_seed;

    /*!
     *  @brief If true, apply offset when generating random points.
     *  - Apply offset to separate points from cube facetsw.
     *  - Apply offset to separate points from cube centers.
     */
    bool flag_gen_isov_apply_offset;

    /*!
     *  @brief If true, generate random points separated by cube centers.
     *  - If more than one isosurface vertex is in a cube, generate the points
     *    separated by cube centers.
     *  - Set to true to generate good distribution of points separated
     *    by cube centers.
     *  - Set to false to test reposition_all_dual_isovert_away_from_cube_facets
     *    and reposition_multi_dual_isovertex_cube.
     */
    bool flag_gen_isov_separated_by_cube_center;

    /*!
     *  @brief Offset for separating isosurface vertices.
     *  - Used in generating isosurface vertex positions away from cube facets
     *    or separating isosurface vertices within cubes.
     */
    COORD_TYPE position_offset;

    /*
     *  @brief If true, generate some random positions on cube boundary.
     *  - Generate coordinates in range [-on_boundary_offset,1+on_boundary_offset]
     *    and then clamp to [0,1].
     */
    bool flag_gen_isov_on_cube_boundary;

    /*
     *  @brief Boundary "width" used for generating some random positions
     *    on cube boundary.
     *  - Generate coordinates in range [-boundary_width,1+boundary_width]
     *    and then clamp to [0,1].
     */
    COORD_TYPE boundary_width;

    /// Default boundary width.
    COORD_TYPE DefaultBoundaryWidth() const
    { return 0.05; }

  };


  // ***************************************************************
  // ISOSURFACE VERTEX POSITION PARAMETERS
  // ***************************************************************

  /*!
   *  @brief Parameters for isosurface vertex positions.
   */
  class ISO_VERTEX_POSITION_PARAM {

  public:

    /// @brief Method to position isosurface vertices.
    VERTEX_POSITION_METHOD vertex_position_method;

    /*!
     *  @brief If true, move all isosurface vertices away from cube facets.
     *  - Move avoids surface self intersections near cube facets.
     *  - Self intersections can be created by degeneracies where
     *    scalar values equal the isovalue.
     */
    bool flag_move_all_isov_away_from_cube_facets;
    
    /*!
     *  @brief If true, move all isosurface vertices away from cube ridges.
     *  - Move avoids surface self intersections near cube ridges.
     *  - Self intersections can be created by degeneracies where
     *    scalar values equal the isovalue.
     *  - In 3D, cube ridges are cube edges.
     */
    bool flag_move_all_isov_away_from_cube_ridges;

    /*!
     *  @brief If true, separate isosurface vertices that
     *    are near cube facets.
     *  - Separate isosurface vertices w0, w1 that are from two different
     *    intersecting grid cubes from intersecting parallel facets.
     */
    bool flag_separate_isov_near_cube_facets;

    /*!
     *  @brief Offset for separating isosurface vertices.
     *  - Used in moving isosurface vertices away from cube facets
     *    or separating isosurface vertices within cubes.
     */
    COORD_TYPE position_offset;

    /*!
     *  @brief Min distance to cube faces.
     *  - Vertices within \a min_distance_to_cube_face to cube faces
     *    are repositioned (when flag_move_all_isov_away_from_cube_facets is true.)
     */
    COORD_TYPE min_distance_to_cube_face;

    /// @brief Minimum number of bins for routine 
    ///   IJKSORT::merge_identical_radix_sort_iter2().
    int min_num_radix_sort_bins;

    /// @brief Separate multiple iso vertices in a cube by cube center.
    bool flag_separate_isov_by_cube_center;

    /// @brief Parameters for assigning random positions.
    RANDOM_POS_PARAM random_pos_param;

    
    // *** Default routines ***
    
    /// @brief Default value of position offset.
    COORD_TYPE DefaultPositionOffset() const
    { return 0.01; }


    /// @brief Default value of min_distance_to_cube_face.
    COORD_TYPE DefaultMinDistanceToCubeFace() const
    { return 0.01; }

    /// @brief Default value of min_num_radix_sort_bins.
    COORD_TYPE DefaultMinNumRadixSortBins() const
    { return 100; }    
    
    /// @brief Default random pos boundary width.
    COORD_TYPE DefaultRandomPosBoundaryWidth() const
    { return random_pos_param.DefaultBoundaryWidth(); }

  };


  // ***************************************************************
  // DUAL CONTOURING INPUT DATA AND DATA STRUCTURES
  // ***************************************************************

  /// Dual contouring parameters
  template <typename NTYPE>
  class DUALISO_DATA_PARAM_BASE {

  protected:
    void Init();


  public:

    /// Number type.
    typedef NTYPE NUMBER_TYPE;
    
    // control parameters
    MESH_TYPE mesh_type;
    INTERPOLATION_TYPE interpolation_type;
    QUAD_TRI_METHOD quad_tri_method;
    ENVELOPE_QUAD_TRI_METHOD envelope_quad_tri_method;

    /// If true, constructing an interval volume.
    bool flag_interval_volume;

    /*!
     *  @brief If true, check whether diagonals are in envelope.
     *  - Used to avoid self-intersections in triangulations.
     */
    bool flag_check_envelope;

    /// @brief If true, allow triangulation of quadrilaterals
    ///   into 4 triangles by adding an extra vertex.
    bool flag_tri4_quad;

    /*!
     *  @brief If true and quadrilateral has only one diagonal 
     *    inside the envelope, may triangulate into 2 triangles.
     *  - If true and one diagonal is inside and one diagonal 
     *    is outside the envelope, then consider triangulation 
     *    with diagonal inside the envelope.
     */
    bool flag_allow_tri2_envelope;

    /*!
     *  @brief If true and quadrilateral has at least one diagonal 
     *    outside the envelope, may triangulate into 4 triangles.
     *  - If true and both diagonals are outside the envelope,
     *    then triangulate into 4 triangles.
     *  - If true and one diagonal is inside and one diagonal is outside the envelope,
     *    then consider triangulation into 4 triangles.
     *    (Triangulation depends on flag_allow_tri2_envelope and envelope_quad_tri_method.)
     */
    bool flag_allow_tri4_envelope;

    /// @brief If true, allow multiple isosurface vertices in a single cube.
    /// - Each isosurface vertex in a cube is on a different isosurface patch.
    bool allow_multiple_iso_vertices;

    /*!
     *  @brief If true, split isosurface vertices that
     *    cause non-manifold conditions.
     *  - Isosurface will be a manifold if isovalue does not
     *    equal any grid scalar value.
     */
    bool flag_split_non_manifold;

    /// @brief If true, split isosurface vertices in cubes
    ///   whose ambiguous facets cause non-manifold conditions.
    /// - Used mainly in testing different split options.
    bool flag_split_non_manifold_ambig;

    /// @brief If true, select the cube with split isosurface patch,
    /// when adjacent cubes share an ambiguous facet and one will have
    /// one isosurface vertex while the other has two isosurface vertices.
    bool flag_select_split;

    /// @brief If true, select configurations of ambiguous cubes
    ///   to increase connectivity of the isosurface.
    bool flag_connect_ambiguous;

    /// If true, isosurfaced patches separate negative vertices.
    bool flag_separate_neg;

    /*!
     *  @brief If true, all isosurface quadrilaterals are dual to grid edges.
     *  - Currently always true, but may be useful in future implementations.
     */
    bool flag_iso_quad_dual_to_grid_edges;

    /// If true, call dual_collapse to merge close isosurface vertices.
    bool flag_dual_collapse;

    /// Parameters for positioning isosurface vertices.
    ISO_VERTEX_POSITION_PARAM iso_vertex_position_param;

    /*!
     *  @brief Parameters for determining coordinates of vertex in interior
     *    of isosurface polytopes.
     *  - Used in adding vertices to triangulate isosurface polytopes.
     */
    INTERIOR_VERTEX_PARAM interior_vertex_param;
    
    /*!
     *  @brief Maximum small magnitude.
     *  - Vectors with magnitude less than or equal to max_small_magnitude
     *    are treated as zero vectors.
     */
    COORD_TYPE max_small_magnitude;

    /*!
     *  @brief Envelope offset.
     *  - Shrink envelope by envelope_offset when determining if
     *    a diagonal is within the envelope.
     *  - Avoids numerical errors that may create intersecting edges.
     *  - Also creates greater separation between edges.
     *  - Must be non-negative.
     */
    COORD_TYPE envelope_offset;

    /*!
     *  @brief Use relative vertex positions to help determine
     *  whether quad diagonal is in envelope.
     *  - Allows diagonals in degenerate case when quad vertices
     *    are co-planar and lie on cube facets.
     */
    bool envelope_flag_use_relative_positions;

    /*!
     *  @brief If true, add all interior vertices before any triangulation.
     *  - Some envelope triangulation routines add interior vertices 
     *    individually, and only when necessary.
     *  - Mainly for testing routines.
     */
    bool flag_force_add_all_interior_vertices;
    
    /// @brief Minimum number of bins used for radix sorting.
    NUMBER_TYPE min_num_radix_sort_bins;
     
    
  public:
    DUALISO_DATA_PARAM_BASE() { Init(); };
    ~DUALISO_DATA_PARAM_BASE() { Init(); };

    // Get functions.

    /// @brief Return flag_separate_isov_near_cube_facets.
    /// - flag_offset_isov_near_cube_facets is in iso_vertex_position_param.
    bool FlagSeparateIsovNearCubeFacets() const
    { return iso_vertex_position_param.flag_separate_isov_near_cube_facets; }

    /// @brief Return flag_move_all_isov_away_from_cube_facets.
    /// - flag_move_all_isov_away_from_cube_facets is 
    ///   in iso_vertex_position_param.
    bool FlagMoveAllIsovAwayFromCubeFacets() const
    { return iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets; }

    /// @brief Return flag_move_all_isov_away_from_cube_ridges.
    /// - flag_move_all_isov_away_from_cube_ridges is 
    ///   in iso_vertex_position_param.
    bool FlagMoveAllIsovAwayFromCubeRidges() const
    { return iso_vertex_position_param.flag_move_all_isov_away_from_cube_ridges; }

    /// @brief Return flag_separate_isov_by_cube_center.
    /// - flag_separate_isov_by_cube_center is
    ///   in isov_vertex_position_param.
    bool FlagSeparateIsovByCubeCenter() const
    { return iso_vertex_position_param.flag_separate_isov_by_cube_center; }

    
    // Set functions.

    /// Set current param from data_flags.
    void Set(const DUALISO_DATA_PARAM_BASE<NUMBER_TYPE> & data_flags);

    /// Set vertex position method.
    void SetVertexPositionMethod
    (const VERTEX_POSITION_METHOD vpos_method)
    {
      iso_vertex_position_param.vertex_position_method = vpos_method;
    }

    /*!
     *  @brief Set flag_move_all_isov_away_from_cube_facets.
     *  @details
     *  Note: Setting this to true automatically implies that
     *    all isosurface vertices are moved away from cube ridges.
     */
    void SetFlagMoveAllIsovAwayFromCubeFacets(const bool flag)
    {
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets = flag;
    }

    /*!
     *  @brief Set flag_move_all_isov_away_from_cube_ridges.
     *  @details
     *  Note: If flag_move_all_isov_away_from_cube_facets is true,
     *    then all isosurface vertices are moved away from cube ridges,
     *    and this flag is ignored.
     */
    void SetFlagMoveAllIsovAwayFromCubeRidges(const bool flag)
    {
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_ridges = flag;
    }        

    /*!
     *  @brief Set flag offset to interior.
     *  @details
     *  Note: Setting this to true automatically implies that
     *    all isosurface vertices are moved away from cube ridges.
     */
    void SetFlagMoveIsovAwayFromCubeFacets(const bool flag)
    {
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets = flag;
    }    
    
    /// @brief Set flag_separate_isov_near_cube_facets.
    void SetFlagSeparateIsovNearCubeFacets(const bool flag)
    {
      iso_vertex_position_param.flag_separate_isov_near_cube_facets = flag;
    }

    /*!
     *  @brief Set min distance to cube face.
     *  - Does not change iso_vertex_position.flag_move_all_isov_away_from_cube_facets.
     *  - Call SetFlagMoveIsovAwayFromCubeFacets() to apply position offset.
     */
    template <typename MIN_DISTANCE_TYPE>
    void SetMinDistanceToCubeFace
    (const MIN_DISTANCE_TYPE min_distance)
    {
      iso_vertex_position_param.min_distance_to_cube_face =
        min_distance;
    }

    /*!
     *  @brief Set position offset.
     *  - Sets position offset both in ISO_VERTEX_POSITION_PARAM
     *    and in ISO_VERTEX_POSITION_PARAM::RANDOM_POS_PARAN,
     *  - Does not change iso_vertex_position.flag_apply_position_offset.
     *  - Call SetFlagApplyPositionOffset(true) to apply position offset.
     */
    template <typename OFFSET_TYPE>
    void SetPositionOffset(const OFFSET_TYPE offset)
    {
      iso_vertex_position_param.position_offset = offset;
      iso_vertex_position_param.random_pos_param.position_offset = offset;
    }

    /// Set flag separate isov by cube center.
    void SetFlagSeparateIsovByCubeCenter(const bool flag)
    { iso_vertex_position_param.flag_separate_isov_by_cube_center = flag; }

    /*!
     *  @brief Set random seed.
     *  - Does not change iso_vertex_position_param.vertex_position_param.
     *  - Call SetVertexPositionMethod(RANDOM_POS) to generate random positions.
     */
    template <typename SEED_TYPE>
    void SetRandomSeed(const SEED_TYPE seed)
    {
      iso_vertex_position_param.random_pos_param.random_seed = seed;
    }

    /*!
     *  @brief Set random position distribution.
     *  - Does not change iso_vertex_position_param.vertex_position_param.
     *  - Call SetVertexPositionMethod(RANDOM_POS) to generate random positions.
     */
    void SetRandomPositionDistribution(const RANDOM_DISTRIBUTION distribution)
    {
      iso_vertex_position_param.random_pos_param.distribution = distribution;
    }

    /*!
     *  @brief Set flag determining if offset is applied in generating random points.
     */
    void SetRPosFlagGenIsovApplyOffset(const bool flag)
    {
      iso_vertex_position_param.random_pos_param.flag_gen_isov_apply_offset = flag;
    }

    /*!
     *  @brief Set flag determining if generated random points 
     *    are separated by cube centers.
     */
    void SetRPosFlagGenIsovSeparatedByCubeCenter(const bool flag)
    {
      iso_vertex_position_param.random_pos_param.flag_gen_isov_separated_by_cube_center = flag;
    }

    /*!
     *  @brief Set flag determining if generate random points are on cube boundary.
     *  - Sets random_pos_param.flag_gen_isov_apply_offset to false.
     */
    void SetRPosFlagGenIsovOnCubeBoundary(const bool flag)
    {
      SetRPosFlagGenIsovApplyOffset(false);
      iso_vertex_position_param.random_pos_param.flag_gen_isov_on_cube_boundary = flag;
    }

    /*!
     *  @brief Set random_pos_param.boundary_width.
     */
    template <typename WTYPE>
    void SetRPosBoundaryWidth(const WTYPE w)
    {
      iso_vertex_position_param.random_pos_param.boundary_width = w;
    }
    

    // Get functions
    bool AllowMultipleIsoVertices() const
      { return(allow_multiple_iso_vertices); }
    bool CollapseFlag() const
      { return flag_dual_collapse; }
    bool SplitNonManifoldFlag() const
      { return flag_split_non_manifold; }
    bool SplitNonManifoldAmbigFlag() const
      { return flag_split_non_manifold_ambig; }
    bool SelectSplitFlag() const
      { return flag_select_split; }
    bool ConnectAmbiguousFlag() const
      { return flag_connect_ambiguous; }
    bool SeparateNegFlag() const
      { return flag_separate_neg; }
    COORD_TYPE MaxSmallMagnitude() const
      { return max_small_magnitude; }
    QUAD_TRI_METHOD QuadTriangulationMethod() const
      { return quad_tri_method; }

    /// Return interpolation type.
    INTERPOLATION_TYPE InterpolationType() const
      { return(interpolation_type); };

    /// Return isosurface vertex position method.
    VERTEX_POSITION_METHOD VertexPositionMethod() const
    { return iso_vertex_position_param.vertex_position_method; }

    /// Return min_distance_to_cube_face.
    COORD_TYPE MinDistanceToCubeFace() const
    { return iso_vertex_position_param.min_distance_to_cube_face; }

    /// Return default min_distance_to_cube_face
    COORD_TYPE DefaultMinDistanceToCubeFace() const
    { return iso_vertex_position_param.DefaultMinDistanceToCubeFace(); }
    
    /// Return position_offset.
    COORD_TYPE PositionOffset() const
    { return iso_vertex_position_param.position_offset; }

    /// Return default position_offset.
    COORD_TYPE DefaultPositionOffset() const
    { return iso_vertex_position_param.DefaultPositionOffset(); }

    /// Return default envelope_offset.
    COORD_TYPE DefaultEnvelopeOffset() const
    { return 0.01; }

    /// Return default interior vertex ratio to cube facet.
    COORD_TYPE DefaultInteriorVertexRatioToCubeFacet() const
    { return interior_vertex_param.DefaultInteriorVertexRatioToCubeFacet(); }

    /// Return random position seed.
    RANDOM_SEED_TYPE RandomPositionSeed() const
    { return iso_vertex_position_param.random_pos_param.random_seed; }

    /// Return random position distribution.
    RANDOM_DISTRIBUTION RandomPositionDistribution() const
    { return iso_vertex_position_param.random_pos_param.distribution; }

    /// Default random pos boundary width.
    COORD_TYPE DefaultRandomPosBoundaryWidth() const
    { return iso_vertex_position_param.DefaultRandomPosBoundaryWidth(); }
    
    /// Return default minimum number of radix sort bins.
    NUMBER_TYPE DefaultMinNumRadixSortBins() const
    { return 100; }
    
    /// Return false if flag_iso_quad_dual_to_grid_edges is false.
    bool CheckIsoQuadDualToGridEdges
    (const char * proc_name, IJK::ERROR & error) const;
  };


  typedef DUALISO_DATA_PARAM_BASE<int> DUALISO_DATA_PARAM;

  /// Input data to Dual Contouring and related algorithms
  template <typename _SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  class DUALISO_DATA_BASE:public DATA_FLAGS_TYPE {

  public:
    typedef typename _SCALAR_GRID_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef _SCALAR_GRID_TYPE SCALAR_GRID_TYPE;


  protected:
    SCALAR_GRID_TYPE scalar_grid;

    // flags
    bool is_scalar_grid_set;

    void Init();
    void FreeAll();


  public:
    DUALISO_DATA_BASE() { Init(); };
    ~DUALISO_DATA_BASE() { FreeAll(); };

    // Set functions
    template <typename SCALAR_GRID_BASE_TYPE>
    void CopyScalarGrid             /// Copy scalar_grid to DUALISO_DATA
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2);
    template <typename SCALAR_GRID_BASE_TYPE>
    void SubsampleScalarGrid        /// Subsample scalar_grid.
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
       const int subsample_resolution);
    template <typename SCALAR_GRID_BASE_TYPE>
    void SupersampleScalarGrid      /// Supersample scalar_grid.
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
       const int supersample_resolution);
    void SetInterpolationType       /// Set type of interpolation.
      (const INTERPOLATION_TYPE interpolation_type);
    void SetVertexPositionMethod    /// Set isosurface vertex position method.
      (const VERTEX_POSITION_METHOD vertex_position_method);
    void SetAllowMultipleIsoVertices /// Set flag for multiple iso vertices.
      (const bool flag);
    void SetConvertQuadToTri
      (const QUAD_TRI_METHOD quad_tri_method);
    void UnsetConvertQuadToTri();

    /// Copy, subsample or supersample scalar grid.
    /// Precondition: flag_subsample and flag_supersample are not both true.
    template <typename SCALAR_GRID_BASE_TYPE>
    void SetScalarGrid
      (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
       const bool flag_subsample, const int subsample_resolution,
       const bool flag_supersample, const int supersample_resolution);

    // Get functions
    bool IsScalarGridSet() const     /// Return true if scalar grid is set.
      { return(is_scalar_grid_set); };
    const DUALISO_SCALAR_GRID_BASE & ScalarGrid() const /// Return scalar_grid
      { return(scalar_grid); };

    /// Check data structure.
    /// Return true if no errors found.
    bool Check(IJK::ERROR & error) const;
  };

  typedef DUALISO_DATA_BASE<DUALISO_SCALAR_GRID, DUALISO_DATA_PARAM> 
  DUALISO_DATA;


  // ***************************************************************
  // CLASSES FOR STORING RUNNING TIMES
  // ***************************************************************

  /*!
   *  @brief Time for positioning dual isosurface vertices.
   *  - Uses system clock() function to determine time.
   *  - System clock() function should return cpu time,
   *    but may return wall time.
   */
  class DUALISO_POSITION_TIME {

  public:
    // all times are in seconds

    /// Time for initial/basic positioning of isosurface vertices.
    float basic;
    float separate_by_cube_center;

    /// @brief Time for separating isosurface vertices near cube facets.
    float separate_isov_near_cube_facets;

    /// @brief Time for moving isosurface vertices away from cube faces (facets or ridges.)
    /// @details
    /// - Note: Parameter does not include "all" and has "faces", not "facets".
    float move_isov_away_from_cube_faces;

    /// Time to create (position) isosurface vertices dual to isosurface polytopes.
    float isov_dual_to_isopoly;
    
    /// Total time for positioning isosurface vertices.
    float total;

    DUALISO_POSITION_TIME();
    void Clear();
    void Add(const DUALISO_POSITION_TIME & dualiso_position_time);
  };

  
  /*!
   *  @brief Dual contouring time.
   *  - Uses system clock() function to determine time.
   *  - System clock() function should return cpu time,
   *    but may return wall time.
   */
  class DUALISO_TIME {

  public:
    // all times are in seconds

    /// time to create data structure for faster isosurface extraction
    float preprocessing;  

    float extract;                    ///< time to extract isosurface mesh
    float merge;                      ///< time to merge identical vertices
    float split_isov;                 ///< time to split isosurface vertices
    DUALISO_POSITION_TIME position;   ///< time to position isosurface vertices
    float collapse;                   ///< time to collapse isosurface vertices
    float envelope;                   ///< time to construct grid edge envelopes
    float triangulation;              ///< time to perform triangulations
    float total;                      ///< total time of all operations

    DUALISO_TIME();
    void Clear();
    void Add(const DUALISO_TIME & dualiso_time);
  };


  // ***************************************************************
  // GRID INFO
  // ***************************************************************

  /// Regular grid information.
  class GRID_INFO {

  public:
    GRID_INFO();                    ///< Constructor.

    VERTEX_INDEX num_cubes;         ///< Number of grid cubes.

    void Clear();                   ///< Clear all data.
  };


  // ***************************************************************
  // SCALAR INFO
  // ***************************************************************

  /// Scalar grid information.
  class SCALAR_INFO {

  protected:
    int dimension;

    void Copy(const SCALAR_INFO & info);  ///< Copy scalar info.
    void Init(const int dimension);       ///< Initialize scalar info.
    void FreeAll();                       ///< Free all memory.

  public:
    SCALAR_INFO() { Init(3); };
    SCALAR_INFO(const int dimension) { Init(dimension); };
    ~SCALAR_INFO();
    SCALAR_INFO(const SCALAR_INFO & info) ///< Copy constructor. 
      { Copy(info); };       
    const SCALAR_INFO & operator =        ///< Copy assignment.
      (const SCALAR_INFO & right);

    VERTEX_INDEX num_non_empty_cubes;     ///< Number of cubes containing an isosurface simplex.

    VERTEX_INDEX num_bipolar_edges; ///< Number of bipolar edges.
    //   (number of edges with some vertex value less than the isovalue
    //    and some vertex value greater than or equal to the isovalue)


    // Set functions.
    void SetDimension(const int dimension); ///< Set dimension.

    // Get functions.
    int Dimension() const { return(dimension); };

    void Clear();     // clear all data
  };


  // ***************************************************************
  // ISOSURFACE VERTEX INFO
  // ***************************************************************

  /*!
   *  @brief Isosurface vertex information.
   *  - Information about number of isosurface vertices repositioned
   *    to avoid isosurface self intersections.
   */
  class ISOV_INFO {

  public:
    ISOV_INFO();

    /// @brief Number of times isosurface vertices moved away
    ///   from grid cube facets.
    /// - Note: An isosurface vertex may be moved away
    ///   from multiple grid cube facets.
    ISO_VERTEX_INDEX num_times_isov_moved_away_from_grid_cube_facets;

    void Clear();     //< Clear all data.
  };

  
  // ***************************************************************
  // MULTI ISOV INFO
  // ***************************************************************

  /// Information on mutiple isosurface vertices per cube
  class MULTI_ISOV_INFO {

  public:
    MULTI_ISOV_INFO();

    VERTEX_INDEX num_cubes_single_isov;
    VERTEX_INDEX num_cubes_multi_isov;
    VERTEX_INDEX num_non_manifold_split;
    VERTEX_INDEX num_1_2_changed;
    VERTEX_INDEX num_connect_changed;
    VERTEX_INDEX num_ambig_ridge_cubes_changed;
    VERTEX_INDEX num_non_ambig_ridge_cubes_changed;

    
    void Clear();     //< Clear all data.
  };


  // ***************************************************************
  // TRIANGULATION INFO
  // ***************************************************************

  /*!
   *  @brief Triangulation information.
   */
  class TRIANGULATION_INFO {

  public:

    /// Total number of triangulated isosurface cubes.
    MESH_SIZE_TYPE num_iso_cubes_tri_total;
    
    /// Number of isosurface cubes triangulated without additional vertices.
    MESH_SIZE_TYPE num_iso_cubes_tri_no_add;

    /// Number of isosurface cubes triangulated with one additional interior vertex.
    MESH_SIZE_TYPE num_iso_cubes_tri_add_interior1;
    
    /// Number of isosurface cubes with at least one diagonal outside the envelope.
    MESH_SIZE_TYPE num_iso_cubes_with_diag_outside_envelope;

    /// @brief Number of isosurface cubes split with one interior vertex
    ///   because some diagonal is outside the envelope.
    MESH_SIZE_TYPE num_tri_add_interior1_with_diag_outside_envelope;

    /// @brief Number of isosurface cubes split without additional vertices
    ///   because some diagonal is outside the envelope.
    MESH_SIZE_TYPE num_tri_no_add_with_diag_outside_envelope;

  public:
    TRIANGULATION_INFO();
    
    void Clear();     //< Clear all data.
  };


  // ***************************************************************
  // DUALISO INFO
  // ***************************************************************

  /// @brief Dual contouring information.
  /// - Statistical and timing information from the Dual Contouring algorithm.
  class DUALISO_INFO {

  public:
    GRID_INFO grid;
    SCALAR_INFO scalar;
    DUALISO_TIME time;
    ISOV_INFO isov;
    MULTI_ISOV_INFO multi_isov;
    TRIANGULATION_INFO triangulation;

    DUALISO_INFO();
    DUALISO_INFO(const int dimension);

    void Clear();     // clear all data
  };


  // ***************************************************************
  // MERGE DATA
  // ***************************************************************

  /// Internal data structure for merge_identical_vertices 
  class MERGE_DATA: public IJK::INTEGER_LIST<MERGE_INDEX, MERGE_INDEX> {

  protected:
    MERGE_INDEX num_edges;             ///< Number of edges.
    MERGE_INDEX num_vertices;          ///< Number of vertices.
    MERGE_INDEX num_obj_per_vertex;    ///< Number of objects per vertex.
    MERGE_INDEX num_obj_per_edge;      ///< Number of objects per edge.
    MERGE_INDEX num_obj_per_grid_vertex; ///< Number of objects per grid vertex.
    MERGE_INDEX vertex_id0;            ///< First vertex identifier.

    /// Initialize.
    void Init(const int dimension, const AXIS_SIZE_TYPE * axis_size,
              const MERGE_INDEX num_obj_per_vertex,
              const MERGE_INDEX num_obj_per_edge);

  public:
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size)
      { Init(dimension, axis_size, 0, 1); };
    MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size,
               const MERGE_INDEX num_obj_per_vertex, 
               const MERGE_INDEX num_obj_per_edge)
      { Init(dimension, axis_size, num_obj_per_vertex, num_obj_per_edge); };

    // get functions
    MERGE_INDEX NumEdges() const        /// Number of edges.
      { return(num_edges); };
    MERGE_INDEX NumVertices() const     /// Number of vertices.
      { return(num_vertices); };
    MERGE_INDEX NumObjPerVertex() const { return(num_obj_per_vertex); };
    MERGE_INDEX NumObjPerEdge() const { return(num_obj_per_edge); };
    MERGE_INDEX NumObjPerGridVertex() const
      { return(num_obj_per_grid_vertex); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv) const { return(vertex_id0 + iv); };
    MERGE_INDEX VertexIdentifier       /// Vertex identifier.
      (const MERGE_INDEX iv, const MERGE_INDEX j) const 
      { return(vertex_id0 + iv + j * num_vertices); };
    MERGE_INDEX EdgeIdentifier         /// Edge identifier.
      (const MERGE_INDEX ie) const { return(ie); };
    MERGE_INDEX EdgeIdentifier         /// edge identifier
      (const MERGE_INDEX ie, const MERGE_INDEX j) const 
      { return(ie + j * num_edges); };

    /// Get first endpoint of edge containing isosurface vertex isov.
    inline VERTEX_INDEX GetFirstEndpoint(const MERGE_INDEX isov) const
      { return(isov/NumObjPerGridVertex()); };

    /// Get direction of edge containing isosurface vertex isov.
    inline MERGE_INDEX GetEdgeDir(const MERGE_INDEX isov) const
      { return(isov%NumObjPerGridVertex()); };

    bool Check(IJK::ERROR & error) const;     ///< Check allocated memory.
  };

  /// @brief Merge data structure for isosurface vertices.
  /// - Vertices are identified by a single integer.
  class ISO_MERGE_DATA: public MERGE_DATA {
  public:
    ISO_MERGE_DATA(const int dimension, const AXIS_SIZE_TYPE * axis_size):
      MERGE_DATA(dimension, axis_size, 1, 0) {};
  };


  // *******************************************************************
  // CLASS DUALISO_DATA_PARAM
  // *******************************************************************

  // Initialize DUALISO_DATA_PARAM
  template <typename NTYPE>
  void DUALISO_DATA_PARAM_BASE<NTYPE>::Init()
  {
    flag_interval_volume = false;
    mesh_type = CUBE_COMPLEX;
    interpolation_type = LINEAR_INTERPOLATION;
    allow_multiple_iso_vertices = true;
    flag_split_non_manifold = true;
    flag_separate_neg = true;
    flag_select_split = false;
    flag_connect_ambiguous = false;
    flag_iso_quad_dual_to_grid_edges = true;
    flag_dual_collapse = false;
    max_small_magnitude = 0.0001;
    quad_tri_method = UNDEFINED_TRI;
    envelope_quad_tri_method = ENVELOPE_UNDEFINED_TRI;
    envelope_offset = DefaultEnvelopeOffset();
    envelope_flag_use_relative_positions = false;
    flag_force_add_all_interior_vertices = false;

    interior_vertex_param.interior_vertex_position_method =
      IJK::INTERIOR_VERTEX_IN_ENVELOPE;

    min_num_radix_sort_bins = DefaultMinNumRadixSortBins();

    flag_check_envelope = false;
    flag_allow_tri2_envelope = true;
    flag_allow_tri4_envelope = true;
    flag_tri4_quad = false;

    // Initialize iso_vertex_position_param.
    iso_vertex_position_param.vertex_position_method =
      CENTROID_EDGE_ISO;

    iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets = false;
    iso_vertex_position_param.flag_move_all_isov_away_from_cube_ridges = false;
    iso_vertex_position_param.flag_separate_isov_near_cube_facets = false;
    iso_vertex_position_param.flag_separate_isov_by_cube_center = false;
    iso_vertex_position_param.min_num_radix_sort_bins =
      iso_vertex_position_param.DefaultMinNumRadixSortBins();
    SetPositionOffset(iso_vertex_position_param.DefaultPositionOffset());
    SetMinDistanceToCubeFace
      (iso_vertex_position_param.DefaultMinDistanceToCubeFace());

    // Initialize random position separation param.
    SetRPosFlagGenIsovApplyOffset(true);
    SetRPosFlagGenIsovSeparatedByCubeCenter(true);

    // Option used for testing.
    flag_split_non_manifold_ambig = false;
  }


  // Set
  template <typename NTYPE>
  void DUALISO_DATA_PARAM_BASE<NTYPE>::Set
  (const DUALISO_DATA_PARAM_BASE<NTYPE> & data_param)
  {
    *this = data_param;
  }


  // Return false if flag_iso_quad_dual_to_grid_edges is false.
  template <typename NTYPE>
  bool DUALISO_DATA_PARAM_BASE<NTYPE>::CheckIsoQuadDualToGridEdges
  (const char * proc_name, IJK::ERROR & error) const
  {
    if (flag_iso_quad_dual_to_grid_edges) {
      return (true);
    }
    else {
      error.AddMessage
        ("Programming error.  Isosurface quadrilaterals may not be dual to grid edges.");
      error.AddMessage
        ("  Procedure ", proc_name, " requires ");
      error.AddMessage
        ("  that all isosurface quadrilaterals are dual to grid edges.");
      return(false);
    }
  }


  // ***************************************************************
  // DUALISO_DATA_BASE CLASS MEMBER FUNCTIONS
  // ***************************************************************

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::Init()
  {
    is_scalar_grid_set = false;
  }

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE, DATA_FLAGS_TYPE>::FreeAll()
  {
    is_scalar_grid_set = false;
  }


  // Copy scalar grid
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE, DATA_FLAGS_TYPE>::
  CopyScalarGrid(const SCALAR_GRID_BASE_TYPE & scalar_grid2)
  {
    scalar_grid.Copy(scalar_grid2);
    scalar_grid.SetSpacing(scalar_grid2.SpacingPtrConst());
    is_scalar_grid_set = true;
  }

  // Subsample scalar grid
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SubsampleScalarGrid(const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
                      const int subsample_resolution)
  {
    const int dimension = scalar_grid2.Dimension();
    IJK::ARRAY<COORD_TYPE> spacing(dimension);

    scalar_grid.Subsample(scalar_grid2, subsample_resolution);

    IJK::copy_coord(dimension, scalar_grid2.SpacingPtrConst(),
                    spacing.Ptr());
    IJK::multiply_coord
      (dimension, subsample_resolution, spacing.PtrConst(), spacing.Ptr());
    scalar_grid.SetSpacing(spacing.PtrConst());

    is_scalar_grid_set = true;
  }

  // Supersample scalar grid
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SupersampleScalarGrid
  (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
   const int supersample_resolution)
  {
    const int dimension = scalar_grid2.Dimension();
    IJK::ARRAY<COORD_TYPE> spacing(dimension);
    scalar_grid.Supersample(scalar_grid2, supersample_resolution);

    IJK::copy_coord(dimension, scalar_grid2.SpacingPtrConst(),
                    spacing.Ptr());
    IJK::divide_coord
      (dimension, supersample_resolution, spacing.PtrConst(), spacing.Ptr());
    scalar_grid.SetSpacing(spacing.PtrConst());

    is_scalar_grid_set = true;
  }

  // Copy, subsample or supersample scalar grid.
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  template <typename SCALAR_GRID_BASE_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::SetScalarGrid
  (const SCALAR_GRID_BASE_TYPE & scalar_grid2, 
   const bool flag_subsample, const int subsample_resolution,
   const bool flag_supersample, const int supersample_resolution)
  {
    IJK::PROCEDURE_ERROR error("DUALISO_DATA_BASE::SetScalarGrid");

    if (flag_subsample && flag_supersample) {
      error.AddMessage
        ("Scalar grid cannot both be subsampled and supersampled.");
      throw error;
    }
  
    if (flag_subsample) {
      // subsample grid
      SubsampleScalarGrid(scalar_grid2, subsample_resolution);
    }
    else if (flag_supersample) {
      // supersample grid
      SupersampleScalarGrid(scalar_grid2, supersample_resolution);
    }
    else {
      CopyScalarGrid(scalar_grid2);
    };
  }

  // Set type of interpolation
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetInterpolationType
  (const INTERPOLATION_TYPE interpolation_type)
  {
    this->interpolation_type = interpolation_type;
  }

  // Set isosurface vertex position method.
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetVertexPositionMethod    
  (const VERTEX_POSITION_METHOD vertex_position_method)
  {
    this->vertex_position_method = vertex_position_method;
  }

  // Set flag for allowing multiple isosurface vertices in a single cube.
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetAllowMultipleIsoVertices(const bool flag)
  {
    this->allow_multiple_iso_vertices = flag;
  }

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  SetConvertQuadToTri(const QUAD_TRI_METHOD method)
  {
    this->use_triangle_mesh = true;
    this->quad_tri_method = method;
  }

  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  void DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  UnsetConvertQuadToTri()
  {
    this->use_triangle_mesh = false;
  }


  /// Check data structure
  template <typename SCALAR_GRID_TYPE, typename DATA_FLAGS_TYPE>
  bool DUALISO_DATA_BASE<SCALAR_GRID_TYPE,DATA_FLAGS_TYPE>::
  Check(IJK::ERROR & error) const
  {
    IJK::ERROR error2;

    if (!IsScalarGridSet()) {
      error.AddMessage("Scalar grid is not set.");
      return(false);
    }

    return(true);
  }

}

#endif
