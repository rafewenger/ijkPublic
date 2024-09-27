/*!
 *  @file ijkMCtable_properties.h
 *  @brief Class defining Marching Cubes isosurface table properties.
 *  - Version 0.5.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2024 Rephael Wenger

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

#ifndef _IJKMCTABLE_PROPERTIES_
#define _IJKMCTABLE_PROPERTIES_

#include "ijk.tpp"
#include "ijkenum.tpp"


namespace IJKMCUBE_TABLE {

  // **************************************************
  // ISOSURFACE/INTERVAL VOLUME LOOKUP TABLE PROPERTIES
  // **************************************************

  /*!
   *  @brief Table type.
   *  - ISOSURFACE: (d-1)-dimensional surface in R^d,
   *  - INTERVAL_VOLUME: d-dimensional volume in R^d,
   */
  typedef enum { ISOSURFACE, INTERVAL_VOLUME, UNDEFINED_TABLE_TYPE }
  LOOKUP_TABLE_TYPE;

  namespace {
    static std::vector<IJK::ENUM_STR<LOOKUP_TABLE_TYPE> > 
      _table_type_pairs =
      { {ISOSURFACE,"ISOSURFACE"}, {INTERVAL_VOLUME,"IVOL"},
        {UNDEFINED_TABLE_TYPE,"Unknown"} };
  };

  /*!
   *  @brief Configuration encodings.
   *  - Standard Marching Cubes lookup table with "-/+" vertex labels
   *    uses binary encoding for "-/+".
   *  - Interval Volume lookup table uses base 3 encoding for "-/ * /+"
   *    where '*' represents a scalar value between the two isovalues.
   *  - NEP (negative-equals-positive) isosurface lookup tables
   *    with "-/=/+" labels uses base 3 encoding for "-/=/+".
   *  - UNKOWN_ENCODING is equivalent to an undefined encoding,
   *    since encoding must always be known/defined.
   */
  typedef enum { BINARY, BASE3, NONSTANDARD_ENCODING, UNKNOWN_ENCODING } 
  ENCODING;

  namespace {
    static std::vector<IJK::ENUM_STR<ENCODING> > _encoding_pairs =
      { {BINARY,"BINARY"}, {BASE3,"BASE3"},
        {NONSTANDARD_ENCODING,"NONSTANDARD"},
        {UNKNOWN_ENCODING,"Unknown"} };
  };

  /*!
   *  @brief Types of vertex labels in this table.
   *  - NEG_POS: Grid vertices are labelled either negative 
   *    (below the isovalue) or positive (above the isovalue). 
   *    Grid vertices with values equal to the isovalue are
   *    given a positive label.
   *  - NEG_EQUALS_POS: Grid vertices are labelled either negative
   *    (below the isovalue), or equals (equals the isovalue),
   *    or positive (above the isovalue).
   */
  typedef enum { NEG_POS, NEG_EQUALS_POS, NEG_STAR_POS,
                 UNKNOWN_GRID_VERTEX_LABEL_TYPE,
                 UNDEFINED_GRID_VERTEX_LABEL_TYPE } 
  GRID_VERTEX_LABEL_TYPE;

  namespace {
    static std::vector<IJK::ENUM_STR<GRID_VERTEX_LABEL_TYPE> > 
      _grid_vertex_label_type_pairs =
      { {NEG_POS,"NegPos"}, {NEG_EQUALS_POS,"NegEqualsPos"},
        {NEG_STAR_POS, "NegStarPos"},
        {UNKNOWN_GRID_VERTEX_LABEL_TYPE, "Unknown"},
        {UNDEFINED_GRID_VERTEX_LABEL_TYPE, "Undefined"} };
  };

  /*!
   *  @brief Indicate how isosurface patches are triangulated.
   *  - CONVEX_HULL: Triangulated through convex hull.
   *  - EDGE_GROUPS: Use edge groups to determine triangulation
   *    in appropriate cases (3D only.)
   *  - (See "Edge Groups: An approach to understanding the mesh
   *    quality of Marching Cubes" by Dietrich et al, IEEE TVCG, 2008.
   */
  typedef enum { CONVEX_HULL, EDGE_GROUPS, 
                 UNKNOWN_ISOSURFACE_TRIANGULATION_TYPE,
                 UNDEFINED_ISOSURFACE_TRIANGULATION_TYPE}
  ISOSURFACE_TRIANGULATION_TYPE;

  namespace {
    static std::vector<IJK::ENUM_STR<ISOSURFACE_TRIANGULATION_TYPE> > 
      _isosurface_triangulation_type_pairs =
      { {CONVEX_HULL,"ConvexHull"}, {EDGE_GROUPS,"EdgeGroups"}, 
        {UNKNOWN_ISOSURFACE_TRIANGULATION_TYPE, "Unknown"},
        {UNDEFINED_ISOSURFACE_TRIANGULATION_TYPE, "Undefined"} };
  };

  /// @brief Indicate if isosurface patches separate positive 
  ///   or negative vertices.
  typedef enum { SEPARATE_NEG, SEPARATE_POS, 
                 UNKNOWN_SEPARATION_TYPE,
                 UNDEFINED_SEPARATION_TYPE}
  ISOSURFACE_SEPARATION_TYPE;

  namespace {
    static std::vector<IJK::ENUM_STR<ISOSURFACE_SEPARATION_TYPE> > 
      _isosurface_separation_type_pairs =
      { {SEPARATE_NEG,"SeparateNeg"}, {SEPARATE_POS,"SeparatePos"}, 
        {UNKNOWN_SEPARATION_TYPE, "Unknown"},
        {UNDEFINED_SEPARATION_TYPE, "Undefined"} };
  };

  /*!
   *  @brief If true, always separate diagonally opposite cube vertices.
   *  - Overrule isosurface_separation_type, and always separate 
   *    diagonally opposite cube vertices.
   *  - Applies only to cubes (hypercubes).
   *  - Created as enum type, not boolean, to include
   *    UNDEFINED_SEP_OPP.
   */
  typedef enum { TRUE_SEPARATE_OPPOSITE, FALSE_SEPARATE_OPPOSITE, 
                 UNDEFINED_SEPARATE_OPPOSITE }
  SEPARATE_OPPOSITE_TYPE;

  namespace {
    static std::vector<IJK::ENUM_STR<SEPARATE_OPPOSITE_TYPE> > 
      _separate_opposite_type_pairs =
      { {TRUE_SEPARATE_OPPOSITE,"True"}, 
        {FALSE_SEPARATE_OPPOSITE,"False"},
        {UNDEFINED_SEPARATE_OPPOSITE, "Undefined"} };
  };

  /*!
   *  @brief Isosurface polytope orientation.
   *  - Orientation of polytopes (usually simplices) forming isosurface.
   *  - POSITIVE_ORIENT: Orient polytopes so that normal points
   *    to positive mesh region.
   *  - NEGATIVE_ORIENT: Orient polytopes so that normal points
   *    to negative mesh region.
   *  - NO_ORIENT: Polytopes not oriented.
   */
  typedef enum { POSITIVE_ORIENT, NEGATIVE_ORIENT, 
                 NO_ORIENT, UNDEFINED_ORIENT }
  ISO_POLY_ORIENTATION;

  namespace {
    static std::vector<IJK::ENUM_STR<ISO_POLY_ORIENTATION> > 
      _iso_poly_orientation_pairs =
      { {POSITIVE_ORIENT,"PositiveOrient"}, 
        {NEGATIVE_ORIENT,"NegativeOrient"},
        {NO_ORIENT, "NoOrient"},
        {UNDEFINED_ORIENT, "Undefined"} };
  };


  class ISOSURFACE_TABLE_PROPERTIES {
  protected:

    /// @brief Strings corresponding to enum LOOKUP_TABLE_TYPE.
    /// - Set table_type_list from _table_type_pairs[] in constructor.
    IJK::ENUM_LIST<LOOKUP_TABLE_TYPE> table_type_list;

    /// @brief Strings corresponding to enum ENCODING.
    /// - Set encoding_list from _encoding_pairs[] in constructor.
    IJK::ENUM_LIST<ENCODING> encoding_list;

    /// @brief Strings corresponding to enum GRID_VERTEX_LABEL_TYPE.
    /// - Set grid_vertex_label_type_list from 
    ///   _grid_vertex_label_type_pairs[] in constructor.
    IJK::ENUM_LIST<GRID_VERTEX_LABEL_TYPE> grid_vertex_label_type_list;

    /// @brief Strings corresponding to enum ISOSURFACE_TRIANGULATION_TYPE.
    /// - Set isosurface_triangluation_type_list from
    ///   _isosurface_triangulation_type_pairs[] in constructor.
    IJK::ENUM_LIST<ISOSURFACE_TRIANGULATION_TYPE>   
      isosurface_triangulation_type_list;

    /// @brief Strings corresponding to enum ISOSURFACE_SEPARATION_TYPE.
    /// - Set isosurface_separation_type_list from 
    ///   _isosurface_separation_type_pairs[] in constructor.
    IJK::ENUM_LIST<ISOSURFACE_SEPARATION_TYPE>   
      isosurface_separation_type_list;

    /// @brief Strings corresponding to enum SEPARATE_OPPOSITE_TYPE.
    /// - Set separate_opposite_type_list from 
    ///   _separte_opposite_type_pairs[] in constructor.
    IJK::ENUM_LIST<SEPARATE_OPPOSITE_TYPE>   
      separate_opposite_type_list;

    /// @brief Strings corresponding to enum ISO_POLY_ORIENTATION.
    /// - Set iso_poly_orientation_list from 
    ///   _iso_poly_orientation_pairs[] in constructor.
    IJK::ENUM_LIST<ISO_POLY_ORIENTATION> iso_poly_orientation_list;

    void Init();


  public:

  ISOSURFACE_TABLE_PROPERTIES():
    table_type_list(UNDEFINED_TABLE_TYPE, _table_type_pairs),
      encoding_list(UNKNOWN_ENCODING, _encoding_pairs),
      grid_vertex_label_type_list
      (UNDEFINED_GRID_VERTEX_LABEL_TYPE, _grid_vertex_label_type_pairs),
      isosurface_triangulation_type_list
      (UNDEFINED_ISOSURFACE_TRIANGULATION_TYPE, 
       _isosurface_triangulation_type_pairs),
      isosurface_separation_type_list
      (UNDEFINED_SEPARATION_TYPE, _isosurface_separation_type_pairs),
      separate_opposite_type_list
      (UNDEFINED_SEPARATE_OPPOSITE, _separate_opposite_type_pairs),
      iso_poly_orientation_list
      (UNDEFINED_ORIENT, _iso_poly_orientation_pairs)
        { Init(); };

    /// @brief Lookup table type.
    LOOKUP_TABLE_TYPE lookup_table_type;

    /// @brief Type of encoding, e.g. binary or base3.
    ENCODING encoding;           

    /// @brief Types of vertex labels in this table:
    /// - negative/positive, or negative/equals/positive.
    GRID_VERTEX_LABEL_TYPE grid_vertex_label_type;

    /// @brief Triangulation type of isosurface patches.
    /// - Indicate how triangulation was generated.
    ISOSURFACE_TRIANGULATION_TYPE isosurface_triangulation_type;

    /// @brief Indicate if isosurface patches separate positive 
    ///   or negative vertices.
    ISOSURFACE_SEPARATION_TYPE isosurface_separation_type;

    /*!
     *  @brief If true, always separate diagonally opposite cube vertices.
     *  - Overrule isosurface_separation_type, and always separate 
     *    diagonally opposite cube vertices.
     *  - Applies only to cubes (hypercubes).
     */
    SEPARATE_OPPOSITE_TYPE separate_opposite;

    /// @brief Orientation of isosurface polytopes (positive or negative).
    ISO_POLY_ORIENTATION iso_poly_orientation;


    /// @name: Get functions.
    ///@{

    /// @brief Return string of enum _table_type.
    std::string TableTypeString
    (const LOOKUP_TABLE_TYPE _table_type) const
    { return table_type_list.String(_table_type); }

    /// @brief Return string of lookup table type.
    std::string TableTypeString() const
    { return TableTypeString(lookup_table_type); }

    /// @brief Return lookup table type.
    LOOKUP_TABLE_TYPE TableType() const
    { return lookup_table_type; }

    /// @brief Return lookup table type defined by string s.
    LOOKUP_TABLE_TYPE TableType(const std::string & s) const
    { return table_type_list.EnumValue(s); }

    /// @brief Return true if lookup table type is undefined.
    bool IsTableTypeUndefined() const
    { return table_type_list.IsUndefined(lookup_table_type); }

    /// @brief Return string of enum _encoding.
    std::string EncodingString(const ENCODING _encoding) const
    { return encoding_list.String(_encoding); };
    
    /// @brief Return string of table encoding.
    std::string EncodingString() const
    { return EncodingString(encoding); }

    /// @brief Return encoding undefined value.
    ENCODING UndefinedEncoding() const
    { return encoding_list.UndefinedValue(); }

    /// @brief Return encoding.
    ENCODING Encoding() const
    { return encoding; }

    /// @brief Return encoding defined by string s.
    ENCODING Encoding(const std::string & s) const
    { return encoding_list.EnumValue(s); }

    /// @brief Return grid vertex label type.
    GRID_VERTEX_LABEL_TYPE GridVertexLabelType() const
    { return grid_vertex_label_type; }

    /// @brief Return grid vertex label type defined by string s.
    GRID_VERTEX_LABEL_TYPE GridVertexLabelType
      (const std::string & s) const
    { return grid_vertex_label_type_list.EnumValue(s); }

    /// @brief Return grid vertex label type undefined value.
    GRID_VERTEX_LABEL_TYPE UndefinedGridVertexLabelType()
    { return grid_vertex_label_type_list.UndefinedValue(); }

    /// @brief Return true if grid vertex label type is defined.
    bool IsGridVertexLabelTypeDefined() const
    { return grid_vertex_label_type_list.IsDefined
        (grid_vertex_label_type); }

    /// @brief Return true if grid vertex label type is undefined.
    bool IsGridVertexLabelTypeUndefined() const
    { return grid_vertex_label_type_list.IsUndefined
        (grid_vertex_label_type); }

    /// @brief Return string of grid vertex label type.
    std::string GridVertexLabelTypeString() const
      { return grid_vertex_label_type_list.String(grid_vertex_label_type); }

    /// @brief Return separation.
    ISOSURFACE_SEPARATION_TYPE SeparationType() const
    { return isosurface_separation_type; }

    /// @brief Return separation type defined by string s.
    ISOSURFACE_SEPARATION_TYPE SeparationType
      (const std::string & s) const
    { return isosurface_separation_type_list.EnumValue(s); }

    /// @brief Return true if isosurface separation type is defined.
    bool IsSeparationTypeDefined() const
    { return isosurface_separation_type_list.IsDefined
        (isosurface_separation_type); }

    /// @brief Return true if isosurface separation type is undefined.
    bool IsSeparationTypeUndefined() const
    { return isosurface_separation_type_list.IsUndefined
        (isosurface_separation_type); }

    /// @brief Return string of enum _sep_type.
    std::string SeparationTypeString
    (const ISOSURFACE_SEPARATION_TYPE _sep_type) const
    { return isosurface_separation_type_list.String(_sep_type); }
    
    /// @brief Return string of separation type.
    std::string SeparationTypeString() const
    { return SeparationTypeString(isosurface_separation_type); }

    /// @brief Return separation type file name label.
    std::string SeparationTypeLabel() const;

    /*!
     *  @brief Return opposite separation type.
     *  - SEPARATE_NEG maps to SEPARATE_POS.
     *  - SEPARATE_POS maps to SEPARATE_NEG.
     *  - Any other values are unchanged.
     */
    ISOSURFACE_SEPARATION_TYPE OppositeSeparationType() const;

    /// @brief Return triangulation type.
    ISOSURFACE_TRIANGULATION_TYPE TriangulationType() const
    { return isosurface_triangulation_type; }

    /// @brief Return triangulation type defined by string s.
    ISOSURFACE_TRIANGULATION_TYPE TriangulationType
      (const std::string & s) const
    { return isosurface_triangulation_type_list.EnumValue(s); }

    /// @brief Return true if triangulation type is defined.
    bool IsTriangulationTypeDefined() const
    { return isosurface_triangulation_type_list.IsDefined
        (isosurface_triangulation_type); }

    /// @brief Return true if triangulation type is undefined.
    bool IsTriangulationTypeUndefined() const
    { return isosurface_triangulation_type_list.IsUndefined
        (isosurface_triangulation_type); }

    /// @brief Return string of triangulation type.
    std::string TriangulationTypeString
    (const ISOSURFACE_TRIANGULATION_TYPE tri_type) const
    { return isosurface_triangulation_type_list.String(tri_type); }

    /// @brief Return string of triangulation type.
    std::string TriangulationTypeString() const
    { return TriangulationTypeString(isosurface_triangulation_type); }

    /// @brief Return triangulation type file name label.
    std::string TriangulationTypeLabel() const;

    /// @brief Return separate opposite.
    SEPARATE_OPPOSITE_TYPE SeparateOpposite() const
    { return separate_opposite; }

    /// @brief Return separate opposite defined by string s.
    SEPARATE_OPPOSITE_TYPE SeparateOpposite
      (const std::string & s) const
    { return separate_opposite_type_list.EnumValue(s); }

    /// @brief Return true if separate opposite is undefined.
    bool IsSeparateOppositeUndefined() const
    { return separate_opposite_type_list.IsUndefined(separate_opposite); }

    /// @brief Return string of enum _sep_opposite.
    std::string SeparateOppositeString
    (const SEPARATE_OPPOSITE_TYPE _sep_opposite) const
      { return separate_opposite_type_list.String(_sep_opposite); }
    
    /// @brief Return string of separate opposite
    std::string SeparateOppositeString() const
    { return SeparateOppositeString(separate_opposite); }

    /// @brief Set separate opposite to true or false.
    void SetSeparateOpposite(const bool flag) {
      if (flag) { separate_opposite = TRUE_SEPARATE_OPPOSITE; }
      else { separate_opposite = FALSE_SEPARATE_OPPOSITE; }
    };

    /// @brief Return isosurface polytope orientation.
    ISO_POLY_ORIENTATION IsoPolyOrientation() const
    { return iso_poly_orientation; }

    /// @brief Return isosurface polytope orientation defined by string s.
    ISO_POLY_ORIENTATION IsoPolyOrientation
      (const std:: string & s) const
    { return iso_poly_orientation_list.EnumValue(s); }

    /// @brief Return true if isosurface polytope orientation is undefined.
    bool IsIsoPolyOrientationUndefined() const
    { return iso_poly_orientation_list.IsUndefined(iso_poly_orientation); }

    /// @brief Return string of enum _orientation.
    std::string IsoPolyOrientationString
    (const ISO_POLY_ORIENTATION _orientation) const
      { return iso_poly_orientation_list.String(_orientation); }
    
    /// @brief Return string of isosurface polytope orientation.
    std::string IsoPolyOrientationString() const
    { return IsoPolyOrientationString(iso_poly_orientation); }

    /// @brief Return isosurface polytope orientation name label.
    std::string IsoPolyOrientationLabel() const;

    /*!
     *  @brief Return opposite orientation type.
     *  - NEGATIVE_ORIENT maps to POSITIVE_ORIENT.
     *  - POSITIVE_ORIENT maps to NEGATIVE_ORIENT.
     *  - Any other values are unchanged.
     */
    ISO_POLY_ORIENTATION OppositeIsoPolyOrientation() const;

    ///@}


    /// @name: Set/Copy functions.
    ///@{

    /// @brief Set encoding based on encoding_str.
    void SetEncoding(const std::string & encoding_str)
    { encoding = Encoding(encoding_str); }

    /// @brief Set grid vertex label type based on string s.
    void SetGridVertexLabelType(const std::string & s)
    { grid_vertex_label_type = GridVertexLabelType(s); }

    /// @brief Set triangulation type based on string s.
    void SetTriangulationType(const std::string & s)
    { isosurface_triangulation_type = TriangulationType(s); }

    /// @brief Set separation type based on string s.
    void SetSeparationType(const std::string & s)
    { isosurface_separation_type = SeparationType(s); }

    /// @brief Set separate opposite based on string s.
    void SetSeparateOpposite(const std::string & s)
    { separate_opposite = SeparateOpposite(s); }

    /// @brief Set isosurface polytope orientation based on string s.
    void SetIsoPolyOrientation(const std::string & s)
    { iso_poly_orientation = IsoPolyOrientation(s); }

    /// @brief Copy properties from isotable_properties.
    void Copy(const ISOSURFACE_TABLE_PROPERTIES & isotable_properties);

    ///@}


    /// @name Check routines.
    //@{

    /// @brief Return true if isotable has properties specified in properties.
    /// - If false, add to error message.
    /// - Ignores any fields in properties that are undefined or unknown.
    bool Check(const ISOSURFACE_TABLE_PROPERTIES & properties,
               IJK::ERROR & error) const;

    /// @brief Return true if lookup_table_type matches table_type.
    /// - Also returns true if table_type is undefined.
    bool CheckTableType
      (const LOOKUP_TABLE_TYPE, IJK::ERROR & error) const;

    /// @brief Return true if isosurface_triangulation_type matches tri_type.
    /// - Also returns true if tri_type is undefined or unknown.
    bool CheckTriangulationType
      (const ISOSURFACE_TRIANGULATION_TYPE tri_type,
       IJK::ERROR & error) const;

    /// @brief Return true if this->encoding matches encoding.
    /// - Also returns true if encoding is undefined.
    bool CheckEncoding(const ENCODING encoding, IJK::ERROR & error) const;

    /// @brief Return true if iso poly orientation matches _orientation.
    /// - Also returns true if _orientation is undefined or unknown.
    bool CheckOrientation(const ISO_POLY_ORIENTATION _orientation,
                          IJK::ERROR & error) const;

    /// @brief Return true if isosurface_separation_type matches sep_type.
    /// - Also returns true if sep_type is undefined or unknown.
    bool CheckSeparationType(const ISOSURFACE_SEPARATION_TYPE sep_type,
                             IJK::ERROR & error) const;

    /// @brief Return true if this->separate_opposite matches sep_opposite.
    /// - Also returns true if sep_opposite is undefined or unknown.
    bool CheckSeparateOpposite
      (const SEPARATE_OPPOSITE_TYPE sep_opposite, IJK::ERROR & error) const;
    
    //@}


    /// @name Print Routines (mainly for debugging).
    //@{

    ///@ brief Print table properties.
    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

    //@}

  };


  // ***************************************************************
  // ISOSURFACE_TABLE_PROPERITES member functions
  // ***************************************************************


  template <typename OSTREAM_TYPE, typename STYPE0>
    void ISOSURFACE_TABLE_PROPERTIES::Print
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    out << line_prefix;
    out << "Encoding: " << EncodingString();
    out << "\n";

    if (!IsGridVertexLabelTypeUndefined()) {
      out << line_prefix;
      out << "Poly vertex label type: " << GridVertexLabelTypeString();
      out << "\n";
    }

    if (!IsTriangulationTypeUndefined()) {
      out << line_prefix;
      out << "Triangulation type: " << TriangulationTypeString();
      out << "\n";
    } 

    if (!IsSeparationTypeUndefined()) {
      out << line_prefix;
      out << "Separation type: " << SeparationTypeString();
      out << "\n";
    }

    if (!IsSeparateOppositeUndefined()) {
      out << line_prefix;
      out << "Separate opposite: " << SeparateOppositeString();
      out << "\n";
    }

    if (!IsIsoPolyOrientationUndefined()) {
      out << line_prefix;
      out << "Isosurface polytope orientation: " 
          << IsoPolyOrientationString();
      out << "\n";
    } 
  }

}

#endif
