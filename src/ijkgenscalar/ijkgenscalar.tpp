/*!
 *  @file ijkgenscalar.tpp
 *  @brief Generate a scalar field.
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2023 Rephael Wenger

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

#ifndef _IJKGENSCALAR_TPP_
#define _IJKGENSCALAR_TPP_

#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "ijk.tpp"
#include "ijkgengeom.tpp"
#include "ijkstring.tpp"


namespace IJKGENSCALAR {

  // ********************************************************
  // TEMPLATE CLASS FIELD_INFO_T
  // ********************************************************

  /// @brief Information about a field.
  /// - Field name, scalar field generator, gradient field generator, etc.
  template <typename OBJECT_PROPERTIES_TYPE,
            typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>
  class FIELD_INFO_T:public IJKGENGEOM::GEOM_INFO_T<OBJECT_PROPERTIES_TYPE> {

  protected:
    void Init();

  public:
    typedef void (*FUNCTION_PTR)
      (const OBJECT_PROPERTIES_TYPE & prop, SCALAR_GRID_TYPE & grid);
    typedef void (*GRADIENT_FUNCTION_PTR)
      (const OBJECT_PROPERTIES_TYPE & prop,
       SCALAR_GRID_TYPE & scalar_grid, GRADIENT_GRID_TYPE & gradient_grid);

  public:
    FUNCTION_PTR function_ptr;
    GRADIENT_FUNCTION_PTR gradient_function_ptr;
    bool flag_maxval;           ///< If true, field has a maximum scalar value.
    bool flag_max_exponent;     ///< If true, field has a maximum exponent.

  public:
    FIELD_INFO_T() { Init(); };

    void Set(const std::string & field_name, FUNCTION_PTR f_ptr)
    {
      this->name = field_name;
      function_ptr = f_ptr;
      gradient_function_ptr = NULL;
      SetAllFlags(false);
    }

    void Set(const std::string & field_name, FUNCTION_PTR f_ptr,
             GRADIENT_FUNCTION_PTR grad_f_ptr)
    {
      this->name = field_name;
      function_ptr = f_ptr;
      gradient_function_ptr = grad_f_ptr;
      SetAllFlags(false);
    }

    void SetAllFlags(const bool flag);
    void SetRandomFlags();
    void SetRandomPowerFlags();
    void SetScalarSetFlags();

    /// Return true if gradient_function_ptr is not NULL.
    bool IsGradientImplemented() const
    {
      if (gradient_function_ptr == NULL) 
        { return(false); }
      else
        { return(true); }
    }

  };


  // ********************************************************
  // TEMPLATE CLASS FIELD_OBJECT_PROPERTIES_T
  // ********************************************************

  /// Object properties.
  template <typename DIM_TYPE, typename COORD_TYPE, typename DIR_TYPE,
            typename RADIUS_TYPE, typename DIFF_TYPE, typename ANGLE_TYPE,
            typename SCALAR_TYPE, typename NUM_TYPE>
  class FIELD_OBJECT_PROPERTIES_T:
    public IJKGENGEOM::OBJECT_PROPERTIES_T
    <DIM_TYPE,COORD_TYPE,DIR_TYPE,RADIUS_TYPE,
     DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE> {
    
  protected:
    void Init();

  public:
    FIELD_OBJECT_PROPERTIES_T(){ Init(); };
    bool gradient_discontinuity_zero;

    /// Copy.
    void Copy(const FIELD_OBJECT_PROPERTIES_T
              <DIM_TYPE, COORD_TYPE, DIR_TYPE, RADIUS_TYPE,
              DIFF_TYPE, ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE> 
              & right)
    {
      if (&right != this) {         // avoid self-assignment
        *this = right;
      }
    };

  };


  // ********************************************************
  // TEMPLATE CLASS SET_SCALAR_DATA_T
  // ********************************************************

  /// Information to set cube.
  template <typename CUBE_INDEX_TYPE, typename ISOTABLE_INDEX_TYPE,
            typename VERTEX_INDEX_TYPE, typename DIRECTION_TYPE>
  class SET_SCALAR_DATA_T {

  public:
    /// Type of data which determines the scalar values.
    typedef enum { SET_CUBE_ISOTABLE_INDEX, 
                   SET_POS_VERTEX,
                   SET_POS_EDGE,
                   SET_POS_FACET,
                   SET_POS_FACET_DIAGONAL,
                   UNDEFINED_SET_CUBE } SET_SCALAR_TYPE;

  public:
    /// Index of cube.
    CUBE_INDEX_TYPE cube_index;

    /// Type of data which determines the scalar values.
    SET_SCALAR_TYPE set_type;

    ISOTABLE_INDEX_TYPE isotable_index;
    VERTEX_INDEX_TYPE vertex_index;
    DIRECTION_TYPE direction;

    /// Index of vertex within facet or cube.
    VERTEX_INDEX_TYPE jth_vertex;

    // Query functions

    bool IsTypeSetIsotableIndex() const
    { return((set_type == SET_CUBE_ISOTABLE_INDEX)); }
    bool IsTypeSetPosVertex() const
    { return((set_type == SET_POS_VERTEX)); }
    bool IsTypeSetPosEdge() const
    { return((set_type == SET_POS_EDGE)); }
    bool IsTypeSetPosFacet() const
    { return((set_type == SET_POS_FACET)); }
    bool IsTypeSetPosFacetDiagonal() const
    { return((set_type == SET_POS_FACET_DIAGONAL)); }

    // Set functions

    /// Set cube configuration to isotable_index.
    void SetCubeIsotableIndex
    (const CUBE_INDEX_TYPE cube_index, 
     const ISOTABLE_INDEX_TYPE isotable_index);

    /// Set vertex to pos.
    void SetPosVertex(const VERTEX_INDEX_TYPE iv);

    /// Set edge to pos.
    void SetPosEdge(const VERTEX_INDEX_TYPE iend, 
                    const DIRECTION_TYPE edge_dir);

    /// Set facet to pos.
    void SetPosFacet(const VERTEX_INDEX_TYPE iv0, 
                     const DIRECTION_TYPE orth_dir);

    /// Set facet diagonal to positive.
    /// @param j Set facet diagonal containing j'th vertex of facet 
    ///          to positive.
    void SetPosFacetDiagonal(const VERTEX_INDEX_TYPE iv0, 
                             const DIRECTION_TYPE orth_dir,
                             const VERTEX_INDEX_TYPE j);
  };


  // ********************************************************
  // TEMPLATE CLASS FIELD_PARAM_T
  // ********************************************************

  /// Input parameters determining field properties.
  template <typename SET_SCALAR_DATA_TYPE,
            typename GEOM_PARAM_TYPE, 
            typename SCALAR_TYPE, 
            typename MIN_MAX_TYPE,
            typename MAX_EXPONENT_TYPE,
            typename ISOTABLE_INDEX_TYPE>
  class FIELD_PARAM_T:public GEOM_PARAM_TYPE {

  public:
    typedef typename GEOM_PARAM_TYPE::NUMBER_TYPE NUMBER_TYPE;

  protected:
    void Init();

  public:
    FIELD_PARAM_T(){ Init(); };

    IJK::SET_VALUE<MIN_MAX_TYPE> maxval;      ///< Max scalar value.

    /// Max exponent.
    IJK::SET_VALUE<MAX_EXPONENT_TYPE> max_exponent;
    
    /// Isosurface lookup table index.
    IJK::SET_VALUE<ISOTABLE_INDEX_TYPE> isotable_index;

    /// List of set instructions.
    std::vector<SET_SCALAR_DATA_TYPE> set_scalar_list;

    /// Scalar value.  (For use with constant scalar field.)
    std::vector<SCALAR_TYPE> scalar_value;

    /// If true, set all boundary vertices to zero.
    bool flag_set_boundary2zero;

    template <typename NTYPE>
    SCALAR_TYPE ScalarValue(const NTYPE i) const
    { return(scalar_value[i]); }

    MIN_MAX_TYPE MaxVal() const
    { return(maxval.Value()); }

    MAX_EXPONENT_TYPE MaxExponent() const
    { return(max_exponent.Value()); }

    ISOTABLE_INDEX_TYPE IsotableIndex() const
    { return(isotable_index.Value()); }

    void SetFieldIndex(const NUMBER_TYPE field_index)
    { this->SetGeomInfoIndex(field_index); }

    NUMBER_TYPE FieldIndex() const
    { return(this->GeomInfoIndex()); }

    template <typename STYPE>
    void SetScalarValue(const STYPE s)
    { scalar_value.push_back(s); }

    template <typename CUBE_INDEX_TYPE2, typename ISOTABLE_INDEX_TYPE2>
    void SetCubeIsotableIndex
    (const CUBE_INDEX_TYPE2 cube_index,
     const ISOTABLE_INDEX_TYPE2 isotable_index);

    template <typename VERTEX_INDEX_TYPE2>
    void SetScalarPosVertex
    (const VERTEX_INDEX_TYPE2 iv);

    template <typename VERTEX_INDEX_TYPE2, typename DIRECTION_TYPE2>
    void SetScalarPosEdge(const VERTEX_INDEX_TYPE2 iv,
                          const DIRECTION_TYPE2 edge_dir);

    template <typename VERTEX_INDEX_TYPE2, typename DIRECTION_TYPE2>
    void SetScalarPosFacet(const VERTEX_INDEX_TYPE2 iv0,
                           const DIRECTION_TYPE2 orth_dir);

    template <typename VERTEX_INDEX_TYPE2, typename DIRECTION_TYPE2,
              typename ITYPE>
    void SetScalarPosFacetDiagonal
    (const VERTEX_INDEX_TYPE2 iv, const DIRECTION_TYPE2 orth_dir,
     const ITYPE j);

    template <typename NTYPE>
    void GetSetCubeIsotableIndexString(const NTYPE i, std::string & s) const;

    template <typename NTYPE>
    void GetSetScalarPosVertexString(const NTYPE i, std::string & s) const;

    template <typename NTYPE>
    void GetSetScalarPosEdgeString(const NTYPE i, std::string & s) const;

    template <typename NTYPE>
    void GetSetScalarPosFacetString(const NTYPE i, std::string & s) const;

    template <typename NTYPE>
    void GetSetScalarPosFacetDiagonalString
    (const NTYPE i, std::string & s) const;


    void GetMaxValString(std::string & s) const;

    void GetRandomBaseString(std::string & s) const;

    void GetMaxExponentString(std::string & s) const;
};

  // ********************************************************
  // MACROS FOR TEMPLATE MEMBER FUNCTIONS
  // ********************************************************

  /// Header for template FIELD_INTO_T
#define _FIELD_INFO_T_HEADER_                                        \
  template <typename OBJECT_PROPERTIES_TYPE,                         \
            typename SCALAR_GRID_TYPE, typename GRADIENT_GRID_TYPE>  \
  void FIELD_INFO_T<OBJECT_PROPERTIES_TYPE,                          \
                    SCALAR_GRID_TYPE,GRADIENT_GRID_TYPE>::


#define _FIELD_OBJECT_PROPERTIES_T_HEADER_                           \
  template <typename DIM_TYPE,                                       \
            typename COORD_TYPE, typename DIR_TYPE,                  \
            typename RADIUS_TYPE, typename DIFF_TYPE,                \
            typename ANGLE_TYPE, typename SCALAR_TYPE,               \
            typename NUM_TYPE>                                       \
  void FIELD_OBJECT_PROPERTIES_T<DIM_TYPE, COORD_TYPE, DIR_TYPE,     \
                           RADIUS_TYPE, DIFF_TYPE, ANGLE_TYPE,       \
                           SCALAR_TYPE, NUM_TYPE>::

#define _SET_SCALAR_DATA_T_HEADER_                                   \
  template <typename CUBE_INDEX_TYPE,                                \
            typename ISOTABLE_INDEX_TYPE,                            \
            typename VERTEX_INDEX_TYPE,                              \
            typename DIRECTION_TYPE>                                 \
  void SET_SCALAR_DATA_T<CUBE_INDEX_TYPE, ISOTABLE_INDEX_TYPE,       \
                         VERTEX_INDEX_TYPE, DIRECTION_TYPE>::

#define _FIELD_PARAM_T_HEADER_A_                                     \
  template <typename SET_SCALAR_DATA_TYPE,                           \
            typename GEOM_PARAM_TYPE,                                \
            typename SCALAR_TYPE,                                    \
            typename MIN_MAX_TYPE,                                   \
	    typename MAX_EXPONENT_TYPE,				     \
            typename ISOTABLE_INDEX_TYPE>

#define _FIELD_PARAM_T_HEADER_B_                                     \
  void FIELD_PARAM_T<SET_SCALAR_DATA_TYPE, GEOM_PARAM_TYPE,          \
                     SCALAR_TYPE, MIN_MAX_TYPE,                      \
		     MAX_EXPONENT_TYPE, 			     \
                     ISOTABLE_INDEX_TYPE>::


#define _FIELD_PARAM_T_HEADER_                                       \
  _FIELD_PARAM_T_HEADER_A_                                           \
  _FIELD_PARAM_T_HEADER_B_


  // ********************************************************
  // TEMPLATE CLASS FIELD_INFO MEMBER FUNCTIONS
  // ********************************************************

  _FIELD_INFO_T_HEADER_
    Init()
    {
      SetAllFlags(false);
    }


  _FIELD_INFO_T_HEADER_
  SetRandomFlags()
  {
    this->flag_random_seed = true;
    flag_maxval = true;
    flag_max_exponent = false;
    this->flag_allow_multi_centers = false;
    this->flag_allow_tilt = false;
  }

  _FIELD_INFO_T_HEADER_
  SetRandomPowerFlags()
  {
    this->flag_random_seed = true;
    flag_maxval = false;
    flag_max_exponent = true;
    this->flag_allow_multi_centers = false;
    this->flag_allow_tilt = false;
  }

  _FIELD_INFO_T_HEADER_
  SetScalarSetFlags()
  {
    flag_maxval = true;
  }

  _FIELD_INFO_T_HEADER_
  SetAllFlags(const bool flag)
  {
    this->IJKGENGEOM::GEOM_INFO_T<OBJECT_PROPERTIES_TYPE>::SetAllFlags(flag);
    this->flag_maxval = false;
    this->flag_max_exponent = false;
  }


  // ********************************************************
  // TEMPLATE CLASS OBJECT_PROPERTIES_T MEMBER FUNCTIONS
  // ********************************************************

  _FIELD_OBJECT_PROPERTIES_T_HEADER_
  Init()
  {
    this->gradient_discontinuity_zero = true;
  }


  // ********************************************************
  // TEMPLATE CLASS SET_SCALAR_DATA_T MEMBER FUNCTIONS
  // ********************************************************

  _SET_SCALAR_DATA_T_HEADER_
  SetCubeIsotableIndex
  (const CUBE_INDEX_TYPE cube_index, 
   const ISOTABLE_INDEX_TYPE isotable_index)
  {
    this->set_type = SET_CUBE_ISOTABLE_INDEX;
    this->cube_index = cube_index;
    this->isotable_index = isotable_index;
  }


  _SET_SCALAR_DATA_T_HEADER_
  SetPosVertex(const VERTEX_INDEX_TYPE iv)
  {
    this->set_type = SET_POS_VERTEX;
    this->vertex_index = iv;
  }


  _SET_SCALAR_DATA_T_HEADER_
  SetPosEdge(const VERTEX_INDEX_TYPE iend,
             const DIRECTION_TYPE edge_dir)
  {
    this->set_type = SET_POS_EDGE;
    this->vertex_index = iend;
    this->direction = edge_dir;
  }


  _SET_SCALAR_DATA_T_HEADER_
  SetPosFacet(const VERTEX_INDEX_TYPE iv0, 
              const DIRECTION_TYPE orth_dir)
  {
    this->set_type = SET_POS_FACET;
    this->vertex_index = iv0;
    this->direction = orth_dir;
  }


  _SET_SCALAR_DATA_T_HEADER_
  SetPosFacetDiagonal(const VERTEX_INDEX_TYPE iv0, 
                      const DIRECTION_TYPE orth_dir,
                      const VERTEX_INDEX_TYPE j)
  {
    this->set_type = SET_POS_FACET_DIAGONAL;
    this->vertex_index = iv0;
    this->direction = orth_dir;
    this->jth_vertex = j;
  }


  // ********************************************************
  // TEMPLATE CLASS FIELD_PARAM_T MEMBER FUNCTIONS
  // ********************************************************

  _FIELD_PARAM_T_HEADER_
  Init()
  {
    this->gradient_discontinuity_zero = true;
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename CUBE_INDEX_TYPE2, typename ISOTABLE_INDEX_TYPE2>
  _FIELD_PARAM_T_HEADER_B_
  SetCubeIsotableIndex
  (const CUBE_INDEX_TYPE2 cube_index,
   const ISOTABLE_INDEX_TYPE2 isotable_index)
  {
    SET_SCALAR_DATA_TYPE set_scalar_data;

    set_scalar_data.SetCubeIsotableIndex(cube_index, isotable_index);
    this->set_scalar_list.push_back(set_scalar_data);
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename VERTEX_INDEX_TYPE2>
  _FIELD_PARAM_T_HEADER_B_
  SetScalarPosVertex(const VERTEX_INDEX_TYPE2 iv)
  {
    SET_SCALAR_DATA_TYPE set_scalar_data;

    set_scalar_data.SetPosVertex(iv);
    this->set_scalar_list.push_back(set_scalar_data);
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename VERTEX_INDEX_TYPE2, typename DIRECTION_TYPE2>
  _FIELD_PARAM_T_HEADER_B_
  SetScalarPosEdge(const VERTEX_INDEX_TYPE2 iv,
                   const DIRECTION_TYPE2 edge_dir)
  {
    SET_SCALAR_DATA_TYPE set_scalar_data;

    set_scalar_data.SetPosEdge(iv, edge_dir);
    this->set_scalar_list.push_back(set_scalar_data);
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename VERTEX_INDEX_TYPE2, typename DIRECTION_TYPE2>
  _FIELD_PARAM_T_HEADER_B_
  SetScalarPosFacet(const VERTEX_INDEX_TYPE2 iv,
                    const DIRECTION_TYPE2 orth_dir)
  {
    SET_SCALAR_DATA_TYPE set_scalar_data;

    set_scalar_data.SetPosFacet(iv, orth_dir);
    this->set_scalar_list.push_back(set_scalar_data);
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename VERTEX_INDEX_TYPE2, 
            typename DIRECTION_TYPE2,
            typename ITYPE>
  _FIELD_PARAM_T_HEADER_B_
  SetScalarPosFacetDiagonal(const VERTEX_INDEX_TYPE2 iv,
                            const DIRECTION_TYPE2 orth_dir,
                            const ITYPE j)
  {
    SET_SCALAR_DATA_TYPE set_scalar_data;

    set_scalar_data.SetPosFacetDiagonal(iv, orth_dir, j);
    this->set_scalar_list.push_back(set_scalar_data);
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename NTYPE>
  _FIELD_PARAM_T_HEADER_B_
  GetSetCubeIsotableIndexString(const NTYPE i, std::string & s) const
  {
    std::string s0, s1;

    if (set_scalar_list[i].IsTypeSetIsotableIndex()) {
      IJK::val2string(set_scalar_list[i].cube_index, s0);
      IJK::val2string(set_scalar_list[i].isotable_index, s1);
    }
    s = s0 + " " + s1;
  }

  _FIELD_PARAM_T_HEADER_A_
  template <typename NTYPE>
  _FIELD_PARAM_T_HEADER_B_
  GetSetScalarPosVertexString(const NTYPE i, std::string & s) const
  {
    std::string s0;

    if (set_scalar_list[i].IsTypeSetPosVertex()) {
      IJK::val2string(set_scalar_list[i].vertex_index, s0);
    }
    s = s0;
  }

  _FIELD_PARAM_T_HEADER_A_
  template <typename NTYPE>
  _FIELD_PARAM_T_HEADER_B_
  GetSetScalarPosEdgeString(const NTYPE i, std::string & s) const
  {
    std::string s0, s1;

    if (set_scalar_list[i].IsTypeSetPosEdge()) {
      IJK::val2string(set_scalar_list[i].vertex_index, s0);
      IJK::val2string(set_scalar_list[i].direction, s1);
    }
    s = s0 + " " + s1;
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename NTYPE>
  _FIELD_PARAM_T_HEADER_B_
  GetSetScalarPosFacetString(const NTYPE i, std::string & s) const
  {
    std::string s0, s1;

    if (set_scalar_list[i].IsTypeSetPosFacet()) {
      IJK::val2string(set_scalar_list[i].vertex_index, s0);
      IJK::val2string(set_scalar_list[i].direction, s1);
    }
    s = s0 + " " + s1;
  }


  _FIELD_PARAM_T_HEADER_A_
  template <typename NTYPE>
  _FIELD_PARAM_T_HEADER_B_
  GetSetScalarPosFacetDiagonalString(const NTYPE i, std::string & s) const
  {
    std::string s0, s1, s2;

    if (set_scalar_list[i].IsTypeSetPosFacetDiagonal()) {
      IJK::val2string(set_scalar_list[i].vertex_index, s0);
      IJK::val2string(set_scalar_list[i].direction, s1);
      IJK::val2string(set_scalar_list[i].jth_vertex, s2);
    }
    s = s0 + " " + s1 + " " + s2;
  }


  _FIELD_PARAM_T_HEADER_
  GetMaxValString(std::string & s) const
  {
    IJK::val2string(MaxVal(), s);
  }

  _FIELD_PARAM_T_HEADER_
  GetRandomBaseString(std::string & s) const
  {
    IJK::val2string(this->RandomBase(), s);
  }

  _FIELD_PARAM_T_HEADER_
  GetMaxExponentString(std::string & s) const
  {
    IJK::val2string(MaxExponent(), s);
  }

};

#endif
