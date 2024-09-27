/*!
 *  @file ijk.tpp
 *  @brief ijk templates defining general ijk objects, i.e. 
 *    classes BOX and ERROR.
 *  - General purpose ijk classes and routines.
 *  - Box, line segments, array and error classes.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008-2023 Rephael Wenger

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

#ifndef _IJK_
#define _IJK_

#include <bitset>
#include <cstddef>
#include <limits>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <utility>


/// @brief IJK namespace for all generic/shared IJK code.
namespace IJK {

  // *****************************************************************
  // CLASSES
  // *****************************************************************

  class ERROR;


  // *****************************************************************
  // TEMPLATE CLASS ARRAY AND ARRAY_L
  // *****************************************************************

  // *** Probably should be replaced by std::array
  /*!
   *  @brief Simple array class for creating static arrays.
   *  - Making this a class guarantees that when a program leaves
   *  the function where this class is declared, the array memory is freed.
   *  - Note: If ARRAY is a an array of pointers, the Free() operation 
   *    and the destructor ARRAY~ DO NOT free any memory pointed
   *    to by this class.
   */
  template <typename ETYPE> class ARRAY {
  protected:
    ETYPE * element;

    template <typename LTYPE>
    void Init(const LTYPE array_length)
    { element = new ETYPE[array_length]; }

    template <typename LTYPE>
    void Init(const LTYPE array_length, const ETYPE init_value)
    { 
      Init(array_length);
      for (LTYPE i = 0; i < array_length; i++) { element[i] = init_value; };
    }

  public:
    template <typename LTYPE>
    ARRAY(const LTYPE array_length) // constructor
    { Init(array_length); }
    template <typename LTYPE>
    ARRAY(const LTYPE array_length, const ETYPE init_value) // constructor
    { Init(array_length, init_value); }
    ~ARRAY();

    // get functions
    ETYPE * Ptr() { return(element); };
    const ETYPE * PtrConst() const { return(element); }
    template <typename ITYPE>
    ETYPE & operator [] (const ITYPE i) { return(*(element+i)); }
    template <typename ITYPE>
    ETYPE operator [] (const ITYPE i) const { return(*(element+i)); }

    // Free function
    void Free();
  };

  /// Class array with length stored
  template <typename ETYPE, typename LTYPE>
  class ARRAY_L:public ARRAY<ETYPE> {
  protected:
    LTYPE length;

    void Init(const LTYPE length) 
    { this->length = length; }

  public:
    ARRAY_L(const LTYPE length):ARRAY<ETYPE>(length) 
      // constructor
    { Init(length); };

    ARRAY_L(const LTYPE length, const ETYPE init_value):
      ARRAY<ETYPE>(length, init_value)
      // constructor
    { Init(length); };

    // get functions
    LTYPE Length() const { return(length); };
    ETYPE * End() { return(this->element+length); };
  };


  // *****************************************************************
  // TEMPLATE CLASS SET_VALUE
  // *****************************************************************

  /// Class SET_VALUE with flag is_set to indicate if value has been set.
  template <typename T>
  class SET_VALUE {

  protected:
    T v;
    bool is_set;

    void Init() 
    { is_set = false; }

    void Init(const T value)
    {
      v = value;
      is_set = false;
    }

  public:
    SET_VALUE(){ Init(); };
    SET_VALUE(const T value) { Init(value); };

    /// Set v to value and is_set to true.
    void Set(const T value)
    {
      v = value;
      is_set = true;
    }

    /// Set is_set to false.
    void Unset()
    { is_set = false; }

    bool IsSet() const { return(is_set); };
    T Value() const { return(v); };    
  };


  /// @brief Template specification SET_VALUE<bool> with additional 
  ///   member function IsSetAndTrue().
  class BOOLEAN_SET_VALUE:public SET_VALUE<bool> {

  public:
    BOOLEAN_SET_VALUE() {};
    BOOLEAN_SET_VALUE(const bool value):SET_VALUE<bool>(value) {};

    bool IsSetAndTrue() const
    { return((this->IsSet() && Value())); }
  };


  // *****************************************************************
  // TEMPLATE CLASS CONSTANT
  // *****************************************************************

  /// Class CONSTANT always returns the same value.
  template <typename ITYPE, typename CTYPE>
  class CONSTANT {
  protected:
    CTYPE c;           ///< The constant value.

  public:
    CONSTANT(const CTYPE c) { this->c = c; };

    // Get functions.
    CTYPE operator [] (const ITYPE i) const { return(c); };
    CTYPE operator () (const ITYPE i) const { return(c); };
  };


  // *****************************************************************
  // TEMPLATE CLASS BOX
  // *****************************************************************

  /// @brief Axis-parallel box data structure.
  /// - Represents box by minimum and maximum coordinates.
  template <typename COORD_TYPE, int INIT_DIM=0> class BOX {
  protected:
    int dimension;
    COORD_TYPE * min_coord;
    COORD_TYPE * max_coord;

    void Init(const int dim);
    void FreeAll();

  public:
    BOX() { Init(INIT_DIM); };
    BOX(const int dimension);
    ~BOX() { FreeAll(); };

    /// Copy constructor.
    BOX(const BOX<COORD_TYPE,INIT_DIM> & box);

    /// Copy assignment.
    const BOX & operator = (const BOX<COORD_TYPE,INIT_DIM> &);

    // get functions

    /// Return box dimension.
    int Dimension() const { return(dimension); };

    /// Return d'th minimum coordinate.
    COORD_TYPE MinCoord(const int d) const
    { return(min_coord[d]); };

    /// Return constant pointer to array of minimum coordinates.
    const COORD_TYPE * MinCoord() const
    { return(min_coord); };

    /// Return d'th maximum coordinate.
    COORD_TYPE MaxCoord(const int d) const
    { return(max_coord[d]); };

    /// Return constant pointer to array of maximum coordinates.
    const COORD_TYPE * MaxCoord() const
    { return(max_coord); };


    /// Return number of vertices along axis d if first vertex starts at min_coord[d].
    COORD_TYPE AxisSize(const int d) const
    { return(max_coord[d]+1-min_coord[d]); }

    /// Return width of box in direction d.
    // *** NOT YET TESTED ***
    COORD_TYPE Width(const int d) const
    { return(max_coord[d]-min_coord[d]); }

    /// Return true if MinCoord(d) <= MaxCoord(d).
    bool IsMinLEMaxCoord(const int d) const
    { return((MinCoord(d) <= MaxCoord(d))); }

    /// If i==0, return MinCoord(d).  If i==1, return MaxCoord(d).
    COORD_TYPE Coord(const int i, const int d) const
    {
      if (i==0) { return(MinCoord(d)); }
      else { return(MaxCoord(d)); }
    }

    /// If i==0, return MinCoord().  If i==1, return MaxCoord().
    const COORD_TYPE * Coord(const int i) const
    {
      if (i==0) { return(MinCoord()); }
      else { return(MaxCoord()); }
    }

    /// Return true if box contains coord[].
    template <typename CTYPE2>
    bool Contains(const CTYPE2 * coord) const;

    // set functions

    /// Set dimensions to d.
    void SetDimension(const int d);

    /// Set d'th minimum coordinate to c.
    void SetMinCoord(const int d, const COORD_TYPE c)
    { min_coord[d] = c; };

    /// Set d'th maximum coordinate to c.
    void SetMaxCoord(const int d, const COORD_TYPE c)
    { max_coord[d] = c; };

    /// Set min_coord[] to coord[].
    template <typename CTYPE2>
    void SetMinCoord(const CTYPE2 * coord);

    /// Set max_coord[] to coord[].
    template <typename CTYPE2>
    void SetMaxCoord(const CTYPE2 * coord);


    /// Set min_coord[] to minc[] and max_coord[] to maxc[].
    template <typename CTYPE2, typename CTYPE3>
    void SetCoord(const CTYPE2 * minc, const CTYPE3 * maxc);

    /// Set all minimum coordinates to c.
    void SetAllMinCoord(const COORD_TYPE c);

    /// Set all maximum coordinates to c.
    void SetAllMaxCoord(const COORD_TYPE c);

    /// Set both min and max coordinate d.
    void SetMinMaxCoord
    (const int d, const COORD_TYPE minc, const COORD_TYPE maxc)
    {
      SetMinCoord(d, minc);
      SetMaxCoord(d, maxc);
    }

    /// Extend the min/max coords so box includes coord.
    template <typename CTYPE2>
    void Extend(const CTYPE2 * coord);

    /// Clip box by increasing min coord.
    template <typename CTYPE2>
    void ClipMin(const int d, const CTYPE2 c)
    { if (min_coord[d] < c) { min_coord[d] = c; } }

    /// Clip box by decreasing max coord.
    template <typename CTYPE2>
    void ClipMax(const int d, const CTYPE2 c)
    { if (max_coord[d] > c) { max_coord[d] = c; } }

    /// Clip box.
    template <typename CTYPE2, typename CTYPE3>
    void Clip(const int d, const CTYPE2 minc, const CTYPE3 maxc)
    {
      ClipMin(d, minc);
      ClipMax(d, maxc);
    }


    // Check functions.

    /*!
     *  @brief Return true if min coordinates are at most max coordinates.
     *  - Return false and add messages to error, if min coordinates
     *    are greater than max coordinates.
     */
    bool CheckCoord(ERROR & error) const;


    // Print routines (mainly for debugging.)

  protected:

    /*!
     *  @brief Print coordinate coord[]. Mainly for debugging.
     *  - Print coordinates coord[].
     *  - Extended version.  Define left/right delimiters and separator.
     *  - Number of coordinates = Dimension().
     *  @param c0 Left delimiter.
     *  @param c1 Separator.
     *  @param c2 Right delimiter.
     */
    template <typename OSTREAM_TYPE>
    void _PrintCoordX
    (OSTREAM_TYPE & out, const COORD_TYPE * coord,
     const char c0, const char c1, const char c2) const;

    /*!
     *  @brief Print coordinate coord[]. Mainly for debugging.
     *  - Version using delimiters '(' and ')' and separator ','.
     *  - Number of coordinates = Dimension().
     */
    template <typename OSTREAM_TYPE>
    void _PrintCoord(OSTREAM_TYPE & out, const COORD_TYPE * coord) const
    { _PrintCoordX(out, coord, '(', ',', ')'); }

    /*!
     *  @brief Print coordinates. (Set prefix and suffix.)
     *  - Mainly for debugging.
     *  - Version adding prefix and suffix strings.
     *  @param s0 Prefix string.
     *  @param s1 Suffix string.
     */
    template <typename OSTREAM_TYPE,
	      typename STYPE0, typename STYPE1>
    void _PrintCoord
    (OSTREAM_TYPE & out,
     const STYPE0 & s0, const COORD_TYPE * coord,
     const STYPE1 & s1) const;

  public:

    /// Print minimum coordinates. Mainly for debugging.
    template <typename OSTREAM_TYPE>
    void PrintMinCoord(OSTREAM_TYPE & out) const
    { _PrintCoord(out, min_coord); }

    /// Print maximum coordinates. Mainly for debugging.
    template <typename OSTREAM_TYPE>
    void PrintMaxCoord(OSTREAM_TYPE & out) const
    { _PrintCoord(out, max_coord); }

    /// @brief Print minimum coordinates. Mainly for debugging.
    /// - Version adding prefix and suffix strings.
    template <typename OSTREAM_TYPE,
	      typename STYPE0, typename STYPE1>
    void PrintMinCoord
    (OSTREAM_TYPE & out,
     const STYPE0 & s0, const STYPE1 & s1) const
    { _PrintCoord(out, s0, min_coord, s1); }

    /// @brief Print maximum coordinates. Mainly for debugging.
    /// - Version adding prefix and suffix strings.
    template <typename OSTREAM_TYPE,
	      typename STYPE0, typename STYPE1>
    void PrintMaxCoord
    (OSTREAM_TYPE & out,
     const STYPE0 & s0, const STYPE1 & s1) const
    { _PrintCoord(out, s0, max_coord, s1); }

  };


  // *****************************************************************
  //! @name CONVERT/INITIALIZE C++ STL VECTORS.
  // *****************************************************************

  //@{

  /// @brief Returns C++ pointer to C array storing in vector data.
  /// - If vector is empty, returns NULL pointer.
  /// - If vector is not empty, returns pointer to a C array.
  template <typename T>
  const T * vector2pointer(const std::vector<T> & v)
  {
    if (v.empty()) { return(NULL); }
    else { return v.data(); }
  }


  /// @brief Returns non const C++ pointer to C array storing in vector data.
  /// - If vector is empty, returns NULL pointer.
  /// - If vector is not empty, returns pointer to a C array.
  /// @pre vector v is not const.
  template <typename T>
  T * vector2pointer(std::vector<T> & v)
  {
    if (v.empty()) { return(NULL); }
    else { return v.data(); }
  }

    
  /// @brief DEPRECATED: Returns non const C++ pointer to C array storing in vector data.
  /// - If vector is empty, returns NULL pointer.
  /// - If vector is not empty, returns pointer to a C array.
  /// @pre vector v is not const.
  template <typename T>
  T * vector2pointerNC(std::vector<T> & v)
  {
    if (v.empty()) { return(NULL); }
    else { return v.data(); }
  }


  /// Initialize a C++ vector to a specific size and value.
  template <typename T>
  void init_vector(std::vector<T> & v, const size_t num, const T value)
  {
    v.clear();
    if (num > 0) { v.resize(num, value); }
  }

  //@}


  // *****************************************************************
  //! @name PUSH MULTIPLE ELEMENTS ON A C++ VECTOR
  // *****************************************************************

  //@{

  /// @brief Push two elements on a C++ vector.
  template <typename T>
  inline void push_back(const T & a0, const T & a1, std::vector<T> & v)
  {
    v.push_back(a0);
    v.push_back(a1);
  }


  /*!
   *  @brief Push three elements on a C++ vector. DEPRECATED.
   *  - Deprecated. Replace by push_backIII().
   */
  template <typename T>
  inline void push_back
  (const T & a0, const T & a1, const T & a2, std::vector<T> & v)
  {
    v.push_back(a0);
    v.push_back(a1);
    v.push_back(a2);
  }

  /*!
   *  @overload
   *  @brief Push four elements on a C++ vector.
   */
  template <typename T>
  inline void push_back
  (const T & a0, const T & a1, const T & a2, const T & a3, std::vector<T> & v)
  {
    v.push_back(a0);
    v.push_back(a1);
    v.push_back(a2);
    v.push_back(a3);
  }


  /// @brief Push three elements on a C++ vector.
  template <typename T>
  inline void push_backIII
  (const T & a0, const T & a1, const T & a2, std::vector<T> & v)
  { 
    v.push_back(a0);
    v.push_back(a1);
    v.push_back(a2);
  }

  
  /// @brief Push array of three elements on a C++ vector.
  template <typename T>
  inline void push_backIII(const T a[3], std::vector<T> & v)
  { push_back(a[0], a[1], a[2], v); }
  
  //@}


  // *****************************************************************
  //! @name SET ARRAY
  // *****************************************************************

  //@{

  /// Set all elements of array a to x.
  template <typename NTYPE, typename ETYPE0, typename ETYPE1>
  void set_c_array(const NTYPE alength, const ETYPE0 x, ETYPE1 a[])
  {
    for (NTYPE i = 0; i < alength; i++)
      { a[i] = x; }
  }

  
  /// Set a[i] to x if flag[i] is true.
  template <typename NTYPE, typename ETYPE0, typename ETYPE1,
            typename FTYPE>
  void set_c_array(const NTYPE alength, const ETYPE0 x, 
                   const FTYPE & flag, ETYPE1 a[])
  {
    for (NTYPE i = 0; i < alength; i++)
      if (flag[i]) 
        { a[i] = x; }
  }

  
  /// Set a[i] to x if flag[i] is true.
  template <typename ETYPE0, typename ETYPE1, typename FTYPE>
  void set_array(const ETYPE0 x, const FTYPE & flag, std::vector<ETYPE1> & a)
  {
    if (a.empty()) { return; };
    set_c_array(a.size(), x, flag, &a.front());
  }

  
  /// Set array a[0..2] to {x0, x1, x2).
  template <typename ETYPE0, typename ETYPE1,
            typename ETYPE2, typename ETYPE3>
  void set_c_array(const ETYPE0 x0, const ETYPE1 x1, const ETYPE2 x2,
                   ETYPE3 a[])
  {
    a[0] = x0;
    a[1] = x1;
    a[2] = x2;
  }
  

  /// Set array a[0..3] to {x0, x1, x2, x3).
  template <typename ETYPE0, typename ETYPE1,
            typename ETYPE2, typename ETYPE3,
            typename ETYPE4>
  void set_c_array(const ETYPE0 x0, const ETYPE1 x1, const ETYPE2 x2,
                   const ETYPE3 x3, ETYPE4 a[])
  {
    a[0] = x0;
    a[1] = x1;
    a[2] = x2;
    a[3] = x3;
  }

  //@}


  // *****************************************************************
  //! @name GET MAX ARRAY ELEMENT
  // *****************************************************************

  //@{

  /*!
   *  @brief Return maximum value of elements in a C array.
   *  - Returns 0 if array has length 0.
   */
  template <typename T, typename NTYPE>
    T get_max_abs_array_value(const T a[], const NTYPE alength)
  {
    if (alength == 0)
      { return 0; }

    T amax = std::abs(a[0]);
    for (NTYPE i = 1; i < alength; i++) {
      const T x = std::abs(a[i]);
      if (x > amax)
	{ amax = x; }
    }

    return amax;
  }


  /*!
   *  @brief Return maximum absolute value of elements in an STL C++ vector.
   *  - Returns 0 if vector has length 0.
   */
  template <typename T>
  T get_max_abs_array_value(const std::vector<T> & v)
  {
    const T * a = IJK::vector2pointer(v);
    return get_max_abs_array_value(a, v.size());
  }

  //@}


  // *****************************************************************
  //! @name BITSET FUNCTIONS
  // *****************************************************************

  //@{

  /*!
   *  @brief return true if flag[i] is non-zero for some i in range [0..n-1].
   *  @param n Number of bits to check.
   *  @pre n <= flagA.size() and n <= flagB.size().
   */
  template <typename BIT_SET_TYPE, typename NTYPE>
  bool is_some_bit_true
  (const BIT_SET_TYPE & flag, const NTYPE n)
  {
    for (NTYPE i = 0; i < n; i++) {
      if (flag[i]) { return true; }
    }

    return false;
  }


  /*!
   *  @brief Return true if flagA[i] = flagB[i] for some i in range [0..n-1].
   *  @param n Number of bits to check.
   *  @pre n <= flagA.size() and n <= flagB.size().
   */
  template <typename BIT_SET_TYPEA, typename BIT_SET_TYPEB,
            typename NTYPE>
  bool does_some_true_bit_match
  (const BIT_SET_TYPEA & flagA,
   const BIT_SET_TYPEB & flagB,
   const NTYPE n)
  {
    for (NTYPE i = 0; i < n; i++) {
      if (flagA[i] && flagB[i])
	{ return true; }
    }

    return false;
  }

  //@}


  // *****************************************************************
  //! @name INTEGER POWER FUNCTIONS
  // *****************************************************************

  //@{


  /// @brief Integer power function.
  /// @param[out] result Equals (base)^p.
  template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
  void int_power(const ITYPE0 base, const ITYPE1 p, 
                 ITYPE2 & result)
  {
    result = 1;
    for (ITYPE1 k = 0; k < p; k++) {
      result *= base;
    }
  }


  /// @brief Integer power function with error checking.
  /// @param[out] result Equals (base)^p.
  template <typename ITYPE0, typename ITYPE1, typename ITYPE2,
            typename ETYPE>
  void int_power(const ITYPE0 base, const ITYPE1 p, 
                 ITYPE2 & result, ETYPE & error)
  {
    const ITYPE2 result_bound = (std::numeric_limits<ITYPE2>::max()/base);
    
    result = 1;
    for (ITYPE1 k = 0; k < p; k++) {

      if (result > result_bound) {
        error.AddMessage
          ("Result out of bounds. ", base, "^", p,
           " is larger than ", std::numeric_limits<ITYPE2>::max(), ".");
        throw error;
      }
      result *= base;
    }
  }


  /// @brief Integer power function. *** DEPRECATED ***
  /// - Return (base)^p.
  template <typename ITYPE0, typename ITYPE1>
  ITYPE0 int_power(const ITYPE0 base, const ITYPE1 p)
  {
    ITYPE0 result = 1;
    for (ITYPE1 k = 0; k < p; k++) {
      result *= base;
    }

    return(result);
  }

  //@}


  // *****************************************************************
  //! @name COUNT TEMPLATE FUNCTIONS
  // *****************************************************************

  //@{

  /// Count number of elements ge X in array a[].
  /// @tparam ATYPE Array type. Could be C array or stl vector.
  template <int X, typename ATYPE, typename NTYPE>
  NTYPE count_ge(const ATYPE a, const NTYPE length)
  {
    NTYPE num_ge = 0;
    for (NTYPE i = 0; i < length; i++) {
      if (a[i] >= X) { num_ge++; }
    }
    return(num_ge);
  }

  /// Count number of elements ge X in stl vector v.
  template <int X, typename T>
  typename std::vector<T>::size_type count_ge(const std::vector<T> & v)
  {
    return(count_ge<X>(v, v.size()));
  }

  //@}


  // *****************************************************************
  //! @name SELECT MIN TEMPLATE FUNCTIONS
  // *****************************************************************

  /*!
   *  @brief Set min_val to minimum value with flag_skip{i} false.
   *  - Set min_val to val0 if flag_skip0 is false and val0 <= val1.
   *  - Set min_val to val0 if flag_skip0 is false and flag_skip1 is true.
   *  - Set min_val to val1 if flag_skip1 is true and val0 > val1.
   *  - Set min_val to val1 if flag_skip1 is true and flag_skip0 is false.
   *  - Set flag_skip_all to true, if both flag_skip0 and flag_skip1 are true.
   */
  template <typename TYPE0, typename TYPE1, typename TYPE2, typename NTYPE>
  void select_min
  (const TYPE0 val0, const bool flag_skip0, 
   const TYPE1 val1, const bool flag_skip1, 
   TYPE2 & min_val, NTYPE & index_selected, bool & flag_skip_all)
  {
    if (flag_skip0) {
      if (flag_skip1) {
        flag_skip_all = true;
        index_selected = 0;
        min_val = val0;
      }
      else {
        flag_skip_all = false;
        index_selected = 1;
        min_val = val1;
      }
    }
    else {
      flag_skip_all = false;

      if (flag_skip1) {
        index_selected = 0;
        min_val = val0;
      }
      else if (val0 <= val1) {
        index_selected = 0;
        min_val = val0;
      }
      else {
        index_selected = 1;
        min_val = val1;
      }
    }
  }


  /// @brief Set min_val to minimum value with flag_skip{i} false.
  /// - Select from 3 values.
  template <typename TYPE0, typename TYPE1, typename TYPE2, 
            typename TYPE3, typename NTYPE>
  void select_minIII
  (const TYPE0 val0, const bool flag_skip0, 
   const TYPE1 val1, const bool flag_skip1, 
   const TYPE2 val2, const bool flag_skip2, 
   TYPE3 & min_val, NTYPE & index_selected, bool & flag_skip_all)
  {
    TYPE3 min_val01;
    bool flag_skip01;
    NTYPE index01_selected, indexB_selected;

    select_min(val0, flag_skip0, val1, flag_skip1,
               min_val01, index01_selected, flag_skip01);
    select_min(min_val01, flag_skip01, val2, flag_skip2,
               min_val, indexB_selected, flag_skip_all);

    if (indexB_selected == 0) 
      { index_selected = index01_selected; }
    else
      { index_selected = 2; }
  }


  /// @brief Set min_val to minimum value with flag_skip{i} false.
  /// - Select from 4 values.
  template <typename TYPE0, typename TYPE1, typename TYPE2, 
            typename TYPE3, typename TYPE4, typename NTYPE>
  void select_minIV
  (const TYPE0 val0, const bool flag_skip0, 
   const TYPE1 val1, const bool flag_skip1, 
   const TYPE2 val2, const bool flag_skip2, 
   const TYPE3 val3, const bool flag_skip3, 
   TYPE4 & min_val, NTYPE & index_selected, bool & flag_skip_all)
  {
    TYPE4 min_val01;
    TYPE4 min_val23;
    bool flag_skip01;
    bool flag_skip23;
    NTYPE index01_selected, index23_selected;
    NTYPE indexB_selected;

    select_min(val0, flag_skip0, val1, flag_skip1,
               min_val01, index01_selected, flag_skip01);
    select_min(val2, flag_skip2, val3, flag_skip3,
               min_val23, index23_selected, flag_skip23);
    select_min(min_val01, flag_skip01, min_val23, flag_skip23,
               min_val, indexB_selected, flag_skip_all);

    if (indexB_selected == 0) 
      { index_selected = index01_selected; }
    else 
      { index_selected = 2 + index23_selected; }
  }


  /// @brief Set min_val to minimum value with flag_skip{i} false.
  /// - Select from 5 values.
  template <typename TYPE0, typename TYPE1, typename TYPE2, 
            typename TYPE3, typename TYPE4, typename TYPE5,
            typename NTYPE>
  void select_minV
  (const TYPE0 val0, const bool flag_skip0, 
   const TYPE1 val1, const bool flag_skip1, 
   const TYPE2 val2, const bool flag_skip2, 
   const TYPE3 val3, const bool flag_skip3, 
   const TYPE4 val4, const bool flag_skip4, 
   TYPE5 & min_val, NTYPE & index_selected, bool & flag_skip_all)
  {
    TYPE5 min_val012;
    TYPE5 min_val34;
    bool flag_skip012;
    bool flag_skip34;
    NTYPE index012_selected, index34_selected;
    NTYPE indexB_selected;

    select_minIII(val0, flag_skip0, val1, flag_skip1, val2, flag_skip2,
                  min_val012, index012_selected, flag_skip012);
    select_min(val3, flag_skip3, val4, flag_skip4,
               min_val34, index34_selected, flag_skip34);
    select_min(min_val012, flag_skip012, min_val34, flag_skip34,
               min_val, indexB_selected, flag_skip_all);

    if (indexB_selected == 0) 
      { index_selected = index012_selected; }
    else 
      { index_selected = 3 + index34_selected; }
  }


  /// @brief Set min_val to minimum value with flag_skip{i} false.
  /// - Select from 6 values.
  template <typename TYPE0, typename TYPE1, typename TYPE2, 
            typename TYPE3, typename TYPE4, typename TYPE5,
            typename TYPE6, typename NTYPE>
  void select_minVI
  (const TYPE0 val0, const bool flag_skip0, 
   const TYPE1 val1, const bool flag_skip1, 
   const TYPE2 val2, const bool flag_skip2, 
   const TYPE3 val3, const bool flag_skip3, 
   const TYPE4 val4, const bool flag_skip4, 
   const TYPE5 val5, const bool flag_skip5, 
   TYPE6 & min_val, NTYPE & index_selected, bool & flag_skip_all)
  {
    TYPE5 min_val0123;
    TYPE5 min_val45;
    bool flag_skip0123;
    bool flag_skip45;
    NTYPE index0123_selected, index45_selected;
    NTYPE indexB_selected;

    select_minIV
      (val0, flag_skip0, val1, flag_skip1, val2, flag_skip2, val3, flag_skip3, 
       min_val0123, index0123_selected, flag_skip0123);
    select_min(val4, flag_skip4, val5, flag_skip5,
               min_val45, index45_selected, flag_skip45);
    select_min(min_val0123, flag_skip0123, min_val45, flag_skip45,
               min_val, indexB_selected, flag_skip_all);

    if (indexB_selected == 0) 
      { index_selected = index0123_selected; }
    else 
      { index_selected = 4 + index45_selected; }
  }


  /// @brief Set min_val to minimum value with flag_skip{i} false.
  /// - Select from 7 values.
  template <typename TYPE0, typename TYPE1, typename TYPE2, 
            typename TYPE3, typename TYPE4, typename TYPE5,
            typename TYPE6, typename TYPE7, typename NTYPE>
  void select_minVII
  (const TYPE0 val0, const bool flag_skip0, 
   const TYPE1 val1, const bool flag_skip1, 
   const TYPE2 val2, const bool flag_skip2, 
   const TYPE3 val3, const bool flag_skip3, 
   const TYPE4 val4, const bool flag_skip4, 
   const TYPE5 val5, const bool flag_skip5, 
   const TYPE6 val6, const bool flag_skip6, 
   TYPE7 & min_val, NTYPE & index_selected, bool & flag_skip_all)
  {
    TYPE6 min_val0123;
    TYPE6 min_val456;
    bool flag_skip0123;
    bool flag_skip456;
    NTYPE index0123_selected, index456_selected;
    NTYPE indexB_selected;

    select_minIV
      (val0, flag_skip0, val1, flag_skip1, val2, flag_skip2, val3, flag_skip3, 
       min_val0123, index0123_selected, flag_skip0123);
    select_minIII(val4, flag_skip4, val5, flag_skip5, val6, flag_skip6,
               min_val456, index456_selected, flag_skip456);
    select_min(min_val0123, flag_skip0123, min_val456, flag_skip456,
               min_val, indexB_selected, flag_skip_all);

    if (indexB_selected == 0) 
      { index_selected = index0123_selected; }
    else 
      { index_selected = 4 + index456_selected; }
  }


  /// @brief Set min_val to minimum value with flag_skip{i} false.
  /// - Select from 7 values.
  template <typename TYPE0, typename TYPE1, typename TYPE2, 
            typename TYPE3, typename TYPE4, typename TYPE5,
            typename TYPE6, typename TYPE7, typename TYPE8,
            typename NTYPE>
  void select_minVIII
  (const TYPE0 val0, const bool flag_skip0, 
   const TYPE1 val1, const bool flag_skip1, 
   const TYPE2 val2, const bool flag_skip2, 
   const TYPE3 val3, const bool flag_skip3, 
   const TYPE4 val4, const bool flag_skip4, 
   const TYPE5 val5, const bool flag_skip5, 
   const TYPE6 val6, const bool flag_skip6, 
   const TYPE7 val7, const bool flag_skip7, 
   TYPE8 & min_val, NTYPE & index_selected, bool & flag_skip_all)
  {
    TYPE6 min_val0123;
    TYPE6 min_val4567;
    bool flag_skip0123;
    bool flag_skip4567;
    NTYPE index0123_selected, index4567_selected;
    NTYPE indexB_selected;

    select_minIV
      (val0, flag_skip0, val1, flag_skip1, val2, flag_skip2, val3, flag_skip3, 
       min_val0123, index0123_selected, flag_skip0123);
    select_minIV
      (val4, flag_skip4, val5, flag_skip5, val6, flag_skip6, val4, flag_skip4,
       min_val4567, index4567_selected, flag_skip4567);
    select_min(min_val0123, flag_skip0123, min_val4567, flag_skip4567,
               min_val, indexB_selected, flag_skip_all);

    if (indexB_selected == 0) 
      { index_selected = index0123_selected; }
    else 
      { index_selected = 4 + index4567_selected; }
  }

  //@}


  // *****************************************************************
  // CLASS LINE_SEGMENT & ASSOCIATED FUNCTIONS
  // *****************************************************************


  /// @brief Line segment between two grid vertices
  template <typename VTYPE>
  class LINE_SEGMENT:public std::pair<VTYPE,VTYPE> {

  public:
    LINE_SEGMENT() {};  // constructor
    LINE_SEGMENT(const VTYPE iv0, const VTYPE iv1)
    { SetEnd(iv0, iv1); };

    // set functions

    /// Set line segment endpoints.
    void SetEnd(const VTYPE iv0, const VTYPE iv1)
    { this->first = iv0; this->second = iv1; Order(); };

    /// Order endpoints so that  V0() <= V1().
    void Order()
    { if (V0() > V1()) { std::swap(this->first, this->second); }; };


    // get functions
    VTYPE V0() const { return(this->first); };  ///< Return endpoint 0.
    VTYPE V1() const { return(this->second); };  ///< Return endpoint 1.
    template <typename ITYPE>
    VTYPE V(const ITYPE i) const       ///< Return endpoint i.
    {
      if (i == 0) { return(V0()); }
      else { return(V1()); }
    }

  };


  //! @name CLASS LINE_SEGMENT ASSOCIATED FUNCTIONS

  //@{

  /// Return true if (V0 < V1) for every pair (V0,V1)
  ///    of line segment endpoints
  template <typename VTYPE>
  bool is_ordered(const std::vector< LINE_SEGMENT<VTYPE> > & list)
  {
    typename std::vector< LINE_SEGMENT<VTYPE> > ::const_iterator pos;
    for (pos = list.begin(); pos != list.end(); ++pos) {
      if (pos->V0() > pos->V1())
        return(false);
    }
    return(true);
  }

  //@}


  // *****************************************************************
  //! @name CREATE OUTPUT FORMATTED STRINGS
  // *****************************************************************

  //@{

  /// @brief Compose string from two elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2> 
  std::string compose_string(T1 a1, T2 a2)
  {
    std::ostringstream os0;

    os0 << a1 << a2;
    return(os0.str());
  }

  
  /// @brief Compose string from three elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3;
    return os0.str();
  }

  
  /// @brief Compose string from four elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3, T4 a4)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4;
    return os0.str();
  }

  
  /// @brief Compose string from five elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5 > 
  std::string compose_string(T1 a1, T2 a2, T3 a3, T4 a4, T5 a5)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5;
    return os0.str();
  }

  
  /// @brief Compose string from seven elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6;
    return os0.str();
  }

  
  /// @brief Compose string from seven elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7;
    return os0.str();
  }

  
  /// @brief Compose string from eight elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8>
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8;
    return os0.str();
  }

  
  /// @brief Compose string from nine elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8, 
             typename T9 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8, T9 a9)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9;
    return os0.str();
  }

  
  /// @brief Compose string from ten elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8, 
             typename T9, typename T10 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8, T9 a9,
   T10 a10)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9 << a10;
    return os0.str();
  }

  
  /// @brief Compose string from ten elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8, 
             typename T9, typename T10, typename T11 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8, T9 a9,
   T10 a10, T11 a11)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9 
        << a10 << a11;
    return os0.str();
  }

  
  /// @brief Compose string from twelve elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8, 
             typename T9, typename T10, typename T11, typename T12 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8, T9 a9,
   T10 a10, T11 a11, T12 a12)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9 
        << a10 << a11 << a12;
    return os0.str();
  }

  
  /// @brief Compose string from twelve elements.
  /// - Format string using output string stream.
  template < typename T1, typename T2, typename T3, typename T4, 
             typename T5, typename T6, typename T7, typename T8, 
             typename T9, typename T10, typename T11, typename T12,
	     typename T13 > 
  std::string compose_string
  (T1 a1, T2 a2, T3 a3, T4 a4, T5 a5, T6 a6, T7 a7, T8 a8, T9 a9,
   T10 a10, T11 a11, T12 a12, T13 a13)
  {
    std::ostringstream os0;

    os0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9 
        << a10 << a11 << a12 << a13;
    return os0.str();
  }

  
  /// @brief Compose string from array.
  /// - Format string using output string stream.
  template <typename ATYPE, typename LTYPE>
  std::string compose_string_from_array
  (const char * prefix, const ATYPE a[], const LTYPE alength,
   const char * suffix)
  {
    std::ostringstream os0;

    os0 << prefix << "(";
    for (LTYPE i = 0; i < alength; i++) {
      if (i > 0) { os0 << ","; };
      os0 << a[i];
    }
    os0 << ")" << suffix;
    
    return os0.str();
  }


  /// @brief Compose string from array using given precision.
  /// - Format string using output string stream.
  template <typename ATYPE, typename LTYPE, typename ITYPE>
  std::string compose_string_from_array_with_precision
  (const char * prefix, const ATYPE a[], const LTYPE alength,
   const char * suffix, const ITYPE precision)
  {
    std::ostringstream os0;

    os0.precision(precision);
    
    os0 << prefix << "(";
    for (LTYPE i = 0; i < alength; i++) {
      if (i > 0) { os0 << ","; };
      os0 << a[i];
    }
    os0 << ")" << suffix;
    
    return os0.str();
  }
  
  //@}


  // *****************************************************************
  // ERROR CLASSES
  // *****************************************************************

  /// Error class.
  class ERROR {

  protected:

    /// Error messages.
    std::vector<std::string> msg;   

    /// Precision of output used in AddArrayMessage().
    SET_VALUE<int> precision;

    
  public:

    /// Constructor.
    ERROR(){};

    /// Constructor and (const char *) error message.
    ERROR(const char * error_msg) { AddMessage(error_msg); };

    /// Constructor and (std::string) error message.
    ERROR(const std::string & error_msg) { AddMessage(error_msg); };


    ERROR(const ERROR & error)                      // copy constructor
    {
      for (int i = 0; i < error.NumMessages(); i++)
        { msg.push_back(error.Message(i)); };
    };

    const ERROR & operator = (const ERROR & right)  // copy assignment
    {
      if (&right != this) {
        for (int i = 0; i < right.NumMessages(); i++)
          { msg.push_back(right.Message(i)); };
      }

      return(*this);
    };

    // get functions
    std::string Message(const int i) const
    {
      if (i >= 0 && i < NumMessages()) { return(msg[i]); }
      else { return(""); }
    }
    int NumMessages() const { return(msg.size()); };

    // get procedure error messages
    std::string ProcMessage(const std::string & procedure_name)
    {
      std::string error_msg = 
        "Error detected in procedure: " +  procedure_name + ".";
      return(error_msg);
    }
    std::string ProcMessage(const char * procedure_name)
    { return(ProcMessage(std::string(procedure_name))); }

    // Set functions

    /// Set precision of output in AddArrayMessage().
    void SetPrecision(const int precision)
    { this->precision.Set(precision); }

    /// Unset output precision. (Use default precision).
    void UnsetPrecision()
    { this->precision.Unset(); }
    
    /// Add error message.
    void AddMessage(const std::string & error_msg)
    { msg.push_back(error_msg); };

    /// Add error message.
    void AddMessage(const char * error_msg)     
    { AddMessage(std::string(error_msg)); };

    /// Add message with two elements.
    template <typename T1, typename T2>
    void AddMessage(T1 m1, T2 m2)
    { AddMessage(compose_string(m1, m2)); };

    /// Add message with three elements.
    template <typename T1, typename T2, typename T3>
    void AddMessage(T1 m1, T2 m2, T3 m3)
    { AddMessage(compose_string(m1, m2, m3)); };

    /// Add message with four elements.
    template <typename T1, typename T2, typename T3, typename T4>
    void AddMessage(T1 m1, T2 m2, T3 m3, T4 m4)
    { AddMessage(compose_string(m1, m2, m3, m4)); };

    /// Add message with five elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5>
    void AddMessage(T1 m1, T2 m2, T3 m3, T4 m4, T5 m5)
    { AddMessage(compose_string(m1, m2, m3, m4, m5)); };

    /// Add message with six elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6>
    void AddMessage(T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6)
    { AddMessage(compose_string(m1, m2, m3, m4, m5, m6)); };

    /// Add message with seven elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7)
    { AddMessage(compose_string(m1, m2, m3, m4, m5, m6, m7)); };

    /// Add message with eight elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8)
    { AddMessage
        (compose_string(m1, m2, m3, m4, m5, m6, m7, m8)); };

    /// Add message with nine elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8,
              typename T9>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8, T9 m9)
    { AddMessage
        (compose_string(m1, m2, m3, m4, m5, m6, m7, m8, m9)); };

    /// Add message with ten elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8,
              typename T9, typename T10>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8, T9 m9,
     T10 m10)
    { AddMessage
        (compose_string(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)); };

    /// Add message with eleven elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8,
              typename T9, typename T10, typename T11>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8, T9 m9,
     T10 m10,  T11 m11)
    { AddMessage
        (compose_string(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11)); };

    /// Add message with twelve elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8,
              typename T9, typename T10, typename T11, typename T12>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8, T9 m9,
     T10 m10,  T11 m11, T12 m12)
    { AddMessage
        (compose_string
	 (m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12)); };

    /// Add message with thirteen elements.
    template <typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename T7, typename T8,
              typename T9, typename T10, typename T11, typename T12,
	      typename T13>
    void AddMessage
    (T1 m1, T2 m2, T3 m3, T4 m4, T5 m5, T6 m6, T7 m7, T8 m8, T9 m9,
     T10 m10,  T11 m11, T12 m12, T13 m13)
    { AddMessage
        (compose_string
         (m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13)); };

    /// Add message with array.
    template <typename ATYPE, typename LTYPE>
    void AddArrayMessage
    (const char * prefix, const ATYPE a[], const LTYPE alength,
     const char * suffix)
    {
      if (precision.IsSet()) {
        AddMessage(compose_string_from_array_with_precision
                   (prefix, a, alength, suffix, precision.Value()));
      }
      else {
        AddMessage(compose_string_from_array(prefix, a, alength, suffix));
      }
    }

    /// Add message with array element
    template <typename ATYPE, typename DTYPE>
    void AddArrayElementMessage
    (const char * prefix, const char * array_name,
     const ATYPE a[], const DTYPE d,
     const char * suffix)
    {
      AddMessage(compose_string
                 (prefix, array_name, "[", d, "] = ",
                  a[d], suffix));
    }

    void AddProcMessage(const char * procname)
    { AddMessage(ProcMessage(procname)); };
    void AddProcMessage(const std::string & procname)
    { AddMessage(ProcMessage(procname)); };
    void SetMessage(const int i, const std::string & error_msg)
    {
      if (i >= 0 && i < NumMessages())
        msg[i] = error_msg;
    };
    void SetMessage(const int i, const char * error_msg)
    { SetMessage(i, std::string(error_msg)); };
    void SetProcMessage(const int i, const char * procname)
    { SetMessage(i, ProcMessage(procname)); };
    void SetProcMessage(const int i, const std::string & procname)
    { SetMessage(i, ProcMessage(procname)); };
    void ClearAll() { msg.clear(); };

    ERROR & operator ()(const std::string & msg1)
    {
      AddMessage(msg1);
      return(*this);
    };

    ERROR & operator ()(const std::string & msg1, const std::string & msg2)
    {
      AddMessage(msg1);
      AddMessage(msg2);
      return(*this);
    }

    // print error messages to ostream out
    void Print(std::ostream & out) const
    {
      for (int i = 0; i < NumMessages(); i++) {
        out << Message(i) << std::endl;
      }
    }

  };

  /// @brief Procedure error class.
  /// Includes procedure name as first error message.
  class PROCEDURE_ERROR:public ERROR {

  protected:
    std::string procedure_name;

  public:
    PROCEDURE_ERROR(const char * procedure_name)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
    };
    PROCEDURE_ERROR(const char * procedure_name, const char * error_msg)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
      AddMessage(error_msg);
    };
    PROCEDURE_ERROR(const std::string & procedure_name, 
                    const std::string & error_msg)
    {
      this->procedure_name = procedure_name;
      AddProcMessage(procedure_name);
      AddMessage(error_msg);
    };

  };


  // *****************************************************************
  //! @name CHECK FUNCTIONS
  // *****************************************************************

  //@{

  /// @brief Return true if array memory is allocated.
  /// - Return false and set error message if array is NULL.
  template <typename T>
  bool check_array_allocated
  (const T * array, const char * array_name, ERROR & error)
  {
    if (array == NULL) {
      error.AddMessage("Programming error. Memory for array ",
                       array_name, "[] not allocated.");
      return(false);
    }

    return(true);
  }


  /// @brief Return true if array is not empty (not NULL).
  /// - Return false and set error message if array is NULL.
  template <typename T>
  bool check_array_non_empty
  (const T * array, const char * array_name, ERROR & error)
  {
    if (array == NULL) {
      error.AddMessage("Programming error. Array ",
                       array_name, "[] is empty.");
      return(false);
    }

    return(true);
  }


  /// @brief Return true if vA.size() = vB.size().
  /// - Return false and set error message
  ///   if vA.size() != vB.size().
  template <typename TA, typename TB>
  bool check_equal_vector_sizes
  (const std::vector<TA> & vA,
   const std::vector<TB> & vB,
   const char * vA_name,
   const char * vB_name,
   ERROR & error)
  {
    if (vA.size() != vB.size()) {
      error.AddMessage
        ("Programming error. Unequal sizes of C++ STL vectors ",
         vA_name, " and ", vB_name, ".");
      error.AddMessage("  ", vA_name, ".size() = ", vA.size(), ".");
      error.AddMessage("  ", vB_name, ".size() = ", vB.size(), ".");
      return(false);
    }

    return(true);
  }


  /// Return true if ptr is NULL
  template <typename T>
  bool check_is_NULL
  (const T * ptr, const char * variable_name, ERROR & error)
  {
    if (ptr != NULL) {
      error.AddMessage
        ("Programming error.  Previously allocated memory for variable ",
         variable_name, " not released.");
      return(false);
    }

    return(true);
  }

  /// Return true if type T has at least num_bits bits.
  template <typename T, typename NTYPE>
  bool check_number_of_bits(const NTYPE num_bits)
  {
    T x = std::numeric_limits<T>::max();
   
    NTYPE n = num_bits;
    while (n > 0) {

      if (x <= 0) { return(false); }

      x = (x >> 1);
      n--;
    }

    return(true);
  }

  //@}


  // *****************************************************************
  // TEMPLATE CLASS BOX MEMBER FUNCTIONS
  // *****************************************************************

  template <typename T, int INIT_DIM> void BOX<T,INIT_DIM>::Init(const int dim)
  {
    dimension = dim;
    min_coord = NULL;
    max_coord = NULL;
    if (dim > 0) { SetDimension(dimension); }
  }

  template <typename T, int INIT_DIM> void BOX<T,INIT_DIM>::FreeAll()
  {
    dimension = 0;
    if (min_coord != NULL) { delete [] min_coord; }
    if (max_coord != NULL) { delete [] max_coord; }
    min_coord = NULL;
    max_coord = NULL;
  }

  template <typename T, int INIT_DIM> BOX<T,INIT_DIM>::
  BOX(const int dimension)
  { Init(dimension); }

  template <typename COORD_TYPE, int INIT_DIM>
  template <typename CTYPE2>
  bool BOX<COORD_TYPE,INIT_DIM>::Contains(const CTYPE2 * coord) const
  {
    for (int d = 0; d < dimension; d++) {
      if (coord[d] < min_coord[d]) { return(false); }
      if (coord[d] > max_coord[d]) { return(false); }
    }
    return(true);
  }

  template <typename COORD_TYPE, int INIT_DIM> 
  void BOX<COORD_TYPE,INIT_DIM>::SetDimension(const int d)
  {
    FreeAll();
    if (d <= 0) { return; }
    min_coord = new COORD_TYPE[d];
    max_coord = new COORD_TYPE[d];
    dimension = d;
  }

  template <typename COORD_TYPE, int INIT_DIM>
  template <typename CTYPE2>
  void BOX<COORD_TYPE,INIT_DIM>::
  SetMinCoord(const CTYPE2 * coord)
  {
    for (int d = 0; d < dimension; d++)
      SetMinCoord(d, coord[d]);
  }

  template <typename COORD_TYPE, int INIT_DIM>
  template <typename CTYPE2>
  void BOX<COORD_TYPE,INIT_DIM>::
  SetMaxCoord(const CTYPE2 * coord)
  {
    for (int d = 0; d < dimension; d++)
      SetMaxCoord(d, coord[d]);
  }

  template <typename COORD_TYPE, int INIT_DIM>
  template <typename CTYPE2, typename CTYPE3>
  void BOX<COORD_TYPE,INIT_DIM>::SetCoord
  (const CTYPE2 * minc, const CTYPE3 * maxc)
  {
    SetMinCoord(minc);
    SetMaxCoord(maxc);
  }

  template <typename COORD_TYPE, int INIT_DIM>
  void BOX<COORD_TYPE,INIT_DIM>::SetAllMinCoord(const COORD_TYPE c)
  {
    for (int d = 0; d < dimension; d++)
      SetMinCoord(d, c);
  }

  template <typename COORD_TYPE,int INIT_DIM>
  void BOX<COORD_TYPE,INIT_DIM>::SetAllMaxCoord(const COORD_TYPE c)
  {
    for (int d = 0; d < dimension; d++)
      SetMaxCoord(d, c);
  }

  // Extend the min/max coords so box includes coord.
  template <typename COORD_TYPE, int INIT_DIM>
  template <typename CTYPE2>
  void BOX<COORD_TYPE,INIT_DIM>::Extend(const CTYPE2 * coord)
  {
    for (int d = 0; d < dimension; d++) {
      if (coord[d] < MinCoord(d))
        { SetMinCoord(d, coord[d]); }
      if (coord[d] > MaxCoord(d))
        { SetMaxCoord(d, coord[d]); }
    }
  }

  // copy constructor template
  template <typename COORD_TYPE, int INIT_DIM> 
  BOX<COORD_TYPE,INIT_DIM>::BOX(const BOX<COORD_TYPE,INIT_DIM> & box)
  {
    Init(box.Dimension());
    for (int d = 0; d < box.Dimension(); d++) {
      SetMinCoord(d, box.MinCoord(d));
      SetMaxCoord(d, box.MaxCoord(d));
    }
  }

  // copy assigment template
  template <typename COORD_TYPE, int INIT_DIM> 
  const BOX<COORD_TYPE,INIT_DIM> & BOX<COORD_TYPE,INIT_DIM>::operator = 
  (const BOX<COORD_TYPE,INIT_DIM> & right)
  {
    if (&right != this) {         // avoid self-assignment
      FreeAll();
      SetDimension(right.Dimension());
      for (int d = 0; d < right.Dimension(); d++) {
        SetMinCoord(d, right.MinCoord(d));
        SetMaxCoord(d, right.MaxCoord(d));
      }
    }
    return *this;
  }

  // Print coord[]. (Extended version.)
  template <typename COORD_TYPE, int INIT_DIM> 
  template <typename OSTREAM_TYPE>
  void BOX<COORD_TYPE,INIT_DIM>::
  _PrintCoordX(OSTREAM_TYPE & out, const COORD_TYPE * coord,
              const char c0, const char c1, const char c2) const
  {
    out << c0;
    for (int d = 0; d < Dimension(); d++) {
      out << coord[d];
      if (d+1 < Dimension()) { out << c1; }
    }
    out << c2;
  }


  // Print coord[]. (Set prefix and suffix.)
  template <typename COORD_TYPE, int INIT_DIM> 
  template <typename OSTREAM_TYPE,
            typename STYPE0, typename STYPE1>
  void BOX<COORD_TYPE,INIT_DIM>::
  _PrintCoord
  (OSTREAM_TYPE & out, 
   const STYPE0 & s0, const COORD_TYPE * coord,
   const STYPE1 & s1) const
  {
    out << s0;
    _PrintCoord(out, coord);
    out << s1;
  }


  // copy constructor template
  template <typename COORD_TYPE, int INIT_DIM> 
  bool BOX<COORD_TYPE,INIT_DIM>::CheckCoord(ERROR & error) const
  {
    for (int d = 0; d < Dimension(); d++) {
      if (!IsMinLEMaxCoord(d)) {
        error.AddMessage
          ("Error.  Minimum coordinate[", d, 
           "] > Maximum coordinate[", d, "].");
        error.AddMessage
          ("  Minimum coordinate[", d, "] = ", MinCoord(d), ".");
        error.AddMessage
          ("  Maximum coordinate[", d, "] = ", MaxCoord(d), ".");

        return(false);
      }
    }

    return(true);
  }


  // *****************************************************************
  // TEMPLATE CLASS ARRAY MEMBER FUNCTIONS
  // *****************************************************************

  template <typename ETYPE> ARRAY<ETYPE>::~ARRAY()
  {
    Free();
  }

  template <typename ETYPE> void ARRAY<ETYPE>::Free()
  {
    delete [] element;
    element = NULL;
  }

}

#endif
