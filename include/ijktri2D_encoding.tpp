/*!
   @file ijktri2D_encoding.tpp
   @brief ijk templates for encoding triangulation of simple polygons.
   - Based on 1,0 Lukasiewicz encoding of convex polygon triangulations 
     described in: 
   - T. Koshy, Catalan Numbers with Applications, 2009, p. 236, and in
   - Saracevic, Masociv, Milosevic, JAVA implementation for
     Triangulation of Convex Polygon Based on Lukasiewicz's Algorithm
     and Binary Trees, Southeast Europe Journal on Soft Computing, 2013.
   - Version 0.4.0

   \par Encoding:
    - Construct the dual graph of the triangulation, including edges crossing
      the polygon edges. (One vertex for each triangle and one vertex outside
      the polygon for each polygon edge.)
    - Tree root is vertex corresponding to (v0, v_{n-1}).
    - Label root as 1.
    - Label each interior vertex 1.
    - Label each leaf as 0.
    - Traverse dual graph (binary tree), reporting vertex labels in "pre-order". 
      (Report label the first time vertex is reached.)
    - Reported labels form the encoding.
    - Reported labels form balanced parentheses (number of 1's always greater than number
      of 0's).
    - Encoding always starts with two 1's and ends with two 0's, so these are implicit and
      not actually represented in the data structure.
*/


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2021-2023 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJKTRI2D_ENCODING_
#define _IJKTRI2D_ENCODING_

#include <bitset>
#include <vector>

#include "ijk.tpp"


namespace IJK {

  // ***************************************************************
  // CLASS FOR POLYGON TRIANGULATION ENCODING
  // ***************************************************************

  template <int BIT_SET_SIZE>
  class POLYGON_TRIANGULATION_ENCODING {

  public:

    /// Return constatnt BIT_SET_SIZE.
    constexpr int BitSetSize() const { return(BIT_SET_SIZE); };

    typedef typename std::bitset<BIT_SET_SIZE> BIT_SET_TYPE;


  protected:

    /// @brief Encoding of a polygon triangulation.
    /// - Note: First two 1's and last two zeros are implicit.
    BIT_SET_TYPE triangulation_encoding;


  protected:

    /// @brief Private function to count actual number of ones in first n bits.
    /// @pre n <= Size().
    template <typename NTYPE>
    NTYPE _CountNumOnes_(const NTYPE n) const;


  public:
    POLYGON_TRIANGULATION_ENCODING() 
    { Clear(); };

    /// Return bit i.
    template <typename ITYPE>
    bool operator [] (const ITYPE i) const
    { return(triangulation_encoding[i]); }

    /// Return bit i.
    template <typename ITYPE>
    bool TriangulationEncoding(const ITYPE i) const
    { return(triangulation_encoding[i]); }

    /// Return bit set size.
    std::size_t Size() const
    { return(triangulation_encoding.size()); }

    /// @brief Return maximum number of polygon vertices supported in this encoding.
    /// - Note: First two 1's and last two zeros are implicit.
    std::size_t MaxNumVertices() const
    { return((Size()/2)+3); }

    /// @brief Number of 1's in triangulation encoding of polygon with n vertices.
    /// - Does not include two implicit 1's at the beginning of the encoding.
    template < typename NTYPE>
    NTYPE NumOnes(const NTYPE numv) const
    { return(numv > 3 ? numv-3: 0); }

    /// Count total number of 1's in the bitset triangulation_encoding.
    int CountNumOnes() const
    { return (triangulation_encoding.count()); }

    /// Return number of bits used to represent triangulation 
    ///   of polygon with numv vertices.
    template <typename NTYPE>
    NTYPE NumBitsUsed(const NTYPE numv) const
    { return(2*NumOnes(numv)); }

    /// Number of triangles represented by the encoding.
    template <typename NTYPE>
    NTYPE NumTriangles(const NTYPE numv) const
    { return(NumOnes(numv)+1); }

    /// Return true if all bits are zero.
    bool IsZero() const
    { return(triangulation_encoding.none()); }

    /*!
     *  @brief Return bit set indicating which vertices are ears.
     *  @param flag_ears[i] True if vertex i is an ear.
     *    - Note: triangulation_encoding[] refers to vertices 
     *      from left to right while flag_ears[] refers to vertices
     *      from right to left.
     *  @pre flag_ears.size() >= numv.
     */
    template <typename NTYPE, typename _BIT_SET_TYPE>
    void GetFlagEars(const NTYPE numv, _BIT_SET_TYPE & flag_ears) const;


    // Set routines.

    /// Set all bits in triangulation_encoding to zero.
    void Clear()
    { triangulation_encoding.reset(); }

    /// set bit it to flag.
    template <typename ITYPE>
    void Set(const ITYPE i, const bool flag)
    { triangulation_encoding[i] = flag; }

    /// @brief Set encoding to represent fan from vertex iv.
    /// @param iv Fan from vertex iv.  0 <= iv <= numv-1.
    template <typename ITYPE, typename NTYPE>
    void SetFan(const ITYPE iv, const NTYPE numv);

    /// Set to "first" encoding, "1111...000".
    template <typename NTYPE>
    void FirstEncoding(const NTYPE numv);

    /// @brief Set to "next" encoding.
    /// - Return false if no next encoding.
    template <typename NTYPE>
    bool NextEncoding(const NTYPE numv);

    /// @brief Print encoding.
    /// - Mainly used for debugging.
    template <typename OSTREAM_TYPE, typename NTYPE>
    void PrintEncodingX(OSTREAM_TYPE & out, const NTYPE numv) const;

    /*!
     *  @brief Print encoding.
     *  - Mainly used for debugging.
     *  - Version that computes number of polygon vertices from number of ones
     *    in triangulation_encoding;
     */
    template <typename OSTREAM_TYPE>
    void PrintEncodingX(OSTREAM_TYPE & out) const;

    /*!
     *  @brief Print encoding.
     *  - Mainly used for debugging.
     *  - Version adding prefix and suffix strings.
     */
    template <typename OSTREAM_TYPE, typename NTYPE,
              typename STYPE0, typename STYPE1>
    void PrintEncoding(OSTREAM_TYPE & out, const NTYPE numv,
                       const STYPE0 & s0, const STYPE1 & s1) const;

    /*!
     *  @brief Print encoding.
     *  - Mainly used for debugging.
     *  - Version adding prefix and suffix strings.
     *  - Version that computes number of polygon vertices from number of ones
     *    in triangulation_encoding;
     */
    template <typename OSTREAM_TYPE,
              typename STYPE0, typename STYPE1>
    void PrintEncoding(OSTREAM_TYPE & out,
                       const STYPE0 & s0, const STYPE1 & s1) const;

    /// @brief Print triangulation representation.
    /// - Mainly used for debugging.
    template <typename OSTREAM_TYPE, typename NTYPE>
    void PrintTriangulationRepresentationX
    (OSTREAM_TYPE & out, const NTYPE numv) const;

    /*!
     *  @brief Print triangulation representation.
     *  - Mainly used for debugging.
     *  - Version that computes number of polygon vertices from number of ones
     *    in triangulation_encoding;
     */
    template <typename OSTREAM_TYPE>
    void PrintTriangulationRepresentationX
    (OSTREAM_TYPE & out) const;

    /*!
     *  @brief Print triangulation representation.
     *  - Mainly used for debugging.
     *  - Version adding prefix and suffix strings.
     */
    template <typename OSTREAM_TYPE, typename NTYPE,
              typename STYPE0, typename STYPE1>
    void PrintTriangulationRepresentation
    (OSTREAM_TYPE & out, const NTYPE numv,
     const STYPE0 & s0, const STYPE1 & s1) const;

    /*!
     *  @brief Print triangulation representation.
     *  - Mainly used for debugging.
     *  - Version that computes number of polygon vertices from number of ones
     *    in triangulation_encoding;
     *  - Version adding prefix and suffix strings.
     */
    template <typename OSTREAM_TYPE,
              typename STYPE0, typename STYPE1>
    void PrintTriangulationRepresentation
    (OSTREAM_TYPE & out,
     const STYPE0 & s0, const STYPE1 & s1) const;


    // Check routines.

    /// Check that (1,triangulation_encoding,0) is balanced parentheses
    ///   with NumOnes(numv) 1's.
    /// - Calls also CheckNumOnes() and CheckNumVertices().
    template <typename NTYPE>
    bool Check(const NTYPE numv, IJK::ERROR & error) const;

    /// Check that BIT_SET_SIZE is large enough to handle the given
    ///   number of vertices.
    template <typename NTYPE>
    bool CheckNumVertices(const NTYPE numv, IJK::ERROR & error) const;

    /// Check that triangulation_encoding has NumOnes(numv) 1's.
    template <typename NTYPE>
    bool CheckNumOnes(const NTYPE numv, IJK::ERROR & error) const;

  };


  // ***************************************************************
  // MEMBER FUNCTIONS FOR POLYGON_TRIANGULATIN_ENCODING.
  // ***************************************************************


  // Return bit set indicating which vertices are ears.
  // @param flag_ears[i] True if vertex i is an ear.
  template <int BIT_SET_SIZE>
  template <typename NTYPE, typename _BIT_SET_TYPE>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  GetFlagEars(const NTYPE numv, _BIT_SET_TYPE & flag_ears) const
  {
    const int THREE(3);

    flag_ears.reset();

    if (numv <= THREE) {
      // All vertices are ears.
      flag_ears.set();
      return;
    }

    // Read encoding backwards, starting at rightmost label 0.

    bool bitC = false;
    bool bitB = false;

    // Initialize vertex (numv-1) to ear.
    flag_ears[numv-1] = true;

    NTYPE iv = numv-2;
    NTYPE num_zeros = 2;
    NTYPE num_ones = 0;
    for (NTYPE i = 0; i < NumBitsUsed(numv); i++) {
      const bool bitA = TriangulationEncoding(i);
      if (bitA) {
	if (!bitB && !bitC)
	  { flag_ears[iv] = true; }
	num_ones++;
      }
      else {
	iv--;
	num_zeros++;
      }

      if ((num_ones+1) == num_zeros) {
	// Triangle incident on (numv-1,0) has right child triangle.
	// Triangle (numv-1,0,1) is not an ear.
	flag_ears[numv-1] = false;
      }

      bitC = bitB;
      bitB = bitA;
    }

    const NTYPE ifirst_bit = NumBitsUsed(numv)-1;
    if (!TriangulationEncoding(ifirst_bit))
      { flag_ears[0] = true; }

    if (iv != 1) {
      IJK::PROCEDURE_ERROR error("POLYGON_TRIANGULATION_ENCODING::GetFlagEars");

      error.AddMessage("Programming error. Incorrect indexing of vertices.");
      error.AddMessage("  iv = ", iv, "  Expected value of iv is 1.");
      throw error;
    }
  }


  // Set encoding to represent fan from vertex iv.
  // @param iv Fan from vertex iv.  0 <= iv <= numv-1.
  template <int BIT_SET_SIZE>
  template <typename ITYPE, typename NTYPE>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  SetFan(const ITYPE iv, const NTYPE numv)
  {
    Clear();

    if (numv < 4) { 
      // Encoding has zero bits.
      return; 
    }

    // Location of first (leftmost bit) representing triangulation.
    const ITYPE ifirst = NumBitsUsed(numv)-1;

    if (iv == 0) {
      for (int i = 0; i < NumOnes(numv); i++)
        { Set(ifirst-i, true); }
      return;
    }
    else if (iv+1 == numv) {
      for (int i = 0; i < NumOnes(numv); i++)
        { Set(2*i, true); }
    }
    else {
      IJK::PROCEDURE_ERROR error("POLYGON_TRIANGULATION_ENCODING::SetFan");

      if (iv > numv || iv < 0) {
        error.AddMessage("Programming error.  Illegal vertex index: ", iv, ".");
        error.AddMessage("  Number of vertices: ", numv, ".");
        error.AddMessage("  Vertex indices must be in range [0:", numv-1, "].");
        throw error;
      }

      for (ITYPE i = 0; i+1 < iv; i++) 
        { Set(ifirst-2*i, true); }

      for (ITYPE i = 0; i+2 < numv-iv; i++) {
        // Note: This only executes if numv-iv >= 3.
        Set(i+numv-iv-3, true); 
      }
    }
  }


  // Set to "first" encoding, "1111...000".
  template <int BIT_SET_SIZE>
  template <typename NTYPE>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  FirstEncoding(const NTYPE numv)
  {
    const NTYPE num1 = NumOnes(numv);

    Clear();

    for (NTYPE i = num1; i < 2*num1; i++)
      { triangulation_encoding[i] = 1; }
  }


  // Set to "next" encoding.
  // - Return false if no next encoding.
  template <int BIT_SET_SIZE>
  template <typename NTYPE>
  bool POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  NextEncoding(const NTYPE numv)
  {
    const NTYPE num1 = NumOnes(numv);
    NTYPE count0 = 0;
    NTYPE count1 = 0;

    for (NTYPE i = 0; i < NumBitsUsed(numv); i++) {

      if (triangulation_encoding[i] == 0) 
        { count0++; }
      else {
        count1++;

        if (count1 <= count0) {
          const NTYPE count0B = count0-1;

          triangulation_encoding[i] = 0;

          // Set [0..(i-1)] to be 111...000.
          for (NTYPE j = 0; j < count0B; j++) 
            { triangulation_encoding[j] = 0; }
          for (NTYPE j = count0B; j < count0B+count1; j++) 
            { triangulation_encoding[j] = 1; }

          return(true);
        }
      }
    }

    return(false);
  }


  // Print encoding.
  // - Mainly used for debugging.
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE, typename NTYPE>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintEncodingX(OSTREAM_TYPE & out, const NTYPE numv) const
  {
    if (numv < 4)
      { return; }

    // Location of first (leftmost bit) representing triangulation.
    const NTYPE ifirst = NumBitsUsed(numv)-1;

    for (NTYPE i = 0; i < NumBitsUsed(numv); i++) 
      { out << int(triangulation_encoding[ifirst-i]); }
  }


  // Print encoding.
  // - Mainly used for debugging.
  // - Version that computes number of polygon vertices from number of ones
  //   in triangulation_encoding;
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintEncodingX(OSTREAM_TYPE & out) const
  {
    const int num_ones = CountNumOnes();
    const int numv = num_ones + 3;

    PrintEncodingX(out, numv);
  }


  // Print encoding.
  // - Mainly used for debugging.
  // - Version adding prefix and suffix strings.
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE, typename NTYPE,
            typename STYPE0, typename STYPE1>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintEncoding(OSTREAM_TYPE & out, const NTYPE numv,
                const STYPE0 & s0, const STYPE1 & s1) const
  {
    out << s0;
    PrintEncodingX(out, numv);
    out << s1;
  }


  // Print encoding.
  // - Mainly used for debugging.
  // - Version that computes number of polygon vertices from number of ones
  //   in triangulation_encoding;
  // - Version adding prefix and suffix strings.
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE,
            typename STYPE0, typename STYPE1>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintEncoding(OSTREAM_TYPE & out,
                const STYPE0 & s0, const STYPE1 & s1) const
  {
    out << s0;
    PrintEncodingX(out);
    out << s1;
  }


  // Print triangulation representation.
  // - Mainly used for debugging.
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE, typename NTYPE>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintTriangulationRepresentationX
  (OSTREAM_TYPE & out, const NTYPE numv) const
  {
    out << "11";
    PrintEncodingX(out, numv);
    out << "00";
  }


  // Print triangulation representation.
  // - Mainly used for debugging.
  // - Version that computes number of polygon vertices from number of ones
  //   in triangulation_encoding;
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintTriangulationRepresentationX(OSTREAM_TYPE & out) const
  {
    out << "11";
    PrintEncodingX(out);
    out << "00";
  }


  // Print triangulation representation.
  // - Mainly used for debugging.
  // - Version adding prefix and suffix strings.
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE, typename NTYPE,
            typename STYPE0, typename STYPE1>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintTriangulationRepresentation
  (OSTREAM_TYPE & out, const NTYPE numv,
   const STYPE0 & s0, const STYPE1 & s1) const
  {
    out << s0;
    PrintTriangulationRepresentationX(out, numv);
    out << s1;
  }


  // Print triangulation representation.
  // - Mainly used for debugging.
  // - Version adding prefix and suffix strings.
  // - Version that computes number of polygon vertices from number of ones
  //   in triangulation_encoding
  template <int BIT_SET_SIZE>
  template <typename OSTREAM_TYPE,
            typename STYPE0, typename STYPE1>
  void POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  PrintTriangulationRepresentation
  (OSTREAM_TYPE & out,
   const STYPE0 & s0, const STYPE1 & s1) const
  {
    out << s0;
    PrintTriangulationRepresentationX(out);
    out << s1;
  }


  // Check that (1,triangulation_encoding,0) is balanced parentheses
  //   with NumOnes(numv) 1's.
  template <int BIT_SET_SIZE>
  template <typename NTYPE>
  bool POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  Check(const NTYPE numv, IJK::ERROR & error) const
  {
    const NTYPE num1 = NumOnes(numv);
    NTYPE count0 = 0;
    NTYPE count1 = 0;

    if (!CheckNumVertices(numv, error)) { return false; }
    if (!CheckNumOnes(numv, error)) { return false; };

    for (NTYPE i = 0; i < NumBitsUsed(numv); i++) {
      if (triangulation_encoding[i] == 0) 
        { count0++; }
      else {
        count1++;

        if (count1 > count0+1) {
          error.AddMessage
            ("Programming error. Unbalanced parentheses.");
          error.AddMessage
            ("  Suffix of last ", i, " bits contains ", count1, 
             " 1's but only ", count0, " 0's.");
          error.AddMessage
            ("  Number of 1's in any suffix should be at most number of 0's plus 1.");
          return false;
        }
      }
    }


    return true;
  }


  // Check that BIT_SET_SIZE is large enough to handle the given
  //   number of vertices.
  template <int BIT_SET_SIZE>
  template <typename NTYPE>
  bool POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  CheckNumVertices(const NTYPE numv, IJK::ERROR & error) const
  {
    if (numv > MaxNumVertices()) {
      error.AddMessage
        ("Programming error. Too many vertices.");
      error.AddMessage
        ("  Compiled instance of triangulation encoding can only support ",
         MaxNumVertices(), " vertices.");
      error.AddMessage
        ("  For polygons with more vertices, increase BIT_SET_SIZE");
      error.AddMessage
        ("  in TRIANGULATION_ENCODING template parameter and recompile.");
      return(false);
    }

    return true;
  }


  // Check that triangulation_encoding has NumOnes(numv) 1's.
  template <int BIT_SET_SIZE>
  template <typename NTYPE>
  bool POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  CheckNumOnes(const NTYPE numv, IJK::ERROR & error) const
  {
    const NTYPE count1 = _CountNumOnes_(NumBitsUsed(numv));
    const NTYPE num1 = NumOnes(numv);

    if (count1 != num1) {
      error.AddMessage("Programming error. Incorrect number of 1's in triangulation encoding.");
      error.AddMessage("  Number of 1's: ", count1, ".");
      error.AddMessage("  Expected number of 1's: ", num1, ".");
      return(false);
    }

    return(true);
  }


  /// Private function to count actual number of ones in first n bits.
  template <int BIT_SET_SIZE>
  template <typename NTYPE>
  NTYPE POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::
  _CountNumOnes_(const NTYPE n) const
  {
    NTYPE count1 = 0;

    for (NTYPE i = 0; i < n; i++) {
      if (triangulation_encoding[i] == 1) 
        { count1++; }
    }

    return(count1);
  }


  // ***************************************************************
  // TRIANGULATE POLYGON BASED ON TRIANGULATION ENCODING
  // ***************************************************************

  /// Generate polygon triangulation based on encoding.
  /// - Polygon vertices are numbered 0 to numv-1 in clockwise order around polygon.
  /// @param[in] gen_triangle Function to generate triangle.
  ///   - Takes 3 parameters, 3 triangle vertices.
  template <typename NTYPE, typename TRI_ENCODING, typename GEN_TRIANGLE>
  void gen_polygon_triangulation_based_on_encoding
  (const NTYPE numv, const TRI_ENCODING & tri_encoding, GEN_TRIANGLE gen_triangle)
  {
    std::vector<NTYPE> vertex_stack(numv);
    NTYPE i, istack, iv;
    IJK::PROCEDURE_ERROR error("TriangulatePolygonBasedOnEncoding");

    if (numv <= 3) {
      if (numv <= 2) {
        // No triangles.
        return; 
      }
      else {
        // One triangle, equivalent to original polygon.
        gen_triangle(0, 1, 2);
        return;
      }
    }

    // Read encoding backwards, starting at rightmost label 0.
    vertex_stack[0] = numv-1;
    vertex_stack[1] = numv-2;
    istack = 1;

    i = 0;
    iv = numv-3;
    for (i = 0; i < 2*(numv-3); i++) {
      if (tri_encoding[i]) {

        if (istack < 1) {

          if (!tri_encoding.CheckNumOnes(numv, error)) {
            // If total number of 1's does not equal number of 0's, report that error.
            throw error;
          }

          error.AddMessage("Error. Triangulation encoding is illegal.");
          error.AddMessage("  Encoding does not represent a balanced 1,0 parenthesization.");
          error.AddMessage("  Suffix [0:", i, "] contains more 1's than 0's.");
          throw error;
        }

        gen_triangle
          (iv, vertex_stack[istack], vertex_stack[istack-1]);
        istack--;
      }
      else {
        istack++;
        vertex_stack[istack] = iv;

        if (iv < 1) {

          if (!tri_encoding.CheckNumOnes(numv, error)) {
            // If total number of 1's does not equal number of 0's, report that error.
            throw error;
          }

          error.AddMessage("Error. Triangulation encoding is illegal.");
          error.AddMessage("  Encoding does not represent a balanced 1,0 parenthesization.");
          error.AddMessage("  Prefix [", i, ":", numv-1, "] contains more 0's than 1's.");
          throw error;
        }

        iv--;
      }
    }

    // Check stack size.
    if (istack < 1) {
      error.AddMessage("Error. Triangulation encoding is illegal.");
      error.AddMessage("  Encoding does not represent a balanced 1,0 parenthesization.");
      error.AddMessage("  Encoding contains more 1's than 0's.");
      throw error;
    }

    if (istack > 1) {
      error.AddMessage("Error. Triangulation encoding is illegal.");
      error.AddMessage("  Encoding does not represent a balanced 1,0 parenthesization.");
      error.AddMessage("  Encoding contains more 0's than 1's.");
      throw error;
    }

    gen_triangle(0, vertex_stack[1], vertex_stack[0]);
  }


  template <typename ITYPE>
  class ADD_TO_TRIANGLE_LIST_TRI2D_ENCODING {

  protected:
    std::vector<ITYPE> * tri_list_ptr;

  public:
    ADD_TO_TRIANGLE_LIST_TRI2D_ENCODING(std::vector<ITYPE> * tri_list_ptrX)
    { tri_list_ptr = tri_list_ptrX; };

  public:
    // Functor definition
    template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
    void operator () 
    (const ITYPE0 iv0, ITYPE1 iv1, const ITYPE2 iv2) const
    {
      tri_list_ptr->push_back(iv0);
      tri_list_ptr->push_back(iv1);
      tri_list_ptr->push_back(iv2);
    };

  };


  /// Generate list of triangles in polygon triangulation based on encoding.
  /// - Polygon vertices are numbered 0 to numv-1 in clockwise order around polygon.
  /// @param[out] Create a list of triangles.
  template <typename NTYPE, typename TRI_ENCODING, typename ITYPE>
  void gen_polygon_triangulation_list_based_on_encoding
  (const NTYPE numv, const TRI_ENCODING & tri_encoding,
   std::vector<ITYPE> & tri_list)
  {
    tri_list.clear();

    gen_polygon_triangulation_based_on_encoding
      (numv, tri_encoding, 
       ADD_TO_TRIANGLE_LIST_TRI2D_ENCODING<NTYPE>(&(tri_list)));
  }


  // ***************************************************************
  // PRINT TRIANGULATION. (MAINLY FOR DEBUGGING.)
  // ***************************************************************

  /// Print vertices of a triangle.
  /// - Mainly for debugging.
  /// @param c0 Left delimiter.
  /// @param c1 Separator.
  /// @param c2 Right delimiter.
  template <typename OSTREAM_TYPE, 
            typename ITYPE0, typename ITYPE1, typename ITYPE2>
  void print_triangle_verticesX_tri2D_encoding
  (OSTREAM_TYPE & out,
   const ITYPE0 iv0, const ITYPE1 iv1, const ITYPE2 iv2,
   const char c0, const char c1, const char c2)
  {
    out << c0 << iv0 << c1 << iv1 << c1 << iv2 << c2;
  }

  /// Print vertices of a triangle.
  /// - Mainly for debugging.
  /// - Version using delimiters '(' and ')' and separator ','.
  template <typename OSTREAM_TYPE, 
            typename ITYPE0, typename ITYPE1, typename ITYPE2>
  void print_triangle_vertices_tri2D_encoding
  (OSTREAM_TYPE & out,
   const ITYPE0 iv0, const ITYPE1 iv1, const ITYPE2 iv2)
  {
    print_triangle_verticesX_tri2D_encoding(out, iv0, iv1, iv2, '(', ',', ')');
  }

  /// Print vertices of a triangle.
  /// - Mainly for debugging.
  /// - Version adding prefix and suffix strings.
  template <typename OSTREAM_TYPE, 
            typename ITYPE0, typename ITYPE1, typename ITYPE2,
            typename STYPE0, typename STYPE1>
  void print_triangle_vertices_tri2D_encoding
  (OSTREAM_TYPE & out, const STYPE0 & s0, 
   const ITYPE0 iv0, const ITYPE1 iv1, const ITYPE2 iv2, 
   const STYPE1 & s1)
  {
    out << s0;
    print_triangle_vertices_tri2D_encoding(out, iv0, iv1, iv2);
    out << s1;
  }

  // Functor class for printing triangle vertices.
  template <typename OSTREAM_TYPE, typename STYPE0, typename STYPE1>
  class PRINT_TRIANGLE_VERTICES_TRI2D_ENCODING {

  protected:
    OSTREAM_TYPE * out;
    std::string s0;
    std::string s1;

  public:
    PRINT_TRIANGLE_VERTICES_TRI2D_ENCODING
    (OSTREAM_TYPE & outX, const STYPE0 & s0X, const STYPE1 & s1X)
    {
      out = &(outX);
      s0 = s0X;
      s1 = s1X;
    };

    // Functor definition
    template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
    void operator () 
    (const ITYPE0 iv0, const ITYPE1 iv1, const ITYPE2 iv2) const
    {
      print_triangle_vertices_tri2D_encoding
        (*out, s0, iv0, iv1, iv2, s1);
    };

  };


  /// Print polygon triangulation based on encoding.
  /// - Polygon vertices are numbered 0 to numv-1.
  /// @pre numv <= TRI_ENCODING::MaxNumVertices().
  /// @param[in] gen_triangle Function to generate triangle.
  ///   - Takes 3 parameters, 3 triangle vertices.
  template <typename OSTREAM_TYPE, 
            typename STYPE0, typename STYPE1, typename STYPE2,
            typename NTYPE, typename TRI_ENCODING>
  void print_polygon_triangulation_based_on_encoding
  (OSTREAM_TYPE & out, 
   const STYPE0 & header_prefix,
   const STYPE1 & header_suffix,
   const STYPE2 & line_prefix,
   const NTYPE numv, const TRI_ENCODING & tri_encoding)
  {
    PRINT_TRIANGLE_VERTICES_TRI2D_ENCODING<OSTREAM_TYPE, STYPE2, const char *>
      print_triangle_vertices(out, line_prefix, "\n");

    tri_encoding.PrintTriangulationRepresentation
      (out, numv, header_prefix, header_suffix);
    out << "\n";
    gen_polygon_triangulation_based_on_encoding
      (numv, tri_encoding, print_triangle_vertices);
  };

}

#endif
