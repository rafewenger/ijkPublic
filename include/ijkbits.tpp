/*!
 *  @file ijkbits.tpp
 *  @brief ijk templates for bit operations.
 *  - Version 0.4.0
 */

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2023 Rephael Wenger

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

#ifndef _IJK_BITS_
#define _IJK_BITS_

// *** DEBUG ***
#include <iostream>

namespace IJK {

  // **************************************************
  // BIT OPERATIONS
  // **************************************************

  /// Convert to an array of digits with the given base.
  template <typename ITYPE, typename BTYPE, typename DTYPE,
            typename NTYPE, typename ETYPE>
  void convert2base(const ITYPE ival, const BTYPE base, 
                    DTYPE * digit, const NTYPE max_num_digits,
                    ETYPE & error)
  {
    ITYPE jval = ival;
    for (NTYPE i = 0; i < max_num_digits; i++) {
      digit[i] = jval % base;
      jval = jval/base;
    }

    if (jval != 0) {
      error.AddMessage("Error converting ", ival, " to base ", base, ".");
      error.AddMessage("Output has more than ", max_num_digits, " digits.");

      throw error;
    };
  }


  /*!
   *  @brief Convert to a string of 0's and 1's.
   *  @param num_bits Number of 0's and 1's in output string.
   *     - If ival has more than num_bits bits, high order bits are ignored.
   */
  template <typename ITYPE, typename NTYPE, typename STYPE>
  void convert2bit_string(const ITYPE ival,
                          const NTYPE num_bits,
                          STYPE & bit_string)
  {
    // Initialize.
    bit_string = "";

    ITYPE jval = ival;
    for (NTYPE i = 0; i < num_bits; i++) {
      const ITYPE bit = jval % 2;
      if (bit == 0) { bit_string = STYPE("0") + bit_string; }
      else { bit_string = STYPE("1") + bit_string; }
      jval = jval/2;
    }
  }


  /// @brief Return true if ival[ipos] equals 1.
  /// @pre ipos < num_bits of ival.
  template <typename ITYPE, typename NTYPE>
  bool is_bit_one(const ITYPE ival, const NTYPE ipos)
  {
    const ITYPE mask = (ITYPE(1) << ipos);
    if ((ival & mask) == 0)
      { return false; }
    else
      { return true; }
  }


  /// Count number of zero and one bits in first num_bits of ival.
  template <typename ITYPE, typename NTYPE, typename NTYPE2, typename NTYPE3>
  void count_bits(const ITYPE ival, const NTYPE num_bits,
                  NTYPE2 & num_zeros, NTYPE3 & num_ones)
  {
    const ITYPE base = 2;

    num_zeros = 0;
    num_ones = 0;

    ITYPE jval = ival;
    for (NTYPE i = 0; i < num_bits; i++) {
      ITYPE bit = jval % base;
      if (bit == 0) { num_zeros++; }
      else { num_ones++; }
      jval = jval/base;
    }
  }


  /// Count number of one bits in first num_bits of ival.
  template <typename ITYPE, typename NTYPE, typename NTYPE2>
  void count_one_bits
  (const ITYPE ival, const NTYPE num_bits, NTYPE2 & num_ones)
  {
    const ITYPE base = 2;

    num_ones = 0;

    ITYPE jval = ival;
    for (NTYPE i = 0; i < num_bits; i++) {
      ITYPE bit = jval % base;
      if (bit == 1) { num_ones++; }
      jval = jval/base;
    }
  }


  /// Count bits with 1 in mask.
  template <typename ITYPE, typename MTYPE, typename NTYPE, 
            typename NTYPE2, typename NTYPE3>
  void count_masked_bits
  (const ITYPE ival, const MTYPE mask, const NTYPE num_bits,
   NTYPE2 & num_zeros, NTYPE3 & num_ones)
  {
    const ITYPE base = 2;

    num_zeros = 0;
    num_ones = 0;

    ITYPE jval = ival;
    MTYPE jmask = mask;
    for (NTYPE i = 0; i < num_bits; i++) {
      ITYPE bit = jval % base;
      MTYPE mask_bit = jmask % base;
      if (mask_bit == 1) {
        if (bit == 0) { num_zeros++; }
        else { num_ones++; }
      }
      jval = jval/base;
      jmask = jmask/base;
    }
  }


  /// Reverse order of bits in ival.
  template <typename ITYPE, typename NTYPE>
  ITYPE reverse_bits(const ITYPE ival, const NTYPE num_bits)
  {
    const ITYPE base = 2;

    ITYPE reverse_val = 0;
    ITYPE jval = ival;
    for (NTYPE i = 0; i < num_bits; i++) {
      ITYPE bit = jval % base;
      reverse_val = base*reverse_val;
      reverse_val = (reverse_val | bit);
      jval = jval/base;
    }

    return(reverse_val);
  }


  /// @brief Copy ival0[i0] equal to ival1[i1].
  /// @pre i0 < num_bits of ival0 and i1 < num_bits of ival1.
  template <typename ITYPE, typename BITS_TYPE0, typename BITS_TYPE1>
  void copy_bit
    (const ITYPE i0, const BITS_TYPE0 & ival0,
     const ITYPE i1, BITS_TYPE1 & ival1)
  {
    const bool flag_bit0 = is_bit_one(ival0, i0);
    const bool flag_bit1 = is_bit_one(ival1, i1);

    if (flag_bit0 == flag_bit1) {
      // ival0[i0] == ival1[i1]. Do nothing.
      return;
    }

    // Flip ival[i1].
    const BITS_TYPE1 mask1 = (BITS_TYPE1(1) << i1);
    ival1 = (ival1 ^ mask1);
  }


  /// Complement first num_bits in ival.
  template <typename ITYPE, typename NTYPE>
  ITYPE complement_bits(const ITYPE ival, const NTYPE num_bits)
  {
    const ITYPE all_ones = (~(ITYPE(0)));
    const ITYPE ival_complement = 
      (~ival) - ((all_ones >> num_bits) << num_bits);

    // Note: (2^{ival}-1 - ival) may not work if num_bits equals number
    //   of bits in type ITYPE.

    return(ival_complement);
  }


  /// Return true if (ival == reverse_bits(ival,num_bits))
  template <typename ITYPE, typename NTYPE>
  bool equals_reverse_bits(const ITYPE ival, const NTYPE num_bits)
  {
    ITYPE reverse_val = reverse_bits(ival, num_bits);
    return((ival == reverse_val));
  }


  /// @brief Return index of first one bit.
  /// - Return num_bits if all bits are zero.
  template <typename ITYPE, typename NTYPE>
  NTYPE get_first_one_bit(const ITYPE val, const NTYPE num_bits)
  {
    ITYPE mask = ITYPE(1);
    for (NTYPE i = 0; i < num_bits; i++) {
      if ((val & mask) != 0) { return(i); }
      mask = (mask << ITYPE(1));
    }
    return(num_bits);
  }


  /// @brief Return index of last one bit.
  /// - Return num_bits if all bits are zero.
  template <typename ITYPE, typename NTYPE>
  NTYPE get_last_one_bit(const ITYPE val, const NTYPE num_bits)
  {
    ITYPE mask = (ITYPE(1) << (num_bits-1));
    for (NTYPE i = 0; i < num_bits; i++) {
      if ((val & mask) != 0) { return(num_bits-1-i); }
      mask = (mask >> ITYPE(1));
    }
    return(num_bits);
  }


  /// Return true if ival has ones in positions k and (num_bits-1-k) and zeros everywhere else.
  template <typename ITYPE, typename NTYPE>
  bool is_two_opposite_ones(const ITYPE ival, const NTYPE num_bits)
  {
    NTYPE num_zeros, num_ones;
    IJK::count_bits(ival, num_bits, num_zeros, num_ones);
 
    if (equals_reverse_bits(ival, num_bits) && (num_ones == 2))
      { return(true); }
    else
      { return(false); }
  }


  /// Return true if ival has zeros in positions k and (num_bits-1-k) and ones everywhere else.
  template <typename ITYPE, typename NTYPE>
  bool is_two_opposite_zeros(const ITYPE ival, const NTYPE num_bits)
  {
    NTYPE num_zeros, num_ones;
    IJK::count_bits(ival, num_bits, num_zeros, num_ones);

    if (equals_reverse_bits(ival, num_bits) && (num_zeros == 2))
      { return(true); }
    else
      { return(false); }

  }


  /// @brief Return true if ival has ones in positions k and (num_bits-1-k).
  /// - Other bits can be zero or one.
  template <typename ITYPE, typename NTYPE>
  bool contains_two_opposite_ones(const ITYPE ival, const NTYPE num_bits)
  {
    const ITYPE ival_reverse_bits = IJK::reverse_bits(ival, num_bits);

    if ((ival & ival_reverse_bits) != 0) { return(true); }

    return(false);
  }


  /// @brief Return true if ival has zeros in positions k and (num_bits-1-k).
  /// - Other bits can be zero or one.
  template <typename ITYPE, typename NTYPE>
  bool contains_two_opposite_zeros(const ITYPE ival, const NTYPE num_bits)
  {
    const ITYPE all_ones = (~(ITYPE(0)));
    const ITYPE ival_complement = 
      (~ival) - ((all_ones >> num_bits) << num_bits);

    return(contains_two_opposite_ones(ival_complement, num_bits));
  }
    
}

#endif
