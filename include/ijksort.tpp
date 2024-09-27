/*!
 *  @file ijksort.tpp
 *  @brief Radix/counting sort integer values.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2023-2024 Rephael Wenger

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


#ifndef _IJKSORT_
#define _IJKSORT_

#include "ijk.tpp"

#include <algorithm>
#include <limits>


namespace IJKSORT {

  // *****************************************************************
  // FORWARD DECLARATIONS
  // *****************************************************************

  template <typename NTYPEB, typename NTYPELB>
  bool check_max_num_bins
    (const NTYPEB max_num_bins, const NTYPELB min_num,
     IJK::ERROR & error);

  template <typename LTYPE, typename NTYPE_IB,
            typename NTYPEB, typename NTYPEL>
  bool check_last_bin_size
  (const LTYPE first_element_in_bin[],
   const NTYPE_IB num_elements_in_bin[],
   const NTYPEB num_bins, const NTYPEL list_length,
   IJK::ERROR & error);

  template <typename LTYPEF, typename LTYPEC, typename NTYPE_IB,
            typename NTYPEB, typename NTYPEL>
  bool check_all_bins_filled
  (const LTYPEF first_element_in_bin[],
   const NTYPE_IB num_elements_in_bin[],
   const NTYPEB num_bins, const LTYPEC current[],
   const NTYPEL list_length, IJK::ERROR & error);

  template <typename LTYPEF, typename LTYPEC, typename NTYPE_IB,
            typename NTYPEB, typename NTYPEL>
  bool check_all_bins_filled
  (const LTYPEF first_element_in_bin[],
   const NTYPE_IB num_elements_in_bin[],
   const NTYPEB num_bins,
   const std::vector<LTYPEC> & current,
   const NTYPEL list_length, IJK::ERROR & error);

  template <typename RTYPE>
  bool check_range
  (const RTYPE rmin, const RTYPE rmax, IJK::ERROR & error);


  // *****************************************************************
  // RADIX SORT INFORMATION
  // *****************************************************************

  /// Information about execution of radix_sort.
  class RADIX_SORT_INFO {
  public:

    /// Number of algorithm iterations.
    unsigned int num_iter;

    /// Number of bins used by each level of radix_sort.
    std::vector<unsigned int> num_bins;

    /// log2_interval used by each level of radix_sort;
    std::vector<unsigned int> log2_interval;

    /// Add iteration.
    template <typename NTYPE, typename LTYPE>
    void AddIter(const NTYPE num_binsX, const LTYPE log2_intervalX)
    {
      num_iter++;
      num_bins.push_back(num_binsX);
      log2_interval.push_back(log2_intervalX);
    };

    /// Clear all information.
    void Clear()
    {
      num_iter = 0;
      num_bins.clear();
      log2_interval.clear();
    };
  };


  /// Empty class to replace RADIX_SORT_INFO when no information is returned.
  class RADIX_SORT_EMPTY_INFO {
  public:

    /// Dummy member function. Does nothing.
    template <typename NTYPE, typename LTYPE>
    void AddIter(const NTYPE num_binX, const LTYPE log2_intervalX)
    {
      // Do nothing.
    };

    /// Dummy member function. Does nothing.
    void Clear()
    {
      // Do nothing.
    };
  };

  // *****************************************************************
  //! @name COMPUTE FUNCTIONS
  // *****************************************************************

  //@{

  /*!
   *  @brief Compute the interval used for each bin.
   *  - Interval must be a power of two.
   *  @param rmin Lower bound on elements in list[].
   *  @param rmax Upper bound on elements in list[].
   *    @pre rmin <= rmax.
   *  @param max_val Maximum value covered by bins.
   *    @pre max_val >= 0.
   *  @param max_num_bins Maximum number of bins to cover [0,max_val].
   *    @pre max_num_bins >= 4.
   *  @param[out] log2_interval log base 2 of the interval.
   *    - log2_interval = recommended_log2_interval
   *      if (2^recommended_log2_interval)*max_num_bins > max_val.
   *  @param[out] interval Interval used for each bin.
   *    - interval = 2^{log2_interval}.
   */
  template <typename RTYPE, typename NTYPE,
            typename LTYPE, typename ITYPE>
  void compute_bin_interval
    (const RTYPE rmin, const RTYPE rmax, const NTYPE max_num_bins,
     LTYPE & log2_interval, ITYPE & interval)
  {
    const int FOUR(4);
    IJK::PROCEDURE_ERROR error("compute_bin_interval");

    if (!check_max_num_bins(max_num_bins, FOUR, error))
      { throw error; }

    const RTYPE range = (rmax-rmin)+1;
    interval = 1;
    log2_interval = 0;
    while (interval*max_num_bins < range) {

      // Extra check to avoid infinite loop.
      if (log2_interval+2 >= std::numeric_limits<ITYPE>::digits) {
        error.AddMessage
          ("Programming error. Not enough bits to store interval.");
        error.AddMessage("  max_num_bins = ", max_num_bins, ".");
        error.AddMessage("  rmin = ", rmin, ".");
        error.AddMessage("  rmax = ", rmax, ".");
        error.AddMessage("  interval = ", interval, ".");
        error.AddMessage
          ("  number of interval bits: ", std::numeric_limits<ITYPE>::digits);
        error.AddMessage
          ("  Need at least ", log2_interval+2, " bits.");
        throw error;
      }

      log2_interval++;
      interval = (interval << 1);
    }
  }


  /*!
   *  @brief Compute the number of bins to cover the range [rmin,rmax].
   *  - Bin ibin contains values in range [ibin*interval,(ibin+1)*interval-1].
   *  @param max_val Maximum value in range.
   *    @pre max_val > 0.
   *  @param interval Interval (number of integer values) covered by each bin.
   *    @pre interval > 0.
   */
  template <typename RTYPE, typename ITYPE, typename NTYPE>
  void compute_num_bins_to_cover_range
    (const RTYPE rmin, const RTYPE rmax, const ITYPE interval, NTYPE & num_bins)
  {
    const RTYPE range = (rmax-rmin)+1;
    IJK::PROCEDURE_ERROR error("compute_num_bins_to_cover_range");

    if (interval < 1) {
      error.AddMessage
        ("Programming error. interval must be at least 1.");
      error.AddMessage("  interval = ", interval, ".");
      throw error;
    }

    num_bins = range/interval;

    // Since max_val >= 0 and interval > 0, the following code
    //   ensures that num_bins>0.
    if (num_bins*interval < range)
      { num_bins++; }
  }


  /*!
   *  @brief Compute index of bin containing ival.
   *  - Bin containing ival is ((ival-rmin)/interval)
   *    where interval = 2^{log2_interval}.
   *  @param ival Integer value.
   *    @pre ival >= rmin.
   */
  template <typename ITYPE, typename RTYPE, typename LTYPE>
  ITYPE compute_bin_index
    (const ITYPE ival, const RTYPE rmin, const LTYPE log2_interval)
  {
    const ITYPE ibin = ((ival-rmin) >> log2_interval);
    return ibin;
  }

  //@}


  // *****************************************************************
  //! @name COUNT NUMBER OF ELEMENTS OR SET FIRST ELEMENT IN BINS
  // *****************************************************************

  //@{

  /// @brief Count number of elements in each bin.
  /// - Elements in each bin are identical.
  template <typename ITYPE, typename RTYPE,
            typename NTYPEL, typename NTYPEB, typename NTYPE_IB>
  void count_num_elements_in_each_bin
    (const ITYPE list[], const NTYPEL list_length,
     const RTYPE rmin, const NTYPEB num_bin, NTYPE_IB num_in_bin[])
  {
    // Initialize num_in_bin[ibin] to zero.
    for (NTYPEB ibin = 0; ibin < num_bin; ibin++)
      { num_in_bin[ibin] = 0; }

    // Count number of elements in each bin.
    for (NTYPEL i = 0; i < list_length; i++) {
      const ITYPE ival = list[i];
      const NTYPEB ibin = ival-rmin;
      num_in_bin[ibin]++;
    }
  }


  /*!
   *  @brief Count number of elements in each bin.
   *  - Bin ibin contains elements in the range
   *    [rmin+interval*ibin,rmin+(interval*ibin)-1]
   *    where interval = 2^log2_interval.
   */
  template <typename ITYPE, typename RTYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename NTYPE_IB>
  void count_num_elements_in_each_radix_sort_bin
    (const ITYPE list[], const NTYPEL list_length,
     const RTYPE rmin, const LTYPE log2_interval,
     const NTYPEB num_bin, NTYPE_IB num_in_bin[])
  {
    // Initialize num_in_bin[ibin] to zero.
    for (NTYPEB ibin = 0; ibin < num_bin; ibin++)
      { num_in_bin[ibin] = 0; }

    // Count number of elements in each bin.
    for (NTYPEL i = 0; i < list_length; i++) {
      const ITYPE ival = list[i];
      const NTYPEB ibin = compute_bin_index(ival, rmin, log2_interval);
      num_in_bin[ibin]++;
    }
  }


  /*!
   *  @brief Set first_element_in_bin[] for each bin.
   *  @param num_bin Number of bins.
   *  @param num_in_bin[ibin] Number of elements in bin ibin.
   *  @param[out] first_element_in_bin[ibin] Location of first
   *    element in bin ibin.
   */
  template <typename NTYPEB, typename NTYPE_IB, typename LTYPE>
  void set_first_element_in_each_bin
    (const NTYPEB num_bins, NTYPE_IB num_in_bin[],
     LTYPE first_element_in_bin[])
  {
    if (num_bins <= 0) {
      // No bins. Nothing to set.
      return;
    }

    first_element_in_bin[0] = 0;
    for (NTYPEB ibin = 1; ibin < num_bins; ibin++) {
      first_element_in_bin[ibin] =
        first_element_in_bin[ibin-1] + num_in_bin[ibin-1];
    }
  }

  //@}


  // *****************************************************************
  //! @name COUNTING SORT
  // *****************************************************************

  //@{

  /*!
   *  @brief Counting sort - Sort integers using bins.
   *  @pre max_num_bins >= (rmax-rmin+1).
   *  @pre rmin <= list[i] <= rmax.
   *    - Warning: counting_sort() does not check these preconditions.
   *    - Procedure may (will) throw a segmentation error
   *      if these preconditions are violated.
   *  @param list[] List of integers.
   *  @param list_length Number of elements in list[].
   *  @param rmin Lower bound on elements in list[].
   *    @pre rmin <= list[i] for every element in list.
   *  @param rmax Upper bound on elements in list[].
   *    @pre rmax >= list[i] for every element in list.
   *    @pre rmin <= rmax.
   *  @param max_num_bins Maximum number of bins.
   *    @pre max_num_bins >= (rmax-rmin+1).
   *  @param[out] new_list[] New list of items containing
   *      a permutation of list[].
   *    - new_list[j], new_list[j+1], ..., new_list[j+k],
   *      represents elements in bin ibin where
   *      j = first_element_in_bin[ibin] and
   *      k = num_elements_in_bin[ibin].
   *    - All elements in bin ibin have values
   *      less than the value of elements in bin (ibin+1).
   *    @pre new_list[] is preallocated to length list_length.
   *  @param[out] new_list_loc[i] Location of list[i] in new_list[].
   *    @pre new_list_loc[i] is preallocated to length list_length.
   *  @param[out] first_element_bin[ibin] Location in new_list[]
   *    of first element in bin ibin.
   *    @pre first_element_ibin[] is preallocated to length max_num_bin.
   *  @param[out] num_elements_in_bin[ibin] Number of elements
   *    in bin ibin.
   *    @pre num_elements_in_bin[] is preallocated
   *    to length max_num_bin.
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename NTYPE_IB>
  void counting_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   ITYPE new_list[],
   LTYPE new_list_loc[],
   LTYPE first_element_in_bin[],
   NTYPE_IB num_elements_in_bin[])
  {
    IJK::PROCEDURE_ERROR error("counting_sort");

    if (list_length == 0) {
      // Nothing to do.
      return;
    }

    if (!check_array_allocated(new_list, "new_list", error))
      { throw error; }

    if (!check_array_allocated(new_list_loc, "new_list_loc", error))
      { throw error; }

    const ITYPE range = (rmax-rmin+1);
    if (max_num_bins < range) {
      error.AddMessage
        ("Programming error. Too few bins for counting sort.");
      error.AddMessage
        ("  Need ", range, " bins to sort range [",
         rmin, ",", rmax, "].");
      error.AddMessage
        ("  Increase the number of bins or use radix_sort.");
      throw error;
    }

    const ITYPE num_bins = range;
    count_num_elements_in_each_bin
      (list, list_length, rmin, num_bins, num_elements_in_bin);

    set_first_element_in_each_bin
      (num_bins, num_elements_in_bin, first_element_in_bin);

    if (!check_last_bin_size
        (first_element_in_bin, num_elements_in_bin,
         num_bins, list_length, error))
      { throw error; }

    std::vector<LTYPE> current(num_bins);
    std::copy(first_element_in_bin, first_element_in_bin+num_bins,
              current.begin());

    for (NTYPEL i = 0; i < list_length; i++) {
      const NTYPEB ibin = list[i] - rmin;
      const LTYPE jloc = current[ibin];
      new_list[jloc] = list[i];
      new_list_loc[i] = jloc;
      current[ibin]++;
    }

    if (!check_all_bins_filled
        (first_element_in_bin, num_elements_in_bin, num_bins,
         current, list_length, error))
      { throw error; }
  }


  /*!
   *  @overload
   *  @brief Counting sort - Sort integers using bins. (Allocate bins.)
   *  - Version that determines num_bins and
   *    allocates first_element_in_bin[] and num_elements_in_bin[].
   *  - @param[out] num_bins Number of bins used in counting_sort().
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB>
  void counting_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   ITYPE new_list[],
   LTYPE new_list_loc[],
   NTYPEB & num_bins)
  {
    IJK::PROCEDURE_ERROR error("counting_sort");

    if (!check_range(rmin, rmax, error))
      { throw error; }

    const ITYPE range = 1 + rmax - rmin;
    num_bins = range;
    std::vector<LTYPE> first_element_in_bin(num_bins);
    std::vector<LTYPE> num_elements_in_bin(num_bins);

    counting_sort(list, list_length, rmin, rmax, num_bins,
                  new_list, new_list_loc,
                  IJK::vector2pointer(first_element_in_bin),
                  IJK::vector2pointer(num_elements_in_bin));
  }


  /*!
   *  @overload
   *  @brief Counting sort - Sort integers using bins. (C++ vector for new list.)
   *  - Version using STL C++ vector for arrays
   *    new_list[] and new_list_loc[].
   *  - This version resizes new_list[] and new_list_loc[]
   *    to size list_length.
   *    (new_list[] and new_list_loc[] DO NOT need to be
   *    preallocated to size list_length.)
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB>
  void counting_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc,
   NTYPEB & num_bins)
  {
    new_list.resize(list_length);
    new_list_loc.resize(list_length);

    counting_sort(list, list_length, rmin, rmax,
                  IJK::vector2pointer(new_list),
                  IJK::vector2pointer(new_list_loc),
                  num_bins);
  }


  /*!
   *  @overload
   *  @brief Counting sort - Sort integers using bins. (C++ vector.)
   *  - Version using STL C++ vector list[], new_list[]
   *    and new_list_loc[].
   *  - This version resizes new_list[] and new_list_loc[]
   *    to size list_length.
   *    (new_list[] and new_list_loc[] DO NOT need to be
   *    preallocated to size list_length.)
   */
  template <typename ITYPE, typename LTYPE, typename NTYPEB>
  void counting_sort
  (const std::vector<ITYPE> & list,
   const ITYPE rmin, const ITYPE rmax,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc,
   NTYPEB & num_bins)
  {
    counting_sort(IJK::vector2pointer(list), list.size(),
                  rmin, rmax, new_list, new_list_loc, num_bins);
  }


  /*!
   *  @overload
   *  @brief Counting sort - Sort integers using bins. (C++ vector.)
   *  - Version that does not return \a num_bins.
   *  - Version using STL C++ vector list[], new_list[]
   *    and new_list_loc[].
   *  - This version resizes new_list[] and new_list_loc[]
   *    to size list_length.
   *    (new_list[] and new_list_loc[] DO NOT need to be
   *    preallocated to size list_length.)
   */
  template <typename ITYPE, typename LTYPE>
  void counting_sort
  (const std::vector<ITYPE> & list,
   const ITYPE rmin, const ITYPE rmax,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc)
  {
    unsigned int num_bins;
    counting_sort(list, rmin, rmax, new_list, new_list_loc, num_bins);
  }


  /*!
   *  @brief Counting sort - Sort integers using bins. 
   *    (Compute range and allocate bins.)
   *  - Compute range from list.
   *  - Memory is proportional to (max_val-min_val) where
   *    max_val and min_val are max and min elements of list[], respectively.
   *  - Running time is proportional to max((max_val-min_val),list.size()).
   *  - WARNING: It is very easy to run out of memory or use very large
   *    amounts of memory if the difference between max_val and min_val
   *    is arbitrarily large.
   *  @param list[] List of integers.
   *  @param[out] new_list[] New list of items containing
   *      a permutation of list[].
   *    - new_list[j], new_list[j+1], ..., new_list[j+k],
   *      represents elements in bin ibin where
   *      j = first_element_in_bin[ibin] and
   *      k = num_elements_in_bin[ibin].
   *    - All elements in bin ibin have values
   *      less than the value of elements in bin (ibin+1).
   *  @param[out] new_list_loc[i] Location of list[i] in new_list[].
   *  @param[out] num_bins Number of bins used in counting_sort_compute_range().
   */
  template <typename ITYPE, typename LTYPE, typename NTYPEB>
  void counting_sort_compute_range
  (const std::vector<ITYPE> & list,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc,
   NTYPEB & num_bins)
  {
    const ITYPE min_val =
      *min_element(list.begin(), list.end());
    const ITYPE max_val =
      *max_element(list.begin(), list.end());

    counting_sort
      (list, min_val, max_val, new_list, new_list_loc, num_bins);
  }


  /*!
   *  @overload
   *  @brief Counting sort - Sort integers using bins. 
   *    (Compute range and allocate bins.)
   *  - Version that does not return num_bins.
   */
  template <typename ITYPE, typename LTYPE>
  void counting_sort_compute_range
  (const std::vector<ITYPE> & list,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc)
  {
    ITYPE num_bins;
    counting_sort_compute_range
      (list, new_list, new_list_loc, num_bins);
  }

  //@}


  // *****************************************************************
  //! @name RADIX SORT
  // *****************************************************************

  //@{

  /*!
   *  @brief Reorder (partially sort) list of integers using bins.
   *  - Algorithm trades memory for running time.
   *  - Memory used is proportional to {max_num_bins}
   *  - Running time is proportional to
   *    {list_length}*{log_2(rmax-rmin)/log2(max_num_bins)}.
   *  @param list[] List of integers.
   *  @param list_length Number of elements in list[].
   *  @param rmin Lower bound on elements in list[].
   *    @pre rmin <= list[i] for every element in list.
   *  @param rmax Upper bound on elements in list[].
   *    @pre rmax >= list[i] for every element in list.
   *    @pre rmin <= rmax.
   *  @param max_num_bins Maximum number of bins.
   *    - Algorithm uses memory proportional to max_num_bins.
   *    - Algorithm runs faster with larger values of max_num_bins.
   *    - Recommendation: max_num_bins should be at least 10,
   *      preferably at least 1000.
   *    @pre max_num_bins >= 4.
   *  @param log2_interval Log base 2 of the interval covered by each bin.
   *    @pre (2^{log2_interval}*max_num_bins >= (rmax-rmin+1).
   *    - Use compute_bin_interval() to compute the correct interval
   *    to use for a given range and number of bins.
   *  @param[out] new_list[] New list of items containing 
   *      a permutation of list[].
   *    - new_list[j], new_list[j+1], ..., new_list[j+k], 
   *      represents elements in bin ibin where 
   *      j = first_element_in_bin[ibin] and
   *      k = num_elements_in_bin[ibin].
   *    - All elements in bin ibin have values
   *      less than the value of elements in bin (ibin+1).
   *    @pre new_list[] is preallocated to length list_length.
   *  @param[out] new_list_loc[i] Location of list[i] in new_list[].
   *    @pre new_list_loc[i] is preallocated to length list_length.
   *  @param[out] first_element_bin[ibin] Location in new_list[]
   *    of first element in bin ibin.
   *    @pre first_element_ibin[] is preallocated to length max_num_bin.
   *  @param[out] num_elements_in_bin[ibin] Number of elements
   *    in bin ibin.
   *    @pre num_elements_in_bin[] is preallocated 
   *    to length max_num_bin.
   *  @param temp_bin[] Temporary array for use by algorithm.
   *    @pre temp_bin[] is preallocated to length max_num_bin.
   *  @param[out] num_bins Actual number of bins used
   *    by split_integer_into_bins().
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename NTYPEB2,
            typename NTYPE_IB>
  void split_integer_list_into_bins
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   const ITYPE log2_interval,
   ITYPE new_list[],
   LTYPE new_list_loc[],
   LTYPE first_element_in_bin[],
   NTYPE_IB num_elements_in_bin[],
   LTYPE temp_bin[],
   NTYPEB2 & num_bins)
  {
    const int FOUR(4);
    IJK::PROCEDURE_ERROR error("split_integer_list_into_bins");

    if (list_length < 1) {
      // Nothing to split.
      return;
    }

    if (!check_max_num_bins(max_num_bins, FOUR, error))
      { throw error; }

    if (!check_range(rmin, rmax, error))
      { throw error; }

    const ITYPE range = rmax-rmin+1;
    const ITYPE interval = (1 << log2_interval);
    if (interval*max_num_bins < range) {
      ITYPE required_log2_interval, required_interval;
      compute_bin_interval
        (rmin, rmax, max_num_bins,
         required_log2_interval, required_interval);
      error.AddMessage("Programming error. Value log2_interval too small.");
      error.AddMessage("  log2_interval = ", log2_interval, ".");
      error.AddMessage
        ("  log2_interval should be ", required_log2_interval,
         " for range [", rmin, ",", rmax, "] and max_num_bins ",
         max_num_bins, ".");
      throw error;
    }

    compute_num_bins_to_cover_range(rmin, rmax, interval, num_bins);

    if (num_bins < 1) {
      // Extra check.
      error.AddMessage("Programming error. num_bins < 1.");
      error.AddMessage("  num_bins = ", num_bins, ".");
      throw error;
    }

    if (num_bins > max_num_bins) {
      error.AddMessage("Programming error. num_bins > max_num_bins");
      error.AddMessage("  num_bins = ", num_bins, ".");
      error.AddMessage("  max_num_bins = ", max_num_bins, ".");
      throw error;
    }

    if (num_bins*interval < range) {
      error.AddMessage("Programming error. num_bins*interval < range.");
      error.AddMessage("  num_bins = ", num_bins, ".");
      error.AddMessage("  interval = ", interval, ".");
      error.AddMessage("  range = ", range, ".");
      throw error;
    }
    
    count_num_elements_in_each_radix_sort_bin
      (list, list_length, rmin, log2_interval, num_bins,
       num_elements_in_bin);

    set_first_element_in_each_bin
      (num_bins, num_elements_in_bin, first_element_in_bin);

    if (!check_last_bin_size
        (first_element_in_bin, num_elements_in_bin,
         num_bins, list_length, error))
      { throw error; }

    LTYPE * current = temp_bin;
    std::copy(first_element_in_bin, first_element_in_bin+num_bins,
              current);

    for (NTYPEL i = 0; i < list_length; i++) {
      const NTYPEB ibin = compute_bin_index(list[i], rmin, log2_interval);
      const LTYPE jloc = current[ibin];
      new_list[jloc] = list[i];
      new_list_loc[i] = jloc;
      current[ibin]++;
    }

    if (!check_all_bins_filled
        (first_element_in_bin, num_elements_in_bin, num_bins,
         current, list_length, error))
      { throw error; }
  };


  /*!
   *  @brief Radix sort order list of integers using bins.
   *  - Algorithm trades memory for running time.
   *  - Memory used is proportional to {max_num_bins}
   *  - Running time is proportional to
   *    {list_length}*log_2(rmax-rmin)/{log2_interval}.
   *  @param list[] List of integers.
   *  @param list_length Number of elements in list[].
   *  @param rmin Lower bound on elements in list[].
   *    @pre rmin <= list[i] for every element in list.
   *  @param rmax Upper bound on elements in list[].
   *    @pre rmax >= list[i] for every element in list.
   *    @pre rmin <= rmax.
   *  @param max_num_bins Maximum number of bins.
   *    - Algorithm uses memory proportional to max_num_bins.
   *    - Algorithm runs faster with larger values of max_num_bins.
   *    - Recommendation: max_num_bins should be at least 10,
   *      preferably at least 1000.
   *    @pre max_num_bins >= 4.
   *  @param[out] new_list[] New list of items containing
   *      a permutation of list[].
   *    - new_list[j], new_list[j+1], ..., new_list[j+k],
   *      represents elements in bin ibin where
   *      j = first_element_in_bin[ibin] and
   *      k = num_elements_in_bin[ibin].
   *    - All elements in bin ibin have values
   *      less than the value of elements in bin (ibin+1).
   *    @pre new_list[] is preallocated to length list_length.
   *  @param[out] new_list_loc[i] Location of list[i] in new_list[].
   *    @pre new_list_loc[i] is preallocated to length list_length.
   *  @param[out] first_element_bin[ibin] Location in new_list[]
   *    of first element in bin ibin.
   *    @pre first_element_ibin[] is preallocated to length max_num_bins.
   *  @param[out] num_elements_in_bin[ibin] Number of elements
   *    in bin ibin.
   *    @pre num_elements_in_bin[] is preallocated
   *    to length max_num_bins.
   *  @param temp_bin[] Temporary array for use by algorithm.
   *    @pre temp_bin[] is preallocated to length max_num_bin.
   *  @param[out] radix_sort_info Information about execution of radix_sort.
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename NTYPE_IB,
            typename INFO_TYPE>
  void radix_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   ITYPE new_list[],
   LTYPE new_list_loc[],
   LTYPE first_element_in_bin[],
   NTYPE_IB num_elements_in_bin[],
   LTYPE temp_bin[],
   INFO_TYPE & radix_sort_info)
  {
    IJK::PROCEDURE_ERROR error("radix_sort");

    // Initialize.
    radix_sort_info.Clear();

    if (list_length == 0) {
      // Nothing to do.
      return;
    }

    if (!check_array_allocated(new_list, "new_list", error))
      { throw error; }

    if (!check_array_allocated(new_list_loc, "new_list_loc", error))
      { throw error; }

    const ITYPE range = 1 + rmax-rmin;
    if (max_num_bins >= range) {
      // Use counting_sort.
      const NTYPEB num_bins = range;
      counting_sort
        (list, list_length, rmin, rmax, num_bins,
         new_list, new_list_loc,
         first_element_in_bin, num_elements_in_bin);
      radix_sort_info.AddIter(num_bins, 0);
      return;
    }

    std::vector<ITYPE> temp_listA(list_length);
    std::vector<LTYPE> temp_locA(list_length);
    ITYPE log2_intervalA, intervalA;
    compute_bin_interval
      (rmin, rmax, max_num_bins, log2_intervalA, intervalA);

    // Actual number of bins used by split_integer_list_into_bins().
    NTYPEB num_bins;
    split_integer_list_into_bins
    (list, list_length, rmin, rmax, max_num_bins, log2_intervalA,
     IJK::vector2pointer(temp_listA),
     IJK::vector2pointer(temp_locA),
     first_element_in_bin, num_elements_in_bin,
     temp_bin, num_bins);
    radix_sort_info.AddIter(num_bins, log2_intervalA);

    std::vector<ITYPE> temp_listC(list_length);
    std::vector<LTYPE> temp_locC(list_length);

    while (intervalA > max_num_bins) {
      ITYPE rminB = rmin;
      ITYPE rmaxB = rminB+intervalA-1;
      NTYPEL iB = 0;
      ITYPE log2_intervalB, intervalB;

      compute_bin_interval
        (rminB, rmaxB, max_num_bins, log2_intervalB, intervalB);

      while (iB < list_length) {
        NTYPEL iendB = iB;
        while ((iendB < list_length) &&
               (temp_listA[iendB] <= rmaxB))
        { iendB++; }

        if (iB < iendB) {
          ITYPE * listB = &(temp_listA[iB]);
          ITYPE * listD = &(temp_listC[iB]);
          LTYPE * list_locD = &(temp_locC[iB]);
          NTYPEL list_lengthB = iendB-iB;
          split_integer_list_into_bins
          (listB, list_lengthB, rminB, rmaxB, max_num_bins,
           log2_intervalB, listD, list_locD,
           first_element_in_bin, num_elements_in_bin, temp_bin,
           num_bins);
          if (iB == 0) {
             // Add radix sort information for level only once.
            radix_sort_info.AddIter(num_bins, log2_intervalB);
          }

          // Adjust temp_locD[] to correctly point to location
          //   in temp_listC[].
          for (NTYPEL iD = 0; iD < list_lengthB; iD++) {
            list_locD[iD] = list_locD[iD] + iB;
          }
        }

        rminB = rmaxB+1;
        rmaxB = rminB+intervalA-1;
        iB = iendB;
      }

      // Copy temp_listC[] into temp_listA[].
      std::copy(temp_listC.begin(), temp_listC.end(),
                temp_listA.begin());

      // Adjust list_locA[] to correctly point to location in temp_listA[].
      for (NTYPEL iA = 0; iA < list_length; iA++) {
        const LTYPE old_locA = temp_locA[iA];
        temp_locA[iA] = temp_locC[old_locA];
      }

      // Update the interval.
      log2_intervalA = log2_intervalB;
      intervalA = intervalB;
    }

    // Apply counting_sort() to complete sort.
    ITYPE rminE = rmin;
    ITYPE rmaxE = rminE+intervalA-1;
    NTYPEL iE = 0;
    num_bins = intervalA;
    radix_sort_info.AddIter(num_bins, 0);

    while (iE < list_length) {
      NTYPEL iendE = iE;
      while ((iendE < list_length) &&
             (temp_listA[iendE] <= rmaxE))
      { iendE++; }

      if (iE < iendE) {
        ITYPE * listE = &(temp_listA[iE]);
        ITYPE * listF = &(new_list[iE]);
        LTYPE * list_locF = &(temp_locC[iE]);
        NTYPEL list_lengthE = iendE-iE;

        counting_sort
          (listE, list_lengthE, rminE, rmaxE, num_bins,
           listF, list_locF,
           first_element_in_bin, num_elements_in_bin);

        // Adjust list_locF[] to correctly point to location
        //   in new_list[].
        for (NTYPEL jF = 0; jF < list_lengthE; jF++) {
          list_locF[jF] = list_locF[jF] + iE;
        }
      }

      rminE = rmaxE+1;
      rmaxE = rminE+intervalA-1;
      iE = iendE;
    }

    // Adjust new_list_loc[] to correctly point to location in new_list[].
    for (NTYPEL i = 0; i < list_length; i++) {
      new_list_loc[i] = temp_locC[temp_locA[i]];
    }
  }


  /*!
   *  @overload
   *  @brief Radix sort order list of integers using bins. (Allocate temp_bin[].)
   *  - Version that allocates array temp_bin[].
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename NTYPE_IB,
            typename INFO_TYPE>
  void radix_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   ITYPE new_list[],
   LTYPE new_list_loc[],
   LTYPE first_element_in_bin[],
   NTYPE_IB num_elements_in_bin[],
   INFO_TYPE & radix_sort_info)
  {
    std::vector<LTYPE> temp_bin(max_num_bins);

    radix_sort(list, list_length, rmin, rmax, max_num_bins,
               new_list, new_list_loc,
               first_element_in_bin, num_elements_in_bin,
               IJK::vector2pointer(temp_bin),
               radix_sort_info);
  }

  /*!
   *  @overload
   *  @brief Radix sort order list of integers using bins. (C++ vector for bins.)
   *  - Version using C++ STL vector for
   *    arrays first_element_in_bin[] and num_elements_in_bin[].
   *  - This version resizes first_element_in_bin[] and
   *    num_elements_in_bin[] to size max_num_bins.
   *    (first_element_in_bin[] and num_elemens_in_bin[] DO NOT need to be
   *    preallocated to size max_num_bins.)
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename NTYPE_IB,
            typename INFO_TYPE>
  void radix_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   ITYPE new_list[],
   LTYPE new_list_loc[],
   std::vector<LTYPE> & first_element_in_bin,
   std::vector<NTYPE_IB> & num_elements_in_bin,
   INFO_TYPE & radix_sort_info)
  {
    first_element_in_bin.resize(max_num_bins);
    num_elements_in_bin.resize(max_num_bins);

    radix_sort
    (list, list_length, rmin, rmax, max_num_bins,
     new_list, new_list_loc,
     IJK::vector2pointer(first_element_in_bin),
     IJK::vector2pointer(num_elements_in_bin),
     radix_sort_info);
  }


  /*!
   *  @overload
   *  @brief Radix sort order list of integers using bins. (Allocate bins.)
   *  - Version that creates arrays
   *    first_element_in_bin[] and num_elements_in_bin[].
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename INFO_TYPE>
  void radix_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   ITYPE new_list[], LTYPE new_list_loc[],
   INFO_TYPE & radix_sort_info)
  {
    std::vector<LTYPE> first_element_in_bin(max_num_bins);
    std::vector<NTYPEL> num_elements_in_bin(max_num_bins);

    radix_sort(list, list_length, rmin, rmax, max_num_bins,
               new_list, new_list_loc,
               IJK::vector2pointer(first_element_in_bin),
               IJK::vector2pointer(num_elements_in_bin),
               radix_sort_info);
  }


  /*!
   *  @overload
   *  @brief Radix sort order list of integers using bins. (C++ vector for new list.)
   *  - Version using C++ STL vector for arrays new_list[]
   *    and new_list_loc[].
   *  - This version resizes new_list[] and new_list_loc[]
   *    to size list_length.
   *    (new_list[] and new_list_loc[] DO NOT need to be
   *    preallocated to size list_length.)
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEL, typename NTYPEB, typename INFO_TYPE>
  void radix_sort
  (const ITYPE list[], const NTYPEL list_length,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc,
   INFO_TYPE & radix_sort_info)
  {
    new_list.resize(list_length);
    new_list_loc.resize(list_length);

    radix_sort(list, list_length, rmin, rmax, max_num_bins,
               IJK::vector2pointer(new_list),
               IJK::vector2pointer(new_list_loc),
               radix_sort_info);
  }


  /*!
   *  @overload
   *  @brief Radix sort order list of integers using bins. (C++ vector.)
   *  - Version using C++ STL vector for all arrays.
   *  - This version resizes new_list[] and new_list_loc[]
   *    to size list_length.
   *    (new_list[] and new_list_loc[] DO NOT need to be
   *    preallocated to size list_length.)
   */
  template <typename ITYPE, typename LTYPE, typename NTYPEB,
            typename INFO_TYPE>
  void radix_sort
  (const std::vector<ITYPE> & list,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc,
   INFO_TYPE & radix_sort_info)
  {
    radix_sort(IJK::vector2pointer(list), list.size(),
               rmin, rmax, max_num_bins, new_list, new_list_loc,
               radix_sort_info);
  }

  /*!
   *  @overload
   *  @brief Radix sort order list of integers using bins. (C++ vector. No info.)
   *  - Version using C++ STL vector for all arrays.
   *  - Version that does not return radix_sort_info.
   *  - This version resizes new_list[] and new_list_loc[]
   *    to size list_length.
   *    (new_list[] and new_list_loc[] DO NOT need to be
   *    preallocated to size list_length.)
   */
  template <typename ITYPE, typename LTYPE,
            typename NTYPEB>
  void radix_sort
  (const std::vector<ITYPE> & list,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB max_num_bins,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc)
  {
    RADIX_SORT_EMPTY_INFO radix_sort_info;
    radix_sort(IJK::vector2pointer(list), list.size(),
               rmin, rmax, max_num_bins, new_list, new_list_loc,
               radix_sort_info);
  }


  /*!
   *  @brief Radix sort order list of integers using bins.
   *  - Allocates between sqrt(range) and 2*sqrt(range) bins
   *    where range = (rmax-rmin)+1.
   *  - Algorithm has at most 2 splitting/sorting iterations.
   *  @param list[] List of integers.
   *  @param list_length Number of elements in list[].
   *  @param rmin Lower bound on elements in list[].
   *    @pre rmin <= list[i] for every element in list.
   *  @param rmax Upper bound on elements in list[].
   *    @pre rmax >= list[i] for every element in list.
   *    @pre rmin <= rmax.
   *  @param min_num_bins MINIMUM (not maximum) number of bins.
   *  - Use at least \a min_num_bins.
   *  - Switches to counting_sort() if (min_num_bins > (rmax-rmin)+1).
   *  @param[out] new_list[] New list of items containing
   *      a permutation of list[].
   *    - new_list[j], new_list[j+1], ..., new_list[j+k],
   *      represents elements in bin ibin where
   *      j = first_element_in_bin[ibin] and
   *      k = num_elements_in_bin[ibin].
   *    - All elements in bin ibin have values
   *      less than the value of elements in bin (ibin+1).
   *    - Resized by radix_sort_iter2() to size list_length.
   *      (Does not need to be preallocated to size list_length.)
   *  @param[out] new_list_loc[i] Location of list[i] in new_list[].
   *    - Resized by radix_sort_iter2() to size list_length.
   *      (Does not need to be preallocated to size list_length.)
   *  @param[out] radix_sort_info Information about execution of radix_sort.
   */
   template <typename ITYPE, typename NTYPEL, typename NTYPEB,
             typename LTYPE, typename INFO_TYPE>
   void radix_sort_iter2
   (const ITYPE list[],
    const NTYPEL list_length,
    const ITYPE rmin, const ITYPE rmax,
    const NTYPEB min_num_bins,
    std::vector<ITYPE> & new_list,
    std::vector<LTYPE> & new_list_loc,
    INFO_TYPE & radix_sort_info)
   {
     const int FOUR(4);
     NTYPEB num_bins;

     radix_sort_info.Clear();

     if (list_length == 0) {
       // Nothing to do;
       return;
     }

     const ITYPE range = (rmax-rmin) + 1;

     if ((min_num_bins >= range) || (range <= 1)) {
       counting_sort
         (list, list_length, rmin, rmax, new_list, new_list_loc, num_bins);
       radix_sort_info.AddIter(num_bins, 0);
     }

     const ITYPE log2_interval = int(floor(std::log2(double(range))/2.0));
     const ITYPE interval = (1 << log2_interval);
     compute_num_bins_to_cover_range(rmin, rmax, interval, num_bins);

     if (num_bins < interval) { num_bins = interval; }

     if (num_bins*interval < FOUR) { num_bins = FOUR; }
     radix_sort
       (list, list_length, rmin, rmax, num_bins,
        new_list, new_list_loc, radix_sort_info);
   }


   /*!
    *  @brief Radix sort order list of integers using bins. (C++ vector.)
    *  - Version using C++ STL vector for array list[].
    */
   template <typename ITYPE, typename LTYPE, typename NTYPEB,
             typename INFO_TYPE>
   void radix_sort_iter2
   (const std::vector<ITYPE> & list,
    const ITYPE rmin, const ITYPE rmax,
    const NTYPEB min_num_bins,
    std::vector<ITYPE> & new_list,
    std::vector<LTYPE> & new_list_loc,
    INFO_TYPE & radix_sort_info)
   {
     radix_sort_iter2
       (IJK::vector2pointer(list), list.size(),
        rmin, rmax, min_num_bins,
        new_list, new_list_loc, radix_sort_info);
   }


  /*!
   *  @overload
   *  @brief Radix sort order list of integers using bins. (No info.)
   *  - Version that does not return radix_sort_info.
   */
  template <typename ITYPE, typename LTYPE, typename NTYPEB>
  void radix_sort_iter2
  (const std::vector<ITYPE> & list,
   const ITYPE rmin, const ITYPE rmax,
   const NTYPEB min_num_bins,
   std::vector<ITYPE> & new_list,
   std::vector<LTYPE> & new_list_loc)
  {
    RADIX_SORT_EMPTY_INFO radix_sort_info;
    radix_sort_iter2
      (list, rmin, rmax, min_num_bins, new_list, new_list_loc,
       radix_sort_info);
  }

  //@}



  // *****************************************************************
  //! @name MERGE IDENTICAL VALUES IN AN INTEGER LIST.
  // *****************************************************************

  /*!
   *  @brief Merge identical values in an integer list using radix_sort_iter2().
   *  @param list0 List of integers in range [rmin,rmax].
   *  @param list0_length Length of list0.
   *  @param[out] list1_nodup List without any duplicate values.
   *  @param[out] list0_map Mapping from list0 to locations in list1_nodup.
   *    @pre Array list0_map is preallocated to length at least list0_length.
   *  @param rmin Lower bound on elements in list[].
   *    @pre rmin <= list[i] for every element in list.
   *  @param rmax Upper bound on elements in list[].
   *    @pre rmax >= list[i] for every element in list.
   *    @pre rmin <= rmax.
   *  @param min_num_bins Minimum number of bins.
   *    - Switches to counting_sort() if (min_num_bins > (rmax-rmin)+1).
   */
  template <typename ITYPE, typename NTYPEL, typename NTYPEB,
            typename RTYPE, typename MTYPE, typename INFO_TYPE>
  void merge_identical_radix_sort_iter2
  (const ITYPE * list0, const NTYPEL list0_length,
   const RTYPE rmin, const RTYPE rmax,
   const NTYPEB min_num_bins,
   std::vector<ITYPE> & list1_nodup,
   MTYPE * list0_map,
   INFO_TYPE & radix_sort_info)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    std::vector<ITYPE> sorted_list;
    std::vector<MTYPE> sorted_list_loc;

    // Initialize
    list1_nodup.clear();

    if (list0_length == 0) {
      // Nothing to do.
      radix_sort_info.Clear();
      return;
    }

    radix_sort_iter2
      (list0, list0_length, rmin, rmax, min_num_bins,
       sorted_list, sorted_list_loc, radix_sort_info);

    // sorted_list_map[i] = location in list1_nodup[] of sorted_list[i].
    std::vector<MTYPE> sorted_list_map(list0_length);
    MTYPE iloc1 = 0;
    ITYPE prev_val = sorted_list[0];
    list1_nodup.push_back(prev_val);
    sorted_list_map[0] = 0;
    for (SIZE_TYPE i = 1; i < sorted_list.size(); i++) {
      if (sorted_list[i] != prev_val) {
        prev_val = sorted_list[i];
        list1_nodup.push_back(prev_val);
        iloc1++;
      }
      sorted_list_map[i] = iloc1;
    }

    // Set list0_map[].
    for (NTYPEL i = 0; i < list0_length; i++) {
      list0_map[i] = sorted_list_map[sorted_list_loc[i]];
    }
  }


  /*!
   *  @overload
   *  @brief Merge identical values in an integer list 
   *    using radix_sort_iter2(). (C++ vector list0_map.)
   *  - Version using C++ STL vector for array list0_map[].
   *  - This version resizes list0_map to size list0_length.
   *    (list0_map[] DOES NOT need to be preallocated to size list0_length.)
   */
  template <typename ITYPE, typename NTYPEL, typename NTYPEB,
            typename RTYPE, typename MTYPE, typename INFO_TYPE>
  void merge_identical_radix_sort_iter2
  (const ITYPE * list0, const NTYPEL list0_length,
   const RTYPE rmin, const RTYPE rmax,
   const NTYPEB min_num_bins,
   std::vector<ITYPE> & list1_nodup,
   std::vector<MTYPE> & list0_map,
   INFO_TYPE & radix_sort_info)
  {
    if (list0_length == 0) {
      // Nothing to do.
      return;
    }

    list0_map.resize(list0_length);

    merge_identical_radix_sort_iter2
    (list0, list0_length, rmin, rmax, min_num_bins, list1_nodup,
     IJK::vector2pointer(list0_map), radix_sort_info);
  }


  /*!
   *  @overload
   *  @brief Merge identical values in an integer list using radix_sort_iter2(). (C++ vector.)
   *  - Version using C++ STL vector for all arrays.
   *  - This version resizes list0_map to size list0.size().
   *    (list0_map[] DOES NOT need to be preallocated to size list0.size().)
   */
  template <typename ITYPE, typename NTYPEB,
            typename RTYPE, typename MTYPE, typename INFO_TYPE>
  void merge_identical_radix_sort_iter2
  (const std::vector<ITYPE> & list0,
   const RTYPE rmin, const RTYPE rmax,
   const NTYPEB min_num_bins,
   std::vector<ITYPE> & list1_nodup,
   std::vector<MTYPE> & list0_map,
   INFO_TYPE & radix_sort_info)
  {
    merge_identical_radix_sort_iter2
    (IJK::vector2pointer(list0), list0.size(),
     rmin, rmax, min_num_bins, list1_nodup,
     list0_map, radix_sort_info);
  }


  /*!
   *  @overload
   *  @brief Merge identical values in an integer list using radix_sort_iter2(). (No info.)
   *  - Version that does not return radix_sort_info.
   *  - Version using C++ STL vector for all arrays.
   *  - This version resizes list0_map to size list0.size().
   *    (list0_map[] DOES NOT need to be preallocated to size list0.size().)
   */
  template <typename ITYPE, typename NTYPEB,
            typename RTYPE, typename MTYPE>
  void merge_identical_radix_sort_iter2
  (const std::vector<ITYPE> & list0,
   const RTYPE rmin, const RTYPE rmax,
   const NTYPEB min_num_bins,
   std::vector<ITYPE> & list1_nodup,
   std::vector<MTYPE> & list0_map)
  {
    RADIX_SORT_EMPTY_INFO radix_sort_info;

    merge_identical_radix_sort_iter2
      (list0, rmin, rmax, min_num_bins, list1_nodup, list0_map,
       radix_sort_info);
  }


  // *****************************************************************
  //! @name CHECK FUNCTIONS
  // *****************************************************************

  //@{

  /// Return false and set error message if rmax < rmin.
  template <typename RTYPE>
  bool check_range
  (const RTYPE rmin, const RTYPE rmax, IJK::ERROR & error)
  {
    if (rmin > rmax) {
      error.AddMessage("Rrror. rmin > rmax.");
      error.AddMessage("  rmin = ", rmin, ".");
      error.AddMessage("  rmax = ", rmax, ".");
      return false;
    }

    return true;
  }


  /// Return false and set error message if max_num_bins < lower_bound.
  template <typename NTYPEB, typename NTYPELB>
  bool check_max_num_bins
    (const NTYPEB max_num_bins, const NTYPELB lower_bound,
     IJK::ERROR & error)
  {
    if (max_num_bins < lower_bound) {
      error.AddMessage
        ("Error. max_num_bins should be at least ",
         lower_bound, ".");
      error.AddMessage
        ("  max_num_bins = ", max_num_bins, ".");
      return false;
    }

    return true;
  }


  /// @brief Check last bin size against list_length.
  /// - Return false and set error if check fails.
  template <typename LTYPE, typename NTYPE_IB,
            typename NTYPEB, typename NTYPEL>
  bool check_last_bin_size
  (const LTYPE first_element_in_bin[],
   const NTYPE_IB num_elements_in_bin[],
   const NTYPEB num_bins,
   const NTYPEL list_length,
   IJK::ERROR & error)
  {
    if (num_bins < 1) {
      // Nothing to check.
      return true;
    }

    const NTYPEB ilast_bin = num_bins-1;
    const NTYPEL iend =
      first_element_in_bin[ilast_bin] + num_elements_in_bin[ilast_bin];
    if (iend != list_length) {
      error.AddMessage
        ("Error in creating bin arrays.");
      error.AddMessage("  num_bins = ", num_bins, ".");
      error.AddMessage
        ("  first_element_in_bin[", ilast_bin, "] = ",
         first_element_in_bin[ilast_bin], ".");
      error.AddMessage
        ("  num_elements_in_bin[", ilast_bin, "] = ",
         num_elements_in_bin[ilast_bin], ".");
      error.AddMessage
        ("  Does not sum to list_length = ", list_length, ".");
      return false;
    }

    return true;
  }


  /// @brief Check that all bins are filled.
  /// - Return false and set error message if check fails.
  template <typename LTYPEF, typename LTYPEC, typename NTYPE_IB,
            typename NTYPEB, typename NTYPEL>
  bool check_all_bins_filled
  (const LTYPEF first_element_in_bin[],
   const NTYPE_IB num_elements_in_bin[],
   const NTYPEB num_bins,
   const LTYPEC current[],
   const NTYPEL list_length,
   IJK::ERROR & error)
  {
    if (num_bins < 1) {
      // Nothing to check.
      return true;
    }

    for (NTYPEB ibin = 1; ibin < num_bins; ibin++) {
      if (current[ibin-1] != first_element_in_bin[ibin]) {
        error.AddMessage
          ("Error. Incorrect insertion of list[] into bins.");
        error.AddMessage
          ("  Inserted ", current[ibin]-first_element_in_bin[ibin-1],
           " elements into bin ", ibin, ".");
        error.AddMessage
          ("  Expected ", num_elements_in_bin[ibin], " elements in bin ",
           ibin, ".");
        return false;
      }
    }

    const NTYPEB ilast_bin = num_bins-1;
    if (NTYPEL(current[ilast_bin]) != list_length) {
      error.AddMessage
        ("Error. Incorrect insertion of list[] into bins.");
      error.AddMessage
        ("  Inserted ", current[ilast_bin]-first_element_in_bin[ilast_bin-1],
         " elements into last bin ", ilast_bin, ".");
      error.AddMessage
        ("  Expected ", num_elements_in_bin[ilast_bin], " elements in last bin ",
         ilast_bin, ".");
      return false;
    }

    return true;
  }


  /// @brief Check that all bins are filled.
  /// - Return false and set error message if check fails.
  /// - Version using C++ STL vector for array current[].
  template <typename LTYPEF, typename LTYPEC, typename NTYPE_IB,
            typename NTYPEB, typename NTYPEL>
  bool check_all_bins_filled
  (const LTYPEF first_element_in_bin[],
   const NTYPE_IB num_elements_in_bin[],
   const NTYPEB num_bins,
   const std::vector<LTYPEC> & current,
   const NTYPEL list_length,
   IJK::ERROR & error)
  {
    return check_all_bins_filled
        (first_element_in_bin, num_elements_in_bin,
         num_bins, IJK::vector2pointer(current),
         list_length, error);
  }


  /*!
   *  @brief Return true if elements are sorted in non-decreasing order.
   *  - Return false and set error message if list[i] > list[i+1]
   *    for some i.
   */
  template <typename ITYPE, typename LTYPE>
  bool check_nondecreasing_order
    (const ITYPE list[], const LTYPE list_length,
     const char * list_name, IJK::ERROR & error)
  {
    for (LTYPE i = 0; i+1 < list_length; i++) {
      if (list[i] > list[i+1]) {
        error.AddMessage
          ("Error. ", list_name, "[", i, "] > ",
           list_name, "[", i+1, "].");
        error.AddMessage
          ("  ", list_name, "[", i, "] = ", list[i], ".");
        error.AddMessage
          ("  ", list_name, "[", i+1, "] = ", list[i+1], ".");
        return false;
      }
    }

    return true;
  }


  /*!
   *  @overload
   *  @brief Return true if elements are sorted in non-decreasing order. (C++ vector.)
   *  - Version using C++ STL vector for list[].
   */
  template <typename ITYPE>
  bool check_nondecreasing_order
    (const std::vector<ITYPE> & list,
     const char * list_name, IJK::ERROR & error)
  {
    return check_nondecreasing_order
      (IJK::vector2pointer(list), list.size(), list_name, error);
  }


  /*!
    *  @brief Return true if elements are sorted in increasing order.
    *  - Return false and set error message if list[i] >= list[i+1]
    *    for some i.
    */
   template <typename ITYPE, typename LTYPE>
   bool check_increasing_order
     (const ITYPE list[], const LTYPE list_length,
      const char * list_name, IJK::ERROR & error)
   {
     for (LTYPE i = 0; i+1 < list_length; i++) {
       if (list[i] >= list[i+1]) {
         error.AddMessage
           ("Error. ", list_name, "[", i, "] > ",
            list_name, "[", i+1, "].");
         error.AddMessage
           ("  ", list_name, "[", i, "] = ", list[i], ".");
         error.AddMessage
           ("  ", list_name, "[", i+1, "] = ", list[i+1], ".");
         return false;
       }
     }

     return true;
   }


   /*!
    *  @overload
    *  @brief Return true if elements are sorted in increasing order. (C++ vector.)
    *  - Version using C++ STL vector for list[].
    */
   template <typename ITYPE>
   bool check_increasing_order
     (const std::vector<ITYPE> & list,
      const char * list_name, IJK::ERROR & error)
   {
     return check_increasing_order
       (IJK::vector2pointer(list), list.size(), list_name, error);
   }


  /*!
   *  @brief Return true if listA[i] = listB[listA_map[i]].
   *  - Return false and set error message if
   *    listA[i] != listB[listA_map[i]] for some i.
   *  @param listA[] List A.
   *  @param listB[] List B. Permutation of listA[].
   *  @param listA_map[i] Location of listA[i] in listB[].
   *    - listA[i] = listB[listA_map[i]].
   */
  template <typename ITYPE, typename LTYPE, typename MTYPE>
  bool check_list_map
    (const ITYPE listA[], const LTYPE listA_length,
     const ITYPE listB[], const LTYPE listB_length,
     const MTYPE listA_map[],
     const char * listA_name,
     const char * listB_name,
     const char * listA_map_name,
     IJK::ERROR & error)
  {
    for (LTYPE i = 0; i+1 < listA_length; i++) {
      const MTYPE iloc = listA_map[i];
      if ((iloc < 0) || (iloc >= listB_length)) {
        error.AddMessage
          ("Error. Illegal value of ", listA_map_name, "[", i, "].");
        error.AddMessage
          ("  ", listA_map_name, "[", i, "] = ", listA_map[i], ".");
        error.AddMessage
          ("  Should be in range [0..", listB_length-1, "].");
        return false;
      }

      if (listA[i] != listB[iloc]) {
        error.AddMessage
          ("Error. Incorrect value ", listA_map_name, "[", i, "].");
        error.AddMessage
          ("  listA_map[", i, "] = ", iloc, ".");
        error.AddMessage
          ("  ", listA_name, "[", i, "] = ", listA[i], ".");
        error.AddMessage
          ("  ", listB_name, "[", iloc, "] = ", listB[iloc], "." );
        error.AddMessage
          ("  ", listA_name, "[", i, "] != ",
           listB_name, "[", iloc, "].");
        return false;
      }
    }

    return true;
  }


  /*!
   *  @overload
   *  @brief Return true if listA[i] = listB[listA_map[i]]. (C++ vector.)
   *  - Version using C++ STL vectors for all lists.
   */
  template <typename ITYPE, typename MTYPE>
  bool check_list_map
    (const std::vector<ITYPE> & listA,
     const std::vector<ITYPE> & listB,
     const std::vector<MTYPE> & listA_map,
     const char * listA_name,
     const char * listB_name,
     const char * listA_map_name,
     IJK::ERROR & error)
  {
    if (!check_equal_vector_sizes
        (listA, listA_map, listA_name, listA_map_name, error))
      { return false; }

    return check_list_map
      (IJK::vector2pointer(listA), listA.size(),
       IJK::vector2pointer(listB), listB.size(),
       IJK::vector2pointer(listA_map),
       listA_name, listB_name, listA_map_name, error);
  }

  //@}

}

#endif
