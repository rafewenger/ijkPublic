/*!
 *  @file ijkdatatable.tpp
 *  @brief ijk templates defining data tables.
 *  - Version 0.4.0
 *  - Used to create tables (rows/columns) of data.
 *  - Can be used to represent histograms.
 *  - Typically, rows are buckets/binscorresponding 
 *    to ranges of values.
 *  - Columns are measurements/values associated
 *    with each bucket.
 *  - Adding a column representing number of items 
 *    in a bucket gives a histogram.
 *  - Columns may also represent sum/average/min/max
 *    of values associated with items in a bin.
 *  - Function compute_bucket() can be used to map items
 *    to a bucket.
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008-2022 Rephael Wenger

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


#ifndef _IJKDATATABLE_
#define _IJKDATATABLE_

#include <string>

#include "ijk.tpp"

/// Classes for tables of data.
namespace IJKDATATABLE {

  using std::string;

  // **************************************************
  // DATA_TABLE
  // **************************************************

  /*!
   *  @brief Base class for data table.
   *  - Data table consists of rows and columns.
   *  - Each column is of type COLUMN_DATA and is
   *    defined within the derived classes.
   *  - Typically, rows are "bins/buckets" and
   *    columns contain information associated
   *    with each bin/bucket.
   */
  template <typename NTYPE> class DATA_TABLE_BASE {

  public:
    typedef NTYPE NUMBER_TYPE;  ///< Number type.

  protected:
    NTYPE num_rows;  ///< Number of rows in table.

    void Init(const NTYPE num_rows) { this->num_rows = num_rows; };

  public:
    DATA_TABLE_BASE<NTYPE>() { Init(0); };
    DATA_TABLE_BASE<NTYPE>(const NTYPE num_rows) { Init(num_rows); };
    ~DATA_TABLE_BASE<NTYPE>() { num_rows = 0; };

    NTYPE NumRows() const { return(num_rows); };
  };


  /*!
   *  @brief Column of data.
   *  - Note that different columns within the same table may
   *    have different types.
   */
  template <typename NTYPE, typename DTYPE> class DATA_COLUMN {

  public:
    typedef NTYPE NUMBER_TYPE;  ///< Number type.

  protected:
    string label;
    DTYPE * data;
    NTYPE num_rows;

    void Init(const string & label, const NTYPE num_rows);
    void FreeAll();
    void SetNumRows(const NTYPE num_rows);

  public:
    DATA_COLUMN<NTYPE,DTYPE>() { Init("", 0); };
    DATA_COLUMN<NTYPE,DTYPE>(const string & label) { Init(label, 0); };
    DATA_COLUMN<NTYPE,DTYPE>(const string & label, const NTYPE num_rows) 
    { Init(label, num_rows); };

    // set functions
    void Set(const NTYPE irow, const DTYPE v) { data[irow] = v; };
    void SetAll(const DTYPE v);
    template<typename ITYPE>
    void SetAtIntervals(const DTYPE v0, const ITYPE interval);
    void Increment(const NTYPE irow) { data[irow]++; };
    void Add(const NTYPE irow, const DTYPE v) { data[irow] += v; };
    void Add(const NTYPE irow0, const NTYPE irow1, const DTYPE v);
    void AddSquare           ///< Add v^2 to row irow
    (const NTYPE irow, const DTYPE v) { data[irow] += (v*v); };

    // get functions
    DTYPE operator [] (const NTYPE irow) const { return(data[irow]); };
    NTYPE NumRows() const { return(num_rows); };
    string Label() const { return(label); };

    // compute functions
    DTYPE Sum() const;  ///< Compute and return sum of all values in column.

    /*!
     *  @brief Return smallest i such that sum of all values in first i rows times x is at least the sum of all values in all rows.
     *  - Return 0 if number of rows is 0.
     */
    template <typename XTYPE>
    NTYPE kthRow(const XTYPE x) const;

    /*!
     *  @brief Return smallest i such that sum of all values in first i rows is at least half the sum of all values in all rows.
     *  - Return 0 if number of rows is 0.
     */
    NTYPE MedianRow() const;

    /*!
     *  @brief Return smallest i such that sum of all values in first i rows is at least (kth_percentile/100) times the sum of all values.
     *  - Return 0 if number of rows is 0.
     */
    template <typename PERCENTILE>
    NTYPE kthPercentileRow(const PERCENTILE kth_percentile) const;
    
    // data type
    typedef DTYPE DATA_TYPE;
  };

  template <typename NTYPE, typename DTYPE> 
  void DATA_COLUMN<NTYPE,DTYPE>::
  Init(const string & label, const NTYPE num_rows)
  {
    this->label = label;
    this->num_rows = num_rows;
    data = NULL;
    if (num_rows > 0) { data = new DTYPE[num_rows]; };
  };

  template <typename NTYPE, typename DTYPE> 
  void DATA_COLUMN<NTYPE,DTYPE>::FreeAll()
  {
    this->label = "";
    delete [] data;
    data = NULL;
    num_rows = 0;
  }

  /// Set all rows.
  template <typename NTYPE, typename DTYPE> 
  void DATA_COLUMN<NTYPE,DTYPE>::SetAll(const DTYPE v)
  {
    for (NTYPE irow = 0; irow < NumRows(); irow++) 
      { data[irow] = v; }
  }

  /// Set all rows at equally spaced intervals.
  template <typename NTYPE, typename DTYPE> 
  template <typename ITYPE>
  void DATA_COLUMN<NTYPE,DTYPE>::SetAtIntervals
  (const DTYPE v0, const ITYPE interval)
  {
    for (NTYPE irow = 0; irow < NumRows(); irow++) 
      { data[irow] = v0 + irow*interval; }
  }

  template <typename NTYPE, typename DTYPE> 
  void DATA_COLUMN<NTYPE,DTYPE>::
  Add(const NTYPE irow0, const NTYPE irow1, const DTYPE v)
  {
    for (NTYPE irow = irow0; irow <= irow1; irow++)
      { Add(irow, v); };
  }

  template <typename NTYPE, typename DTYPE> 
  void DATA_COLUMN<NTYPE,DTYPE>::SetNumRows(const NTYPE num_rows)
    // set number of rows
  {
    FreeAll();
    this->num_rows = num_rows;
    data = NULL;
    if (num_rows > 0) { data = new DTYPE[num_rows]; };
  }

  /// Compute and return sum of all values in column.
  template <typename NTYPE, typename DTYPE>
  DTYPE DATA_COLUMN<NTYPE,DTYPE>::Sum() const
  {
    DTYPE sum = 0;
    for (NTYPE irow = 0; irow < NumRows(); irow++) 
      { sum += data[irow]; }

    return(sum);
  }

  // Return smallest i such that sum of all values in first i rows times x
  //   is at least the sum of all values in all rows.
  // Return 0 if number of rows is 0.
  template <typename NTYPE, typename DTYPE>
  template <typename XTYPE>
  NTYPE DATA_COLUMN<NTYPE,DTYPE>::kthRow(const XTYPE x) const
  {
    if (num_rows == 0) { return(0); }

    DTYPE sum = Sum();
    NTYPE irow = 0;
    DTYPE sum2 = data[0];
    while (x*sum2 < sum && irow < num_rows) {
      irow++;
      sum2 += data[irow];
    }

    return(irow);
  }

  // Return smallest i such that sum of all values in first i rows 
  //   is at least half the sum of all values in all rows.
  // Return 0 if number of rows is 0.
  template <typename NTYPE, typename DTYPE>
  NTYPE DATA_COLUMN<NTYPE,DTYPE>::MedianRow() const
  {
    return(kthRow(2));
  }

  // Return smallest i such that sum of all values in first (i-1) rows
  //   is less than (kth_percentile/100) times the sum of all values.
  template <typename NTYPE, typename DTYPE>
  template <typename PERCENTILE>
  NTYPE DATA_COLUMN<NTYPE,DTYPE>::
  kthPercentileRow(const PERCENTILE kth_percentile) const
  {
    return(kthRow(100.0/kth_percentile));
  }


// **************************************************
// Compute Bucket
// **************************************************

  /// Return the bucket containing the input value.
  template <typename STYPE0, typename STYPE1, 
            typename ITYPE, typename NTYPE>
  NTYPE compute_bucket(const STYPE0 value, const STYPE1 min_value,
                       const ITYPE interval, const NTYPE num_buckets)
  {
    NTYPE ibucket = NTYPE((value-min_value)/interval);
    if (ibucket >= num_buckets) { ibucket = num_buckets-1; };
    if (ibucket < 0) { ibucket = 0; };

    return(ibucket);
  }

  /// Return the number of buckets in the given range.
  template <typename NTYPE, typename STYPE, typename ITYPE>
  NTYPE compute_num_buckets(const STYPE min_value, const STYPE max_value,
			    const ITYPE interval)
  {
    NTYPE max_bucket = NTYPE((max_value-min_value)/interval);
    NTYPE num_buckets = max_bucket+1;
    if (num_buckets <= 0) { num_buckets = 1; };

    return(num_buckets);
  }

// **************************************************
// Check
// **************************************************

  /// Check data columns
  /// Return true if columns have same number of rows.
  template <typename NTYPE0, typename DTYPE0, 
            typename NTYPE1, typename DTYPE1>
  bool check_data_columns (const DATA_COLUMN<NTYPE0,DTYPE0> & col0,
			   const DATA_COLUMN<NTYPE1,DTYPE1> & col1,
			   IJK::ERROR & error)
  {
    if (col0.NumRows() != col1.NumRows()) {
      error.AddMessage("Data columns have different numbers of rows.");
      error.AddMessage("One column has ", col0.NumRows(), " rows.");
      error.AddMessage("Another column has ", col1.NumRows(), " rows.");
      return(false);
    }
    return(true);
  }
}

#endif
