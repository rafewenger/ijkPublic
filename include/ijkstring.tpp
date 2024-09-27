/*!
 *  @file ijkstring.tpp
 *  @brief ijk templates for converting strings to values/arrays
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2016 Rephael Wenger

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

#ifndef _IJKSTRING_
#define _IJKSTRING_

#include <cstring>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>


namespace IJK {


  // *****************************************************************
  //! @name CONVERT STRING TO NUMERIC VALUES/ARRAYS
  // *****************************************************************

  namespace {

    /// Remove trailing blanks from string s
    inline void remove_trailing_blanks(std::string & s)
    {
      size_t pos = 0;
      for (size_t i = 0; i < s.length(); i++) {
        if (!isspace(s[i])) { pos = i+1; }
      }
      if (pos < s.length()) { s.erase(pos); };
        
    }
  };

  ///@{

  /// Convert string to value
  template <typename VTYPE>
  bool string2val(const char * s, VTYPE & val)
  {
    std::istringstream v_string;

    std::string s2 = s;
    remove_trailing_blanks(s2);

    v_string.str(s2);

    v_string >> val;

    if (!v_string.eof()) 
      { return(false); }
    else
      { return(true); }
  }


  /// Convert string to array of elements and append to vector.
  template <typename ETYPE>
  bool string2vector_append(const char * s, std::vector<ETYPE> & v)
  {
    std::istringstream v_string;

    std::string s2 = s;
    remove_trailing_blanks(s2);
    if (s2.size() == 0) {
      // Empty string. Vector v has length 0.
      return(true);
    }

    v_string.str(s2);
    while (v_string.good()) {
      ETYPE x;
      v_string >> x;
      v.push_back(x);
    }

    if (!v_string.eof()) 
      { return(false); }
    else
      { return(true); }
  }

  /// Convert string to array of elements
  template <typename ETYPE>
  bool string2vector(const char * s, std::vector<ETYPE> & v)
  {
    v.clear();
    return(string2vector_append(s, v));
  }

  ///@}


  // *****************************************************************
  //! @name CONVERT NUMERIC VALUE OR BOOL TO STRING
  // *****************************************************************

  /// @{

  /// @brief Convert value to string.
  /// - Return false if error in converting to string.
  template <typename T>
  bool val2string(const T x, std::string & s_out)
  {
    std::ostringstream s_stream;

    s_stream << x;

    if (s_stream.bad()) { 
      s_out.clear();
      return(false); 
    }

    s_out = s_stream.str();
    return(true);
  }

  /// @brief Convert array to string.
  /// - Return false if error in converting to string.
  template <typename T, typename I>
  bool array2string
  (const T x[], const I length, const char * separator, std::string & s_out)
  {
    std::ostringstream s_stream;

    if (length <= 0) { return(true); };

    s_stream << x[0];
    if (s_stream.bad()) { 
      s_out.clear();
      return(false); 
    };
    
    for (I i = 1; i < length; i++) {
      s_stream << separator << x[i];
      if (s_stream.bad()) {
        s_out.clear();
        return(false); 
      }
    }

    s_out = s_stream.str();
    return(true);
  }

  /// @brief Convert C++ vector to string.
  /// - Return false if error in converting to string.
  template <typename T, typename I>
  bool vector2string
  (const std::vector<T> & x, const char * separator, std::string & s_out)
  {
    if (x.size() <= 0) { return(true); };

    return(array2string(&(x.front()), x.size(), separator, s_out));
  }


  /// @brief Convert bool to string true or false.
  /// - Return string.
  /// - Mainly useful for printing bool value.
  inline std::string bool2string(const bool flag)
  {
    if (flag) { return std::string("true"); }
    else { return std::string("false"); }
  }


  /// @brief Create string from prefix and value.
  /// - Return false if error in creating string.
  template <typename T>
  bool create_string(const std::string & prefix, const T x, std::string & s_out)
  {
    const bool flag = IJK::val2string(x, s_out);
    if (flag)
      { s_out = prefix + s_out; }
    else {
      // Set s_out to prefix.
      s_out = prefix;
    }

    return flag;
  }


  /// @brief Create string from prefix and value.
  /// - Return false if error in creating string.
  /// - Version with const char * prefix.
  template <typename T>
  bool create_string(const char * prefix, const T x, std::string & s_out)
  {
    return create_string(std::string(prefix), x, s_out);
  }


  ///@}


  // *****************************************************************
  //! @name SPLIT STRING INTO PREFIX AND SUFFIX.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Split string at last occurrence of character c into prefix and suffix.
   *  @param s Input string.
   *  @param c Split at character c.
   *  @param[out] prefix Prefix. All characters before last occurrence of c.
   *    - Does not include last occurrence of character c.
   *  @param[out] suffix Suffix. All characters after last occurrence of c.
   *    - Does not include last occurrence of character c.
   */
  template <typename STRING_TYPE>
  void split_string(const STRING_TYPE & s, const char c,
                    STRING_TYPE & prefix, STRING_TYPE & suffix)
  {
    typename STRING_TYPE::size_type i = s.rfind(c);
    if (i == STRING_TYPE::npos) {
      prefix = s;
      suffix = "";
    }
    else {
      if (i > 0) { prefix = s.substr(0,i); }
      else { prefix = ""; };

      if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
      else { suffix = ""; };
    }
  }


  /*!
   *  @brief Split string into prefix and suffix.
   *  @param s Input string. (Type char *.)
   *  - Version with input string of type char *.
   */
  template <typename STRING_TYPE>
  void split_string(const char * s, const char c,
                    STRING_TYPE & prefix, STRING_TYPE & suffix)
  {
    split_string(STRING_TYPE(s), c, prefix, suffix);
  }

  ///@}


  // *****************************************************************
  //! @name GET FILENAME
  // *****************************************************************

  /*!
   *  @brief Get filename from pathname.
   *  - If pathname ends in ".{suffix}", remove ".{suffix}.
   */
  inline void get_filename_remove_suffix
  (const char * pathname, const char * suffix, std::string & filename)
  {
    // remove path from file name
    std::string prefix0, suffix0;
    split_string(pathname, '/', prefix0, suffix0);
    if (suffix0 == "") 
      { filename = pathname; }
    else
      { filename = suffix0; }

    split_string(filename, '.', prefix0, suffix0);
    if (suffix == suffix0) 
      { filename = prefix0; }
  }


  // *****************************************************************
  /// @name QUERY STRING
  // *****************************************************************

  ///@{

  /// Return true if string represents object of given type.
  template <typename T>
  bool is_type(const std::string s) 
  {
    std::istringstream s_stream(s);
    T x;

    s_stream >> x;
    return (s_stream.eof() && !s_stream.fail());
  }


  /// Return true if suffix is a suffix of s.
  inline bool is_suffix(const std::string & s, const std::string & suffix)
  {
    typedef std::string::size_type SIZE_TYPE;

    const SIZE_TYPE s_len = s.size();
    const SIZE_TYPE suffix_len = suffix.size();

    if (s_len < suffix_len)
      { return(false); }

    SIZE_TYPE k0 = s_len - suffix_len;
    for (SIZE_TYPE k = k0; k < s_len; k++) {
      if (s[k] != suffix[k-k0]) { return(false); };
    }
  
    return(true);
  }


  /*!
   *  @brief Return true if suffix is a suffix of s.
   *  - Version where s has type char *.
   */
  inline bool is_suffix(const char * s, const std::string & suffix)
  {
    if (s == NULL) { return(false); }
    return(is_suffix(std::string(s), suffix));
  }

  /// Return true if suffix is a suffix of s.
  /// Version where both s and suffix have type char *.
  inline bool is_suffix(const char * s, const char * suffix)
  {
    if (s == NULL) { return(false); }
    return(is_suffix(std::string(s), std::string(suffix)));
  }

  ///@}


  // *****************************************************************
  //! @name MAP ENUM VALUE TO STRING
  // *****************************************************************


  /// Map of enum types to (char *) strings.
  template <typename ENUM_TYPE>
  class ENUM_STRING_MAP:
    public std::unordered_map<ENUM_TYPE, const char *>
  {};


  // *****************************************************************
  //! @name *** DEPRECATED *** GET SUFFIX TYPE
  // *****************************************************************

  /// @brief *** DEPRECATED *** List of pairs of suffix types and suffix strings.
  template <typename SUFFIX_TYPE>
  class SUFFIX_TYPE_LIST:
    public std::vector< std::pair<SUFFIX_TYPE, const char *> >
  {};


  /*!
   *  @brief Return type of suffix of string s.
   *  - Version where string s is a C++ STL string.
   *  @tparam SUFFIX_TYPE An enum type.
   *  @param suffix_list[] List of pairs of suffix types and suffix strings.
   *  @param unknown_suffix_type Return type unknown_suffix_type if
   *    string s does not have suffix in suffix_list.
   */
  template <typename SUFFIX_TYPE>
  SUFFIX_TYPE get_suffix_type
  (const std::string & s, const SUFFIX_TYPE_LIST<SUFFIX_TYPE> & suffix_list,
   const SUFFIX_TYPE unknown_suffix_type)
  {
    typedef typename SUFFIX_TYPE_LIST<SUFFIX_TYPE>::size_type SIZE_TYPE;

    for (SIZE_TYPE i = 0; i < suffix_list.size(); i++) {
      if (is_suffix(s, suffix_list[i].second))
        { return(suffix_list[i].first); }    
    }

    return(unknown_suffix_type);
  }


  /*!
   *  @brief Return type of suffix of string s.
   *  - Version where string s has type (char *).
   *  @tparam SUFFIX_TYPE An enum type.
   *  @param s String of type (char *).
   *  @param suffix_list[] List of pairs of suffix types and suffix strings.
   *  @param unknown_suffix_type Return type unknown_suffix_type if
   *    string s does not have suffix in suffix_list.
   */
  template <typename SUFFIX_TYPE>
  SUFFIX_TYPE get_suffix_type
  (const char * s, const SUFFIX_TYPE_LIST<SUFFIX_TYPE> & suffix_list,
   const SUFFIX_TYPE unknown_suffix_type)
  {
    if (s == NULL) { return(unknown_suffix_type); }

    return(get_suffix_type(std::string(s), suffix_list, unknown_suffix_type));
  }


  /*!
   *  @brief *** DEPRECATED *** Return suffix string for type suffix_type.
   *  @tparam SUFFIX_TYPE An enum type.
   *  @param suffix_list[] List of pairs of suffix types and suffix strings.
   *  @param unknown_type_string Return string unknown_type_string if
   *    suffix_type is not in suffix_list.
   */
  template <typename SUFFIX_TYPE>
  const char * get_suffix_string
  (const SUFFIX_TYPE suffix_type, 
   const SUFFIX_TYPE_LIST<SUFFIX_TYPE> & suffix_list,
   const char * unknown_type_string)
  {
    typedef typename SUFFIX_TYPE_LIST<SUFFIX_TYPE>::size_type SIZE_TYPE;

    for (SIZE_TYPE i = 0; i < suffix_list.size(); i++) {
      if (suffix_type == suffix_list[i].first)
        { return(suffix_list[i].second); }    
    }

    return(unknown_type_string);
  }


}

#endif
