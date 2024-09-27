/*!
 *  @file ijkenum.tpp
 *  @brief Support for enum types.
 *  - Support for strings for printing enum type.
 *  - Version 0.1.0
 */


#ifndef _IJKENUM_
#define _IJKENUM_

#include <string>

#include "ijk.tpp"

namespace IJK {

  // *****************************************************************
  // CLASSES
  // *****************************************************************

  /*!
   *  @brief Structure storing enum values and string reprentation.
   *  - Intended for:
   *    static ENUM_STR<ENUM_TYPE>> enum_str_pair[] = 
   *      {{ {A,"AAA"}, {B,"BBB"}... }}
   */
  template <typename _ENUM_TYPE>
  class ENUM_STR {
    
  public:
    typedef _ENUM_TYPE ENUM_TYPE;

  public:
    /// @brief Enum value.
    const ENUM_TYPE value;

    /// @brief String representing enum value.
    const char * str;
  };


  /// @brief List of enum values and associated strings.
  template <typename _ENUM_TYPE>
  class ENUM_LIST {

  public:
    typedef _ENUM_TYPE ENUM_TYPE;

  protected:
    std::vector<std::string> str;

    /// @brief Value returned if can't match shape.
    ENUM_TYPE undefined_value;

    void Init(const ENUM_TYPE undefined,
	      const ENUM_STR<ENUM_TYPE> enum_str_pair[],
	      const int num_pairs);

  public:

    /// @brief Constructor with C array and array length.
    /// @param undefined Return undefined if can't match shape.
    ENUM_LIST(const ENUM_TYPE undefined,
	      const ENUM_STR<ENUM_TYPE> enum_str_pair[],
	      const int num_pairs)
    { Init(undefined, enum_str_pair, num_pairs); }
    
    /// @brief Constructor with C++ container that supports arrays,
    ///   either std::array or std::vector.
    /// @param undefined Return undefined if can't match shape.
    template <typename ATYPE>
    ENUM_LIST(const ENUM_TYPE undefined,
	      const ATYPE & enum_str_pair)
    { Init(undefined, enum_str_pair.data(), enum_str_pair.size()); }
    

    /// @brief Return number of enum values that are represented
    ///   in this data structure.
    /// - Note: Some enum values may have empty strings.
    int NumEnum() const
    { return str.size(); }

    /// @brief Return enum string.
    /// @pre int(enum_value) < NumEnum().
    std::string String(const ENUM_TYPE enum_value) const
    { return str[int(enum_value)]; }

    /// @brief Return enum undefined value.
    ENUM_TYPE UndefinedValue() const
    { return undefined_value; }

    /// @brief Return true if enum_value is undefined.
    bool IsUndefined(const ENUM_TYPE enum_value) const
    { return (enum_value == UndefinedValue()); }

    /// @brief Return true if enum_value is defined.
    /// - Equivalent to "!IsUndefined(enum_value)".
    bool IsDefined(const ENUM_TYPE enum_value) const
    { return !IsUndefined(enum_value); }

    /// @brief Return enum value named by enum_str.
    /// - All cases ignored (converted to lower case).
    ENUM_TYPE EnumValue(const std::string & enum_str) const;


    // *** Helpful functions ***

    /// @brief Return string in lower case.
    static std::string ToLower(const std::string & s)
    {
      std::string s2 = s;
      for (int i = 0; i < s2. size(); i++)
        { s2[i] = std::tolower(s2[i]); }
      return s2;
    };

  };


  // *****************************************************************
  // MEMBER FUNCTIONS
  // *****************************************************************

  template <typename _ENUM_TYPE>
  void ENUM_LIST<_ENUM_TYPE>::Init
  (const ENUM_TYPE undefined,
   const ENUM_STR<ENUM_TYPE> enum_str_pair[],
   const int num_pairs)
  {
    IJK::PROCEDURE_ERROR error("ENUM_LIST::Init().");

    for (int i = 0; i < num_pairs; i++) {
      const ENUM_TYPE enum_value  = enum_str_pair[i].value;
      const int ivalue = int(enum_value);

      // Resize C++ vector, so that str[undefined] is allocated.
      str.resize(int(undefined)+1);

      // Store undefined value.
      undefined_value = undefined;

      // Resize C++ vector str, if necessary.
      if (ivalue >= str.size()) 
        { str.resize(ivalue+1); }

      if (str[ivalue] != "") {
        error.AddMessage
          ("Programming error. Attempt to initialize an enum string twice.");
        error.AddMessage
          ("  Attempting to initialize string for enum value ", ivalue,
           " with string ", enum_str_pair[i].str, ".");
        error.AddMessage
          ("  String for enum value already set to ", str[ivalue], ".");
        throw error;
      }

      str[ivalue] = enum_str_pair[i].str;
    }
  }


  // Return enum value named by enum_str.
  template <typename _ENUM_TYPE>
  _ENUM_TYPE ENUM_LIST<_ENUM_TYPE>::
  EnumValue(const std::string & enum_str) const
  {
    const std::string strL = ToLower(enum_str);
    for (int i = 0; i < NumEnum(); i++) {
      const std::string str2L = ToLower(String(ENUM_TYPE(i)));
      if (strL == str2L) 
        { return _ENUM_TYPE(i); }
    }
    return undefined_value;
  };

};

#endif
