/*!
 *  @file ijkhash.tpp
 *  @brief ijk templates for hashing.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2018-2022 Rephael Wenger

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

#ifndef _IJKHASH_
#define _IJKHASH_

#include <unordered_map>
#include <utility>


namespace IJK {

  template <typename T>
  inline std::size_t change_hash
  (const T key, std::size_t & hash_value)
  {
    hash_value ^= key + 0x9e3779b9 + (hash_value<<6) + (hash_value>>2);

    return(hash_value);
  }

  /// Hash two arguments (a,b).
  template <typename TYPEA, typename TYPEB>
  std::size_t hash_func(const TYPEA a, const TYPEB b)
  {
    std::size_t hash_value = a + 0x9e3779b9;
    change_hash(b, hash_value);

    return(hash_value);
  }

  /// Hash pair (a,b).
  /// @tparam TYPEA Type of first element in key.
  /// @tparam TYPEB Type of second element in key.
  template <typename TYPEA, typename TYPEB>
  std::size_t hash_func(const std::pair<TYPEA,TYPEB> & key)
  {
    return(hash_func(key.first, key.second));
  }

  template <typename TYPEA, typename TYPEB>
  struct HASH_PAIR {

    HASH_PAIR(){};

    size_t operator() (const std::pair<TYPEA,TYPEB> & key) const
    { return(hash_func(key)); }
  };

  template <typename TYPEA, typename TYPEB, typename DATA_TYPE>
  class unordered_map_key_pair:
    public std::unordered_map
  <std::pair<TYPEA,TYPEB>, DATA_TYPE, HASH_PAIR<TYPEA,TYPEB> > 
  {
    
  protected:
    typedef typename std::unordered_map
    <std::pair<TYPEA,TYPEB>, int, HASH_PAIR<TYPEA,TYPEB> > BASE_CLASS;

  public:

    typedef typename BASE_CLASS::iterator iterator;
    typedef typename BASE_CLASS::const_iterator const_iterator;

    unordered_map_key_pair(){};

    std::pair<const_iterator,bool> 
    insert (const TYPEA a, const TYPEB b, const DATA_TYPE data)
    {
      return(BASE_CLASS::insert
             (std::make_pair( std::make_pair(a,b), data)));
    }

    iterator find(const TYPEA a, const TYPEB b)
    { return(BASE_CLASS::find(std::make_pair(a,b))); }

    const_iterator find(const TYPEA a, const TYPEB b) const
    { return(BASE_CLASS::find(std::make_pair(a,b))); }
  };

}

#endif

