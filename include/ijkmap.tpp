/*!
 *  @file ijkmap.tpp
 *  @brief ijk templates for handling mappings and permutations.
 *  - Version 0.4.0
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

#ifndef _IJKMAP_
#define _IJKMAP_

#include "ijk.tpp"

namespace IJK {


  /*!
   *  @brief Invert the permutation given by mapA[].
   *  @param mapA[] Element i goes to mapA[i].
   *  @param[out] mapB[] Mapping such that mapB[mapA[i]] = i.
   *     @pre Array mapB[] is preallocated to length 
   *        at least num_elements.
   *  @param num_elements Number of elements in mapA[] and mapB[].
   */
  template <typename TYPEA, typename TYPEB, typename NTYPE>
  void invert_permutation
  (const TYPEA * mapA, TYPEB * mapB, const NTYPE num_elements)
  {
    for (NTYPE i = 0; i < num_elements; i++)
      { mapB[mapA[i]] = i; }
  }


  /*!
   *  @overload
   *  @brief Invert the permutation given by mapA[]. (C++ STL vector.)
   *  - Version using C++ STL vector for mapB[].
   *  @param mapA[] Element i goes to mapA[i].
   *  @param[out] mapB[] Mapping such that mapB[mapA[i]] = i.
   *  @param num_elements Number of elements in mapA[] and mapB[].
   */
  template <typename TYPEA, typename TYPEB, typename NTYPE>
  void invert_permutation
  (const TYPEA * mapA, std::vector<TYPEB> & mapB, const NTYPE num_elements)    
  {
    mapB.resize(num_elements);

    invert_permtutation(mapA, IJK::vector2pointer(mapB), num_elements);
  }

  
  /*!
   *  @overload
   *  @brief Invert the permutation given by mapA[]. (C++ STL vector.)
   *  - Version using C++ STL vector for mapA[] and mapB[].
   *  @param mapA[] Element i goes to mapA[i].
   *  @param[out] mapB[] Mapping such that mapB[mapA[i]] = i.
   */  
  template <typename TYPEA, typename TYPEB, typename NTYPE>
  void invert_permutation
  (const std::vector<TYPEA> & mapA, std::vector<TYPEB> & mapB)
  {
    invert_permtutation
      (IJK::vector2pointer(mapA), mapB, mapA.size());
  }
}

#endif

