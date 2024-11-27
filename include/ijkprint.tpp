/*!
 *  @file ijkprint.tpp
 *  @brief Templates for printing lists, coordinates, etc.
 *  - Useful for debugging or for routines printing information
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2023 Rephael Wenger

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

#ifndef _IJKPRINT_
#define _IJKPRINT_

#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>


#include "ijk.tpp"

namespace IJK {

  // *****************************************************************
  // PRINT LIST
  // *****************************************************************

  /// Print values in list.  No parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list_values
  (std::ostream & out, const ETYPE * list, const NTYPE length,
   const char separator)
  {
    for (NTYPE i = 0; i < length; i++) {
      if (i > 0) { out << separator; };
      out << list[i];
    }
  }

  
  /// Print list with left and right delimiters.
  template <typename ETYPE, typename NTYPE>
  void print_list(std::ostream & out, const ETYPE * list, const NTYPE length,
                  const char separator,
                  const char left_delim, const char right_delim)
  {
    out << left_delim;
    print_list_values(out, list, length, separator);
    out << right_delim;
  }

  
  /// @overload
  /// Print list enclosed in parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list(std::ostream & out, const ETYPE * list, const NTYPE length,
                  const char separator)
  {
    print_list(out, list, length, separator, '(', ')');
  }


  /// @overload
  /// Print list separated by commas and enclosed in parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list(std::ostream & out, const ETYPE * list, const NTYPE length)
  {
    print_list(out, list, length, ',');
  }


  /// @overload
  /// Print list preceded and followed by char strings.
  template <typename ETYPE, typename NTYPE>
  void print_list
  (std::ostream & out, const char * s1, 
   const ETYPE * list, const NTYPE length, const char * s2)
  {
    out << s1;
    print_list(out, list, length, ',');
    out << s2;
  }


  /// @overload
  /// @brief Print list separated by commas and enclosed in parentheses.
  /// - Version using C++ STL vector list[].
  template <typename ETYPE>
  void print_list(std::ostream & out, const std::vector<ETYPE> & list)
  {
    print_list(out, IJK::vector2pointer(list), list.size());
  }


  /// @overload
  /// @brief Print list preceded and followed by char strings.
  /// - Version using C++ STL vector list[].
  template <typename ETYPE>
  void print_list(std::ostream & out, const char * s1,
                  const std::vector<ETYPE> & list, const char * s2)
  {
    print_list(out, s1, IJK::vector2pointer(list), list.size(), s2);
  }


  // *****************************************************************
  // PRINT LIST VALUES AS INTEGERS
  // *****************************************************************

  /// Print values in list as integers. (Cast to integers.)  No parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list_values_as_int
  (std::ostream & out, const ETYPE * list, const NTYPE length,
   const char separator)
  {
    for (NTYPE i = 0; i < length; i++) {
      if (i > 0) { out << separator; };
      out << int(list[i]);
    }
  }

  
  /// Print list values as integers with left and right delimiters.
  template <typename ETYPE, typename NTYPE>
  void print_list_as_int
  (std::ostream & out, const ETYPE * list, const NTYPE length,
   const char separator,
   const char left_delim, const char right_delim)
  {
    out << left_delim;
    print_list_values_as_int(out, list, length, separator);
    out << right_delim;
  }


  /// @overload
  /// @brief Print list values as integers enclosed in parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list_as_int
  (std::ostream & out, const ETYPE * list, const NTYPE length,
   const char separator)
  {
    print_list_as_int(out, list, length, separator, '(', ')');
  }


  /// @overload
  /// @brief Print list values as int separated by commas
  ///   and enclosed in parentheses.
  template <typename ETYPE, typename NTYPE>
  void print_list_as_int
  (std::ostream & out, const ETYPE * list, const NTYPE length)
  {
    print_list_as_int(out, list, length, ',');
  }

  
  /// @overload
  /// @brief Print list values as int separated by commas
  ///   and enclosed in parentheses. (C++ STL vector.)
  /// - Version using C++ STL vector list[].
  template <typename ETYPE>
  void print_list_as_int(std::ostream & out, const std::vector<ETYPE> & list)
  {
    print_list_as_int(out, IJK::vector2pointer(list), list.size());
  }


  /// @overload
  /// @brief Print list values as int, preceded and
  ///   followed by char strings.
  template <typename ETYPE, typename NTYPE>
  void print_list_as_int
  (std::ostream & out, const char * s1, 
   const ETYPE * list, const NTYPE length, const char * s2)
  {
    out << s1;
    print_list_as_int(out, list, length, ',');
    out << s2;
  }

  
  /// @overload
  /// @brief Print list values as int, precede and
  ///   followed by char strings. (C++ STL vector.)
  /// - Version using C++ STL vector list[].
  template <typename ETYPE>
  void print_list_as_int
  (std::ostream & out, const char * s1,
   const std::vector<ETYPE> & list, const char * s2)
  {
    print_list_as_int
      (out, s1, IJK::vector2pointer(list), list.size(), s2);
  }

  
  // *****************************************************************
  // PRINT LIST OF TUPLES
  // *****************************************************************

  /// Print list of tuples.
  template <typename ETYPE, typename N0TYPE, typename N1TYPE>
  void print_list_of_tuples
  (std::ostream & out, const ETYPE * list, 
   const N0TYPE tuple_size, const N1TYPE num_tuples,
   const char separator0=',', const char separator1=' ',
   const char left_delim='(', const char right_delim=')')
  {
    for (N1TYPE i = 0; i < num_tuples; i++) {
      N1TYPE k = i*tuple_size;
      print_list
        (out, list+k, tuple_size, separator0, left_delim, right_delim);
      if (i+1 < num_tuples)
        { out << separator1; }
    }
  }

  /// Print list of tuples.
  template <typename ETYPE, typename NTYPE>
  void print_list_of_tuples
  (std::ostream & out, const std::vector<ETYPE> & list, 
   const NTYPE tuple_size)
  {
    typedef typename std::vector<ETYPE>::size_type SIZE_TYPE;

    if (list.size() <= 0) { return; }

    SIZE_TYPE num_tuples = list.size()/tuple_size;
    print_list_of_tuples(out, &(list[0]), tuple_size, num_tuples);
  }


  // *****************************************************************
  // PRINT 3D COORDINATES
  // *****************************************************************

  /// Print 3 coordinates.
  template <typename CTYPE>
  void print_coord3D(std::ostream & out, const CTYPE * coord)
  {
    print_list(out, coord, 3);
  }

  /// Print 3 coordinates.
  /// - Version using C++ STL vector coord[].
  template <typename CTYPE>
  void print_coord3D(std::ostream & out, 
                     const std::vector<CTYPE> & coord)
  {
    print_coord3D(out, IJK::vector2pointer(coord));
  }

  /// Print 3 coordinates followed by a string.
  template <typename CTYPE>
  void print_coord3D(std::ostream & out, const CTYPE * coord, const char * s)
  {
    print_coord3D(out, coord);
    out << s;
  }

  /// Print 3 coordinates preceded and followed by a string.
  template <typename CTYPE>
  void print_coord3D
  (std::ostream & out, const char * s0, const CTYPE * coord, const char * s1)
  {
    out << s0;
    print_coord3D(out, coord);
    out << s1;
  }

  /// Print 3 coordinates preceded and followed by a string.
  /// - Version using C++ STL vector coord[].
  template <typename CTYPE>
  void print_coord3D
  (std::ostream & out, const char * s0, 
   const std::vector<CTYPE> & coord, const char * s1)
  {
    out << s0;
    print_coord3D(out, coord);
    out << s1;
  }

  /// Print two pairs of 3 coordinates preceded, separated and 
  ///   followed by a string.
  template <typename CTYPE>
  void print_coord3D
  (std::ostream & out, const char * s0, const CTYPE * coord0, const char * s1,
   const CTYPE * coord1, const char * s2)
  {
    out << s0;
    print_coord3D(out, coord0);
    out << s1;
    print_coord3D(out, coord1);
    out << s2;
  }


  // *****************************************************************
  // PRINT GRID COORDINATES
  // *****************************************************************

  template <typename GTYPE, typename VTYPE>
  void print_grid_coord(std::ostream & out, const GTYPE & grid, const VTYPE iv)
  {
    typename GTYPE::DIMENSION_TYPE dimension = grid.Dimension();
    VTYPE coord[dimension];

    grid.ComputeCoord(iv, coord);
    print_list(out, coord, dimension);
  }

  template <typename GTYPE, typename VTYPE>
  void print_grid_coord(std::ostream & out, const GTYPE & grid, const VTYPE iv,
                        const char * s)
  {
    typename GTYPE::DIMENSION_TYPE dimension = grid.Dimension();
    VTYPE coord[dimension];

    grid.ComputeCoord(iv, coord);
    print_list(out, coord, dimension);
    out << s;
  }


  // *****************************************************************
  // PRINT GRID SIZE
  // *****************************************************************

  template <typename GTYPE>
  void print_grid_size(std::ostream & out, const GTYPE & grid, const char * s)
  {
    if (grid.Dimension() < 1) 
      { out << "Grid dimension: " << grid.Dimension(); }
    else {
      out << "Grid ";
      for (int d = 0; d < grid.Dimension(); d++) {
        out << grid.AxisSize(d);
        if (d+1 < grid.Dimension()) 
          { out << "x"; }
      };
    }
    out << s; 
  }


  // *****************************************************************
  // PRINT POLYGON VERTICES
  // *****************************************************************

  /// Print triangle vertices, preceded and followed by a string.
  template <typename VTYPE>
  void print_triangle_vertices
  (std::ostream & out, const char * s0, const VTYPE tri_vert[],
   const char * s1)
  {
    const int NUM_VERT_PER_TRI(3);

    out << s0;
    print_list(out, tri_vert, NUM_VERT_PER_TRI);
    out << s1;
  }

  /// Print quadrilateral vertices, preceded and followed by a string.
  template <typename VTYPE>
  void print_quad_vertices
  (std::ostream & out, const char * s0, const VTYPE quad_vert[],
   const char * s1)
  {
    const int NUM_VERT_PER_QUAD(4);

    out << s0;
    print_list(out, quad_vert, NUM_VERT_PER_QUAD);
    out << s1;
  }

  /// Print pentagon vertices, preceded and followed by a string.
  template <typename VTYPE>
  void print_pentagon_vertices
  (std::ostream & out, const char * s0, const VTYPE pentagon_vert[],
   const char * s1)
  {
    const int NUM_VERT_PER_PENTAGON(5);

    out << s0;
    print_list(out, pentagon_vert, NUM_VERT_PER_PENTAGON);
    out << s1;
  }


  /// Print hexahedron vertices, preceded and followed by a string.
  template <typename VTYPE>
  void print_hexahedron_vertices
  (std::ostream & out, const char * s0, const VTYPE hexahedron_vert[],
   const char * s1)
  {
    const int NUM_VERT_PER_HEXAHEDRON(8);

    out << s0;
    print_list(out, hexahedron_vert, NUM_VERT_PER_HEXAHEDRON);
    out << s1;
  }


  // *****************************************************************
  // PRINT TIME
  // *****************************************************************

  inline void print_time(std::ostream & out, const char * s, 
                         const clock_t & t_start, const clock_t & t_end)
  {
    out << s;
    out << float(t_end-t_start)/CLOCKS_PER_SEC;
    out << std::endl;
  }


  // *****************************************************************
  // PRINT PERCENTAGE
  // *****************************************************************

  /// Convert (numerator/denominator) to percentage.
  /// @pre denominator is not zero.
  template <typename T0, typename T1, typename T2>
  void compute_percent
  (const T0 numerator, const T1 denominator, T2 & percent)
  {
    if (denominator == 0) {
      percent = 100*numerator;
    }
    else {
      percent = 100*(double(numerator)/denominator);
    }
  }

  /// Print percentage(if denominator is not zero.)
  /// Percentage is rounded down to nearest integer.
  template <typename T0, typename T1>
  void print_percent
  (std::ostream & out, const T0 numerator, const T1 denominator)
  {
    if (denominator != 0) {
      float percent;
      compute_percent(numerator, denominator, percent);
      percent = std::floor(percent);
      out << "(" << percent << "%)";
    }
  }

}

#endif
