/*!
 *  @file ijkrandom.tpp
 *  @brief ijk template classes for random number generation.
 *  - Mainly for testing by generating random isosurface
 *    vertex positions.
 *  - Version 0.4.0
*/


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2022-2023 Rephael Wenger

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


#ifndef _IJKRANDOM_TPP_
#define _IJKRANDOM_TPP_

#include <random>

#include "ijk.tpp"
#include "ijkcoord.tpp"

namespace IJK {

  // *****************************************************************
  //! @name Generate random number in range [0,1].
  // *****************************************************************

  ///@{

  /*!
   *  @brief Abstract base class for generating double in range [0,1].
   */
  class GENERATE_RANDOM01_BASE {

  public:

    /*!
     *  Return random double in range [0,1].
     *  - Random01() must be defined in derived classes.
     *  - Distribution is determined in derived classes.
     *  @param random_engine Pseudorandom number generator.
     *    - Simplest type is std::minstd_rand.
     */
    template <typename RANDOM_ENGINE>
    double Random01(RANDOM_ENGINE & random_engine);
  };


  /*!
   *  @brief Generate double in range [0,1] with uniform distribution.
   */
  class GENERATE_UNIFORM_RANDOM01:
    public GENERATE_RANDOM01_BASE {

  protected:
    std::uniform_real_distribution<double> random_uniform_real;

  public:

    /// Constructor
    GENERATE_UNIFORM_RANDOM01():random_uniform_real(0.0,1.0)
    {}


    /*!
     *  Return random double in range [0,1] with uniform distribution.
     *  @param random_engine Pseudorandom number generator.
     *    - Simplest type is std::minstd_rand.
     */
    template <typename RANDOM_ENGINE>
    double Random01(RANDOM_ENGINE & random_engine)
    {
      return random_uniform_real(random_engine);
    }

  };



  /*!
   *  @brief Generating random double in range [0,1] with U-quadratic distribution.
   *  - Generating function: f(u) = ((u/4) - (1/8))^(1/3) + (1/2)
   *    where u is generated with the uniform distribution.
   *  - The generating function is the inverse of the CDF:
   *       F(x) = 4((x - (1/2))^3 + (1/8)).
   */
  class GENERATE_U_QUADRATIC_RANDOM01: 
    public GENERATE_RANDOM01_BASE {

  protected:
    std::uniform_real_distribution<double> random_uniform_real;

  public:

    /// Constructor
    GENERATE_U_QUADRATIC_RANDOM01():random_uniform_real(0.0,1.0)
    {}

    /*!
     *  Return random double in range [0,1] with U-quadratic distribution.
     *  @param random_engine Pseudorandom number generator.
     *    - Simplest type is std::minstd_rand.
     */
    template <typename RANDOM_ENGINE>
    double Random01(RANDOM_ENGINE & random_engine)
    {
      const double xu = random_uniform_real(random_engine);
      const double xq = std::cbrt((xu/4.0) - (0.125)) + 0.5;

      // Clamp to [0,1] just in case.
      const double xq2 = IJK::clamp_coord_to_range(xq, 0.0, 1.0);

      return xq2;
    }

  };

  ///@}


  // *****************************************************************
  //! @name Generate random number in arbitrary range.
  // *****************************************************************

  /*!
   *  @brief Abstract base class for generating random positions.
   *  @tparam GENERATE_RANDOM01_TYPE 
   *    Class with definition of function Random01.
   */
  template <typename GENERATE_RANDOM01_TYPE>
  class GENERATE_RANDOM_BASE:public GENERATE_RANDOM01_TYPE {

  protected:

    /*!
     *  @brief Map x in range [0,1] to range [x0,x1].
     *  @param x Float in range [0,1].
     *  @pre x is in range [0,1].
     */
    template <typename XTYPE, typename XTYPE0, typename XTYPE1,
	      typename OFFSET_TYPE>
    double Map01ToRange
    (const XTYPE x, const XTYPE0 x0, const XTYPE1 x1,
     const OFFSET_TYPE offset) const
    {
      const double xdiff = x1 - x0;

      if (xdiff <= 2*offset)
 	  { return ((x0+x1)/2.0); }
      else
 	  { return (x0 + offset + x*(xdiff-2*offset)); }
    }


  public:

    /*!
     *  @brief Return random double in range [x0+offset,x1-offset].
     *  - Random(x0,x1,offset) must be defined in derived classes.
     *  - Distribution is determined by distribution
     *    of GENERATE_RANDOM01_TYPE::Random01.
     */
    template <typename XTYPE0, typename XTYPE1, typename OFFSET_TYPE,
	      typename RANDOM_ENGINE>
    double Random
    (const XTYPE0 x0, const XTYPE1 x1, const OFFSET_TYPE offset,
     RANDOM_ENGINE & random_engine)
    {
      const double x = Random01(random_engine);
      return Map01ToRange(x, x0, x1, offset);
    }


    /*!
     *  Call base class Random01().
     */
    template <typename RANDOM_ENGINE>
    double Random01(RANDOM_ENGINE & random_engine)
    { return GENERATE_RANDOM01_TYPE::Random01(random_engine); }


    /*!
     *  @brief Return random double in range [offset, 1-offset].
     */
    template <typename OFFSET_TYPE, typename RANDOM_ENGINE>
    double Random01
    (const OFFSET_TYPE offset, 
     RANDOM_ENGINE & random_engine)
    { return Random(0, 1, offset, random_engine); }


    /*!
     *  @brief Return random double in range [x0,x1].
     *  - Generate random number in range [x0-boundary_width,x1+boundary_width]
     *    and then clamp to [x0,x1].
     *  - Distribution is determined by distribution
     *    of RANDOM01_TYPE::Random01.
     *  @param boundary_width Width of boundary.
     *    - Values in range [x0-boundary_width, boundary_width] and
     *      [x1, x1+boundary_width] are clamped to the boundary (x0 or x1).
     *    - Note: boundary_width enlarges random number range
     *      to [x0-boundary_width,x1+boundary_width]
     *      in contrast with offset argument to Random() that shrinks
     *      random number range to [x0+offset,x1-offset].
     */
    template <typename XTYPE0, typename XTYPE1, typename OFFSET_TYPE,
              typename RANDOM_ENGINE>
    double RandomAndClamp
    (const XTYPE0 x0, const XTYPE1 x1, const OFFSET_TYPE boundary_width,
     RANDOM_ENGINE & random_engine)
    {
      double x = Random(x0, x1, -boundary_width, random_engine);
      if (x < x0) { x = x0; }
      if (x > x1) { x = x1; }
      return x;
    }

  };


  typedef GENERATE_RANDOM_BASE<GENERATE_UNIFORM_RANDOM01> 
  GENERATE_UNIFORM_RANDOM;

  typedef GENERATE_RANDOM_BASE<GENERATE_U_QUADRATIC_RANDOM01>
  GENERATE_U_QUADRATIC_RANDOM;

  ///@}

}


#endif
