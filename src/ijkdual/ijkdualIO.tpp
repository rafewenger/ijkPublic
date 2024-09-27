/*!
 *  @file ijkdualIO.tpp
 *  @brief IO templates for ijkdual.
 * - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2024 Rephael Wenger

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

#ifndef _IJKDUALIO_TPP_
#define _IJKDUALIO_TPP_

#include <vector>

#include "ijkdual_types.h"

namespace IJKDUAL {

  // ******************************************************************
  //! @name OUTPUT ISOSURFACE
  // ******************************************************************

  //@{

  /*!
   *  @brief Output dual isosurface.
   *  - Use \a output_info.mesh_type to determine type of isosurface.
   */
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    const int THREE(3);
    const int dimension = output_info.Dimension();
    
    if ((dimension == THREE) &&
        (output_info.mesh_type == SIMPLICIAL_COMPLEX)) {
      output_dual_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.simplex_vert, dualiso_info, io_time);
    }
    else if ((dimension == THREE) &&
             (output_info.mesh_type == MIXED_MESH)) {
      output_dual_quad_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dual_isosurface.simplex_vert, 
         dualiso_info, io_time);
    }
    else if ((dimension == THREE) && output_info.flag_dual_collapse) {
      output_dual_quad_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dual_isosurface.simplex_vert,
         dualiso_info, io_time);
    }
    else {
      output_dual_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dualiso_info, io_time);
    }
  }


  /*!
   *  @brief Output dual isosurface and report isosurface information.
   *  - Output file(s) and format(s) depend on flags in \a output_info.
   *  - May write to more than one file.
   *  - Output (report) isosurface information to stdout.
   */
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & slist,
   const DUALISO_INFO_TYPE & dualiso_info, 
   IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info(output_info, dualiso_data, 
                      vertex_coord, slist, dualiso_info);
    }

    if (!output_info.flag_nowrite) 
      { write_dual_mesh(output_info, vertex_coord, slist, io_time); }
  }


  /*!
   *  @brief Output isosurface of  cubes (line segments, quads, hexahedra, ...)
   *  - Isosurface contains only cubes.
   */
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_cube_complex_isosurface
  (const OUTPUT_INFO_TYPE & output_info,
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    output_dual_isosurface
      (output_info, dualiso_data, dual_isosurface.vertex_coord,
       dual_isosurface.isopoly_vert, dualiso_info, io_time);
  }


  /// Output isosurface of triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_tri_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info(output_info, dualiso_data,
                      vertex_coord, tri_vert, dualiso_info);
    }

    if (!output_info.flag_nowrite) {
      write_dual_mesh(output_info, vertex_coord, tri_vert, io_time);
    }
  }


  /// Output isosurface of quadrilaterals and triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_quad_tri_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_quad_tri_iso_info
        (output_info, dualiso_data, vertex_coord, quad_vert, tri_vert, 
         dualiso_info);
    }

    if (!output_info.flag_nowrite) {
      write_dual_quad_tri_mesh
        (output_info, vertex_coord, quad_vert, tri_vert, io_time);
    }
  }


  /// Output isosurface with colored facets.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename COLOR_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_isosurface_color
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & slist,
   const COLOR_TYPE * front_color, 
   const COLOR_TYPE * back_color,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info
        (output_info, dualiso_data, vertex_coord, slist, dualiso_info);
    }
  
    if (!output_info.flag_nowrite) {
      write_dual_mesh_color
        (output_info, vertex_coord, slist, front_color, back_color, io_time);
    }
  }


  /// Output isosurface with colored facets.
  /// - Version with input DUAL_ISOSURFACE_TYPE.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename COLOR_TYPE, typename IO_TIME_TYPE>
  void output_dual_isosurface_color
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    output_dual_isosurface_color
      (output_info, dualiso_data, 
       dual_isosurface.vertex_coord, dual_isosurface.isopoly_vert,
       front_color, back_color, dualiso_info, io_time);
  }

  //@}


  // ******************************************************************
  //! @name REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // ******************************************************************

  //@{

  /// Report information about cubes with multiple vertices.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_INFO_TYPE>  
  void report_isov_info
  (const OUTPUT_INFO_TYPE & output_info, const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int DIM3(3);
    const char * indent4 = "    ";

    using namespace std;

    if (output_info.Dimension() >= DIM3) {
      if (output_info.FlagSeparateIsovNearCubeFacets() ||
          output_info.FlagMoveAllIsovAwayFromCubeFacets()) {
        cout << indent4 << "# times isosurface vertices moved away from cube facets: "
             << dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets
             << endl;
      }
    }

  }

  
  /// Report information about cubes with multiple vertices.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_INFO_TYPE>  
  void report_multi_isov_info
  (const OUTPUT_INFO_TYPE & output_info, const DUALISO_INFO_TYPE & dualiso_info)
  {
    const char * indent4 = "    ";

    using namespace std;

    if (output_info.allow_multiple_iso_vertices) {
      cout << indent4 << "# active (non-empty) cubes: "
           << dualiso_info.scalar.num_non_empty_cubes << endl;
      cout << indent4 << "# cubes with single isosurface vertex: "
           << dualiso_info.multi_isov.num_cubes_single_isov << endl;
      cout << indent4 << "# cubes with multiple isosurface vertices: "
           << dualiso_info.multi_isov.num_cubes_multi_isov << endl;

      if (output_info.flag_split_non_manifold) {
        cout << indent4 << "# cubes changed to 2 iso vertices to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_non_manifold_split << endl;
        cout << indent4 
             << "# ambiguous ridge cubes changed to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_ambig_ridge_cubes_changed << endl;
        cout << indent4 
             << "# non-ambiguous ridge cubes changed to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_non_ambig_ridge_cubes_changed << endl;
      }

      if (output_info.flag_select_split) {
        cout << indent4 << "# cubes changed in selecting isosurface patch splits: "
             << dualiso_info.multi_isov.num_1_2_changed << endl;
      }

      if (output_info.flag_connect_ambiguous) {
        cout << indent4 << "# ambiguous cubes changed to connect isosurface: "
             << dualiso_info.multi_isov.num_connect_changed << endl;
      }
    }
  }


  /// Report triangulation information.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_INFO_TYPE>  
  void report_triangulation_info
  (const OUTPUT_INFO_TYPE & output_info, const DUALISO_INFO_TYPE & dualiso_info)
  {
    const char * indent4 = "    ";
    const char * indent6 = "      ";
    const char * indent8 = "        ";

    using namespace std;

    if (output_info.flag_check_envelope) {
      cout << indent4 << "# quads with some diagonal outside envelope: "
           << dualiso_info.triangulation.num_iso_cubes_with_diag_outside_envelope
           << endl;
    }
    
    if (output_info.flag_trimesh ||
        (output_info.flag_check_envelope &&
         (output_info.flag_allow_tri2_envelope ||
          output_info.flag_allow_tri4_envelope))) {

      cout << indent4 << "# triangulated quads: "
           << dualiso_info.triangulation.num_iso_cubes_tri_total
           << endl;
      
      if (output_info.quad_tri_method != TRI4_ALL_QUADS ||
          (output_info.flag_check_envelope &&
           output_info.flag_allow_tri2_envelope)) {
        
        cout << indent6 << "# quads split into 2 triangles: "
           << dualiso_info.triangulation.num_iso_cubes_tri_no_add
           << endl;

        if (output_info.flag_check_envelope &&
            output_info.flag_allow_tri2_envelope) {

          cout << indent8 << "# quads split into 2 triangles bcuz diagonal outside envelope: "
               << dualiso_info.triangulation.num_tri_no_add_with_diag_outside_envelope
               << endl;
        }
      }

      if (output_info.flag_tri4_quad ||
          (output_info.flag_check_envelope &&
           output_info.flag_allow_tri4_envelope)) {
        
        cout << indent6 << "# quads split into 4 triangles: "
           << dualiso_info.triangulation.num_iso_cubes_tri_add_interior1
           << endl;

        if (output_info.flag_check_envelope &&
            output_info.flag_allow_tri4_envelope) {
          cout << indent8 << "# quads split into 4 triangles bcuz diagonal outside envelope: "
               << dualiso_info.triangulation.num_tri_add_interior1_with_diag_outside_envelope
               << endl;
        }
      }
    }
  }

  /// Report isosurface information.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUALISO_INFO_TYPE>
  void report_iso_info
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & plist, 
   const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int TWO(2);
    const int THREE(3);
    const int FOUR(4);
    const int dimension = output_info.dimension;
    const int numv_per_isopoly = output_info.num_vertices_per_isopoly;
    const char * polytopes_name = "polytopes";
    
    using namespace std;

    VERTEX_INDEX numv = (vertex_coord.size())/dimension;
    VERTEX_INDEX num_poly = (plist.size())/numv_per_isopoly;
    IJK::PROCEDURE_ERROR error("report_iso_info");

    if (output_info.flag_interval_volume) {
      cout << "  Interval volume [" 
           << output_info.isovalue[0] << ":"
           << output_info.isovalue[1] << "].  "
           << numv << " ivol vertices.  "
           << num_poly << " ivol polytopes." << endl;
    }
    else {

      if (output_info.Dimension() == TWO) {
        polytopes_name = "line segments";
      }
      else if (output_info.Dimension() == THREE) {

        // Default in 3D.
        polytopes_name = "polygons";

        if (output_info.mesh_type == SIMPLICIAL_COMPLEX) {
          polytopes_name = "triangles";
        }
        else if (output_info.mesh_type == CUBE_COMPLEX) {
          polytopes_name = "quadrilaterals";
        };
      }
      else if (output_info.Dimension() == FOUR) {
        polytopes_name = "hexahedra";

        // Triangulation not implemented.
      }
      
      cout << "  Isovalue " << output_info.isovalue[0] << ".  "
           << numv << " isosurface vertices.  "
           << num_poly << " isosurface " << polytopes_name
           << "." << endl;
    }

    if (!output_info.flag_use_stdout && !output_info.flag_silent &&
        !output_info.flag_terse) {             
      report_multi_isov_info(output_info, dualiso_info);
      report_isov_info(output_info, dualiso_info);
      report_triangulation_info(output_info, dualiso_info);
    }
  }


  /// Report information about isosurface quadrilaterals and triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUALISO_INFO_TYPE>
  void report_quad_tri_iso_info
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & quad_list, 
   const std::vector<VERTEX_INDEX> & tri_list, 
   const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int dimension = output_info.dimension;
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_TRI(3);
    const char * indent4 = "    ";

    using namespace std;

    const VERTEX_INDEX numv = (vertex_coord.size())/dimension;
    const VERTEX_INDEX num_quad = (quad_list.size())/NUM_VERT_PER_QUAD;
    const VERTEX_INDEX num_tri = (tri_list.size())/NUM_VERT_PER_TRI;

    cout << "  Isovalue " << output_info.isovalue[0] << ".  " 
         << numv << " isosurface vertices.  "
         << num_quad+num_tri << " isosurface polygons." << endl;
    cout << indent4 << num_quad << " quadrilaterals.  "
         << num_tri << " triangles." << endl;

    report_multi_isov_info(output_info, dualiso_info);
    report_triangulation_info(output_info, dualiso_info);
  }

  //@}

}

#endif
