/*!
 *  \file ijkdualIO.cpp
 *  @brief IO routines for ijkdual
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2024 Rephael Wenger

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


#include <assert.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


#include "ijkcommand_line.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkisoIO.tpp"
#include "ijkIO.tpp"
#include "ijkmesh.tpp"
#include "ijkstring.tpp"
#include "ijkprint.tpp"

#include "ijkdualIO.h"

using namespace IJK;
using namespace IJKDUAL;

using namespace std;

// ******************************************************************
// PARSE COMMAND LINE
// ******************************************************************

// local namespace
namespace {

  typedef typename IJK::COMMAND_LINE_OPTIONS
    <COMMAND_LINE_OPTION_TYPE,COMMAND_LINE_OPTION_GROUP>
  COMMAND_LINE_OPTIONS;

  COMMAND_LINE_OPTIONS options;
};

// local namespace
namespace {

  /// @brief Convert string s into interpolation type.
  INTERPOLATION_TYPE get_interpolation_type(char * s)
  {
    INTERPOLATION_TYPE type = LINEAR_INTERPOLATION;

    if (strcmp(s, "linear") == 0) 
      { type = LINEAR_INTERPOLATION; }
    else if (strcmp(s, "multilinear") == 0) 
      { type = MULTILINEAR_INTERPOLATION; }
    else {
      cerr << "Error in input parameter -interpolate.  Illegal interpolation type: " 
           << s << "." << endl;
      exit(64);
    }

    return(type);
  }

  // Convert string s into vertex position method.
  VERTEX_POSITION_METHOD get_vertex_position_method(char * s)
  {
    VERTEX_POSITION_METHOD method = CENTROID_EDGE_ISO;
    string str = s;

    if (str == "cube_center") {
      method = CUBE_CENTER;
    }
    else if (str == "centroid") {
      method = CENTROID_EDGE_ISO;
    }
    else if (str == "ivol_lifted02") {
      method = IVOL_LIFTED02;
    }
    else if (str == "random") {
      method = RANDOM_ISOV_POS;
    }
    else {
      cerr << "Error in input parameter -position.  Illegal position method: " 
           << str << "." << endl;
      exit(64);
    }

    return(method);
  }


  /// @brief Convert string s into random distribution.
  RANDOM_DISTRIBUTION get_random_distribution(char * s)
  {
    RANDOM_DISTRIBUTION distribution = UNDEFINED_DISTRIBUTION;
    string str = s;

    if (str == "uniform") {
      distribution = UNIFORM_DISTRIBUTION;
    }
    else if (str == "Uquadratic" || str == "quadratic") {
      distribution = U_QUADRATIC_DISTRIBUTION;
    }
    else {
      cerr << "Error in input parameter -rpos_distribution.  Illegal distribution: "
           << str << "." << endl;
      exit(64);
    }

    return(distribution);
  }

}

namespace {

  void set_command_line_options(const IO_INFO & io_info)
  {
    int iarg;

    options.SetLabelWidth(15);
    options.AddLabelTab(20);

    // Regular options.
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (TRIMESH_OPT, "TRIMESH_OPT", REGULAR_OPTG, "-trimesh", 
       "Triangulate all isosurface quadrilaterals.");
    options.AddToHelpMessage
      (TRIMESH_OPT,
       "Default triangulation method: Maximize minimum angle.",
       "Option only works with 3D scalar data.");

    options.AddOptionNoArg
      (TRIMESH_MAX_MIN_ANGLE_OPT, "TRIMESH_MAX_MIN_ANGLE_OPT", REGULAR_OPTG,
       "-trimesh_max_min_angle",
       "Triangulate all isosurface quadrilaterals.");
    options.AddToHelpMessage
    (TRIMESH_MAX_MIN_ANGLE_OPT,
     "Split quadrilateral into two triangles, choosing split",
     "that maximizes the minimum triangle angle.",
     "Option only works with 3D scalar data.");

    options.AddOptionNoArg
      (TRIMESH_UNIFORM_OPT, "TRIMESH_UNIFORM_OPT", REGULAR_OPTG,
       "-trimesh_uniform",
       "Triangulate all isosurface quadrilaterals.");
    options.AddToHelpMessage
      (TRIMESH_UNIFORM_OPT, "Uniformly split quadrilaterals into two triangle.",
       "Option only works with 3D scalar data.");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOptionNoArg
       (TRIMESH_TRI4_MAX_MIN_ANGLE_OPT, "TRIMESH_TRI4_MAX_MIN_ANGLE_OPT", REGULAR_OPTG,
        "-trimesh_tri4_max_min_angle",
        "Triangulate all isosurface quadrilaterals.");
    options.AddToHelpMessage
      (TRIMESH_TRI4_MAX_MIN_ANGLE_OPT,
       "Split quadrilaterals into two or four triangles.",
       "Add interior vertex to triangulate into four triangles.",
       "Choose triangulation which maximizes the minimum angle.",
       "Option only works with 3D scalar data.");
    options.AddSynonym(TRIMESH_TRI4_MAX_MIN_ANGLE_OPT, "-tri4");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (ENVELOPE_TRI2_OR_TRI4_OPT, "ENVELOPE_TRI2_OR_TRI4_OPT",
       REGULAR_OPTG, "-envelope",
       "Check envelope when triangulating.");
    options.AddToHelpMessage
      (ENVELOPE_TRI2_OR_TRI4_OPT,
       "Envelope is union of four tetrahedron where each",
       "tetrahedron is the convex hull of each quad edge",
       "and the grid edge dual to the quad.");
    options.AddToHelpMessage
      (ENVELOPE_TRI2_OR_TRI4_OPT,
       "See \"Intersection free contouring on an octree grid\"",
       "by Ju and Udeshi, Proceedings of Pacific Graphics, 2006.",
       "Triangulate if one or both diagonals are outside the envelope.");
    options.AddToHelpMessage
      (ENVELOPE_TRI2_OR_TRI4_OPT,
       "If only one diagonal is outside the envelope,",
       "triangulate using other diagonal or into four triangles,",
       "choosing triangulation that maximizes minimum angle.",
       "If both diagonals are outside_the envelope,",
       "triangulate into four triangles.");
    options.AddSynonym(ENVELOPE_TRI2_OR_TRI4_OPT, "-envelope_tri2_or_tri4");

    options.AddOptionNoArg
      (ENVELOPE_ONLY_TRI2_OPT, "ENVELOPE_ONLY_TRI2_OPT",
       REGULAR_OPTG, "-envelope_only_tri2",
       "Check envelope when triangulating.");
    options.AddToHelpMessage
      (ENVELOPE_ONLY_TRI2_OPT,
       "Triangulate only into two triangles based on envelope.",
       "If one diagonal is outside the envelope and",
       "one is inside the envelope, triangulate with the diagonal",
       "inside the envelope.",
       "If both diagonals are outside the envelope, ignore the envelope.",
       "Do not triangulate into four triangles based on envelope.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOption1Arg
      (ENVELOPE_OFFSET_OPT, "ENVELOPE_OFFSET_OPT",
       REGULAR_OPTG, "-envelope_offset", "{offset}",
       "Set envelope offset.");
    options.AddToHelpMessageDefaultValue
      (ENVELOPE_OFFSET_OPT, "-envelope_offset", io_info.DefaultEnvelopeOffset());

    options.AddOptionNoArg
      (ENVELOPE_USE_RELATIVE_POSITIONS_OPT,
       "ENVELOPE_USE_RELATIVE_POSITIONS_OPT", REGULAR_OPTG,
       "-envelope_rel_pos",
       "Use relative positions to help determine diagonals in envelope.");

    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (TRI4_IN_ENVELOPE_OPT, "TRI4_IN_ENVELOPE_OPT", REGULAR_OPTG,
       "-tri4_in_envelope",
       "Position additional triangulation vertex in envelope",
       "formed by isosurface quadrilateral edges and dual grid edges.");
    options.AddToHelpMessage(TRI4_IN_ENVELOPE_OPT, "(Default.)");

    options.AddOptionNoArg
      (TRI4_ON_GRID_EDGE_OPT, "TRI4_ON_GRID_EDGE_OPT", REGULAR_OPTG,
       "-tri4_on_grid_edge",
       "Position additional triangulation vertex on grid edge.");
       
    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);
    
    options.AddOptionNoArg
     (SEPARATE_ISOV_OPT, "SEPARATE_ISOV_OPT", REGULAR_OPTG,
      "-separate_isov",
      "Separate isosurface vertices to avoid mesh self intersections.");
    options.AddToHelpMessage
      (SEPARATE_ISOV_OPT,
       "Separate isosurface vertices from cube faces by offset.",
       "Separate multiple isosurface vertices within a cube",
       "by cube centers.",
       "Offset is applied both to separation from cube faces and",
       "from cube center coordinates to avoid self intersections",
       "caused by collinearities and coplanarities.");
    options.AddToHelpMessageDefaultValue
      (SEPARATE_ISOV_OPT, "offset", io_info.DefaultPositionOffset());
    options.AddToHelpMessage
      (SEPARATE_ISOV_OPT,
       "Equivalent to \"-sep_isov_near_cube_facets -separate_by_cc\".",
       "(Use -help_more for explanations of those options.)");

    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (SEP_ISOV_NEAR_CUBE_FACETS_OPT, "SEP_ISOV_NEAR_CUBE_FACETS_OPT",
       REGULAR_OPTG, "-sep_isov_near_cube_facets",
       "Separate isosurface vertices that are near shared cube facets.");
    options.AddToHelpMessage
      (SEP_ISOV_NEAR_CUBE_FACETS_OPT,
       "Separate isosurface vertices if two vertices from different",
       "grid cubes are both near the same cube facet.",
       "Move vertices that are within {min_distance} of cube facet.",
       "Move distance {position_offset} away from cube facet.");
    options.AddToHelpMessage
      (SEP_ISOV_NEAR_CUBE_FACETS_OPT,
       "Also moves all isosurface vertices away from cube ridges.",
       "Use option \"-min_distance_to_cube_face\" to set {min_distance}.",
       "Use option \"-position_offset\" to set {position_offset}.");

    options.AddOptionNoArg
    (MOVE_ALL_ISOV_AWAY_FROM_CUBE_FACETS_OPT,
     "MOVE_ALL_ISOV_AWAY_FROM_CUBE_FACETS_OPT",
     REGULAR_OPTG, "-move_isov_away_from_facets",
     "Move all isosurface vertices away from cube facets.");
    options.AddToHelpMessage
      (MOVE_ALL_ISOV_AWAY_FROM_CUBE_FACETS_OPT,
     "Move distance {position_offset} away from cube facets.",
     "Use option \"-position_offset\" to set {position_offset}.",
     "Supercedes -sep_isov_near_cube_facets if both options",
     "are set.");
    options.AddSynonym
      (MOVE_ALL_ISOV_AWAY_FROM_CUBE_FACETS_OPT, "-move_to_interior");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    
    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOption1Arg
      (SUBSAMPLE_OPT, "SUBSAMPLE_OPT", REGULAR_OPTG,
       "-subsample", "S", "Subsample grid at every S vertices.");
    options.AddToHelpMessage
      (SUBSAMPLE_OPT, "S must be an integer greater than 1.");

    options.AddOption1Arg
      (SUPERSAMPLE_OPT, "SUPERSAMPLE_OPT", REGULAR_OPTG,
       "-supersample", "S",
       "Supersample grid at every S vertices.");
    options.AddToHelpMessage
      (SUPERSAMPLE_OPT, "S must be an integer greater than 1.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOption1Arg
      (POSITION_OPT, "POSITION_OPT", REGULAR_OPTG,
       "-position", "{centroid|cube_center}",
       "Isosurface vertex position method.");
    iarg = options.AddArgChoice
      (POSITION_OPT, "centroid",
       "Position isosurface vertices at centroid of");
    options.AddToHelpArgMessage
      (POSITION_OPT, iarg,
       "intersection of grid edges and isosurface.");
    iarg = options.AddArgChoice
      (POSITION_OPT, "cube_center",
       "Position isosurface vertices at cube_centers.");

    options.AddOptionNoArg
      (CUBE_CENTER_OPT, "CUBE_CENTER_OPT", REGULAR_OPTG, "-cube_center",
       "Position isosurface vertices at cube_centers.");
    options.AddToHelpMessage
      (CUBE_CENTER_OPT, "Equivalent to \"-position cube_center\".");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddUsageOptionBeginOr(REGULAR_OPTG);
    options.AddOptionNoArg
      (MANIFOLD_OPT, "MANIFOLD_OPT", REGULAR_OPTG, "-manifold",
       "Isosurface mesh has a manifold representation. (Default.)");
    options.AddToHelpMessage
      (MANIFOLD_OPT, "Allows multiple isosurface vertices per cube and",
       "splits isosurface vertices to avoid non-manifold edges.");

    options.AddOptionNoArg
      (MULTI_ISOV_OPT, "MULTI_ISOV_OPT", REGULAR_OPTG, "-multi_isov",
       "Allow multiple isosurface vertices per cube but");
    options.AddToHelpMessage
      (MULTI_ISOV_OPT,
       "don't change isosurface mesh to avoid non-manifold edges.",
       "Isosurface may have non-manifold edges/vertices.");

    options.AddOptionNoArg
      (SINGLE_ISOV_OPT, "SINGLE_ISOV_OPT", REGULAR_OPTG, "-single_isov",
       "Single isosurface vertex per cube.");
    options.AddToHelpMessage
      (SINGLE_ISOV_OPT, "Isosurface may have non-manifold edges/vertices.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOptionNoArg
      (OFF_OPT, "OFF_OPT", REGULAR_OPTG, "-off", 
       "Output in geomview OFF format. (Default.)");

    options.AddOptionNoArg
      (PLY_OPT, "PLY_OPT", REGULAR_OPTG, "-ply", 
       "Output in Stanford Polygon File Format (PLY).");
    options.AddToHelpMessage
      (PLY_OPT, "(Allowed only with 3D scalar data.)");

    options.AddOptionNoArg
      (FIG_OPT, "FIG_OPT", REGULAR_OPTG, "-fig",
       "Output in FIG (xfig) format.");
    options.AddToHelpMessage
      (FIG_OPT, "(Allowed only with 2D scalar data.)");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddUsageOptionBeginOr(REGULAR_OPTG);
    options.AddOption1Arg
      (OUTPUT_FILENAME_OPT, "OUTPUT_FILENAME_OPT", REGULAR_OPTG, 
       "-o", "{output_filename}", 
       "Write isosurface to file {output_filename}.");

    options.AddOption1Arg
      (OUTPUT_FILENAME_PREFIX_OPT, "OUTPUT_FILENAME_PREFIX_OPT", REGULAR_OPTG, 
       "-prefix", "{prefix}", 
       "Use {prefix} as filename prefix in constructing",
       "output file name.");

    options.AddOptionNoArg
      (STDOUT_OPT, "STDOUT_OPT", REGULAR_OPTG, 
       "-stdout", "Write isosurface to standard output.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (LABEL_WITH_ISOVALUE_OPT, "LABEL_WITH_ISOVALUE_OPT", REGULAR_OPTG,
       "-label_with_isovalue", "Include isovalue in output file name.");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddUsageOptionBeginOr(REGULAR_OPTG);
    
    options.AddOptionNoArg
      (SILENT_OPT, "SILENT_OPT", REGULAR_OPTG, "-silent", 
       "Silent mode.  No output messages.");
    options.AddSynonym(SILENT_OPT, "-s");

    options.AddOptionNoArg
      (TERSE_OPT, "TERSE_OPT", REGULAR_OPTG, "-terse",
       "Terse mode. Terse output message.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    
    options.AddOptionNoArg
      (NO_WARN_OPT, "NO_WARN_OPT", REGULAR_OPTG, "-no_warn", 
       "Do not output non-manifold warning.");

    options.AddOptionNoArg
      (NO_WRITE_OPT, "NO_WRITE_OPT", REGULAR_OPTG, "-no_write", 
       "Don't wite isosurface.");

    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (TIME_OPT, "TIME_OPT", REGULAR_OPTG, "-time", "Output running time.");

    options.AddOptionNoArg
      (DETAILED_TIME_OPT, "DETAILED_TIME_OPT", REGULAR_OPTG, "-detailed_time",
       "Output detailed running time.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);

    /* NOT CURRENTLY USED
    options.AddOptionNoArg
      (INFO_OPT, "INFO_OPT", REGULAR_OPTG, "-info", 
       "Output more information about isosurface.");
    */

    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (USAGE_OPT, "USAGE_OPT", REGULAR_OPTG, "-usage", 
       "Print usage message.");

    options.AddOptionNoArg
       (MORE_OPTIONS_OPT, "MORE_OPTIONS_OPT", REGULAR_OPTG, "-more_options",
        "Print additional (more) options.");

    options.AddOptionNoArg
       (TESTING_OPTIONS_OPT, "TESTING_OPTIONS_OPT", REGULAR_OPTG, "-testing_options",
        "Print options used for testing.");

    options.AddOptionNoArg
      (ALL_OPTIONS_OPT, "ALL_OPTIONS_OPT", REGULAR_OPTG, "-all_options", 
       "Print usage message with all options.");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOptionNoArg
      (HELP_OPT, "HELP_OPT", REGULAR_OPTG, "-help",
       "Print help message for standard options.");
    options.AddSynonym(HELP_OPT, "-h");

    options.AddOptionNoArg
       (HELP_MORE_OPT, "HELP_MORE_OPT", REGULAR_OPTG, "-help_more",
        "Print help message for additional (more) options.");

    options.AddOptionNoArg
        (HELP_TESTING_OPT, "HELP_TESTING_OPT", REGULAR_OPTG, "-help_testing",
         "Print help message for testing options.");

    options.AddOptionNoArg
      (HELP_ALL_OPT, "HELP_ALL_OPT", REGULAR_OPTG, "-help_all", 
       "Print help message for all options.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);


    // Extended options.

    options.AddUsageOptionBeginOr(EXTENDED_OPTG);

    options.AddOptionNoArg
      (SEP_NEG_OPT, "SEP_NEG_OPT", EXTENDED_OPTG, "-sep_neg", 
       "Isosurface patches separate negative grid vertices.");

    options.AddOptionNoArg
      (SEP_POS_OPT, "SEP_POS_OPT", EXTENDED_OPTG, "-sep_pos", 
       "Isosurface patches separate positive grid vertices.");

    options.AddUsageOptionEndOr(EXTENDED_OPTG);

    options.AddOptionNoArg
      (SELECT_SPLIT_OPT, "SELECT_SPLIT_OPT", EXTENDED_OPTG, "-select_split", 
       "Select which cube has configuration of split isosurface");
    options.AddToHelpMessage
      (SELECT_SPLIT_OPT, 
       "vertices where adjacent cubes share an ambiguous facet.");

    options.AddOptionNoArg
      (CONNECT_AMBIG_OPT, "CONNECT_AMBIG_OPT", EXTENDED_OPTG, 
       "-connect_ambig", "Select connection of ambiguous cubes to increase");
    options.AddToHelpMessage
      (CONNECT_AMBIG_OPT, "connectivity of the isosurface."); 

    options.AddUsageOptionNewline(EXTENDED_OPTG);
    options.AddUsageOptionBeginOr(EXTENDED_OPTG);

    options.AddOptionNoArg
      (TRIMESH_SPLIT_MAX_OPT, "TRIMESH_SPLIT_MAX_OPT", EXTENDED_OPTG,
       "-trimesh_split_max", 
       "Triangulate all isosurface quadrilaterals.");
    options.AddToHelpMessage
      (TRIMESH_SPLIT_MAX_OPT, "Split max quadrilateral angle.",
       "Option only works with 3D scalar data.");

    options.AddOptionNoArg
      (TRIMESH_TRI4_ALL_QUADS_OPT, "TRIMESH_TRI4_ALL_QUADS_OPT", EXTENDED_OPTG, 
       "-trimesh_tri4_all",
       "Triangulate all isosurface quadrilaterals");
    options.AddToHelpMessage
      (TRIMESH_TRI4_ALL_QUADS_OPT, "into four triangles.",
       "Option only works with 3D scalar data.");

    options.AddUsageOptionEndOr(EXTENDED_OPTG);
    options.AddUsageOptionNewline(EXTENDED_OPTG);

    options.AddOptionNoArg
      (ENVELOPE_TRI2_OR_TRI4_PREFER_TRI2_OPT, "ENVELOPE_TRI2_OR_TRI4_PREFER_TRI2_OPT",
       EXTENDED_OPTG, "-envelope_tri2_or_tri4_prefer_tri2",
       "Check envelope when triangulating.");
    options.AddToHelpMessage
      (ENVELOPE_TRI2_OR_TRI4_PREFER_TRI2_OPT,
       "Triangulate if one or both diagonals are outside the envelope.",
       "If only one diagonal is outside the envelope,",
       "triangulate using other diagonal.",
       "If both diagonals are outside_the envelope,",
       "triangulate into four triangles.");

    options.AddOptionNoArg
      (ENVELOPE_ONLY_TRI4_OPT, "ENVELOPE_ONLY_TRI4_OPT",
       EXTENDED_OPTG, "-envelope_only_tri4",
       "Triangulate if some diagonal is outside envelope.");
    options.AddToHelpMessage
      (ENVELOPE_ONLY_TRI4_OPT,
       "Triangulate into four triangles, even if one diagonal",
       "is inside the envelope and one diagonal is outside.");

    options.AddUsageOptionNewline(EXTENDED_OPTG);

    options.AddOptionNoArg
    (SEPARATE_ISOV_BY_CUBE_CENTER_OPT,
     "SEPARATE_ISOV_BY_CUBE_CENTER_OPT", EXTENDED_OPTG,
     "-separate_isov_by_cube_center",
     "Separate multiple isosurface vertices in a cube.");
    options.AddToHelpMessage
    (SEPARATE_ISOV_BY_CUBE_CENTER_OPT,
     "Separate isosurface vertices by the cube center.",
     "Apply position offset to further separate vertices",
     "if position offset flag is true.");
    options.AddSynonym
      (SEPARATE_ISOV_BY_CUBE_CENTER_OPT, "-separate_by_cc");

    options.AddUsageOptionNewline(EXTENDED_OPTG);
    
    options.AddOption1Arg
      (MIN_DISTANCE_TO_CUBE_FACE_OPT,
       "MIN_DISTANCE_TO_CUBE_FACE_OPT", EXTENDED_OPTG,
       "-min_distance_to_cube_face", "{x}",
       "Set position offset and apply to isosurface vertices.");
    options.AddToHelpMessage
      (MIN_DISTANCE_TO_CUBE_FACE_OPT,
       "Apply offset to isosurface vertex positions",
       "when vertices are within (Linf) distance {x} of cube faces.",
       "Necessary to avoid surface self intersections",
       "in some degenerate cases.");
    options.AddSynonym
      (MIN_DISTANCE_TO_CUBE_FACE_OPT, "-min_distance2face");
    options.AddToHelpMessageDefaultValue
      (MIN_DISTANCE_TO_CUBE_FACE_OPT, "minium distance",
       io_info.DefaultMinDistanceToCubeFace());
    
    options.AddUsageOptionNewline(EXTENDED_OPTG);
    options.AddUsageOptionBeginOr(EXTENDED_OPTG);

    options.AddOption1Arg
      (POSITION_OFFSET_OPT, "POSITION_OFFSET_OPT", EXTENDED_OPTG,
       "-position_offset", "{x}",
       "Set position offset and apply to isosurface vertices.");
    options.AddToHelpMessage
      (POSITION_OFFSET_OPT, "Apply offset {x} to isosurface vertex positions",
       "to move isosurface vertices away from cube facets",
       "and/or away from each other.",
       "Necessary to avoid surface self intersections",
       "in some degenerate cases.");
    options.AddToHelpMessageDefaultValue
      (POSITION_OFFSET_OPT, "offset", io_info.DefaultPositionOffset());

    options.AddUsageOptionNewline(EXTENDED_OPTG);

    options.AddOptionNoArg
       (NO_POSITION_OFFSET_OPT, "NO_OFFSET_OPT", EXTENDED_OPTG,
        "-no_position_offset", "Do not apply offset to isosurface vertex positions.");
    options.AddToHelpMessage
      (NO_POSITION_OFFSET_OPT,
       "Do not reposition isosurface vertices to separate",
       "from facets or from other isosurface vertices.");
    options.AddSynonym(NO_POSITION_OFFSET_OPT, "-no_offset");

    options.AddUsageOptionEndOr(EXTENDED_OPTG);
    options.AddUsageOptionNewline(EXTENDED_OPTG);

    // Quadrilateral edge intersection (qei) options.
    options.AddUsageOptionBeginOr(EXTENDED_OPTG);

    options.AddOptionNoArg
      (QEI_INTERPOLATE_SCALAR_OPT, "QEI_INTERPOLATE_SCALAR_OPT", 
       EXTENDED_OPTG, "-qei_interpolate_scalar", 
       "Compute quadrilateral edge intersection using scalar interpolation.");
    options.AddToHelpMessage
      (QEI_INTERPOLATE_SCALAR_OPT,
       "Computer intersection of isosurface quadrilaterals and",
       "grid edges using linear interpolation of scalar values",
       "at grid edge endpoints.");

    options.AddOptionNoArg
      (QEI_BILINEAR_OPT, "QEI_BILINEAR_OPT", 
       EXTENDED_OPTG, "-qei_bilinear", 
       "Compute quadrilateral edge intersection using bilinear interpolation.");
    options.AddToHelpMessage
      (QEI_BILINEAR_OPT,
       "Compute intersection of isosurface quadrilaterals and",
       "grid edges by intersecting the bilinear surface patch",
       "defined by the quadrilateral and the grid edge.");

    options.AddOptionNoArg
      (QEI_AVERAGE_OPT, "QEI_AVERAGE_OPT", EXTENDED_OPTG,
       "-qei_average", 
       "Compute quadrilateral edge intersection using projection average.");
    options.AddToHelpMessage
      (QEI_AVERAGE_OPT,
       "Compute intersection of isosurface quadrilaterals and",
       "grid edges using the average of the projection",
       "of quadrilateral vertices on the grid edge.");

    options.AddUsageOptionEndOr(EXTENDED_OPTG);
    options.AddUsageOptionNewline(EXTENDED_OPTG);

    options.AddOption1Arg
      (INTERIOR_VERTEX_RATIO_OPT, "INTERIOR_VERTEX_RATIO_OPT",
       EXTENDED_OPTG, "-interior_vertex_ratio", "{R}",
       "Set ratio of distance of interior vertex to cube facet",
       "over distance of isosurface polytope vertex to cube facet.");
    options.AddToHelpMessage
      (INTERIOR_VERTEX_RATIO_OPT,
       "Used in positioning additional vertex to split a quadrilateral",
       "into four triangles.",
       "Must be non-negative and less than 0.5.");
    options.AddToHelpMessageDefaultValue
      (INTERIOR_VERTEX_RATIO_OPT, "-interior_vertex_ratio",
       io_info.DefaultInteriorVertexRatioToCubeFacet());

    options.AddUsageOptionNewline(EXTENDED_OPTG);
    
    options.AddOption1Arg
      (OUT_ISOV_OPT, "OUT_ISOV_OPT", EXTENDED_OPTG, 
       "-out_isov", "{output_filename}", 
       "Write information about isosurface vertices");
    options.AddToHelpMessage
      (OUT_ISOV_OPT, "to file {output_filename}.");

    options.AddOptionNoArg
      (NO_COMMENTS_OPT, "NO_COMMMENTS_OPT", EXTENDED_OPTG,
       "-no_comments",
       "Do not put any comments in isosurface output file.");
    options.AddToHelpMessage
      (NO_COMMENTS_OPT,
       "(Mainly used for testing/comparing outputs",
       "without comparing differences in comments.)");


    // Testing options.

    options.AddOption1Arg
      (POSITION_TESTING_OPT, "POSITION_TESTING_OPT", TESTING_OPTG,
       "-position", "{random}",
       "Isosurface vertex position method for testing.");
    iarg = options.AddArgChoice
      (POSITION_TESTING_OPT, "random", 
       "Randomly position isosurface vertex in cube.");
    options.AddToHelpArgMessage
      (POSITION_TESTING_OPT, iarg,
       "Requires option \"-random_seed {S}\".");

    options.AddOption1Arg
      (RANDOM_SEED_OPT, "RANDOM_SEED_OPT", TESTING_OPTG, 
       "-random_seed", "{S}", 
       "Set random generator seed to {S}.");

    options.AddUsageOptionNewline(TESTING_OPTG);
    options.AddUsageOptionBeginOr(TESTING_OPTG);

    options.AddOption1Arg
      (RANDOM_POS_OPT, "RANDOM_POS_OPT", TESTING_OPTG, 
       "-random_pos", "{S}", 
       "Randomly position isosurface vertices.");
    options.AddToHelpMessage
      (RANDOM_POS_OPT, "S is the random generator seed.",
       "Use non-uniform distribution with more coordinates",
       "on boundary of permitted region.");

    options.AddOption1Arg
      (RANDOM_POS_UNIFORM_OPT, "RANDOM_POS_UNIFORM_OPT", TESTING_OPTG, 
       "-random_posU", "{S}",
       "Randomly position isosurface vertices.");
    options.AddToHelpMessage
      (RANDOM_POS_UNIFORM_OPT, 
       "Uniformly distribute vertices in permitted region.",
       "S is the random generator seed.");

    options.AddUsageOptionNewline(TESTING_OPTG);

    options.AddOption1Arg
      (RANDOM_POS_DISTRIBUTION_OPT, "RANDOM_POS_DISTRIBUTION_OPT", TESTING_OPTG,
       "-rpos_distribution", "{uniform|Uquadratic}",
       "Random position distribution.");
    iarg = options.AddArgChoice
      (RANDOM_POS_DISTRIBUTION_OPT, "uniform",
       "Randomly position using uniform distribution.");
    iarg = options.AddArgChoice
      (RANDOM_POS_DISTRIBUTION_OPT, "Uquadratic",
       "Randomly position using U-quadratic distribution.");
    options.AddToHelpArgMessage
      (RANDOM_POS_DISTRIBUTION_OPT, iarg,
       "U-quadratic distribution positions more vertices",
       "near boundaries.");

    options.AddUsageOptionEndOr(TESTING_OPTG);
    options.AddUsageOptionNewline(TESTING_OPTG);
    options.AddUsageOptionBeginOr(TESTING_OPTG);

    options.AddOptionNoArg
      (RANDOM_SAMPLE_ENTIRE_CUBE_OPT, "RANDOM_SAMPLE_ENTIRE_CUBE",
       TESTING_OPTG, "-rpos_sample_entire_cube",
       "Generate random points sampling entire cube.");
    options.AddToHelpMessage
      (RANDOM_SAMPLE_ENTIRE_CUBE_OPT,
       "Useful for testing options -apply_position_offset",
       "and -separate_isov_by_cube_centers.");

    options.AddOptionNoArg
      (RANDOM_SAMPLE_CUBE_BOUNDARY_OPT, "RANDOM_SAMPLE_CUBE_BOUNDARY",
       TESTING_OPTG, "-rpos_sample_boundary",
       "Generate random points sampling entire cube and cube boundary.");
    options.AddToHelpMessage
       (RANDOM_SAMPLE_CUBE_BOUNDARY_OPT,
        "Useful for testing options -apply_position_offset",
        "and -separate_isov_by_cube_centers.");

    options.AddUsageOptionEndOr(TESTING_OPTG);
    options.AddUsageOptionNewline(TESTING_OPTG);

    options.AddOption1Arg
      (RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT,
       "RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT", TESTING_OPTG,
       "-rpos_boundary_width", "{w}",
       "Set width to be used for random sampling of boundary.");
    options.AddToHelpMessage
      (RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT,
       "Width of size 0.05 puts approximately 15%",
       "of the isosurface vertices on 3D cube boundaries.");
    options.AddToHelpMessage
      (RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT,
       "Width of size 0.1 puts approximately 25%",
       "of the isosurface vertices on 3D cube boundaries.");
    options.AddToHelpMessage
      (RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT,
       "Width of size w puts (1 - (1/(1+w))^3",
       "of the isosurface vertices on 3D cube boundaries.");
    options.AddToHelpMessageDefaultValue
      (RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT, "boundary width",
       io_info.DefaultRandomPosBoundaryWidth());

    options.AddUsageOptionNewline(TESTING_OPTG);

    options.AddOptionNoArg
      (TRI4_CENTROID_OPT, "TRI4_CENTROID_OPT", TESTING_OPTG, "-tri4_centroid",
       "Position additional triangulation vertex at centroid",
       "of isosurface quadrilateral vertices.");

    options.AddUsageOptionNewline(TESTING_OPTG);

    options.AddOptionNoArg
      (SPLIT_NON_MANIFOLD_AMBIG_OPT, "SPLIT_NON_MANIFOLD_AMBIG_OPT",
       TESTING_OPTG, "-split_non_manifold_ambig",
       "Split non-manifold isosurface vertices",
       "caused by ambiguous facets.");
    options.AddToHelpMessage
      (SPLIT_NON_MANIFOLD_AMBIG_OPT,
       "Used with -multi_isov option to test routines",
       "that prevent non-manifold conditions.",
       "To ensure manifold, use -manifold option, not -multi_isov.",
       "The -manifold option implicitly includes this option.");

    options.AddOptionNoArg
      (FORCE_ADD_ALL_INTERIOR_VERTICES_OPT,
       "FORCE_ADD_ALL_INTERIOR_VERTICES_OPT",
       TESTING_OPTG, "-force_add_all_interior_vertices",
       "Always add all interior vertices before any triangulation.");
  }

};


/*!
 *  @brief Store option information in io_info.
 *  - Return true if optA is processed, i.e. optA matches one
 *    of the cases in the switch statement.
 */
template <typename IO_INFO_TYPE>
bool process_option
(const COMMAND_LINE_OPTION_TYPE optA, const int argc, char **argv,
 int & iarg, IO_INFO_TYPE & io_info)
{
  VERTEX_POSITION_METHOD vpos_method;
  RANDOM_SEED_TYPE random_seed;
  RANDOM_DISTRIBUTION random_distribution;
  IJK::ERROR error;

  switch(optA) {

  case SUBSAMPLE_OPT:
    io_info.subsample_resolution = get_arg_int(iarg, argc, argv, error);
    io_info.flag_subsample = true;
    iarg++;
    break;

  case SUPERSAMPLE_OPT:
    io_info.supersample_resolution = get_arg_int(iarg, argc, argv, error);
    io_info.flag_supersample = true;
    iarg++;
    break;

  case POSITION_OPT:
    iarg++;
    if (iarg >= argc) usage_error();
    vpos_method = get_vertex_position_method(argv[iarg]);
    io_info.SetVertexPositionMethod(vpos_method);
    break;

  case CUBE_CENTER_OPT:
    io_info.SetVertexPositionMethod(CUBE_CENTER);
    break;

  case MANIFOLD_OPT:
    io_info.allow_multiple_iso_vertices = true;
    io_info.flag_split_non_manifold = true;
    break;

  case SPLIT_NON_MANIFOLD_AMBIG_OPT:
    // Used with -multi_isov for testing.
    // No need to set this to true if flag_split_non_manifold is true;
    io_info.flag_split_non_manifold_ambig = true;
    break;

  case MULTI_ISOV_OPT:
    io_info.allow_multiple_iso_vertices = true;
    io_info.flag_split_non_manifold = false;
    break;

  case SINGLE_ISOV_OPT:
    io_info.allow_multiple_iso_vertices = false;
    break;

  case SELECT_SPLIT_OPT:
    io_info.allow_multiple_iso_vertices = true;
    io_info.flag_select_split = true;
    break;

  case CONNECT_AMBIG_OPT:
    io_info.flag_connect_ambiguous = true;
    io_info.allow_multiple_iso_vertices = true;
    break;

  case SEP_NEG_OPT:
    io_info.allow_multiple_iso_vertices = true;
    io_info.flag_separate_neg = true;
    break;

  case SEP_POS_OPT:
    io_info.allow_multiple_iso_vertices = true;
    io_info.flag_separate_neg = false;
    break;

  case TRIMESH_OPT:
    io_info.mesh_type = SIMPLICIAL_COMPLEX;
    io_info.flag_trimesh = true;
    break;

  case TRIMESH_SPLIT_MAX_OPT:
    io_info.mesh_type = SIMPLICIAL_COMPLEX;
    io_info.quad_tri_method = SPLIT_MAX_ANGLE;
    io_info.flag_trimesh = true;
    break;

  case TRIMESH_MAX_MIN_ANGLE_OPT:
    io_info.mesh_type = SIMPLICIAL_COMPLEX;
    io_info.quad_tri_method = MAX_MIN_ANGLE;
    io_info.flag_trimesh = true;
    break;

  case TRIMESH_TRI4_ALL_QUADS_OPT:
    io_info.mesh_type = SIMPLICIAL_COMPLEX;
    io_info.quad_tri_method = TRI4_ALL_QUADS;
    io_info.flag_tri4_quad = true;
    io_info.flag_trimesh = true;
    break;

  case TRIMESH_TRI4_MAX_MIN_ANGLE_OPT:
    io_info.mesh_type = SIMPLICIAL_COMPLEX;
    io_info.quad_tri_method = TRI4_MAX_MIN_ANGLE;
    io_info.flag_tri4_quad = true;
    io_info.flag_trimesh = true;
    break;

  case ENVELOPE_TRI2_OR_TRI4_OPT:
    io_info.flag_check_envelope = true;
    io_info.flag_allow_tri2_envelope = true;
    io_info.flag_allow_tri4_envelope = true;
    io_info.envelope_quad_tri_method = ENVELOPE_MAX_MIN_ANGLE;
    break;

  case ENVELOPE_ONLY_TRI2_OPT:
    io_info.flag_check_envelope = true;
    io_info.flag_allow_tri2_envelope = true;
    io_info.flag_allow_tri4_envelope = false;
    io_info.envelope_quad_tri_method = ENVELOPE_ONLY_TRI2;
    break;

  case ENVELOPE_TRI2_OR_TRI4_PREFER_TRI2_OPT:
    io_info.flag_check_envelope = true;
    io_info.flag_allow_tri2_envelope = true;
    io_info.flag_allow_tri4_envelope = true;
    io_info.envelope_quad_tri_method = ENVELOPE_PREFER_TRI2;
    break;

  case ENVELOPE_ONLY_TRI4_OPT:
    io_info.flag_check_envelope = true;
    io_info.flag_allow_tri2_envelope = false;
    io_info.flag_allow_tri4_envelope = true;
    io_info.envelope_quad_tri_method = ENVELOPE_ONLY_TRI4;
    break;

  case ENVELOPE_OFFSET_OPT:
    {
      const COORD_TYPE envelope_offset =
        get_arg_float(iarg, argc, argv, error);
      io_info.envelope_offset = envelope_offset;
      iarg++;
    }
    break;
    
  case ENVELOPE_USE_RELATIVE_POSITIONS_OPT:
    io_info.envelope_flag_use_relative_positions = true;
    break;
    
  case SEPARATE_ISOV_OPT:
    io_info.SetFlagSeparateIsovByCubeCenter(true);
    io_info.SetFlagSeparateIsovNearCubeFacets(true);
    io_info.SetFlagMoveAllIsovAwayFromCubeRidges(true);
    break;

  case MIN_DISTANCE_TO_CUBE_FACE_OPT:
    {
      const COORD_TYPE min_distance =
        get_arg_float(iarg, argc, argv, error);
      io_info.SetMinDistanceToCubeFace(min_distance);
      io_info.SetFlagMoveAllIsovAwayFromCubeFacets(true);
      
      iarg++;
    }
    break;
    
  case POSITION_OFFSET_OPT:
    {
      const COORD_TYPE position_offset =
        get_arg_float(iarg, argc, argv, error);
      io_info.SetPositionOffset(position_offset);

      // *** SHOULD THIS BE SET HERE TO FALSE??? ***
      io_info.SetFlagMoveAllIsovAwayFromCubeFacets(false);
      
      iarg++;
    }
    break;

  case NO_POSITION_OFFSET_OPT:
    io_info.SetFlagMoveAllIsovAwayFromCubeFacets(false);
    io_info.SetFlagMoveAllIsovAwayFromCubeRidges(false);
    io_info.SetFlagSeparateIsovNearCubeFacets(false);
    io_info.SetFlagMoveAllIsovAwayFromCubeRidges(false);
    break;
    
  case SEP_ISOV_NEAR_CUBE_FACETS_OPT:
    io_info.SetFlagSeparateIsovNearCubeFacets(true);
    io_info.SetFlagMoveAllIsovAwayFromCubeRidges(true);
    break;

  case MOVE_ALL_ISOV_AWAY_FROM_CUBE_FACETS_OPT:
    io_info.SetFlagMoveAllIsovAwayFromCubeFacets(true);
    break;

  case SEPARATE_ISOV_BY_CUBE_CENTER_OPT:
    io_info.SetFlagSeparateIsovByCubeCenter(true);
    break;

  case TRIMESH_UNIFORM_OPT:
    io_info.mesh_type = SIMPLICIAL_COMPLEX;
    io_info.quad_tri_method = UNIFORM_TRI;
    io_info.flag_trimesh = true;
    break;

  case TRI4_IN_ENVELOPE_OPT:
    io_info.interior_vertex_param.interior_vertex_position_method =
      INTERIOR_VERTEX_IN_ENVELOPE;
    io_info.is_tri4_position_method_set = true;
    break;

  case TRI4_ON_GRID_EDGE_OPT:
    io_info.interior_vertex_param.interior_vertex_position_method =
      INTERIOR_VERTEX_ON_GRID_EDGE;
    io_info.is_tri4_position_method_set = true;
    break;
    
  case TRI4_CENTROID_OPT:
    io_info.interior_vertex_param.interior_vertex_position_method =
      INTERIOR_VERTEX_AT_CENTROID;
    io_info.is_tri4_position_method_set = true;
    break;

  case QEI_INTERPOLATE_SCALAR_OPT:
    io_info.interior_vertex_param.interior_vertex_position_method =
      INTERIOR_VERTEX_ON_GRID_EDGE;
    io_info.interior_vertex_param.poly_edge_intersection_method =
      INTERPOLATE_EDGE_ENDPOINT_SCALARS;
    io_info.is_tri4_position_method_set = true;    
    io_info.is_qei_method_set = true;
    break;

  case QEI_BILINEAR_OPT:
    io_info.interior_vertex_param.interior_vertex_position_method =
      INTERIOR_VERTEX_ON_GRID_EDGE;    
    io_info.interior_vertex_param.poly_edge_intersection_method =
      POLY_EDGE_MULTILINEAR_INTERPOLATION;
    io_info.is_tri4_position_method_set = true;    
    io_info.is_qei_method_set = true;
    break;

  case QEI_AVERAGE_OPT:
    io_info.interior_vertex_param.interior_vertex_position_method =
      INTERIOR_VERTEX_ON_GRID_EDGE;    
    io_info.interior_vertex_param.poly_edge_intersection_method =
      POLY_EDGE_AVERAGE_PROJECTION;
    io_info.is_tri4_position_method_set = true;    
    io_info.is_qei_method_set = true;
    break;

  case INTERIOR_VERTEX_RATIO_OPT:
    {
      io_info.interior_vertex_param.interior_vertex_position_method =
        INTERIOR_VERTEX_IN_ENVELOPE;    
      io_info.is_tri4_position_method_set = true;    
      const float R = get_arg_float(iarg, argc, argv, error);
      io_info.interior_vertex_param.interior_vertex_ratio_to_cube_facet = R;
      iarg++;
      break;
    }
    
  case FORCE_ADD_ALL_INTERIOR_VERTICES_OPT:
    io_info.flag_force_add_all_interior_vertices = true;
    break;

#ifdef INCLUDE_RANDOM

  case RANDOM_SEED_OPT:
    random_seed =
      get_arg_int(iarg, argc, argv, error);
    io_info.SetRandomSeed(random_seed);
    io_info.is_random_seed_set = true;
    iarg++;
    break;

  case RANDOM_POS_OPT:
    random_seed =
       get_arg_int(iarg, argc, argv, error);
    io_info.SetRandomSeed(random_seed);
    io_info.SetVertexPositionMethod(RANDOM_ISOV_POS);
    io_info.is_random_seed_set = true;
    iarg++;
    break;

  case RANDOM_POS_DISTRIBUTION_OPT:
    iarg++;
    if (iarg >= argc) usage_error();
    random_distribution = get_random_distribution(argv[iarg]);
    io_info.SetRandomPositionDistribution(random_distribution);
    break;

  case RANDOM_SAMPLE_ENTIRE_CUBE_OPT:
    io_info.SetRPosFlagGenIsovApplyOffset(false);
    io_info.SetRPosFlagGenIsovSeparatedByCubeCenter(false);
    break;

  case RANDOM_SAMPLE_CUBE_BOUNDARY_OPT:
    io_info.SetRPosFlagGenIsovApplyOffset(false);
    io_info.SetRPosFlagGenIsovSeparatedByCubeCenter(false);
    io_info.SetRPosFlagGenIsovOnCubeBoundary(true);
    break;


  case RANDOM_SAMPLE_BOUNDARY_WIDTH_OPT:
    io_info.SetRPosFlagGenIsovApplyOffset(false);
    io_info.SetRPosFlagGenIsovSeparatedByCubeCenter(false);
    io_info.SetRPosFlagGenIsovOnCubeBoundary(true);
    {
      const COORD_TYPE width =
        get_arg_float(iarg, argc, argv, error);
      io_info.SetRPosBoundaryWidth(width);
      iarg++;
    }
    break;

#else

  case RANDOM_SEED_OPT:
    cerr << "Error. Flag -random_seed not implemented in this compiled version of ijkdual." << endl;
    usage_error();
    break;

  case RANDOM_POS_OPT:
    cerr << "Error. Flag -random_pos not implemented in this compiled version of ijkdual." << endl;
    usage_error();
    break;

#endif

  case OFF_OPT:
    io_info.flag_output_off = true;
    io_info.is_file_format_set = true;
    break;

  case PLY_OPT:
    io_info.flag_output_ply = true;
    io_info.is_file_format_set = true;
    break;

  case FIG_OPT:
    io_info.flag_output_fig = true;
    io_info.is_file_format_set = true;
    break;

  case OUTPUT_FILENAME_OPT:
    iarg++;
    if (iarg >= argc) usage_error();
    io_info.output_filename = argv[iarg];
    break;

  case OUTPUT_FILENAME_PREFIX_OPT:
    iarg++;
    if (iarg >= argc) usage_error();
    io_info.output_filename_prefix = argv[iarg];
    break;

  case STDOUT_OPT:
    io_info.flag_use_stdout = true;
    break;

  case LABEL_WITH_ISOVALUE_OPT:
    io_info.flag_label_with_isovalue = true;
    break;

  case NO_WRITE_OPT:
    io_info.flag_nowrite = true;
    break;

  case SILENT_OPT:
    io_info.flag_silent = true;
    break;

  case TERSE_OPT:
    io_info.flag_terse = true;
    break;    

  case NO_WARN_OPT:
    io_info.flag_no_warn = true;
    break;

  case TIME_OPT:
    io_info.flag_report_time = true;
    break;

  case DETAILED_TIME_OPT:
    io_info.flag_report_time = true;
    io_info.flag_report_time_detailed = true;
    break;

  case INFO_OPT:
    io_info.flag_report_info = true;
    break;

  case USAGE_OPT:
    usage(cout, 0);
    break;

  case MORE_OPTIONS_OPT:
    print_options(cout, EXTENDED_OPTG, 0);
    break;

  case TESTING_OPTIONS_OPT:
    print_options(cout, TESTING_OPTG, 0);
    break;

  case ALL_OPTIONS_OPT:
    usage_all(cout, 0);
    break;

  case HELP_OPT:
    help();
    break;

  case HELP_MORE_OPT:
    print_options_help(EXTENDED_OPTG);
    break;

  case HELP_TESTING_OPT:
    print_options_help(TESTING_OPTG);
    break;

  case HELP_ALL_OPT:
    help_all();
    break;

  case OUT_ISOV_OPT:
    iarg++;
    if (iarg >= argc) usage_error();
    io_info.report_isov_filename = argv[iarg];
    io_info.flag_report_all_isov = true;
    break;

  case NO_COMMENTS_OPT:
    io_info.flag_no_comments = true;
    break;

  default:
    return(false);
  };

  return(true);
}


/*!
 *  @brief Store isovalues and input filename in io_info.
 */
template <typename IO_INFO_TYPE>
void process_isovalues_and_input_filename
(const int argc, char **argv, const int iarg, IO_INFO_TYPE & io_info)
{
  // check for more parameter tokens
  for (int j = iarg; j+1 < argc; j++) {
    COMMAND_LINE_OPTION_TYPE optA;
    if (options.GetOption(argv[j], optA) ||
        (argv[j][0] == '-' && !is_type<float>(argv[j]))) {
      // argv[iarg] is not an isovalue
      cerr << "Usage error. Illegal parameter: " << argv[iarg] << endl;
      cerr << endl;
      usage_error();
    }
  }

  if (iarg+2 > argc) {
    cerr << "Error.  Missing input isovalue or input file name." << endl;
    cerr << endl;
    usage_error();
  };

  // store isovalues
  for (int j = iarg; j+1 < argc; j++) {
    SCALAR_TYPE value;
    if (!IJK::string2val(argv[j], value)) {
      cerr << "Error. \"" << argv[j] << "\" is not a valid input isovalue." 
           << endl;
      usage_error();
    }

    io_info.isovalue_string.push_back(argv[j]);
    io_info.isovalue.push_back(value);
  }

  io_info.input_filename = argv[argc-1];
}


/// Check io_info for parameter conflicts.
template <typename IO_INFO_TYPE>
void check_io_info(IO_INFO_TYPE & io_info)
{

  if (io_info.flag_subsample && io_info.flag_supersample) {
    cerr << "Usage error.  Can't use both -subsample and -supersample parameters."
         << endl;
    exit(64);
  }

  if (io_info.interior_vertex_param.interior_vertex_ratio_to_cube_facet >= 0.5) {
    cerr << "Usage error. Interior vertex ratio must be less than 0.5."
         << endl;
    exit(64);   
  }

  if (io_info.interior_vertex_param.interior_vertex_ratio_to_cube_facet < 0.0) {
    cerr << "Usage error. Interior vertex ratio must be non-negative."
         << endl;
    exit(64);     
  }
  
}


/*!
 *  @brief Process io_info.
 *  - Check consistency of options.
 *  - Some options are set based on others.
 *  - Some inconsistent options generate usage error messages and exit().
 */
template <typename IO_INFO_TYPE>
void process_io_info(IO_INFO_TYPE & io_info)
{
  if (io_info.flag_subsample && io_info.subsample_resolution <= 1) {
    cerr << "Usage error.  Subsample resolution must be an integer greater than 1."
         << endl;
    exit(64);
  };

  if (io_info.output_filename != "" && io_info.flag_use_stdout) {
    cerr << "Usage error.  Can't use both -o and -stdout parameters."
         << endl;
    exit(64);
  };

  if (io_info.isovalue.size() > 1 && io_info.output_filename != "") {
    cerr << "Usage error.  Cannot specify output file when input contains"
         << endl;
    cerr << "  more than one isovalue." << endl;
    exit(64);
  }

  if (!io_info.is_file_format_set) {
    if (io_info.IsOutputFilenameSet()) {
      const IO_INFO::OUTPUT_FORMAT output_format =
        io_info.GetFileType
          (io_info.output_filename, OUTPUT_DEFINITIONS::OFF);
      io_info.SetOutputFormat(output_format);
      io_info.SetOutputFilename(output_format, io_info.output_filename);
    }
    else {
      io_info.flag_output_off = true;
    }
  }
  else {
    if (io_info.IsOutputFilenameSet()) {
      if (io_info.NumOutputFormats() > 1) {
        cerr << "Usage error. Cannot set output filename when more than one"
             << endl
             << "  output format is selected." << endl;
        exit(64);
      }

      io_info.SetOutputFilename(io_info.output_filename);
    }
  }

  if (io_info.flag_split_non_manifold) {
    // Not really necessary, but...
    io_info.flag_split_non_manifold_ambig = true;
  }

  check_io_info(io_info);
}


// parse command line
// control parameters, followed by one or more isovalues, 
// followed by input file name
void IJKDUAL::parse_command_line(int argc, char **argv, IO_INFO & io_info)
{
  IJK::ERROR error;

  set_command_line_options(io_info);

  // Set defaults.
  io_info.quad_tri_method = MAX_MIN_ANGLE;

  if (argc == 1) { usage_error(); };

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    COMMAND_LINE_OPTION_TYPE optA;
    if (!options.GetOption(argv[iarg], optA)) {
      // Unknown parameter.  Possibly negative scalar value.
      break;
    }

    if (!process_option(optA, argc, argv, iarg, io_info)) {
      error.AddMessage
        ("Programming error. Option ", options.Option(optA).option_name, 
         " not implemented.");
      throw error;
    }

    iarg++;
  }

  // remaining parameters should be list of isovalues followed
  // by input file name

  process_isovalues_and_input_filename(argc, argv, iarg, io_info);
  process_io_info(io_info);

  if (io_info.mesh_type == CUBE_COMPLEX) {
    if (io_info.flag_check_envelope) {
      io_info.mesh_type = MIXED_MESH;
      io_info.quad_tri_method = TRI_DIAGONALS_OUTSIDE_ENVELOPE;
    }
  }


  if (io_info.VertexPositionMethod() == RANDOM_ISOV_POS) {
    if (io_info.is_qei_method_set) {
      if (io_info.interior_vertex_param.poly_edge_intersection_method ==
          IJK::INTERPOLATE_EDGE_ENDPOINT_SCALARS) {
        cerr << "Error.  Can't use both -random_pos (or variations) and -qei_interpolate_scalar."
             << endl;
        exit(64);
      }
    }
    else if (io_info.interior_vertex_param.poly_edge_intersection_method ==
             IJK::INTERPOLATE_EDGE_ENDPOINT_SCALARS) {
      // Change quad_edge_intersection_method.
      io_info.interior_vertex_param.poly_edge_intersection_method =
        IJK::POLY_EDGE_MULTILINEAR_INTERPOLATION;
    }
  }

  if (io_info.VertexPositionMethod() == RANDOM_ISOV_POS) {
    if (!io_info.is_random_seed_set) {
      cerr << "Usage error. Option \"-random_seed <seed>\" required" << endl;
      cerr << "  with random vertex position method." << endl;
      usage_error();
    }

    if (io_info.RandomPositionDistribution() == UNDEFINED_DISTRIBUTION) {
      // Default to uniform.
      io_info.SetRandomPositionDistribution(UNIFORM_DISTRIBUTION);
    }

  }

  if (io_info.is_random_seed_set) {
    if (io_info.VertexPositionMethod() != RANDOM_ISOV_POS) {
      cerr << "Warning. Vertex position method is not random." << endl;
      cerr << "  Ignoring option \"-random_seed <seed>\"." << endl;
      cerr << endl;
    }
  }

  if (io_info.is_tri4_position_method_set) {
    if (!io_info.flag_tri4_quad &&
        (!io_info.flag_check_envelope ||
         !io_info.flag_allow_tri4_envelope)) {
      cerr << "Warning. Vertex position method for interior vertex is set" << endl
           << "  but triangulation into four triangles is not enabled." << endl;
      cerr << "  Interior vertex is used to triangulate into 4 triangles." << endl;
      cerr << "  Use option -tri4 or -envelope to enable triangulation" << endl
           << "  into four triangles under various circumstances." << endl;
      cerr << endl;
    }
  }

  if (io_info.MinDistanceToCubeFace() > io_info.PositionOffset()) {
    cerr << "Warning. Minimum distance to grid cube face is greater than position offset."
         << endl;
    cerr << "  Resetting position offset to equal minimum distance to grid cube face."
         << endl;
    cerr << endl;
    io_info.SetPositionOffset(io_info.MinDistanceToCubeFace());
  }

}


// Check input information/flags.
void IJKDUAL::check_and_reset_input
(const DUALISO_SCALAR_GRID_BASE & scalar_grid,
 IO_INFO & io_info)
{
  const int DIM3(3);
  
  if (io_info.isovalue.size() > 1 && io_info.flag_use_stdout) {
    cerr << "Usage error. Cannot use stdout for more than one isovalue.";
    exit(64);
  }

  if (io_info.flag_output_ply) {
    if (scalar_grid.Dimension() != DIM3) {
      cerr << "Usage error. Cannot use option -ply on scalar grids with dimension other than 3.";
      exit(64);
    }
  }

  if (io_info.flag_output_fig) {
    if (scalar_grid.Dimension() != DIM2) {
      cerr << "Usage error. Cannot use option -fig on scalar grids with dimension other than 2.";
      exit(64);
    }
  }

  if (!io_info.flag_no_warn) {
    if (io_info.flag_trimesh) {
      if (scalar_grid.Dimension() == DIM2) {
        cerr << "Warning. Cannot triangulate 1 dimensional isocontour." << endl;
        cerr << "  Ignoring triangulation flags." << endl;
        cerr << endl;

        io_info.mesh_type = CUBE_COMPLEX;
        io_info.flag_trimesh = false;
      }
      else if (scalar_grid.Dimension() != DIM3) {
        cerr << "Warning. Triangulation not implemented on scalar grids with dimension other than 3." << endl;
        cerr << "  Ignoring triangulation flags." << endl;
        cerr << endl;

        io_info.mesh_type = CUBE_COMPLEX;
        io_info.flag_trimesh = false;
      }
    }

    if (io_info.flag_check_envelope) {
      if (scalar_grid.Dimension() == DIM2) {
        cerr << "Warning. No need to check envelope on 1 dimensional isocontours."
             << endl;
        cerr << "  Ignoring envelope flags." << endl;
        cerr << endl;

        io_info.mesh_type = CUBE_COMPLEX;
        io_info.flag_check_envelope = false;
      }
      else if (scalar_grid.Dimension() != DIM3) {
        cerr << "Warning. Check envelope not implemented on scalar grids with dimension other than 3."
             << endl;
        cerr << "  Ignoring envelope flags." << endl;
        cerr << endl;

        io_info.mesh_type = CUBE_COMPLEX;
        io_info.flag_check_envelope = false;
      }
    }

    if (scalar_grid.Dimension() > DIM3) {
      if (io_info.iso_vertex_position_param.flag_separate_isov_near_cube_facets) {
        cerr << "Warning. Dimension > 3." << endl;
        cerr << "  Separate isov from cube faces not fully implemented"
             << endl;
        cerr << "  for dimension > 3." << endl;
        cerr << endl;
      }
    }
  }
}


// ******************************************************************
// WRITE_DUAL_MESH
// ******************************************************************

// Write dual mesh.
void IJKDUAL::write_dual_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist)
{
  const int DIM1(1);
  const int DIM2(2);
  const int DIM3(3);
  std::vector<std::string> comment;

  if (!output_info.flag_no_comments) {
    if (output_info.mesh_dimension == DIM1) {
      comment.push_back(string("Dual contouring isocontour."));
    }
    else if (output_info.mesh_type == SIMPLICIAL_COMPLEX) {
      comment.push_back(string("Dual contouring triangulated isosurface."));
    }
    else if (output_info.mesh_type == CUBE_COMPLEX) {
      if (output_info.mesh_dimension == DIM2) {
        comment.push_back(string("Dual contouring (quadrilateral) isosurface."));
      }
      else if (output_info.mesh_dimension == DIM3) {
        comment.push_back(string("Dual contouring (hexahedral) isosurface."));
      }
      else {
        comment.push_back(string("Dual contouring (hypercube) isosurface."));
      }
    }
    else {
      comment.push_back(string("Dual contouring isosurface."));
    }
    add_meshfile_comments(output_info, comment);
  }

  write_mesh(output_info, vertex_coord, plist, comment);
}


// Write dual mesh and record output time.
void IJKDUAL::write_dual_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & plist,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_dual_mesh(output_info, vertex_coord, plist);

  io_time.write_time += wall_time.getElapsed();
}


// Write dual isosurface mesh of quad and triangles.
// @param output_info Output information.
// @param vertex_coord List of vertex coordinates.
void IJKDUAL::write_dual_quad_tri_mesh
(const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const std::vector<VERTEX_INDEX> & tri_vert)
{
  const int NUMV_PER_QUAD = 4;
  const int NUMV_PER_TRI = 3;
  const int dimension = output_info.dimension;
  const bool flag_use_stdout = output_info.flag_use_stdout;
  std::vector<VERTEX_INDEX> quad_vert2;
  std::vector<std::string> comment;
  ofstream output_file;
  string ofilename;
  PROCEDURE_ERROR error("write_dual_quad_tri_mesh");

  if (dimension != 3) {
    error.AddMessage("Programming error.  Illegal dimension ", dimension, ".");
    error.AddMessage("   Routine only allowed for dimension 3.");
    throw error;
  }


  quad_vert2.resize(quad_vert.size());
  std::copy(quad_vert.begin(), quad_vert.end(), quad_vert2.begin());

  IJK::reorder_quad_vertices(quad_vert2);

  comment.push_back(string("Dual contouring isosurface."));
  add_meshfile_comments(output_info, comment);

  switch (output_format) {

    case OUTPUT_DEFINITIONS::OFF:
      if (!flag_use_stdout) {
        ofilename = output_info.output_off_filename;
        output_file.open(ofilename.c_str(), ios::out);
        ijkoutOFF(output_file, dimension, vertex_coord, 
                  quad_vert2, NUMV_PER_QUAD, tri_vert, NUMV_PER_TRI,
                  comment);

        output_file.close();
      }
      else {
        throw error("Output stdout not implemented.");
      };
      break;

    case OUTPUT_DEFINITIONS::PLY:
      if (dimension == 3) {
        if (!flag_use_stdout) {
          ofilename = output_info.output_ply_filename;
          output_file.open(ofilename.c_str(), ios::out);
          ijkoutPLY(output_file, dimension, vertex_coord,
                    quad_vert2, NUMV_PER_QUAD, tri_vert, NUMV_PER_TRI);
          output_file.close();
        }
        else {
          ijkoutPLY(cout, dimension, vertex_coord,
                    quad_vert2, NUMV_PER_QUAD, tri_vert, NUMV_PER_TRI);
        }
      }
      else throw error
        ("Illegal dimension. PLY format is only for dimension 3.");
    break;

    default:
      throw error("Illegal output format.");
      break;
  }

  if (!flag_use_stdout && !output_info.flag_silent && ofilename != "")
    cout << "Wrote output to file: " << ofilename << endl;
}


// Write dual isosurface mesh of quad and triangles.
// @param output_info Output information.
// @param vertex_coord List of vertex coordinates.
void IJKDUAL::write_dual_quad_tri_mesh
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const std::vector<VERTEX_INDEX> & tri_vert)
{
  IJK::PROCEDURE_ERROR error("write_dual_quad_tri_mesh");

  if (output_info.flag_output_off) {
    if (output_info.output_off_filename != "") {
      write_dual_quad_tri_mesh
        (output_info, OUTPUT_DEFINITIONS::OFF, vertex_coord, quad_vert, tri_vert);
    }
    else {
      error.AddMessage("Programming error. Geomview OFF file name not set.");
      throw error;
    }
  }

  if (output_info.flag_output_ply) {
    if (output_info.output_ply_filename != "") {
      write_dual_quad_tri_mesh
        (output_info, OUTPUT_DEFINITIONS::PLY, vertex_coord, quad_vert, tri_vert);
    }
    else {
      error.AddMessage("Programming error. PLY file name not set.");
      throw error;
    }
  }

}


void IJKDUAL::write_dual_quad_tri_mesh
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const std::vector<VERTEX_INDEX> & tri_vert,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_dual_quad_tri_mesh(output_info, vertex_coord, quad_vert, tri_vert);

  io_time.write_time += wall_time.getElapsed();
}


// Write dual isosurface triangular mesh, color vertices
// @param output_info Output information.
// @param vertex_coord List of vertex coordinates.
// @param tri_vert[] List of triangle vertices.
//        tri_vert[3*i+k] is k'th vertex of triangle i.
void IJKDUAL::write_dual_tri_mesh_color_vertices
(const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & tri_vert,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
  const int NUMV_PER_TRI = 3;
  const int dimension = output_info.dimension;
  const bool flag_use_stdout = output_info.flag_use_stdout;
  ofstream output_file;
  string ofilename;
  PROCEDURE_ERROR error("write_dual_tri_mesh");

  if (dimension != 3) {
    error.AddMessage("Programming error.  Illegal dimension ", dimension, ".");
    error.AddMessage("   Routine only allowed for dimension 3.");
    throw error;
  }

  switch (output_format) {

    case OUTPUT_DEFINITIONS::OFF:
      if (!flag_use_stdout) {
        ofilename = output_info.output_off_filename;
        output_file.open(ofilename.c_str(), ios::out);
        ijkoutColorVertOFF(output_file, dimension, NUMV_PER_TRI,
                           vertex_coord, tri_vert, front_color, back_color);
        output_file.close();
      }
      else {
        ijkoutColorVertOFF(std::cout, dimension, NUMV_PER_TRI,
                           vertex_coord, tri_vert, front_color, back_color);
      };
      break;

    default:
      throw error("Illegal output format.");
      break;
  }

  if (!flag_use_stdout && !output_info.flag_silent && ofilename != "")
    cout << "Wrote output to file: " << ofilename << endl;
}


// Write dual isosurface triangular mesh, color vertices
// @param output_info Output information.
// @param vertex_coord List of vertex coordinates.
// @param tri_vert[] List of triangle vertices.
//        tri_vert[3*i+k] is k'th vertex of triangle i.
void IJKDUAL::write_dual_tri_mesh_color_vertices
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & tri_vert,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
  IJK::PROCEDURE_ERROR error("write_dual_tri_mesh_color_vertices");

  if (output_info.flag_output_off) {
    if (output_info.output_off_filename != "") {
      write_dual_tri_mesh_color_vertices
        (output_info, OUTPUT_DEFINITIONS::OFF, vertex_coord, tri_vert, front_color, back_color);
    }
    else {
      error.AddMessage("Programming error. Geomview OFF file name not set.");
      throw error;
    }
  }

  if (output_info.flag_output_ply) {
    cerr << "Warning: Writing tri mesh with colored vertices not"
         << endl
         << "  implemented for .ply format." 
         << endl;
    cerr << "Skipping output of .ply format." << endl;
  }

}

void IJKDUAL::write_dual_tri_mesh_color_vertices
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & tri_vert,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_dual_tri_mesh_color_vertices
    (output_info, vertex_coord, tri_vert, front_color, back_color);

  io_time.write_time += wall_time.getElapsed();
}


// Write dual isosurface mesh of quad and triangles.  Color vertices.
// @param output_info Output information.
// @param vertex_coord List of vertex coordinates.
void IJKDUAL::write_dual_quad_tri_mesh_color_vertices
(const OUTPUT_INFO & output_info, const OUTPUT_FORMAT output_format,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const std::vector<VERTEX_INDEX> & tri_vert,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
  const int NUMV_PER_QUAD = 4;
  const int NUMV_PER_TRI = 3;
  const int dimension = output_info.dimension;
  const bool flag_use_stdout = output_info.flag_use_stdout;
  std::vector<VERTEX_INDEX> quad_vert2;
  ofstream output_file;
  string ofilename;
  PROCEDURE_ERROR error("write_dual_quad_tri_mesh_color_vertices");

  if (dimension != 3) {
    error.AddMessage("Programming error.  Illegal dimension ", dimension, ".");
    error.AddMessage("   Routine only allowed for dimension 3.");
    throw error;
  }


  quad_vert2.resize(quad_vert.size());
  std::copy(quad_vert.begin(), quad_vert.end(), quad_vert2.begin());

  IJK::reorder_quad_vertices(quad_vert2);

  switch (output_format) {

    case OUTPUT_DEFINITIONS::OFF:
      if (!flag_use_stdout) {
        ofilename = output_info.output_off_filename;
        output_file.open(ofilename.c_str(), ios::out);
        ijkoutColorVertOFF
          (output_file, dimension, vertex_coord, 
           quad_vert2, NUMV_PER_QUAD, tri_vert, NUMV_PER_TRI,
           front_color, back_color);
        output_file.close();
      }
      else {
        throw error("Output stdout not implemented.");
      };
      break;

    default:
      throw error("Illegal output format.");
      break;
  }

  if (!flag_use_stdout && !output_info.flag_silent && ofilename != "")
    cout << "Wrote output to file: " << ofilename << endl;
}


// Write dual isosurface mesh of quad and triangles.  Color vertices.
// @param output_info Output information.
// @param vertex_coord List of vertex coordinates.
void IJKDUAL::write_dual_quad_tri_mesh_color_vertices
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const std::vector<VERTEX_INDEX> & tri_vert,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
  IJK::PROCEDURE_ERROR error("write_dual_quad_tri_mesh_color_vertices");

  if (output_info.flag_output_off) {
    if (output_info.output_off_filename != "") {
      write_dual_quad_tri_mesh_color_vertices
        (output_info, OUTPUT_DEFINITIONS::OFF, vertex_coord, quad_vert, tri_vert,
         front_color, back_color);
    }
    else {
      error.AddMessage("Programming error. Geomview OFF file name not set.");
      throw error;
    }
  }

  if (output_info.flag_output_ply) {
    cerr << "Warning: Writing quad/tri mesh with colored vertices not"
         << endl
         << "  implemented for .ply format." 
         << endl;
    cerr << "Skipping output of .ply format." << endl;
  }

}


void IJKDUAL::write_dual_quad_tri_mesh_color_vertices
(const OUTPUT_INFO & output_info,
 const std::vector<COORD_TYPE> & vertex_coord,
 const std::vector<VERTEX_INDEX> & quad_vert,
 const std::vector<VERTEX_INDEX> & tri_vert,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  write_dual_quad_tri_mesh_color_vertices
    (output_info, vertex_coord, quad_vert, tri_vert, front_color, back_color);

  io_time.write_time += wall_time.getElapsed();
}


// ******************************************************************
// RESCALE ROUTINES
// ******************************************************************

namespace {

  void grow_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
  {
    for (unsigned int i = 0; i < vertex_coord.size(); i++) {
      vertex_coord[i] = scale * vertex_coord[i];
    };
  }

  void shrink_coord(const int scale, vector<COORD_TYPE> & vertex_coord)
  {
    for (unsigned int i = 0; i < vertex_coord.size(); i++) {
      vertex_coord[i] = vertex_coord[i]/scale;
    };
  }

  bool unit_spacing(const int dimension, const COORD_TYPE * spacing)
  {
    for (int d = 0; d < dimension; d++) {
      if (!AIR_EXISTS(spacing[d])) { return(true); }
      else if (spacing[d] != 1.0) { return(false); };
    }

    return(true);
  }

}


void IJKDUAL::rescale_vertex_coord
(const int dimension, const COORD_TYPE * grid_spacing,
 std::vector<COORD_TYPE> & vertex_coord)
{
  if (unit_spacing(dimension, grid_spacing)) { return; }

  if (vertex_coord.size() == 0) { return; };

  const VERTEX_INDEX numv = vertex_coord.size()/dimension;
  for (int iv = 0; iv < numv; iv++) {
    for (int d = 0; d < dimension; d++) {
      vertex_coord[iv*dimension+d] *= grid_spacing[d];
    }
  }
}


void IJKDUAL::rescale_vertex_coord(const std::vector<COORD_TYPE> & grid_spacing,
                                   std::vector<COORD_TYPE> & vertex_coord)
{
  const int dimension = grid_spacing.size();
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (grid_spacing.size() < 1) {
    error.AddMessage("Illegal size ", grid_spacing.size(), 
                     " of array grid spacing.");
    error.AddMessage("Size must equal vertex dimension.");
    throw error;
  }

  rescale_vertex_coord(dimension, vector2pointer(grid_spacing),
                       vertex_coord);
}


// Rescale subsampled/supersampled vertex coordinates.
void IJKDUAL::rescale_vertex_coord
(const int subsample_resolution, const int supersample_resolution, COORD_ARRAY & vertex_coord)
{
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (subsample_resolution <= 0) {
    error.AddMessage("Illegal subsample resolution ", subsample_resolution, ".");
    error.AddMessage("  Subsample resolution must be a positive integer");
  }

  if (supersample_resolution <= 0) {
    error.AddMessage("Illegal supersample resolution ", supersample_resolution, ".");
    error.AddMessage("  Supersample resolution must be a positive integer");
  }

  if (vertex_coord.size() == 0) { return; };

  if (subsample_resolution != 1)
    { grow_coord(subsample_resolution, vertex_coord); };

  if (supersample_resolution != 1)
    { shrink_coord(supersample_resolution, vertex_coord); };
}


// Rescale subsampled/supersampled vertex coordinates.
// Also rescale to reflect grid spacing.
void IJKDUAL::rescale_vertex_coord
(const int subsample_resolution, const int supersample_resolution,
 const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord)
{
  rescale_vertex_coord(subsample_resolution, supersample_resolution, vertex_coord);
  rescale_vertex_coord(grid_spacing, vertex_coord);
}


// ******************************************************************
// CREATE MESH FILE COMMENTS
// ******************************************************************

void IJKDUAL::add_meshfile_comments
(const OUTPUT_INFO & output_info,
 std::vector<std::string> & comment)
{
  std::string str;

  if (output_info.flag_no_comments) { return; }

  output_info.AddIsovalueComment(comment);
  output_info.AddSubsampleResolutionComment(comment);
  output_info.AddSupersampleResolutionComment(comment);

  vpos2string(output_info.VertexPositionMethod(), str);
  str = "VertexPositionMethod: " + str;
  comment.push_back(str);

  if (output_info.allow_multiple_iso_vertices)
    { comment.push_back("IsoVerticesPerCube: Multiple"); }
  else
    { comment.push_back("IsoVerticesPerCube: Single"); }
  if (output_info.flag_check_envelope)
    { comment.push_back("CheckEnvelope: True"); }

  if (output_info.flag_separate_neg) {
    comment.push_back("IsosurfaceSeparates: NegativeGridVertices");
  }
  else {
    comment.push_back("IsosurfaceSeparates: PositiveGridVertices");
  }

  if (output_info.mesh_type == SIMPLICIAL_COMPLEX ||
      output_info.mesh_type == MIXED_MESH) {
    tri_method2string(output_info.quad_tri_method, str);
    str = "TriangulationMethod: " + str;
    comment.push_back(str);

    if (output_info.quad_tri_method == TRI_DIAGONALS_OUTSIDE_ENVELOPE) {
      if (output_info.flag_allow_tri4_envelope) {
        if (output_info.flag_allow_tri2_envelope)
          { comment.push_back("TriangulateQuads: TwoOrFourTriangles"); }
        else
          { comment.push_back("TriangulateQuads: OnlyFourTriangles"); }
      }
      else if (output_info.flag_allow_tri2_envelope)
        { comment.push_back("TriangulateQuads: OnlyTwoTriangles"); }
    }
    else if (output_info.flag_tri4_quad) {

      if (output_info.quad_tri_method == TRI4_ALL_QUADS)
        { comment.push_back("TriangulateQuads: OnlyFourTriangles"); }
      else
        { comment.push_back("TriangulateQuads: TwoOrFourTriangles"); }
    }
    else
      { comment.push_back("TriangulateQuads: OnlyTwoTriangles"); }
  }

  if (output_info.mesh_type == SIMPLICIAL_COMPLEX ||
      output_info.mesh_type == MIXED_MESH) {

    if (output_info.quad_tri_method == TRI4_ALL_QUADS ||
        output_info.quad_tri_method == TRI4_MAX_MIN_ANGLE ||
        (output_info.quad_tri_method == TRI_DIAGONALS_OUTSIDE_ENVELOPE &&
         output_info.flag_allow_tri4_envelope)) {

      splitv_pos2string(output_info.interior_vertex_param.interior_vertex_position_method, str);
      str = "SplitVertexPosition: " + str;
      comment.push_back(str);

      if (output_info.interior_vertex_param.PositionOnGridEdge()) {
        qei2string(output_info.interior_vertex_param.poly_edge_intersection_method, str);
        str = "QuadEdgeIntersection: " + str;
        comment.push_back(str);
      }
      else if (output_info.interior_vertex_param.PositionInEnvelope()) {
        const COORD_TYPE ratio =
          output_info.interior_vertex_param.interior_vertex_ratio_to_cube_facet;
        IJK::create_string("SplitVertexRatioToCubeFacet: ", ratio, str);
        comment.push_back(str);
      }
    }
  }
  

#ifdef INCLUDE_RANDOM
  if (output_info.VertexPositionMethod() == RANDOM_ISOV_POS) {
    IJK::create_string
      ("RandomSeed: ", output_info.RandomPositionSeed(), str);
    comment.push_back(str);
  }
#endif

  if (output_info.VertexPositionMethod() == RANDOM_ISOV_POS) {
    IJK::create_string
      ("PositionOffset: ", output_info.PositionOffset(), str);
    comment.push_back(str);
  }

}


// ******************************************************************
// CREATE STRING FROM enum TYPE
// ******************************************************************

void IJKDUAL::vpos2string(const VERTEX_POSITION_METHOD vpos,
                          std::string & vpos_string)
{
  switch (vpos) {

    case CUBE_CENTER:
      vpos_string = "CubeCenter";
      break;

    case CENTROID_EDGE_ISO:
      vpos_string = "CentroidEdge";
      break;

    case DIAGONAL_INTERPOLATION:
      vpos_string = "DiagonalInterpolation";
      break;

    case IVOL_LIFTED02:
      vpos_string = "IntervalVolumeLifted";
      break;

    case QDUAL_INTERPOLATION:
      vpos_string = "QualityDualInterpolation";
      break;

    case RANDOM_ISOV_POS:
      vpos_string = "Random";
      break;

    default:
        vpos_string = "Undefined";
        break;
  }

}


void IJKDUAL::tri_method2string
(const QUAD_TRI_METHOD quad_tri_method,
 std::string & method_string)
{
  switch(quad_tri_method) {

    case UNIFORM_TRI:
      method_string = "Uniform";
      break;

    case SPLIT_MAX_ANGLE:
      method_string = "SplitMaxAngle";
      break;

    case MAX_MIN_ANGLE:
    case TRI4_MAX_MIN_ANGLE:
      method_string = "MaxMinAngle";
      break;

    case TRI4_ALL_QUADS:
      method_string = "Tri4AllQuads";
      break;

    case TRI_DIAGONALS_OUTSIDE_ENVELOPE:
      method_string = "TriDiagonalsOutsideEnvelope";
      break;

    default:
      method_string = "Undefined";
  }

}


void IJKDUAL::splitv_pos2string
(const INTERIOR_VERTEX_POSITION_METHOD position_method,
 std::string & splitv_pos_string)
{
  switch (position_method) {

  case IJK::INTERIOR_VERTEX_ON_GRID_EDGE:
    splitv_pos_string = "OnGridEdge";
    break;

  case IJK::INTERIOR_VERTEX_IN_ENVELOPE:
    splitv_pos_string = "InEnvelope";
    break;

  case IJK::INTERIOR_VERTEX_AT_CENTROID:
    splitv_pos_string = "AtCentroid";
    break;

  default:
    splitv_pos_string = "Undefined";
  }

}


void IJKDUAL::qei2string
(const IJK::POLY_EDGE_INTERSECTION_METHOD qei_method,
 std::string & qei_string)
{
  switch (qei_method) {

  case IJK::POLY_EDGE_MULTILINEAR_INTERPOLATION:
    qei_string = "BilinearSurface";
    break;

  case INTERPOLATE_EDGE_ENDPOINT_SCALARS:
    qei_string = "InterpolateScalar";
    break;

  case POLY_EDGE_AVERAGE_PROJECTION:
    qei_string = "AverageProjection";
    break;

  default:
    qei_string = "Undefined";
  }
  
}

  
// ******************************************************************
// REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
// ******************************************************************

void IJKDUAL::report_num_cubes
(const DUALISO_GRID & full_scalar_grid, const IO_INFO & io_info, 
 const DUALISO_DATA & dualiso_data)
{
  report_num_cubes(full_scalar_grid, io_info, dualiso_data.ScalarGrid());
}


void IJKDUAL::report_num_cubes
(const DUALISO_GRID & full_scalar_grid, const IO_INFO & io_info, 
 const DUALISO_GRID & dualiso_data_grid)
{
  const int num_grid_cubes = full_scalar_grid.ComputeNumCubes();
  const int num_cubes_in_dualiso_data = dualiso_data_grid.ComputeNumCubes();

  if (!io_info.flag_use_stdout && !io_info.flag_silent &&
      !io_info.flag_terse) {

    if (io_info.flag_subsample) {
      // subsampled grid
      cout << num_grid_cubes << " grid cubes.  "
           << num_cubes_in_dualiso_data << " subsampled grid cubes." << endl;
    }
    else if (io_info.flag_supersample) {
      // supersample grid
      cout << num_grid_cubes << " grid cubes.  "
           << num_cubes_in_dualiso_data << " supersampled grid cubes." << endl;
    }
    else {
      // use full_scalar_grid
      cout << num_grid_cubes << " grid cubes." << endl;
    }
  }

}


void IJKDUAL::warn_non_manifold(const IO_INFO & io_info)
{
  const char * mesh_str = "mesh";

  if (io_info.flag_use_stdout || io_info.flag_silent ||
      io_info.flag_no_warn) { return; }

  if (io_info.isovalue.size() > 1) { mesh_str = "meshes"; }

  if (!io_info.allow_multiple_iso_vertices ||
      !io_info.flag_split_non_manifold) {

    cout << "*** Warning: Isosurface " << mesh_str
         << " may have non-manifold vertices or edges."   << endl;
    cout << endl;
  }
}


// ******************************************************************
// REPORT TIMING INFORMATION
// ******************************************************************

void IJKDUAL::report_dualiso_time
(const IO_INFO & io_info, const DUALISO_TIME & dualiso_time, 
 const char * mesh_type_string)
{
  cout << "CPU time to run ijkdual: " 
       << dualiso_time.total << " seconds." << endl;

  cout << "    Time to extract " << mesh_type_string << " polytopes: "
       << dualiso_time.extract << " seconds." << endl;
  cout << "    Time to merge identical "
       << mesh_type_string << " vertices: " 
       << dualiso_time.merge << " seconds." << endl;
  
  if (io_info.flag_report_time_detailed) {
    cout << "    Total time to position "
         << mesh_type_string << " vertices: "
         << dualiso_time.position.total << " seconds." << endl;
    cout << "      Basic (initial) time to position "
         << mesh_type_string << " vertices: "
         << dualiso_time.position.basic << " seconds." << endl;

    if (io_info.FlagSeparateIsovNearCubeFacets()) {
      cout << "      Time to separate vertices near cube facets: "
           << dualiso_time.position.separate_isov_near_cube_facets
           << " seconds." << endl;
    }

    if (io_info.FlagMoveAllIsovAwayFromCubeFacets()) {
      cout << "      Time to move vertices away from cube facets: "
           << dualiso_time.position.move_isov_away_from_cube_faces
           << " seconds." << endl;
    }

    if (io_info.FlagSeparateIsovByCubeCenter()) {
      cout << "      Time to separate vertices by cube centers: "
         << dualiso_time.position.separate_by_cube_center << " seconds." << endl;
    }

    if (io_info.flag_tri4_quad ||
        (io_info.flag_check_envelope &&
         io_info.flag_allow_tri4_envelope)) {
      cout << "      Time to create (position) vertices dual to isosurface polytopes: "
           << dualiso_time.position.isov_dual_to_isopoly << " seconds." << endl;
    }
  }
  else {
    cout << "    Time to position "
         << mesh_type_string << " vertices: "
         << dualiso_time.position.total << " seconds." << endl;
  }
  
  if (io_info.allow_multiple_iso_vertices) {
    cout << "    Time to split iso vertices: "
         << dualiso_time.split_isov << " seconds." << endl;
  }
  
  if (io_info.flag_check_envelope) {
    cout << "    Time to determine diagonals in envelopes: "
         << dualiso_time.envelope << " seconds." << endl;
  }

  if ((io_info.mesh_type == SIMPLICIAL_COMPLEX) ||
      (io_info.flag_check_envelope &&
       (io_info.flag_allow_tri2_envelope ||
        io_info.flag_allow_tri4_envelope))) {
      cout << "    Time to triangulate quadrilaterals: "
           << dualiso_time.triangulation << " seconds." << endl;
  }
  
  if (io_info.flag_dual_collapse) {
    cout << "    Time to merge " << mesh_type_string << " vertices: "
         << dualiso_time.collapse << " seconds." << endl;
  }
}


void IJKDUAL::report_time
(const IO_INFO & io_info, const IO_TIME & io_time, 
 const DUALISO_TIME & dualiso_time, const double total_elapsed_time)
{
  const char * ISOSURFACE_STRING = "isosurface";
  const char * mesh_type_string = NULL;
  
  mesh_type_string = ISOSURFACE_STRING;

  cout << "Time to read file " << io_info.input_filename << ": "
       << io_time.read_nrrd_time << " seconds." << endl;

  report_dualiso_time(io_info, dualiso_time, mesh_type_string);
  if (!io_info.flag_nowrite) {
    cout << "Time to write "
         << mesh_type_string << ": " 
         << io_time.write_time << " seconds." << endl;
  };
  cout << "Total elapsed time: " << total_elapsed_time
       << " seconds." << endl;
}


// ******************************************************************
// USAGE/HELP MESSAGES
// ******************************************************************

// local namespace
namespace {

  void usage_msg(std::ostream & out)
  {
    out << "Usage: ijkdual [OPTIONS] {isovalue1 isovalue2 ...} {input filename}" << endl;
  }

  void print_options_title
  (std::ostream & out,
   const COMMAND_LINE_OPTION_GROUP group)
  {
    switch(group) {
      
    case REGULAR_OPTG:
      out << "OPTIONS:" << endl;
      break;

    case EXTENDED_OPTG:
      out << "MORE OPTIONS:" << endl;
      break;

    case TESTING_OPTG:
      out << "TESTING OPTIONS:" << endl;
      break;

    default:
      out << "OTHER OPTIONS:" << endl;
      break;
    }
  }

  void options_msg(std::ostream & out,
		   const COMMAND_LINE_OPTION_GROUP group)
  {
    print_options_title(out, group);
    options.PrintUsageOptions(out, group);
  };

  void help_msg()
  {
    usage_msg(cout);
    cout << endl;
    cout << "ijkdual - Dual contouring isosurface generation algorithm." 
         << endl;
  }

  void print_help_options(const COMMAND_LINE_OPTION_GROUP group)
  {
    print_options_title(cout, group);

    for (std::size_t i = 0; i < options.list.size(); i++) {
      if (options.list[i].Group() == group)
        { options.PrintHelpMessage(cout, i); }
    }
  }
}

void IJKDUAL::usage(std::ostream & out, const int return_code)
{
  usage_msg(out);
  options_msg(out, REGULAR_OPTG);
  exit(return_code);
}

void IJKDUAL::usage_error()
{
  usage(cerr, 10);
}

void IJKDUAL::print_options
  (std::ostream & out, const COMMAND_LINE_OPTION_GROUP group,
   const int return_code)
{
  options_msg(out, group);
  exit(return_code);
}

void IJKDUAL::usage_all(std::ostream & out, const int return_code)
{
  usage_msg(out);
  options_msg(out, REGULAR_OPTG);
  options_msg(out, EXTENDED_OPTG);
  options_msg(out, TESTING_OPTG);
  exit(return_code);
}

void IJKDUAL::help_all()
{
  help_msg();
  cout << endl;

  print_help_options(REGULAR_OPTG);
  cout << endl;
  print_help_options(EXTENDED_OPTG);
  cout << endl;
  print_help_options(TESTING_OPTG);
  cout << endl;

  exit(0);
}

void IJKDUAL::print_options_help
  (const COMMAND_LINE_OPTION_GROUP group)
{
  print_help_options(group);
  exit(0);
}

void IJKDUAL::help()
{
  help_msg();
  cout << endl;
  print_help_options(REGULAR_OPTG);
  exit(0);
}


// ******************************************************************
// CLASS OUTPUT_INFO
// ******************************************************************

void IJKDUAL::OUTPUT_INFO::SetDimension(const int d)
{
  IJK::PROCEDURE_ERROR error("OUTPUT_INFO::SetDimension");

  dimension = d;

  if (flag_interval_volume)
    { mesh_dimension = dimension; }
  else
    { mesh_dimension = dimension - 1; }

  if (mesh_type == SIMPLICIAL_COMPLEX) {
    num_vertices_per_isopoly = mesh_dimension+1;
  }
  else {
    IJK::int_power(2, mesh_dimension, num_vertices_per_isopoly, error);
  }

}

namespace {

  void split_string(const string & s, const char c,
                    string & prefix, string & suffix)
  // split string at last occurrence of character c into prefix and suffix
  {
    string::size_type i = s.rfind(c);
    if (i == string::npos) {
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

}


// ******************************************************************
// SET ROUTINES
// ******************************************************************

// Set output information, copying information from io_info.
void IJKDUAL::set_output_info
(const IO_INFO & io_info, const int dimension, const int i,
 OUTPUT_INFO & output_info)
{
  // Note: SetDimension() should always be called AFTER output_info.Set().
  output_info.Set(io_info);
  output_info.SetDimension(dimension);

  output_info.output_isovalue.clear();
  output_info.output_isovalue.push_back(io_info.isovalue[i]);

  if (io_info.flag_interval_volume) {
    if (i+1 < int(io_info.isovalue.size()))
      { output_info.output_isovalue.push_back(io_info.isovalue[i+1]); }
  }

  if (!io_info.are_output_filenames_set) {
    output_info.ConstructOutputFilenames(i);
  }
}
