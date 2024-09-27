/*!
 *  @file ijkmcubeIO.cxx
 *  @brief IO routines for ijkmcube.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009-2024 Rephael Wenger

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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ijkmcubeIO.h"
#include "ijkmcube_util.h"
#include "ijkMCtable_xitIO.h"

#include "ijkcommand_line.tpp"
#include "ijkIO.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkstring.tpp"

// *** DEBUG ***
#include "ijkprint.tpp"

using namespace IJK;
using namespace IJKMCUBE;
using namespace IJKMCUBE_TABLE;

using namespace std;

// **************************************************
// PARSE COMMAND LINE
// **************************************************

// local namespace
namespace {

  typedef typename IJK::COMMAND_LINE_OPTIONS
    <COMMAND_LINE_OPTION_TYPE,COMMAND_LINE_OPTION_GROUP>
  COMMAND_LINE_OPTIONS;

  COMMAND_LINE_OPTIONS options;


  /// @brief Get high resolution coordinates.
  void get_high_resolution_arg
    (const int iarg, const int argc, char **argv, IO_INFO & io_info)
  {
    IJK::BOX<GRID_COORD_TYPE> box;
    vector<GRID_COORD_TYPE> coord;

    IJK::get_arg_multiple_arguments(iarg, argc, argv, coord);

    if (coord.size() == 1) {
      cerr << "Usage error in -highres coordinates: "
          << "\"" << argv[iarg+1] << "\"" << endl;
      cerr << "  Coordinate list must be contained in quotation marks." << endl;
      cerr << "  Number of coordinates must be at least two." << endl;
      exit(-1);
    }

    int box_dim = coord.size()/2;
    box.SetDimension(box_dim);

    for (int d = 0; d < box_dim; d++) {
      const GRID_COORD_TYPE minc = coord[d];
      const GRID_COORD_TYPE maxc = coord[d+box_dim];

      box.SetMinCoord(d, coord[d]);
      box.SetMaxCoord(d, coord[d+box_dim]);
    }

    io_info.high_resolution_region.push_back(box);
    io_info.high_resolution_string.push_back(argv[iarg+1]);

  }


  ISOSURFACE_TOPOLOGY get_topology(char * s)
    // convert string s into parameter token
  {
    ISOSURFACE_TOPOLOGY topology = ISOTABLE_TOPOLOGY;

    string str = s;

    if (str == "isotable") 
      { topology = ISOTABLE_TOPOLOGY; }
    else if (str == "adecider") 
      { topology = ASYMPTOTIC_DECIDER_TOPOLOGY; }
    else if (str == "cube_decider") 
      { topology = CUBE_DECIDER_TOPOLOGY; }
    else if (str == "linear") 
      { topology = LINEAR_TOPOLOGY; }
    else {
      cerr << "Error in input parameter -topology.  Illegal isosurface topology: " 
           << s << "." << endl;
      exit(10);
    }

    return(topology);
  }


  INTERPOLATION_TYPE get_interpolation_type(char * s)
    // convert string s into parameter token
  {
    INTERPOLATION_TYPE type = LINEAR_INTERPOLATION;

    if (strcmp(s, "linear") == 0) 
      { type = LINEAR_INTERPOLATION; }
    else if (strcmp(s, "multilinear") == 0) 
      { type = MULTILINEAR_INTERPOLATION; }
    else {
      cerr << "Error in input parameter -interpolate.  Illegal interpolation type: " 
           << s << "." << endl;
      exit(10);
    }

    return(type);
  }


  /// @brief Check that highres options have correct number of parameters.
  /// @param dimension Grid dimension.
  void check_highres
  (const int dimension, const IO_INFO & io_info)

  {
    using std::string;

    for (unsigned int i = 0; i < io_info.high_resolution_region.size(); i++) {
      if (io_info.high_resolution_region[i].Dimension() != dimension) {
        const string msg =
          string("Usage error in option: -highres \"") +
            io_info.high_resolution_string[i] + "\"";
        cerr << msg << endl;
        cerr << "Option -highres should be followed by a string containing "
             << 2*dimension << " coordinates." << endl;
        exit(-1);
      }
    }

    for (unsigned int i = 0; i < io_info.high_resolution_region.size(); i++) {
      for (int d = 0; d < dimension; d++) {
        const GRID_COORD_TYPE minc =
          io_info.high_resolution_region[i].MinCoord(d);
        const GRID_COORD_TYPE maxc =
          io_info.high_resolution_region[i].MaxCoord(d);
        if (minc > maxc) {
          cerr << "Usage error in -highres coordinates: "
               << "\"" << io_info.high_resolution_string[i]
               << "\"" << endl;
          cerr << "  Minimum coordinate " << minc
               << " (position " << d+1 << ")"
               << " is greater than maximum coordinate " << maxc
               << " (position " << d+dimension+1 << ")." << endl;
          cerr << "  List all minimum coordinates before maximum coordinates."
              << endl;
          exit(-1);
        }
      }
    }
  }

}


namespace {

  void set_command_line_options(const IO_INFO & io_info)
  {
    int iarg;

    options.SetLabelWidth(15);
    options.AddLabelTab(20);

    // Regular options.

    options.AddOptionNoArg
      (OCTREE_OPT, "OCTREE_OPT", REGULAR_OPTG,
       "-octree", "Create octree for faster isosurface extraction.");

    options.AddOption1Arg
      (REGION_OPT, "REGION_OPT", REGULAR_OPTG,
       "-region", "-L", "Preprocess into regions of dimension LxLx...");
    options.AddToHelpMessage
      (REGION_OPT, "for fast isosurface extraction.",
       "Each region has a minimum and maximum scalar value.");

    options.AddOptionNoArg
      (USE_CUBE_LIST_OPT, "USE_CUBE_LIST_OPT", REGULAR_OPTG,
       "-use_cube_list", "Preprocess by creating list of mixed cubes.");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOptionNoArg
      (NEP_OPT, "NEP_OPT", REGULAR_OPTG, "-nep",
       "Use negative, equals, positive isosurface lookup table.");
    options.AddToHelpMessage
      (NEP_OPT, "Use isosurface lookup table that differentiates scalar values",
       "less than (negative), equal to, and greater than (positive) the isovalue.");

    options.AddOption1Arg
      (NEP_DUP_OPT, "NEP_DUP_OPT", REGULAR_OPTG, "-nep_dup", "N",
       "Set parameter for processing duplicate isosurface triangles.");
    options.AddToHelpMessage
      (NEP_DUP_OPT,
       "The -nep option can cause the creation of isosurface patches",
       "lying entirely in cube facets.",
       "Two cubes sharing a facet could generate duplicate isosurface patches",
       "and duplicate simplices on that facet.");
    options.AddToHelpMessage
      (NEP_DUP_OPT,
       "-nep_dup 0: All duplicate isosurface triangles are deleted.",
       "(Duplicate triangles are replace by 0 triangles.)",
       "-nep_dup 2: Isosurface may contain duplicate isosurface triangles.",
       "(Allow 2 identical triangles in isosurface.)");
    options.AddToHelpMessageDefaultValue
      (NEP_DUP_OPT, "nep_dup", io_info.DefaultNEPnumDup());


    options.AddOption1Arg
      (SNAP_OPT, "SNAP_OPT", REGULAR_OPTG,
       "-snap", "D", "Snap isosurface vertices within distance D of grid vertices.");
    options.AddToHelpMessage
      (SNAP_OPT, "Value D must be in range [0.0, 0.5].");

    options.AddOptionNoArg
      (IVOL_OPT, "IVOL_OPT", REGULAR_OPTG, "-ivol",
       "Generate interval volume.");
    options.AddToHelpMessage
      (IVOL_OPT, "Requires at least two isovalues.");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOption1Arg
      (HIGHRES_OPT, "HIGHRES_OPT", REGULAR_OPTG,
       "-highres", "\"region coordinates\"",
       "High resolution region in multiresolution isosurface.");
    options.AddToHelpMessage
      (HIGHRES_OPT, "Coordinates must be surrounded by quotation marks.",
       "List coordinates of lowest point in region",
       "followed by coordinates of highest point.",
       "Option may be used multiple times.");

    options.AddOption1Arg
      (NUMRES_LEVELS_OPT, "NUMRES_LEVELS_OPT", REGULAR_OPTG,
       "-numres_levels", "N", "Number of resolution levels.");
    options.AddToHelpMessage
      (NUMRES_LEVELS_OPT,
       "Determines maximum resolution of high resolution regions.",
       "Each level represents a further subdivision of a (hyper) cube",
       "into 2^d subcubes.",
       "Number of resolution levels must be at least 2.");
    options.AddToHelpMessageDefaultValue
      (NUMRES_LEVELS_OPT, "number of resolution levels",
       io_info.num_resolution_levels);

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOption1Arg
      (TOPOLOGY_OPT, "TOPOLOGY_OPT", REGULAR_OPTG,
       "-topology", "{isotable|cube_decider|adecider}",
       "Construct isosurface with given topology type.");
    options.AddArgChoice
      (TOPOLOGY_OPT, "isotable",
       "Topology based on isosurface lookup table");
    options.AddArgChoice
      (TOPOLOGY_OPT, "cube_decider",
       "Topology based on resolving cube ambiguities.");
    options.AddArgChoice
      (TOPOLOGY_OPT, "adecider", "Asymptotic decider topology.");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOption1Arg
      (INTERPOLATE_OPT, "INTERPOLATE_OPT", REGULAR_OPTG,
       "-interpolate", "{linear|multilinear}",
       "Use linear or multilinear interpolation.");
    options.AddArgChoice
      (INTERPOLATE_OPT, "linear",
       "Use linear interpolation on scalar values of edge endpoints",
       "to determine isosurface edge intersections.");
    options.AddArgChoice
      (INTERPOLATE_OPT, "multilinear",
       "Use multilinear interpolation on scalar values of cube vertices",
       "to determine scalar values of new mesh vertices.",
       "Multilinear interpolation only works with asymptotic decider",
       "or linear topology.");

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

    /* Suppress for now.
     * ijkIO.tpp does not produce proper colored mesh files.
     */
    /*
    options.AddOptionNoArg
      (COLOR_ALTERNATING_OPT, "COLOR_ALTERNATING_OPT", REGULAR_OPTG,
       "-color_alternating", "Color alternating cubes. Works only with -off.");
    */

    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (EDGE_GROUPS_OPT, "EDGE_GROUPS_OPT", REGULAR_OPTG,
       "-edge_groups", "Use isosurface lookup table based on edge groups.");
    options.AddToHelpMessage
      (EDGE_GROUPS_OPT, 
       "Isosurface patch triangulation based on edge groups.",
       "(Default for cubes in dimension 3.)");
    options.AddToHelpMessage
      (EDGE_GROUPS_OPT, "See \"Edge Groups:",
       "An approach to understanding the mesh quality",
       "of Marching Cubes\" by Dietrich et al, IEEE TVCG, 2008.");

    options.AddOptionNoArg
      (CHULL_OPT, "CHULL_OPT", REGULAR_OPTG,
       "-chull", "Use isosurface lookup table based on convex hull projection.");
    options.AddToHelpMessage
      (CHULL_OPT, "(Default for all polygons/polytopes other than cubes",
       "or for all dimensions other than 3.)");
    options.AddToHelpMessage
      (CHULL_OPT, "See \"Isosurface construction in any dimension\"",
       "by Bhaniramka et al, IEEE TVCG, 2004.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (SEP_NEG_OPT, "SEP_NEG_OPT", REGULAR_OPTG,
       "-sep_neg", 
       "Use isosurface lookup table separating negative grid vertices.",
       "(Default.)");

    options.AddOptionNoArg
      (SEP_POS_OPT, "SEP_POS_OPT", REGULAR_OPTG,
       "-sep_pos", 
       "Use isosurface lookup table separating positive grid vertices.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (SEP_OPP_OPT, "SEP_OPP_OPT", REGULAR_OPTG,
       "-sep_opp", 
       "Use isosurface lookup table separating opposite cube vertices.",
       "(Default.)");

    options.AddOptionNoArg
      (NO_SEP_OPP_OPT, "NO_SEP_OPP_OPT", REGULAR_OPTG,
       "-no_sep_opp", 
       "Do not require separation of opposite cube vertices.");
    options.AddToHelpMessage
      (NO_SEP_OPP_OPT, "Opposite cube vertices may not be separated",
       "in isosurface lookup table.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (NEG_ORIENT_OPT, "NEG_ORIENT_OPT", REGULAR_OPTG,
       "-neg_orient",
       "Use isosurface lookup table with negative orientation. (Default.)");
    options.AddToHelpMessage
      (NEG_ORIENT_OPT,
       "Isosurface simplices/polygons/polytopes are oriented",
       "so that normal points to negative region.");

    options.AddOptionNoArg
      (POS_ORIENT_OPT, "POS_ORIENT_OPT", REGULAR_OPTG,
       "-pos_orient",
       "Use isosurface lookup table with positive orientation.");
    options.AddToHelpMessage
      (POS_ORIENT_OPT,
       "Isosurface simplices/polygons/polytopes are oriented",
       "so that normal points to positive region.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);
    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOptionNoArg
      (OFF_OPT, "OFF_OPT", REGULAR_OPTG,
       "-off", "Output in Geomview .off format. (Default.)");

    options.AddOptionNoArg
      (PLY_OPT, "PLY_OPT", REGULAR_OPTG,
       "-ply", "Output in Standford Polygon File Format (PLY).");
    options.AddToHelpMessage
      (PLY_OPT, "(Allowed only with 3D scalar data.)");

    options.AddOptionNoArg
      (FIG_OPT, "FIG_OPT", REGULAR_OPTG,
       "-fig", "Output in FIG (xfig) format.");
    options.AddToHelpMessage
      (FIG_OPT, "(Allowed only with 2D scalar data.)");

    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOption1Arg
      (ISOTABLE_DIR_OPT, "ISOTABLE_DIR_OPT", REGULAR_OPTG,
       "-isotable_dir", "{isotable_directory}",
       "Specify directory containing appropriate isosurface lookup table(s).");

    options.AddOption1Arg
      (ISOTABLE_PATH_OPT, "ISOTABLE_PATH_OPT", REGULAR_OPTG,
       "-isotable", "{isotable_path}",
       "Specify isosurface isosurface lookup table.",
       "(Does not work with option -highres.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);

    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOption1Arg
      (OUTPUT_FILENAME_OPT, "OUTPUT_FILENAME_OPT", REGULAR_OPTG,
       "-o", "{output_filename}",
       "Write isosurface to file {output_filename}.");

    options.AddOption1Arg
      (OUTPUT_FILENAME_PREFIX_OPT, "OUTPUT_FILENAME_PREFIX_OPT", REGULAR_OPTG,
       "-prefix", "{prefix}",
       "Use {prefix} as filename prefix in constructing output file name.");

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
      (SILENT_OPT, "SILENT_OPT", REGULAR_OPTG,
       "-silent", "Silent mode. Do not print any output.");
    options.AddSynonym(SILENT_OPT, "-s");

    options.AddOptionNoArg
      (VERBOSE_OPT, "VERBOSE_OPT", REGULAR_OPTG,
       "-verbose", "Verbose mode. Print additional information.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (NO_WRITE_OPT, "NO_WRITE_OPT", REGULAR_OPTG,
       "-no_write", "Don't write isosurface. (Used mainly in testing.)");

    options.AddOptionNoArg
      (TIME_OPT, "TIME_OPT", REGULAR_OPTG,
       "-time", "Output running time of various operations.");

    options.AddUsageOptionNewline(REGULAR_OPTG);
    options.AddUsageOptionBeginOr(REGULAR_OPTG);

    options.AddOptionNoArg
      (USAGE_OPT, "USAGE_OPT", REGULAR_OPTG, "-usage",
       "Print usage message.");

    options.AddOptionNoArg
       (MORE_OPTIONS_OPT, "MORE_OPTIONS_OPT", REGULAR_OPTG, "-more_options",
        "Print additional (more) options.");

    options.AddOptionNoArg
      (ALL_OPTIONS_OPT, "ALL_OPTIONS_OPT", REGULAR_OPTG, "-all_options",
       "Print usage message with all options.");

    options.AddUsageOptionNewline(REGULAR_OPTG);

    options.AddOptionNoArg
      (HELP_OPT, "HELP_OPT", REGULAR_OPTG,
       "-help", "Print this help message.");
    options.AddSynonym(HELP_OPT, "-h");

    options.AddOptionNoArg
        (HELP_MORE_OPT, "HELP_MORE_OPT", REGULAR_OPTG, "-help_more",
         "Print help message for additional (more) options.");

     options.AddOptionNoArg
       (HELP_ALL_OPT, "HELP_ALL_OPT", REGULAR_OPTG, "-help_all",
        "Print help message for all options.");

    options.AddUsageOptionEndOr(REGULAR_OPTG);


    // More options.

    options.AddUsageOptionBeginOr(EXTENDED_OPTG);

    options.AddOptionNoArg
      (EDGE1_OPT, "EDGE1_OPT", EXTENDED_OPTG,
       "-edge1", "Represent grid edges by a single identifier.");

    options.AddOptionNoArg
      (EDGE2_OPT, "EDGE2_OPT", EXTENDED_OPTG,
       "-edge2", "Represent grid/mesh edges as a pair of endpoints.");

    options.AddUsageOptionEndOr(EXTENDED_OPTG);

    options.AddOptionNoArg
      (NO_COMMENTS_OPT, "NO_COMMENTS_OPT", EXTENDED_OPTG,
       "-no_comments",
       "Do not put any comments in isosurface output file.");
    options.AddToHelpMessage
      (NO_COMMENTS_OPT,
       "(Mainly used for testing/comparing outputs",
       "without comparing differences in comments.)");
  }


  // Process option optA.
  // - Return true if option OPTA is processed (i.e., matches some case in
  //   switch statement.)
  template <typename IO_INFO_TYPE>
  bool process_option
  (const COMMAND_LINE_OPTION_TYPE optA, const int argc, char **argv,
   int & iarg, IO_INFO_TYPE & io_info)
  {
    switch(optA) {

    case REGION_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.region_length);
      io_info.use_minmax = true;
      break;

    case OCTREE_OPT:
      io_info.use_octree = true;
      break;

    case USE_CUBE_LIST_OPT:
      io_info.use_cube_list = true;
      break;

    case NEP_OPT:
      io_info.use_nep = true;
      io_info.isotable_properties.grid_vertex_label_type =
        NEG_EQUALS_POS;
      io_info.isotable_properties.encoding = BASE3;
      break;

    case NEP_DUP_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.nep_num_dup);
      break;

    case IVOL_OPT:
      io_info.flag_interval_volume = true;
      io_info.isotable_properties.grid_vertex_label_type =
        NEG_STAR_POS;
      io_info.isotable_properties.lookup_table_type =
        INTERVAL_VOLUME;
      io_info.isotable_properties.encoding = BASE3;
      break;

    case SNAP_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%f", &io_info.snap_value);
      io_info.flag_snap = true;
      break;

    case TOPOLOGY_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.isosurface_topology = get_topology(argv[iarg]);
      break;

    case INTERPOLATE_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.interpolation_type = get_interpolation_type(argv[iarg]);
      break;

    case SUBSAMPLE_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.subsample_resolution);
      io_info.flag_subsample = true;
      break;

    case SUPERSAMPLE_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.supersample_resolution);
      io_info.flag_supersample = true;
      break;

    case HIGHRES_OPT:
      get_high_resolution_arg(iarg, argc, argv, io_info);
      io_info.use_multires = true;
      iarg++;
      break;

    case NUMRES_LEVELS_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &io_info.num_resolution_levels);
      break;

    case EDGE_GROUPS_OPT:
      io_info.isotable_properties.isosurface_triangulation_type = 
        IJKMCUBE_TABLE::EDGE_GROUPS;
      break;

    case CHULL_OPT:
      io_info.isotable_properties.isosurface_triangulation_type = 
        IJKMCUBE_TABLE::CONVEX_HULL;
      break;

    case SEP_NEG_OPT:
      io_info.isotable_properties.isosurface_separation_type = 
        SEPARATE_NEG;
      break;

    case SEP_POS_OPT:
      io_info.isotable_properties.isosurface_separation_type = 
        SEPARATE_POS;
      break;

    case SEP_OPP_OPT:
      io_info.isotable_properties.separate_opposite =
        TRUE_SEPARATE_OPPOSITE;
      break;

    case NO_SEP_OPP_OPT:
      io_info.isotable_properties.separate_opposite =
        FALSE_SEPARATE_OPPOSITE;
      break;

    case NEG_ORIENT_OPT:
      io_info.isotable_properties.iso_poly_orientation =
        IJKMCUBE_TABLE::NEGATIVE_ORIENT;
      break;

    case POS_ORIENT_OPT:
      io_info.isotable_properties.iso_poly_orientation =
        IJKMCUBE_TABLE::POSITIVE_ORIENT;
      break;

    case EDGE1_OPT:
      io_info.edge_representation = EDGE_ID;
      break;

    case EDGE2_OPT:
      io_info.edge_representation = EDGE_ENDPOINT_PAIR;
      break;

    case LOG_INTERVAL_OPT:
      int log2_interval;
      iarg++;
      if (iarg >= argc) usage_error();
      sscanf(argv[iarg], "%d", &log2_interval);
      io_info.merge_edges_parameters.SetLog2Interval(log2_interval);
      break;

    case ISOTABLE_DIR_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.isotable_directory = argv[iarg];
      break;

    case ISOTABLE_PATH_OPT:
      iarg++;
      if (iarg >= argc) usage_error();
      io_info.isotable_path = argv[iarg];
      break;      

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

    case VERBOSE_OPT:
      io_info.flag_verbose = true;
      break;

    case TIME_OPT:
      io_info.flag_report_time = true;
      break;

    case USAGE_OPT:
      usage(cout, 0);
      break;

    case MORE_OPTIONS_OPT:
      print_options(cout, EXTENDED_OPTG, 0);
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

    case HELP_ALL_OPT:
      help_all();
      break;

    case NO_COMMENTS_OPT:
      io_info.flag_no_comments = true;
      break;

    default:
      return false;
    };

    return true;
  }


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
      exit(-1);
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
    if (!io_info.is_file_format_set) {
      if (io_info.output_filename != "") {
        const IO_INFO::OUTPUT_FORMAT output_format =
          io_info.GetFileType(io_info.output_filename, OUTPUT_DEFINITIONS::OFF);

        io_info.SetOutputFormat(output_format);
        io_info.SetOutputFilename(output_format, io_info.output_filename);
      }
      else {
        io_info.flag_output_off = true;
      }
    }
    else {
      if (io_info.output_filename != "") {
        if (io_info.NumOutputFormats() > 1) {
          cerr << "Usage error. Cannot set output filename when more than one"
               << endl
               << "  output format is selected." << endl;
          exit(-1);
        }

        io_info.SetOutputFilename(io_info.output_filename);
      }
    }

    if (io_info.use_octree && io_info.use_minmax) {
      cerr << "Usage error.  Can't use both -region and -octree parameters.";
      usage_error();
    }

    if (io_info.flag_interval_volume && io_info.isovalue.size() < 2) {
      cerr << "Usage error.  Need at least two isovalues for interval volume generation." << endl;
      exit(-1);
    }

    if (io_info.flag_subsample && io_info.subsample_resolution <= 1) {
      cerr << "Usage error.  Subsample resolution must be an integer greater than 1."
           << endl;
      exit(-1);
    };

    if (io_info.IsOutputFilenameSet() && io_info.flag_use_stdout) {
      cerr << "Usage error.  Can't use both -o and -stdout parameters."
           << endl;
      exit(-1);
    };

    if (io_info.flag_snap &&
        (io_info.snap_value < 0.0 || io_info.snap_value > 0.5)) {
      cerr << "Usage error.  Illegal snap value " << io_info.snap_value << "."
           << endl;
      cerr << "        Snap value must be in range [0.0, 0.5]." << endl;
      exit(-1);
    };

    if (io_info.flag_subsample && io_info.use_multires) {
      cerr << "Usage error.  Can't use both -subsample and -highres parameters."
           << endl;
      exit(-1);
    }

    if (io_info.flag_subsample && io_info.flag_supersample) {
      cerr << "Usage error.  Can't use both -subsample and -supersample parameters."
           << endl;
      exit(-1);
    }

    if (io_info.flag_subsample && io_info.subsample_resolution <= 1) {
      cerr << "Usage error.  Subsample resolution must be an integer greater than 1."
           << endl;
      exit(-1);
    };

    if (io_info.IsOutputFilenameSet() && io_info.flag_use_stdout) {
      cerr << "Usage error.  Can't use both -o and -stdout parameters."
           << endl;
      exit(-1);
    };

    if (io_info.IsOutputFilenameSet()) {

      if (io_info.flag_interval_volume) {
        if (io_info.isovalue.size() > 2) {
          cerr << "Usage error.  Cannot specify -ivol output file when input contains"
                << endl;
          cerr << "  more than two isovalues." << endl;
          exit(-1);
        }
      }
      else {
        if (io_info.isovalue.size() > 1) {
          cerr << "Usage error.  Cannot specify isosurface output file when input contains"
              << endl;
          cerr << "  more than one isovalue." << endl;
          exit(-1);
        }
      }
    }

    if (io_info.flag_snap) {
      io_info.isotable_properties.grid_vertex_label_type =
        NEG_EQUALS_POS;
      io_info.isotable_properties.encoding = IJKMCUBE_TABLE::BASE3;
    }

    if (io_info.isotable_properties.IsGridVertexLabelTypeUndefined()) {
      io_info.isotable_properties.grid_vertex_label_type = NEG_POS;
      io_info.isotable_properties.encoding = IJKMCUBE_TABLE::BINARY;
    }

    if (io_info.isotable_properties.IsTableTypeUndefined()) {
      io_info.isotable_properties.lookup_table_type = ISOSURFACE;
    }

    check_io_info(io_info);
  }

}


// parse command line
// control parameters, followed by one or more isovalues,
// followed by input file name
void IJKMCUBE::parse_command_line(int argc, char **argv, IO_INFO & io_info)
{
  IJK::ERROR error;

  set_command_line_options(io_info);

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
}


// Check input information/flags.
void IJKMCUBE::check_input
(const IO_INFO & io_info, 
 const MC_SCALAR_GRID_BASE & scalar_grid)
{
  if (io_info.isotable_directory == "") {
    cerr << "Usage error.  Unknown isotable directory." << endl;
    cerr << "  Use -dir {isotable_directory} argument or set environment variable IJK_ISOTABLE_DIR." << endl;
    exit(-1);
  }

  if (io_info.flag_interval_volume) {
    // Construct interval volume
    if (io_info.isovalue.size() > 2 && io_info.flag_use_stdout) {
      cerr << "Usage error.  Cannot use stdout for more than one interval volume." << endl;
      exit(-1);
    }

    if (io_info.isovalue.size() > 2 && io_info.IsOutputFilenameSet()) {
      cerr << "Usage error.  Cannot specify output file for more than one interval volume." << endl;
      exit(-1);
    }
  }
  else {
    // Construct isosurface
    if (io_info.isovalue.size() > 1 && io_info.flag_use_stdout) {
      cerr << "Usage error.  Cannot use stdout for more than one isovalue." << endl;
      exit(-1);
    }

    if (io_info.isovalue.size() > 1 && io_info.IsOutputFilenameSet()) {
      cerr << "Usage error.  Cannot specify output file for more than one isovalue." << endl;
      exit(-1);
    }
  }

  check_highres(scalar_grid.Dimension(), io_info);
}


// **************************************************
// PATH_DELIMITER
// **************************************************

namespace {

#ifdef _WIN32
  const char PATH_DELIMITER = '\\';
#else
  const char PATH_DELIMITER = '/';
#endif
}

// **************************************************
// READ ISOSURFACE LOOKUP TABLE(S)
// **************************************************

/// Read cube isosurface lookup table.
void IJKMCUBE::read_cube_isotable
(const int dimension, const IO_INFO & io_info, 
 ISOSURFACE_TABLE & cube_isotable, IO_TIME & io_time)
{
  read_isosurface_table
    (dimension, "cube", io_info, cube_isotable, io_time);
}


/// Read polytope isosurface lookup table.
void IJKMCUBE::read_poly_isotable
(const IO_INFO & io_info, MC_DATA & mc_data, IO_TIME & io_time)
{
  const int dimension = mc_data.ScalarGrid().Dimension();
  PROCEDURE_ERROR error("read_poly_isotable");

  if (mc_data.UseNEP() || mc_data.Snap()) {
    read_cube_isotable
      (dimension, io_info, mc_data.isotable.cube_nep, io_time);
    set_in_facet_iso_patches(mc_data.isotable.cube_nep);

    if (!check_dimension(mc_data.isotable.cube_nep, mc_data.ScalarGrid(), 
                         error))
      { throw(error); };
  }
  else {
    read_cube_isotable
      (dimension, io_info, mc_data.isotable.cube, io_time);

    if (!check_dimension(mc_data.isotable.cube, mc_data.ScalarGrid(), error))
      { throw(error); };
  }

  ISOSURFACE_TOPOLOGY isosurface_topology = mc_data.IsosurfaceTopology();
  bool flag_multires = mc_data.UseMultires();

  if (isosurface_topology != ISOTABLE_TOPOLOGY) {
    IJK::ERROR warning;

    if (!check_cube_isotable_fits_topology
        (mc_data.isotable.cube, isosurface_topology, warning)) {
      cerr << endl;
      cerr << "*** Warning: ";
      warning.Print(cerr);
      cerr << "  Continuing..." << endl;
      cerr << endl;
    }
  }

  if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY ||
      isosurface_topology == LINEAR_TOPOLOGY ||
      flag_multires) {

    read_isosurface_table
      (dimension, "pyramid", io_info, mc_data.isotable.pyramid, io_time);

    if (!check_dimension(mc_data.isotable.pyramid, 
                         mc_data.ScalarGrid(), error))
      { throw(error); };
  }

  if (flag_multires) {
    read_isosurface_table
      (dimension, "simplex", io_info, mc_data.isotable.simplex, io_time);

    if (!check_dimension(mc_data.isotable.simplex, 
                         mc_data.ScalarGrid(), error))
      { throw(error); };
  }

  if (isosurface_topology != ISOTABLE_TOPOLOGY) 
    { mc_data.isotable.ComputeAmbiguityInformation(); }

  if (io_info.flag_verbose && !io_info.flag_use_stdout) {
    // Print newline after printing read isotable information.
    cout << endl;
  }
}


namespace {
  
  class ISOTABLE_FILENAME_STR {
  public:
    std::string prefix;
    std::string middle;
    std::string suffix;
    std::string pathname_prefix;    
  };

  void open_file_local
  (const ISOTABLE_FILENAME_STR & filename_str,
   std::string & filename, 
   std::string & pathname,
   std::ifstream & isotable_file)
  {
    filename = filename_str.prefix + filename_str.middle + 
      filename_str.suffix;
    pathname = filename_str.pathname_prefix + filename;

    isotable_file.open(pathname, ios::in);
  }

}


void IJKMCUBE::open_isotable_file
(const int dimension, const char * poly_name,
 const IO_INFO & io_info, const bool flag_use_local_directory,
 std::string & isotable_filename,
 std::string & isotable_pathname,
 std::ifstream & isotable_file,
 std::vector<std::string> & path_list)
{
  using namespace IJKMCUBE_TABLE;
  const ISOTABLE_TYPE isotable_type =
    get_isotable_type(io_info.isotable_properties);
  const std::string isotable_directory = io_info.isotable_directory;
  ISOTABLE_FILENAME_STR filename_str;

  const std::string filename_prefix =
    get_isotable_filename_prefix(dimension, poly_name, isotable_type);
  const std::string filename_suffix =
    get_isotable_filename_suffix(dimension, poly_name, isotable_type);
  std::string isotable_pathname_prefix;

  if (flag_use_local_directory)
    { isotable_pathname_prefix = ""; }
  else
    { isotable_pathname_prefix = isotable_directory + PATH_DELIMITER; }
    
  filename_str.prefix = filename_prefix;
  filename_str.suffix = filename_suffix;
  filename_str.pathname_prefix = isotable_pathname_prefix;

  const std::string triangulation_label =
    io_info.isotable_properties.TriangulationTypeLabel();
  const std::string separation_label =
    io_info.isotable_properties.SeparationTypeLabel();
  const std::string orientation_label =
    io_info.isotable_properties.IsoPolyOrientationLabel();

  if (triangulation_label != "") {
    if (separation_label != "") {
      if (orientation_label != "") {
        filename_str.middle =
          "." + triangulation_label + "." + separation_label +
          "." + orientation_label + ".";

        open_file_local
          (filename_str, isotable_filename, isotable_pathname, isotable_file);
        path_list.push_back(isotable_pathname);
        if (isotable_file) { return; }
      }

      filename_str.middle =
        "." + triangulation_label + "." + separation_label + ".";

      open_file_local
        (filename_str, isotable_filename, isotable_pathname, isotable_file);
      path_list.push_back(isotable_pathname);
      if (isotable_file) { return; }
    }
    else if (orientation_label != "") {
        filename_str.middle =
          "." + triangulation_label + "." + orientation_label + ".";

        open_file_local
          (filename_str, isotable_filename, isotable_pathname, isotable_file);
        path_list.push_back(isotable_pathname);
        if (isotable_file) { return; }
    }

    filename_str.middle =
      "." + triangulation_label + ".";

    open_file_local
      (filename_str, isotable_filename, isotable_pathname, isotable_file);
    path_list.push_back(isotable_pathname);
    if (isotable_file) { return; }
  }
  else if (separation_label != "") {
    if (orientation_label != "") {
      filename_str.middle =
        "." + separation_label + "." + orientation_label + ".";

      open_file_local
        (filename_str, isotable_filename, isotable_pathname, isotable_file);
      path_list.push_back(isotable_pathname);
      if (isotable_file) { return; }
    }

    filename_str.middle = "." + separation_label + ".";

    open_file_local
      (filename_str, isotable_filename, isotable_pathname, isotable_file);
    path_list.push_back(isotable_pathname);
    if (isotable_file) { return; }
  }
  else if (orientation_label != "") {
    filename_str.middle = "." + orientation_label + ".";
    
    open_file_local
      (filename_str, isotable_filename, isotable_pathname, isotable_file);
    path_list.push_back(isotable_pathname);
    if (isotable_file) { return; }
  }


  filename_str.middle = ".";
  open_file_local
    (filename_str, isotable_filename, isotable_pathname, isotable_file);
  path_list.push_back(isotable_pathname);
}


void IJKMCUBE::read_isosurface_table
(const int dimension, const char * poly_name,
 const IO_INFO & io_info, ISOSURFACE_TABLE & isotable)
{
  std::string isotable_filename, isotable_pathname, isotable_pathnameB;
  std::vector<std::string> path_listA;
  std::vector<std::string> path_listB;

  ifstream isotable_file;

  if (io_info.isotable_path == "") {
    open_isotable_file
      (dimension, poly_name, io_info, false,
       isotable_filename, isotable_pathname, isotable_file,
       path_listA);

    if (!isotable_file) {
      isotable_file.clear();
      open_isotable_file
        (dimension, poly_name, io_info, true,
         isotable_filename, isotable_pathnameB, isotable_file,
         path_listB);
    };

    if (!isotable_file) {
      cerr << "Unable to find/open isosurface table file. Tried: "
           << endl;
      for (int i = 0; i < path_listA.size(); i++)
        {  cerr << "  " << path_listA[i] << endl; }
      for (int i = 0; i < path_listB.size(); i++)
        {  cerr << "  " << path_listB[i] << endl; }        
      exit(30);
    };
  }
  else {
    isotable_pathname = io_info.isotable_path;
    isotable_file.open(io_info.isotable_path, ios::in);

    if (!isotable_file) {
      cerr << "Unable to open isosurface table "
           << io_info.isotable_path << "." << endl;
      exit(30);
    }
  }

  try {
    if (io_info.flag_verbose && !io_info.flag_use_stdout) {
      
      if (io_info.isotable_properties.TableType() ==
          INTERVAL_VOLUME)
        { cout << "Reading interval volume table in file: "; }
      else
        { cout << "Reading isosurface table in file: "; }
      cout << isotable_pathname << endl;
    }

    IJKXIO::read_xit(isotable_file, isotable);
  }
  catch(...) {
    cerr << "Error reading file: " << isotable_pathname << "." << endl;
    throw;
  };

  isotable_file.close();

  // Check for incorrect isosurface lookup table.
  ERROR error;

  try {
    const LOOKUP_TABLE_TYPE table_type =
      io_info.isotable_properties.TableType();
    if (!isotable.Properties().CheckTableType(table_type, error)) 
      { throw error; }

    const ENCODING encoding = io_info.isotable_properties.Encoding();
    if (!isotable.Properties().CheckEncoding(encoding, error)) 
      { throw error; }
  }
  catch (ERROR & error) {
    cerr << "Usage error. Incorrect isosurface lookup table: "
         << isotable_pathname << endl;
    error.Print(cerr);
    cerr << "Exiting." << endl;
    exit(-1);
  }

  // Check for inconsistencies in isosurface lookup table.
  if (!isotable.Check(error)) {
    cerr << "Warning: Data structure inconsistency in isosurface table "
         << isotable_pathname << "." << endl;
    error.Print(cerr);
    cerr << "  Attempting to continue..." << endl << endl;
  }

  // Check for inconsistencies in isosurface lookup table.
  ERROR error2;
  if (!isotable.Properties().Check(io_info.isotable_properties, error2)) {
    cerr << "Usage error. Incorrect isosurface lookup table: "
         << isotable_pathname << endl;
    error2.Print(cerr);
    cerr << "Exiting." << endl;
    exit(-1);
  }

}


void IJKMCUBE::read_isosurface_table
(const int dimension, const char * poly_name, const IO_INFO & io_info, 
 ISOSURFACE_TABLE & isotable, IO_TIME & io_time)
{
  ELAPSED_TIME wall_time;

  wall_time.getElapsed();

  read_isosurface_table
    (dimension, poly_name, io_info, isotable);

  io_time.read_table_time = wall_time.getElapsed();
}


ISOTABLE_TYPE IJKMCUBE::get_isotable_type
(const ISOSURFACE_TABLE_PROPERTIES & properties)
{
  ISOTABLE_TYPE isotable_type;
  if (properties.GridVertexLabelType() == NEG_EQUALS_POS)
    { isotable_type = NEP; }
  else if (properties.TableType() == INTERVAL_VOLUME)
    { isotable_type = IVOL; }
  else
    { isotable_type = BINARY; }

  return isotable_type;
}


// Get default isotable directory and get directory from environment.
void IJKMCUBE::get_isotable_directory(std::string & isotable_directory)
{
#ifdef IJK_ISOTABLE_DIR
  isotable_directory = std::string(IJK_ISOTABLE_DIR);
#endif

  const char * envir_isotable_dir = getenv("IJK_ISOTABLE_DIR");
  if (envir_isotable_dir != NULL) {
    isotable_directory = std::string(envir_isotable_dir);
  };
}

// **************************************************
// OUTPUT ISOSURFACE
// **************************************************

void IJKMCUBE::output_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const MCUBE_INFO & mcube_info)
{
  output_isosurface(output_info, mc_data, 
                    mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
                    mcube_info);
}


void IJKMCUBE::output_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const MCUBE_INFO & mcube_info)
{
  if (!output_info.flag_use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, mc_data, vertex_coord, slist, mcube_info);
  }

  if (!output_info.flag_nowrite) {
    std::vector<std::string> comment;

    if (!output_info.flag_no_comments) {

      if (output_info.flag_interval_volume)
        { comment.push_back("Marching Cubes interval volume."); }
      else
        { comment.push_back("Marching Cubes isosurface."); }

      add_meshfile_comments(output_info, comment);
    }

    IJK::write_mesh(output_info, vertex_coord, slist, comment);
  }
}


void IJKMCUBE::output_nep_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const MCUBE_INFO & mcube_info)
{
  output_nep_isosurface
    (output_info, mc_data, 
     mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
     mcube_info);
}


void IJKMCUBE::output_nep_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const MCUBE_INFO & mcube_info)
{
  if (!output_info.flag_use_stdout && !output_info.flag_silent) {
    report_nep_info(output_info, mc_data, vertex_coord, slist, mcube_info);
  }

  if (!output_info.flag_nowrite) {
    std::vector<std::string> comment;

    if (!output_info.flag_no_comments) {
      comment.push_back("Marching Cubes (NEP) isosurface.");
      comment.push_back
        ("- Constructed with NEP (negative, equals, positive) isosurface lookup table.");

      add_meshfile_comments(output_info, comment);
    };

    IJK::write_mesh(output_info, vertex_coord, slist, comment);

  }
}


void IJKMCUBE::output_snapmc_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const SNAP_INFO & snap_info)
{
  output_snapmc_isosurface
    (output_info, mc_data, 
     mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
     snap_info);
}


void IJKMCUBE::output_snapmc_isosurface
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const SNAP_INFO & snap_info)
{
  if (!output_info.flag_use_stdout && !output_info.flag_silent) {
    report_snap_info(output_info, mc_data, vertex_coord, slist, snap_info);
  }

  if (!output_info.flag_nowrite) {
    std::vector<std::string> comment;

    if (!output_info.flag_no_comments) {
      comment.push_back("Marching Cubes SnapMC isosurface.");

      add_meshfile_comments(output_info, comment);
    };

    write_mesh(output_info, vertex_coord, slist, comment);
  }
}


void IJKMCUBE::output_isosurface_color
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 const MCUBE_INFO & mcube_info)
{
  output_isosurface_color
    (output_info, mc_data, 
     mc_isosurface.vertex_coord, mc_isosurface.simplex_vert,
     front_color, back_color, mcube_info);
}


void IJKMCUBE::output_isosurface_color
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
 const MCUBE_INFO & mcube_info)
{
  if (!output_info.flag_use_stdout && !output_info.flag_silent) {
    report_iso_info(output_info, mc_data, vertex_coord, slist, mcube_info);
  }
  
  if (!output_info.flag_nowrite) {
    write_mcube_mesh_color
      (output_info, vertex_coord, slist, front_color, back_color);
  }
}


void IJKMCUBE::output_isosurface_color_alternating
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const MC_ISOSURFACE & mc_isosurface,
 const MCUBE_INFO & mcube_info)
{
  const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
  const VERTEX_INDEX nums = 
    mc_isosurface.simplex_vert.size()/numv_per_simplex;

  IJK::ARRAY<COLOR_TYPE> front_color(4*nums);
  IJK::ARRAY<COLOR_TYPE> back_color(4*nums);
  set_color_alternating
    (mc_data.ScalarGrid(), mc_isosurface.cube_containing_simplex,
     front_color.Ptr());
  set_color_alternating
    (mc_data.ScalarGrid(), mc_isosurface.cube_containing_simplex,
     back_color.Ptr());

  output_isosurface_color
    (output_info, mc_data, mc_isosurface.vertex_coord,
     mc_isosurface.simplex_vert, 
     front_color.PtrConst(), back_color.PtrConst(),
     mcube_info);
}


// **************************************************
// WRITE_MESH
// **************************************************

void IJKMCUBE::write_mcube_mesh
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist)
{
  IJK::write_mesh(output_info, vertex_coord, slist);
}


void IJKMCUBE::write_mcube_mesh_color
(const OUTPUT_INFO & output_info,
 const vector<COORD_TYPE> & vertex_coord, const vector<VERTEX_INDEX> & slist,
 const COLOR_TYPE * front_color, const COLOR_TYPE * back_color)
{
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.NumVerticesPerSimplex();
  const bool flag_use_stdout = output_info.flag_use_stdout;

  ofstream output_file;
  ERROR error_mcube("write_mcube_mesh_color");

  std::string ofilename;

  if (output_info.flag_output_off) {
    if (!flag_use_stdout) {
      ofilename = output_info.output_off_filename;
      output_file.open(ofilename.c_str(), ios::out);
      ijkoutColorFacesOFF(output_file, dimension, numv_per_simplex,
                          vertex_coord, slist, front_color, back_color);
      output_file.close();
    }
    else {
      ijkoutColorFacesOFF(std::cout, dimension, numv_per_simplex,
                          vertex_coord, slist, front_color, back_color);
    };
  }
  else {
    throw error_mcube
      ("Writing color mesh currently only implemented for .off output format.");
  }

  if (!flag_use_stdout && !output_info.flag_silent)
    cout << "Wrote output to file: " << ofilename << endl;
}


// **************************************************
// RESCALE ROUTINES
// **************************************************

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

  bool unit_spacing(const std::vector<COORD_TYPE> & spacing)
    // return true if spacing not defined or spacing along all axes equals 1.0
  {
    for (unsigned int d = 0; d < spacing.size(); d++) {
      if (!AIR_EXISTS(spacing[d])) { return(true); }
      else if (spacing[d] != 1.0) { return(false); };
    }

    return(true);
  }

  void rescale_coord(const std::vector<COORD_TYPE> & grid_spacing,
                     std::vector<COORD_TYPE> & vertex_coord)
  {
    const int dimension = grid_spacing.size();

    if (unit_spacing(grid_spacing)) { return; }

    const VERTEX_INDEX numv = vertex_coord.size()/dimension;
    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dimension; d++) {
        vertex_coord[iv*dimension+d] *= grid_spacing[d];
      }
    };
  }

}


// Rescale subsampled/supersampled vertex coordinates.
// - Also rescale to reflect grid spacing.
void IJKMCUBE::rescale_vertex_coord
(const OUTPUT_INFO & output_info, vector<COORD_TYPE> & vertex_coord)
{
  const int subsample_resolution = output_info.subsample_resolution;
  const int supersample_resolution = output_info.supersample_resolution;
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (subsample_resolution <= 0) {
    error.AddMessage("Illegal subsample resolution ", subsample_resolution, ".");
    error.AddMessage("  Subsample resolution must be a positive integer");
  }

  if (supersample_resolution <= 0) {
    error.AddMessage("Illegal supersample resolution ", supersample_resolution, ".");
    error.AddMessage("  Supersample resolution must be a positive integer");
  }

  if (output_info.dimension != output_info.grid_spacing.size()) {
    error.AddMessage("Size of grid spacing array does not equal volume dimension.");
    error.AddMessage("  Grid spacing array has ", 
                     output_info.grid_spacing.size(), " elements.");
    error.AddMessage("  Volume dimension = ", output_info.dimension, ".");
  }

  if (output_info.subsample_resolution != 1)
    { grow_coord(output_info.subsample_resolution, vertex_coord); };

  if (output_info.supersample_resolution != 1)
    { shrink_coord(output_info.supersample_resolution, vertex_coord); };

  rescale_coord(output_info.grid_spacing, vertex_coord);
}


// Rescale subsampled/supersampled vertex coordinates.
// - Also rescale to reflect grid spacing.
void IJKMCUBE::rescale_vertex_coord
(const int subsample_resolution, const int supersample_resolution,
 const COORD_ARRAY & grid_spacing, COORD_ARRAY & vertex_coord)
{
  PROCEDURE_ERROR error("rescale_vertex_coord");

  if (subsample_resolution <= 0) {
    error.AddMessage("Illegal grow factor ", subsample_resolution, ".");
    error.AddMessage("  Grow factor must be a positive integer");
  }

  if (supersample_resolution <= 0) {
    error.AddMessage("Illegal shrink factor ", supersample_resolution, ".");
    error.AddMessage("  Shrink factor must be a positive integer");
  }

  if (vertex_coord.size() == 0) { return; };

  if (grid_spacing.size() < 1) {
    error.AddMessage("Illegal size ", grid_spacing.size(), 
                     " of array grid spacing.");
    error.AddMessage("Size must equal vertex dimension.");
    throw error;
  }

  if (subsample_resolution != 1)
    { grow_coord(subsample_resolution, vertex_coord); };

  if (supersample_resolution != 1)
    { shrink_coord(supersample_resolution, vertex_coord); };

  rescale_coord(grid_spacing, vertex_coord);
}


// **************************************************
// CREATE MESH FILE COMMENTS
// **************************************************

void IJKMCUBE::add_meshfile_comments
(const OUTPUT_INFO & output_info,
 std::vector<std::string> & comment)
{
  std::string str;

  if (output_info.flag_no_comments) { return; }

  output_info.AddIsovalueComment(comment);
  output_info.AddSubsampleResolutionComment(comment);
  output_info.AddSupersampleResolutionComment(comment);

  topology2string(output_info.isosurface_topology, str);
  str = "IsosurfaceTopology: " + str;
  comment.push_back(str);

  if (output_info.use_nep)
    { comment.push_back("IsosurfaceTableType: NegativeEqualsPositive"); }

  if (output_info.Snap()) {
    IJK::create_string
      ("SnapValue: ", output_info.SnapValue(), str);
    comment.push_back(str);
  }

  output_info.AddHighResolutionComments(comment);
};


// **************************************************
// CREATE STRING FROM enum TYPE
// **************************************************

void IJKMCUBE::topology2string
(const ISOSURFACE_TOPOLOGY & isosurface_topology,
 std::string & topology_string)
{
  switch(isosurface_topology) {

    case ISOTABLE_TOPOLOGY:
      topology_string = "Isotable";
      break;

    case CUBE_DECIDER_TOPOLOGY:
      topology_string = "CubeDecider";
      break;

    case ASYMPTOTIC_DECIDER_TOPOLOGY:
      topology_string = "AsymptoticDecider";
      break;

    default:
      topology_string = "Unknown";
      break;
  }
}


// **************************************************
// REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
// **************************************************

void IJKMCUBE::report_num_cubes
(const MC_GRID & full_scalar_grid, const IO_INFO & io_info, 
 const MC_DATA & mc_data)
{
  const int num_grid_cubes = full_scalar_grid.ComputeNumCubes();
  const int num_cubes_in_mc_data = 
    mc_data.ScalarGrid().ComputeNumCubes();

  if (!io_info.flag_use_stdout && !io_info.flag_silent) {

    if (io_info.flag_subsample) {
      // subsampled grid
      cout << num_grid_cubes << " grid cubes.  "
           << num_cubes_in_mc_data << " subsampled grid cubes." << endl;
    }
    else if (io_info.flag_supersample) {
      // supersample grid
      cout << num_grid_cubes << " grid cubes.  "
           << num_cubes_in_mc_data << " supersampled grid cubes." << endl;
    }
    else {
      // use full_scalar_grid
      cout << num_grid_cubes << " grid cubes." << endl;
    }
  }

}


void IJKMCUBE::report_iso_info
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, 
 const vector<VERTEX_INDEX> & slist, 
 const MCUBE_INFO & mcube_info)
{
  const int DIM2(2);
  const int dimension = output_info.dimension;
  const int numv_per_simplex = output_info.NumVerticesPerSimplex();
  const int isosurface_topology = mc_data.IsosurfaceTopology();

  const char * indent4 = "    ";
  string grid_element_name = "cubes";
  if (dimension == DIM2) { grid_element_name = "squares"; };

  VERTEX_INDEX numv = (vertex_coord.size())/dimension;
  VERTEX_INDEX nums = (slist.size())/numv_per_simplex;
  VERTEX_INDEX num_grid_cubes = mcube_info.grid.num_cubes;
  VERTEX_INDEX num_non_empty_cubes = mcube_info.scalar.num_non_empty_cubes;

  float percent = 0.0;
  if (num_grid_cubes > 0)
    { percent = float(num_non_empty_cubes)/float(num_grid_cubes); }
  int ipercent = int(100*percent);
  if (output_info.output_isovalue.size() > 0) {
    cout << "  Isovalue " << output_info.output_isovalue[0] << ".  "
        << numv << " isosurface vertices.  "
        << nums << " isosurface simplices." << endl;
    cout << indent4 << num_non_empty_cubes
        << " (" << ipercent << "%) non-empty " << grid_element_name << "." << endl;
  }

  if (isosurface_topology == ASYMPTOTIC_DECIDER_TOPOLOGY) {
    cout << indent4 << mcube_info.scalar.num_ambiguous_cubes 
         << " ambiguous " << grid_element_name << "." << endl;
    if (dimension > 2) {
      cout << indent4 << mcube_info.scalar.num_non_empty_pyramids << " non-empty pyramids." << endl;
      cout << indent4 << mcube_info.scalar.num_ambiguous_pyramids << " ambiguous_pyramids." << endl;
    }
  }

  /* saddle topology not implemented
     if (mcube_info.scalar.Dimension() == 3 && 
     isosurface_topology == SADDLE_TOPOLOGY) {
     cout << indent4 << mcube_info.scalar.NumCubesWithSaddle(1)
     << " cubes with 1 saddle." << endl;
     cout << indent4 << mcube_info.scalar.NumCubesWithSaddle(2)
     << " cubes with 2 saddles." << endl;
     }
  */

}


void IJKMCUBE::report_nep_info
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, 
 const vector<VERTEX_INDEX> & slist, 
 const MCUBE_INFO & mcube_info)
{
  report_iso_info(output_info, mc_data, vertex_coord, slist, mcube_info);

  if (!output_info.flag_use_stdout && !output_info.flag_silent) {
    cout << "    " << mcube_info.nep.num_non_empty_boundary_facets
         << " non-empty boundary cube facets." << endl;

    if (mc_data.NEPNumDup() != 2) {
      cout << "    " << mcube_info.nep.num_in_facet_cubes
           << " cubes with isosurface patches contained in a facet." << endl;
      cout << "    " << mcube_info.nep.num_dup_iso_patches
           << " duplicate isosurface patches eliminated." << endl;
    }
  }
}


void IJKMCUBE::report_snap_info
(const OUTPUT_INFO & output_info, const MC_DATA & mc_data,
 const vector<COORD_TYPE> & vertex_coord, 
 const vector<VERTEX_INDEX> & slist, 
 const SNAP_INFO & snap_info)
{
  report_nep_info(output_info, mc_data, vertex_coord, slist, snap_info);

  if (!output_info.flag_use_stdout && !output_info.flag_silent) {
    cout << "    " << snap_info.num_snapped_iso_vertices
         << " snapped isosurface vertices." << endl;
  }
}


// **************************************************
// REPORT TIMING INFORMATION
// **************************************************

void IJKMCUBE::report_mcube_time
(const IO_INFO & io_info, const MCUBE_TIME & mcube_time, 
 const char * mesh_type_string)
{
  cout << "CPU time to run Marching Cubes: " 
       << mcube_time.total << " seconds." << endl;
  if (io_info.use_octree) {
    cout << "    Time to create octree: "
         << mcube_time.preprocessing << " seconds." << endl;
  }
  if (io_info.use_minmax) {
    cout << "    Time to create min/max regions: "
         << mcube_time.preprocessing << " seconds." << endl;
  }
  if (io_info.flag_snap) {
    cout << "    Time to snap grid scalar values: "
         << mcube_time.snap << " seconds." << endl;
  }
  if (io_info.use_cube_list) {
    cout << "    Time to create list of non-empty cubes: "
         << mcube_time.preprocessing << " seconds." << endl;
  }
  if (io_info.use_multires) {
    cout << "    Time to preprocess multi-resolution mesh: "
         << mcube_time.process_multires << " seconds." << endl;
  }
  cout << "    Time to extract " << mesh_type_string << " triangles: "
       << mcube_time.extract << " seconds." << endl;
  cout << "    Time to merge identical "
       << mesh_type_string << " vertices: " 
       << mcube_time.merge << " seconds." << endl;
  cout << "    Time to position "
       << mesh_type_string << " vertices: "
       << mcube_time.position << " seconds." << endl;
}


void IJKMCUBE::report_time
(const IO_INFO & io_info, const IO_TIME & io_time, 
 const MCUBE_TIME & mcube_time, const double total_elapsed_time)
{
  const char * ISOSURFACE_STRING = "isosurface";
  const char * INTERVAL_VOLUME_STRING = "interval volume";
  const char * mesh_type_string = NULL;
  
  if (!io_info.flag_interval_volume) {
    mesh_type_string = ISOSURFACE_STRING;
  }
  else {
    mesh_type_string = INTERVAL_VOLUME_STRING;
  };

  cout << "Time to read file " << io_info.input_filename << ": "
       << io_time.read_nrrd_time << " seconds." << endl;

  cout << "Time to read " << mesh_type_string << " lookup tables: "
       << io_time.read_table_time << " seconds." << endl;

  report_mcube_time(io_info, mcube_time, mesh_type_string);
  if (!io_info.flag_nowrite) {
    cout << "Time to write "
         << mesh_type_string << ": " 
         << io_time.write_time << " seconds." << endl;
  };
  cout << "Total elapsed time: " << total_elapsed_time
       << " seconds." << endl;
}


// **************************************************
// USAGE/HELP MESSAGES
// **************************************************

// local namespace
namespace {

  void usage_msg(std::ostream & out)
  {
    out << "Usage: ijkmcube [OPTIONS] {isovalue1 isovalue2 ...} {input filename}" << endl;
  }


  void print_options_title(std::ostream & out,
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
    cout << "ijkmcube - Marching Cubes algorithm and variants."
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


void IJKMCUBE::usage(std::ostream & out, const int return_code)
{
  usage_msg(out);
  options_msg(out, REGULAR_OPTG);
  exit(return_code);
}


void IJKMCUBE::usage_error()
{
  usage(cerr, -1);
}


void IJKMCUBE::usage_all(std::ostream & out, const int return_code)
{
  usage_msg(out);
  options_msg(out, REGULAR_OPTG);
  options_msg(out, EXTENDED_OPTG);
  exit(return_code);
}


void IJKMCUBE::help()
{
  help_msg();
  cout << endl;
  print_help_options(REGULAR_OPTG);
  exit(0);
}


void IJKMCUBE::print_options
  (std::ostream & out, const COMMAND_LINE_OPTION_GROUP group,
   const int return_code)
{
  options_msg(out, group);
  exit(return_code);
}


void IJKMCUBE::help_all()
{
  help_msg();
  cout << endl;

  print_help_options(REGULAR_OPTG);
  cout << endl;
  print_help_options(EXTENDED_OPTG);
  cout << endl;

  exit(0);
}


void IJKMCUBE::print_options_help
  (const COMMAND_LINE_OPTION_GROUP group)
{
  print_help_options(group);
  exit(0);
}


// **************************************************
// CLASS IO_INFO
// **************************************************

void IJKMCUBE::IO_INFO::Init()
{
  flag_color_alternating = false;
}


// **************************************************
// SET ROUTINES
// **************************************************

void IJKMCUBE::set_mc_data
(const IO_INFO & io_info, MC_DATA & mc_data, MCUBE_TIME & mcube_time)
{
  const bool flag_use_stdout = io_info.flag_use_stdout;
  const bool flag_silent = io_info.flag_silent;
  PROCEDURE_ERROR error("set_mc_data");

  if (!mc_data.IsScalarGridSet()) {
    error.AddMessage("Programming error. Scalar field must be set before set_mc_data is called.");
    throw error;
  }
  
  // Set data structures in mc data  
  clock_t t0 = clock();

  if (io_info.use_octree) {
    if (io_info.flag_snap)
      { mc_data.SetSnapOctree(); }
    else
      { mc_data.SetOctree(); }
  }
  else if (io_info.use_minmax) {
    if (io_info.flag_snap)
      { mc_data.SetSnapMinmaxRegions(io_info.region_length); }
    else
      { mc_data.SetMinmaxRegions(io_info.region_length); };
  }

  // Set flags in mc data  
  if (io_info.use_nep) {
    mc_data.SetNEPOn(io_info.nep_num_dup);
  };

  if (io_info.flag_interval_volume)
    { mc_data.SetFlagIntervalVolume(true); };

  if (io_info.flag_snap)
    { mc_data.SetSnapOn(io_info.snap_value, io_info.nep_num_dup); }

  if (io_info.use_multires) {
    mc_data.SetHighResolutionRegions(io_info.high_resolution_region);
    mc_data.SetMultiresGrid(io_info.num_resolution_levels);
  }

  if (io_info.use_octree || io_info.use_minmax || io_info.use_multires) {
    clock_t t1 = clock();
    mcube_time.preprocessing += float(t1-t0)/CLOCKS_PER_SEC;
    mcube_time.total += mcube_time.preprocessing;
  }

  if (io_info.use_cube_list) { mc_data.SetUseCubeList(true); }

  if (io_info.flag_color_alternating) 
    { mc_data.SetCubeContainingSimplexFlag(true); };

  mc_data.SetEdgeRepresentation(io_info.edge_representation);
  mc_data.SetIsosurfaceTopology(io_info.isosurface_topology);
  mc_data.SetInterpolationType(io_info.interpolation_type);

  mc_data.SetMergeEdgesParameters(io_info.merge_edges_parameters);
}

void IJKMCUBE::set_output_info
(const ISOSURFACE_TABLE & isotable, const IO_INFO & io_info, 
 const int i, OUTPUT_INFO & output_info)
{
  output_info.Set(io_info);

  output_info.dimension = isotable.Dimension();
  output_info.num_vertices_per_isopoly = isotable.NumVerticesPerSimplex();

  output_info.output_isovalue.clear();
  output_info.output_isovalue.push_back(io_info.isovalue[i]);

  if (io_info.flag_interval_volume) {
    if (i+1 < int(io_info.isovalue.size()))
      { output_info.output_isovalue.push_back(io_info.isovalue[i+1]); }
    output_info.mesh_dimension = output_info.dimension;
  }
  else {
    output_info.mesh_dimension = output_info.dimension-1;
  }

  if (io_info.IsOutputFilenameSet()) {
    output_info.output_filename = io_info.output_filename;
  }
  else {
    output_info.ConstructOutputFilenames(i, io_info.flag_interval_volume);
  }
}


void IJKMCUBE::set_color_alternating
(const MC_GRID & grid, const vector<VERTEX_INDEX> & cube_list, 
 COLOR_TYPE * color)
{
  const int dimension = grid.Dimension();
  IJK::ARRAY<GRID_COORD_TYPE> coord(dimension);

  const COLOR_TYPE red[3] = { 1.0, 0.0, 0.0 };
  const COLOR_TYPE blue[3] = { 0.0, 0.0, 1.0 };
  const COLOR_TYPE green[3] = { 0.0, 1.0, 0.0 };

  VERTEX_INDEX icube = 0;
  int parity = 0;
  COLOR_TYPE * color_ptr = color;
  for (unsigned int i = 0; i < cube_list.size(); i++) {
    int new_cube = cube_list[i];
    if (icube != new_cube) {
      icube = new_cube;
      grid.ComputeCoord(icube, coord.Ptr());
      int sum = 0;
      for (int d = 0; d < dimension; d++) 
        { sum += coord[d]; }
      parity = sum%2;
    }

    if (parity == 0) 
      { std::copy(red, red+3, color_ptr); }
    else
      { std::copy(blue, blue+3, color_ptr); }

    // set opacity
    color_ptr[3] = 1.0;
    color_ptr += 4;
  }
  
}

