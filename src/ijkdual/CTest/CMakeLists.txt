PROJECT(TESTIJKDUAL)

#---------------------------------------------------------

ENABLE_TESTING()

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED IJK_DIR)
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED IJK_DIR)

IF (NOT DEFINED IJK_DATA_DIR)
  SET(IJK_DATA_DIRA "${IJK_DIR}/data")
  SET(IJK_DATA_DIRB "$ENV{HOME}/Data/Nrrd")
  SET(IJK_DATA_DIRC "/c/users/$ENV{USERNAME}/Data/Nrrd")
  IF (EXISTS ${IJK_DATA_DIRA})
    SET(IJK_DATA_DIR ${IJK_DATA_DIRA})
  ELSEIF (EXISTS ${IJK_DATA_DIRB})
    SET(IJK_DATA_DIR ${IJK_DATA_DIRB})
    ELSEIF (EXISTS ${IJK_DATA_DIRC})
    SET(IJK_DATA_DIR ${IJK_DATA_DIRC})
  ELSE()
    MESSAGE( FATAL_ERROR "Unable to find data directory. Directories ${IJK_DATA_DIRA} and ${IJK_DATA_DIRB} and ${IJK_DATA_DIRC} do not exist.")
  ENDIF()
  
ENDIF (NOT DEFINED IJK_DATA_DIR)


# Compare two programs, ${NEW_PROGRAM} and ${PREV_PROGRAM}.
IF (NOT DEFINED NEW_PROGRAM)
  SET(NEW_PROGRAM "../ijkdual")
ENDIF()

IF (NOT DEFINED PREV_PROGRAM)
  SET(PREV_PROGRAM "../ijkdual-ctest-prototype")
ENDIF()

MESSAGE(STATUS "This ctest compares new program \"${NEW_PROGRAM}\" against existing program \"${PREV_PROGRAM}\".\n")


#---------------------------------------------------------
# MACROS
#---------------------------------------------------------

# Add comparison test. Specify output file suffix.
MACRO(ADD_COMPARISON_TEST_X tname cargs isovalue infile outfile_suffix)

  IF (DEFINED NO_COMMENTS)
    # Set and use NO_COMMENTS to avoid unused variable warning.
    SET(NO_COMMENTS "-no_comments")
    SET(arg_list "${cargs} ${NO_COMMENTS} ${isovalue} ${infile}")
  ELSE()
    SET(arg_list "${cargs} ${isovalue} ${infile}")
  ENDIF()

  ADD_TEST( ${tname} ${CMAKE_COMMAND} -Dprog0=${NEW_PROGRAM} -Dprog1=${PREV_PROGRAM} -Dtest_arg=${arg_list} -Doutput0=out0.${outfile_suffix} -Doutput1=out1.${outfile_suffix} -Dflag_write_to_file=true -P ${CMAKE_SOURCE_DIR}/compare_prog.cmake )

ENDMACRO()


# Add standard comparison test. Output to .off files.
MACRO(ADD_COMPARISON_TEST tname cargs isovalue infile)
  ADD_COMPARISON_TEST_X(${tname} "${cargs}" ${isovalue} ${infile} "off")
ENDMACRO()


# Add comparison test on tdata.
MACRO(ADD_COMPARISON_TEST_ON_DATA_X tdata test_suffix cargs output_suffix)

  IF ( "${test_suffix}" STREQUAL "")
    SET( tname "${tdata}" )
  ELSE()
    SET( tname ${tdata}_${test_suffix} )
  ENDIF()

  ADD_COMPARISON_TEST_X( ${tname} "${cargs}" ${${tdata}_isovalue} ${${tdata}_file}
         	         ${output_suffix} )
  SET_PROPERTY( TEST ${tname} PROPERTY LABELS ${${tdata}_speed} )

  IF ( "${test_suffix}" STREQUAL "")
    SET_PROPERTY( TEST ${tname} APPEND PROPERTY LABELS "basic" )
  ENDIF()

ENDMACRO()


# Set FLAG_SKIP_TEST to TRUE for certain combinations of tests and cargs. 
FUNCTION(QUERY_SKIP_TEST tdata cargs)

  # Default
  SET(FLAG_SKIP_TEST FALSE PARENT_SCOPE)

   IF ( "${cargs}" MATCHES ".*-supersample.*" )
     IF ( NOT ("${${tdata}_speed}" MATCHES "^veryfast" ) )
       SET(FLAG_SKIP_TEST TRUE PARENT_SCOPE)
    ENDIF()
  ENDIF()

ENDFUNCTION()


# Add comparison teston all data. Specify output file suffix.
MACRO(ADD_COMPARISON_TEST_ON_ALL_DATA_X test_suffix cargs output_suffix)

  FOREACH (tdata ${testdata_list})

    QUERY_SKIP_TEST("${tdata}" "${cargs}") 
    IF (${FLAG_SKIP_TEST})
      CONTINUE()
    ENDIF()

    ADD_COMPARISON_TEST_ON_DATA_X( "${tdata}" "${test_suffix}" 
    				   "${cargs}" "${output_suffix}" )

  ENDFOREACH(tdata)

ENDMACRO()


MACRO(ADD_COMPARISON_TEST_ON_ALL_DATA test_suffix cargs)
  ADD_COMPARISON_TEST_ON_ALL_DATA_X( "${test_suffix}" "${cargs}" "off" )
ENDMACRO()


# Add comparison test on all 2D data. Specify output file suffix.
MACRO(ADD_COMPARISON_TEST_ON_ALL_2D_DATA_X test_suffix cargs output_suffix)

  FOREACH (tdata ${test2Ddata_list})

    QUERY_SKIP_TEST("${tdata}" "${cargs}") 
    IF (${FLAG_SKIP_TEST})
      CONTINUE()
    ENDIF()

    IF (NOT(${${tdata}_dimension} EQUAL 2))
      CONTINUE()
    ENDIF()

    ADD_COMPARISON_TEST_ON_DATA_X( "${tdata}" "${test_suffix}" 
    				   "${cargs}" "${output_suffix}" )

  ENDFOREACH(tdata)

ENDMACRO()


# Add comparison test on all 3D data. Specify output file suffix.
MACRO(ADD_COMPARISON_TEST_ON_ALL_3D_DATA_X test_suffix cargs output_suffix)

  FOREACH (tdata ${test3Ddata_list})
    QUERY_SKIP_TEST("${tdata}" "${cargs}") 
    IF (${FLAG_SKIP_TEST})
      CONTINUE()
    ENDIF()

    IF (NOT(${${tdata}_dimension} EQUAL 3))
      CONTINUE()
    ENDIF()

    ADD_COMPARISON_TEST_ON_DATA_X( "${tdata}" "${test_suffix}" 
    				   "${cargs}" "${output_suffix}" )

  ENDFOREACH(tdata)

ENDMACRO()


# Add comparison test on all 3D data.
MACRO(ADD_COMPARISON_TEST_ON_ALL_3D_DATA test_suffix cargs)
  ADD_COMPARISON_TEST_ON_ALL_3D_DATA_X("${test_suffix}" "${cargs}" "off")
ENDMACRO()


MACRO(SET_DATA_PROP tdata isovalue dimension tfile tspeed)
  SET(${tdata}_isovalue ${isovalue})
  SET(${tdata}_file ${tfile})
  SET(${tdata}_dimension ${dimension})
  SET(${tdata}_speed ${tspeed})
ENDMACRO()

MACRO(SET_2D_DATA_PROP tdata isovalue tfile)
  SET_DATA_PROP("${tdata}" ${isovalue} 2 "${tfile}" "veryfastII")
ENDMACRO()

MACRO(SET_3D_DATA_PROP tdata isovalue tfile tspeed)
  SET_DATA_PROP("${tdata}" ${isovalue} 3 "${tfile}" ${tspeed})
ENDMACRO()


#---------------------------------------------------------
# SET VARIABLES
#---------------------------------------------------------

SET(volvis_slices_dir "${IJK_DATA_DIR}/2D/volvis_slices")
SET(volvis_dir "${IJK_DATA_DIR}/3D/volvis")
SET(vollib_dir "${IJK_DATA_DIR}/3D/vollib")

# 2D test sets
SET_2D_DATA_PROP( Silicium-Slice16z 120 ${volvis_slices_dir}/silicium-slice16z.nrrd )
SET_2D_DATA_PROP( Fuel-Slice32z 70 ${volvis_slices_dir}/fuel-slice32z.nrrd )
SET_2D_DATA_PROP( Neghip-Slice32z 100 ${volvis_slices_dir}/neghip-slice32z.nrrd )
SET_2D_DATA_PROP( MarschnerLobb-Slice21z 100 ${volvis_slices_dir}/marschnerlobb-slice21z.nrrd )
SET_2D_DATA_PROP( MarschnerLobb-Slice22z 100 ${volvis_slices_dir}/marschnerlobb-slice22z.nrrd )
SET_2D_DATA_PROP( MarschnerLobb-Slice26z 100 ${volvis_slices_dir}/marschnerlobb-slice26z.nrrd )
SET(test2Ddata_list Silicium-Slice16z Fuel-Slice32z Neghip-Slice32z 
    MarschnerLobb-Slice21z MarschnerLobb-Slice22z 
    MarschnerLobb-Slice26z )  

# 3D test sets, integer isovalues
SET_3D_DATA_PROP( Silicium120 120 ${volvis_dir}/silicium.nhdr "veryfastII" )
SET_3D_DATA_PROP( Fuel70 70 ${volvis_dir}/fuel.nhdr "veryfastII" )
SET_3D_DATA_PROP( Neghip100 100 ${volvis_dir}/neghip.nhdr "veryfastII" )
SET_3D_DATA_PROP( MarschnerLobb100 100 ${volvis_dir}/marschnerlobb.nhdr "veryfastII" )
SET_3D_DATA_PROP( Engine80 80 "${volvis_dir}/engine.nhdr" "fast" )
SET_3D_DATA_PROP( Aneurysm30 30 "${volvis_dir}/aneurysm.nhdr" "fastII" )
SET_3D_DATA_PROP( Frog45 45 "${vollib_dir}/Frog.nhdr" "fastII" )
SET_3D_DATA_PROP( Lobster20 20 "${volvis_dir}/lobster.nhdr" "fast" )
SET_3D_DATA_PROP( Bonsai20 20 "${volvis_dir}/bonsai.nhdr" "slow" )
SET(test3DdataI_list Silicium120 Fuel70 Neghip100 MarschnerLobb100
    Engine80 Aneurysm30 Frog45 Lobster20 Bonsai20 )
  
# 3D test sets, non-integer isovalues
SET_3D_DATA_PROP( Silicium120p5 120.5 ${volvis_dir}/silicium.nhdr "veryfastII" )
SET_3D_DATA_PROP( Fuel70p5 70.5 ${volvis_dir}/fuel.nhdr "veryfastII" )
SET_3D_DATA_PROP( Neghip100p5 100.5 ${volvis_dir}/neghip.nhdr "veryfastII" )
SET_3D_DATA_PROP( MarschnerLobb100p5 100.5 ${volvis_dir}/marschnerlobb.nhdr "veryfastII" )
SET_3D_DATA_PROP( Engine80p5 80.5 "${volvis_dir}/engine.nhdr" "fastII" )
SET_3D_DATA_PROP( Aneurysm30p5 30.5 "${volvis_dir}/aneurysm.nhdr" "fastII" )
SET_3D_DATA_PROP( Frog45p5 45.5 "${vollib_dir}/Frog.nhdr" "fastII" )
SET_3D_DATA_PROP( Lobster20p5 20.5 "${volvis_dir}/lobster.nhdr" "fast" )
SET_3D_DATA_PROP( Bonsai20p5 20.5 "${volvis_dir}/bonsai.nhdr" "slow" )
SET(test3DdataNI_list Silicium120p5 Fuel70p5 Neghip100p5 MarschnerLobb100p5
    Engine80p5 Aneurysm30p5 Frog45p5 Lobster20p5 Bonsai20p5 )  


# All test sets
SET(test3Ddata_list ${test3DdataI_list} ${test3DdataNI_list})
SET(testdata_list ${test2Ddata_list} ${test3Ddata_list})


#---------------------------------------------------------
# ADD TESTS
#---------------------------------------------------------

ADD_COMPARISON_TEST_ON_ALL_DATA("" "")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA_X("ply" "-ply" "ply")
ADD_COMPARISON_TEST_ON_ALL_2D_DATA_X("fig" "-fig" "fig")
ADD_COMPARISON_TEST_ON_ALL_DATA("centroid" "-position centroid")
ADD_COMPARISON_TEST_ON_ALL_DATA("center" "-position cube_center")
ADD_COMPARISON_TEST_ON_ALL_DATA("sep_neg" "-sep_neg")
ADD_COMPARISON_TEST_ON_ALL_DATA("sep_pos" "-sep_pos")
ADD_COMPARISON_TEST_ON_ALL_DATA("manifold" "-manifold")
ADD_COMPARISON_TEST_ON_ALL_DATA("multi_isov" "-multi_isov")
ADD_COMPARISON_TEST_ON_ALL_DATA("single_isov" "-single_isov")
ADD_COMPARISON_TEST_ON_ALL_DATA("single_isov_centroid" "-single_isov -position centroid")
ADD_COMPARISON_TEST_ON_ALL_DATA("single_isov_cube_center" "-single_isov -position cube_center")
ADD_COMPARISON_TEST_ON_ALL_DATA("select_split" "-select_split")
ADD_COMPARISON_TEST_ON_ALL_DATA("connect_ambig" "-connect_ambig")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh" "-trimesh")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA_X("trimesh_ply" "-trimesh -ply" "ply")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_uniform" "-trimesh_uniform")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_split_max" "-trimesh_split_max")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_single_isov" "-trimesh -single_isov")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_uniform_single_isov" "-trimesh_uniform -single_isov")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4" "-tri4")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_max_min_angle" "-trimesh_tri4_max_min_angle")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_single_isov" "-tri4 -single_isov")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_all" "-trimesh_tri4_all")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_all_single_isov" "-trimesh_tri4_all -single_isov")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_all_qei_bilinear" "-trimesh_tri4_all -tri4_on_grid_edge -qei_bilinear")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_all_qei_average" "-trimesh_tri4_all -tri4_on_grid_edge -qei_average")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_all_qei_interpolateS" "-trimesh_tri4_all -tri4_on_grid_edge -qei_interpolate_scalar")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope" "-envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope_single_isov" "-envelope -single_isov")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope_only_tri2" "-envelope_only_tri2")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope_only_tri4" "-envelope_only_tri4")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope_tri2_or_tri4_prefer_tri2" "-envelope_tri2_or_tri4_prefer_tri2")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope_prefer_tri2_qei_bilinear" "-envelope_tri2_or_tri4_prefer_tri2 -tri4_on_grid_edge -qei_bilinear")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope_prefer_tri2_qei_average" "-envelope_tri2_or_tri4_prefer_tri2 -tri4_on_grid_edge -qei_average")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("envelope_prefer_tri2_qei_interpolateS" "-envelope_tri2_or_tri4_prefer_tri2 -tri4_on_grid_edge -qei_interpolate_scalar")

ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_envelope" "-trimesh -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_envelope_only_tri2" "-trimesh -envelope_only_tri2")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_envelope_only_tri4" "-trimesh -envelope_only_tri4")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_envelope_tri2_or_tri4_prefer_tri2" "-trimesh -envelope_tri2_or_tri4_prefer_tri2")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("trimesh_envelope_single_isov" "-trimesh -envelope -single_isov")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_envelope" "-tri4 -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_envelope_only_tri2" "-tri4 -envelope_only_tri2")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_envelope_only_tri4" "-tri4 -envelope_only_tri4")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_envelope_tri2_or_tri4_prefer_tri2" "-tri4 -envelope_tri2_or_tri4_prefer_tri2")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("tri4_envelope" "-tri4 -envelope -single_isov")
ADD_COMPARISON_TEST_ON_ALL_DATA("subsample" "-subsample 2")
ADD_COMPARISON_TEST_ON_ALL_DATA("supersample" "-supersample 2")
ADD_COMPARISON_TEST_ON_ALL_DATA("separate_isov" "-separate_isov")
ADD_COMPARISON_TEST_ON_ALL_DATA("trimesh_separate_isov" "-trimesh -separate_isov")
ADD_COMPARISON_TEST_ON_ALL_DATA("separate_by_cc" "-separate_by_cc")
ADD_COMPARISON_TEST_ON_ALL_DATA("move_isov_away_from_facets" "-move_isov_away_from_facets")
ADD_COMPARISON_TEST_ON_ALL_DATA("trimesh_move_isov_away_from_facets" "-trimesh -move_isov_away_from_facets")  
ADD_COMPARISON_TEST_ON_ALL_DATA("random_pos" "-position random -random_seed 100 -rpos_sample_boundary")
ADD_COMPARISON_TEST_ON_ALL_DATA("random_pos_single_isov" "-single_isov -position random -random_seed 100 -rpos_sample_boundary")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_envelope" "-position random -random_seed 100 -envelope -rpos_sample_boundary")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_envelope_single_isov" "-position random -random_seed 100 -rpos_sample_boundary -envelope -single_isov")  
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_trimesh_envelope" "-position random -random_seed 100 -rpos_sample_boundary -trimesh -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_trimesh_envelope_single_isov" "-position random -random_seed 100 -rpos_sample_boundary -trimesh -envelope -single_isov")  
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_tri4_envelope" "-position random -random_seed 100 -rpos_sample_boundary -tri4 -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_tri4_envelope_single_isov" "-position random -random_seed 100 -rpos_sample_boundary -tri4 -envelope -single_isov")
ADD_COMPARISON_TEST_ON_ALL_DATA("rpos_move" "-position random -random_seed 100 -rpos_sample_boundary -move_isov_away_from_facets")
ADD_COMPARISON_TEST_ON_ALL_DATA("rpos_move_single_isov" "-single_isov -position random -random_seed 100 -rpos_sample_boundary -move_isov_away_from_facets")
ADD_COMPARISON_TEST_ON_ALL_DATA("rpos_Uquadratric" "-position random -rpos_distribution Uquadratic -random_seed 100")
ADD_COMPARISON_TEST_ON_ALL_DATA("rpos_Uquadratic_single_isov" "-single_isov -position random -rpos_distribution Uquadratic -random_seed 100")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rposUQ_envelope" "-position random -rpos_distribution Uquadratic -random_seed 100 -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rposUQ_envelope_single_isov" "-single_isov -position random -rpos_distribution Uquadratic -random_seed 100 -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rposUQ_trimesh_envelope" "-position random -rpos_distribution Uquadratic -random_seed 100 -trimesh -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rposUQ_trimesh_envelope_single_isov" "-single_isov -position random -rpos_distribution Uquadratic -random_seed 100 -trimesh -envelope")
  ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rposUQ_tri4_envelope" "-position random -rpos_distribution Uquadratic -random_seed 100 -tri4 -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rposUQ_tri4_envelope_single_isov" "-single_isov -position random -rpos_distribution Uquadratic -random_seed 100 -tri4 -envelope")
ADD_COMPARISON_TEST_ON_ALL_DATA("rpos_separate_isov" "-position random -random_seed 100 -rpos_sample_boundary -separate_isov")
ADD_COMPARISON_TEST_ON_ALL_DATA("rpos_separate_by_cc" "-position random -random_seed 100 -rpos_sample_boundary -separate_by_cc")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_trimesh_separate_isov" "-position random -random_seed 100 -rpos_sample_boundary -separate_isov -trimesh -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_trimesh_separate_by_cc" "-position random -random_seed 100 -rpos_sample_boundary -separate_by_cc -trimesh -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_trimesh_move" "-position random -random_seed 100 -rpos_sample_boundary -move_isov_away_from_facets -trimesh -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_trimesh_sep_nearF" "-position random -random_seed 100 -rpos_sample_boundary -sep_isov_near_cube_facets -trimesh -envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_trimesh_sep_nearF_tri4inE" "-position random -random_seed 100 -rpos_sample_boundary -sep_isov_near_cube_facets -trimesh -envelope -tri4_in_envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_tri4_sep_nearF_tri4inE" "-position random -random_seed 100 -rpos_sample_boundary -sep_isov_near_cube_facets -tri4 -envelope -tri4_in_envelope")
ADD_COMPARISON_TEST_ON_ALL_3D_DATA("rpos_tri4_sep_nearB2_tri4inE" "-position random -random_seed 100 -rpos_sample_boundary -rpos_boundary_width 0.2 -sep_isov_near_cube_facets -tri4 -envelope -tri4_in_envelope")


 
