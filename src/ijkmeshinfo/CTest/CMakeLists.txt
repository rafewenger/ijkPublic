PROJECT(TESTIJKMESHINFO)

#---------------------------------------------------------

ENABLE_TESTING()

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})


#---------------------------------------------------------
# TEST MACROS
#---------------------------------------------------------

# Add "-mesh dim <mdim>" to arg_list if mesh_dim is defined.
MACRO(ADD_MESH_DIM arg_list tdata)

  SET(mesh_dim ${${tdata}_mesh_dim})

  IF (DEFINED mesh_dim)
    SET(arg_list "${${arg_list}} -mesh_dim ${mesh_dim}")
  ENDIF(DEFINED mesh_dim)

ENDMACRO()

# Append arg1 and ${tdata}_${arg2} to listX if ${tdata}_${arg2} is defined.
MACRO(APPEND_IF_DEFINED listX arg1 arg2 tdata)

  SET(X ${${tdata}_${arg2}})

  IF (DEFINED X)
    LIST(APPEND ${listX} "${arg1}" ${X})
  ENDIF(DEFINED X)
  
ENDMACRO()


# Append arg1 and ${tdata}_${arg2} to listX if ${tdata}_${arg2} is defined.
# Otherwise, append arg3.
MACRO(APPEND_SELECT listX arg1 arg2 arg3 tdata)

  SET(X ${${tdata}_${arg2}})

  IF (DEFINED X)
    LIST(APPEND ${listX} "${arg1}" ${X})
  ELSE (DEFINED X)
    LIST(APPEND ${listX} "${arg3}")
  ENDIF(DEFINED X)
  
ENDMACRO()


MACRO(ADD_RE_TEST tname commandX args infile reg_exp)

  SET(arg_list "${args} ${infile}")

  # Eliminate extra spaces
  STRING(REGEX REPLACE "[ ][ ]+" " " arg_list "${arg_list}")
  STRING(REGEX REPLACE "^[ ]+" "" arg_list "${arg_list}")

  SEPARATE_ARGUMENTS(arg_list)
  ADD_TEST( ${tname} ${commandX} ${arg_list} )
  SET_TESTS_PROPERTIES(${tname} PROPERTIES PASS_REGULAR_EXPRESSION "${reg_exp}")

ENDMACRO()


# Note: Strings in re_list should not contain ";".
MACRO(ADD_RE_LIST_TEST tname commandX args infile re_list)

  STRING(REGEX REPLACE ";" ".*" reg_exp "${re_list}")
  ADD_RE_TEST(${tname} ${commandX} "${args}" "${infile}" "${reg_exp}")

ENDMACRO()


MACRO(ADD_GENERAL_INFO_TEST tname commandX args tdata)

  SET(infile ${${tdata}_file})
  SET(vdim ${${tdata}_vdim})
  SET(mdim ${${tdata}_mdim})
  SET(numv ${${tdata}_numv})
  SET(nump ${${tdata}_nump})
  SET(poly_type ${${tdata}_poly_type})
  SET(lower_coord ${${tdata}_lower_coord})
  SET(upper_coord ${${tdata}_upper_coord})

  SET(re_list "Volume dimension" ${vdim} 
              "Mesh dimension" ${mdim} 
              "Number of mesh vertices" ${numv} )

  APPEND_IF_DEFINED(re_list "Number of mesh edges" nume ${tdata} )
  LIST(APPEND re_list "Number of mesh ${poly_type}" ${nump} )
  APPEND_IF_DEFINED(re_list "#V - #E [+] #F" euler_num ${tdata} )
  APPEND_IF_DEFINED(re_list "Number of.*with duplicate vertices"
                    num_poly_with_dup_vert ${tdata} )
  APPEND_IF_DEFINED(re_list "Number of duplicate" num_dup_poly ${tdata} )

  IF ((${poly_type} STREQUAL "polygons") OR
      (${poly_type} STREQUAL "triangles"))
    APPEND_IF_DEFINED(re_list "Min polygon angle" min_angle ${tdata} )
    APPEND_IF_DEFINED(re_list "Max polygon angle" max_angle ${tdata} )
    APPEND_IF_DEFINED(re_list "Min polygon edge length" 
                              min_edge_length ${tdata} )
    APPEND_IF_DEFINED(re_list "Max polygon edge length" 
                              max_edge_length ${tdata} )
  ELSE()
    APPEND_IF_DEFINED(re_list "Min edge length" 
                              min_edge_length ${tdata} )
    APPEND_IF_DEFINED(re_list "Max edge length" 
                              max_edge_length ${tdata} )
  ENDIF()

  IF (DEFINED lower_coord)
    LIST(APPEND re_list 
                "Bounding box" "${lower_coord}" "${upper_coord}" )
  ENDIF(DEFINED lower_coord)

  SET(arg_list "${args}")
  ADD_MESH_DIM(arg_list ${tdata})
  ADD_RE_LIST_TEST(
        ${tname} ${commandX} "${arg_list}" "${infile}" "${re_list}" )

ENDMACRO()


MACRO (ADD_GENERAL_INFO_TEST_ON_DATA commandX args suffix tdata_list)

  FOREACH (tdata ${tdata_list})
    ADD_GENERAL_INFO_TEST(${tdata}_${suffix} ${commandX} "${args}" ${tdata})
  ENDFOREACH(tdata)

ENDMACRO()


MACRO(ADD_MANIFOLD_TEST tname commandX args tdata)

  SET(infile ${${tdata}_file})
  SET(num_degenerate_poly ${${tname}_num_degenerate_poly})

  IF (${${tdata}_mdim} LESS 3)
    SET(poly_name "polygons")
  ELSE (${${tdata}_mdim} LESS 3)
    SET(poly_name "polytopes")
  ENDIF(${${tdata}_mdim} LESS 3)

  SET(re_list "")
  APPEND_SELECT(re_list "Num degenerate ${poly_name}" 
                        num_degenerate_poly 
                        "No degenerate ${poly_name}" ${tdata} )
  APPEND_SELECT(re_list "Num duplicate ${poly_name}"
                        num_dup_poly "No duplicate ${poly_name}" ${tdata} )

  IF (${${tdata}_mdim} GREATER 2)
    APPEND_SELECT(re_list "Num non-manifold facets"
                          num_nonmanifold_facets
                          "No non-manifold facets" ${tdata} )
  ENDIF(${${tdata}_mdim} GREATER 2)

  IF (${${tdata}_mdim} EQUAL 2)
    APPEND_SELECT(re_list "Num non-manifold edges"
                          num_nonmanifold_edges 
                          "No non-manifold edges" ${tdata} )
  ENDIF(${${tdata}_mdim} EQUAL 2)

  APPEND_SELECT(re_list "Num non-manifold vertices"
                num_nonmanifold_vertices "No non-manifold vertices" ${tdata} )


  IF (${${tdata}_mdim} LESS ${${tdata}_vdim})

    IF (${${tdata}_flag_boundary})

      IF (DEFINED ${tdata}_num_facets_inside_BB)
        LIST(APPEND re_list "Number of surface boundary facets inside bounding box" ${${tdata}_num_facets_inside_BB})
        LIST(APPEND re_list "Number of surface boundary facets at least 1 from bounding box boundary" ${${tdata}_num_facets_deep_inside_BB})
      ELSE (DEFINED ${tdata}_num_facets_inside_BB)
        LIST(APPEND re_list "Surface boundary lies on bounding box boundary")
      ENDIF(DEFINED ${tdata}_num_facets_inside_BB)

    ELSE (${${tdata}_flag_boundary})
      LIST(APPEND re_list "Surface has no boundary")
    ENDIF(${${tdata}_flag_boundary})

  ENDIF(${${tdata}_mdim} LESS ${${tdata}_vdim})


  SET(arg_list "-manifold ${args}")
  ADD_MESH_DIM(arg_list ${tdata})
  ADD_RE_LIST_TEST(
        ${tname} ${commandX} "${arg_list}" "${infile}" "${re_list}" )

ENDMACRO()


MACRO (ADD_MANIFOLD_TEST_ON_DATA commandX args suffix tdata_list)

  FOREACH (tdata ${tdata_list})
    ADD_MANIFOLD_TEST(${tdata}_${suffix} ${commandX} "${args}" ${tdata})
  ENDFOREACH(tdata)

ENDMACRO()


MACRO(ADD_TERSE_MANIFOLD_TEST tname commandX args tdata)

  SET(infile ${${tdata}_file})
  SET(num_degenerate_poly ${${tname}_num_degenerate_poly})

  SET(arg_list "-manifold ${args}")
  ADD_MESH_DIM(arg_list ${tdata})
  ADD_RE_LIST_TEST(
        ${tname} ${commandX} "${arg_list}" "${infile}" "${re_list}" )

ENDMACRO()


#---------------------------------------------------------
# COMPARISON TEST MACROS
#---------------------------------------------------------

MACRO(ADD_COMPARISON_TEST tname cargs infile)

  SET(arg_list "${cargs} ${infile}")

  ADD_TEST( ${tname} ${CMAKE_COMMAND} -Dprog0=${prog0} -Dprog1=${prog1} -Dtest_arg=${arg_list} -Doutput0=out0.off -Doutput1=out1.off -Dno_fail_on_return1=1 -P ${CMAKE_SOURCE_DIR}/compare_prog.cmake )
ENDMACRO()

MACRO(ADD_COMPARISON_TEST_ON_DATA test_suffix cargs tdata_list)

  FOREACH (tdata ${tdata_list})

    IF ( "${test_suffix}" STREQUAL "")
      SET( tname "${tdata}" )
    ELSE ( "${test_suffix}" STREQUAL "") 
      SET( tname ${tdata}_${test_suffix} )
    ENDIF ( "${test_suffix}" STREQUAL "") 

    SET(arg_list "${cargs}")
    ADD_MESH_DIM(arg_list ${tdata})

    ADD_COMPARISON_TEST( ${tname} "${arg_list}" ${${tdata}_file} )

  ENDFOREACH(tdata)

ENDMACRO()


MACRO(ADD_COMPARE_OUTPUT_FILES_TEST tname cargs infile suffix_list)

  # Change ";" to " " to avoid list expansion with ";"
  STRING(REPLACE ";" " " suffix_listB "${suffix_list}")

  SET(arg_list "${cargs}")

  ADD_TEST( ${tname} ${CMAKE_COMMAND} -Dprog0=${prog0} -Dprog1=${prog1} -Dinfile=${infile} -Dtest_arg=${arg_list} -Dsuffix_list=${suffix_listB} -Dno_fail_on_return1=1 -P ${CMAKE_SOURCE_DIR}/compare_output_files.cmake )

ENDMACRO()


MACRO(ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA 
      test_suffix cargs suffix_list tdata_list)

  FOREACH (tdata ${tdata_list})

    IF ( "${test_suffix}" STREQUAL "")
      SET( tname "${tdata}" )
    ELSE ( "${test_suffix}" STREQUAL "") 
      SET( tname ${tdata}_${test_suffix} )
    ENDIF ( "${test_suffix}" STREQUAL "") 

    SET(arg_list "${cargs}")
    ADD_MESH_DIM(arg_list ${tdata})

    ADD_COMPARE_OUTPUT_FILES_TEST( 
      ${tname} "${arg_list}" ${${tdata}_file} "${suffix_list}" )

  ENDFOREACH(tdata)

ENDMACRO()


#---------------------------------------------------------
# SET MACROS
#---------------------------------------------------------


MACRO(SET_EULER_NUM tdata euler_num)
  SET(${tdata}_euler_num ${euler_num})
ENDMACRO()

MACRO(SET_MESH_DIM tdata mesh_dim)
  SET(${tdata}_mesh_dim ${mesh_dim})
ENDMACRO()


MACRO(SET_TRIANGLE_PROP 
      tdata tfile vdim mdim numv nume nump euler_num)

  SET(${tdata}_file ${tfile})
  SET(${tdata}_vdim ${vdim})
  SET(${tdata}_mdim ${mdim})
  SET(${tdata}_numv ${numv})
  SET(${tdata}_nume ${nume})
  SET(${tdata}_nump ${nump})
  SET(${tdata}_poly_type "triangles")
  SET_EULER_NUM(${tdata} ${euler_num})
ENDMACRO()


MACRO(SET_SIMPLEX_PROP 
      tdata tfile vdim mdim numv nume nump)

  SET(${tdata}_file ${tfile})
  SET(${tdata}_vdim ${vdim})
  SET(${tdata}_mdim ${mdim})
  SET(${tdata}_numv ${numv})
  SET(${tdata}_nume ${nume})
  SET(${tdata}_nump ${nump})
  SET(${tdata}_poly_type "simplices")
ENDMACRO()


MACRO(SET_TETRAHEDRON_PROP 
      tdata tfile vdim numv nume nump euler_num)

  SET(${tdata}_file ${tfile})
  SET(${tdata}_vdim ${vdim})
  SET(${tdata}_mdim 3)
  SET(${tdata}_numv ${numv})
  SET(${tdata}_nume ${nume})
  SET(${tdata}_nump ${nump})
  SET(${tdata}_poly_type "simplices")
  SET_EULER_NUM(${tdata} ${euler_num})

  SET_MESH_DIM(${tdata} 3)

ENDMACRO()


MACRO(SET_QUAD_PROP 
      tdata tfile vdim mdim numv nump )

  SET(${tdata}_file ${tfile})
  SET(${tdata}_vdim ${vdim})
  SET(${tdata}_mdim ${mdim})
  SET(${tdata}_numv ${numv})
  SET(${tdata}_nump ${nump})
  SET(${tdata}_poly_type "polygons")

  SET_MESH_DIM(${tdata} 2)

ENDMACRO()


MACRO(SET_POLYGON_PROP 
      tdata tfile vdim mdim numv nump )

  SET(${tdata}_file ${tfile})
  SET(${tdata}_vdim ${vdim})
  SET(${tdata}_mdim ${mdim})
  SET(${tdata}_numv ${numv})
  SET(${tdata}_nump ${nump})
  SET(${tdata}_poly_type "polygons")

  SET_MESH_DIM(${tdata} ${mdim})

ENDMACRO()


MACRO(SET_HEXAHEDRA_PROP 
      tdata tfile vdim mdim numv nump )

  SET(${tdata}_file ${tfile})
  SET(${tdata}_vdim ${vdim})
  SET(${tdata}_mdim ${mdim})
  SET(${tdata}_numv ${numv})
  SET(${tdata}_nump ${nump})
  SET(${tdata}_poly_type "polytopes")

ENDMACRO()


MACRO(SET_ANGLE_PROP tdata amin amax)
  SET(${tdata}_min_angle ${amin})
  SET(${tdata}_max_angle ${amax})
ENDMACRO()


MACRO(SET_ELENGTH_PROP tdata lmin lmax)
  SET(${tdata}_min_edge_length ${lmin})
  SET(${tdata}_max_edge_length ${lmax})
ENDMACRO()


MACRO(SET_BBOX_PROP tdata lower_coord upper_coord)
  SET(${tdata}_lower_coord "${lower_coord}")
  SET(${tdata}_upper_coord "${upper_coord}")
ENDMACRO()


# Set angle, edge length and bounding box properties
MACRO(SET_ALB_PROP tdata amin amax lmin lmax lower_coord upper_coord)
  SET_ANGLE_PROP(${tdata} ${amin} ${amax})
  SET_ELENGTH_PROP(${tdata} ${lmin} ${lmax})
  SET_BBOX_PROP(${tdata} ${lower_coord} ${upper_coord})
ENDMACRO()

# Set edge length and bounding box properties
MACRO(SET_ELENGTH_BBOX_PROP tdata lmin lmax lower_coord upper_coord)
  SET_ELENGTH_PROP(${tdata} ${lmin} ${lmax})
  SET_BBOX_PROP(${tdata} ${lower_coord} ${upper_coord})
ENDMACRO()

MACRO(SET_NUM_POLY_WITH_DUP_VERT tdata nump)
  SET(${tdata}_num_poly_with_dup_vert ${nump})
ENDMACRO()

MACRO(SET_NUM_DUP_POLY tdata nump)
  SET(${tdata}_num_dup_poly ${nump})
ENDMACRO()


# Set ${tdata}_${field} to ${value} if ${value} is non-zero.
MACRO(SET_IF_NONZERO tdata field value)

   IF (NOT ${value} EQUAL 0)
     SET(${tdata}_${field} ${value})
   ENDIF(NOT ${value} EQUAL 0)

ENDMACRO()

MACRO(SET_NON_MANIFOLD2D tdata num_degenerate_poly num_dup_poly 
                         num_nonmanifold_edges num_nonmanifold_vertices)

  SET_IF_NONZERO(${tdata} num_degenerate_poly ${num_degenerate_poly})
  SET_IF_NONZERO(${tdata} num_dup_poly ${num_dup_poly})
  SET_IF_NONZERO(${tdata} num_nonmanifold_edges ${num_nonmanifold_edges})
  SET_IF_NONZERO(${tdata} num_nonmanifold_vertices ${num_nonmanifold_vertices})

ENDMACRO()

MACRO(SET_NON_MANIFOLD3D tdata num_degenerate_poly num_dup_poly 
                         num_nonmanifold_facets 
                         num_nonmanifold_edges 
                         num_nonmanifold_vertices)

  SET_IF_NONZERO(${tdata} num_degenerate_poly ${num_degenerate_poly})
  SET_IF_NONZERO(${tdata} num_dup_poly ${num_dup_poly})
  SET_IF_NONZERO(${tdata} num_nonmanifold_facets ${num_nonmanifold_facets})
  SET_IF_NONZERO(${tdata} num_nonmanifold_edges ${num_nonmanifold_edges})
  SET_IF_NONZERO(${tdata} num_nonmanifold_vertices ${num_nonmanifold_vertices})

ENDMACRO()


MACRO(SET_NON_MANIFOLD_ND tdata num_degenerate_poly num_dup_poly 
                          num_nonmanifold_facets 
                          num_nonmanifold_edges 
                          num_nonmanifold_vertices)

  SET_IF_NONZERO(${tdata} num_degenerate_poly ${num_degenerate_poly})
  SET_IF_NONZERO(${tdata} num_dup_poly ${num_dup_poly})
  SET_IF_NONZERO(${tdata} num_nonmanifold_facets ${num_nonmanifold_facets})
  SET_IF_NONZERO(${tdata} num_nonmanifold_edges ${num_nonmanifold_edges})
  SET_IF_NONZERO(${tdata} num_nonmanifold_vertices ${num_nonmanifold_vertices})

ENDMACRO()


MACRO(SET_NON_MANIFOLD_HEX 
      tdata num_degenerate_poly num_dup_poly num_nonmanifold_facets
      num_nonmanifold_edges num_nonmanifold_vertices)

  SET_NON_MANIFOLD3D(${tdata} ${num_degenerate_poly} ${num_dup_poly}
                              ${num_nonmanifold_facets} 
                              ${num_nonmanifold_edges} 
                              ${num_nonmanifold_vertices} )
ENDMACRO()

MACRO(SET_BOUNDARY tdata num_facets_inside_BB num_facets_deep_inside_BB)

  SET(${tdata}_flag_boundary true)

  IF (NOT ${num_facets_inside_BB} EQUAL 0)
    SET(${tdata}_num_facets_inside_BB  ${num_facets_inside_BB})
    SET(${tdata}_num_facets_deep_inside_BB  ${num_facets_deep_inside_BB})
  ENDIF(NOT ${num_facets_inside_BB} EQUAL 0)

ENDMACRO()

MACRO(SET_BOUNDARY_FLAG tdata flag)
  SET(${tdata}_flag_boundary ${flag})
ENDMACRO()


#---------------------------------------------------------
# SET VARIABLES
#---------------------------------------------------------

SET(prog0 "../ijkmeshinfo")
SET(prog1 "../ijkmeshinfo-ctest-prototype")
SET(test_dir "../TestData")

SET(min_max_angle_suffixes "min_poly_angle_freq.gplt" 
                           "max_poly_angle_freq.gplt")
SET(edge_length_suffixes "edge_length_freq.gplt")
SET(edge_length_internalE_suffixes "edge_length_freq.internalE.gplt")
SET(min_max_jacobian_suffixes "min_jdet_freq.gplt" "max_jdet_freq.gplt"
                              "min_normalized_jdet_freq.gplt"
                              "max_normalized_jdet_freq.gplt")
SET(min_max_jdet_vJacobian_suffixes 
       "min_jdet_freq.V.gplt" "max_jdet_freq.V.gplt"
       "min_normalized_jdet_freq.V.gplt" "max_normalized_jdet_freq.V.gplt")

SET(min_max_jdet_internalP_suffixes 
       "min_jdet_freq.internalP.gplt" "max_jdet_freq.internalP.gplt"
       "min_normalized_jdet_freq.internalP.gplt"
       "max_normalized_jdet_freq.internalP.gplt")
SET(min_max_jdet_internalV_suffixes 
       "min_jdet_freq.internalV.gplt" "max_jdet_freq.internalV.gplt"
       "min_normalized_jdet_freq.internalV.gplt"
       "max_normalized_jdet_freq.internalV.gplt")
SET(min_max_jshape_suffixes "min_jshape_freq.gplt" "max_jshape_freq.gplt")
SET(min_max_jshape_internalV_suffixes 
       "min_jshape_freq.internalV.gplt" "max_jshape_freq.internalV.gplt")
SET(min_max_jshape_vJacobian_suffixes 
       "min_jshape_freq.V.gplt" "max_jshape_freq.V.gplt")


SET_TRIANGLE_PROP(CubeT "${test_dir}/cubeT.off" 3 2 8 18 12 2)
SET_ALB_PROP(CubeT 45 90 1 1.41421 "(0,0,0)" "(1,1,1)")

SET_QUAD_PROP(CubeQ "${test_dir}/cubeQ.off" 3 2 8 6)
SET_ALB_PROP(CubeQ 90 90 1 1 "(0,0,0)" "(1,1,1)")

SET_QUAD_PROP(InternalQA "${test_dir}/internalQ.A.off" 3 2 8 5)
SET_ALB_PROP(InternalQA 78.9042 101.819 0.824621 1.03923 "(0,0,0)" "(1,1,1)")

SET_POLYGON_PROP(PolygonA "${test_dir}/polygon.A.off" 3 2 23 6)
SET_ALB_PROP(PolygonA 23.0215 165.964 1 5.09902 "(0,-4,0)" "(12,3,2)")

SET_TRIANGLE_PROP(DupVertTA "${test_dir}/dup_vertT.A.off" 3 2 8 18 12 2)
SET_ALB_PROP(DupVertTA 45 90 1 1.41421 "(0,0,0)" "(1,1,1)")
SET_NUM_POLY_WITH_DUP_VERT(DupVertTA 3)
SET_NON_MANIFOLD2D(DupVertTA 3 0 0 3)
SET_BOUNDARY_FLAG(DupVertTA true)

SET_QUAD_PROP(DupVertQA "${test_dir}/dup_vertQ.A.off" 3 2 8 6)
SET_ALB_PROP(DupVertQA 90 90 1 1 "(0,0,0)" "(1,1,1)")
SET_NUM_POLY_WITH_DUP_VERT(DupVertQA 2)
SET_NON_MANIFOLD2D(DupVertQA 2 0 0 0)
SET_BOUNDARY_FLAG(DupVertQA true)

SET_QUAD_PROP(DupPolyQA "${test_dir}/dup_polyQ.A.off" 3 2 8 6)
SET_ALB_PROP(dupPolyQQ 90 90 1 1 "(0,0,0)" "(1,1,1)")
SET_NUM_POLY_WITH_DUP_VERT(DupPolyQA 1)
SET_NUM_DUP_POLY(DupPolyQA 2)
SET_NON_MANIFOLD2D(DupPolyQA 1 2 0 0)
SET_BOUNDARY_FLAG(DupPolyQA true)

SET_QUAD_PROP(NonmanifoldEdgeQA 
              "${test_dir}/nonmanifold_edgeQ.A.off" 3 2 8 7)
SET_ALB_PROP(NonmanifoldEdgeQA 90 90 1 1.41421 "(0,0,0)" "(1,1,1)")
SET_NON_MANIFOLD2D(NonmanifoldEdgeQA 0 0 2 0)
SET_BOUNDARY_FLAG(NonmanifoldEdgeQA true)

SET_QUAD_PROP(NonmanifoldEdgeQB
              "${test_dir}/nonmanifold_edgeQ.B.off" 3 2 10 7)
SET_ALB_PROP(NonmanifoldEdgeQB 90 90 1 1 "(0,0,0)" "(2,1,1)")
SET_NON_MANIFOLD2D(NonmanifoldEdgeQB 0 0 1 0)
SET_BOUNDARY_FLAG(NonmanifoldEdgeQB true)

SET_TRIANGLE_PROP(NonmanifoldTwocubeTA
                   "${test_dir}/nonmanifold_twocubeT.A.off" 3 2 15 36 24 3)
SET_ALB_PROP(NonmanifoldTwocubeTA 45 90 1 1.41421 "(0,0,0)" "(2,2,2)")
SET_NON_MANIFOLD2D(NonmanifoldTwocubeTA 0 0 0 1)

SET_QUAD_PROP(NonmanifoldTwocubeQA 
              "${test_dir}/nonmanifold_twocubeQ.A.off" 3 2 15 12)
SET_ALB_PROP(NonmanifoldTwocubeQA 90 90 1 1 "(0,0,0)" "(2,2,2)")
SET_NON_MANIFOLD2D(NonmanifoldTwocubeQA 0 0 0 1)

SET_POLYGON_PROP(NonmanifoldVertPoly2DA 
                 "${test_dir}/nonmanifold_vertPoly2D.A.off" 3 2 6 2)
SET_NON_MANIFOLD2D(NonmanifoldVertPoly2DA 0 0 0 1)
SET_BOUNDARY_FLAG(NonmanifoldVertPoly2DA true)

SET_QUAD_PROP(IncorrectOrientQA
              "${test_dir}/incorrect_orientQ.A.off" 3 2 8 6)
SET_ALB_PROP(IncorrectOrientQA 90 90 1 1 "(0,0,0)" "(1,1,1)")

SET_TRIANGLE_PROP(BoundTA
                   "${test_dir}/boundT.A.off" 3 2 8 18 11 1)
SET_ALB_PROP(BoundTA 45 90 1 1.41421 "(0,0,0)" "(1,1,1)" )
SET_BOUNDARY_FLAG(BoundTA true)

SET_TRIANGLE_PROP(BoundTB
                   "${test_dir}/boundT.B.off" 3 2 260 648 384 -4)
SET_ALB_PROP(BoundTB 37.39 89.63 0.71 1.25 
                     "(0.143867,0.143867,0.143867)" 
                     "(8.85613,8.85613,8.85613)" )
SET_BOUNDARY(BoundTB 144 0)

SET_TRIANGLE_PROP(BoundTC
                   "${test_dir}/boundT.C.off" 3 2 1248 3732 2484 0)
SET_ALB_PROP(BoundTC 11.85 154.76 0.42 1.62
                     "(0.0463917,0.0463917,5.73368)"
                     "(18.9536,18.9536,13.2663)" )
SET_BOUNDARY(BoundTC 12 6)

SET_QUAD_PROP(BoundQA
              "${test_dir}/boundQ.A.off" 3 2 8 5)
SET_ALB_PROP(BoundQA 90 90 1 1 "(0,0,0)" "(1,1,1)" )
SET_BOUNDARY_FLAG(BoundQA true)

SET_QUAD_PROP(BoundQB
              "${test_dir}/boundQ.B.off" 3 2 260 192)
SET_ALB_PROP(BoundQB 52.8543 133.16 0.72 1.005 
                     "(0.143867,0.143867,0.143867)" 
                     "(8.85613,8.85613,8.85613)" )
SET_BOUNDARY(BoundQB 144 0)

SET_QUAD_PROP(BoundQC
              "${test_dir}/boundQ.C.off" 3 2 1248 1240)
SET_ALB_PROP(BoundQC 51.81 154.76 0.42 1.19
                     "(0.0463917,0.0463917,5.73368)"
                     "(18.9536,18.9536,13.2663)" )
SET_BOUNDARY(BoundQC 16 8)

SET_TRIANGLE_PROP(SelfITA
                   "${test_dir}/selfI.A.off" 3 2 6 6 2 2)
SET_ALB_PROP(SelfIA 47.12 79.48 1.11 1.5 "(0,0,0)" "(1,1,1)")
SET_BOUNDARY(SelfITA 2 0)

SET_TRIANGLE_PROP(SelfITB
                   "${test_dir}/selfI.B.off" 3 2 6 6 2 2)
SET_ALB_PROP(SelfITB 60 60 1.41 1.41 "(0,0,0)" "(1,1,1)")

SET_TRIANGLE_PROP(SelfITC
                   "${test_dir}/selfI.C.off" 3 2 6 6 2 2)
SET_ALB_PROP(SelfITC 36.86 90 1.5 2.82 "(0,0,0)" "(2,2,2)")

SET_TRIANGLE_PROP(SelfITD
                   "${test_dir}/selfI.D.off" 3 2 6 6 2 2)
SET_ALB_PROP(SelfITD 0 79.56 0 1.41 "(123.934,162,37.5)" "(124,163.5,39)")

SET_TRIANGLE_PROP(SelfITE
	          "${test_dir}/selfI.E.off" 3 2 8 12 6 2)
SET_ALB_PROP(SelfITE 60 60 1.41 1.41 "(0,0,0)" "(2,2,2)")

SET_TRIANGLE_PROP(Fuel35
            "${test_dir}/fuel.isov=35.off" 3 2 3143 9348 6224 19)
SET_ALB_PROP(Fuel35 0 134.732 0 1.52238 "(0.212121,20,20)" "(57.7031,43,43)")

SET_TRIANGLE_PROP(Neghip20
            "${test_dir}/neghip.isov=20.off" 3 2 22340 66722 44420 38)
SET_ALB_PROP(Neghip20 0 178.938 0 1.59844 "(0.0666667,5.2,0.251515)" "(63,57.1857,62.7485)")
  
SET_TRIANGLE_PROP(SphereCapTA
                   "${test_dir}/sphere_capT.A.off" 3 2 257 752 496 1)
SET_ALB_PROP(SphereCapTA 19.57 100.918 0.063 0.42 
                         "(0.133975,0.00214106,0.00214106)" 
                         "(2,1.99786,1.99786)")

SET_POLYGON_PROP(SphereCapPA
              "${test_dir}/sphere_capP.A.off" 3 2 257 256)
SET_ALB_PROP(SphereCapPA 22.42 100.91 0.063 0.389
                         "(0.133975,0.00214106,0.00214106)"
                         "(2,1.99786,1.99786)")

SET_TETRAHEDRON_PROP(CubeVolTetA
                     "${test_dir}/cubeVolTet.A.off" 3 8 19 6 -5)
SET_TETRAHEDRON_PROP(BallOctantTetA
                     "${test_dir}/ball_octantTet.A.off" 
                     3 19 63 28 -16)
SET_TETRAHEDRON_PROP(SlabTetA
                     "${test_dir}/slabTet.A.off" 3 18 57 24 -15)

SET_TETRAHEDRON_PROP(NonmanifoldFacetTetA
                     "${test_dir}/nonmanifold_facetTet.A.off"
                     3 9 12 3 0)
SET_NON_MANIFOLD3D(NonmanifoldFacetTetA 0 0 1 0 0)

SET_TETRAHEDRON_PROP(NonmanifoldVertTetA
                     "${test_dir}/nonmanifold_vertTet.A.off"
                     3 7 12 2 -3)
SET_NON_MANIFOLD3D(NonmanifoldVertTetA 0 0 0 0 1)

SET_TETRAHEDRON_PROP(NonmanifoldVertTetB
                     "${test_dir}/nonmanifold_vertTet.B.off"
                     3 9 24 8 -7)
SET_NON_MANIFOLD3D(NonmanifoldVertTetB 0 0 0 0 1)

SET_TETRAHEDRON_PROP(IncorrectOrientTetA
                     "${test_dir}/incorrect_orientTet.A.off" 
                     3 8 19 6 -5)

SET_TETRAHEDRON_PROP(IncorrectOrientTetB
                     "${test_dir}/incorrect_orientTet.B.off" 
                     3 8 19 6 -5)

SET_SIMPLEX_PROP(Simplex4D "${test_dir}/simplex4D.off" 4 4 5 10 1)
SET_ELENGTH_BBOX_PROP(Simplex4D 1 1.41421 "(0,0,0,0)" "(1,1,1,1)")

SET_SIMPLEX_PROP(Simp4DA "${test_dir}/simp4D.A.off" 4 4 6 14 2)
SET_ELENGTH_BBOX_PROP(Simp4DA 1 1.73205 "(0,0,0,0)" "(1,1,1,1)")

SET_SIMPLEX_PROP(Simplex5D "${test_dir}/simplex5D.off" 5 5 6 15 1)
SET_ELENGTH_BBOX_PROP(Simplex5D 1 1.41421 "(0,0,0,0)" "(1,1,1,1)")

SET_SIMPLEX_PROP(Simp5DA "${test_dir}/simp5D.A.off" 5 5 7 20 2)
SET_ELENGTH_BBOX_PROP(Simp5DA 1 2 "(0,0,0,0,0)" "(1,1,1,1,1)")

SET_SIMPLEX_PROP(NonmanifoldFacetSimp4DA 
                 "${test_dir}/nonmanifold_facetSimp4D.A.off"
                  4 4 7 18 3)
SET_NON_MANIFOLD_ND(NonmanifoldFacetSimp4DA 0 0 1 0 0)

SET_SIMPLEX_PROP(NonmanifoldVertSimp4DA 
                 "${test_dir}/nonmanifold_vertSimp4D.A.off"
                  4 4 9 20 2)
SET_NON_MANIFOLD_ND(NonmanifoldVertSimp4DA 0 0 0 0 1)

SET_SIMPLEX_PROP(NonmanifoldFacetSimp5DA 
                 "${test_dir}/nonmanifold_facetSimp5D.A.off"
                 5 5 8 25 3)
SET_NON_MANIFOLD_ND(NonmanifoldFacetSimp5DA 0 0 1 0 0)

SET_SIMPLEX_PROP(NonmanifoldVertSimp5DA 
                 "${test_dir}/nonmanifold_vertSimp5D.A.off"
                  5 5 11 30 2)
SET_NON_MANIFOLD_ND(NonmanifoldVertSimp5DA 0 0 0 0 1)

SET_SIMPLEX_PROP(IncorrectOrientSimp4DA
                 "${test_dir}/incorrect_orientSimp4D.A.off"
                  4 4 6 14 2)

SET_SIMPLEX_PROP(IncorrectOrientSimp5DA
                 "${test_dir}/incorrect_orientSimp5D.A.off"
                  5 5 7 20 2)

SET_HEXAHEDRA_PROP(SphereHexA
                   "${test_dir}/sphereHex.A.off" 3 3 342 200)

SET_HEXAHEDRA_PROP(Hexahedra4DA
                   "${test_dir}/hexahedra4D.A.off" 4 3 22 5)
SET_BOUNDARY(Hexahedra4DA 5 0)

SET_HEXAHEDRA_PROP(CubeVolHexA
                   "${test_dir}/cubeVolHex.A.off" 3 3 64 27)
SET_BBOX_PROP(CubeVolHexA "(0.5,0.5,0.5)" "(3.5,3.5,3.5)")

SET_HEXAHEDRA_PROP(CubeVolHexB
                   "${test_dir}/cubeVolHex.B.off" 3 3 252 150)
SET_BBOX_PROP(CubeVolHexB "(1,0.5,0.5)" "(5,5.5,5.5)")

SET_HEXAHEDRA_PROP(NonmanifoldVertHex3DA
                   "${test_dir}/nonmanifold_vertHex3D.A.off" 3 3 35 9)
SET_BBOX_PROP(NonmanifoldVertHex3DA "(0.5,0.5,0.5)" "(2.5,2.5,2.5)")
SET_NON_MANIFOLD_HEX(NonmanifoldVertHex3DA 0 0 0 0 1)

SET_HEXAHEDRA_PROP(NonmanifoldVertHex4DA
                   "${test_dir}/nonmanifold_vertHex4D.A.off" 4 3 28 2)
SET_BBOX_PROP(NonmanifoldVertHex4DA "(0.5,0.5,0.5,0.75)" "(2.5,2.5,2.5,1.5)")
SET_NON_MANIFOLD_HEX(NonmanifoldVertHex4DA 0 0 0 0 1)
SET_BOUNDARY(NonmanifoldVertHex4DA 6 0)

SET_HEXAHEDRA_PROP(NonmanifoldEdgeHex3DA
                   "${test_dir}/nonmanifold_edgeHex3D.A.off" 3 3 14 2)
SET_BBOX_PROP(NonmanifoldEdgeHex3DA "(0,0,0)" "(2,1,2)")
SET_NON_MANIFOLD_HEX(NonmanifoldEdgeHex3DA 0 0 0 1 0)

SET_HEXAHEDRA_PROP(NonmanifoldEdgeHex4DA
                   "${test_dir}/nonmanifold_edgeHex4D.A.off" 4 3 14 2)
SET_BBOX_PROP(NonmanifoldEdgeHex4DA 
              "(0.5,0.5,0.5,1.25)" "(2.5,2.5,1.375.1.8125)")
SET_NON_MANIFOLD_HEX(NonmanifoldEdgeHex4DA 0 0 0 1 0)
SET_BOUNDARY(NonmanifoldEdgeHex4DA 6 0)

SET_HEXAHEDRA_PROP(NonmanifoldEdgeHex4DB
                   "${test_dir}/nonmanifold_edgeHex4D.B.off" 4 3 28 6)
SET_BBOX_PROP(NonmanifoldEdgeHex4DB "(0.5,0.5,0.5,0.75)" "(2.5,2.5,2.5,1.5)")
SET_NON_MANIFOLD_HEX(NonmanifoldEdgeHex4DB 0 0 0 1 0)
SET_BOUNDARY(NonmanifoldEdgeHex4DB 6 0)

SET_HEXAHEDRA_PROP(NonmanifoldFacetHex3DA
                   "${test_dir}/nonmanifold_facetHex3D.A.off" 3 3 16 3)
SET_BBOX_PROP(NonmanifoldFacetHex4DA "(1,1,1)" "(2,4,3)" )
SET_NON_MANIFOLD_HEX(NonmanifoldFacetHex3DA 0 0 1 0 0)

SET_HEXAHEDRA_PROP(NonmanifoldFacetHex4DA
                   "${test_dir}/nonmanifold_facetHex4D.A.off" 4 3 22 5)
SET_BBOX_PROP(NonmanifoldFacetHex4DA "(0.5,0.5,0.5,0.25)" "(2.5,2.5,1.75.1.5)")
SET_NON_MANIFOLD_HEX(NonmanifoldFacetHex4DA 0 0 1 0 0)
SET_BOUNDARY(NonmanifoldFacetHex4DA 5 0)

SET_HEXAHEDRA_PROP(IncorrectOrientHexA
                     "${test_dir}/incorrect_orientHex.A.off" 3 3 16 3)
SET_HEXAHEDRA_PROP(IncorrectOrientHexB
                     "${test_dir}/incorrect_orientHex.B.off" 3 3 14 3)

SET_HEXAHEDRA_PROP(IrregularHexFacetIntersectionA
        "${test_dir}/irregular_hex_facet_intersection.A.off"
        3 3 13 2)        
SET_HEXAHEDRA_PROP(IrregularHexFacetIntersectionB
        "${test_dir}/irregular_hex_facet_intersection.B.off"
        3 3 16 4)

SET(cube_data_list CubeT CubeQ)
SET(polygon_data_list PolygonA)
SET(dup_data_list DupVertTA DupVertQA DupPolyQA)
SET(nonmanifold_data_list 
    NonmanifoldEdgeQA NonmanifoldEdgeQB
    NonmanifoldTwocubeTA NonmanifoldTwocubeQA
    NonmanifoldVertPoly2DA
    NonmanifoldFacetTetA NonmanifoldVertTetA NonmanifoldVertTetB
    NonmanifoldFacetSimp4DA NonmanifoldFacetSimp5DA
    NonmanifoldVertSimp4DA NonmanifoldVertSimp5DA
    NonmanifoldVertHex3DA NonmanifoldVertHex4DA 
    NonmanifoldEdgeHex3DA NonmanifoldEdgeHex4DA NonmanifoldEdgeHex4DB
    NonmanifoldFacetHex3DA NonmanifoldFacetHex4DA
    IrregularHexFacetIntersectionA IrregularHexFacetIntersectionB)
SET(boundary_data_list
    BoundTA BoundTB BoundTC BoundQA BoundQB BoundQC)
SET(selfI_no_grid_data_list SelfITA SelfITB SelfITC SelfITD SelfITE)
SET(selfI_data_list ${selfI_no_grid_data_list} Fuel35 Neghip20)
SET(sphere_data_list SphereCapTA SphereCapPA)
SET(tet_data_list CubeVolTetA BallOctantTetA SlabTetA)
SET(simplices_data_list Simplex4D Simp4DA Simplex5D Simp5DA)
SET(hex3D_data_list SphereHexA CubeVolHexA CubeVolHexB)
SET(hex_data_list ${hex3D_data_list} Hexahedra4DA)
SET(incorrect_orient_list IncorrectOrientHexA IncorrectOrientHexB
                          IncorrectOrientQA 
                          IncorrectOrientTetA IncorrectOrientTetB
                          IncorrectOrientSimp4DA IncorrectOrientSimp5DA)
SET(misc_data_list InternalQA)

SET(testdata_list ${cube_data_list} ${polygon_data_list} 
                  ${dup_data_list} 
                  ${nonmanifold_data_list} ${boundary_data_list}
                  ${selfI_data_list} ${sphere_data_list} 
                  ${tet_data_list} ${simplices_data_list} ${hex_data_list} 
                  ${misc_data_list})

SET(internal_test_list InternalQA BallOctantTetA SlabTetA)

SET(manifold_test_list ${cube_data_list} ${nonmanifold_data_list}
                       "SelfITA" ${dup_data_list} ${tet_data_list}
                       ${hex_data_list} ${boundary_data_list})

SET(orientation_test_list ${cube_data_list} 
                          ${incorrect_orient_list} ${nonmanifold_data_list}
                          ${hex_data_list})

SET(plot_angles_test_list ${cube_data_list} ${sphere_data_list})
SET(plot_edge_lengths_test_list ${hex3D_data_list})
SET(plot_jacobian_test_list ${hex3D_data_list})
SET(plot_jshape_test_list ${hex3D_data_list})


#---------------------------------------------------------
# ADD TESTS
#---------------------------------------------------------

IF (NOT no_re_test)

ADD_GENERAL_INFO_TEST_ON_DATA(${prog0} "" "info" "${testdata_list}")

# Manifold test
ADD_MANIFOLD_TEST_ON_DATA(${prog0} "" manifold "${manifold_test_list}")

ENDIF (NOT no_re_test)


#---------------------------------------------------------
# ADD COMPARISON TESTS
#---------------------------------------------------------

IF (NOT no_compare)
 
  ADD_COMPARISON_TEST_ON_DATA("info_compare" "" "${testdata_list}")
  ADD_COMPARISON_TEST_ON_DATA("reverse_orient_compare" "-reverse_orient" 
                              "${hex_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("for_each_compare" "-for_each_type" 
                              "${polygon_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("internal_compare" "-internal" 
                              "${internal_test_list}")
  ADD_COMPARISON_TEST_ON_DATA("internal_vert_compare" "-internal_vert" 
                              "${hex_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("manifold_compare" 
                              "-manifold" "${manifold_test_list}")
  ADD_COMPARISON_TEST_ON_DATA("manifold_terse_compare" 
                              "-manifold -terse" "${manifold_test_list}")
  ADD_COMPARISON_TEST_ON_DATA("manifold_deep_compare" 
                              "-manifold -report_deep" "${boundary_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("manifold_deep_terse_compare" 
                              "-manifold -report_deep -terse" 
                              "${boundary_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("orient_compare" 
                              "-oriented_manifold" "${orientation_test_list}")

  ADD_COMPARISON_TEST_ON_DATA("selfI_compare" "-selfI" "${selfI_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("selfI_terse_compare" "-selfI -terse" "${selfI_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("selfI_no_grid_compare" "-selfI_no_grid" "${selfI_no_grid_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("min_max_compare" "-out_min_angle -out_max_angle"
                              "${cube_data_list};${sphere_data_list};${tet_data_list};${hex_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("min_max_list_compare" 
                               "-out_min_angle -out_max_angle -plist"
                              "${cube_data_list};${sphere_data_list};${tet_data_list};${hex_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("min_jdet_compare" "-out_min_jdet -vlist"
                              "${hex3D_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("vlist_min" "-vlist_min"
                              "${hex3D_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("vlist_min_internal_vert" 
                              "-vlist_min -internal_vert"
                              "${hex3D_data_list}")

  ADD_COMPARISON_TEST_ON_DATA("plist_compare" "-plist" 
                              "${sphere_data_list};${tet_data_list};${hex_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("angleLE40_compare" "-angle_le 40" 
                              "${sphere_data_list};${tet_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("angleGE80_compare" " -angle_ge 80" 
                              "${sphere_data_list};${tet_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("angleLE40_plist_compare" "-plist -angle_le 40" 
                              "${sphere_data_list};${tet_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("angleGE80_plist_compare" "-plist -angle_ge 80" 
                              "${sphere_data_list};${tet_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("angleLE40_internal_compare" 
                              "-internal -angle_le 40" 
                              "${sphere_data_list};${tet_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("angleGE80_internal_compare" 
                              "-internal -angle_ge 80" 
                              "${sphere_data_list};${tet_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("plist_max_out_2_compare" 
                              "-plist -max_out 2" "${cube_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("plist_max_out_10_compare" 
                              "-plist -max_out 10" "${cube_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("vlist_compare" "-vlist"
                              "${cube_data_list};${dup_data_list};${sphere_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("vertex3_compare" "-vertex 3"
                              "${cube_data_list};${sphere_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("containsv5_compare" "-containsv 5"
                              "${cube_data_list};${sphere_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("containse2_5_compare" "-containse 2 5" "CubeT")
  ADD_COMPARISON_TEST_ON_DATA("containse2_6_compare" "-containse 2 6" "CubeT")
  ADD_COMPARISON_TEST_ON_DATA("containse9_25_compare" 
                              "-containse 9 25" "SphereCapTA")

  ADD_COMPARISON_TEST_ON_DATA("elist_compare" "-elist" "${hex_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("einfo_compare" "-edge_info" "${hex_data_list}")

  ADD_COMPARISON_TEST_ON_DATA("vJacobian" "-vJacobian" "${hex_data_list}")
  ADD_COMPARISON_TEST_ON_DATA("pAndvJacobian" "-pJacobian -vJacobian" 
                                              "${hex_data_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_angles" "-plot_angles" "${min_max_angle_suffixes}"
    "${plot_angles_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_edge_lengths" "-plot_edge_lengths" "${edge_length_suffixes}"
    "${plot_edge_lengths_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_edge_lengths_internalE" "-plot_edge_lengths -internal_edge" 
    "${edge_length_internalE_suffixes}" "${plot_edge_lengths_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_jacobian" "-plot_jacobian" "${min_max_jacobian_suffixes}"
    "${plot_jacobian_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_jdet_internalP" "-plot_jacobian -internal" 
    "${min_max_jdet_internalP_suffixes}" "${plot_jacobian_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_jdet_internalV" "-plot_jacobian -internal_vert" 
    "${min_max_jdet_internalV_suffixes}" "${plot_jacobian_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_vJacobian" "-plot_jacobian -vJacobian" 
    "${min_max_jdet_vJacobian_suffixes}" "${plot_jacobian_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_jshape" "-plot_jshape" 
    "${min_max_jshape_suffixes}" "${plot_jshape_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_jshape_internalV" "-plot_jshape -internal_vert" 
    "${min_max_jshape_internalV_suffixes}" "${plot_jshape_test_list}")

  ADD_COMPARE_OUTPUT_FILES_TEST_ON_DATA(
    "plot_jshape_vJacobian" "-plot_jshape -vJacobian" 
    "${min_max_jshape_vJacobian_suffixes}" "${plot_jshape_test_list}")

ENDIF(NOT no_compare)
