/*!
 *  @file ijkgrid_macros.h
 *  @brief ijk macros for grid classes.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2022 Rephael Wenger

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

#ifndef _IJKGRID_MACROS_H_
#define _IJKGRID_MACROS_H_

// ***********************************************************************
//! @name GENERAL MACROS
// ***********************************************************************

///@{

/// Macro for defining a local variable. <br>
/// Uses "for loop" to create local variable.
#define IJK_LOCAL(_init_statement,_counter0,_counter1)          \
  for (int _counter0=0, _counter1=0; _counter0<1; _counter0++)  \
    for (_init_statement; _counter1<1; _counter1++)

///@}


// ***********************************************************************
//! @name GRID CUBE MACROS
// ***********************************************************************

///@{

/// Local version of macro.
#define IJK_FOR_EACH_GRID_CUBE_LOCAL(_iv,_grid,_VTYPE,_cube_list,_i,_endv) \
  if (_grid.Dimension() > 0 && _grid.AxisSize(0) > 1)                \
    IJK_LOCAL(IJK::FACET0_CUBE_LIST<_VTYPE> _cube_list(_grid),       \
              _k0 ## __LINE__, _k1 ## __LINE__)                      \
      for (_VTYPE _i = 0; _i < _cube_list.NumCubes(); _i++)          \
        for (_VTYPE _iv = _cube_list.CubeIndex(_i),                  \
               _endv = _iv + _grid.AxisSize(0)-1;                    \
             _iv < _endv; _iv++)


/// For each grid cube.
#define IJK_FOR_EACH_GRID_CUBE(_iv,_grid,_VTYPE)                     \
  IJK_FOR_EACH_GRID_CUBE_LOCAL(_iv,_grid,_VTYPE,                     \
                               _cube_list ## __LINE__,               \
                               _i ## __LINE__,                       \
                               _endv ## __LINE__)

///@}


// ***********************************************************************
//! @name GRID FACET MACROS
// ***********************************************************************

///@{

/// Macro for defining local GRID_FACET_LIST.<br>
/// Uses "for loop" to create local variable.
#define IJK_LOCAL_GRID_FACET_LIST(_orth_dir,_grid,_DTYPE,_VTYPE,_facet_list,_counter0,counter1) \
   for (int _counter0=0, _counter1=0; _counter0<1; _counter0++)  \
     for (IJK::GRID_FACET_LIST<_DTYPE,_VTYPE> _facet_list(_grid,_orth_dir); _counter1<1; _counter1++)

/// Local version of macro.
#define IJK_FOR_EACH_GRID_FACET_ORTHOGONAL_TO_LOCAL(_iv,_orth_dir,_grid,_VTYPE,_facet_list,_i,_endv) \
  if (_grid.Dimension() > 0 && _grid.AxisSize(_orth_dir) > 0)        \
    IJK_LOCAL_GRID_FACET_LIST(_orth_dir,_grid,_VTYPE,_VTYPE,         \
                              _facet_list,                           \
                              _k0 ## __LINE__, _k1 ## __LINE__)      \
      for (_VTYPE _i = 0,                                            \
             _axis_inc = _grid.AxisIncrement(_orth_dir),             \
             _axis_size = _grid.AxisSize(_orth_dir);                 \
           _i < _facet_list.NumFacets(); _i++)                       \
        for (_VTYPE _iv = _facet_list.FacetIndex(_i),                \
               _endv = _iv + _axis_size*_axis_inc;                   \
             _iv < _endv; _iv += _axis_inc)

/// For each grid facet orthogonal to direction _orth_dir.
#define IJK_FOR_EACH_GRID_FACET_ORTHOGONAL_TO(_iv,_orth_dir,_grid,_VTYPE) \
  IJK_FOR_EACH_GRID_FACET_ORTHOGONAL_TO_LOCAL(_iv,_orth_dir,              \
                                              _grid,_VTYPE,               \
                                              _facet_list ## __LINE__,    \
                                              _i ## __LINE__,             \
                                              _endv ## __LINE__)

/// For each grid facet.
#define IJK_FOR_EACH_GRID_FACET(_iv,_orth_dir,_grid,_VTYPE)               \
  for (_VTYPE _orth_dir = 0; _orth_dir < _grid.Dimension(); _orth_dir++)  \
    IJK_FOR_EACH_GRID_FACET_ORTHOGONAL_TO(_iv,_orth_dir,_grid,_VTYPE)

/// Local version of macro.
#define IJK_FOR_EACH_INTERIOR_GRID_FACET_ORTHOGONAL_TO_LOCAL(_iv,_orth_dir,_grid,_VTYPE,_facet_list,_i,_endv) \
  if (_grid.Dimension() > 0 && _grid.AxisSize(_orth_dir) > 2)        \
    IJK_LOCAL_GRID_FACET_LIST(_orth_dir,_grid,_VTYPE,_VTYPE,         \
                              _facet_list,                           \
                              _k0 ## __LINE__, _k1 ## __LINE__)      \
      for (_VTYPE _i = 0,                                            \
             _axis_inc = _grid.AxisIncrement(_orth_dir),             \
             _axis_size = _grid.AxisSize(_orth_dir);                 \
           _i < _facet_list.NumFacets(); _i++)                       \
        for (_VTYPE _iv = _facet_list.FacetIndex(_i)+_axis_inc,      \
               _endv = _facet_list.FacetIndex(_i)+(_axis_size-1)*_axis_inc; \
             _iv < _endv; _iv += _axis_inc)

/// For each interior grid facet orthogonal to direction _orth_dir.
#define IJK_FOR_EACH_INTERIOR_GRID_FACET_ORTHOGONAL_TO(_iv,_orth_dir,_grid,_VTYPE) \
  IJK_FOR_EACH_INTERIOR_GRID_FACET_ORTHOGONAL_TO_LOCAL(_iv,_orth_dir,     \
                                              _grid,_VTYPE,               \
                                              _facet_list ## __LINE__,    \
                                              _i ## __LINE__,             \
                                              _endv ## __LINE__)

/// For each interior grid facet.
#define IJK_FOR_EACH_INTERIOR_GRID_FACET(_iv,_orth_dir,_grid,_VTYPE)      \
  for (_VTYPE _orth_dir = 0; _orth_dir < _grid.Dimension(); _orth_dir++)  \
    IJK_FOR_EACH_INTERIOR_GRID_FACET_ORTHOGONAL_TO(_iv,_orth_dir,_grid,_VTYPE)

/// Local version of macro.
#define IJK_FOR_EACH_BOUNDARY_GRID_FACET_ORTHOGONAL_TO_LOCAL(_iv,_orth_dir,_grid,_VTYPE,_facet_list,_i,_endv) \
  if (_grid.Dimension() > 0 && _grid.AxisSize(_orth_dir) > 0)        \
    IJK_LOCAL_GRID_FACET_LIST(_orth_dir,_grid,_VTYPE,_VTYPE,         \
                              _facet_list,                           \
                              _k0 ## __LINE__, _k1 ## __LINE__)      \
      for (_VTYPE _i = 0,                                            \
             _axis_inc = _grid.AxisIncrement(_orth_dir),             \
             _axis_size = _grid.AxisSize(_orth_dir);                 \
           _i < _facet_list.NumFacets(); _i++)                       \
        for (_VTYPE _iv = _facet_list.FacetIndex(_i),                \
               _endv = _iv + _axis_size*_axis_inc;                   \
             _iv < _endv;                                            \
             _iv += std::max((_axis_size-1),_VTYPE(1))*_axis_inc)

/// For each boundary grid facet orthogonal to direction _orth_dir.
#define IJK_FOR_EACH_BOUNDARY_GRID_FACET_ORTHOGONAL_TO(_iv,_orth_dir,_grid,_VTYPE) \
  IJK_FOR_EACH_BOUNDARY_GRID_FACET_ORTHOGONAL_TO_LOCAL(_iv,_orth_dir,     \
                                              _grid,_VTYPE,               \
                                              _facet_list ## __LINE__,    \
                                              _i ## __LINE__,             \
                                              _endv ## __LINE__)

/// For each boundary grid facet.
#define IJK_FOR_EACH_BOUNDARY_GRID_FACET(_iv,_orth_dir,_grid,_VTYPE)      \
  for (_VTYPE _orth_dir = 0; _orth_dir < _grid.Dimension(); _orth_dir++)  \
    IJK_FOR_EACH_BOUNDARY_GRID_FACET_ORTHOGONAL_TO(_iv,_orth_dir,_grid,_VTYPE)

///@}


// ***********************************************************************
//! @name GRID EDGE MACROS
// ***********************************************************************

///@{

/// Local version of macro.
#define IJK_FOR_EACH_GRID_EDGE_LOCAL(_iend0,_edge_dir,_grid,_VTYPE,_vlist,_i,_axis_inc,_axis_size,_endv) \
  if (_grid.Dimension() > 0)                                        \
    IJK_LOCAL(IJK::FACET_VERTEX_LIST<_VTYPE>               \
              _vlist(_grid, 0, true),                               \
              _k0 ## __LINE__, _k1 ## __LINE__)                     \
      for (_VTYPE _edge_dir = 0; _edge_dir < _grid.Dimension();     \
           _edge_dir++)                                             \
        if (_grid.AxisSize(_edge_dir) > 0)                          \
          IJK_LOCAL(_vlist.GetVertices(_grid, _edge_dir),           \
                   _k2 ## __LINE__, _k3 ## __LINE__)                \
            for (_VTYPE _i = 0,                                     \
                   _axis_inc = _grid.AxisIncrement(_edge_dir),      \
                   _axis_size = _grid.AxisSize(_edge_dir);          \
                 _i < _vlist.NumVertices(); _i++)                   \
              for (_VTYPE _iend0 = _vlist.VertexIndex(_i),          \
                     _endv = _iend0 + (_axis_size-1)*_axis_inc;     \
                   _iend0 < _endv;                                  \
                   _iend0 +=  _axis_inc)

/// For each grid edge.
#define IJK_FOR_EACH_GRID_EDGE(_iend0,_edge_dir,_grid,_VTYPE) \
  IJK_FOR_EACH_GRID_EDGE_LOCAL(_iend0,_edge_dir,_grid,_VTYPE, \
                                        _vlist ## __LINE__,         \
                                        _i ## __LINE__,             \
                                        _axis_inc ## __LINE__,      \
                                        _axis_size ## __LINE__,     \
                                        _endv ## __LINE__)

/// Local version of macro.
#define IJK_FOR_EACH_INTERIOR_GRID_EDGE_IN_DIRECTION_LOCAL(_iend0,_edge_dir,_grid,_VTYPE,_vlist,_i,_axis_inc,_axis_size,_endv) \
  if (_grid.AxisSize(_edge_dir) > 0)                                \
    IJK_LOCAL(IJK::FACET_INTERIOR_VERTEX_LIST<_VTYPE>               \
              _vlist(_grid, _edge_dir, false, true),                \
            _k0 ## __LINE__, _k1 ## __LINE__)                       \
      for (_VTYPE _i = 0,                                           \
             _axis_inc = _grid.AxisIncrement(_edge_dir),            \
             _axis_size = _grid.AxisSize(_edge_dir);                \
           _i < _vlist.NumVertices(); _i++)                         \
        for (_VTYPE _iend0 = _vlist.VertexIndex(_i),                \
               _endv = _iend0 + (_axis_size-1)*_axis_inc;           \
             _iend0 < _endv;                                        \
             _iend0 +=  _axis_inc)

/// For each grid edge in direction _edge_dir.
/// - If grid dimension is 1, this processes all grid edges.
#define IJK_FOR_EACH_INTERIOR_GRID_EDGE_IN_DIRECTION(_iend0,_edge_dir,_grid,_VTYPE) \
       IJK_FOR_EACH_INTERIOR_GRID_EDGE_IN_DIRECTION_LOCAL           \
       (_iend0,_edge_dir,_grid,_VTYPE,                              \
        _vlist ## __LINE__,                                         \
        _i ## __LINE__,                                             \
        _axis_inc ## __LINE__,                                      \
        _axis_size ## __LINE__,                                     \
        _endv ## __LINE__)

/// For each grid edge which does not lie on any grid boundary facet.
#define IJK_FOR_EACH_INTERIOR_GRID_EDGE(_iend0,_edge_dir,_grid,_VTYPE) \
  if (_grid.Dimension() > 0)                                        \
    for (_VTYPE _edge_dir = 0; _edge_dir < _grid.Dimension();       \
         _edge_dir++)                                               \
      IJK_FOR_EACH_INTERIOR_GRID_EDGE_IN_DIRECTION                  \
        (_iend0,_edge_dir,_grid,_VTYPE)

///@}


// ***********************************************************************
//! @name GRID VERTEX MACROS
// ***********************************************************************

///@{

/// Local version of macro.
#define IJK_FOR_EACH_INTERIOR_GRID_VERTEX_LOCAL(_iv,_grid,_VTYPE,_vlist,_i,_endv) \
  if (_grid.Dimension() > 0)                                        \
    if (_grid.AxisSize(0) > 1)                                      \
      IJK_LOCAL(IJK::FACET_INTERIOR_VERTEX_LIST<_VTYPE>             \
                _vlist(_grid, 0, false, true),                      \
                _k0 ## __LINE__, _k1 ## __LINE__)                   \
        for (_VTYPE _i = 0;                                         \
             _i < _vlist.NumVertices(); _i++)                       \
          for (_VTYPE _iv = _vlist.VertexIndex(_i)+1,               \
                 _endv = _iv + (_grid.AxisSize(0)-2);               \
               _iv < _endv;                                         \
               _iv++)

/// For each grid vertex in the interior of the grid,
/// i.e. which does not lie on a boundary grid facet.
/// - Processes 0 vertices if grid dimension is 0.
#define IJK_FOR_EACH_INTERIOR_GRID_VERTEX(_iv,_grid,_VTYPE)         \
  IJK_FOR_EACH_INTERIOR_GRID_VERTEX_LOCAL                           \
  (_iv,_grid,_VTYPE, _vertex_list ## __LINE__, _i ## __LINE,        \
   _endv ## __LINE__)

/// Local version of macro.
#define IJK_FOR_EACH_BOUNDARY_GRID_VERTEX_LOCAL(_iv,_grid,_VTYPE,_vertex_list,_i) \
  IJK_LOCAL(IJK::GRID_BOUNDARY_VERTEX_LIST<_VTYPE> _vertex_list(_grid), \
            _k0 ## __LINE__, _k1 ## __LINE__)                       \
    for (_VTYPE _i = 0, _iv = _vertex_list.VertexIndex(_i);         \
         _i < _vertex_list.NumVertices();                           \
         _i++, _iv = _vertex_list.VertexIndex(_i))

/// For each boundary grid vertex, i.e., each grid vertex which
/// lies on the a boundary grid facet.
#define IJK_FOR_EACH_BOUNDARY_GRID_VERTEX(_iv,_grid,_VTYPE)         \
  IJK_FOR_EACH_BOUNDARY_GRID_VERTEX_LOCAL(_iv,_grid,_VTYPE,         \
                                          _vertex_list ## __LINE__, \
                                          _i ## __LINE__)

/// Local version of macro.
#define IJK_FOR_EACH_VERTEX_IN_GRID_FACET_INTERIOR_LOCAL(_iv,_orth_dir,_grid,_VTYPE,_vlist,_i) \
  IJK_LOCAL(IJK::FACET_INTERIOR_VERTEX_LIST<_VTYPE>                 \
            _vlist(_grid, _orth_dir, false),                        \
            _k0 ## __LINE__, _k1 ## __LINE__)                       \
    if (_vlist.NumVertices() > 0)                                   \
      for (_VTYPE _i = 0, _iv = _vlist.VertexIndex(_i);       \
           _i < _vlist.NumVertices();                               \
           _i++, _iv = _vlist.VertexIndex(_i))

// Direction orthogonal to facet (_orth_dir) is an input parameter.
// Computes vertices for grid facets incident on origin (0,0,...,0).
/// For each grid vertex on the inteior of the grid boundary facet
/// orthogonal to direction _orth_dir 
/// and incident on the origin (0,0,...,0).
#define IJK_FOR_EACH_VERTEX_IN_GRID_FACET_INTERIOR(_iv,_orth_dir,_grid,_VTYPE) \
  IJK_FOR_EACH_VERTEX_IN_GRID_FACET_INTERIOR_LOCAL                  \
  (_iv,_orth_dir,_grid,_VTYPE,_vertex_list ## __LINE__,_i ## __LINE__)

///@}

#endif

