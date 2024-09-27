/*!
 *  @file ijkgenscalar.cxx
 *  @brief Generate a scalar field. (Store in .nrrd format.)
 *  - Version v0.3.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2011-2023 Rephael Wenger

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

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "ijkNrrd.h"
#include "ijkgenscalar.h"
#include "ijkgenscalarIO.h"

#include "ijkcube.tpp"
#include "ijkcommand_line.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkgradientfield.tpp"
#include "ijkscalarfield.tpp"
#include "ijkstring.tpp"

// Types
typedef IJK::NRRD_DATA<int, int> NRRD_HEADER;

using namespace IJKSCALARFIELD;
using namespace IJKGENSCALAR;
using namespace IJKGENGEOM;
using namespace std;

// global constants
const int DIM2(2);
const int DIM3(3);
const int DIM4(4);

// global variables
char * ofilename = NULL;
std::vector<AXIS_SIZE_TYPE> grid_axis_size;
std::vector<FIELD_INFO> field_info;
FIELD_PARAM field_param;
bool output_gradients = false;
bool flag_silent = false;
bool flag_gzip = false;
double max_small_magnitude(0.0001);

// functions to read field parameters
void init_field_info(vector<FIELD_INFO> & field_info);
void read_dimension();
void read_num_vert_per_axis();
void read_field_name(const vector<FIELD_INFO> & field_info);
void prompt_param(const std::string & field_name);
void set_field_param(const vector<FIELD_INFO> & field_info);

// global functions
void intersect_with_wedge
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid);
void intersect_with_wedge
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & grid, GRADIENT_GRID & gradient);
void generate_field(SCALAR_GRID & scalar_grid);
void generate_gradient_field
(SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid);
void set_random_param
(const int dimension, const int num_objects);

// misc functions
void construct_gradient_filename
(const char * scalar_filename, string & gradient_filename);
void parse_command_line(const int argc, char **argv);
void usage_error(), help_msg();
void check_input_param(const int dimension);
void check_param(const SCALAR_GRID & grid, const std::string & field_name);
void check_field();
void check_field_options();
bool check_dimension(const int dimension, const int dimension2,
                     const char * text, IJK::ERROR & error);

int main(int argc, char **argv)
{
  IJK::PROCEDURE_ERROR error("ijkgenscalar");

  try {

    init_field_info(field_info);

    parse_command_line(argc, argv);

    if (ofilename == NULL) {
      error.AddMessage("Missing output file name");
      throw error;
    }

    if (!field_param.geom_info_index.IsSet()) 
      { read_field_name(field_info); }

    // Set parameters related to field.
    set_field_param(field_info);
    check_field();

    int ifield = field_param.FieldIndex();
    std::string field_name = field_info[ifield].name;

    if (!field_param.dimension.IsSet()) 
      { read_dimension(); }

    check_field_options();
    check_input_param(field_param.Dimension());

    if (field_name == "isotable_entry") {
      // Fix axis size to 2. (Grid is a single cube.)
      grid_axis_size.clear();
      grid_axis_size.push_back(2);
    }
    else if (field_name == "set_xcube") {
      // Fix axis size to 4. (Grid is 4x4x4x...)
      grid_axis_size.clear();
      grid_axis_size.push_back(4);
    }
    else if (field_name == "set_x2cube") {
      // Set grid axis size to 4x4x...x5. (Two interior cubes.)
      grid_axis_size.clear();

      for (int d = 0; d+1 < field_param.Dimension(); d++) 
        { grid_axis_size.push_back(4); }
      grid_axis_size.push_back(5);
    }
    
    if (grid_axis_size.size() == 0)
      { read_num_vert_per_axis(); }

    if (grid_axis_size.size() == field_param.Dimension()) 
      { field_param.SetAxisSize(grid_axis_size); }
    else 
      { field_param.SetAxisSize(grid_axis_size[0]); }

    SCALAR_GRID scalar_grid
      (field_param.Dimension(), &(field_param.axis_size[0]));

    if (field_info[ifield].name == "dist2planes")
      { field_param.flag_multi_normals = true; }

    if (field_param.spacing.size() >= scalar_grid.Dimension())
      { scalar_grid.SetSpacing(field_param.SpacingPtrConst()); };

    if (field_param.flag_grid_center && field_param.center.size() == 0) {
      for (int d = 0; d < field_param.Dimension(); d++) {
        COORD_TYPE x = (field_param.axis_size[d]-1)/2.0;
        x *= scalar_grid.Spacing(d);
        field_param.center.push_back(x);
      }
    }

    // Prompt for any missing parameters.
    prompt_param(field_name);

    check_param(scalar_grid, field_name);

    int dimension = field_param.Dimension();
    if (output_gradients) {
      GRADIENT_GRID gradient_grid;

      gradient_grid.SetSize
        (dimension, scalar_grid.AxisSize(), dimension);

      if (!field_info[ifield].IsGradientImplemented()) {
        error.AddMessage
          ("Programming error. Gradient not implemented for function ",
           field_info[ifield].name, ".");
        throw error;
      }

      if (field_param.spacing.size() >= gradient_grid.Dimension()) {
        { gradient_grid.SetSpacing(field_param.SpacingPtrConst()); };
      }

      generate_gradient_field(scalar_grid, gradient_grid);

      if (!flag_silent) {
        cout << "Writing scalar field to " << ofilename << endl;
      }
      write_scalar_grid
        (ofilename, scalar_grid, field_param, field_info, flag_gzip);

      string gradient_filename;
      construct_gradient_filename(ofilename, gradient_filename);
      if (!flag_silent) {
        cout << "Writing gradient field to " << gradient_filename << endl;
      }
      write_gradient_grid
        (gradient_filename, gradient_grid, field_param, field_info, flag_gzip);
    }
    else {

      generate_field(scalar_grid);

      if (field_param.flag_set_boundary2zero) {
        set_boundary_values(0, scalar_grid);
      }

      if (!flag_silent) {
        cout << "Writing scalar field to " << ofilename << endl;
      }
      write_scalar_grid
        (ofilename, scalar_grid, field_param, field_info, flag_gzip);
    }

  } 
  catch (IJK::ERROR & error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}


// **************************************************
// Compute functions
// **************************************************

void compute_center_translation_vector
(const int dimension, COORD_TYPE translate[])
{
  const int ifield = field_param.FieldIndex();
  const std::string field_name = field_info[ifield].name;
  IJK::PROCEDURE_ERROR error("compute_center_translation_vector");

  if (field_param.flag_tilt) {

    if (field_param.NumDirections() <= 0) {
      error.AddMessage("Programming error. Missing tilt/axis direction.");
      throw error;
    };

    if (field_name == "cube") {
      IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);

      if (field_param.NumSideDirections() <= 0) {
        error.AddMessage("Programming error. Missing side direction.");
        throw error;
      };

      compute_orthogonal_basis3D
        (field_param.DirectionPtrConst(0), field_param.SideDirectionPtrConst(0),
         zdir.Ptr(), xdir.Ptr(), ydir.Ptr());
      transform_basis3D
        (field_param.TranslatePtrConst(0), 
         xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), translate);
    }
    else {
      rotate_vector_ed
        (dimension, max_small_magnitude, field_param.DirectionPtrConst(0),
         field_param.TranslatePtrConst(0), translate);
    }
  }
  else {
    IJK::copy_coord
      (dimension, field_param.TranslatePtrConst(0), translate);
  }
}

template <typename T>
void compute_num_table_entries(const int dimension, T & num_table_entries)
{
  T num_cube_vertices;
  IJK::PROCEDURE_ERROR error("compute_num_table_entries");

  try {
    num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
    IJK::int_power(2, num_cube_vertices, num_table_entries, error);
  } catch (IJK::ERROR & error2) {
    cerr << "Input error. Dimension " << dimension
         << " is too large for field isotable_entry."
         << endl;
    cerr << "Isosurface lookup table has more than "
         << std::numeric_limits<ISOTABLE_INDEX_TYPE>::max()
         << " entries." << endl;
    exit(50);
  }

}


// **************************************************
// Prompt for scalar field parameters
// **************************************************

template <typename CTYPE>
void prompt_coordinates
(const char * prompt, const int num_coord,
 std::vector<CTYPE> & coord)
{
  cout << prompt << " (" << num_coord << " values) : ";
  for (int ic = 0; ic < num_coord; ic++) {
    CTYPE c;
    cin >> c;
    coord.push_back(c);
  }
}

template <typename CTYPE>
void prompt_coordinates
(const std::ostringstream & prompt, const int num_coord,
 std::vector<CTYPE> & coord)
{
  prompt_coordinates(prompt.str().c_str(), num_coord, coord);
}

/// Prompt for a single scalar value.
template <typename STYPE>
void prompt_scalar
(const char * prompt, STYPE & s)
{
  cout << prompt << " : ";
  cin >> s;
}

/// Prompt for a single scalar value and add to list
template <typename STYPE>
void prompt_scalar
(const char * prompt, std::vector<STYPE> & list)
{
  cout << prompt << " : ";
  STYPE s;
  cin >> s;
  list.push_back(s);
}

template <typename STYPE>
void prompt_scalar
(const std::ostringstream & prompt, std::vector<STYPE> & list)
{
  prompt_scalar(prompt.str().c_str(), list);
}

void prompt_num_objects(const int field_index)
{
  string field_name = field_info[field_index].name;

  int n(0);
  while (true) {

    if (field_name == "cube") {
      cout << "Enter number of cubes: ";
    }
    else if (field_name == "cylinder") {
      cout << "Enter number of cylinders: ";
    }
    else if (field_name == "annulus") {
      cout << "Enter number of annuli: ";
    }
    else if (field_name == "dist2planes") {
      cout << "Enter number of planes: ";
    }
    else {
      cout << "Enter number of objects: ";
    }

    cin >> n;
    if (n > 0) { break; };

    cout << "Illegal number.  Enter a positive integer." << endl;
  }

  field_param.num_objects.Set(n);
}

void prompt_center_translate(const int dimension)
{
  if (field_param.translate.size() == 0) {
    prompt_coordinates("Enter center translation vector",
                       dimension, field_param.translate);
  }
}

void prompt_location
(const int dimension, const int num_obj, 
 const char * obj_name, const char * loc_name)
{
  int k = field_param.center.size()/dimension;

  // Clear any stray coordinates
  if (field_param.center.size()%dimension != 0) 
    { field_param.center.resize(k*dimension); }

  if (k >= num_obj) { return; };

  if (num_obj == 1) {
    ostringstream prompt;
    prompt << "Enter coordinates of " << obj_name << " " << loc_name;
    prompt_coordinates(prompt, dimension, field_param.center);
  }
  else {

    if (field_param.flag_stack) {

      if (k < 1) {
        ostringstream prompt;
        prompt << "Enter coordinates of first " 
               << obj_name << " " << loc_name;
        prompt_coordinates(prompt, dimension, field_param.center);
      }

      prompt_center_translate(dimension);
    }
    else {
      for (int i = k; i < num_obj; i++) {
        if (field_param.center.size() <= i*dimension) {
          ostringstream prompt;
          prompt << "Enter coordinates of " << obj_name << " " 
                 << i << " " << loc_name;
          prompt_coordinates(prompt, dimension, field_param.center);
        }
      }
    }
  }
}

void prompt_center
(const int dimension, const int num_obj, const char * obj_name)
{
  prompt_location(dimension, num_obj, obj_name, "center");
}

void prompt_apex
(const int dimension, const int num_obj, const char * obj_name)
{
  prompt_location(dimension, num_obj, obj_name, "apex");
}

void prompt_two_axis_directions(const int dimension)
{
  if (field_param.direction.size() == 0) {
    prompt_coordinates("Enter first axis direction",
                       dimension, field_param.direction);
  }

  if (field_param.side_direction.size() == 0) {
    prompt_coordinates("Enter second axis direction",
                       dimension, field_param.side_direction);
  }
}

void prompt_wedge(const int dimension)
{
  const int field_index = field_param.FieldIndex();
  string field_name = field_info[field_index].name;

  if (field_param.wedge_angle.size() == 0) {
    prompt_scalar("Enter wedge angle", field_param.wedge_angle);
  }

  if (field_param.side_direction.size() == 0) {
    prompt_coordinates
      ("Enter coordinates of side direction",
       dimension, field_param.side_direction);
  }

  if (field_param.NumWedgeIsovalues() < 1) {
    prompt_scalar("Enter isovalue where wedge passes through the center", 
                  field_param.wedge_isovalue);
  }

}

void prompt_flange()
{
  if (field_param.flange_width.size() == 0) {
    prompt_scalar("Enter flange width", field_param.flange_width);
  }

  if (field_param.flange_height.size() == 0) {
    prompt_scalar("Enter flange height", field_param.flange_height);
  }
}

void prompt_cube
(const int dimension, const int num_cubes)
{
  prompt_center(dimension, num_cubes, "cube");

  if (field_param.flag_tilt) 
    { prompt_two_axis_directions(dimension); }
}

void prompt_cylinder
(const int dimension, const int num_cylinders, const bool flag_closed)
{
  prompt_center(dimension, num_cylinders, "cylinder");

  if (field_param.flag_tilt) {
    if (field_param.direction.size() == 0) {
      prompt_coordinates("Enter cylinder axis direction",
                         dimension, field_param.direction);
    }
  }
  else {
    // Set axis direction to (1,0,0,...)
    field_param.direction.clear();
    field_param.direction.resize(dimension, 0);
    if (dimension > 0)
      { field_param.direction[0] = 1; }
  }

  if (flag_closed) {
    if (field_param.length_difference.size() == 0) {
      prompt_scalar("Enter difference between cylinder length and diameter", 
                    field_param.length_difference);
    }
  }

}

void prompt_cone(const int dimension, const bool flag_closed)
{
  const int num_cones = 1;

  prompt_center(dimension, num_cones, "cone");

  if (field_param.flag_tilt) {
    if (field_param.direction.size() == 0) {
      prompt_coordinates("Enter cone axis direction",
                         dimension, field_param.direction);
    }
  }
  else {
    // Set axis direction to (1,0,0,...)
    field_param.direction.clear();
    field_param.direction.resize(dimension, 0);
    if (dimension > 0)
      { field_param.direction[0] = 1; }
  }

  if (field_param.angle.size() == 0) {
    prompt_scalar("Enter cone angle", field_param.angle);
  }

}

void prompt_frustrum(const int dimension)
{
  const int num_frustra = 1;

  prompt_apex(dimension, num_frustra, "frustrum");

  if (field_param.flag_tilt) {
    if (field_param.direction.size() == 0) {
      prompt_coordinates("Enter frustrum axis direction",
                         dimension, field_param.direction);
    }
  }
  else {
    // Set axis direction to (1,0,0,...)
    field_param.direction.clear();
    field_param.direction.resize(dimension, 0);
    if (dimension > 0)
      { field_param.direction[0] = 1; }
  }

  if (field_param.angle.size() == 0) {
    prompt_scalar("Enter frustrum angle", field_param.angle);
  }

  if (field_param.dist2near0.size() == 0) {
    prompt_scalar("Enter distance from apex to near plane", 
                  field_param.dist2near0);
  }

  if (field_param.dist2far0.size() == 0) {
    prompt_scalar("Enter distance from apex to far plane", 
                  field_param.dist2far0);
  }

}

void prompt_cannon(const int dimension)
{
  const int num_cannon = 1;

  prompt_apex(dimension, num_cannon, "cannon");

  if (field_param.flag_tilt) {
    if (field_param.direction.size() == 0) {
      prompt_coordinates("Enter cannon axis direction",
                         dimension, field_param.direction);
    }
  }
  else {
    // Set axis direction to (1,0,0,...)
    field_param.direction.clear();
    field_param.direction.resize(dimension, 0);
    if (dimension > 0)
      { field_param.direction[0] = 1; }
  }

  if (field_param.angle.size() == 0) {
    prompt_scalar("Enter cannon angle", field_param.angle);
  }

  if (field_param.dist2near0.size() == 0) {
    prompt_scalar("Enter distance from apex to near plane", 
                  field_param.dist2near0);
  }

  if (field_param.dist2ball_center.size() == 0) {
    prompt_scalar("Enter distance from apex to ball center", 
                  field_param.dist2ball_center);
  }
}

void prompt_octahedron(const int dimension, const int num_objects)
{
  prompt_center(dimension, num_objects, "octahedron");
}

void prompt_dist2planes(const int dimension, const int num_planes)
{
  if (field_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of point on all planes",
                       dimension, field_param.center);
  }

  int k = field_param.normal.size()/dimension;

  // Clear any stray coordinates
  if (field_param.normal.size()%dimension != 0) 
    { field_param.normal.resize(k*dimension); }

  if (k >= num_planes) { return; };

  for (int i = k; i < num_planes; i++) {
    ostringstream prompt;
    prompt << "Enter coordinates of normal to plane " << i;

    prompt_coordinates(prompt, dimension, field_param.normal);
  }
}

void prompt_edge(const int dimension)
{
  if (field_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of point on edge",
                       dimension, field_param.center);
  }

  if (field_param.angle.size() == 0) {
    prompt_scalar
      ("Enter dihedral angle", field_param.angle);
  }

  if (field_param.normal.size() == 0) {
    prompt_coordinates
      ("Enter coordinates of normal to first plane", 
       dimension, field_param.normal);
  }

  if (field_param.rotation_direction.size() == 0) {
    prompt_coordinates
      ("Enter rotation direction", dimension, field_param.rotation_direction);
  }
}

void prompt_square_cylinder(const int dimension)
{
  if (field_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of cylinder center",
                       dimension, field_param.center);
  }

  if (field_param.direction.size() == 0) {
    prompt_coordinates("Enter cylinder axis direction",
                       dimension, field_param.direction);
  }

  if (field_param.side_direction.size() == 0) {
    prompt_coordinates("Enter direction orthogonal to square side",
                       dimension, field_param.side_direction);
  }

  if (field_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between cylinder height and width", 
                  field_param.length_difference);
  }
}

void prompt_square_cylinder_rot45(const int dimension)
{
  if (field_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of cylinder center",
                       dimension, field_param.center);
  }

  if (field_param.direction.size() == 0) {
    prompt_coordinates("Enter cylinder axis direction",
                       dimension, field_param.direction);
  }

  if (field_param.side_direction.size() == 0) {
    prompt_coordinates("Enter square corner direction",
                       dimension, field_param.side_direction);
  }

  if (field_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between cylinder height and width", 
                  field_param.length_difference);
  }
}

void prompt_crossing_square_cylinders(const int dimension)
{
  if (field_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of crossing center",
                       dimension, field_param.center);
  }

  if (field_param.direction.size() == 0) {
    prompt_coordinates("Enter direction of first cylinder",
                       dimension, field_param.direction);
  }

  if (field_param.side_direction.size() == 0) {
    prompt_coordinates("Enter direction of second cylinder",
                       dimension, field_param.side_direction);
  }

  if (field_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between cylinder height and width", 
                  field_param.length_difference);
  }
}

void prompt_annulus(const int dimension)
{
  if (field_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of annulus center",
                       dimension, field_param.center);
  }

  if (field_param.direction.size() == 0) {
    prompt_coordinates("Enter annulus axis direction",
                       dimension, field_param.direction);
  }

  if (field_param.radius.size() == 0) {
    prompt_scalar
      ("Enter annulus radius", field_param.radius);
  }

  if (field_param.length_difference.size() == 0) {
    prompt_scalar("Enter difference between annulus height and width", 
                  field_param.length_difference);
  }

}

void prompt_torus(const int dimension, const int num_tori)
{
  prompt_center(dimension, num_tori, "torus");

  if (field_param.direction.size() == 0) {
    prompt_coordinates("Enter torus axis direction",
                       dimension, field_param.direction);
  }

  if (field_param.radius.size() == 0) {
    prompt_scalar("Enter torus radius", field_param.radius);
  }
}

void prompt_constant_scalar()
{
  if (field_param.scalar_value.size() == 0) {
    SCALAR_TYPE s;
    prompt_scalar("Enter scalar value", s);
    field_param.SetScalarValue(s);
  }
}


void prompt_constant_unit_gradient(const int dimension)
{
  if (field_param.center.size() == 0) {
    prompt_coordinates("Enter coordinates of point with scalar value zero",
                       dimension, field_param.center);
  }

  if (field_param.direction.size() == 0) {
    prompt_coordinates("Enter gradient direction",
                       dimension, field_param.direction);
  }
}

template <typename T>
void prompt_maxval(IJK::SET_VALUE<T> & maxval)
{
  T x;

  cout << "Enter maximum scalar value: ";
  cin >> x;
  maxval.Set(x);
}

template <typename T>
void prompt_random_base(IJK::SET_VALUE<T> & random_base)
{
  T x;

  cout << "Enter base for random number generation: ";
  cin >> x;
  random_base.Set(x);
}

template <typename T>
void prompt_max_exponent(IJK::SET_VALUE<T> & max_exponent)
{
  T x;

  cout << "Enter maximum exponent (Maxval = (base)^(max_exponent)): ";
  cin >> x;
  max_exponent.Set(x);
}

template <typename T>
void prompt_random_seed(IJK::SET_VALUE<T> & random_seed)
{
  T x;

  cout << "Enter random seed: ";
  cin >> x;
  random_seed.Set(x);
}

/// Prompt for generation of random integers.
void prompt_random_int()
{
  if (!field_param.maxval.IsSet()) 
    { prompt_maxval(field_param.maxval); }

  if (!field_param.random_seed.IsSet())
    { prompt_random_seed(field_param.random_seed); }
}

/// Prompt for generation of base raised to a random power.
void prompt_random_power()
{
  if (!field_param.random_base.IsSet()) 
    { prompt_random_base(field_param.random_base); }

  if (!field_param.max_exponent.IsSet()) 
    { prompt_max_exponent(field_param.max_exponent); }

  if (!field_param.random_seed.IsSet())
    { prompt_random_seed(field_param.random_seed); }

}

template <typename T>
void prompt_isotable_index
(const int dimension, IJK::SET_VALUE<T> & isotable_index)
{
  T num_table_entries;
  IJK::PROCEDURE_ERROR error("prompt_isotable_index");

  compute_num_table_entries(dimension, num_table_entries);

  if (num_table_entries < 1) {
    isotable_index.Set(0);
    return;
  }

  T x;
  do {
    cout << "Enter isosurface lookup table index [0," << num_table_entries-1
         << "]: ";
    cin >> x;

    if (x < 0) {
      cout << "Error.  Isosurface lookup table index must be non-negative."
           << endl;
    }
    else if (x >= num_table_entries) {
      cout << "Error.  Isosurface lookup table index must be less than "
           << num_table_entries << "." << endl;
    }

  } while (x < 0 || x >= num_table_entries);

  isotable_index.Set(x);
}


void prompt_set_scalar_field()
{
  if (!field_param.maxval.IsSet()) 
    { prompt_maxval(field_param.maxval); }

  if (field_param.set_scalar_list.size() == 0) {
    int n;
    prompt_scalar("Enter number of cubes to set", n);
    for (int i = 0; i < n; i++) {

      int cube_index;
      ISOTABLE_INDEX_TYPE isotable_index;
      SCALAR_TYPE s;
      cout << "Enter cube " << i << " index and isotable index: ";
      cin >> cube_index;
      cin >> s;
      field_param.SetCubeIsotableIndex(cube_index, isotable_index);
    }
  }

}


void prompt_set_xcube(const int dimension)
{
  const int num_cube_vertices = IJK::compute_num_cube_vertices(dimension);

  if (field_param.scalar_value.size() == 0) {

    cout << "Enter scalar values at cube vertices (" 
         << num_cube_vertices << " values) : ";

    for (int ic = 0; ic < num_cube_vertices; ic++) {
      COORD_TYPE c;
      cin >> c;
      field_param.scalar_value.push_back(c);
    }
  }

}


void prompt_set_x2cube(const int dimension)
{
  const int num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
  const int num_facet_vertices = num_cube_vertices/2;

  if (field_param.scalar_value.size() == 0) {

    cout << "Enter scalar values at two adjacent cube vertices (" 
         << num_cube_vertices+num_facet_vertices << " values) : ";

    for (int ic = 0; ic < num_cube_vertices+num_facet_vertices; ic++) {
      COORD_TYPE c;
      cin >> c;
      field_param.scalar_value.push_back(c);
    }
  }

}


/// Prompt for generation of isosurface lookup table entry.
void prompt_isotable_index(const int dimension)
{
  if (!field_param.isotable_index.IsSet())
    { prompt_isotable_index(dimension, field_param.isotable_index); }
}

/// Prompt for any missing parameters.
void prompt_param(const std::string & field_name)
{
  const int dimension = field_param.Dimension();

  if (!field_param.num_objects.IsSet()) {
    if (field_param.flag_multi_centers || field_param.flag_multi_normals) {
      prompt_num_objects(field_param.FieldIndex());
    }
    else {
      field_param.num_objects.Set(1); 
    }
  }

  const int num_objects = field_param.NumObjects();

  set_random_param(dimension, num_objects);

  if (field_name == "cube") {
    prompt_cube(dimension, num_objects);
  }
  else if (field_name == "sphere") {
    prompt_center(dimension, num_objects, "sphere");
  }
  else if (field_name == "cylinder") {
    prompt_cylinder(dimension, num_objects, true);
  }
  else if (field_name == "open_cylinder") {
    prompt_cylinder(dimension, num_objects, false);
  }
  else if (field_name == "cone") {
    prompt_cone(dimension, true);
  }
  else if (field_name == "frustrum") {
    prompt_frustrum(dimension);
  }
  else if (field_name == "cannon") {
    prompt_cannon(dimension);
  }
  else if (field_name == "square") {
    prompt_square_cylinder(dimension);
  }
  else if (field_name == "square45") {
    prompt_square_cylinder_rot45(dimension);
  }
  else if (field_name == "squareX" ||
           field_name == "square45X") {
    prompt_crossing_square_cylinders(dimension);
  }
  else if (field_name == "annulus") {
    prompt_annulus(dimension);
  }
  else if (field_name == "torus") {
    prompt_torus(dimension, num_objects);
  }
  else if (field_name == "dist2planes") {
    prompt_dist2planes(dimension, num_objects);
  }
  else if (field_name == "edge") {
    prompt_edge(dimension);
  }
  else if (field_name == "octahedron") {
    prompt_octahedron(dimension, num_objects);
  }
  else if (field_name == "constant") {
    prompt_constant_scalar();
  }
  else if (field_name == "constant_unit_gradient") {
    prompt_constant_unit_gradient(dimension);
  }
  else if (field_name == "randomint") {
    prompt_random_int();
  }
  else if (field_name == "randomP") {
    prompt_random_power();
  }
  else if (field_name == "isotable_entry") {
    prompt_isotable_index(dimension);
  }

  if (field_param.flag_flange) 
    { prompt_flange(); }

  if (field_param.flag_wedge)
    { prompt_wedge(dimension); }
}


// **************************************************
// Generate scalar fields
// **************************************************

/// Generate scalar field whose isosurfaces have flanges.
void generate_field_with_flange
(const OBJECT_PROPERTIES & object_properties, SCALAR_GRID & scalar_grid)
{
  const int ifield = field_param.FieldIndex();
  OBJECT_PROPERTIES prop;
  SCALAR_GRID gridB;
  IJK::PROCEDURE_ERROR error("generate_field_with_flange");

  gridB.SetSize(scalar_grid);
  gridB.SetSpacing(scalar_grid);
  prop.Copy(object_properties);

  if (field_info[ifield].flag_length_difference) {

    if (prop.length_difference.size() == 0) 
      { prop.length_difference.push_back(0); }

    DIFF_TYPE length_diff = prop.length_difference[0];

    // Create field whose isosurfaces are narrow cylinders/annuli.
    prop.length_difference[0] = length_diff - 2*prop.flange_width[0];
    field_info[ifield].function_ptr(prop, scalar_grid);
    scalar_grid.Add((-1)*prop.flange_width[0]);

    // Create field whose isosurfaces are short cylinders/annuli.
    prop.length_difference[0] = length_diff + 2*prop.flange_height[0];
    field_info[ifield].function_ptr(prop, gridB);

    // Combine two fields.
    min_scalar(scalar_grid, gridB, scalar_grid);
  }
  else if (field_info[ifield].flag_clip) {

    if (prop.dist2near0.size() == 0) {
      error.AddMessage
        ("Programming error.  Missing distance to near clipping plane.");
      throw error;
    }

    if (prop.dist2far0.size() == 0) {
      error.AddMessage
        ("Programming error.  Missing distance to far clipping plane.");
      throw error;
    }

    // Create field whose isosurfaces are frustrums.
    field_info[ifield].function_ptr(prop, scalar_grid);

    // Create field whose isosurfaces are fatter, shorter frustrums.
    prop.dist2near0[0] += (prop.flange_height[0] + prop.flange_width[0]);
    prop.dist2far0[0] -= (prop.flange_height[0] + prop.flange_width[0]);
    field_info[ifield].function_ptr(prop, gridB);
    gridB.Add((-1)*prop.flange_width[0]);

    // Combine two fields.
    min_scalar(scalar_grid, gridB, scalar_grid);
  }
}

/// Generate scalar field and gradients whose isosurfaces have flanges.
void generate_gradient_field_with_flange
(const OBJECT_PROPERTIES & object_properties, 
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid )
{
  const int ifield = field_param.FieldIndex();
  const bool gradient_discontinuity_zero =
    field_param.gradient_discontinuity_zero;
  OBJECT_PROPERTIES prop;
  SCALAR_GRID gridB;
  GRADIENT_GRID gradientB;
  IJK::PROCEDURE_ERROR error("generate_gradient_field_with_flange");

  gridB.SetSize(scalar_grid);
  gridB.SetSpacing(scalar_grid);
  gradientB.SetSize(gradient_grid);
  gradientB.SetSpacing(gradient_grid);
  prop.Copy(object_properties);

  if (field_info[ifield].flag_length_difference) {

    if (prop.length_difference.size() == 0) 
      { prop.length_difference.push_back(0); }

    DIFF_TYPE length_diff = prop.length_difference[0];

    // Create field whose isosurfaces are narrow cylinders/annuli.
    prop.length_difference[0] = length_diff - 2*prop.flange_width[0];
    field_info[ifield].gradient_function_ptr(prop, scalar_grid, gradient_grid);
    scalar_grid.Add((-1)*prop.flange_width[0]);

    // Create field whose isosurfaces are short cylinders/annuli.
    prop.length_difference[0] = length_diff + 2*prop.flange_height[0];
    field_info[ifield].gradient_function_ptr(prop, gridB, gradientB);

    // Combine two surfaces.
    min_scalar_select_gradient
      (scalar_grid, gradient_grid, gridB, gradientB, 
       gradient_discontinuity_zero, scalar_grid, gradient_grid);
  }
  else if (field_info[ifield].flag_clip) {

    if (prop.dist2near0.size() == 0) {
      error.AddMessage
        ("Programming error.  Missing distance to near clipping plane.");
      throw error;
    }

    if (prop.dist2far0.size() == 0) {
      error.AddMessage
        ("Programming error.  Missing distance to far clipping plane.");
      throw error;
    }

    // Create field whose isosurfaces are frustrums.
    field_info[ifield].gradient_function_ptr(prop, scalar_grid, gradient_grid);

    // Create field whose isosurfaces are fatter, shorter frustrums.
    prop.dist2near0[0] += (prop.flange_height[0] + prop.flange_width[0]);
    prop.dist2far0[0] -= (prop.flange_height[0] + prop.flange_width[0]);
    field_info[ifield].gradient_function_ptr(prop, gridB, gradientB);
    gridB.Add((-1)*prop.flange_width[0]);

    // Combine two fields.
    min_scalar_select_gradient
      (scalar_grid, gradient_grid, gridB, gradientB, 
       gradient_discontinuity_zero, scalar_grid, gradient_grid);  }
}

/// Generate scalar field representing a single object.
void generate_single_object_field
(const OBJECT_PROPERTIES & object_properties, SCALAR_GRID & scalar_grid)
{
  const int ifield = field_param.FieldIndex();

  if (field_param.flag_flange) {
    generate_field_with_flange(object_properties, scalar_grid);
  }
  else {
    field_info[ifield].function_ptr(object_properties, scalar_grid);
  }

  if (field_param.flag_wedge) 
    { intersect_with_wedge(object_properties, scalar_grid); }
}

/// Generate scalar field and gradients representing a single object.
void generate_single_object_gradient_field
(const OBJECT_PROPERTIES & object_properties, 
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  const int ifield = field_param.FieldIndex();

  if (field_param.flag_flange) {
    generate_gradient_field_with_flange
      (object_properties, scalar_grid, gradient_grid);
  }
  else {
    field_info[ifield].gradient_function_ptr
      (object_properties, scalar_grid, gradient_grid);
  }

  if (field_param.flag_wedge) 
    { intersect_with_wedge(object_properties, scalar_grid, gradient_grid); }
}


/// Generate scalar field.
void generate_field(SCALAR_GRID & scalar_grid)
{
  const int dimension = scalar_grid.Dimension();
  const int num_objects = field_param.NumObjects();
  const int ifield = field_param.FieldIndex();
  const std::string field_name = field_info[ifield].name;
  OBJECT_PROPERTIES object_properties;

  object_properties.Copy(field_param);

  if (num_objects == 1 || !field_param.flag_multi_centers) {
    generate_single_object_field(object_properties, scalar_grid);
  }
  else {
    IJK::ARRAY<COORD_TYPE> translate(dimension);

    SCALAR_GRID gridB;
    gridB.SetSize(scalar_grid);
    gridB.SetSpacing(scalar_grid);

    if (field_param.flag_stack) {
      compute_center_translation_vector(dimension, translate.Ptr());

      generate_single_object_field(object_properties, scalar_grid);

      for (int i = 1; i < num_objects; i++) {
        IJK::add_coord
          (dimension, translate.PtrConst(), object_properties.CenterPtr(0),
           object_properties.CenterPtr(0));

        generate_single_object_field(object_properties, gridB);
        min_scalar(scalar_grid, gridB, scalar_grid);
      }
    }
    else {
      generate_single_object_field(object_properties, scalar_grid);

      for (int i = 1; i < num_objects; i++) {
        IJK::copy_coord(dimension, field_param.CenterPtrConst(i), 
                        object_properties.CenterPtr(0));
        if (field_param.NumDirections() > i) {
          IJK::copy_coord(dimension, field_param.DirectionPtrConst(i), 
                          object_properties.DirectionPtr(0));
        }

        generate_single_object_field(object_properties, gridB);
        min_scalar(scalar_grid, gridB, scalar_grid);
      }
    }

  }

}                    

/// Generate scalar field and gradients.
void generate_gradient_field
(SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  const int dimension = scalar_grid.Dimension();
  const int num_objects = field_param.NumObjects();
  const int ifield = field_param.FieldIndex();
  const std::string field_name = field_info[ifield].name;
  const bool gradient_discontinuity_zero =
    field_param.gradient_discontinuity_zero;
  OBJECT_PROPERTIES object_properties;

  object_properties.Copy(field_param);

  if (num_objects == 1 || !field_param.flag_multi_centers) {
    generate_single_object_gradient_field
      (object_properties, scalar_grid, gradient_grid);
  }
  else {
    IJK::ARRAY<COORD_TYPE> translate(dimension);

    SCALAR_GRID gridB;
    GRADIENT_GRID gradientB;
    gridB.SetSize(scalar_grid);
    gridB.SetSpacing(scalar_grid);
    gradientB.SetSize(gradient_grid);
    gradientB.SetSpacing(gradient_grid);

    if (field_param.flag_stack) {
      compute_center_translation_vector(dimension, translate.Ptr());

      generate_single_object_gradient_field
        (object_properties, scalar_grid, gradient_grid);

      for (int i = 1; i < num_objects; i++) {
        IJK::add_coord
          (dimension, translate.PtrConst(), object_properties.CenterPtr(0),
           object_properties.CenterPtr(0));

        generate_single_object_gradient_field
          (object_properties, gridB, gradientB);
        min_scalar_select_gradient
          (scalar_grid, gradient_grid, gridB, gradientB, 
           gradient_discontinuity_zero, scalar_grid, gradient_grid);
      }
    }
    else {
      generate_single_object_gradient_field
        (object_properties, scalar_grid, gradient_grid);

      for (int i = 1; i < num_objects; i++) {
        IJK::copy_coord(dimension, field_param.CenterPtrConst(i), 
                        object_properties.CenterPtr(0));
        if (field_param.NumDirections() > i) {
          IJK::copy_coord(dimension, field_param.DirectionPtrConst(i), 
                          object_properties.DirectionPtr(0));
        }

        generate_single_object_gradient_field
          (object_properties, gridB, gradientB);
        min_scalar_select_gradient
          (scalar_grid, gradient_grid, gridB, gradientB, 
           gradient_discontinuity_zero, scalar_grid, gradient_grid);
      }
    }

  }

}                    


// Set scalar values.
void set_scalar_values(const FIELD_PARAM & field_param, SCALAR_GRID & grid)
{
  const int dimension = field_param.Dimension();
  const SCALAR_TYPE neg_value = 0;
  const SCALAR_TYPE pos_value = field_param.MaxVal();
  const CUBE_TYPE cube(dimension);
  IJK::PROCEDURE_ERROR error("set_scalar_values");

  // Set any cubes in field_param.set_scalar_list.
  for (int i = 0; i < field_param.set_scalar_list.size(); i++) {

    if (field_param.set_scalar_list[i].IsTypeSetIsotableIndex()) {
      const CUBE_INDEX_TYPE cube_index = 
        field_param.set_scalar_list[i].cube_index;
      const ISOTABLE_INDEX_TYPE isotable_index =
        field_param.set_scalar_list[i].isotable_index;
    
      if (!grid.CheckCubeIndex(cube_index, error)) { throw error; }

      set_cube_scalar
        (cube_index, isotable_index, neg_value, pos_value, grid);
    }
    else if (field_param.set_scalar_list[i].IsTypeSetPosVertex()) {
      const VERTEX_INDEX_TYPE iv =
        field_param.set_scalar_list[i].vertex_index;

      grid.Set(iv, pos_value);
    }
    else if (field_param.set_scalar_list[i].IsTypeSetPosEdge()) {
      const VERTEX_INDEX_TYPE iv =
        field_param.set_scalar_list[i].vertex_index;
      const int edge_dir =
        field_param.set_scalar_list[i].direction;
      const VERTEX_INDEX_TYPE iv2 = grid.NextVertex(iv, edge_dir);

      grid.Set(iv, pos_value);
      grid.Set(iv2, pos_value);
    }
    else if (field_param.set_scalar_list[i].IsTypeSetPosFacet()) {
      const VERTEX_INDEX_TYPE iv0 =
        field_param.set_scalar_list[i].vertex_index;
      const int orth_dir =
        field_param.set_scalar_list[i].direction;

      for (int j = 0; j < grid.NumFacetVertices(); j++) {
        const VERTEX_INDEX_TYPE jv = grid.FacetVertex(iv0, orth_dir, j);
        grid.Set(jv, pos_value);
      }
    }
    else if (field_param.set_scalar_list[i].IsTypeSetPosFacetDiagonal()) {
      const VERTEX_INDEX_TYPE iv0 =
        field_param.set_scalar_list[i].vertex_index;
      const int orth_dir =
        field_param.set_scalar_list[i].direction;
      const int jth_vertex =
        field_param.set_scalar_list[i].jth_vertex;
      const int j2 = cube.OppositeFacetVertex(orth_dir, jth_vertex);

      const VERTEX_INDEX_TYPE jv = 
        grid.FacetVertex(iv0, orth_dir, jth_vertex);
      const VERTEX_INDEX_TYPE jv2 = 
        grid.FacetVertex(iv0, orth_dir, j2);

      grid.Set(jv, pos_value);
      grid.Set(jv2, pos_value);
    }

  }
}


// Set scalar values of cube from list.
void set_cube_from_list
(const FIELD_PARAM & field_param, const CUBE_INDEX_TYPE icube, 
 SCALAR_GRID & grid)
{
  IJK::PROCEDURE_ERROR error("set_cube_from_list");

  if (grid.NumCubeVertices() > field_param.scalar_value.size()) {
    error.AddMessage("Error.  Missing scalar values.");
    throw error;
  }

  for (int i = 0; i < grid.NumCubeVertices(); i++) {
    const VERTEX_INDEX_TYPE iv = grid.CubeVertex(icube, i);
    grid.Set(iv, field_param.scalar_value[i]);
  }
}


// Set scalar values of two adjacent cubes from list.
// Second cube is further +(0,0,...,1) translate of icube.
void set_two_adjacent_cubes_from_list
(const FIELD_PARAM & field_param, const CUBE_INDEX_TYPE icube, 
 SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  IJK::PROCEDURE_ERROR error("set_two_adjacent_cubes_from_list");

  if (dimension < 1) {
    error.AddMessage("Error.  Dimension must be positive.");
    throw error;
  }

  if (grid.NumCubeVertices() + grid.NumFacetVertices() > 
      field_param.scalar_value.size()) {
    error.AddMessage("Error.  Missing scalar values.");
    error.AddMessage
      ("  Requires ", grid.NumCubeVertices() + grid.NumFacetVertices(),
       " scalar values.");
    throw error;
  }

  set_cube_from_list(field_param, icube, grid);

  const CUBE_INDEX_TYPE icube2 = icube + grid.AxisIncrement(dimension-1);

  for (int i = 0; i < grid.NumFacetVertices(); i++) {
    const int j = i + grid.NumFacetVertices();
    const VERTEX_INDEX_TYPE jv = grid.CubeVertex(icube2, j);
    grid.Set(jv, field_param.scalar_value[i+grid.NumCubeVertices()]);
  }
}


// Set scalar values of interior vertices from list.
void set_interior_vertices_from_list
(const std::vector<SCALAR_TYPE> & scalar_value, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  const int num_interior_vertices = 
    grid.ComputeNumInteriorVertices();
  IJK::PROCEDURE_ERROR error("set_interior_vertices_from_list");

  if (dimension < 1) {
    // Nothing to set.
    return;
  }

  if (scalar_value.size() < num_interior_vertices) {
    error.AddMessage("Error.  Missing scalar values.");
    error.AddMessage
      ("  Requires ", num_interior_vertices, " scalar values.");
    throw error;
  }

  IJK::FACET_INTERIOR_VERTEX_LIST<VERTEX_INDEX_TYPE> 
    facet_vlist(grid, dimension-1, false, true);

  VERTEX_INDEX_TYPE k = 0;
  for (AXIS_SIZE_TYPE i = 1; i+1 < grid.AxisSize(dimension-1); i++) {
    const VERTEX_INDEX_TYPE facet_increment = i*grid.AxisIncrement(dimension-1);
    for (VERTEX_INDEX_TYPE j = 0; j < facet_vlist.NumVertices(); j++) {
      const VERTEX_INDEX_TYPE jv = 
	facet_vlist.VertexIndex(j) + facet_increment;
      grid.Set(jv, scalar_value[k]);
      k++;
    }
  }

}


void gen_cube(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();

  if (field_param.flag_tilt) {
    IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);

    compute_orthogonal_basis3D
      (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
       zdir.Ptr(), xdir.Ptr(), ydir.Ptr());
    gen_cube
      (prop.CenterPtrConst(0), 
       xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), grid);
  }
  else {
    gen_dist2point_Linf(prop.CenterPtrConst(0), grid);
  }
}

void gen_gradient_cube
(const OBJECT_PROPERTIES & prop, 
SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  const int dimension = scalar_grid.Dimension();
  const bool gradient_discontinuity_zero =
    field_param.gradient_discontinuity_zero;

  if (field_param.flag_tilt) {
    IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);

    compute_orthogonal_basis3D
      (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
       zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

    gen_gradient_cube
      (prop.CenterPtrConst(0), 
       xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), 
       scalar_grid, gradient_grid);
  }
  else {
    gen_gradient_dist2point_Linf
      (prop.CenterPtrConst(0), gradient_discontinuity_zero, 
       scalar_grid, gradient_grid);
  }
}

void gen_sphere(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  gen_dist2point_L2(prop.CenterPtrConst(0), grid);
}

void gen_gradient_sphere
(const OBJECT_PROPERTIES & prop, 
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  gen_gradient_dist2point_L2
    (prop.CenterPtrConst(0), scalar_grid, gradient_grid);
}

void gen_closed_cylinder
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & scalar_grid)
{
  gen_closed_cylinder
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.length_difference[0], scalar_grid);
}

void gen_gradient_closed_cylinder
(const OBJECT_PROPERTIES & prop, 
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  gen_gradient_closed_cylinder
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.length_difference[0], scalar_grid, gradient_grid);
}

void gen_closed_cone
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  const double height0 = 0;

  if (dimension < 1) { return; }

  if (prop.flag_smooth_tip) {
    gen_closed_cone_smooth_tip
      (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
       prop.angle[0], height0, grid);
  }
  else {
    gen_closed_cone
      (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
       prop.angle[0], height0, grid);
  }
}

void gen_gradient_closed_cone
(const OBJECT_PROPERTIES & prop, 
 SCALAR_GRID & grid, GRADIENT_GRID & gradient)
{
  const int dimension = grid.Dimension();
  const double height0 = 0;
  SCALAR_GRID gridB;

  if (dimension < 1) { return; }

  if (prop.flag_smooth_tip) {
    gen_gradient_closed_cone_smooth_tip
      (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
       prop.angle[0], height0, grid, gradient);
  }
  else {
    gen_gradient_closed_cone
      (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
       prop.angle[0], height0, grid, gradient);
  }
}

void gen_frustrum
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  const double height0 = 0;

  if (dimension < 1) { return; }

  gen_frustrum
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.angle[0], 
     prop.dist2near0[0], prop.dist2far0[0],
     grid);
}

void gen_gradient_frustrum
(const OBJECT_PROPERTIES & prop, 
 SCALAR_GRID & grid, GRADIENT_GRID & gradient)
{
  const int dimension = grid.Dimension();
  const double height0 = 0;

  if (dimension < 1) { return; }

  gen_gradient_frustrum
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.angle[0], 
     prop.dist2near0[0], prop.dist2far0[0],
     grid, gradient);
}

// Generate a field with cannon shaped isosurfaces.
void gen_cannon
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  const double height0 = 0;
  IJK::PROCEDURE_ERROR error("gen_cannon");

  if (dimension < 1) { return; }

  if (!prop.CheckNumCenters(1, error)) { throw error; }
  if (!prop.CheckNumDirections(1, error)) { throw error; }
  if (!prop.CheckNumAngles(1, error)) { throw error; }
  if (!prop.CheckNumDist2BallCenters(1, error)) { throw error; }

  gen_cannon
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.angle[0], 
     prop.dist2near0[0], prop.dist2ball_center[0],
     grid);
}

// Generate a field with cannon shaped isosurfaces.
void gen_gradient_cannon
(const OBJECT_PROPERTIES & prop, 
 SCALAR_GRID & grid, GRADIENT_GRID & gradient)
{
  const int dimension = grid.Dimension();
  const double height0 = 0;
  IJK::PROCEDURE_ERROR error("gen_cannon");

  if (dimension < 1) { return; }

  if (!prop.CheckNumCenters(1, error)) { throw error; }
  if (!prop.CheckNumDirections(1, error)) { throw error; }
  if (!prop.CheckNumAngles(1, error)) { throw error; }
  if (!prop.CheckNumDist2BallCenters(1, error)) { throw error; }

  gen_gradient_cannon
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.angle[0], 
     prop.dist2near0[0], prop.dist2ball_center[0],
     grid, gradient);
}

/// Compute a direction at a given angle from \a dir1.
/// @param dir1[] Direction.
/// @param dir2[] Second direction.
/// @param angle Angle in degrees.
/// @param[out] new_dir[] Direction at given with \a dir1 in plane
///                  determined by \a dir1 and \a dir2.
/// @pre dir1[] is a unit vector.
/// @pre dir2[] is a unit vector.
template <typename COORD1_TYPE, typename COORD2_TYPE, typename COORD3_TYPE>
void compute_dir_at_angle
(const int dimension, const COORD1_TYPE dir1[], 
 const COORD2_TYPE dir2[], const ANGLE_TYPE angle,
 COORD3_TYPE new_dir[])
{
  const double radians = M_PI*angle/180.0;

  IJK::multiply_coord(dimension, cos(radians), dir1, new_dir);
  IJK::add_scaled_coord
    (dimension, sin(radians), dir2, new_dir, new_dir);
}

template <typename T>
T ijk_mod(const T x, const int y)
{
  if (x >= 0) {
    return(fmod(x, y));
  }
  else {
    return(y-fmod(-x, y));
  }
}


/// Compute two normal vectors to wedge planes.
/// @param[out] normal[] Normals to wedge planes.
template <typename CTYPE>
void compute_wedge_plane_normals
(const OBJECT_PROPERTIES & prop, CTYPE normal[2*DIM3])
{
  IJK::ARRAY<COORD_TYPE> xdir(DIM3), ydir(DIM3), zdir(DIM3);
  IJK::PROCEDURE_ERROR error("compute_wedge_plane_normals");

  if (prop.NumWedgeAngles() < 1) {
    error.AddMessage("Programming error.  Missing wedge angle.");
    throw error;
  }

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0), 
     zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  IJK::copy_coord(DIM3, xdir.PtrConst(), normal);

  compute_dir_at_angle
    (DIM3, xdir.PtrConst(), ydir.PtrConst(),
     prop.wedge_angle[0], normal+DIM3);
  IJK::multiply_coord(DIM3, -1, normal+DIM3, normal+DIM3);
}

/// Intersect scalar field with wedge field.
void intersect_with_wedge
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int num_planes = 2;
  double normal[2*DIM3];
  SCALAR_GRID gridB;
  IJK::PROCEDURE_ERROR error("intersect_with_wedge");

  if (!check_dimension(grid.Dimension(), DIM3, "Wedge", error))
    { throw error; }

  gridB.SetSize(grid);
  gridB.SetSpacing(grid);

  compute_wedge_plane_normals(prop, normal);

  if (ijk_mod(field_param.wedge_angle[0], 360) <= 180) {
    gen_max_dist2planes
      (prop.CenterPtrConst(0), normal, num_planes, gridB);
  }
  else {
    gen_min_dist2planes
      (prop.CenterPtrConst(0), normal, num_planes, gridB);
  }
  gridB.Add(field_param.wedge_isovalue[0]);

  max_scalar(grid, gridB, grid);
}

/// Intersect scalar field with wedge field.
void intersect_with_wedge
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & grid, GRADIENT_GRID & gradient)
{
  const int num_planes = 2;
  double normal[2*DIM3];
  SCALAR_GRID gridB;
  GRADIENT_GRID gradientB;
  IJK::PROCEDURE_ERROR error("intersect_with_wedge");

  if (!check_dimension(grid.Dimension(), DIM3, "Wedge", error))
    { throw error; }

  gridB.SetSize(grid);
  gridB.SetSpacing(grid);
  gradientB.SetSize(gradient);
  gradientB.SetSpacing(gradient);

  compute_wedge_plane_normals(prop, normal);

  if (ijk_mod(field_param.wedge_angle[0], 360) <= 180) {
    gen_gradient_max_dist2planes
      (prop.CenterPtrConst(0), normal, num_planes, gridB, gradientB);
  }
  else {
    gen_gradient_min_dist2planes
      (prop.CenterPtrConst(0), normal, num_planes, gridB, gradientB);
  }
  gridB.Add(field_param.wedge_isovalue[0]);

  max_scalar_select_gradient
    (grid, gradient, gridB, gradientB, 
     field_param.gradient_discontinuity_zero, grid, gradient);
}

void gen_octahedron(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  gen_dist2point_L1(prop.CenterPtrConst(0), grid);
}

void gen_gradient_octahedron
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  if (field_param.gradient_discontinuity_zero) {
    gen_gradient_dist2point_L1_zero
      (prop.CenterPtrConst(0), scalar_grid, gradient_grid);
  }
  else {
    gen_gradient_dist2point_L1_select
      (prop.CenterPtrConst(0), scalar_grid, gradient_grid);
  }
}

void gen_open_cylinder(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  gen_dist2line_L2
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0), grid);

}

void gen_gradient_open_cylinder
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  gen_gradient_dist2line_L2
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0), 
     scalar_grid, gradient_grid);
}

void gen_square_cylinder
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);
  IJK::PROCEDURE_ERROR error("gen_square_cylinder");

  if (!check_dimension(dimension, DIM3, "Square cylinder", error))
    { throw error; }

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  gen_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), 
     prop.length_difference[0], grid);
}

void gen_gradient_square_cylinder
(const OBJECT_PROPERTIES & prop, 
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  const int dimension = scalar_grid.Dimension();
  IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);
  IJK::PROCEDURE_ERROR error("gen_gradient_square_cylinder");

  if (!check_dimension(dimension, DIM3, "Square cylinder", error))
    { throw error; }

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  gen_gradient_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), 
     prop.length_difference[0], scalar_grid, gradient_grid);
}

void gen_square_cylinder_rot45
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  gen_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), 
     prop.length_difference[0], grid);
}

void gen_gradient_square_cylinder_rot45
(const OBJECT_PROPERTIES & prop, 
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  const int dimension = scalar_grid.Dimension();
  IJK::ARRAY<COORD_TYPE> xdir(dimension), ydir(dimension), zdir(dimension);

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     zdir.Ptr(), xdir.Ptr(), ydir.Ptr());

  gen_gradient_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     xdir.PtrConst(), ydir.PtrConst(), zdir.PtrConst(), 
     prop.length_difference[0], scalar_grid, gradient_grid);
}

void gen_crossing_square_cylinders
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  IJK::ARRAY<COORD_TYPE> dir0(dimension);
  IJK::ARRAY<COORD_TYPE> dir1(dimension);
  IJK::ARRAY<COORD_TYPE> dir2(dimension);

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     dir0.Ptr(), dir1.Ptr(), dir2.Ptr());
	
  SCALAR_GRID gridA, gridB, gridC;
  gridA.SetSize(grid);
  gridA.SetSpacing(grid);
  gridB.SetSize(grid);
  gridB.SetSpacing(grid);
  gridC.SetSize(grid);
  gridC.SetSpacing(grid);

  gen_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir1.PtrConst(), dir2.PtrConst(),
     prop.length_difference[0], gridA);
  gen_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     dir1.PtrConst(), dir2.PtrConst(), dir0.PtrConst(),
     prop.length_difference[0], gridB);
  gen_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir2.PtrConst(), dir1.PtrConst(),
     prop.length_difference[0], gridC);

  min_scalar(gridA, gridB, grid);
  min_scalar(gridC, grid, grid);
}

void gen_gradient_crossing_square_cylinders
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  const int dimension = scalar_grid.Dimension();
  IJK::ARRAY<COORD_TYPE> dir0(dimension);
  IJK::ARRAY<COORD_TYPE> dir1(dimension);
  IJK::ARRAY<COORD_TYPE> dir2(dimension);

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     dir0.Ptr(), dir1.Ptr(), dir2.Ptr());
	
  SCALAR_GRID gridA, gridB, gridC;
  GRADIENT_GRID gradientA, gradientB, gradientC;
  gridA.SetSize(scalar_grid);
  gridA.SetSpacing(scalar_grid);
  gridB.SetSize(scalar_grid);
  gridB.SetSpacing(scalar_grid);
  gridC.SetSize(scalar_grid);
  gridC.SetSpacing(scalar_grid);
  gradientA.SetSize(gradient_grid);
  gradientA.SetSpacing(gradient_grid);
  gradientB.SetSize(gradient_grid);
  gradientB.SetSpacing(gradient_grid);
  gradientC.SetSize(gradient_grid);
  gradientC.SetSpacing(gradient_grid);

  gen_gradient_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir1.PtrConst(), dir2.PtrConst(),
     prop.length_difference[0], gridA, gradientA);
  gen_gradient_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     dir1.PtrConst(), dir2.PtrConst(), dir0.PtrConst(),
     prop.length_difference[0], gridB, gradientB);
  gen_gradient_closed_square_cylinder
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir2.PtrConst(), dir1.PtrConst(),
     prop.length_difference[0], gridC, gradientC);

  min_scalar_select_gradient
    (gridA, gradientA, gridB, gradientB, 
     field_param.gradient_discontinuity_zero, scalar_grid, gradient_grid);
  min_scalar_select_gradient
    (scalar_grid, gradient_grid, gridC, gradientC,
     field_param.gradient_discontinuity_zero, scalar_grid, gradient_grid);
}

void gen_crossing_square45_cylinders
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  IJK::ARRAY<COORD_TYPE> origin(dimension);
  IJK::ARRAY<COORD_TYPE> dir0(dimension);
  IJK::ARRAY<COORD_TYPE> dir1(dimension);
  IJK::ARRAY<COORD_TYPE> dir2(dimension);

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     dir0.Ptr(), dir1.Ptr(), dir2.Ptr());
	
  SCALAR_GRID gridA, gridB, gridC;
  gridA.SetSize(grid);
  gridA.SetSpacing(grid);
  gridB.SetSize(grid);
  gridB.SetSpacing(grid);
  gridC.SetSize(grid);
  gridC.SetSpacing(grid);

  gen_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir1.PtrConst(), dir2.PtrConst(),
     prop.length_difference[0], gridA);
  gen_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     dir1.PtrConst(), dir2.PtrConst(), dir0.PtrConst(),
     prop.length_difference[0], gridB);
  gen_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir2.PtrConst(), dir1.PtrConst(),
     prop.length_difference[0], gridC);

  min_scalar(gridA, gridB, grid);
  min_scalar(gridC, grid, grid);
}

void gen_gradient_crossing_square45_cylinders
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  const int dimension = scalar_grid.Dimension();
  IJK::ARRAY<COORD_TYPE> dir0(dimension);
  IJK::ARRAY<COORD_TYPE> dir1(dimension);
  IJK::ARRAY<COORD_TYPE> dir2(dimension);

  compute_orthogonal_basis3D
    (prop.DirectionPtrConst(0), prop.SideDirectionPtrConst(0),
     dir0.Ptr(), dir1.Ptr(), dir2.Ptr());
	
  SCALAR_GRID gridA, gridB, gridC;
  GRADIENT_GRID gradientA, gradientB, gradientC;
  gridA.SetSize(scalar_grid);
  gridA.SetSpacing(scalar_grid);
  gridB.SetSize(scalar_grid);
  gridB.SetSpacing(scalar_grid);
  gridC.SetSize(scalar_grid);
  gridC.SetSpacing(scalar_grid);
  gradientA.SetSize(gradient_grid);
  gradientA.SetSpacing(gradient_grid);
  gradientB.SetSize(gradient_grid);
  gradientB.SetSpacing(gradient_grid);
  gradientC.SetSize(gradient_grid);
  gradientC.SetSpacing(gradient_grid);

  gen_gradient_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir1.PtrConst(), dir2.PtrConst(),
     prop.length_difference[0], gridA, gradientA);
  gen_gradient_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     dir1.PtrConst(), dir2.PtrConst(), dir0.PtrConst(),
     prop.length_difference[0], gridB, gradientB);
  gen_gradient_closed_square_cylinder_rot45
    (prop.CenterPtrConst(0), 
     dir0.PtrConst(), dir2.PtrConst(), dir1.PtrConst(),
     prop.length_difference[0], gridC, gradientC);

  min_scalar_select_gradient
    (gridA, gradientA, gridB, gradientB, 
     field_param.gradient_discontinuity_zero, scalar_grid, gradient_grid);
  min_scalar_select_gradient
    (scalar_grid, gradient_grid, gridC, gradientC,
     field_param.gradient_discontinuity_zero, scalar_grid, gradient_grid);
}

void gen_annulus(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  gen_annulus
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.radius[0], prop.length_difference[0], grid);
}

void gen_gradient_annulus
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid, GRADIENT_GRID & gradient)
{
  gen_gradient_annulus
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.radius[0], prop.length_difference[0], 
     grid, gradient);
}

void gen_torus(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  gen_torus
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.radius[0], grid);
}

void gen_gradient_torus
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid, GRADIENT_GRID & gradient)
{
  gen_gradient_torus
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0),
     prop.radius[0], grid, gradient);
}

void gen_dist2planes(const OBJECT_PROPERTIES & prop, SCALAR_GRID & scalar_grid)
{
  const int num_planes = field_param.NumObjects(); 
  gen_max_dist2planes(prop.CenterPtrConst(0),
                      prop.NormalPtrConst(), num_planes,
                      scalar_grid);
}

void gen_gradient_dist2planes
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient)
{
  const int num_planes = field_param.NumObjects();

  gen_gradient_max_dist2planes
    (prop.CenterPtrConst(0), prop.NormalPtrConst(), num_planes,
     scalar_grid, gradient);
}

void gen_edge(const OBJECT_PROPERTIES & prop, SCALAR_GRID & scalar_grid)
{
  const int dimension = scalar_grid.Dimension();
  const int num_planes = 2;
  IJK::ARRAY<double> normal(2*dimension);
  IJK::ARRAY<double> v2(dimension);

  IJK::copy_coord(dimension, prop.NormalPtrConst(), normal.Ptr());
  normalize_vector(dimension, normal.Ptr());

  compute_normalized_orthogonal_vector
    (dimension, prop.RotationDirectionPtrConst(), 
     normal.PtrConst(), v2.Ptr());

  compute_dir_at_angle
    (dimension, normal.PtrConst(), v2.PtrConst(), prop.angle[0], 
     normal.Ptr()+dimension);
  IJK::multiply_coord
    (DIM3, -1, normal.Ptr()+dimension, normal.Ptr()+dimension);

  if (ijk_mod(prop.angle[0], 360) <= 180) {
    gen_max_dist2planes
      (prop.CenterPtrConst(0), normal.PtrConst(), num_planes, scalar_grid);
  }
  else {
    gen_min_dist2planes
      (prop.CenterPtrConst(0), normal.PtrConst(), num_planes, scalar_grid);
  }
}

void gen_gradient_edge
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient)
{
  const int dimension = scalar_grid.Dimension();
  const int num_planes = 2;
  IJK::ARRAY<double> normal(2*dimension);
  IJK::ARRAY<double> v2(dimension);

  IJK::copy_coord(dimension, prop.NormalPtrConst(), normal.Ptr());
  normalize_vector(dimension, normal.Ptr());

  compute_normalized_orthogonal_vector
    (dimension, prop.RotationDirectionPtrConst(), normal.PtrConst(), v2.Ptr());

  compute_dir_at_angle
    (dimension, normal.PtrConst(), v2.PtrConst(), prop.angle[0], 
     normal.Ptr()+dimension);
  IJK::multiply_coord(DIM3, -1, normal.Ptr()+dimension, normal.Ptr()+dimension);

 if (ijk_mod(prop.angle[0], 360) <= 180) {
   gen_gradient_max_dist2planes
     (prop.CenterPtrConst(0), normal.PtrConst(), num_planes, 
      scalar_grid, gradient);
 }
 else {
   gen_gradient_min_dist2planes
     (prop.CenterPtrConst(0), normal.PtrConst(), num_planes, 
      scalar_grid, gradient);
 }

}

void gen_constant_scalar_field
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  gen_constant(field_param.ScalarValue(0), grid);
}

void gen_sfield_with_constant_unit_gradient
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  gen_constant_unit_gradient
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0), grid);
}

void gen_gfield_with_constant_unit_gradient
(const OBJECT_PROPERTIES & prop,
 SCALAR_GRID & scalar_grid, GRADIENT_GRID & gradient_grid)
{
  gen_constant_unit_gradient
    (prop.CenterPtrConst(0), prop.DirectionPtrConst(0), 
     scalar_grid, gradient_grid);
}


void set_scalar_field(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  prompt_set_scalar_field();
  gen_constant(0, grid);
  set_scalar_values(field_param, grid);
}


void set_xcube(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  CUBE_INDEX_TYPE icube;

  // Generate cube index of cube (1,1,1,...)
  icube = 0;
  for (int d = 0; d < dimension; d++)
    { icube += grid.AxisIncrement(d); }

  prompt_set_xcube(grid.Dimension());
  gen_constant(0, grid);
  set_cube_from_list(field_param, icube, grid);
  grid.SetBoundaryFromInteriorNeighbors();
}


void set_x2cube(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();
  CUBE_INDEX_TYPE icube;

  // Generate cube index of cube (1,1,1,...)
  icube = 0;
  for (int d = 0; d < dimension; d++)
    { icube += grid.AxisIncrement(d); }

  prompt_set_x2cube(grid.Dimension());
  gen_constant(0, grid);
  set_two_adjacent_cubes_from_list(field_param, icube, grid);
  grid.SetBoundaryFromInteriorNeighbors();
}


void set_interior_grid_vertices
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int dimension = grid.Dimension();

  /* DEBUG
  prompt_set_interior_vertices(grid.Dimension());
  */
  gen_constant(0, grid);
  set_interior_vertices_from_list(field_param.scalar_value, grid);
  grid.SetBoundaryFromInteriorNeighbors();
}


void gen_random_int(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  IJK::PROCEDURE_ERROR error("gen_random_int");

  gen_random_int(field_param.MaxVal(), field_param.RandomSeed(), grid);
  set_scalar_values(field_param, grid);
}

void gen_random_power_signed
(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  const int base = field_param.RandomBase();

  const SCALAR_TYPE s = pow(base, field_param.MaxExponent());
  field_param.maxval.Set(s);

  gen_random_power_signed
    (base, field_param.MaxExponent(), field_param.RandomSeed(), grid);
}

void gen_isotable_entry(const OBJECT_PROPERTIES & prop, SCALAR_GRID & grid)
{
  IJK::PROCEDURE_ERROR error("gen_isotable_entry");

  for (int d = 0; d < grid.Dimension(); d++) {
    if (grid.AxisSize(d) != 2) {
      error.AddMessage("Programming error. Incorrect axis size.");
      error.AddMessage("Axis size for axis ", d, " is ", grid.AxisSize(d),
                       ".  Should be 2.");
      throw error;
    }
  }

  set_cube_scalar(0, field_param.IsotableIndex(), -1, 1, grid);
}


// **************************************************
// Read functions
// **************************************************

void read_dimension()
{
  int dimension;

  cout << "Enter domain dimension: ";
  cin >> dimension;

  field_param.SetDimension(dimension);

  if (dimension < 1) {
    cerr << "Illegal dimension.  Dimension must be at least 1." << endl;
    exit(20);
  };
}

void read_num_vert_per_axis()
{
  int size;

  cout << "Enter # grid vertices per axis: ";
  cin >> size;

  grid_axis_size.clear();
  grid_axis_size.push_back(size);

  if (size < 2) {
    cerr << "Illegal grid axis size.  Must be at least 2." << endl;
    exit(30);
  };

}

// **************************************************
// FIELD_INFO functions
// **************************************************

void init_field_info(vector<FIELD_INFO> & field_info)
{
  FIELD_INFO info;

  info.Set("cube", gen_cube, gen_gradient_cube);
  info.SetCubeFlags();
  info.SetDescription
    ("L-infinity distance to a point. Isosurfaces are cubes.");
  field_info.push_back(info);

  info.Set("sphere", gen_sphere, gen_gradient_sphere);
  info.SetSphereFlags();
  info.SetDescription
    ("L2 distance to a point. Isosurfaces are spheres.");
  field_info.push_back(info);

  info.Set("cylinder", gen_closed_cylinder, gen_gradient_closed_cylinder);
  info.SetCylinderFlags();
  info.SetDescription
    ("Minimum of L2 distance to a line and an orthogonal plane.");
  info.AddToDescription("Isosurfaces are closed cylinders.");
  field_info.push_back(info);

  info.Set("cone", gen_closed_cone, gen_gradient_closed_cone);
  info.SetConeFlags();
  info.SetDescription
    ("Max signed distance (L2) to a cone and a clipping plane.");
  info.AddToDescription("Isosurfaces are cones with a circle base.");
  field_info.push_back(info);

  info.Set("frustrum", gen_frustrum, gen_gradient_frustrum);
  info.SetFrustrumFlags();
  info.SetDescription
    ("Max signed distance (L2) to a cone and two clipping planes.");
  info.AddToDescription("Isosurfaces are frustra (truncated cones.)");
  field_info.push_back(info);

  info.Set("cannon", gen_cannon, gen_gradient_cannon);
  info.SetFrustrumFlags();
  info.SetDescription("Cannon shaped isosurfaces.");
  field_info.push_back(info);

  info.Set("square", gen_square_cylinder, gen_gradient_square_cylinder);
  info.SetCylinderFlags();
  info.SetDescription("Isosurfaces are closed square cylinders.");
  field_info.push_back(info);

  info.Set("square45", gen_square_cylinder_rot45,
           gen_gradient_square_cylinder_rot45);
  info.SetCylinderFlags();
  info.SetDescription("Isosurfaces are closed square cylinders");
  info.AddToDescription("  rotated 45 degrees around their central axes.");
  field_info.push_back(info);

  info.Set("squareX", gen_crossing_square_cylinders, 
           gen_gradient_crossing_square_cylinders);
  info.SetCrossingCylinderFlags();
  info.SetDescription
    ("Isosurfaces are 3 intersecting closed square cylinders.");
  field_info.push_back(info);

  info.Set("square45X", gen_crossing_square45_cylinders,
           gen_gradient_crossing_square45_cylinders);
  info.SetCrossingCylinderFlags();
  info.SetDescription
    ("Isosurfaces are 3 intersecting closed square cylinders.");
  info.AddToDescription("Square cylinders are rotated 45 degrees"); 
  info.AddToDescription("  around their central axes.");
  field_info.push_back(info);

  info.Set("annulus", gen_annulus, gen_gradient_annulus);
  info.SetAnnulusFlags();
  info.SetDescription
    ("Minimum of L2 distance to a cylinder and an orthogonal plane.");
  info.AddToDescription("Isosurfaces are thickened annuli.");
  field_info.push_back(info);

  info.Set("torus", gen_torus, gen_gradient_torus);
  info.SetTorusFlags();
  info.SetDescription("L2 distance to a circle.  Isosurfaces are torii.");
  field_info.push_back(info);

  info.Set("edge", gen_edge, gen_gradient_edge);
  info.SetFlagCenter(true);
  info.SetDescription
    ("Scalar field is the maximum signed distance to two planes.");
  info.AddToDescription
    ("Isosurfaces have a single edge with the given dihedral angle.");
  field_info.push_back(info);

  info.Set("dist2planes", gen_dist2planes, gen_gradient_dist2planes);
  info.SetFlagCenter(true);
  info.SetDescription
    ("Scalar field is the maximum signed distance to a set of planes.");
  field_info.push_back(info);

  info.Set("octahedron", gen_octahedron, gen_gradient_octahedron);
  info.SetOctahedronFlags();
  info.SetDescription("Isosurfaces are octahedra.");
  field_info.push_back(info);

  info.Set("constant", gen_constant_scalar_field);
  info.SetDescription("Field with constant scalar value.");
  field_info.push_back(info);

  info.Set("constant_unit_gradient", 
           gen_sfield_with_constant_unit_gradient,
           gen_gfield_with_constant_unit_gradient);
  info.SetFlagCenter(true);
  info.SetDescription("Field with constant unit gradient.");
  info.AddToDescription("Isosurfaces are planes orthogonal to gradient.");
  field_info.push_back(info);

  info.Set("open_cylinder", gen_open_cylinder, gen_gradient_open_cylinder);
  info.SetCylinderFlags();
  info.SetDescription
    ("L2 distance to a line.  Isosurfaces are open cylinders.");
  field_info.push_back(info);

  info.Set("set", set_scalar_field);
  info.SetDescription("Set scalar values in individual cubes.");
  info.AddToDescription("Extend scalar values to grid boundary.");
  info.AddToDescription("Creates grid with axis size 4.");
  info.SetScalarSetFlags();
  field_info.push_back(info);

  info.Set("setv_interior", set_interior_grid_vertices);
  info.SetDescription("Set interior grid vertices.");
  info.AddToDescription
    ("Sets boundary vertices from neighboring interior vertices.");
  field_info.push_back(info);

  info.Set("set_xcube", set_xcube);
  info.SetDescription("Set scalar values in a single cube. (Extend values to boundary.)");
  info.AddToDescription("Grid is 4x4x4 (or 4x...x4.)");
  info.AddToDescription
    ("Sets boundary vertices from neighboring interior vertices.");
  field_info.push_back(info);

  info.Set("set_x2cube", set_x2cube);
  info.SetDescription("Set scalar values in two cubes. (Extend values to boundary.)");
  info.AddToDescription("Grid is 4x4x5 (or 4x...x4x5.)");
  info.AddToDescription
    ("Sets boundary vertices from neighboring interior vertices.");
  field_info.push_back(info);

  info.Set("randomint", gen_random_int);
  info.SetRandomFlags();
  info.SetDescription
    ("Assign random integers as scalar values at grid vertices.");
  field_info.push_back(info);

  info.Set("randomP", gen_random_power_signed);
  info.SetRandomPowerFlags();
  info.SetDescription
    ("Assign base^{random power} as scalar values at grid vertices.");
  field_info.push_back(info);

  info.Set("isotable_entry", gen_isotable_entry);
  info.SetDescription
    ("Assign +1/-1 scalar values corresponding to +/- configuration");
  info.AddToDescription("  in the isosurface lookup table.");
  field_info.push_back(info);
}

/// Output field names.
/// @param with_gradients 
///           If true, output only fields with gradients implemented.
void out_field_names
(const vector<FIELD_INFO> & field_info, 
 const bool with_gradients)
{
  bool newline_flag(false);

  int iout = 0;
  for (int i = 0; i < field_info.size(); i++) {
    if (!with_gradients || 
        field_info[i].IsGradientImplemented()) {
      if (iout > 0) { 
        cout << ","; 
        newline_flag = false;
        if (iout%5 == 0) { 
          newline_flag = true;
          cout << endl; 
        }
      }
      cout << field_info[i].name;
      iout++;
    }
  }

  if (!newline_flag) { cout << endl; }
}

/// Output field descriptions.
void out_field_descriptions
(const vector<FIELD_INFO> & field_info)
{
  const int COL1_WIDTH = 15;

  for (int i = 0; i < field_info.size(); i++) {
    int num_description_lines = field_info[i].description.size();
    string s = field_info[i].name + ":";

    cout << left << setw(COL1_WIDTH) << s;

    int j = 0;
    if ((num_description_lines > 0) && (s.size()+1 <= COL1_WIDTH)) {
      cout << field_info[i].description[0];
      j++;
    }
    cout << endl;

    while (j < num_description_lines) {
      cout << setw(COL1_WIDTH) << " " 
           << field_info[i].description[j] << endl;
      j++;
    }
  }

}

void read_field_name(const vector<FIELD_INFO> & field_info)
{
  string s;
  int ifield;

  cout << "Enter scalar field name (h for help): ";
  cin >> s;

  while (true) {
    if (!find_geom_info_name(field_info, s, ifield)) {
      if (s != "h" && s != "help") {
        cout << "Illegal field name: " << s << endl;
      }
    }
    else if (output_gradients &&
             !field_info[ifield].IsGradientImplemented()) {
      cout << "Computing gradients not implemented for function " 
             << field_info[ifield].name << "." << endl;
    }
    else {
      break;
    }

    cout << endl;
    cout << "Valid field names: ";
    out_field_names(field_info, output_gradients);
    cout << endl;
    cout << "Enter scalar field name (h for help): ";
    cin >> s;
  }

  field_param.SetFieldIndex(ifield);
}

// **************************************************
// FIELD_PARAM functions
// **************************************************

void set_field_param(const vector<FIELD_INFO> & field_info)
{
  int ifield = field_param.FieldIndex();

  if (field_info[ifield].flag_allow_multi_centers) {
    if (field_param.num_objects.IsSet()) {
      if (field_param.NumObjects() > 1) 
        { field_param.flag_multi_centers = true; }
    }
  }
  else
    { field_param.flag_multi_centers = false; }

  if (!field_info[ifield].flag_allow_tilt) 
    { field_param.flag_tilt = false; }

  if (!field_info[ifield].flag_allow_flange)
    { field_param.flag_flange = false; }

  if (!field_info[ifield].flag_allow_wedge)
    { field_param.flag_wedge = false; }
}

void set_random_param
(const int dimension, const int num_objects)
{
  // Local types.
  typedef std::uniform_int_distribution<long>
    _UNIFORM_INT_DISTRIBUTION;
  typedef _UNIFORM_INT_DISTRIBUTION *
    _UNIFORM_INT_DISTRIBUTION_PTR;

  const int MAX_GRAD_COORD = 5;
  const int MIN_GRAD_COORD = -MAX_GRAD_COORD;
  IJK::PROCEDURE_ERROR error("set_random_param");

  if (field_param.flag_random_centers) {

    if (!field_param.randompos_seed.IsSet()) {
      error.AddMessage("Programming error. Missing random position seed.");
      error.AddMessage("Random position seed must be set to generate random centers.");
      throw error;
    }

    IJK::ARRAY<_UNIFORM_INT_DISTRIBUTION_PTR> generate_randompos(dimension);
    std::minstd_rand randompos_engine(field_param.randompos_seed.Value());

    for (int d = 0; d < dimension; d++) {
      generate_randompos[d] = 
	new std::uniform_int_distribution<long>
	(0, field_param.axis_size[d]);
    }


    for (int i = field_param.NumCenters(); i < num_objects; i++) {
      for (int d = 0; d < dimension; d++) {
	long x = (*generate_randompos[d])(randompos_engine);

	// Not necessary, but just in case.
        x = x%(field_param.axis_size[d]);

        field_param.center.push_back(x);
      }
    }

    for (int d = 0; d < dimension; d++) {
      delete generate_randompos[d]; 
      generate_randompos[d] = NULL;
    }
  }


  if (field_param.flag_random_directions) {

    if (!field_param.randomdir_seed.IsSet()) {
      error.AddMessage("Programming error. Missing random direction seed.");
      error.AddMessage("Random direction seed must be set to generate random directions.");
      throw error;
    }

    _UNIFORM_INT_DISTRIBUTION
      generate_randomdir(MIN_GRAD_COORD, MAX_GRAD_COORD);
    std::minstd_rand randomdir_engine(field_param.randomdir_seed.Value());

    for (int i = field_param.NumDirections(); i < num_objects; i++) {
      for (int d = 0; d < dimension; d++) {
	long x = generate_randomdir(randomdir_engine);
        field_param.direction.push_back(x);
      }
    }
  }

}


// **************************************************
// STRING PROCESSING
// **************************************************

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

void construct_gradient_filename
(const char * scalar_filename, string & gradient_filename)
{
  string prefix, suffix;

  gradient_filename = scalar_filename;

  split_string(gradient_filename, '.', prefix, suffix);
  if (suffix == "nrrd" || suffix == "nhdr") { 
    gradient_filename = prefix; 
    gradient_filename = gradient_filename + ".grad." + suffix;
  }
  else {
    gradient_filename = gradient_filename + ".grad";
  }
}


// **************************************************
// COMMAND LINE PARSER
// **************************************************

void check_input_param(const int dimension)
{
  int ifield = field_param.FieldIndex();
  std::string field_name = field_info[ifield].name;

  if (field_param.center.size()%dimension != 0) {
    cerr << "Input error."
         << " Number of arguments in option -center does not match dimension "
         << dimension << ".";
    exit(50);
  }

  if (field_param.direction.size()%dimension != 0) {
    cerr << "Input error."
         << " Number of arguments in option -dir does not match dimension "
         << dimension << ".";
    exit(50);
  }

  if (field_param.normal.size()%dimension != 0) {
    cerr << "Input error."
         << " Number of arguments in option -dir does not match dimension "
         << dimension << ".";
    exit(50);
  }

  if (field_name == "isotable_entry") {

    ISOTABLE_INDEX_TYPE num_table_entries;
    compute_num_table_entries(dimension, num_table_entries);

    if (field_param.isotable_index.IsSet()) {

      if (field_param.IsotableIndex() < 0) {
        cerr << "Input error.  Isosurface table index must be non-negative."
             << endl;
        exit(55);
      }

      if (field_param.IsotableIndex() >= num_table_entries) {
        cerr << "Input error.  Isosurface table index must be less than "
             << num_table_entries << "." << endl;
        exit(57);
      }
    }
  }

  if (field_name == "set_xcube") {
    const int num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
    const int num_scalar_values = field_param.scalar_value.size();

    if ((num_scalar_values > 0) && (num_scalar_values < num_cube_vertices)) {
      cerr << "Input error.  Too few scalar values in option -scalar." << endl;
      cerr << "  In dimension " << dimension
           << ", field set_xcube requires " << num_cube_vertices
           << " scalar values." << endl;
      cerr << "  (Use \"\" to form string containing multiple values.)" << endl;
      exit(58);
    }
  }

  if (field_name == "set_x2cube") {
    const int num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
    const int num_facet_vertices = num_cube_vertices/2;
    const int num_scalar_values = field_param.scalar_value.size();

    if ((num_scalar_values > 0) && 
	(num_scalar_values < num_cube_vertices+num_facet_vertices)) {
      cerr << "Input error.  Too few scalar values in option -scalar." << endl;
      cerr << "  In dimension " << dimension
           << ", field set_x2cube requires " 
	   << num_cube_vertices+num_facet_vertices
           << " scalar values." << endl;
      cerr << "  (Use \"\" to form string containing multiple values.)" << endl;
      exit(58);
    }
  }

  if ((grid_axis_size.size() > 1) && (grid_axis_size.size() != dimension)) {
    cerr << "Input error. Incorrect number of values in -axis \"...\" option." << endl;
    cerr << "Either use -axis <N> ";
    if (dimension == DIM2) {
      cerr << "or -axis \"<N0> <N1>\" when dimension is 2." << endl;
    }
    else if (dimension == DIM3) {
      cerr << "or -axis \"<N0> <N1> <N2>\" when dimension is 3." << endl;
    }
    else if (dimension == DIM4) {
      cerr << "or -axis \"<N0> <N1> <N2> <N3>\" when dimension is 4." << endl;
    }
    else {
      cerr << "or -axis \"<N0> <N1> ...<Nk>\" where k=(dimension-1)." << endl;
    }
    exit(58);
  }
}


void check_param(const SCALAR_GRID & grid, const std::string & field_name)
{
  if (field_name == "cannon") {
    const COORD_TYPE dist2center = field_param.dist2ball_center[0];
    const COORD_TYPE dist2near0 = field_param.dist2near0[0];
    const ANGLE_TYPE angle = field_param.angle[0];
    const double angle_radians = (angle*M_PI)/180.0;

    if ((dist2center - dist2near0) < sin(angle_radians)*dist2center) {
      COORD_TYPE min_dist = dist2near0/(1-sin(angle_radians));
      cerr << "Usage error.  Ball center is too close to near clipping plane."
           << endl;
      cerr << "  Distance to ball center should be greater than " 
           << min_dist << endl;
      exit(71);
    }
  }

  if (field_name == "randomP") {
    const SCALAR_TYPE base = field_param.RandomBase();
    const SCALAR_TYPE s = pow(base, field_param.MaxExponent());
    const SCALAR_TYPE s_neg = -s;

    if (s == HUGE_VAL || s_neg == -HUGE_VAL) {
      cerr << "Usage error. Max exponent "
	   << field_param.MaxExponent() << " is too large." << endl;
      cerr << "  " << base << "^" << field_param.MaxExponent()
	   << " is too large for memory." << endl;
      exit(-1);
    }
  }


  if (field_param.set_scalar_list.size() > 0) {
    IJK::ERROR error;

    for (int i = 0; i < field_param.set_scalar_list.size(); i++) {

      if (field_param.set_scalar_list[i].IsTypeSetIsotableIndex()) {

        const CUBE_INDEX_TYPE cube_index = 
          field_param.set_scalar_list[i].cube_index;

        if (!grid.CheckCubeIndex(cube_index, error)) {
          cerr << "Usage error.  Error in -set_cube. Cube "
               << cube_index << " is not a grid cube." << endl;
          exit(70);
        }
      }
      else if (field_param.set_scalar_list[i].IsTypeSetPosVertex()) {

        const VERTEX_INDEX_TYPE vertex_index = 
          field_param.set_scalar_list[i].vertex_index;

        if (!grid.CheckVertexIndex(vertex_index, error)) {
          cerr << "Usage error.  Error in -set_pos_vertex. Vertex "
               << vertex_index << " is out of range." << endl;
          exit(70);
        }
      }
      else if (field_param.set_scalar_list[i].IsTypeSetPosEdge()) {

        const VERTEX_INDEX_TYPE vertex_index = 
          field_param.set_scalar_list[i].vertex_index;

        if (!grid.CheckVertexIndex(vertex_index, error)) {
          cerr << "Usage error.  Error in -set_pos_edge. Edge endpoint "
               << vertex_index << " is out of range." << endl;
          exit(70);
        }
      }
      else if (field_param.set_scalar_list[i].IsTypeSetPosFacet()) {

        const VERTEX_INDEX_TYPE vertex_index = 
          field_param.set_scalar_list[i].vertex_index;

        if (!grid.CheckVertexIndex(vertex_index, error)) {
          cerr << "Usage error.  Error in -set_pos_facet. Facet vertex "
               << vertex_index << " is out of range." << endl;
          exit(70);
        }
      }
      else if (field_param.set_scalar_list[i].IsTypeSetPosFacetDiagonal()) {

        const VERTEX_INDEX_TYPE vertex_index = 
          field_param.set_scalar_list[i].vertex_index;

        if (!grid.CheckVertexIndex(vertex_index, error)) {
          cerr << "Usage error.  Error in -set_pos_facet_diagonal. Facet vertex "
               << vertex_index << " is out of range." << endl;
          exit(70);
        }
      }
    }
  }
}

void check_field()
{
  int ifield = field_param.FieldIndex();

  // Currently, nothing to check.
}

void check_field_options()
{
  int ifield = field_param.FieldIndex();

  if (output_gradients) {

    if (!field_info[ifield].IsGradientImplemented()) {
      cerr << "Usage error.  Option -grad cannot be used with field "
           << field_info[ifield].name << "." << endl;
      usage_error();
    }
  }

  if (field_info[ifield].name == "cube") {
    if (field_param.flag_tilt && field_param.Dimension() != DIM3) {
      cerr << "Error."
           << " Illegal dimension " << field_param.Dimension() 
           << " used with field cube and option tilt."
           << endl;
      cerr << "  Field cube and option tilt can only be used with dimension " 
           << DIM3 << "." << endl;
      exit(50);
    }
  }


  if (field_param.flag_wedge && field_param.Dimension() != DIM3) {
    cerr << "Usage error."
         << " Illegal dimension " << field_param.Dimension() 
         << " used with field " << field_info[ifield].name << "." << endl;
    cerr << "Wedge can only be used with dimension " << DIM3 << "." << endl;
    exit(55);
  }

}

bool check_dimension(const int dimension, const int dimension2,
                     const char * text, IJK::ERROR & error)
{
  if (dimension == dimension2) {
    return(true);
  }
  else {
    error.AddMessage("Programming error.  Illegal dimension ", dimension, ".");
    error.AddMessage
      ("  ", text, " can only be used with dimension ", dimension2, ".");
    return(false);
  }
}


int get_field_index
(const int iarg, const int argc, char **argv)
{
  int ifield(0);

  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error();
  }

  if (!find_geom_info_name(field_info, argv[iarg+1], ifield)) {
    cerr << "Illegal field name: " << argv[iarg+1] << endl;
    exit(51);
  }

  return(ifield);
}

void set_cube_isotable_index
(const int iarg, const int argc, char **argv,
 FIELD_PARAM & field_param)
{
  int isotable_index;
  std::vector<int> cube_index;
  IJK::ERROR error;

  throw_error_on_missing_argument(iarg, argc, argv, 2, error);
  IJK::get_arg_multiple_arguments(iarg, argc, argv, cube_index);

  if (cube_index.size() < 1) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Missing cube indices in string: \"", argv[iarg+1], "\"");
    throw error;
  }

  if (!IJK::string2val(argv[iarg+2], isotable_index)) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Non-integer character in string: ", argv[iarg+2], "");
    throw error;
  }

  for (int i = 0; i < cube_index.size(); i++) {
    field_param.SetCubeIsotableIndex(cube_index[i], isotable_index);
  }
}


void set_pos_vertex
(const int iarg, const int argc, char **argv,
 FIELD_PARAM & field_param)
{
  std::vector<int> vertex_index;
  IJK::ERROR error;

  throw_error_on_missing_argument(iarg, argc, argv, error);
  IJK::get_arg_multiple_arguments(iarg, argc, argv, vertex_index);

  if (vertex_index.size() < 1) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Missing vertex indices in string: \"", argv[iarg+1], "\"");
    throw error;
  }

  for (int i = 0; i < vertex_index.size(); i++) {
    field_param.SetScalarPosVertex(vertex_index[i]);
  }
}


void set_pos_edge
(const int iarg, const int argc, char **argv,
 FIELD_PARAM & field_param)
{
  std::vector<int> vertex_index;
  int direction;
  IJK::ERROR error;

  throw_error_on_missing_argument(iarg, argc, argv, 2, error);
  IJK::get_arg_multiple_arguments(iarg, argc, argv, vertex_index);

  if (vertex_index.size() < 1) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Missing vertex indices in string: \"", argv[iarg+1], "\"");
    throw error;
  }

  if (!IJK::string2val(argv[iarg+2], direction)) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Non-integer character in string: ", argv[iarg+2], "");
    throw error;
  }

  for (int i = 0; i < vertex_index.size(); i++) {
    field_param.SetScalarPosEdge(vertex_index[i], direction);
  }
}


void set_pos_facet
(const int iarg, const int argc, char **argv,
 FIELD_PARAM & field_param)
{
  std::vector<int> vertex_index;
  int direction;
  IJK::ERROR error;

  throw_error_on_missing_argument(iarg, argc, argv, 2, error);
  IJK::get_arg_multiple_arguments(iarg, argc, argv, vertex_index);

  if (vertex_index.size() < 1) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Missing vertex indices in string: \"", argv[iarg+1], "\"");
    throw error;
  }

  if (!IJK::string2val(argv[iarg+2], direction)) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Non-integer character in string: ", argv[iarg+2], "");
    throw error;
  }

  for (int i = 0; i < vertex_index.size(); i++) {
    field_param.SetScalarPosFacet(vertex_index[i], direction);
  }
}


void set_pos_facet_diagonal
(const int iarg, const int argc, char **argv,
 FIELD_PARAM & field_param)
{
  std::vector<int> vertex_index;
  int direction;
  int jth_vertex;
  IJK::ERROR error;

  throw_error_on_missing_argument(iarg, argc, argv, 3, error);
  IJK::get_arg_multiple_arguments(iarg, argc, argv, vertex_index);

  if (vertex_index.size() < 1) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Missing vertex indices in string: \"", argv[iarg+1], "\"");
    throw error;
  }

  if (!IJK::string2val(argv[iarg+2], direction)) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Non-integer character in string: ", argv[iarg+2], "");
    throw error;
  }

  if (!IJK::string2val(argv[iarg+3], jth_vertex)) {
    error.AddMessage
      ("Usage error.  Error in argument for option: ", argv[iarg], "");
    error.AddMessage
      ("Non-integer character in string: ", argv[iarg+3], "");
    throw error;
  }

  for (int i = 0; i < vertex_index.size(); i++) {
    field_param.SetScalarPosFacetDiagonal
      (vertex_index[i], direction, jth_vertex);
  }
}

void parse_command_line(const int argc, char **argv)
{
  using IJK::get_arg_int;
  using IJK::get_arg_float;
  using IJK::get_arg2_float;
  using IJK::get_arg_multiple_arguments;

  if (argc == 2 && string(argv[1]) == "-version") {
    cout << "Version: " << ProgramVersion() << endl;
    exit(0);
  }

  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {

    string s = argv[iarg];

    if (s == "-dim") {
      int dimension = get_arg_int(iarg, argc, argv);
      field_param.SetDimension(dimension);
      iarg++;
    }
    else if (s == "-asize") {
      grid_axis_size.clear();
      get_arg_multiple_arguments(iarg, argc, argv, grid_axis_size);
      iarg++;
    }
    else if (s == "-field") {
      int ifield = get_field_index(iarg, argc, argv);
      field_param.SetFieldIndex(ifield);
      iarg++;
    }
    else if (s == "-grad") {
      output_gradients = true;
    }
    else if (s == "-disc_zero") {
      field_param.gradient_discontinuity_zero = true;
    }
    else if (s == "-disc_select") {
      field_param.gradient_discontinuity_zero = false;
    }
    else if (s == "-tilt") {
      field_param.flag_tilt = true;
    }
    else if (s == "-no_tilt") {
      field_param.flag_tilt = false;
    }
    else if (s == "-multi") {
      field_param.flag_multi_centers = true;
    }
    else if (s == "-single") {
      field_param.flag_multi_centers = false;
    }
    else if (s == "-two") {
      field_param.flag_multi_centers = true;
      field_param.num_objects.Set(2);
    }
    else if (s == "-n") {
      int num_obj = get_arg_int(iarg, argc, argv);
      field_param.num_objects.Set(num_obj);
      iarg++;
    }
    else if (s == "-stack") {
      field_param.flag_stack = true;
      field_param.flag_multi_centers = true;
    }
    else if (s == "-flange") {
      field_param.flag_flange = true;
    }
    else if (s == "-flange_wh") {
      field_param.flag_flange = true;

      float w, h;
      get_arg2_float(iarg, argc, argv, w, h);
      field_param.flange_width.push_back(w);
      field_param.flange_height.push_back(h);
      iarg = iarg+2;
    }
    else if (s == "-wedge") {
      field_param.flag_wedge = true;
    }
    else if (s == "-center") {
      get_arg_multiple_arguments(iarg, argc, argv, field_param.center);
      iarg++;
    }
    else if (s == "-grid_center") {
      field_param.flag_grid_center = true;
    }
    else if (s == "-randompos") {
      field_param.flag_random_centers = true;
      int seed = get_arg_int(iarg, argc, argv);
      field_param.randompos_seed.Set(seed);
      iarg++;
    }
    else if (s == "-randomdir") {
      field_param.flag_random_directions = true;
      int seed = get_arg_int(iarg, argc, argv);
      field_param.randomdir_seed.Set(seed);
      iarg++;
    }
    else if (s == "-base") {
      int base = get_arg_int(iarg, argc, argv);
      field_param.random_base.Set(base);
      iarg++;
    }
    else if (s == "-dir") {
      get_arg_multiple_arguments(iarg, argc, argv, field_param.direction);
      iarg++;
    }
    else if (s == "-side_dir") {
      get_arg_multiple_arguments(iarg, argc, argv, field_param.side_direction);
      iarg++;
    }
    else if (s == "-normal") {
      get_arg_multiple_arguments(iarg, argc, argv, field_param.normal);
      iarg++;
    }
    else if (s == "-rot_dir") {
      get_arg_multiple_arguments
        (iarg, argc, argv, field_param.rotation_direction);
      iarg++;
    }
    else if (s == "-spacing") {
      get_arg_multiple_arguments(iarg, argc, argv, field_param.spacing);
      iarg++;
    }
    else if (s == "-translate") {
      get_arg_multiple_arguments(iarg, argc, argv, field_param.translate);
      iarg++;
    }
    else if (s == "-radius") {
      float x = get_arg_float(iarg, argc, argv);
      field_param.radius.push_back(x);
      iarg++;
    }
    else if (s == "-length_diff") {
      float x = get_arg_float(iarg, argc, argv);
      field_param.length_difference.push_back(x);
      iarg++;
    }
    else if (s == "-wedge_angle") {
      ANGLE_TYPE angle = get_arg_float(iarg, argc, argv);
      field_param.wedge_angle.push_back(angle);
      field_param.flag_wedge = true;
      iarg++;
    }
    else if (s == "-angle") {
      ANGLE_TYPE angle = get_arg_float(iarg, argc, argv);
      field_param.angle.push_back(angle);
      iarg++;
    }
    else if (s == "-wedge_isovalue") {
      SCALAR_TYPE wedge_isovalue = get_arg_float(iarg, argc, argv);
      field_param.wedge_isovalue.push_back(wedge_isovalue);
      iarg++;
    }
    else if (s == "-clip0") {
      float x1, x2;
      get_arg2_float(iarg, argc, argv, x1, x2);
      field_param.dist2near0.push_back(x1);
      field_param.dist2far0.push_back(x2);
      iarg = iarg+2;
    }
    else if (s == "-clip_near") {
      float x = get_arg_float(iarg, argc, argv);
      field_param.dist2near0.push_back(x);
      iarg++;
    }
    else if (s == "-dist2center") {
      float x = get_arg_float(iarg, argc, argv);
      field_param.dist2ball_center.push_back(x);
      iarg++;
    }
    else if (s == "-smooth_tip") {
      field_param.flag_smooth_tip = true;
    }
    else if (s == "-seed") {
      int x = get_arg_int(iarg, argc, argv);
      field_param.random_seed.Set(x);
      iarg++;
    }
    else if (s == "-scalar") {
      get_arg_multiple_arguments(iarg, argc, argv, field_param.scalar_value);
      iarg++;
    }
    else if (s == "-maxval") {
      int x = get_arg_int(iarg, argc, argv);
      field_param.maxval.Set(x);
      iarg++;
    }
    else if (s == "-maxexp") {
      int x = get_arg_int(iarg, argc, argv);
      field_param.max_exponent.Set(x);
      iarg++;
    }
    else if (s == "-isotable_index") {
      ISOTABLE_INDEX_TYPE x = get_arg_int(iarg, argc, argv);
      field_param.isotable_index.Set(x);
      iarg++;
    }
    else if (s == "-set_cube") {
      set_cube_isotable_index(iarg, argc, argv, field_param);
      iarg = iarg+2;
    }
    else if (s == "-set_pos_vertex") {
      set_pos_vertex(iarg, argc, argv, field_param);
      iarg = iarg+1;
    } 
    else if (s == "-set_pos_edge") {
      set_pos_edge(iarg, argc, argv, field_param);
      iarg = iarg+2;
    } 
    else if (s == "-set_pos_facet") {
      set_pos_facet(iarg, argc, argv, field_param);
      iarg = iarg+2;
    } 
    else if (s == "-set_pos_facet_diagonal") {
      set_pos_facet_diagonal(iarg, argc, argv, field_param);
      iarg = iarg+3;
    } 
    else if (s == "-bzero") {
      field_param.flag_set_boundary2zero = true;
    }
    else if (s == "-gzip") {
      flag_gzip = true;
    }
    else if (s == "-list") {
      cout << "Valid field names: " << endl;
      out_field_descriptions(field_info);
      exit(20);
    }
    else if (s == "-s") {
      flag_silent = true;
    }
    else if (s == "-version") {
      cout << "Version: " << ProgramVersion() << endl;
    }
    else if (s == "-help") {
      help_msg();
    }
    else {
      cerr << "Usage error.  Illegal parameter: " << s << endl;
      usage_error();
    };
    iarg++;
  };

  if (iarg+1 > argc) { 
    cerr << "Error.  Missing file name." << endl << endl;
    usage_error(); 
  }
  else if (iarg+1 < argc) 
    { usage_error(); }

  ofilename = argv[iarg];
}

void usage_msg()
{
  cerr << "Usage: ijkgenscalar [OPTIONS] <filename>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  [-field <field name>] [-grad]" << endl;
  cerr << "  [-dim <dimension>] [-asize <N> | -asize \"<N0> <N1>...\" ]" << endl;
  cerr << "  [-tilt | -no_tilt]" << endl;
  cerr << "  [-single | -two | -multi | -n <num_obj>] [-stack]" << endl;
  cerr << "  [-center \"<coord>\" | -grid_center | -randompos <S>]" << endl;
  cerr << "  [-dir \"<coord>\" | -randomdir <S>]" << endl;
  cerr << "  [-side_dir \"<coord>\"] [-rot_dir \"<coord>\"]" << endl;
  cerr << "  [-radius <x>] [-length_diff <x>] [-angle <A>]" << endl;
  cerr << "  [-clip0 <N> <F>] [-clip_near <N>] [-dist2center <D>]" << endl;
  cerr << "  [-flange | -flange_wh <W> <H>] [-wedge | wedge_angle <A>]" 
       << endl;
  cerr << "  [-wedge_isovalue <S>]" << endl;
  cerr << "  [-smooth_tip]" << endl;
  cerr << "  [-normal \"<coord>\"] [-translate \"<coord>\"]" << endl;
  cerr << "  [-disc_zero] [-disc_select]" << endl;
  cerr << "  [-spacing \"<space_coord>\"]" << endl;
  cerr << "  [-seed <S>] [-maxval <M>] [-base <B>] [-maxexp <M>]" << endl; 
  cerr << "  [-scalar <S> | -scalar \"<scalar_values>\"]" << endl;
  cerr << "  [-bzero]" << endl;
  cerr << "  [-isotable_index <I>]" << endl;
  cerr << "  [-set_cube <C> <I> | -set_cube \"<cube_list\" <I> ]" << endl;
  cerr << "  [-set_pos_vertex <V>] [-set_pos_edge <V> <dir>]" << endl;
  cerr << "  [-set_pos_facet <V> <orth_dir>]" << endl;
  cerr << "  [-set_pos_facet_diagonal <V> <orth_dir> <j>]" << endl;
  cerr << "  [-gzip] [-s] [-list] [-version] [-help]" << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help_msg()
{
  cerr << "Usage: ijkgenscalar [OPTIONS] <filename>" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -field <field name>: Set field."
       << "  (Use -list to list fields.)" << endl;
  cerr << "  -grad:      Write gradients to file." << endl;
  cerr << "  -dim <D>:   Set dimension to <D>." << endl;
  cerr << "  -asize <N>: Set number of vertices along each grid axis to <N>."
       << endl;
  cerr << "  -asize \"<N0> <N1>...\": Set number of vertices along grid axis 0 to <N0>,"
       << endl
       << "              along grid axis 1 to <N1>, etc." << endl;
  cerr << "  -tilt:      Use directions to tilt objects." << endl;
  cerr << "  -no_tilt:   Objects are not tilted." << endl; 
  cerr << "  -single:    Single object (default)." << endl;
  cerr << "  -two:       Two objects." << endl;
  cerr << "  -multi:     Multiple objects." << endl;
  cerr << "  -n <N>:     Number of objects is <N>." << endl;
  cerr << "  -stack:     Stack of translated objects." << endl;
  cerr << "  -center \"<coord>\": Set center to given coordinates." << endl;
  cerr << "  -grid_center: Set center to center of grid." << endl;
  cerr << "  -randompos <S>: Assign random positions to centers." << endl
       << "      Use <S> as random seed for positions." << endl;
  cerr << "  -dir \"<coord>\": Set direction to given coordinates." << endl;
  cerr << "  -randomdir <S>: Assign random directions." << endl
       << "      Use <S> as random seed for positions." << endl;
  cerr << "  -start_dir \"<coord>\": Set start direction to given coordinates." 
       << endl;
  cerr << "  -side_dir \"<coord>\": Set side direction to given coordinates." 
       << endl;
  cerr << "  -rot_dir \"<coord>\": Set rotation direction to given coordinates."
       << endl;
  cerr << "  -radius <x>: Set radius to <x>." << endl;
  cerr << "  -length_diff <x>:  Set (length-width) or (height-width) to <x>."
       << endl;
  cerr << "            Negative <x> sets width greater than length or height."
       << endl;
  cerr << "  -angle <A>:  Set cone, frustrum or edge angle to <A>." << endl;
  cerr << "  -clip0 <N> <F>: Set distance from isovalue 0 frustrum apex" << endl
       << "                    to near clipping plane to <N>." << endl
       << "                  Set distance to far clipping plane to <F>." 
       << endl;
  cerr << "  -clip_near <N>: Set distance from isovalue 0 cannon apex" << endl
       << "                    to near clipping plane to <N>." << endl;
  cerr << "  -dist2center <D>: Set distance from cannon apex to ball center."
       << endl;
  cerr << "  -flange: Add flange to object." << endl;
  cerr << "  -flange_wh <W> <H>: Set flange width to <W> and height to <H>."
       << endl;
  cerr << "  -wedge: Intersect object with wedge." << endl;
  cerr << "  -wedge_angle <A>: Set wedge angle to <A>." << endl;
  cerr << "  -wedge_isovalue <S>: Set scalar value where wedge passes" << endl
       << "                       through the center to <S>." << endl;
  cerr << "     Note: This determines how the wedge intersects objects."
       << endl;
  cerr << "  -smooth_tip : Smooth the cone tip." << endl;
  cerr << "  -normal \"<coord>\": Set plane normals to given coordinates." 
       << endl;
  cerr << "  -translate \"<coord>\": Set stack translation vector to given coordinates." 
       << endl;
  cerr << "  -disc_zero:  Set gradient to zero at gradient discontinuity." 
       << endl;
  cerr << "  -disc_select: Select gradient from one incident surface" << endl;
  cerr << "                at gradient discontinuity." << endl;
  cerr << "  -spacing \"<space_coord>\": Set spacing between grid vertices"
       << endl
       << "                 along each axis." << endl;
  cerr << "  -seed <S>:   Set random seed to <S>." << endl;
  cerr << "  -scalar <S>: Set a single scalar value." << endl;
  cerr << "               Used with field constant." << endl;
  cerr << "  -scalar \"<scalar_values>\": Set list of scalar values." << endl;
  cerr << "               Used with field set_xcube or set_x2cube or setv_interior." << endl;
  cerr << "  -maxval <M>: Set maximum random scalar value to <M>." << endl;
  cerr << "  -base <B>:   Set base for randomP number generation to <B>." 
       << endl
       << "               Generate random int x and +/- sign S," << endl
       << "               and then compute S*B^x."
       << endl;
  cerr << "  -maxexp <M>: Set maximum exponent for randomP number generation to <E>." << endl
       << "               Generate random int x in range [0,<M>] and" << endl
       << "               random +/- sign S, and then compute S*B^x."
       << endl;
  cerr << "  -bzero:      Set values at all boundary vertices to zero." << endl;
  cerr << "  -isotable_index <I>: Set isosurface table index to <I>." << endl
       << "               Used only with field isotable_entry." << endl;
  cerr << "  -set_cube <C> <I>: Set cube <C> to isosurface table configuration <I>." << endl
       << "               Used with fields set or randomint." << endl;
  cerr << "  -set_cube \"<cube_list\" <I>: Set all cubes in cube_list" << endl
       << "                 to isosurface configuration <I>." << endl;
  cerr << "  -set_pos_vertex <V>: Set vertex <V> to be positive." << endl;
  cerr << "  -set_pos_edge <V> <dir>: Set all vertices of edge to be positive."
       << endl;
  cerr << "     Edge has lower endpoint <V> and direction <dir>." << endl;
  cerr << "  -set_facet <V> <orth_dir>: Set all vertices of facet" << endl;
  cerr << "     to be positive.  Facet has lower/leftmost vertex <V> and"
       << endl
       << "     orthogonal direction <orth_dir>." << endl;
  cerr << "  -set_pos_facet_diagonal <V> <orth_dir> <j>:" << endl;
  cerr << "     Set vertices of facet diagonal to be positive." << endl;
  cerr << "     Facet has lower/leftmost vertex <V> and orthogonal direction <orth_dir>." << endl;
  cerr << "     Diagonal contains j'th vertex of facet." << endl;
  cerr << "  -gzip:       Compress output in gzip format." << endl;
  cerr << "  -list:       List fields." << endl;
  cerr << "  -s:          Silent mode." << endl;
  cerr << "  -version:    Print program version." << endl;
  cerr << "  -help:       Print this help message." << endl;

  exit(20);
}
