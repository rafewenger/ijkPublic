/*!
 *  @file ijkMCtable_xitIO.cxx
 *  @brief I/O routines for .xit (XML isosurface table) file
 *  - .xit is an xml format 
 *  - Version 0.6
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2024 Rephael Wenger

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

#include <string.h>

#include <iostream>
#include <sstream>
#include <ctime>

#include "expat.h"

#include "ijkMCtable_xitIO.h"

using namespace IJK;
using namespace IJKXIO;

//**************************************************
// read_xit
//**************************************************

namespace {
  // global error class
  ERROR read_error;
}

// types
typedef enum { START_TAG, END_TAG } XML_TAG_TYPE;
typedef enum { NO_ELEMENT,
               ISOTABLE, VERSION, CREATION_DATE, TABLE_TYPE,
               DIMENSION, 
               POLY, POLY_SHAPE,
               VERTICES, NUM_VERTICES, C,
               EDGES, NUM_EDGES, V, W,
               FACETS, NUM_FACETS, F,
               ISO_VERTICES, IN_E, IN_V, IN_F, L,
               TABLE, ENCODING, TABLE_PROPERTIES,
               POLY_VERTEX_LABEL_TYPE, SEPARATION_TYPE,
               TRIANGULATION_TYPE, SEPARATE_OPPOSITE,
               ISO_POLY_ORIENTATION,
               NUM_ENTRIES, S,
               UNKNOWN_ELEMENT } XML_ELEMENT;

/// IO user data
typedef struct {
  XML_TAG_TYPE tag_type;
  XML_ELEMENT element;
  bool flag_eof;
  /* NOT YET USED
     string char_data; */
  } USER_DATA;

const char * xml_strings[] = 
  { "No element", "isotable",
    "version", "creationDate", "tableType",
    "dimension",
    "poly", "polyShape",
    "vertices", "numVertices", "c",
    "edges", "numEdges", "v", "w",
    "facets", "numFacets", "f",
    "isoVertices", "inE", "inV", "inF", "L",
    "table", "encoding", "tableProperties",
    "polyVertexLabelType", "separationType",
    "triangulationType", "separateOpposite", 
    "isoPolyOrientation",
    "numEntries", "s",
    "Unknown element" };

// local functions
namespace {

  // *** Parsing functions ***
  
  void read_tag
  (std::istream & in, XML_Parser parser, USER_DATA & user_data);
  void read_start_tag(std::istream & in, XML_Parser parser, 
                      const XML_ELEMENT element, USER_DATA & user_data);
  void read_end_tag(std::istream & in, XML_Parser parser, 
                    const XML_ELEMENT element, USER_DATA & user_data);
  std::string read_string(std::istream & in, XML_Parser parser, 
                          USER_DATA & user_data);
  void XMLCALL startElement(void *data, const char *name, 
                            const char **attr);
  void XMLCALL endElement(void *data, const char *ename);
  void initUserData(USER_DATA & user_data);

  
  // *** Read data functions ***

  /// @brief Read all data from version 1.0 xit file.
  /// - Starts at reading dimension.
  void read_xit_V1_data
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  /// @brief Read all data from version 2.0 xit file.
  /// - Starts at reading dimension.
  void read_xit_V2_data
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_xit_version
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   std::string & version_str);
  
  XIT_VERSION_TYPE determine_xit_version
  (const std::string & version_str);

  void read_table_type
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  
  void read_dimension
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_poly_shape
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_polytope_vertices
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_polytope_edges
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_polytope_facets
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  
  void read_isosurface_vertex
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table, const int isov);

  void read_isosurface_vertices
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_table_encoding
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_table_properties
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  void read_mcube_lookup_table
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  
  void XMLCALL startPolyElement(void *data, const char *name, 
                                const char **attr);
  void XMLCALL startIsoVerticesElement(void *data, const char *name,
                                       const char **attr);
  void XMLCALL endIsoVerticesElement(void *data, const char *ename);
  void XMLCALL startTableElement(void *data, const char *name, 
                                 const char **attr);
  void XMLCALL startTableEntriesElement
  (void *data, const char *name, const char **attr);
  void XMLCALL startTablePropertiesElement
  (void *data, const char *name, const char **attr);
  void XMLCALL endTablePropertiesElement(void *data, const char *name);
}


// Read xit file.
void IJKXIO::read_xit
(std::istream & in, IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
{
  XML_Parser parser = NULL;
  USER_DATA user_data;

  read_error.AddProcMessage("read_xit");

  try {

    initUserData(user_data);

    parser = XML_ParserCreate(NULL);

    if (!parser) 
      throw read_error("Couldn't allocate memory for XML parser.");

    XML_SetElementHandler(parser, startElement, endElement);
    XML_SetUserData(parser, (void *) &user_data);

    // read isotable start tag
    read_start_tag(in, parser, ISOTABLE, user_data);

    // read in version, creationDate
    std::string version_str;
    read_xit_version(in, parser, user_data, version_str);
    read_start_tag(in, parser, CREATION_DATE, user_data);

    const XIT_VERSION_TYPE xit_version = 
      determine_xit_version(version_str);

    if (xit_version == XIT_VERSION_1_0) {
      read_xit_V1_data(in, parser, user_data, table);
    }
    else if (xit_version == XIT_VERSION_2_0) {
      read_xit_V2_data(in, parser, user_data, table);
    }
    else {
      const std::string msg = "Unknown .xit version: " + version_str;
      throw read_error(msg);
    }

  } catch(...) {
    if (parser != NULL)
      XML_ParserFree(parser);

    // re-throw the exception
    throw;
  }

  XML_ParserFree(parser);
}


// Read (old) version 1.0 xit file.
void IJKXIO::read_xit_V1
(std::istream & in, IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
{
  XML_Parser parser = NULL;
  USER_DATA user_data;

  read_error.AddProcMessage("read_xit");

  try {

    initUserData(user_data);

    parser = XML_ParserCreate(NULL);

    if (!parser) 
      throw read_error("Couldn't allocate memory for XML parser.");

    XML_SetElementHandler(parser, startElement, endElement);
    XML_SetUserData(parser, (void *) &user_data);

    // read in version, creationDate
    read_start_tag(in, parser, VERSION, user_data);
    read_start_tag(in, parser, CREATION_DATE, user_data);

    read_xit_V1_data(in, parser, user_data, table);
  } catch(...) {
    if (parser != NULL)
      XML_ParserFree(parser);

    // re-throw the exception
    throw;
  }

  XML_ParserFree(parser);
}


namespace {

  // *** Parsing functions ***

  void read_tag
  (std::istream & in, XML_Parser parser, USER_DATA & user_data)
  {
    const int BUFFSIZE = 1024;
    static char buffer[BUFFSIZE];

    char s;
    int done = 0;
    int length = 0;
    while (length < BUFFSIZE && in.get(s) && s != '\n' && s != '>') {
      buffer[length] = s;
      length++;
    };
    if (length < BUFFSIZE && s == '>' && in.good()) {
      buffer[length] = s;
      length++;
    };

    if (in.eof()) {
      user_data.flag_eof = 1;
      done = 1;
    }

    if (XML_Parse(parser, buffer, length, done) == XML_STATUS_ERROR) {
      std::string msg = "Parse error " + 
        std::string(XML_ErrorString(XML_GetErrorCode(parser)));
      throw read_error(msg);
    }
  }

  void read_start_tag(std::istream & in, XML_Parser parser, 
                      const XML_ELEMENT element, USER_DATA & user_data)
  {
    user_data.element = NO_ELEMENT;
    while(in.good() &&
          (user_data.element != element || 
           user_data.tag_type != START_TAG)) {

      read_tag(in, parser, user_data);
    }

    if (in.eof() || in.fail()) {
      std::string msg = 
        "Missing start tag: <" + std::string(xml_strings[element]) + ">";
      throw read_error(msg);
    };
  }

  void read_end_tag(std::istream & in, XML_Parser parser, 
                    const XML_ELEMENT element, USER_DATA & user_data)
  {
    user_data.element = NO_ELEMENT;
    while(in.good() &&
          (user_data.element != element || user_data.tag_type != END_TAG)) {
      read_tag(in, parser, user_data);
    }

    if (in.eof() || in.fail()) {
      std::string msg = 
        "Missing end tag: <" + std::string(xml_strings[element]) + ">";
      throw read_error(msg);
    };
  }

  std::string read_string(std::istream & in, XML_Parser parser, 
                          USER_DATA & user_data)
    // read string, starting with first non-blank character and ending
    // with next blank character
  {
    std::string s;
    char c;

    // read white space before string
    in.get(c);
    while (in.good() && isspace(c)) { in.get(c); }

    // read string
    while (in.good() && !isspace(c) && c != '<') {
      s = s + c;
      in.get(c);
    }
    if (in.good() && c == '<')
      in.putback(c);

    return(s);
  }

  void XMLCALL startElement
  (void *data, const char *name, const char **attr)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = START_TAG;

    if (strcmp(name, xml_strings[VERSION]) == 0) {
      user_data->element = VERSION;
    }
    else if (strcmp(name, xml_strings[ISOTABLE]) == 0) {
      user_data->element = ISOTABLE;
    }
    else if (strcmp(name, xml_strings[CREATION_DATE]) == 0) {
      user_data->element = CREATION_DATE;
    }
    else if (strcmp(name, xml_strings[TABLE_TYPE]) == 0) {
      user_data->element = TABLE_TYPE;
    }
    else if (strcmp(name, xml_strings[DIMENSION]) == 0) {
      user_data->element = DIMENSION;
    }
    else {
      user_data->element = UNKNOWN_ELEMENT;
    };
  }

  void XMLCALL endElement(void *data, const char *ename)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = END_TAG;
  }

  /* NOT YET USED
     static void XMLCALL
     charDataHandler(void *data, const XML_Char *s, int length)
     {
     USER_DATA * user_data = (USER_DATA *) data;

     if (user_data == NULL) return;

     user_data->char_data += s;
     }
  */

  void initUserData(USER_DATA & user_data)
  {
    user_data.tag_type = START_TAG;
    user_data.element = NO_ELEMENT;
    user_data.flag_eof = false;
  }

  
  // *** Read data functions ***

  /// Read all data from version 2.0 xit file.
  void read_xit_V2_data
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    // read table type
    read_table_type(in, parser, user_data, table);

    // read dimension
    read_dimension(in, parser, user_data, table);
                   
    // read in polytope description
    XML_SetElementHandler(parser, startPolyElement, endElement);
    read_start_tag(in, parser, POLY, user_data);
    read_poly_shape(in, parser, user_data, table);
    read_polytope_vertices(in, parser, user_data, table);    
    read_polytope_edges(in, parser, user_data, table);
    read_polytope_facets(in, parser, user_data, table);

    // read isosurface vertices
    read_isosurface_vertices(in, parser, user_data, table);
      
    // read isosurface lookup table
    XML_SetElementHandler(parser, startTableElement, endElement);
    read_start_tag(in, parser, TABLE, user_data);
    read_table_encoding(in, parser, user_data, table);
    read_table_properties(in, parser, user_data, table);
    read_mcube_lookup_table(in, parser, user_data, table);
  }

  
  /// Read all data from version 1.0 xit file.
  void read_xit_V1_data
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    // read dimension
    read_dimension(in, parser, user_data, table);
                   
    // read in polytope description
    XML_SetElementHandler(parser, startPolyElement, endElement);
    read_start_tag(in, parser, POLY, user_data);
    read_polytope_vertices(in, parser, user_data, table);    
    read_polytope_edges(in, parser, user_data, table);
    read_polytope_facets(in, parser, user_data, table);

    // read isosurface vertices
    read_isosurface_vertices(in, parser, user_data, table);
      
    // read isosurface lookup table
    XML_SetElementHandler(parser, startTableElement, endElement);
    read_start_tag(in, parser, TABLE, user_data);
    read_table_encoding(in, parser, user_data, table);
    read_mcube_lookup_table(in, parser, user_data, table);
  }

  void read_xit_version
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   std::string & version_str)
  {
    read_start_tag(in, parser, VERSION, user_data);
    version_str = read_string(in, parser, user_data);
  }

  XIT_VERSION_TYPE determine_xit_version
  (const std::string & version_str)
  {
    if ((version_str == "1.0") ||
        (version_str == "1")) {
      return XIT_VERSION_1_0;
    }
    else if ((version_str == "2.0") ||
	     (version_str == "2")) {
      return XIT_VERSION_2_0;
    }
    else {
      const std::string version_prefix =
        version_str.substr(0, version_str.find('.'));
      if (version_prefix == "2") {
        return XIT_VERSION_2_x;
      }
      else {
        return UNKNOWN_XIT_VERSION;
      }
    }
  }


  void read_table_type
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using IJKMCUBE_TABLE::LOOKUP_TABLE_TYPE;

    read_start_tag(in, parser, TABLE_TYPE, user_data);
    const std::string table_type_str = read_string(in, parser, user_data);
    const LOOKUP_TABLE_TYPE table_type =
      table.Properties().TableType(table_type_str);
    table.SetTableType(table_type);
  }
  
  
  void read_dimension
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using IJKMCUBE_TABLE::ISOSURFACE;
    using IJKMCUBE_TABLE::INTERVAL_VOLUME;

    int polyDimension = 0;
    int surfaceDimension = 0;

    read_start_tag(in, parser, DIMENSION, user_data);
    in >> polyDimension;
    in >> surfaceDimension;
    if (in.fail()) throw read_error("Error reading dimension.");
    table.SetDimension(polyDimension);
    table.SetSimplexDimension(surfaceDimension);

    if (polyDimension == surfaceDimension) 
      { table.SetTableType(INTERVAL_VOLUME); }
    else if (polyDimension == (surfaceDimension+1))
      { table.SetTableType(ISOSURFACE); }
  }


  void read_poly_shape
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using IJKMCUBE_TABLE::POLYTOPE_SHAPE;

    read_start_tag(in, parser, POLY_SHAPE, user_data);
    const std::string shape_str = read_string(in, parser, user_data);
    const POLYTOPE_SHAPE shape = table.Polytope().Shape(shape_str);
    table.SetPolyShape(shape);
  }

  
  void read_polytope_vertices
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    const int polyDimension = table.Dimension();
    int numVertices = 0;

    read_start_tag(in, parser, VERTICES, user_data);
    read_start_tag(in, parser, NUM_VERTICES, user_data);
    in >> numVertices;
    if (in.fail()) throw read_error("Error reading number of vertices.");
    table.SetNumPolyVertices(numVertices);

    for (int i = 0; i < numVertices; i++) {
      read_start_tag(in, parser, C, user_data);
      if (in.fail()) {
        read_error.AddMessage
          ("Error reading vertex coordinates of polyhedron vertex ", i, ".");
        throw read_error;
      };
    
      int coord;
      for (int d = 0; d < polyDimension; d++) {
        in >> coord;
        if (in.fail()) {
          read_error.AddMessage
            ("Error reading polyhedron vertex ", i, ", coordinate ", d, ".");
          throw read_error;
        }
        table.SetPolyVertexCoord(i, d, coord);
      }
    }
  }

  void read_polytope_edges
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    const int numVertices = table.Polytope().NumVertices();
    int numEdges;
    std::ostringstream msg;
    
    read_start_tag(in, parser, EDGES, user_data);
    read_start_tag(in, parser, NUM_EDGES, user_data);
    in >> numEdges;
    if (in.fail()) throw read_error("Error reading number of edges.");
    table.SetNumPolyEdges(numEdges);

    for (int i = 0; i < numEdges; i++) {
      read_start_tag(in, parser, V, user_data);
      int end0, end1;
      in >> end0;
      in >> end1;
      if (in.fail()) {
        read_error.AddMessage("Error reading polyhedron edge ", i, ".");
        throw read_error;
      };
      if (end0 < 0 || end1 < 0 || end0 >= numVertices || 
	  end1 >= numVertices) {
        msg.str("");
        msg << "Error reading polyhedron edge " << i 
            << ".  Illegal vertex indices: " << end0 << " " << end1 << "."
            << std::ends;
        throw read_error(msg.str());
      };
      table.SetPolyEdge(i, end0, end1);
    }
  }

  void read_polytope_facets
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    const int numVertices = table.Polytope().NumVertices();
    int numFacets;
    std::ostringstream msg;

    read_start_tag(in, parser, FACETS, user_data);
    read_start_tag(in, parser, NUM_FACETS, user_data);
    in >> numFacets;
    if (in.fail()) throw read_error("Error reading number of facets.");
    table.SetNumPolyFacets(numFacets);

    for (int i = 0; i < numFacets; i++) {
      read_start_tag(in, parser, F, user_data);
      int numv;
      in >> numv;
      if (in.fail()) {
        read_error.AddMessage("Error reading polyhedron facet ", i, ".");
        throw read_error;
      }

      table.SetPolyNumFacetVertices(i, numv);
    
      for (int j = 0; j < numv; j++) {
        int v;
        in >> v;
        if (in.fail()) {
          msg.str("");
          msg << "Error reading polyhedron facet " << i << ", vertex " 
              << j << std::ends;
          throw read_error(msg.str());
        }
        if (v < 0 || v >= numVertices) {
          msg.str("");
          msg << "Error reading polyhedron facet " << i << ", vertex " 
              << j << ".  Illegal vertex index: " << v << "." << std::ends;
          throw read_error(msg.str());
        }
        table.SetPolyFacetVertex(i, j, v);
      }
    }
  }
  
  void read_isosurface_vertex
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table, const int isov)
  {
    user_data.element = NO_ELEMENT;
    read_start_tag(in, parser, W, user_data);

    while (in.good() && 
           ((user_data.element != IN_V && user_data.element != IN_E &&
             user_data.element != IN_F && user_data.element != C) ||
            user_data.tag_type != START_TAG)) {
      read_tag(in, parser, user_data);
    }

    switch(user_data.element) {
    case IN_V:
      int iv;
      in >> iv;
      table.SetIsoVertexType
        (isov, IJKMCUBE_TABLE::ISOSURFACE_VERTEX::VERTEX);
      table.SetIsoVertexFace(isov, iv);
      read_end_tag(in, parser, IN_V, user_data);
      break;

    case IN_E:
      int ie;
      in >> ie;
      table.SetIsoVertexType
        (isov, IJKMCUBE_TABLE::ISOSURFACE_VERTEX::EDGE);
      table.SetIsoVertexFace(isov, ie);
      read_end_tag(in, parser, IN_E, user_data);
      break;

    case IN_F:
      int jf;
      in >> jf;
      table.SetIsoVertexType(isov, IJKMCUBE_TABLE::ISOSURFACE_VERTEX::FACET);
      table.SetIsoVertexFace(isov, jf);
      read_end_tag(in, parser, IN_F, user_data);
      break;

    case C:
      int num_coord = table.Dimension();
      table.SetIsoVertexNumCoord(isov, num_coord);
      for (int d = 0; d < num_coord; d++) {
        IJKMCUBE_TABLE::ISOSURFACE_VERTEX::COORD_TYPE coord;
        in >> coord;
        table.SetIsoVertexCoord(isov, d, coord);
      }
      read_end_tag(in, parser, C, user_data);
      break;
    }

    read_tag(in, parser, user_data);
    while (in.good() && 
           ((user_data.element != L || user_data.tag_type != START_TAG) &&
            (user_data.element != W || user_data.tag_type != END_TAG))) {
      read_tag(in, parser, user_data);
    }

    if (in.eof() || in.fail()) {
      std::string msg = "Missing end tag for isosurface vertex tag <w>.";
      throw read_error(msg);
    };

    if (user_data.element == L && user_data.tag_type == START_TAG) {
      std::string label = read_string(in, parser, user_data);
      table.SetIsoVertexLabel(isov, label);
    }
  }

  void read_isosurface_vertices
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    int numIsosurfaceVertices;

    XML_SetElementHandler
      (parser,  startIsoVerticesElement, endIsoVerticesElement);
    read_start_tag(in, parser, ISO_VERTICES, user_data);
    read_start_tag(in, parser, NUM_VERTICES, user_data);
    in >> numIsosurfaceVertices;
    table.SetNumIsosurfaceVertices(numIsosurfaceVertices);
    for (int iv = 0; iv < numIsosurfaceVertices; iv++) {
      read_isosurface_vertex(in, parser, user_data, table, iv);
    }
  }


  void read_table_encoding
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    read_start_tag(in, parser, ENCODING, user_data);
    const std::string encoding_str = read_string(in, parser, user_data);
    table.SetEncoding(encoding_str);
  }


  /// @brief Read poly vertex label type.
  /// - Already read POLY_VERTEX_LABEL_TYPE start tag.
  void read_poly_vertex_label_type
  (std::istream & in, XML_Parser parser, USER_DATA & user_data)
  {
    const std::string s = read_string(in, parser, user_data);

  }

  void read_table_properties
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using IJKMCUBE_TABLE::GRID_VERTEX_LABEL_TYPE;
    using IJKMCUBE_TABLE::ISOSURFACE_SEPARATION_TYPE;
    using IJKMCUBE_TABLE::ISOSURFACE_TRIANGULATION_TYPE;

    XML_SetElementHandler
      (parser,  startTablePropertiesElement, endTablePropertiesElement);

    read_start_tag(in, parser, TABLE_PROPERTIES, user_data);

    // Clear user_data.element.
    user_data.element = NO_ELEMENT;
    read_tag(in, parser, user_data);

    while(in.good() &&
          (user_data.element != TABLE_PROPERTIES)) {
      if (user_data.element == POLY_VERTEX_LABEL_TYPE) {
        const std::string label_str = 
          read_string(in, parser, user_data);
        table.SetGridVertexLabelType(label_str);
      }
      else if (user_data.element == SEPARATION_TYPE) {
        const std::string label_str = 
          read_string(in, parser, user_data);
        table.SetSeparationType(label_str);
      }
      else if (user_data.element == TRIANGULATION_TYPE) {        
        const std::string label_str = 
          read_string(in, parser, user_data);
        table.SetTriangulationType(label_str);
      }
      else if (user_data.element == SEPARATE_OPPOSITE) {
        const std::string label_str = 
          read_string(in, parser, user_data);
        table.SetSeparateOpposite(label_str);
      }
      else if (user_data.element == ISO_POLY_ORIENTATION) {
        const std::string label_str = 
          read_string(in, parser, user_data);
        table.SetIsoPolyOrientation(label_str);
      }
      else {
        // Ignore element.
      }

      read_tag(in, parser, user_data);
    }

    if (in.fail()) {
      read_error.AddMessage("Error reading table properties.");
      read_error.AddMessage("  Missing end tag </tableProperties>.");
      throw read_error;
    }

  }

  void read_mcube_lookup_table
  (std::istream & in, XML_Parser parser, USER_DATA & user_data,
   IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    const int surfaceDimension = table.SimplexDimension();
    int numEntries;
    std::ostringstream msg;

    XML_SetElementHandler(parser, startTableEntriesElement, endElement);

    read_start_tag(in, parser, NUM_ENTRIES, user_data);
    in >> numEntries;
    if (in.fail())
      throw read_error("Error reading number of lookup table entries.");
    table.SetNumTableEntries(numEntries);

    if (table.NumTableEntries() != numEntries)
      throw read_error("Error in number of table entries.");

    for (int it = 0; it < numEntries; it++) {
      read_start_tag(in, parser, S, user_data);
      int nums;
      in >> nums;
      if (in.fail()) {
        read_error.AddMessage("Error reading table entry ", it, ".");
        throw read_error;
      }

      table.SetNumSimplices(it, nums);

      for (int js = 0; js < nums; js++) {
        for (int k = 0; k <= surfaceDimension; k++) {
          int ie;
          in >> ie;
          if (in.fail()) {
            msg.str("");
            msg << "Error reading table entry " << it
                << ", simplex " << js << ", "
                << " vertex/edge " << k << std::ends;
            throw read_error(msg.str());
          }

          // NOTE: There is no check on the value of ie

          table.SetSimplexVertex(it, js, k, ie);
        }
      };
    }
  }

  void XMLCALL startPolyElement(void *data, const char *name, 
                                const char **attr)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = START_TAG;

    if (strcmp(name, xml_strings[POLY_SHAPE]) == 0) {
      user_data->element = POLY_SHAPE;
    }
    else if (strcmp(name, xml_strings[VERTICES]) == 0) {
      user_data->element = VERTICES;
    }
    else if (strcmp(name, xml_strings[NUM_VERTICES]) == 0) {
      user_data->element = NUM_VERTICES;
    }
    else if (strcmp(name, xml_strings[C]) == 0) {
      user_data->element = C;
    }
    else if (strcmp(name, xml_strings[EDGES]) == 0) {
      user_data->element = EDGES;
    }
    else if (strcmp(name, xml_strings[NUM_EDGES]) == 0) {
      user_data->element = NUM_EDGES;
    }
    else if (strcmp(name, xml_strings[V]) == 0) {
      user_data->element = V;
    }
    else if (strcmp(name, xml_strings[FACETS]) == 0) {
      user_data->element = FACETS;
    }
    else if (strcmp(name, xml_strings[NUM_FACETS]) == 0) {
      user_data->element = NUM_FACETS;
    }
    else if (strcmp(name, xml_strings[F]) == 0) {
      user_data->element = F;
    }
    else if (strcmp(name, xml_strings[POLY]) == 0) {
      user_data->element = POLY;
    }
    else {
      user_data->element = UNKNOWN_ELEMENT;
    };

  }

  void setIsoVerticesElement(void *data, const char *name)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (strcmp(name, xml_strings[ISO_VERTICES]) == 0) {
      user_data->element = ISO_VERTICES;
    }
    else if (strcmp(name, xml_strings[NUM_VERTICES]) == 0) {
      user_data->element = NUM_VERTICES;
    }
    else if (strcmp(name, xml_strings[W]) == 0) {
      user_data->element = W;
    }
    else if (strcmp(name, xml_strings[IN_V]) == 0) {
      user_data->element = IN_V;
    }
    else if (strcmp(name, xml_strings[IN_E]) == 0) {
      user_data->element = IN_E;
    }
    else if (strcmp(name, xml_strings[IN_F]) == 0) {
      user_data->element = IN_F;
    }  
    else if (strcmp(name, xml_strings[C]) == 0) {
      user_data->element = C;
    }
    else if (strcmp(name, xml_strings[L]) == 0) {
      user_data->element = L;
    }
    else {
      user_data->element = UNKNOWN_ELEMENT;
    };
  }


  void XMLCALL startIsoVerticesElement(void *data, const char *name,
                                       const char **attr)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = START_TAG;

    setIsoVerticesElement(data, name);
  }

  void XMLCALL endIsoVerticesElement(void *data, const char *name)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = END_TAG;

    setIsoVerticesElement(data, name);
  }
						      

  void XMLCALL startTableElement(void *data, const char *name, 
                                 const char **attr)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = START_TAG;
    if (strcmp(name, xml_strings[TABLE]) == 0) {
      user_data->element = TABLE;
    }
    else if (strcmp(name, xml_strings[ENCODING]) == 0) {
      user_data->element = ENCODING;
    }
    else if (strcmp(name, xml_strings[TABLE_PROPERTIES]) == 0) {
      user_data->element = TABLE_PROPERTIES;
    }
    else if (strcmp(name, xml_strings[NUM_ENTRIES]) == 0) {
      user_data->element = NUM_ENTRIES;
    }
    else if (strcmp(name, xml_strings[S]) == 0) {
      user_data->element = S;
    }
    else {
      user_data->element = UNKNOWN_ELEMENT;
    };
  }


  void XMLCALL startTableEntriesElement
  (void *data, const char *name, const char **attr)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = START_TAG;
    if (strcmp(name, xml_strings[NUM_ENTRIES]) == 0) {
      user_data->element = NUM_ENTRIES;
    }
    else if (strcmp(name, xml_strings[S]) == 0) {
      user_data->element = S;
    }
    else {
      user_data->element = UNKNOWN_ELEMENT;
    };
  }


  void XMLCALL startTablePropertiesElement
  (void *data, const char *name, const char **attr)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    user_data->tag_type = START_TAG;
    if (strcmp(name, xml_strings[TABLE_PROPERTIES]) == 0) {
      user_data->element = TABLE_PROPERTIES;
    }
    else if (strcmp(name, xml_strings[POLY_VERTEX_LABEL_TYPE]) == 0) {
      user_data->element = POLY_VERTEX_LABEL_TYPE;
    }
    else if (strcmp(name, xml_strings[TRIANGULATION_TYPE]) == 0) {
      user_data->element = TRIANGULATION_TYPE;
    }
    else if (strcmp(name, xml_strings[SEPARATION_TYPE]) == 0) {
      user_data->element = SEPARATION_TYPE;
    }
    else if (strcmp(name, xml_strings[SEPARATE_OPPOSITE]) == 0) {
      user_data->element = SEPARATE_OPPOSITE;
    }
    else if (strcmp(name, xml_strings[ISO_POLY_ORIENTATION]) == 0) {
      user_data->element = ISO_POLY_ORIENTATION;
    }
    else {
      user_data->element = UNKNOWN_ELEMENT;
    }
  }


  void XMLCALL endTablePropertiesElement(void *data, const char *name)
  {
    USER_DATA * user_data = (USER_DATA *) data;

    if (user_data == NULL) return;

    if (strcmp(name, xml_strings[TABLE_PROPERTIES]) == 0) {
      user_data->element = TABLE_PROPERTIES;
    }
    else {
      user_data->element = UNKNOWN_ELEMENT;
    }
  }

}


// ***************************************************************
// write_xit
// ***************************************************************

// local functions
namespace {
  void write_header_info
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table,
   const char * version_str);
  void write_start_tag(std::ostream & out, const XML_ELEMENT xml_element);
  void write_end_tag(std::ostream & out, const XML_ELEMENT xml_element);
  void writeln_start_tag(std::ostream & out, const XML_ELEMENT xml_element);
  void writeln_end_tag(std::ostream & out, const XML_ELEMENT xml_element);
  template <typename T> void write_xml_element
  (std::ostream & out, const XML_ELEMENT xml_element, T x);
  template <typename T> void writeln_xml_element
  (std::ostream & out, const XML_ELEMENT xml_element, T x);  
  void write_table_type
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  void write_dimension
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  void write_polytope
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  void write_polytope_V1
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  void write_isosurface_vertices
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  void write_table_encoding
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  void write_table_properties
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
  void write_mcube_lookup_table_entries
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
}

// Write xit file.
void IJKXIO::write_xit
(std::ostream & out, const XIT_VERSION_TYPE version,
 const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
{
  IJK::PROCEDURE_ERROR error("write_xit");

  if (version == XIT_VERSION_1_0) {
    write_xit_V1(out, table);
  }
  else if ((version == XIT_VERSION_2_0) || (version == XIT_VERSION_2_x)) {
    write_xit_V2(out, table);
  }
  else {
    error.AddMessage("Unable to write unknown .xit version.");
    throw error;
  }
}


// Write xit file.
void IJKXIO::write_xit_V2
(std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
{
  using std::endl;

  // Write header information.
  write_header_info(out, table, "2.0");

  // Write table type.
  write_table_type(out, table);

  // Write dimension.
  write_dimension(out, table);

  // Write polytope.
  write_polytope(out, table);

  // Write isosurface vertices.
  write_isosurface_vertices(out, table);

  // Write mcube lookup table.
  writeln_start_tag(out, TABLE);
  write_table_encoding(out, table);
  write_table_properties(out, table);
  write_mcube_lookup_table_entries(out, table);
  writeln_end_tag(out, TABLE);

  writeln_end_tag(out, ISOTABLE);
}


// Write (old) version 1.0 xit file.
void IJKXIO::write_xit_V1
(std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
{
  using std::endl;

  // Write header information.
  write_header_info(out, table, "1.0");

  // Write dimension.
  write_dimension(out, table);

  // Write polytope.
  write_polytope_V1(out, table);

  // Write isosurface vertices.
  write_isosurface_vertices(out, table);

  // Write mcube lookup table.
  out << "<table>" << endl;
  write_table_encoding(out, table);
  write_mcube_lookup_table_entries(out, table);
  out << "</table>" << endl;

  out << "</isotable>" << endl;
}

namespace {

  void write_start_tag(std::ostream & out, const XML_ELEMENT xml_element)
  {
    out << "<" << xml_strings[xml_element] << ">";
  }

  // Write start tag followed by new line.
  void writeln_start_tag(std::ostream & out, const XML_ELEMENT xml_element)
  {
    write_start_tag(out, xml_element);
    out << std::endl;
  }

  void write_end_tag(std::ostream & out, const XML_ELEMENT xml_element)
  {
    out << "</" << xml_strings[xml_element] << ">";
  }


  // Write end tag followed by new line.
  void writeln_end_tag(std::ostream & out, const XML_ELEMENT xml_element)
  {
    write_end_tag(out, xml_element);
    out << std::endl;
  }


  template <typename T>
  void write_xml_element
  (std::ostream & out, const XML_ELEMENT xml_element, T x)
  {
    write_start_tag(out, xml_element);
    out << " " << x << " ";
    write_end_tag(out, xml_element);
  }


  // Specialization to pass string by reference.
  void write_xml_element
  (std::ostream & out, const XML_ELEMENT xml_element, 
   const std::string & str)
  {
    write_xml_element(out, xml_element, str.c_str());
  }


  template <typename T>
  void writeln_xml_element
  (std::ostream & out, const XML_ELEMENT xml_element, T x)
  {
    write_xml_element(out, xml_element, x);
    out << std:: endl;
  }


  void write_xml_header(std::ostream & out)
  {
    out << "<?xml version=\"1.0\"?>" << std::endl;
  }

  void write_isotable_header
  (std::ostream & out, const char * version)
  {
    writeln_start_tag(out, ISOTABLE);
    out << "<!-- Isosurface lookup table -->" << std::endl;
    writeln_xml_element(out, VERSION, version);
  }

  std::string get_date()
  {
    const size_t BUFFER_LENGTH = 16;
    char buffer[BUFFER_LENGTH];
    time_t rawtime;
    struct tm * timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, BUFFER_LENGTH, "%Y-%m-%d", timeinfo);
    std::string date_str = buffer;
    return date_str;
  }


  void write_creation_date(std::ostream & out)
  {
    std::string date_str = get_date();
    writeln_xml_element(out, CREATION_DATE, date_str);
  }


  void write_table_type
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    out << "<tableType> " << table.Properties().TableTypeString()
        << " </tableType>" << endl;
  }

  
  void write_dimension
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    out << "<dimension> " << table.Dimension() << "  "
        << table.SimplexDimension() << " </dimension>" << endl;
  }


  void write_header_info
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table,
   const char * version_str)
  {
    write_xml_header(out);
    write_isotable_header(out, version_str);
    write_creation_date(out);
  }


  void write_polytope_shape
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    writeln_xml_element(out, POLY_SHAPE, table.Polytope().ShapeString());
  }


  void write_polytope_vertices
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    writeln_start_tag(out, VERTICES);
    writeln_xml_element
      (out, NUM_VERTICES, table.Polytope().NumVertices());
    for (int i = 0; i < table.Polytope().NumVertices(); i++) {
      out << "<c> ";
      for (int d = 0; d < table.Dimension(); d++) {
        out << table.Polytope().VertexCoord(i, d) << " ";
      };
      out << "</c>" << endl;
    }
    writeln_end_tag(out, VERTICES);
  }

  void write_polytope_edges
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    writeln_start_tag(out, EDGES);
    writeln_xml_element(out, NUM_EDGES, table.Polytope().NumEdges());
    for (int i = 0; i < table.Polytope().NumEdges(); i++) {
      out << "<v> ";
      for (int j = 0; j < 2; j++) {
        out << table.Polytope().EdgeEndpoint(i, j) << " ";
      };
      out << "</v>" << endl;
    };
    writeln_end_tag(out, EDGES);
  }


  void write_polytope_facets
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    writeln_start_tag(out, FACETS);
    writeln_xml_element(out, NUM_FACETS, table.Polytope().NumFacets());
    for (int i = 0; i < table.Polytope().NumFacets(); i++) {
      out << "<f> ";
      out << table.Polytope().NumFacetVertices(i) << " ";
      for (int jv = 0; jv < table.Polytope().NumVertices(); jv++) {
        if (table.Polytope().IsVertexInFacet(i, jv)) {
          out << jv << " ";
        };
      };
      out << "</f>" << endl;
    };
    writeln_end_tag(out, FACETS);
  }


  void write_polytope
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    writeln_start_tag(out, POLY);
    write_polytope_shape(out, table);
    write_polytope_vertices(out, table);
    write_polytope_edges(out, table);
    write_polytope_facets(out, table);
    writeln_end_tag(out, POLY);
  }


  void write_polytope_V1
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    out << "<poly>" << endl;
    write_polytope_vertices(out, table);
    write_polytope_edges(out, table);
    write_polytope_facets(out, table);
    out << "</poly>" << endl;
  }


  void write_isosurface_vertices
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    writeln_start_tag(out, ISO_VERTICES);
    writeln_xml_element(out, NUM_VERTICES, table.NumIsosurfaceVertices());
    for (int i = 0; i < table.NumIsosurfaceVertices(); i++) {
      out << "<w> ";
      IJKMCUBE_TABLE::ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE vtype =
        table.IsosurfaceVertex(i).Type();
      int iv, ie, jf, numc;
      switch(vtype) {
      case IJKMCUBE_TABLE::ISOSURFACE_VERTEX::VERTEX:
        iv = table.IsosurfaceVertex(i).Face();
        write_xml_element(out, IN_V, iv);
        break;

      case IJKMCUBE_TABLE::ISOSURFACE_VERTEX::EDGE:
        ie = table.IsosurfaceVertex(i).Face();
        write_xml_element(out, IN_E, ie);
        break;

      case IJKMCUBE_TABLE::ISOSURFACE_VERTEX::FACET:
        jf = table.IsosurfaceVertex(i).Face();
        write_xml_element(out, IN_F, jf);
        break;

      case IJKMCUBE_TABLE::ISOSURFACE_VERTEX::POINT:
        numc = table.IsosurfaceVertex(i).NumCoord();
        out << "<c> ";
        if (numc > 0) {
          out << "<c> ";
          for (int ic = 0; ic < numc; ic++) {
            out << table.IsosurfaceVertex(i).Coord(ic) << " ";
          }
          out  << "</c>";
        }
        break;
      }

      if (table.IsosurfaceVertex(i).IsLabelSet()) {
        out << " ";
        write_xml_element(out, L, table.IsosurfaceVertex(i).Label());
      }
      out << " </w>" << endl;
    }
    writeln_end_tag(out, ISO_VERTICES);
  }

  void write_table_encoding
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    writeln_xml_element(out, ENCODING, table.EncodingName());
  }


  void write_table_properties
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    writeln_start_tag(out, TABLE_PROPERTIES);
    if (!table.Properties().IsGridVertexLabelTypeUndefined()) {
      writeln_xml_element(out, POLY_VERTEX_LABEL_TYPE, 
                          table.Properties().GridVertexLabelTypeString());
    }

    if (!table.Properties().IsTriangulationTypeUndefined()) {
      writeln_xml_element(out, TRIANGULATION_TYPE,
                          table.Properties().TriangulationTypeString());
    }

    if (!table.Properties().IsSeparationTypeUndefined()) {
      writeln_xml_element(out, SEPARATION_TYPE,
                          table.Properties().SeparationTypeString());
    }

    if (table.Polytope().Shape() == IJKMCUBE_TABLE::CUBE) {
      if (!table.Properties().IsSeparateOppositeUndefined()) {
        writeln_xml_element(out, SEPARATE_OPPOSITE,
                            table.Properties().SeparateOppositeString());
      }
    }

    if (!table.Properties().IsIsoPolyOrientationUndefined()) {
      writeln_xml_element(out, ISO_POLY_ORIENTATION,
                          table.Properties().IsoPolyOrientationString());
    }
    writeln_end_tag(out, TABLE_PROPERTIES);
  }


  void write_mcube_lookup_table_entries
  (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table)
  {
    using std::endl;

    out << "<numEntries> " 
        << table.NumTableEntries() 
        << " </numEntries>" << endl;
    for (int it = 0; it < table.NumTableEntries(); it++) {
      const int nums = table.NumSimplices(it);
      out << "<s> ";
      out << nums << " ";
      for (int js = 0; js < nums; js++) {
        for (int kv = 0; kv < table.SimplexDimension()+1; kv++) {
          out << int(table.SimplexVertex(it, js, kv)) << " ";
        }
      }
      out << "</s>" << endl;
    };
  }

}

