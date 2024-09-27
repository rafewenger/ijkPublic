/*!
 *  @file ijkcommand_line.tpp
 *  @brief Templates for command line options.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2023 Rephael Wenger

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

#ifndef _IJKCOMMAND_LINE_
#define _IJKCOMMAND_LINE_

#include <cstddef>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkstring.tpp"


namespace IJK {

  // forward declarations
  inline void throw_error_on_missing_argument
  (const int iarg, const int argc, char **argv, IJK::ERROR & error);
  inline void throw_error_on_missing_argument
  (const int iarg, const int argc, char **argv, const int num_arg,
   IJK::ERROR & error);

  // ***************************************************************
  // GET ARGUMENT OPTIONS
  // ***************************************************************

  /// Get float argument argv[iarg+1].
  inline float get_arg_float(const int iarg, const int argc, char **argv,
                           IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    float x;
    if (!IJK::string2val(argv[iarg+1], x)) {
      error.AddMessage
        ("Usage error.  Error in argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-numeric character in string: ", argv[iarg+1], "");
      throw error;
    }

    return(x);
  }


  /// @brief Get float argument argv[iarg+1].
  /// - Version without parameter IJK::ERROR.
  inline float get_arg_float(const int iarg, const int argc, char **argv)
  {
    IJK::ERROR error;
    return(get_arg_float(iarg, argc, argv, error));
  }

  /// Get two float arguments argv[iarg+1] and argv[iarg+2].
  inline void get_arg2_float
  (const int iarg, const int argc, char **argv,
   float & x1, float & x2, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, 2, error);

    if (!IJK::string2val(argv[iarg+1], x1)) {
      error.AddMessage
        ("Usage error.  Error in first argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-numeric character in string: ", argv[iarg+1], "");
      throw error;
    }

    if (!IJK::string2val(argv[iarg+2], x2)) {
      error.AddMessage
        ("Usage error.  Error in second argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-numeric character in string: ", argv[iarg+2], "");
      throw error;
    }
  }


  /// @brief Get two float arguments argv[iarg+1] and argv[iarg+2].
  /// - Version without IJK::ERROR.
  inline void get_arg2_float
  (const int iarg, const int argc, char **argv,
   float & x1, float & x2)
  {
    IJK::ERROR error;
    get_arg2_float(iarg, argc, argv, x1, x2, error);
  }


  /// Get integer argument argv[iarg+1].
  inline int get_arg_int(const int iarg, const int argc, char **argv,
                         IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    int x;
    if (!IJK::string2val(argv[iarg+1], x)) {
      error.AddMessage
        ("Usage error.  Error in argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-integer character in string: ", argv[iarg+1], "");
      throw error;
    }

    return(x);
  }


  /// @brief Get integer argument argv[iarg+1].
  /// - Version without parameter IJK::ERROR.
  inline int get_arg_int(const int iarg, const int argc, char **argv)
  {
    IJK::ERROR error;
    return(get_arg_int(iarg, argc, argv, error));
  }


  /// Get two integer arguments argv[iarg+1] and argv[iarg+2].
  inline void get_arg2_int
  (const int iarg, const int argc, char **argv,
   int & x1, int & x2, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, 2, error);

    if (!IJK::string2val(argv[iarg+1], x1)) {
      error.AddMessage
        ("Usage error.  Error in first argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-integer character in string: ", argv[iarg+1], "");
      throw error;
    }

    if (!IJK::string2val(argv[iarg+2], x2)) {
      error.AddMessage
        ("Usage error.  Error in second argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-integer character in string: ", argv[iarg+2], "");
      throw error;
    }
  }


  /// @brief Get two integer arguments argv[iarg+1] and argv[iarg+2].
  /// - Version without IJK::ERROR.
  inline void get_arg2_int
  (const int iarg, const int argc, char **argv,
   int & x1, int & x2)
  {
    IJK::ERROR error;
    get_arg2_int(iarg, argc, argv, x1, x2, error);
  }


  /// Get integer argument argv[iarg+1] and float arguments argv[iarg+2] and argv[iarg+3].
  inline void get_arg3_int_float_float
  (const int iarg, const int argc, char **argv,
   int & i1, float & x1, float & x2, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, 3, error);

    if (!IJK::string2val(argv[iarg+1], i1)) {
      error.AddMessage
        ("Usage error.  Error in first argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-integer character in string: ", argv[iarg+1], "");
      throw error;
    }

    if (!IJK::string2val(argv[iarg+2], x1)) {
      error.AddMessage
        ("Usage error.  Error in second argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-numeric character in string: ", argv[iarg+2], "");
      throw error;
    }

    if (!IJK::string2val(argv[iarg+3], x2)) {
      error.AddMessage
        ("Usage error.  Error in third argument for option: ", argv[iarg], "");
      error.AddMessage
        ("Non-numeric character in string: ", argv[iarg+3], "");
      throw error;
    }
  }


  /// @brief Get integer argument argv[iarg+1] and float arguments argv[iarg+2] and argv[iarg+3].
  /// - Version without IJK::ERROR.
  inline void get_arg3_int_float_float
  (const int iarg, const int argc, char **argv,
   int & i1, float & x1, float & x2)
  {
    IJK::ERROR error;
    get_arg3_int_float_float(iarg, argc, argv, i1, x1, x2, error);
  }


  /// @brief Get boolean (true or false) argument.
  /// - Accepts true, false, True, False, t, f, T, F, yes, no, y, n.
  inline bool get_boolean_argument
  (const int iarg, const int argc, char **argv, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    std::string arg = argv[iarg+1];
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);

    if (arg == "true" || arg == "t") { return(true); }
    else if (arg == "false" || arg == "f") { return(false); }
    else if (arg == "yes" || arg == "y") { return(true); }
    else if (arg == "no" || arg == "n") { return(false); }

    // Illegal argument.
    error.AddMessage
      ("Usage error. Argument for option ", argv[iarg],
       " must be \"true\" or \"false\" or \"yes\" or \"no\".");
    throw error;
  }


  /// Get string argument argv[iarg+1] and convert to list of arguments.
  template <typename ETYPE>
  inline void get_arg_multiple_arguments
  (const int iarg, const int argc, char **argv,
   std::vector<ETYPE> & v, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    if (!IJK::string2vector(argv[iarg+1], v)) {
      error.AddMessage
        ("Usage error.  Error in argument for option: ", argv[iarg], ".");
      error.AddMessage("Illegal character in string: \"", argv[iarg+1], "\"");
      throw error;
    }
  }


  /// @brief Get string argument argv[iarg+1] and convert to list of arguments.
  /// - Version without IJK::ERROR.
  template <typename ETYPE>
  inline void get_arg_multiple_arguments
  (const int iarg, const int argc, char **argv,
   std::vector<ETYPE> & v)
  {
    IJK::ERROR error;
    get_arg_multiple_arguments(iarg, argc, argv, v, error);
  }


  /// @brief Get string argument argv[iarg+1] and convert to list of (int) arguments.
  /// - Specialization for type int.
  template <>
  inline void get_arg_multiple_arguments
  (const int iarg, const int argc, char **argv,
   std::vector<int> & v, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    if (!IJK::string2vector(argv[iarg+1], v)) {
      error.AddMessage
        ("Usage error.  Error in argument for option: ", argv[iarg], ".");
      error.AddMessage
        ("Non-integer character in string: \"", argv[iarg+1], "\"");
      throw error;
    }
  }


  /// @brief Get string argument argv[iarg+1] and convert to list of arguments.
  /// - Specialization for type float.
  template <>
  inline void get_arg_multiple_arguments
  (const int iarg, const int argc, char **argv,
   std::vector<float> & v, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    if (!IJK::string2vector(argv[iarg+1], v)) {
      error.AddMessage
        ("Usage error.  Error in argument for option: ", argv[iarg], ".");
      error.AddMessage
        ("Non-numeric character in string: \"", argv[iarg+1], "\"");
      throw error;
    }
  }


  /// @brief Get string argument argv[iarg+1] and convert to list of (float) arguments.
  /// - Specialization for type float.
  /// - Version without IJK::ERROR.
  template <>
  inline void get_arg_multiple_arguments
  (const int iarg, const int argc, char **argv,
   std::vector<float> & v)
  {
    IJK::ERROR error;
    get_arg_multiple_arguments(iarg, argc, argv, v, error);
  }


  /// @brief Get string argument argv[iarg+1] and convert to list of (double) arguments.
  /// - Specialization for type double.
  template <>
  inline void get_arg_multiple_arguments
  (const int iarg, const int argc, char **argv,
   std::vector<double> & v, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    if (!IJK::string2vector(argv[iarg+1], v)) {
      error.AddMessage
        ("Usage error.  Error in argument for option: ", argv[iarg], ".");
      error.AddMessage
        ("Non-numeric character in string: \"", argv[iarg+1], "\"");
      throw error;
    }
  }


  /// Get string argument argv[iarg+1], append to list of arguments.
  template <typename ETYPE>
  inline void get_arg_multiple_append
  (const int iarg, const int argc, char **argv,
   std::vector<ETYPE> & v, IJK::ERROR & error)
  {
    throw_error_on_missing_argument(iarg, argc, argv, error);

    if (!IJK::string2vector_append(argv[iarg+1], v)) {
      error.AddMessage
        ("Usage error.  Error in argument for option: ", argv[iarg], ".");
      error.AddMessage("Illegal character in string: \"", argv[iarg+1], "\"");
      throw error;
    }
  }


  // **************************************************
  // CHECK NUMBER OF ARGUMENTS
  // **************************************************

  /// Throw an error on a missing argument.
  inline void throw_error_on_missing_argument
  (const int iarg, const int argc, char **argv, IJK::ERROR & error)
  {
    if (iarg+1 >= argc) { 
      error.AddMessage
        ("Usage error. Missing argument for option ", argv[iarg], ".");
      throw error;
    }
  }

  /// @brief Throw an error on a missing argument.
  /// - Version for more than one argument.
  inline void throw_error_on_missing_argument
  (const int iarg, const int argc, char **argv, const int num_arg,
   IJK::ERROR & error)
  {
    if (iarg+num_arg >= argc) { 
      error.AddMessage
        ("Usage error. Missing arguments for option ", argv[iarg], ".");
      error.AddMessage
        ("  Option ", argv[iarg], " requires ", num_arg, " arguments.");
      throw error;
    }
  }


  // **************************************************
  // CLASS COMMAND_LINE_OPTIONS
  // **************************************************

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE=int>
  class COMMAND_LINE_OPTIONS {

  public:

    /// Permitted argument.
    class ARG_INFO {
    public:
      std::vector<std::string> str;
      std::vector<std::string> help_message;

      ARG_INFO(const char * arg_str, const char * arg_help_message);
    };

    /// Option information.
    class OPTION_INFO {

    public:
      OPTION_TYPE option;
      std::string option_name;
      OPTION_GROUP option_group;
      NTYPE num_arguments;
      std::vector<std::string> str;      

      /// Argument string for usage and help messages.
      std::string arg_str;

      /// Help message.
      std::vector<std::string> help_message;

      /// If true, argument must be from a list of predefined strings.
      bool is_arg_choice;

      /// List of permitted arguments.
      std::vector<ARG_INFO> permitted_arg_list;

    public:
      OPTION_INFO(const OPTION_TYPE optA, const char * optA_name,
                  const OPTION_GROUP optA_group,
                  const char * optA_str,
                  const NTYPE optA_num_arguments,
                  const char * optA_help_label, const char * optA_help_message);

      // get commands
      NTYPE NumArguments() const
      { return(num_arguments); }
      OPTION_GROUP Group() const
      { return(option_group); }
      bool IsArgChoice() const
      { return(is_arg_choice); }

      // set commands

      /// @brief Add argument to option choices.
      /// - Return argument choice number.
      std::size_t AddArgChoice
      (const char * arg_str, const char * arg_help_message);

    };


  public:
    typedef enum 
      { OPTION_ITEM, NEWLINE_ITEM, MESSAGE_ITEM, 
        BEGIN_OR_ITEM, END_OR_ITEM } USAGE_ITEM_TYPE;

    /// Usage message item.
    class USAGE_ITEM {
    protected:
      USAGE_ITEM_TYPE type;
      OPTION_TYPE option;
      std::string message;

      void Init();

    public:
      USAGE_ITEM();
      USAGE_ITEM(const USAGE_ITEM_TYPE type);
      USAGE_ITEM(const OPTION_TYPE option);

      OPTION_TYPE Option() const
      { return(option); }
      USAGE_ITEM_TYPE Type() const
      { return(type); }
      const std::string & Message() const
      { return(message); }
      NTYPE MessageLength() const
      { return(message.length()); }

      void SetOption(const OPTION_TYPE option);
      void SetType(const USAGE_ITEM_TYPE type)
        { this->type = type; }
      void SetNewline(const bool flag);
      void SetMessage(const std::string & message);
    };

    typedef std::vector<USAGE_ITEM> USAGE_ITEM_VECTOR;

  protected:
    std::vector<bool> is_option_defined;
    std::vector<NTYPE> option_list_index;

  public:
    typedef std::vector<std::string> STRING_VECTOR;

  public:
    std::vector<OPTION_INFO> list;

    static const NTYPE DEFAULT_LABEL_WIDTH = 15;
    static const NTYPE DEFAULT_USAGE_OPTION_INDENT = 5;
    static const NTYPE DEFAULT_HELP_INDENT = 2;
    static const NTYPE DEFAULT_HELP_TEXT_INDENT = 20;
    static const NTYPE DEFAULT_ARG_LABEL_WIDTH = 15;
    static const NTYPE DEFAULT_HELP_ARG_INDENT = 6;
    static const NTYPE DEFAULT_HELP_ARG_TEXT_INDENT = 24;

    NTYPE usage_option_indent;
    std::string usage_message;
    std::vector<USAGE_ITEM_VECTOR> usage_options_message;
    NTYPE help_label_width;
    NTYPE help_indent;
    NTYPE help_text_indent;
    std::vector<NTYPE> help_label_tab;
    NTYPE help_arg_label_width;
    NTYPE help_arg_indent;
    NTYPE help_arg_text_indent;
    std::vector<NTYPE> help_arg_label_tab;

  public:
    COMMAND_LINE_OPTIONS();
    COMMAND_LINE_OPTIONS(const NTYPE label_width);

    // get functions

    /// Return true if option "optA" is defined.
    bool IsOptionDefined(const OPTION_TYPE optA) const;

    /// Return true if option group "group" is defined.
    bool IsGroupDefined(const OPTION_GROUP group) const;

    /// Return location (index) of option in C++ vector list.
    std::size_t OptionListIndex(const OPTION_TYPE optA) const;

    /// Return integer associated with option optA.
    std::size_t OptionNum(const OPTION_TYPE optA) const
    { return(std::size_t(optA)); }

    /// Return integer associated with optA_group
    std::size_t OptionGroupNum(const OPTION_GROUP optA_group) const
    { return(std::size_t(optA_group)); }

    /// Return reference to option optA.
    const OPTION_INFO & Option(const OPTION_TYPE optA) const;

    /// Get option.  Return false if option not found.
    bool GetOption(const char * s, OPTION_TYPE & option);

    // set functions

    /// Add an option.
    void AddOption
    (const OPTION_TYPE optA, const char * optA_name,
     const OPTION_GROUP optA_group,
     const char * opt_str,
     const NTYPE num_arguments,
     const char * help_label, const char * help_message);

    /// Add an option with no arguments.
    void AddOptionNoArg
    (const OPTION_TYPE optA, const char * optA_name,
     const OPTION_GROUP optA_group,
     const char * opt_str,
     const char * help_message)
    { AddOption
        (optA, optA_name, optA_group, opt_str, 0, "", help_message); }

    /// @brief Add an option with no arguments. (2 help lines.)
    /// - Version that creates help message with two lines.
    void AddOptionNoArg
    (const OPTION_TYPE optA, const char * optA_name,
     const OPTION_GROUP optA_group,
     const char * opt_str,
     const char * help_message1,
     const char * help_message2);


    /// Add an option with 1 argument.
    void AddOption1Arg
    (const OPTION_TYPE optA, const char * optA_name,
     const OPTION_GROUP optA_group,
     const char * opt_str,
     const char * arg_str, const char * help_message)
    { AddOption(optA, optA_name, optA_group, opt_str, 1,
                arg_str, help_message); }

    /// @overload
    /// @brief Add an option with 1 argument. (2 help lines.)
    /// - Version that creates help message with two lines.
    void AddOption1Arg
    (const OPTION_TYPE optA, const char * optA_name,
     const OPTION_GROUP optA_group,
     const char * opt_str,
     const char * arg_str,
     const char * help_message1,
     const char * help_message2);

    /// Add an option with 2 arguments.
    void AddOption2Arg
    (const OPTION_TYPE optA, const char * optA_name,
     const OPTION_GROUP optA_group,
     const char * opt_str,
     const char * arg_str, const char * help_message)
    { AddOption(optA, optA_name, optA_group, opt_str, 2,
                arg_str, help_message); }

    /// Add a synonym for option A, i.e., another string which invokes
    ///   option A.
    void AddSynonym(const OPTION_TYPE optA,
                    const char * optA_str);

    /// @brief Add permitted argument to option.
    /// - Return argument choice number.
    std::size_t AddArgChoice
    (const OPTION_TYPE optA, 
     const char * arg_str, 
     const char * arg_help_message);

    /*!
       *  @brief Add permitted argument to option.
       *  - Return argument choice number.
       *  - Version that creates help message with two lines.
       */
     std::size_t AddArgChoice
     (const OPTION_TYPE optA,
      const char * arg_str,
      const char * arg_help_message1,
      const char * arg_help_message2);

    /*!
      *  @brief Add permitted argument to option.
      *  - Return argument choice number.
      *  - Version that creates help message with four lines.
      */
    std::size_t AddArgChoice
    (const OPTION_TYPE optA,
     const char * arg_str,
     const char * arg_help_message1,
     const char * arg_help_message2,
     const char * arg_help_message3,
     const char * arg_help_message4);


    /// Set usage message.
    void SetUsageMessage(const char * msg);
    void SetUsageMessage(const char * s0, const char * s1);
    void SetUsageMessage(const char * s0, const char * s1, const char * s2);

    /// Add usage option.
    void AddUsageOption 
    (const OPTION_TYPE option, const OPTION_GROUP opt_group);
    void AddUsageOption
    (const OPTION_GROUP opt_group, const USAGE_ITEM_TYPE type);
    void AddUsageOptionNewline(const OPTION_GROUP opt_group)
    { AddUsageOption(opt_group, NEWLINE_ITEM); }
    void AddUsageOptionBeginOr(const OPTION_GROUP opt_group)
    { AddUsageOption(opt_group, BEGIN_OR_ITEM); }
    void AddUsageOptionEndOr(const OPTION_GROUP opt_group)
    { AddUsageOption(opt_group, END_OR_ITEM); }
    void AddUsageOptionMessage
    (const OPTION_GROUP opt_group, const std::string & message);

    /// Remove last usage option item.
    void PopUsageOption(const OPTION_GROUP opt_group);

    void SetUsageOptionIndent(const NTYPE usage_option_indent);

    // Help message set functions.
    void MakeHelpLabel
    (const std::size_t index, std::string & help_label) const;

    /// @brief Add line to help message.
    void AddToHelpMessage(const OPTION_TYPE optA, const char * str);

    /// @brief Add line to help  message. (std::string version.)
    void AddToHelpMessage(const OPTION_TYPE optA, const std::string & str)
      { AddToHelpMessage(optA, str.c_str()); }

    /// @brief Add two lines to help message.
    void AddToHelpMessage
    (const OPTION_TYPE optA, const char * str0, const char * str1);

    /// @brief Add three lines to help message.
    void AddToHelpMessage
    (const OPTION_TYPE optA, const char * str0, const char * str1,
     const char * str2);

    /// @brief Add four lines to help message.
    void AddToHelpMessage
    (const OPTION_TYPE optA, const char * str0, const char * str1,
     const char * str2, const char * str3);

    /// @brief Add five lines to help message.
    void AddToHelpMessage
    (const OPTION_TYPE optA, const char * str0, const char * str1,
     const char * str2, const char * str3, const char * str4);

    /// @brief Add six lines to help message.
    void AddToHelpMessage
    (const OPTION_TYPE optA, const char * str0, const char * str1,
     const char * str2, const char * str3, const char * str4,
     const char * str5);

    /// @brief Add default value message to help message.
    template <typename _VALUE_TYPE>
    void AddToHelpMessageDefaultValue
    (const OPTION_TYPE optA, const std::string & variable_name,
     const _VALUE_TYPE value);

    /// @brief Add default value message to help message.
    /// - Version using (char *) for variable_name.
    template <typename _VALUE_TYPE>
    void AddToHelpMessageDefaultValue
    (const OPTION_TYPE optA, const char * variable_name,
     const _VALUE_TYPE value)
    {
      AddToHelpMessageDefaultValue
        (optA, std::string(variable_name), value);
    }

    /// @brief Add one line to help arg message.
    void AddToHelpArgMessage
    (const OPTION_TYPE optA, const std::size_t iarg, const char * str);

    /// @brief Add two lines to help arg message.
    void AddToHelpArgMessage
    (const OPTION_TYPE optA, const std::size_t iarg, 
     const char * str0, const char * st1);

    void SetLabelWidth(const NTYPE label_width);
    void SetHelpIndent(const NTYPE help_indent);
    void SetHelpIndent
    (const NTYPE help_indent, const NTYPE help_text_indent);
    void SetHelpTextIndent(const NTYPE indent);
    void AddLabelTab(const NTYPE tab);
    void SetArgLabelWidth(const NTYPE arg_label_width);
    void SetHelpArgIndent(const NTYPE arg_indent);
    void SetHelpArgIndent
    (const NTYPE arg_indent, const NTYPE arg_text_indent);
    void SetHelpArgTextIndent(const NTYPE indent);
    void AddArgLabelTab(const NTYPE tab);


    /// Print usage options message.
    template <typename OSTREAM_TYPE>
    void PrintUsageOptions
    (OSTREAM_TYPE & out, const OPTION_GROUP opt_group) const;

    /// Print help message help_message.
    template <typename OSTREAM_TYPE>
    void PrintHelpMessage
    (OSTREAM_TYPE & out, const std::string & help_label,
     const std::vector<std::string> & help_message,
     const NTYPE help_label_width,
     const NTYPE help_indent,
     const NTYPE help_text_indent,
     const std::vector<NTYPE> & help_label_tab) const;

    /// Print help message of option in list[index].
    template <typename OSTREAM_TYPE>
    void PrintHelpMessage
    (OSTREAM_TYPE & out, const std::size_t list_index) const;

    /// Print help message of permitted argument iarg for option list[index].
    template <typename OSTREAM_TYPE>
    void PrintArgHelpMessage(OSTREAM_TYPE & out, const std::size_t index,
                             const std::size_t iarg) const;
  };


  // ***************************************************************
  // CLASS COMMAND_LINE_OPTIONS MEMBER FUNCTIONS
  // ***************************************************************

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::COMMAND_LINE_OPTIONS()
  {
    SetUsageOptionIndent(DEFAULT_USAGE_OPTION_INDENT);
    SetLabelWidth(DEFAULT_LABEL_WIDTH);
    SetHelpIndent(DEFAULT_HELP_INDENT, DEFAULT_HELP_TEXT_INDENT);
    SetArgLabelWidth(DEFAULT_ARG_LABEL_WIDTH);
    SetHelpArgIndent(DEFAULT_HELP_ARG_INDENT, DEFAULT_HELP_ARG_TEXT_INDENT);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  COMMAND_LINE_OPTIONS(const NTYPE label_width)
  {
    SetLabelWidth(label_width);
    SetArgLabelWidth(label_width);

    SetUsageOptionIndent(DEFAULT_USAGE_OPTION_INDENT);
    SetHelpIndent(DEFAULT_HELP_INDENT);
    SetHelpArgIndent(DEFAULT_HELP_ARG_INDENT, DEFAULT_HELP_ARG_TEXT_INDENT);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  bool COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  IsOptionDefined(const OPTION_TYPE optA) const
  {
    const std::size_t ioptA = OptionNum(optA);
    if (ioptA >= is_option_defined.size()) { return(false); }
    else {
      return(is_option_defined[ioptA]);
    }
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  std::size_t COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  OptionListIndex(const OPTION_TYPE optA) const
  {
    const std::size_t ioptA = OptionNum(optA);

    if (!IsOptionDefined(optA)) { return(0); }
    else if (ioptA >= option_list_index.size()) { return(0); }
    else {
      return(option_list_index[ioptA]);
    }
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  const typename 
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::OPTION_INFO & 
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  Option(const OPTION_TYPE optA) const
  {
    const std::size_t ioptA = OptionNum(optA);
    IJK::PROCEDURE_ERROR error("COMMAND_LINE_OPTIONS::Option");

    if (!IsOptionDefined(optA)) {
      error.AddMessage("Programming error. Option ", ioptA, " not defined.");
      throw error;
    }

    const std::size_t index = OptionListIndex(optA);
    return(list[index]);
  }

  // Get option.  Return false if option not found.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  bool COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  GetOption(const char * s, OPTION_TYPE & option)
  {
    for (std::size_t i = 0; i < list.size(); i++) {
      for (std::size_t j = 0; j < list[i].str.size(); j++) {
        if (list[i].str[j] == s) {
          option = list[i].option;
          return(true);
        }
      }
    }

    return(false);
  }


  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddOption
  (const OPTION_TYPE optA, const char * optA_name,
   const OPTION_GROUP optA_group,
   const char * optA_str,
   const NTYPE num_arguments,
   const char * arg_str, const char * help_message)
  {
    IJK::PROCEDURE_ERROR error("COMMAND_LINE_OPTIONS::AddOption");
    const std::size_t ioptA = OptionNum(optA);

    if (IsOptionDefined(optA)) {
      std::size_t index = OptionListIndex(optA);
      error.AddMessage
        ("Programming error.  Option ", list[index].option_name,
         " is already defined.");
      throw error;
    }

    if (ioptA >= is_option_defined.size()) {
      is_option_defined.resize(ioptA+1, false);
      option_list_index.resize(ioptA+1, 0);
    }

    is_option_defined[ioptA] = true;
    option_list_index[ioptA] = list.size();

    list.push_back
      (OPTION_INFO(optA, optA_name, optA_group, optA_str, num_arguments,
                   arg_str, help_message));

    AddUsageOption(optA, optA_group);
  }


  // @brief Add an option with no arguments.
  // - Version that creates help message with two lines.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddOptionNoArg
  (const OPTION_TYPE optA, const char * optA_name,
   const OPTION_GROUP optA_group,
   const char * opt_str,
   const char * help_message1,
   const char * help_message2)
  {
    AddOptionNoArg
      (optA, optA_name, optA_group, opt_str, help_message1);
    AddToHelpMessage(optA, help_message2);
  }


  // @brief Add an option with 1 argument.
  // - Version that creates help message with two lines.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddOption1Arg
  (const OPTION_TYPE optA, const char * optA_name,
   const OPTION_GROUP optA_group,
   const char * opt_str,
   const char * arg_str,
   const char * help_message1,
   const char * help_message2)
  {
    AddOption1Arg
      (optA, optA_name, optA_group, opt_str, arg_str, help_message1);
    AddToHelpMessage(optA, help_message2);
  }


  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddSynonym(const OPTION_TYPE optA, const char * optA_str)
  {
    IJK::PROCEDURE_ERROR error("COMMAND_LINE_OPTIONS::AddSynonym");
    const std::size_t ioptA = OptionNum(optA);

    if (!IsOptionDefined(optA)) {
      error.AddMessage
        ("Programming error.  Option ", ioptA, " not defined.");
      error.AddMessage("  Unable to add option synonym ", optA_str, ".");
      throw error;
    }

    const std::size_t index = option_list_index[ioptA];
    list[index].str.push_back(std::string(optA_str));
  }


  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  std::size_t COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddArgChoice
  (const OPTION_TYPE optA, const char * arg_str,
   const char * arg_help_message)
  {
    IJK::PROCEDURE_ERROR
      error("COMMAND_LINE_OPTIONS::AddArgChoice");
    const std::size_t ioptA = OptionNum(optA);

    if (!IsOptionDefined(optA)) {
      error.AddMessage
        ("Programming error.  Option ", ioptA, " not defined.");
      error.AddMessage("  Unable to add option argument ", arg_str, ".");
      throw error;
    }

    const std::size_t index = option_list_index[ioptA];
    return(list[index].AddArgChoice(arg_str, arg_help_message));
  }


  // Add argument to optA.
  // - Create help message with two lines.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  std::size_t COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddArgChoice
  (const OPTION_TYPE optA, const char * arg_str,
   const char * arg_help_message1,
   const char * arg_help_message2)
  {
    IJK::PROCEDURE_ERROR
      error("COMMAND_LINE_OPTIONS::AddArgChoice");
    const std::size_t ioptA = OptionNum(optA);

    if (!IsOptionDefined(optA)) {
      error.AddMessage
        ("Programming error.  Option ", ioptA, " not defined.");
      error.AddMessage("  Unable to add option argument ", arg_str, ".");
      throw error;
    }

    const std::size_t index = option_list_index[ioptA];
    const std::size_t iarg =
      list[index].AddArgChoice(arg_str, arg_help_message1);
    list[index].permitted_arg_list[iarg].help_message.push_back
      (std::string(arg_help_message2));

    return iarg;
  }

  // Add argument to optA.
  // - Create help message with four lines.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  std::size_t COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddArgChoice
  (const OPTION_TYPE optA, const char * arg_str,
   const char * arg_help_message1,
   const char * arg_help_message2,
   const char * arg_help_message3,
   const char * arg_help_message4)
  {
    IJK::PROCEDURE_ERROR
      error("COMMAND_LINE_OPTIONS::AddArgChoice");
    const std::size_t ioptA = OptionNum(optA);

    if (!IsOptionDefined(optA)) {
      error.AddMessage
        ("Programming error.  Option ", ioptA, " not defined.");
      error.AddMessage("  Unable to add option argument ", arg_str, ".");
      throw error;
    }

    const std::size_t index = option_list_index[ioptA];
    const std::size_t iarg =
      list[index].AddArgChoice(arg_str, arg_help_message1);
    list[index].permitted_arg_list[iarg].help_message.push_back
      (std::string(arg_help_message2));
    list[index].permitted_arg_list[iarg].help_message.push_back
       (std::string(arg_help_message3));
    list[index].permitted_arg_list[iarg].help_message.push_back
       (std::string(arg_help_message4));

    return iarg;
  }


  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddUsageOption
  (const OPTION_TYPE option, const OPTION_GROUP opt_group)
  {
    const std::size_t igroup = OptionGroupNum(opt_group);

    if (igroup >= usage_options_message.size()) {
      usage_options_message.resize(igroup+1);
    }

    usage_options_message[igroup].push_back(USAGE_ITEM(option));
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddUsageOption
  (const OPTION_GROUP opt_group, const USAGE_ITEM_TYPE type)
  {
    const std::size_t igroup = OptionGroupNum(opt_group);

    if (igroup >= usage_options_message.size()) {
      usage_options_message.resize(igroup+1);
    }

    const NTYPE k = usage_options_message[igroup].size();
    usage_options_message[igroup].resize(k+1);
    usage_options_message[igroup][k].SetType(type);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddUsageOptionMessage
  (const OPTION_GROUP opt_group, const std::string & message)
  {
    const std::size_t igroup = OptionGroupNum(opt_group);
    USAGE_ITEM usage_item;

    if (igroup >= usage_options_message.size()) {
      usage_options_message.resize(igroup+1);
    }

    usage_item.SetMessage(message);
    usage_options_message[igroup].push_back(usage_item);
  }

  // Remove last usage option item.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  PopUsageOption(const OPTION_GROUP opt_group)
  {
    const std::size_t igroup = OptionGroupNum(opt_group);

    if (igroup >= usage_options_message.size()) {
      // Do nothing.
      return;
    }

    if (usage_options_message[igroup].size() == 0) {
      // Do nothing.
      return;
    }

    usage_options_message[igroup].pop_back();
  }


  // Add line to help message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpMessage(const OPTION_TYPE optA, const char * str)
  {
    IJK::PROCEDURE_ERROR error("COMMAND_LINE_OPTIONS::AddToHelpMessage");
    const std::size_t ioptA = OptionNum(optA);

    if (!IsOptionDefined(optA)) {
      error.AddMessage
        ("Programming error.  Option ", ioptA, " not defined.");
      error.AddMessage
        ("  Unable to add help message \"", str, "\"");
      throw error;
    }

    const std::size_t index = option_list_index[ioptA];
    list[index].help_message.push_back(std::string(str));
  }


  // Add two lines to help message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpMessage
  (const OPTION_TYPE optA, const char * str0, const char *str1)
  {
    AddToHelpMessage(optA, str0);
    AddToHelpMessage(optA, str1);
  }


  // Add three lines to help message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpMessage
  (const OPTION_TYPE optA, const char * str0, const char *str1,
   const char * str2)
  {
    AddToHelpMessage(optA, str0);
    AddToHelpMessage(optA, str1);
    AddToHelpMessage(optA, str2);
  }


  // Add four lines to help message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpMessage
  (const OPTION_TYPE optA, const char * str0, const char *str1,
   const char * str2, const char * str3)
  {
    AddToHelpMessage(optA, str0);
    AddToHelpMessage(optA, str1);
    AddToHelpMessage(optA, str2);
    AddToHelpMessage(optA, str3);
  }


  // Add five lines to help message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpMessage
  (const OPTION_TYPE optA, const char * str0, const char *str1,
   const char * str2, const char * str3, const char * str4)
  {
    AddToHelpMessage(optA, str0);
    AddToHelpMessage(optA, str1);
    AddToHelpMessage(optA, str2);
    AddToHelpMessage(optA, str3);
    AddToHelpMessage(optA, str4);
  }

  // Add six lines to help message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpMessage
  (const OPTION_TYPE optA, const char * str0, const char *str1,
   const char * str2, const char * str3, const char * str4,
   const char * str5)
  {
    AddToHelpMessage(optA, str0);
    AddToHelpMessage(optA, str1);
    AddToHelpMessage(optA, str2);
    AddToHelpMessage(optA, str3);
    AddToHelpMessage(optA, str4);
    AddToHelpMessage(optA, str5);
  }


  // Add default value message to help message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  template <typename _VALUE_TYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpMessageDefaultValue
  (const OPTION_TYPE optA, const std::string & variable_name,
   const _VALUE_TYPE value)
  {
    std::string value_str;
    IJK::val2string(value, value_str);
    std::string msg =
      "Default value of " + variable_name + ": " + value_str + ".";
    AddToHelpMessage(optA, msg);
  }


  // Add line to help arg message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpArgMessage
  (const OPTION_TYPE optA, const std::size_t iarg, const char * str)
  {
    IJK::PROCEDURE_ERROR error("COMMAND_LINE_OPTIONS::AddToHelpArgMessage");
    const std::size_t ioptA = OptionNum(optA);

    if (!IsOptionDefined(optA)) {
      error.AddMessage
        ("Programming error.  Option ", ioptA, " not defined.");
      error.AddMessage
        ("  Unable to add help message \"", str, "\"");
      throw error;
    }

    const std::size_t index = option_list_index[ioptA];
    const std::size_t num_arg = list[index].permitted_arg_list.size();
    if (iarg >= num_arg) {
      error.AddMessage
        ("Programming error.  Option ", list[index].option_name, 
         " only has ", num_arg, " permitted arguments.");
      error.AddMessage("  Argument ", iarg, " is not defined.");
      error.AddMessage
        ("  Unable to add help message \"", str, "\"");
      throw error;
    }

    list[index].permitted_arg_list[iarg].help_message.push_back
      (std::string(str));
  }


  // Add two lines to help arg message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddToHelpArgMessage
  (const OPTION_TYPE optA, const std::size_t iarg,
   const char * str0, const char * str1)
  {
    AddToHelpArgMessage(optA, iarg, str0);
    AddToHelpArgMessage(optA, iarg, str1);
  }


  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetLabelWidth(const NTYPE label_width)
  {
    help_label_width = label_width;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetHelpIndent(const NTYPE help_indent)
  {
    this->help_indent = help_indent;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetHelpIndent
  (const NTYPE help_indent, const NTYPE help_text_indent)
  {
    SetHelpIndent(help_indent);
    SetHelpTextIndent(help_text_indent);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetHelpTextIndent(const NTYPE help_text_indent)
  {
    this->help_text_indent = help_text_indent;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetArgLabelWidth(const NTYPE arg_label_width)
  {
    help_arg_label_width = arg_label_width;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetHelpArgIndent(const NTYPE help_arg_indent)
  {
    this->help_arg_indent = help_arg_indent;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetHelpArgIndent
  (const NTYPE arg_indent, const NTYPE arg_text_indent)
  {
    SetHelpArgIndent(arg_indent);
    SetHelpArgTextIndent(arg_text_indent);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetHelpArgTextIndent(const NTYPE indent)
  {
    this->help_arg_text_indent = indent;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddLabelTab(const NTYPE tab)
  {
    help_label_tab.push_back(tab);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  AddArgLabelTab(const NTYPE tab)
  {
    help_arg_label_tab.push_back(tab);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetUsageMessage(const char * msg)
  {
    this->usage_message = msg;
  }


  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetUsageMessage
  (const char * s0, const char * s1)
  {
    this->usage_message = std::string(s0) + std::string(s1);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetUsageMessage
  (const char * s0, const char * s1, const char * s2)
  {
    this->usage_message = 
      std::string(s0) + std::string(s1) + std::string(s2);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  SetUsageOptionIndent(const NTYPE usage_option_indent)
  {
    this->usage_option_indent = usage_option_indent;
  }

  // Print usage options message.
  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  template <typename OSTREAM_TYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  PrintUsageOptions
  (OSTREAM_TYPE & out, const OPTION_GROUP opt_group) const
  {
    const std::size_t igroup = OptionGroupNum(opt_group);

    if (igroup >= usage_options_message.size()) {
      // No usage options messages. Print nothing.
      return;
    }

    bool flag_at_beginning_of_line = true;
    bool flag_or = false;
    NTYPE first_item_in_clause = true;
    for (std::size_t i = 0; i < usage_options_message[igroup].size(); i++) {
      const USAGE_ITEM_TYPE item_type =
        usage_options_message[igroup][i].Type();

      if (item_type == OPTION_ITEM) {
        if (flag_at_beginning_of_line) 
          { out << std::setw(usage_option_indent) << " "; }
        else if (!flag_or)
          { out << " "; }

        const OPTION_TYPE optA = usage_options_message[igroup][i].Option();
        const std::size_t index = OptionListIndex(optA);

        if (!flag_or) { out << "["; }

        for (std::size_t j = 0; j < list[index].str.size(); j++) {
          if (!first_item_in_clause) { 
            if (flag_or && flag_at_beginning_of_line) {
              out << " ";
            }
            else {
              out << " | "; 
            }
          }
          out << list[index].str[j];
          if (list[index].num_arguments > 0) {
            if (flag_or || j+1 == list[index].str.size())
              { out << " " << list[index].arg_str; }
          }
          first_item_in_clause = false;
          flag_at_beginning_of_line = false;
        }

        if (!flag_or) { 
          out << "]"; 
          first_item_in_clause = true;
        }

        flag_at_beginning_of_line = false;
      }
      else if (item_type == MESSAGE_ITEM) {
        if (flag_at_beginning_of_line) {
          out << std::setw(usage_option_indent) << " ";
        }
        else
          { out << " "; }

        out << usage_options_message[igroup][i].Message();
        flag_at_beginning_of_line = false;
      }
      else if (item_type == NEWLINE_ITEM) {
        if (flag_or) { out << " |"; }
        out << std::endl;
        flag_at_beginning_of_line = true;
      }
      else if (item_type == BEGIN_OR_ITEM) {
        flag_or = true;
        first_item_in_clause = true;

        if (flag_at_beginning_of_line) 
          { out << std::setw(usage_option_indent) << " "; }
        else 
          { out << " "; }

        out << "[";
        flag_at_beginning_of_line = false;
      }
      else if (item_type == END_OR_ITEM) {
        out << "]";
        flag_or = false;
        first_item_in_clause = true;
        flag_at_beginning_of_line = false;
      }
    }

    if (!flag_at_beginning_of_line) {
      if (flag_or) { out << "]"; }
      out << std::endl;
    }

  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  MakeHelpLabel(const std::size_t index, std::string & help_label) const
  {
    help_label = "";
    for (std::size_t j = 0; j < list[index].str.size(); j++) {
      if (j > 0) { help_label += ", "; }
      help_label += list[index].str[j];
    }

    if (list[index].num_arguments > 0) {
      if (list[index].str.size() > 1) {
        help_label += "  ";
      }
      else {
        help_label += " ";
      }
      help_label += list[index].arg_str;
    }
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  template <typename OSTREAM_TYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  PrintHelpMessage(OSTREAM_TYPE & out, 
                   const std::string & help_label,
                   const std::vector<std::string> & help_message,
                   const NTYPE help_label_width,
                   const NTYPE help_indent,
                   const NTYPE help_text_indent,
                   const std::vector<NTYPE> & help_label_tab) const
  {
    using std::endl;

    if (help_message.size() < 1) {
      // No help message.
      return;
    }

    std::string::size_type label_length = help_label.length();
    size_t imsg = 0;

    if (help_indent > 0) { out << std::setw(help_indent) << " "; }

    if (label_length <= std::size_t(help_label_width)) {
      out << std::setw(help_label_width) << std::left 
          << help_label
          << " : " << help_message[imsg] << endl;
      imsg++;
    }
    else {
      for (std::size_t j = 0; j < help_label_tab.size(); j++) {
        if (label_length <= std::size_t(help_label_tab[j])) {
          out << std::setw(help_label_tab[j]) << std::left 
              << help_label
              << " : " << help_message[0] << endl;
          imsg++;
          break;
        }
      }
    }

    if (imsg == 0) { out << help_label << " :" << endl; }

    for (std::size_t j = imsg; j < help_message.size(); j++) {
      out << std::setw(help_text_indent) << " " 
          << help_message[j] << endl;
    }
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  template <typename OSTREAM_TYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  PrintHelpMessage(OSTREAM_TYPE & out, const std::size_t index) const
  {
    std::string help_label;

    MakeHelpLabel(index, help_label);

    PrintHelpMessage
      (out, help_label, list[index].help_message,
       help_label_width, help_indent, help_text_indent, help_label_tab);

    if (list[index].IsArgChoice()) {
      for (std::size_t iarg = 0; 
           iarg < list[index].permitted_arg_list.size(); iarg++) {
        PrintArgHelpMessage(out, index, iarg);
      }
    }

  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  template <typename OSTREAM_TYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  PrintArgHelpMessage(OSTREAM_TYPE & out, const std::size_t index,
                      const std::size_t iarg) const
  {
    std::string help_label;

    for (std::size_t j = 0; j < list[index].permitted_arg_list[iarg].str.size();
         j++) {
      if (j > 0) { help_label += ", "; }
      help_label += list[index].permitted_arg_list[iarg].str[j];
    }

    PrintHelpMessage
           (out, help_label,
       list[index].permitted_arg_list[iarg].help_message,
       help_arg_label_width, help_arg_indent, help_arg_text_indent, 
       help_arg_label_tab);
  }


  // ***************************************************************
  // CLASS COMMAND_LINE_OPTIONS::OPTION_INFO MEMBER FUNCTIONS
  // ***************************************************************

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  OPTION_INFO::OPTION_INFO
  (const OPTION_TYPE optA, const char * optA_name,
   const OPTION_GROUP optA_group,
   const char * optA_str,
   const NTYPE optA_num_arguments,
   const char * optA_arg_str, const char * optA_help_message)
  {
    option = optA;
    option_name = std::string(optA_name);
    option_group = optA_group;
    str.push_back(std::string(optA_str));
    num_arguments = optA_num_arguments;
    arg_str = std::string(optA_arg_str);
    help_message.push_back(std::string(optA_help_message));
    is_arg_choice = false;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  std::size_t COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  OPTION_INFO::AddArgChoice
  (const char * arg_str, const char * arg_help_message)
  {
    const std::size_t iarg = permitted_arg_list.size();
    is_arg_choice = true;
    permitted_arg_list.push_back(ARG_INFO(arg_str, arg_help_message));
    return(iarg);
  }


  // ***************************************************************
  // CLASS COMMAND_LINE_OPTIONS::ARG_INFO MEMBER FUNCTIONS
  // ***************************************************************

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  ARG_INFO::ARG_INFO(const char * arg_str,
                     const char * arg_help_message)
  {
    str.push_back(arg_str);
    help_message.push_back(arg_help_message);
  }


  // ***************************************************************
  // CLASS COMMAND_LINE_OPTIONS::USAGE_ITEM MEMBER FUNCTIONS
  // ***************************************************************

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  USAGE_ITEM::USAGE_ITEM()
  {
    Init();
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  USAGE_ITEM::USAGE_ITEM(const USAGE_ITEM_TYPE type)
  {
    Init();
    this->type = type;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  USAGE_ITEM::USAGE_ITEM(const OPTION_TYPE option)
  {
    Init();
    SetOption(option);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  USAGE_ITEM::Init()
  {
    // this->option = OPTION_TYPE(0);
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  USAGE_ITEM::SetOption(const OPTION_TYPE option)
  {
    this->option = option;
    type = OPTION_ITEM;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  USAGE_ITEM::SetNewline(const bool flag)
  {
    type = NEWLINE_ITEM;
  }

  template <typename OPTION_TYPE, typename OPTION_GROUP, typename NTYPE>
  void COMMAND_LINE_OPTIONS<OPTION_TYPE,OPTION_GROUP,NTYPE>::
  USAGE_ITEM::SetMessage(const std::string & message)
  {
    this->message = message;
    type = MESSAGE_ITEM;
  }

}

#endif
