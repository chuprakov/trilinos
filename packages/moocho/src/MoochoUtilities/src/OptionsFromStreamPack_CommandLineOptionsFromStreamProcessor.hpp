// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef OPTIONS_FORM_STEAM_PACK_COMMANDLINE_OPTIONS_FROM_STREAM_PROCESSOR_HPP
#define OPTIONS_FORM_STEAM_PACK_COMMANDLINE_OPTIONS_FROM_STREAM_PROCESSOR_HPP

#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace OptionsFromStreamPack {

/** \brief Reads from a file and/or parses from the commandline to initalize
 * an <tt>OptionsFromStream</tt> object.
 *
 * ToDo: Finish documentation!
 */
class CommandLineOptionsFromStreamProcessor {
public:

  /** \brief Construct with default values for the options file name and extra
   * options.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->get_options_file_name()==options_file_name</tt>
   * <li><tt>this->get_extra_options_str()==extra_options_str</tt>
   * </ul>
   */
  CommandLineOptionsFromStreamProcessor(
    const std::string &options_file_name = ""
    ,const std::string &extra_options_str = ""
    );

  /** \brief Set the <tt>OptionsFromStream</tt> object that will be used to
   * fill options in to.
   *
   * Postconditions:<ul>
   * <li><tt>this->get_options()=options.get()</tt>
   * </ul>
   */
  void set_options(
    Teuchos::RefCountPtr<OptionsFromStream> const& options
    );

  /** \brief Just return the <tt>OptionsFromStream</tt> object in its current
   * state.
   *
   * This function does not force the options to be processed.
   */
  Teuchos::RefCountPtr<OptionsFromStream> get_options() const;

  /** \brief Set the options file name manually (can be used for default value
   * for --ofs-options-file commandline option).
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->get_options_file_name()==options_file_name</tt>
   * </ul>
   */
  void set_options_file_name( 
    const std::string &options_file_name
    );

  /** \brief Get the name of the file that will be used to read options.
   */
  std::string get_options_file_name() const;

  /** \brief Set the extra commandline options string (can be used for default value
   * for --ofs-extra-options commandline option).
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->get_extra_options_str()==extra_options_str</tt>
   * </ul>
   */
  void set_extra_options_str( 
    const std::string &extra_options_str
    );

  /** \brief Get the current value of the extra options.
   */
  std::string get_extra_options_str() const;

  /** \brief Setup a comandline processor before it processes commandline
   * options or reads form a file.
   *
   * Sets up two options on the commandline processor:<ul>
   * <li><tt>--ofs-options-file</tt>: Sets the value returned from <tt>this->get_options_file_name()</tt>
   * <li><tt>--ofs-extra-options</tt>: Sets the value returned from <tt>this->get_extra_options_str()</tt>
   * </ul>
   */
  void setup_commandline_processor(
    Teuchos::CommandLineProcessor *clp
    );

  /** \brief Read the options file and/or process the commandline options.
   *
   * If <tt>this->get_options_file_name().length()>0</tt>, then this file is
   * read and processed first.  Then, if
   * <tt>this->get_extra_options_str().length()>0</tt>, this is followed by
   * the processing of these extra options which may overide options set in
   * the input options file.
   */
  void process_options();

  /** \brief Calls <tt>process_options()</tt> and returns
   * <tt>get_options()</tt>
   */
  Teuchos::RefCountPtr<OptionsFromStream> process_and_get_options();

private:

  bool                                        options_are_processed_;
  std::string                                 options_file_name_;
  std::string                                 extra_options_str_;
  Teuchos::RefCountPtr<OptionsFromStream>     options_;

};

}	// end namespace OptionsFromStreamPack

#endif	// OPTIONS_FORM_STEAM_PACK_COMMANDLINE_OPTIONS_FROM_STREAM_PROCESSOR_HPP
