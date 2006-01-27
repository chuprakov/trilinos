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

#include "OptionsFromStreamPack_CommandLineOptionsFromStreamProcessor.hpp"

// Define this if you want to debug the parser
//#define PRINT_COMMAND_LINE_OPTIONS_FROM_STREAM_PROCESSOR_TRACE

namespace OptionsFromStreamPack {

CommandLineOptionsFromStreamProcessor::CommandLineOptionsFromStreamProcessor(
  const std::string &options_file_name
  ,const std::string &extra_options_str
  )
  :options_are_processed_(false)
  ,options_file_name_(options_file_name)
  ,extra_options_str_(extra_options_str)
{}

void CommandLineOptionsFromStreamProcessor::set_options(
  Teuchos::RefCountPtr<OptionsFromStream> const& options
  )
{
  options_ = options;
}

Teuchos::RefCountPtr<OptionsFromStream>
CommandLineOptionsFromStreamProcessor::get_options() const
{
  return options_;
}

void CommandLineOptionsFromStreamProcessor::set_options_file_name( 
  const std::string &options_file_name
  )
{
  options_file_name_ = options_file_name;
}

std::string CommandLineOptionsFromStreamProcessor::get_options_file_name() const
{
  return options_file_name_;
}

void CommandLineOptionsFromStreamProcessor::set_extra_options_str( 
  const std::string &extra_options_str
  )
{
  extra_options_str_ = extra_options_str;
}

std::string CommandLineOptionsFromStreamProcessor::get_extra_options_str() const
{
  return extra_options_str_;
}

void CommandLineOptionsFromStreamProcessor::setup_commandline_processor(
  Teuchos::CommandLineProcessor *clp
  )
{
  clp->setOption("ofs-options-file",&options_file_name_,"The name of the file containing input options for OptionsFromStream object.");
  clp->setOption("ofs-extra-options",&extra_options_str_,"Extra options in format \"OptGroup1{name1=val1,...,namen=valn}:OptGroup2{name1=val1,...,namen=valn}:...\"");
  // RAB: 2006/01/27: Note, this value contains no semi-columns since this
  // conflicts with Trilinos' runtests script.  Therefore, I have to replace
  // the ',' separators with ';' below!
  // Note: we can leave off the last ',' since it turns out that the
  // way the new OptionsFromStream::parse_options(...) is written that
  // the last semicolon in an options group is not necessary!
}

void CommandLineOptionsFromStreamProcessor::process_options()
{
  if(options_are_processed_) return;
  // Process the file options first
  if(options_file_name_.length()) {
    std::ifstream options_in(options_file_name_.c_str());
    if(options_in) {
      if(!options_.get())
        options_ = Teuchos::rcp(new OptionsFromStream());
      options_->read_options(options_in);
    }
  }
  // Process the extra commandline options
  const int len = extra_options_str_.length();
  if(len) {
    if(!options_.get())
      options_ = Teuchos::rcp(new OptionsFromStream());
    const char colon = ':';
    const std::string::size_type npos = std::string::npos;
    std::ostringstream ooptsstream;
    ooptsstream << "\nbegin_options\n\n";
    std::string::size_type start_i = 0, last_i = 0;
    while(true) {
      last_i = extra_options_str_.find(colon,start_i);
      std::string optgroup = extra_options_str_.substr(start_i,last_i);
      std::replace( optgroup.begin(), optgroup.end(), ',',';' ); // See above!
      ooptsstream << "options_group " << optgroup << "\n";
      if(last_i == npos) break;
      start_i = last_i + 1;
    }
    ooptsstream << "\nend_options\n";
    const std::string options_str = ooptsstream.str();
#ifdef PRINT_COMMAND_LINE_OPTIONS_FROM_STREAM_PROCESSOR_TRACE
    std::cout << "options_str:\n" << options_str;
#endif
    std::istringstream ioptsstream(options_str);
    options_->read_options(ioptsstream);
  }
  options_are_processed_ = true;
}

Teuchos::RefCountPtr<OptionsFromStream>
CommandLineOptionsFromStreamProcessor::process_and_get_options()
{
  process_options();
  return get_options();
}

}	// end namespace OptionsFromStreamPack
