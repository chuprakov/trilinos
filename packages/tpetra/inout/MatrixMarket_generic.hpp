//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __MatrixMarket_generic_hpp
#define __MatrixMarket_generic_hpp

#include "Teuchos_ScalarTraits.hpp"
#ifdef HAVE_TEUCHOS_COMPLEX
#  include <complex>
#endif // HAVE_TEUCHOS_COMPLEX
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace Tpetra {
  namespace MatrixMarket {

    /// \brief Maximum line length for a Matrix Market file.
    ///
    /// Matrix Market files are only allowed to have this many ASCII
    /// characters per line.  In "tolerant" parsing mode, this
    /// restriction is relaxed.
    int maxLineLength();

    /// \brief True if line is a comment line, false otherwise.  
    ///
    /// If false is returned, line.substr(start,size) is a string
    /// containing valid, possibly noncomment data.  (size could be
    /// std::string::npos, indicating that there is no comment
    /// character in the line.  In that case, the substr() call with
    /// size==npos does the right thing.)  In tolerant mode, empty or
    /// whitespace-only lines are considered comment lines.
    bool 
    checkCommentLine (const std::string& line, 
		      size_t& start, 
		      size_t& size,
		      const size_t lineNumber,
		      const bool tolerant);

    /// \fn readPatternData
    /// \brief Read "<rowIndex> <colIndex>" from a line.
    ///
    /// Matrix Market files that store a sparse matrix pattern
    /// (a.k.a. an unweighted graph) do so with one pattern entry per
    /// line, stored as space-delimited ASCII text with the row index
    /// followed by the column index.  Both the row and column indices
    /// are 1-based.  This function attempts to read one line from the
    /// given input stream istr, extract the row and column indices,
    /// and write them to the corresponding output variables.
    ///
    /// \param istr [in/out] Input stream from which to attempt to
    ///   read one line.
    ///
    /// \param rowIndex [out] On output: if successful, the row index
    ///   read from the line.
    ///
    /// \param colIndex [out] On output: if successful, the column
    ///   index read from the line.
    ///
    /// \param lineNumber [in] The current line number.  Used only for
    ///   diagnostic error messages.
    ///
    /// \param tolerant [in] Whether to parse tolerantly.  In tolerant
    ///   mode, if this function fails in any way to read any of the
    ///   data, it will return false without throwing an exception.
    ///   Otherwise, this function will either throw an exception or
    ///   return true.
    ///
    /// \return True if this function successfully read the line from
    ///   istr and extracted both the row and column index, false
    ///   otherwise.  If tolerant==false, this function never returns
    ///   false; it either returns true or throws an exception.
    template<class Ordinal>
    bool
    readPatternData (std::istream& istr, 
		     Ordinal& rowIndex, 
		     Ordinal& colIndex, 
		     const size_t lineNumber,
		     const bool tolerant)
    {
      Ordinal the_rowIndex, the_colIndex;

      if (istr.eof() || istr.fail())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "Unable to read any data from line " << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      istr >> the_rowIndex;
      if (istr.fail())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "Failed to get row index from line " << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      else if (istr.eof())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "No more data after row index on line " << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      istr >> the_colIndex;
      if (istr.fail())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "Failed to get column index from line " << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      rowIndex = the_rowIndex;
      colIndex = the_colIndex;
      return true;
    }

    /// \fn readRealData
    /// \brief Read "<rowIndex> <colIndex> <realValue>" from a line.
    ///
    /// Matrix Market files that store a sparse matrix with real
    /// values do so with one sparse matrix entry per line.  It is
    /// stored as space-delimited ASCII text: the row index, the
    /// column index, and the real value, in that order.  Both the row
    /// and column indices are 1-based.  This function attempts to
    /// read one line from the given input stream istr, extract the
    /// row and column indices and the real value, and write them to
    /// the corresponding output variables.
    ///
    /// \param istr [in/out] Input stream from which to attempt to
    ///   read one line.
    ///
    /// \param rowIndex [out] On output: if successful, the row index
    ///   read from the line.
    ///
    /// \param colIndex [out] On output: if successful, the column
    ///   index read from the line.
    ///
    /// \param realValue [out] On output: if successful, the real
    ///   value read from the line.
    ///
    /// \param lineNumber [in] The current line number.  Used only for
    ///   diagnostic error messages.
    ///
    /// \param tolerant [in] Whether to parse tolerantly.  In tolerant
    ///   mode, if this function fails in any way to read any of the
    ///   data, it will return false without throwing an exception.
    ///   Otherwise, this function will either throw an exception or
    ///   return true.
    ///
    /// \return True if this function successfully read the line from
    ///   istr and extracted the row and column index and the nonzero
    ///   value, false otherwise.  If tolerant==false, this function
    ///   never returns false; it either returns true or throws an
    ///   exception.
    template<class Ordinal, class Real>
    bool
    readRealData (std::istream& istr, 
		  Ordinal& rowIndex, 
		  Ordinal& colIndex, 
		  Real& realValue, 
		  const size_t lineNumber,
		  const bool tolerant)
    {
      Real the_realValue;
      if (! readPatternData (istr, rowIndex, colIndex, lineNumber, tolerant))
	{
	  if (tolerant) 
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "Failed to read pattern data from line " << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      if (istr.eof())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "No more data after pattern data on line " << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      istr >> the_realValue;
      if (istr.fail())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "Failed to get real value from line " << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      realValue = the_realValue;
      return true;
    }

#ifdef HAVE_TEUCHOS_COMPLEX

    /// \fn readComplexData
    /// \brief Read "<rowIndex> <colIndex> <realPart> <imagPart>" from a line.
    ///
    /// Matrix Market files that store a sparse matrix with complex
    /// values do so with one sparse matrix entry per line.  It is
    /// stored as space-delimited ASCII text: the row index, the
    /// column index, the real part, and the imaginary part, in that
    /// order.  Both the row and column indices are 1-based.  This
    /// function attempts to read one line from the given input stream
    /// istr, extract the row and column indices and the real and
    /// imaginary parts, and write them to the corresponding output
    /// variables.
    ///
    /// \param istr [in/out] Input stream from which to attempt to
    ///   read one line.
    ///
    /// \param rowIndex [out] On output: if successful, the row index
    ///   read from the line.
    ///
    /// \param colIndex [out] On output: if successful, the column
    ///   index read from the line.
    ///
    /// \param realPart [out] On output: if successful, the real
    ///   part of the matrix entry's value read from the line.
    ///
    /// \param imagPart [out] On output: if successful, the imaginary
    ///   part of the matrix entry's value read from the line.
    ///
    /// \param lineNumber [in] The current line number.  Used only for
    ///   diagnostic error messages.
    ///
    /// \param tolerant [in] Whether to parse tolerantly.  In tolerant
    ///   mode, if this function fails in any way to read any of the
    ///   data, it will return false without throwing an exception.
    ///   Otherwise, this function will either throw an exception or
    ///   return true.
    ///
    /// \return True if this function successfully read the line from
    ///   istr and extracted all the output data, false otherwise.  If
    ///   tolerant==false, this function never returns false; it
    ///   either returns true or throws an exception.
    template<class Ordinal, class Real>
    bool
    readComplexData (std::istream& istr, 
		     Ordinal& rowIndex, 
		     Ordinal& colIndex, 
		     Real& realPart,
		     Real& imagPart,
		     const size_t lineNumber,
		     const bool tolerant)
    {
      Real the_realPart, the_imagPart;
      if (! readRealData (istr, rowIndex, colIndex, the_realPart, lineNumber, tolerant))
	{
	  if (tolerant) 
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "Failed to read pattern data and/or real value from line " 
		 << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      if (istr.eof())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "No more data after real value on line " 
		 << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      istr >> the_imagPart;
      if (istr.fail())
	{
	  if (tolerant)
	    return false;
	  else
	    {
	      std::ostringstream os;
	      os << "Failed to get imaginary value from line " 
		 << lineNumber << " of input";
	      throw std::invalid_argument(os.str());
	    }
	}
      realPart = the_realPart;
      imagPart = the_imagPart;
      return true;
    }
#endif // HAVE_TEUCHOS_COMPLEX

    template<class Ordinal>
    bool
    readPatternLine (const std::string& line, // only the data-containing part
		     Ordinal& rowIndex, 
		     Ordinal& colIndex, 
		     const size_t lineNumber,
		     const bool tolerant)
    {
      // The part of the line that contains data
      std::istringstream istr (line);
      return readPatternData (istr, rowIndex, colIndex, lineNumber, tolerant);
    }

    template<class Ordinal, class Real>
    bool
    readRealLine (const std::string& line, 
		  Ordinal& rowIndex, 
		  Ordinal& colIndex, 
		  Real& realValue, 
		  const size_t lineNumber,
		  const bool tolerant)
    {
      size_t start, size;
      if (checkCommentLine (line, start, size, lineNumber, tolerant)) 
	{
	  // {
	  //   using std::cerr;
	  //   using std::endl;
	  //   cerr << "Comment line: " << lineNumber << endl;
	  // }
	  return false; // It's a comment line
	}
      // If it's an empty line, checkCommentLine() will throw an
      // exception if non-tolerant parsing is being performed, so
      // we need only return false otherwise.
      if (size == 0) 
	{
	  if (! tolerant)
	    throw std::logic_error("Should never get here! checkCommentLine() "
				   "is supposed to catch empty lines.");
	  else
	    {
	      // {
	      //   using std::cerr;
	      //   using std::endl;
	      //   cerr << "Empty line: " << lineNumber << endl;
	      // }
	      return false;
	    }
	}
      // The part of the line that contains data
      std::istringstream istr (line.substr (start, size));
      return readRealData (istr, rowIndex, colIndex, realValue, lineNumber, tolerant);
    }

#ifdef HAVE_TEUCHOS_COMPLEX
    template<class Ordinal, class Real>
    bool
    readComplexLine (const std::string& line, 
		     Ordinal& rowIndex, 
		     Ordinal& colIndex,
		     Real& realPart,
		     Real& imagPart,
		     const size_t lineNumber,
		     const bool tolerant)
    {
      size_t start, end;
      if (checkCommentLine (line, start, end, lineNumber, tolerant)) 
	return false; // It's a comment line
      // If it's an empty line, checkCommentLine() will throw an
      // exception if non-tolerant parsing is being performed, so
      // we need only return false otherwise.
      if (end == 0) 
	{
	  if (tolerant)
	    throw std::logic_error("Should never get here! checkCommentLine() "
				   "is supposed to catch empty lines.");
	  else
	    return false;
	}
      // The part of the line that contains data
      std::istringstream istr (line.substr (start, end));
      // Read the data
      Real the_realPart, the_imagPart;
      const bool success = 
	readComplexData (istr, rowIndex, colIndex, the_realPart, the_imagPart, 
			 lineNumber, tolerant);
      if (success)
	{
	  realPart = the_realPart;
	  imagPart = the_imagPart;
	}
      return success;
    }
#endif // HAVE_TEUCHOS_COMPLEX

    template<class Ordinal, class PatternCallback>
    std::pair<bool, std::vector<size_t> >
    readPatternCoordinateData (std::istream& in, 
			       PatternCallback add,
			       const size_t startingLineNumber, 
			       const bool tolerant)
    {
      std::string line;
      size_t lineNumber = startingLineNumber;
      bool anySucceeded = false;
      bool allSucceeded = true;
      std::vector<size_t> badLineNumbers; 
      size_t validDataLines = 0;
      while (getline (in, line))
	{
	  size_t start, size;
	  if (checkCommentLine (line, start, size, lineNumber, tolerant))
	    continue; // it's a comment line
	  const std::string theLine = line.substr (start, size);

	  Ordinal rowIndex, colIndex;
	  const bool localSuccess = readPatternLine (theLine, rowIndex, colIndex, lineNumber++, tolerant);
	  anySucceeded = anySucceeded || localSuccess;
	  allSucceeded = allSucceeded && localSuccess;
	  if (! localSuccess)
	    badLineNumbers.push_back (lineNumber);
	  else
	    {
	      add (rowIndex, colIndex);
	      validDataLines++;
	    }
	}
      if (lineNumber == startingLineNumber)
	anySucceeded = true; // Trivially true
    
      return std::make_pair (allSucceeded, badLineNumbers);
    }


    template<class Ordinal>
    bool
    readCoordinateDimensions (std::istream& in, 
			      Ordinal& numRows, 
			      Ordinal& numCols,
			      Ordinal& numNonzeros,
			      const size_t lineNumber,
			      const bool tolerant)
    {
      Ordinal the_numRows, the_numCols, the_numNonzeros;
      std::string line;

      if (! getline(in, line))
	{
	  std::ostringstream os;
	  os << "Failed to get data from line " << lineNumber << " of input stream";
	  throw std::invalid_argument (os.str());
	}
      std::istringstream istr (line);

      if (istr.eof() || istr.fail())
	{
	  std::ostringstream os;
	  os << "Unable to read any data from line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      istr >> the_numRows;
      if (istr.fail())
	{
	  std::ostringstream os;
	  os << "Failed to get number of rows from line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      else if (istr.eof())
	{
	  std::ostringstream os;
	  os << "No more data after number of rows on line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      istr >> the_numCols;
      if (istr.fail())
	{
	  std::ostringstream os;
	  os << "Failed to get number of columns from line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      else if (istr.eof())
	{
	  std::ostringstream os;
	  os << "No more data after number of columns on line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      istr >> the_numNonzeros;
      if (istr.fail())
	{
	  std::ostringstream os;
	  os << "Failed to get number of nonzeros from line " << lineNumber << " of input";
	  throw std::invalid_argument(os.str());
	}
      numRows = the_numRows;
      numCols = the_numCols;
      numNonzeros = the_numNonzeros;
      return true;
    }

  } // namespace MatrixMarket
} // namespace Tpetra

#endif // __MatrixMarket_generic_hpp
