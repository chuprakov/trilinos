#ifndef MATLABREADER_H
#define MATLABREADER_H

#include "TSF.h"
#include "StrUtils.h"
#include <string>
#include <iostream>
#include <fstream>

namespace TSF
{
  using std::string;
  using std::istream;

  /** 
   *  Dup keeps the data for duplicate entries. Specifically, it keeps
   *  the i and j location and the two values.  If an entry is repeated
   *  more than two times, there will be multiple entries.  
   */
  class Dup
    {

    public:
      /**
       *  Constructor for Matrix
       */
      Dup(const int& i, const int& j, const double& d, const double& dNew)
        :
        i_(i),
        j_(j),
        d_(d),
        dNew_(dNew)
        {
          ;
        }

      Dup(){;}

      /**
       *  Constructor for Vector
       */
      Dup(const int& i, const double& d, const double& dNew)
        :
        i_(i),
        d_(d),
        dNew_(dNew)
        {
          j_ = 0;
        }

      /**
       *  Methods to return the data
       */
      int getI() const  {return i_;}
      int getJ() const {return j_;}
      double getD() const {return d_;}
      double getDNew() const {return dNew_;}


    private:
      int i_;
      int j_;
      double d_;
      double dNew_;



    };





  /**
   *  Class for the data associated with a matrix or vec
   */

  class MatVec
    {

    public:

      /**
       * Empty Ctor.
       */
      MatVec() : isVector_(false), isMatrix_(false) {;}

      /**
       * Ctor.
       */
      MatVec(const string& name)
        :
        name_(name), isVector_(false), isMatrix_(false)
        {
          numRows_ = 0;
          numCols_ = 0;
        }


      /**
       *  add Entry for matrix
       */
      void addEntry(const int& i, const int& j, const double& d);
    
      /**
       *  add Entry for vector
       */
      void addEntry(const int& i, const double& d);

      /**
       *  Return the number of rows
       */
      int numRows() const {return numRows_;}


      /**
       *  Return the number of cols
       */
      int numCols() const {return numCols_;}

      /**
       *  Is this MatVec a vector?
       */
      bool isVector() const {return isVector_;}
      
      /**
       *  Is this MatVec a matrix?
       */
      bool isMatrix() const {return isMatrix_;}

      /**
       *  Return the TSFHashtable associated with row i
       */
      const TSFHashtable<int, double>& getRowHash(int i) const {return rows_[i];}

      /**
       *  Return the TSFHashtable associated with col j
       */
      const TSFHashtable<int, double>& getColHash(int j) const {return cols_[j];}


      /**
       *  Return the TSFHashtable associated with a vector
       */
      const TSFHashtable<int, double>& getVecHash() const {return vec_;}


      /**
       *  Return the TSFArray containing the duplicates information
       */
      const TSFArray<Dup>& getDups() const  {return dups_;}

      /**
       *  Add new information if an entry for matrix location is
       *  not already present.  Addd to the list of duplicates if
       *  necessary.  
       *  NOTE: The first value is NOT replaced.  This keeps
       *  the row and column information consistent 
       */
      bool putWDups(TSFHashtable<int, double>& hash, const int& i, 
                    const double& d);

    private:
      TSFArray<TSFHashtable<int, double> > rows_;
      TSFArray<TSFHashtable<int, double> > cols_;

      TSFHashtable<int, double> vec_;
      string name_;

      TSFArray<Dup> dups_;

      int numRows_;
      int numCols_;

      bool isVector_;
      bool isMatrix_;

    };





  /**
   *  This reads a matlab .m file and makes a sparse problem
   *  The expected format is:
   *    a set of matlab statements of the form
   *    arrayName[i,j] = ...; or
   *    vectorName[i] = ...;
   */
  class MatlabReader : public TSFMatrixReaderBase
    {
    public:

      /**
       *  Constructor reads the file
       */
      /*
      MatlabReader(istream& ios);
      */

      /**
       *  Constructor reads the file
       */
      MatlabReader(const string& filename);

      /**
       * List names and dimensions of matrices and vector in the file
       */
      void summary() const;

      /**
       * Obtain a specific matrix from the file
       */
      TSFLinearOperator readMatrix(const TSFVectorType& vectortype,
					   const string& name) const;

      /**
       * Obtain a specific vector from the file
       */
      TSFVector readVector(const TSFVectorType& vectorType,
			   const string& name) const;

      /**
       * Obtain matrix named "A" from the file
       */
      TSFLinearOperator read(const TSFVectorType& vectorType) const
	{ return readMatrix(vectorType, "A"); }

      /** 
       * Main method to process the file.
       */
      void doRead();

      /**
       *  Fills a TSFVector containing the contents of a vector
       */

      void createVector(const string name, TSFVector& vector) const;

      /**
       *  Creates a sparse representation of both the row and col
       *   oriented information
       */
      void createMatrix(const string name, TSFLinearOperator& mat) const;

      /**
       *  Creates a double[] containing the contents of a vector
       */

      /**
       *  Returns the number of rows in the array name
       */
      int getNumRows(const string& name) const
	{
	  TSFSmartPtr<MatVec> mv = arrays_.get(name);
	  return mv->numRows();
	}

      /**
       *  Returns the number of cols in the array name
       */

      int getNumCols(const string& name) const
        {
          TSFSmartPtr<MatVec> mv = getMV(name);
          return mv->numCols();
        }

      /**
       *  Reurns the TSFArray of Dups containing the duplicates information
       */
      const TSFArray<Dup>& getDups(const string name) const
        {
          return getMV(name)->getDups();
        }


    private:


      TSFHashtable<string, TSFSmartPtr<MatVec> >  arrays_;
      TSFArray<string> lines_;


      /**
       *  Method to check that name is valid.  Quits if name is not valid.
       */
      const TSFSmartPtr<MatVec> getMV(string name) const;


      /**
       * String tokenizer.
       */
      TSFArray<string> stringTokenizer(const string& str, 
                                    const string& delim);

      int findNextDelimiter(const string& str, int offset, 
                            const string& delim);

      int findNextNonDelimiter(const string& str, int offset, 
                               const string& delim);

    };

}

#endif
