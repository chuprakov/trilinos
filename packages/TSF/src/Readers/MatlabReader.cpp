#include "MatlabReader.h"

using namespace TSF;

/**
 *  MatlabReader methods
 */

/**
 * Constructor
 */
/*
MatlabReader::MatlabReader(istream& ios)
{
  //cout << "In MatlabReader: going to readFile \n";
  lines_ = StrUtils::readFile(ios, '#');
  doRead();
}
*/

/**
 * Constructor
 */
MatlabReader::MatlabReader(const string& filename)
  : TSFMatrixReaderBase(filename)
{
  // filebuf fb;
  // fb.open(filename.c_str(), ios::in);
  // ifstream is(&fb);
  ifstream is(filename.c_str());
  lines_ = StrUtils::readFile(is, '#');
  is.close();
  // fb.close();
  doRead();
}

void MatlabReader::summary() const
{
  TSFArray<string> names;
  TSFArray<TSFSmartPtr<MatVec> > mvs;
  arrays_.arrayify(names, mvs);

  for(int i=0; i<names.length(); i++) {
    string name = names[i];
    TSFSmartPtr<MatVec> mv = mvs[i];

    if(mv->isVector()) {
      cout << "Vector " << name << " is " << mv->numRows() << " x 1" <<  endl;
    }
    else if(mv->isMatrix()) {
      cout << "Matrix " << name << " is "
	   << mv->numRows() << " x " << mv->numCols() << endl;
    }
    else {
      cout << name << " is undefined." << endl;
    }
  }
}

TSFLinearOperator MatlabReader::readMatrix(const TSFVectorType& vectorType,
					   const string& name) const
{                     
  TSFVectorSpace rowSpace = vectorType.createSpace(getNumRows(name));
  TSFVectorSpace colSpace = vectorType.createSpace(getNumCols(name));
  TSFLinearOperator rtn = vectorType.createMatrix(rowSpace, colSpace);
    
  createMatrix(name, rtn);
  return rtn;
}

TSFVector MatlabReader::readVector(const TSFVectorType& vectorType,
				   const string& name) const
{
  int size = arrays_.get(name)->getVecHash().size();
  TSFVector rtn = vectorType.createSpace(size).createMember();
  createVector(name, rtn);
  
  return rtn;
}

void MatlabReader::doRead()
{
  string delim = "(),=;[] ";
  TSFSmartPtr<MatVec> matVec;
  for (int i = 0; i < lines_.length(); i++)
	{
      cout << lines_[i] << "\n";
	  TSFArray<string> token = stringTokenizer(lines_[i], delim);
      for (int k = 0; k < token.length(); k++)
        {
          cout << token[k] << "\n";
        } 
	  if (!arrays_.containsKey(token[0]))
		{
	  cout << "Adding " << token[0] << " to arrays_ \n";
		  string id = token[0];
		  matVec = new MatVec(id);
		  arrays_.put(id, matVec);
		}

	  if (token.length() == 4) 
		{
          //cout << "   adding entry \n";
		  int ii = StrUtils::atoi(token[1]) - 1;
                  // int ii = atoi(token[1].c_str()) - 1;
		  int jj = StrUtils::atoi(token[2]) - 1;
		  // int jj = atoi(token[2].c_str()) - 1;
		  double d = StrUtils::atof(token[3]);
		  // double d = atof(token[3].c_str());

          cout << "here \n";
		  matVec->addEntry(ii, jj, d);
		}
	  else // must be a one-dimensional array
		{
          //cout << "shouldn't be here";
		  int ii = StrUtils::atoi(token[1]) - 1;
		  // int ii = atoi(token[1].c_str()) - 1;
		  double d = StrUtils::atof(token[2]);
		  // double d = atof(token[2].c_str());
		  matVec->addEntry(ii, d);
		}
	} 
}

    /**
     *  This converts a matlab number with D-05 to one with E-05 etc.
     *  so that it is legal for Java
     */

/*     public String eForD(String d) */
/*     { */
/*         StringBuffer sb = new StringBuffer(d); */
/*         for (int i = 0; i < sb.length(); i++) */
/*             { */
/*                 if ( (new Character(sb.charAt(i))).equals(new Character('D')) ) */
/*                     { */
/*                         sb.setCharAt(i, 'E'); */
/*                         break; */
/*                     } */
/*             } */
/*         return sb.toString(); */
/*     } */

    /**
     *  Creates a double[] containing the contents of a vector
     */

 void MatlabReader::createVector(const string name, TSFVector& vector) const
{
  cout << "\n In create vector: name = " << name << "\n";
  cout << "vector space: " << vector.space() << "\n";
  TSFSmartPtr<MatVec> vec = arrays_.get(name);
  TSFHashtable<int, double> hash = vec->getVecHash();
  cout << "  From CreateVector: hash size = " << hash.size() << "\n";
  TSFArray<int> indices;
  TSFArray<double> values;
  hash.arrayify(indices, values);
  for (int i = 0; i < hash.size(); i++)
{
  cout << "indices[" << i << "] " << indices[i] << " =  " << values[i] << "\n";
}
//cout << indices.toString() << "\n";
//  cout << values.toString() << "\n";
  for (int i = 0; i < hash.size(); i++)
	{
	  //vector.setElement(indices[i], values[i]);
      vector[indices[i]] =  values[i];
	}
  cerr << "out of create vector: name = " << name << "\n";
}


/**
 *  Method to check that name is valid.  Quits if name is not valid.
 */
const TSFSmartPtr<MatVec> MatlabReader::getMV(string name) const
{
  if (!arrays_.containsKey(name))
	{
	  cout << "In MatlabReader: vector" << name  <<  " not present" ;
	  exit(2);
	}
  return arrays_.get(name);


}


void MatlabReader::createMatrix(const string name, TSFLinearOperator& mat) const
{
  cout << " domain = " << mat.domain().dim() << " range = " << 
    mat.range().dim();
  TSFSmartPtr<MatVec> mv = getMV(name);
  TSFMatrixView matView(mat);
  
  int numRows = mv->numRows();
  int numIndicesPerRow [numRows];
  for(int i=0; i<numRows; i++) numIndicesPerRow[i] = mv->getRowHash(i).size();
  cerr << "rows: " << mv->numRows() << "\nnumIndicesPerRow: " << numIndicesPerRow << "\n";
  matView.setBandwidth(mv->numRows(), numIndicesPerRow);

  matView.freezeStructure();
  for (int i = 0; i < mv->numRows(); i++)
	{
      cout << "Row " << i << "\n";
	  TSFArray<int> indices(0);
	  TSFArray<double> values(0);
	  TSFHashtable<int, double> hash = mv->getRowHash(i);
	  hash.arrayify(indices, values);
      matView.setRowStructure(i, indices.length(), &(indices[0]));
      for (int j = 0; j < indices.length(); j++)
        {
          cout << "   " << indices[j] << "  " << values[j] << "\n";
          matView.setElement(i, indices[j], values[j]);
          cout << "     added element \n";
        }
      //cout << "   about to addToRow \n";
      //matView.addToRow(i, indices.length(), &(indices[0]), &(values[0])); 
      //cout << "   back from addToRow \n";
	}
  matView.freezeValues();
}

TSFArray<string> MatlabReader::stringTokenizer(const string& str, 
										const string& delim)
{
  TSFArray<string> rtn(0);
  unsigned int start = 0;
        
  while(start < str.length())
    {
      start =  findNextNonDelimiter(str, start, delim);
      int stop = findNextDelimiter(str, start, delim);
      if (start-stop == 0) return rtn;
      string sub = StrUtils::subString(str, start, stop);
      //cout << "In stringTokenizer: new token = " << sub << "\n";
      rtn.append(sub);
      start =  findNextNonDelimiter(str, stop, delim);
    }
  return rtn;
}

int MatlabReader::findNextDelimiter(const string& str, int offset, 
                                 const string& delim)
{
  //cout << "In findNextDelimiter \n";
  for (unsigned int i=0; i<(str.length()-offset); i++)
    {
      for (unsigned int j = 0; j < delim.length(); j++)
        {
          if (str[i + offset] == delim[j])
            {
              //cout << "returning " << i + offset << "\n";
              return i + offset;
            }
        }
    }
  return str.length();
}

int  MatlabReader::findNextNonDelimiter(const string& str, int offset,
                                     const string& delim)
{
  //cout << "In findNextNonDelimiter \n";
  bool isDelim = false;
  for (unsigned int i=0; i<(str.length()-offset); i++)
    {
      isDelim = false;
      //cout << "str[i+offset] = " << str[i+offset] << "\n";
      for (unsigned int j = 0; j < delim.length(); j++)
        {
          //cout <<  "delim[" << j << "] = " << delim[j] << "\n";
          if (str[i + offset] == delim[j])
            {
              //cout << "    breaking out... \n";
              isDelim = true;
              break;
            }
        }
      if (!isDelim)
        {
          //cout << "returning " << i + offset << "\n";
          return i + offset;
        }

    }
  return str.length();
}

/*
 *  MatVec methods
 */



void MatVec::addEntry(const int& i, const int& j, const double& d)
{
  isMatrix_ = true;
  cout << "In addEntry\n";
  bool dup = false;
  if (i > numRows_ - 1) 
    {
      numRows_ = i + 1;
      //TSFHashtable<int, double> ht(0);
      //cout << "In addEntry2 rows_.length() = " << rows_.length() << "\n ";
      rows_.resize(numRows_);
  cout << "In addEntry2.4\n";
  //rows_[i] = ht;
  cout << "In addEntry2.5\n";
    }
  if (j > numCols_ - 1) 
    {
      numCols_ = j + 1;
      //TSFHashtable<int, double> ht(0);
      cols_.resize(numCols_);
  cout << "In addEntry2.6\n";
  //cols_[i] = ht;
  cout << "In addEntry2.7\n";
    }
  dup = putWDups(rows_[i], j, d);
  cout << "In addEntry3\n";
  if (!dup)
	{
	  dup = putWDups(cols_[j], i, d);
	}
}



void MatVec::addEntry(const int& i, const double& d)
{
  isVector_ = true;
  if(i >= numRows_) numRows_ = i + 1;
  putWDups(vec_, i, d);
  cout <<"vec_.size() = " << vec_.size();
}


bool MatVec::putWDups(TSFHashtable<int, double>& hash, const int& i, 
					  const double& d)
{
  cout << "In MatVec putWDups: hash.size() ="<< hash.size() << "\n";
  if (hash.size() == 0) 
    {
      cout << "putting in" << i << " and " << d << "\n";  
      hash.put(i, d);
    }
  //bool hasKey = hash.containsKey(i);
  cout << "In MatVec putWDups2: \n";
  if (hash.containsKey(i))
    //if (hasKey)
	{
	  double dOld = hash.get(i);
	  Dup dup(i, 0, dOld, d);
	  dups_.append(dup);
	  return true;
	}
  else
    {
      cout << "putting in" << i << " and " << d << "\n";  
      hash.put(i, d);
    }
  return false;
}
