#include "TSFIntHashtable.h"

using namespace TSF;

TSFIntHashtable::TSFIntHashtable(int capacity)
	: data_(), count_(0), capacity_(TSFHashUtils::nextPrime(capacity))
{
	data_.resize(capacity_);
}

bool TSFIntHashtable::containsKey(int key) const 
{
	const TSFArray<IntPair>& candidates = data_[key % capacity_];
	
	for (int i=0; i<candidates.length(); i++)
		{
			if (candidates[i].key_ == key) 
				{
					mostRecentKey_ = key;
					mostRecentValue_ = candidates[i].value_;
					return true;
				}
		}
	return false;
}

int TSFIntHashtable::get(int key) const 
{
	if (!containsKey(key)) 
		{
			TSFError::raise("key=" + TSF::toString(key) + " not found in "
											"hashtable " + toString());
		}
	
	return mostRecentValue_;
}

void TSFIntHashtable::put(int key, int value) 
{
	int index = key % capacity_;

	// check for duplicate key
	for (int i=0; i<data_[index].length(); i++)
		{
			if (data_[index][i].key_ == key)
				{
					data_[index][i].value_ = value;
					return;
				}
		}
	
	// no duplicate key, so increment element count by one.
	count_++;

	// check for need to resize.
	if (count_ > capacity_) 
		{
			capacity_ = TSFHashUtils::nextPrime(capacity_+1);
			rehash();
			// recaluate index
			index = key % capacity_;
		}

	data_[index].append(IntPair(key, value));
}


void TSFIntHashtable::rehash()
{
	TSFArray<TSFArray<IntPair> > tmp(capacity_);

	for (int i=0; i<data_.length(); i++)
		{
			for (int j=0; j<data_[i].length(); j++)
				{
					int newIndex = data_[i][j].key_ % capacity_;
					tmp[newIndex].append(data_[i][j]);
				}
		}

	data_ = tmp;
}


string TSFIntHashtable::toString() const 
{
	string str = "[";
	int count = 0;
	for (int i=0; i<data_.length(); i++)
		{
			for (int j=0; j<data_[i].length(); j++)
				{
					str = str + "(" + TSF::toString(data_[i][j].key_)
						+ ", " + TSF::toString(data_[i][j].value_) + ")";
					if (count < count_-1) str = str + ",";
					count++;
				}
		}
	str = str + "]";
	return str;
}



	
	

	
		
	
	


	






			
