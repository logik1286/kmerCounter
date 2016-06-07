#include <kchashdb.h>
#include <string>
#include "kmerError.h"
#include <vector>
#include <array>
#include <limits>
#include <cstdio>
#include <iostream>
#include <list>

#define ENABLE_CACHE 

//Base template for databasing mechanism to support large data files
//
//This class implements a caching technique to keep track
//of a certain number (specified in topCacheSize) of the most
//frequent kmers on the fly. This elimiates the need to
//search the database for top kmers after processing.
//There is some overhead in computing the cache. There isn't
//a major improvement for small files, but it's much more
//noticable for larger cases that fall out of Kyoto's RAM cache
//(controlled by db.tune_map)... around 25% improvement in total runtime.
template <class KeyType, typename CounterType> class BaseKeyCache{

    typedef std::pair<KeyType, CounterType> EntryType;
	size_t topCacheSize;
	public:

    //typedef datum KeyCache<std:;array<KeyT, NKeys>>::datum;

	BaseKeyCache(std::string _cacheName, bool _cleanUp=true, size_t _topCacheSize=50):cacheName(_cacheName),cleanUp(_cleanUp),topCacheSize(_topCacheSize){

		//in case we didn't clean up before
		remove(cacheName.c_str());

		db.tune_options(kyotocabinet::HashDB::TLINEAR);
		db.tune_buckets(1000*1000*50);
		db.tune_map((size_t)1024*1024*1024*4); //1GB


		if(!db.open(cacheName.c_str(), kyotocabinet::HashDB::OWRITER | kyotocabinet::HashDB::OCREATE))
			throw KmerError(1, "Unable to create database");


	}

	KmerError getTopKeys(std::vector<EntryType> & entries, size_t total, CounterType minValueFilter = 0){

		entries.clear();
#ifdef ENABLE_CACHE
		if(total <= topCache.size()){
			//pull from cache
			entries.reserve(total);
			auto it = topCache.begin();
			for(auto i = 0; i < total; i++){
				entries.push_back(*it);
				it++;
			};

			return KmerError();
		}
#endif
		std::list<EntryType> temp;

		kyotocabinet::DB::Cursor * cur = db.cursor();

		cur->jump();

		datum key, data;
		loadNextEntry(cur, key, data);

		KmerError err;

		while (key.dptr!=nullptr){


			err = updateTopCache(key, data, temp, total, false);
			if(err.isError()) return err;

			delete key.dptr; key.dptr = nullptr;
			delete data.dptr; data.dptr = nullptr;
			loadNextEntry(cur, key, data);

		}

		std::cout << "Top keys got:" << temp.size() << " keys" << std::endl;

		delete cur; delete key.dptr; delete data.dptr;

		auto sortFunction = [](const EntryType & first, const EntryType & second){ return first.second > second.second;};
		if(temp.size() != total){
			//need to sort, sort hasn't been called yet
			temp.sort(sortFunction);
		}

		//copy to rtn vector
		entries.reserve(temp.size());
		for(auto & i : temp){
			entries.push_back(i);
		}

		return KmerError();
	}


	KmerError incrementKey(EntryType key){

		auto dKey = createDatum(key);

		auto checkRtn = db.check(dKey.first.dptr,dKey.first.dsize);

		KmerError err;

		if(checkRtn==-1){

			if(!db.add(dKey.first.dptr, dKey.first.dsize, dKey.second.dptr, dKey.second.dsize))
				return KmerError(1, "Error when adding key");

#ifdef ENABLE_CACHE
			err = updateTopCache(dKey.first, dKey.second, topCache, topCacheSize);
#endif

		}else{

			datum current;

			//it does exist. updatea
			current.dptr = db.get(dKey.first.dptr, dKey.first.dsize, &current.dsize);

			//if(!db.get(dKey.first.dptr, dKey.first.dsize, current.dptr, current.dsize))
			if(current.dptr == nullptr)
				return KmerError(2, "Failure during fetch. Odd, GDBM said this entry existed.");

			CounterType * cPtr = (CounterType *)current.dptr;

			//add this counter, but don't allow overflows
			if(std::numeric_limits<CounterType>::max() - *cPtr > key.second)//no overflow
				*cPtr += key.second;
			else //overflow will happen
				*cPtr = std::numeric_limits<CounterType>::max();

			if(!db.set(dKey.first.dptr, dKey.first.dsize, current.dptr, current.dsize))
				err = KmerError(1, "Error when updating key ");

#ifdef ENABLE_CACHE
			err |= updateTopCache(dKey.first, current, topCache, topCacheSize);
#endif

			delete current.dptr;

		}


		return err;


	}


	protected:

	std::string cacheName;
	kyotocabinet::HashDB db;
	bool cleanUp;
	struct datum{
		char * dptr;
		size_t dsize;
	};

	virtual KeyType getKeyFromDatum(datum & data) = 0;
	virtual std::pair<datum, datum> createDatum(EntryType & input) = 0;

	public:
	void dumpAllKeys(){

		kyotocabinet::DB::Cursor * cur = db.cursor();
		cur->jump();


		datum key, data;


		loadNextEntry(cur, key, data);

		while (key.dptr!=nullptr){


			KeyType k = getKeyFromDatum(key);

			CounterType * c = (CounterType *)data.dptr;

			for(auto & i : k)
				std::cout << +i << ",";

			std::cout << " = " << +*c << std::endl;


			delete key.dptr; key.dptr = nullptr;
			delete data.dptr; data.dptr = nullptr;
			loadNextEntry(cur, key, data);

		};

		delete key.dptr;
		delete data.dptr;

		delete cur;

		return;


	}
	protected:


	std::list<EntryType> topCache;

	KmerError updateTopCache(datum & key, datum & data, std::list<EntryType> & temp, size_t total, bool verbose = false){

		auto sortFunction = [](const EntryType & first, const EntryType & second){ return first.second > second.second;};


		CounterType c = *((CounterType *)data.dptr);

		KeyType k = getKeyFromDatum(key);

		if(temp.size()==0){
			//Cache is empty, just add it
			if(verbose) std::cout << "Cache primed" << std::endl;
			temp.push_back(std::make_pair(k,c));
			return KmerError();
		};

		if(c > temp.back().second || temp.size() < total){

			if(verbose)
				std::cout << "Updating... key = " << k[0] << " count = " << +c << " current size = " << temp.size() << " total = " << total << std::endl;

			auto currentPos = temp.end();
			auto nextPos = temp.end();
			bool found = false;
			for(auto i = temp.begin(); i != temp.end(); i++){

				if(c > i->second && !found){
					nextPos = i;
					found = true;
				};


				if(i->first == k){
					//Found my current entry. Flag, move on to next
					currentPos = i;
					break;
				};

			};

			temp.insert(nextPos, std::make_pair(k,c));

			if(currentPos!=temp.end())
				temp.erase(currentPos);


			while(temp.size() > total) temp.pop_back();

		};


		return KmerError();

	};

	void loadNextEntry(kyotocabinet::DB::Cursor * cur, datum & key, datum & data){
		key.dptr = cur->get_key(&key.dsize, false);
		if(key.dptr!=nullptr)
			data.dptr = cur->get_value(&data.dsize, true);
	};


	virtual ~BaseKeyCache(){

		//close cache
		db.close();

		//now delete database
		if(cleanUp)
			remove(cacheName.c_str());

	};

};


//TODO: Remove cacheName as input parameter, have it autogenerated using a unique ID.
//Right now the application uses the same name for its temporary files, which
//could be bad when trying to run multiple instances.
template <class KeyType, typename CounterType> class KeyCache : public BaseKeyCache<KeyType, CounterType>{

	typedef typename BaseKeyCache<KeyType, CounterType>::datum datum;


	public:
	typedef std::pair<KeyType, CounterType> EntryType;
	KeyCache(std::string cacheName, bool cleanUp=false, size_t topCache=100);

	~KeyCache();

	protected:




	private:
	std::string cacheName;
	kyotocabinet::HashDB db;
	bool cleanUp;

	KeyType getKeyFromDatum(datum & data);
	std::pair<datum, datum> createDatum(EntryType & input);

};

//array type
template<> template <typename KeyT, size_t Nkeys, typename CounterT> class KeyCache<std::array<KeyT, Nkeys>, CounterT> : public BaseKeyCache<std::array<KeyT, Nkeys>, CounterT>{

	//ensure key using unsigned counters and key types
	static_assert(std::is_unsigned<KeyT>::value == true, "KeyCache array key precision should be unsigned.");
	static_assert(std::is_unsigned<CounterT>::value == true, "KeyCache counter precision should be unsigned.");


	typedef typename BaseKeyCache<std::array<KeyT, Nkeys>, CounterT>::datum datum;

	public:
	typedef KeyT KeyPrecision;
	typedef std::array<KeyPrecision, Nkeys> KeyType;
	typedef CounterT CounterType;
	typedef std::pair<KeyType, CounterType> EntryType;

	//typedef datum KeyCache<std:;array<KeyT, NKeys>>::datum;

	KeyCache(std::string _cacheName, bool _cleanUp=true, size_t _topCacheSize=100):BaseKeyCache<std::array<KeyT, Nkeys>, CounterT>(_cacheName,_cleanUp, _topCacheSize){};

	~KeyCache(){}

	protected:
	KeyType getKeyFromDatum(datum & data){
		KeyType k;
		KeyPrecision * kPtr = (KeyPrecision *)data.dptr;
		for(auto i = 0; i < k.size(); i++)
			k[i] = kPtr[i];

		return k;

	};

	std::pair<datum, datum> createDatum(EntryType & input){

		std::pair<datum,datum> rtn;
		rtn.first.dptr = (char *)input.first.data();
		rtn.first.dsize = sizeof(KeyPrecision) * input.first.size();

		rtn.second.dsize = sizeof(CounterType);
		rtn.second.dptr = (char *)&input.second;

		return rtn;

	}



};

template<> template <typename CounterT> class KeyCache<std::vector<char>, CounterT> : public BaseKeyCache<std::vector<char>, CounterT>{

	typedef typename BaseKeyCache<std::vector<char>, CounterT>::datum datum;



	public:

	//typedef datum KeyCache<std:;array<KeyT, NKeys>>::datum;
	typedef char KeyPrecision;
	typedef std::vector<KeyPrecision> KeyType;
	typedef CounterT CounterType;
	typedef std::pair<KeyType, CounterType> EntryType;


	KeyCache(std::string _cacheName, bool _cleanUp=true, size_t _topCacheSize=100):BaseKeyCache<std::vector<char>, CounterT>(_cacheName,_cleanUp, _topCacheSize){};

	~KeyCache(){}

	protected:
	KeyType getKeyFromDatum(datum & data){
		KeyType k;
		KeyPrecision * kPtr = (KeyPrecision *)data.dptr; k.resize(data.dsize);
		for(auto i = 0; i < k.size(); i++)
			k[i] = kPtr[i];

		return k;

	};

	std::pair<datum, datum> createDatum(EntryType & input){

		std::pair<datum,datum> rtn;
		rtn.first.dptr = (char *)input.first.data();
		rtn.first.dsize = sizeof(KeyPrecision) * input.first.size();

		rtn.second.dsize = sizeof(CounterType);
		rtn.second.dptr = (char *)&input.second;

		return rtn;

	}


};
