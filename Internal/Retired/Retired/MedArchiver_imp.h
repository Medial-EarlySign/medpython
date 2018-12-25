		
// Templated add
template <class T> size_t MedArchiver::add (T& newData) {

	size_t size = newData.get_size() ;
	data.resize(index + size) ;
	
	if (newData.serialize(&(data[index])) != size)
		return -1 ;
	else {
		index += size ;
		return size ;
	}
}


// map add
template <class T, class S> size_t MedArchiver::add (map <T,S> & newData) {
	size_t size = newData.size() ;
	size_t tot_size = size ;
	add((unsigned char *) &size, sizeof(size_t)) ;
	for (auto it=newData.begin(); it!=newData.end(); ++it) {
		T temp_t = it->first ;
		S temp_s = it->second ;
		tot_size += add(temp_t);
		tot_size += add(temp_s);

	}
	return tot_size ;
}



// vector add
template <class T> size_t MedArchiver::add (vector<T>& newData) {
	size_t size = newData.size() ;
	size_t tot_size = size ;
	add((unsigned char *) &size, sizeof(size_t)) ;
	for (size_t i=0; i<size; i++)
		tot_size += add(newData[i]) ;

	return tot_size ;
}

// Extract
template <class T> size_t MedArchiver::get (T& newData) {

	size_t size = newData.deserialize(&(data[index])) ;
	if (size < 0)
		return -1 ;
	else {
		index += size ;
		return size ;
	}
}

// vector extract
template <class T> size_t MedArchiver::get (vector<T>& newData) {
	size_t size ;
	get((unsigned char *)&size, sizeof(size_t)) ;
	newData.resize(size) ;
	size_t tot_size = sizeof(size_t);
	for (size_t i=0; i<size; i++) {
		tot_size += get(newData[i]);
	}
	index+=tot_size;
	return tot_size ;
}



// map extract
template <class T, class S> size_t  MedArchiver::get (map <T,S> & newData) {
	size_t size ;
	get((unsigned char *)&size, sizeof(size_t)) ;
	size_t tot_size = sizeof(size_t);
	for (size_t i=0; i<size; i++) {
		
		T temp_t;
		tot_size += get( (T) temp_t);

		S temp_s;
		tot_size += get( (S) temp_s);

		newData[temp_t] = temp_s;
	}
	index+=tot_size;
	return tot_size ;
}