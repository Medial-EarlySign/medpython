//
// Implementation file for creating Embedding features
//

#include <MedProcessTools/MedProcessTools/EmbeddingGenerator.h>
#include <Logger/Logger/Logger.h>
#include <MedMat/MedMat/MedMat.h>


#define LOCAL_SECTION LOG_FTRGNRTR
#define LOCAL_LEVEL	LOG_DEF_LEVEL

//-------------------------------------------------------------------------------------------------------------------------
int EmbeddingGenerator::init(map<string, string>& mapper)
{
	for (auto entry : mapper) {
		string field = entry.first;
		if (field == "f_scheme") f_scheme = entry.second;
		else if (field == "f_layers") f_layers = entry.second;
		else if (field == "name_prefix") name_prefix = entry.second;
		else if (field == "to_layer") to_layer = stoi(entry.second);
		else if (field == "tags") boost::split(tags, entry.second, boost::is_any_of(","));
		else if (field != "fg_type")
			MLOG("Unknown parameter \'%s\' for BasicFeatGenerator\n", field.c_str());
	}

	// read must have scheme and layers

	if (emc.read_from_file(f_scheme) < 0) {
		MTHROW_AND_ERR("EmbeddingGenerator ERROR: Can't open or parse scheme file %s . A valid scheme file is a must.\n", f_scheme.c_str());
		return -1;
	}

	if (embedder.init_from_text_file(f_layers) < 0) {
		MTHROW_AND_ERR("EmbeddingGenerator ERROR: Can't open or parse layers file %s . A valid layers file is a must.\n", f_layers.c_str());
		return -1;
	}

	if (to_layer < 0 || to_layer > embedder.layers.size()) {
		MTHROW_AND_ERR("EmbeddingGenerator ERROR: to_layer %d is not in allowed range of %d-%d.\n", to_layer, 0, (int)embedder.layers.size());
		return -1;
	}

	e_dim = embedder.layers[to_layer].out_dim;
	//MLOG("==> 1 ==> to_layer %d e_dim %d\n", to_layer, e_dim);

	names.clear();
	set_names();

	req_signals = emc.sigs_to_load;

	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------
void EmbeddingGenerator::set_names()
{
	names.clear();
	e_dim = embedder.layers[to_layer].out_dim;
	for (int j = 0; j < e_dim; j++) {
		string name = "FTR_" + int_to_string_digits(serial_id, 6) + "." + name_prefix + ".col_" + int_to_string_digits(j, 4);
		names.push_back(name);
	}
}

//-------------------------------------------------------------------------------------------------------------------------
int EmbeddingGenerator::_learn(MedPidRepository& rep, const MedSamples& samples, vector<RepProcessor *> processors)
{
	// currently nothing done ... BUT we may need to once Model option is inside.

	e_dim = embedder.layers[to_layer].out_dim;
	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------
void EmbeddingGenerator::set_signal_ids(MedSignals& sigs)
{
	emc.init_sids(sigs);
}

//-------------------------------------------------------------------------------------------------------------------------
void EmbeddingGenerator::init_tables(MedDictionarySections& dict)
{
	emc.init_tables(dict);
}

//-------------------------------------------------------------------------------------------------------------------------
int EmbeddingGenerator::_generate(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data)
{
	// currently calling the generate_by_rec version
	return generate_by_rec(rec, features, index, num, _p_data);
}

//-------------------------------------------------------------------------------------------------------------------------
int EmbeddingGenerator::generate_by_rec(PidDynamicRec& rec, MedFeatures& features, int index, int num, vector<float *> &_p_data)
{
	int use_shrink = 1; // atm : always do this. We can optionally export to params.


	// prep times vector
	vector<int> times;
	for (int i = 0; i < num; i++)
		times.push_back(features.samples[index + i].time);

	map<int, float> out_line;
	vector<pair<int, float>> line;
	vector<float> out;
	for (int ver = 0; ver < times.size(); ver++) {
		emc.get_pid_out_line(rec, ver, times[ver], use_shrink, out_line);

		MedSparseMat::convert_map_to_line(out_line, line);
		embedder.apply_sparse(line, out, to_layer);
		for (int j = 0; j < e_dim; j++) {
			float *pfeat_j = _p_data[j] + index;
			pfeat_j[ver] = out[j];
		}

	}
	
	return 0;

	//MLOG("=====> pid %d index %d num %d times (%d) : %d :: to_layer %d :: e_dim %d\n", rec.pid, index, num, times.size(), times[0], to_layer, e_dim);

	// first : generating a sparse matrix
	MedSparseMat smat;
	emc.add_pid_lines(rec, smat, times, use_shrink); // not the most efficient by far... but a start

	//MLOG("Got smat : %d lines\n", smat.lines.size());
	//for (auto e : smat.lines[0]) MLOG("%d,%f : ", e.first, e.second);
	//MLOG("\n");
	// second : generate embedding
	MedMat<float> emat;
	embedder.get_all_embeddings(smat, to_layer, emat);
	//MLOG("Got emat of size %d x %d names %d\n", emat.nrows, emat.ncols, names.size());


	// at the moment assuming the order in emat is the order we need. This is true when the times are sorted to begin with
	for (int j = 0; j < e_dim; j++) {
		float *pfeat_j = _p_data[j] + index;
		for (int i = 0; i < num; i++)
			pfeat_j[i] = emat(i, j);
	}

	return 0;
}