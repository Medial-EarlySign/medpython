#include "micNet2.h"

//=================================================================================
// NodeInfo
//=================================================================================

// format to read a node-info: 
// <id:[int]|type:[LR,SoftMax,Norm,Regression,MaxOut,Chooser,Noise]|n_hidden:[int]|sources:[int:int:...]|lr:[float]|drop:[float]|...
// ...keep|A:[float]|B:[float]|lambda:[float]|momentum:[float]|decay:[float]|max_width:[int]|noise_std:[float]|...
// ...wgt_std:[float]>
//.........................................................................................................................................
void NodeInfo::init_from_string(const string &in_str)
{
	string str_in = in_str;
	boost::replace_all(str_in, ">", "");

	vector<string> fields;
	boost::split(fields, str_in, boost::is_any_of("|"));

	for (auto fitem : fields) {
		vector<string> f;
		boost::split(f, fitem, boost::is_any_of(":"));

		if (f[0] == "id") id = stoi(f[1]);
		if (f[0] == "type") type = f[1];
		if (f[0] == "n_hidden") n_hidden = stoi(f[1]);
		if (f[0] == "sources") {
			for (int i=1; i<f.size(); i++)
				sources.push_back(stoi(f[1]));
		}
		if (f[0] == "lr") learn_rate = stof(f[1]);
		if (f[0] == "drop_in") drop_in = stof(f[1]);
		if (f[0] == "keep") keep_node = 1;
		if (f[0] == "A") A = stof(f[1]);
		if (f[0] == "B") B = stof(f[1]);
		if (f[0] == "lambda") lambda = stof(f[1]);
		if (f[0] == "momentum") momentum = stof(f[1]);
		if (f[0] == "decay") rate_decay = stof(f[1]);
		if (f[0] == "max_width") max_width = stoi(f[1]);
		if (f[0] == "noise_std") noise_std = stof(f[1]);
		if (f[0] == "wgt_std") wgt_std = stof(f[1]);
		if (f[0] == "train") to_train = stoi(f[1]);

	}
}

//=================================================================================
// micNetNode
//=================================================================================
//.........................................................................................................................................
int micNetNode::Forward()
{
	if (get_S_flag) {
		fast_multiply_medmat_(batch_in, wgt, S);
	}

	int rc = forward();

	return rc;
}

//=================================================================================
// micNodeLeakyReLU
//=================================================================================
//.........................................................................................................................................
int micNodeLeakyReLU::forward()
{
	for (int i=0; i<n_in; i++)
		for (int j=0; j<k_out; j++) {
			if (S(i, j) >= 0) {
				G(i, j) = ni.A;
				batch_out(i, j) = S(i, j)*ni.A;
			}
			else {
				G(i, j) = ni.B;
				batch_out(i, j) = S(i, j)*ni.B;
			}
		}

	return 0;
}
