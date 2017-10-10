#include "MedPlot.h"
#include <fstream>
#include <ctime>
#include <iostream>
#include <boost/filesystem/path.hpp>

#ifdef _WIN32
//define something for Windows
string BaseResourcePath = "W:\\Graph_Infra";
#else
string BaseResourcePath = "/nas1/Work/Graph_Infra";
#endif

using namespace std;

inline char separator()
{
#if defined _WIN32 || defined __CYGWIN__
	return '\\';
#else
	return '/';
#endif
}

string float2Str(float num) {
	//return to_string((round(num * 1000) / 1000));
	char res[50];
	snprintf(res, sizeof(res), "%2.4f", num);
	return string(res);
}

map<float, float> BuildHist(vector<float> featNums) {
	map<float, float> hist;
	for (float n : featNums)
	{
		if (hist.find(n) == hist.end()) {
			hist[n] = 0;
		}
		++hist[n];
	}
	return hist;
}

map<float, float> BuildAggeration(const vector<vector<float>> &vec_x, const vector<float> &y,
	float(*aggFunction)(const vector<float> &),
	float(*combineFeat)(const vector<float>&)) {
	map<float, float> res;
	map<float, vector<float>> data;

	for (size_t i = 0; i < vec_x.size(); ++i)
	{
		if (y.size() != vec_x[i].size()) {
			std::cout << "y_size=" << y.size() << " vec_" << i << "_size=" << vec_x[i].size() << endl;
			throw logic_error("not all feature vectors has same length of target. there is no matching between them");
		}
	}

	vector<float> rowData(vec_x.size());
	float bucket;
	for (int i = 0; i < y.size(); ++i)
	{
		for (int k = 0; k < vec_x.size(); ++k) {
			rowData[k] = vec_x[k][i];
		}
		if (combineFeat == NULL && vec_x.size() == 1) {
			bucket = vec_x[0][i];
		}
		else {
			bucket = combineFeat(rowData);
		}

		data[bucket].push_back(y[i]);
	}

	for (auto it = data.begin(); it != data.end(); ++it) {
		res[it->first] = aggFunction(it->second);
	}

	return res;
}

void Build3Data(const vector<float> &x1, const vector<float> &x2,
	const vector<float> &y,
	float(*aggFunction)(const vector<float> &), vector<vector<float>> &data, int min_filter_cnt) {
	//aggregate for each tuples of x1,x2 aggFucntion on y list results
	if (x1.size() != x2.size() || x1.size() != y.size()) {
		throw invalid_argument("arrays must have same size");
	}
	data = vector<vector<float>>(3);
	map<float, map<float, vector<float>>> d;
	for (size_t i = 0; i < x1.size(); ++i)
	{
		d[x1[i]][x2[i]].push_back(y[i]);
	}
	for (auto it = d.begin(); it != d.end(); ++it)
	{
		for (auto jt = it->second.begin(); jt != it->second.end(); ++jt)
		{
			if (jt->second.size() < min_filter_cnt)
				continue; //filtered out
			data[0].push_back(it->first);
			data[1].push_back(jt->first);
			data[2].push_back(aggFunction(jt->second));
		}
	}

	if (data[0].size() == 0)
		throw invalid_argument("filtered all points - min_filter_cnt is too high or axis bining is needed");
}

void createHtmlGraph(string outPath, vector<map<float, float>> data, string title, string xName, string yName, 
	vector<string> seriesNames, int refreshTime, string chart_type)
{
	string x_name = "x";
	string y_name = "y";
	if (chart_type == "pie") {
		x_name = "labels";
		y_name = "values";
	}

	/*ifstream jsFile;
	jsFile.open(BaseResourcePath + separator() + "plotly-latest.min.js");
	if (!jsFile.is_open()) {
		throw logic_error("Unable to open js file");
	}
	string jsData((istreambuf_iterator<char>(jsFile)),
		istreambuf_iterator<char>());
	ofstream jsOut;*/
	size_t lastDirPos = outPath.find_last_of("/\\");
	boost::filesystem::path p(outPath);
	boost::filesystem::path outDir = p.parent_path();

	/*cerr << "writing: [" << outDir.string() + "/plotly-latest.min.js" << "]\n";
	jsOut.open(outDir.string() + "/plotly-latest.min.js");
	jsOut << jsData;
	jsOut.close();*/

	ifstream file(BaseResourcePath + separator() + "Graph_HTML.txt");
	if (!file.is_open()) {
		throw logic_error("Unable to open file");
	}
	string ln;
	string content = "";
	while (getline(file, ln)) {
		content += "\n" + ln;
	}
	file.close();

	size_t ind = content.find("{0}");
	if (ind == string::npos) {
		throw invalid_argument("Not Found in template");
	}

	string rep = "";

	for (size_t i = 0; i < data.size(); ++i)
	{
		map<float, float> dmap = data[i];

		rep += "var series" + to_string(i) + " = {\n type: '" + chart_type + "',\n mode: 'lines+markers',\n " + x_name +": [";
		for (auto it = dmap.begin(); it != dmap.end(); ++it) {
			rep += float2Str(it->first) + ", ";
		}
		if (rep[rep.size() - 2] == ',') {
			rep = rep.substr(0, rep.size() - 2);
		}
		rep += "], \n";

		rep += y_name + ": [";
		for (auto it = dmap.begin(); it != dmap.end(); ++it) {
			rep += float2Str(it->second) + ", ";
		}
		if (rep[rep.size() - 2] == ',') {
			rep = rep.substr(0, rep.size() - 2);
		}
		rep += "] \n";
		if (seriesNames.size() > 0) {
			if (seriesNames.size() != data.size()) {
				throw invalid_argument("wrong number of series names passed with data graphs to plot");
			}
			rep += ", name: '";
			rep += seriesNames[i];
			rep += "' \n";
		}
		rep += "};\n";
	}

	if (refreshTime > 0) {
		char buf[100];
		snprintf(buf, 100, "setTimeout(function() { window.location.reload(1); }, %d);", refreshTime);
		rep += buf;
	}


	rep += "var data = [";
	for (size_t i = 0; i < data.size(); ++i)
	{
		rep += " series" + to_string(i) + ", ";
	}
	rep = rep.substr(0, rep.size() - 2);
	rep += " ]; \n";

	rep += "var layout = { \n  title:'";
	rep += title + "',\n ";
	if (chart_type != "pie") {
		rep += "xaxis: { title : '";
		rep += xName;
		rep += "'}, \n yaxis: { title: '";
		rep += yName + "'},\n ";
	}
	
	rep += "height: 800, \n    width: 1200 \n }; ";

	content.replace(ind, 3, rep);
	content.replace(content.find("\"plotly-latest.min.js\""), 22, "\"W:\\Graph_Infra\\plotly-latest.min.js\"");

	ofstream myfile;
	cerr << "writing: [" << outPath << "]\n";
	myfile.open(outPath);
	myfile << content;
	myfile.close();
}

void createHtml3D(string outPath, const vector<vector<float>> &vec3d, bool heatmap, string title, string xName, string yName, string zName) {
	if (vec3d.size() != 3) {
		throw invalid_argument("please pass 3 signal vectors as input");
	}
	vector<string> ind2axis = { "x", "y", "z" };

	ifstream jsFile;
	jsFile.open(BaseResourcePath + separator() + "plotly-latest.min.js");
	if (!jsFile.is_open()) {
		throw logic_error("Unable to open js file");
	}
	string jsData((istreambuf_iterator<char>(jsFile)),
		istreambuf_iterator<char>());
	ofstream jsOut;
	size_t lastDirPos = outPath.find_last_of("/\\");
	string outDir = outPath.substr(0, lastDirPos) + separator();
	if (lastDirPos == string::npos)
	{
		outDir = "";
	}

	jsOut.open(outDir + "plotly-latest.min.js");
	jsOut << jsData;
	jsOut.close();

	ifstream file(BaseResourcePath + separator() + "Graph_HTML.txt");
	if (!file.is_open()) {
		throw logic_error("Unable to open file");
	}
	string ln;
	string content = "";
	while (getline(file, ln)) {
		content += "\n" + ln;
	}
	file.close();

	size_t ind = content.find("{0}");
	if (ind == string::npos) {
		throw invalid_argument("Not Found in template");
	}

	string type = "scatter3d";
	if (heatmap)
		type = "heatmap";

	string rep = "";
	rep += "var series" + to_string(0) + " = {\n type: '" + type + "', \n mode: 'markers'";
	for (size_t i = 0; i < vec3d.size(); ++i) {
		rep += ",\n" + ind2axis[i] + ": [";
		rep += float2Str(vec3d[i][0]);
		for (size_t j = 1; j < vec3d[i].size(); ++j)
		{
			rep += ", " + float2Str(vec3d[i][j]);
		}

		rep += "]";
	}
	rep += "\n";
	rep += ", name: '";
	rep += "series0";
	rep += "' \n";
	rep += "};\n";

	rep += "var data = [";
	for (size_t i = 0; i < 1; ++i)
	{
		rep += " series" + to_string(i) + ", ";
	}
	rep = rep.substr(0, rep.size() - 2);
	rep += " ]; \n";

	rep += "var layout = { \n  title:'";
	rep += title;
	if (!heatmap)
		rep += "', \n scene: {\n xaxis: { title : '";
	else
		rep += "', \n xaxis: { title : '";
	rep += xName;
	rep += "'}, \n yaxis: { title: '";
	rep += yName;
	rep += "'}, \n zaxis: { title: '";
	rep += zName;
	if (!heatmap)
		rep += "' }\n }, \n height: 800, \n    width: 1200 \n }; ";
	else
		rep += "'}, \n height: 800, \n    width: 1200 \n }; ";

	content.replace(ind, 3, rep);

	ofstream myfile;
	myfile.open(outPath);
	myfile << content;
	myfile.close();
}

string createCsvFile(map<float, float> data)
{
	string out = "X, Y\n";
	for (auto it = data.begin(); it != data.end(); ++it)
	{
		out += to_string(it->first) + ", " + to_string(it->second) + "\n";
	}

	return out;
}

string createCsvFile(vector<vector<float>> &data, const vector<string> &headers)
{
	string out = "ROWS";
	for (size_t i = 0; i < headers.size(); ++i)
	{
		out += ", " + headers[i];
	}
	out += "\n";

	for (size_t i = 0; i < data.size(); ++i)
	{
		out += "ROW_" + to_string(i);
		for (size_t k = 0; k < data[i].size(); ++k)
		{
			out += ", " + float2Str(data[i][k]);
		}
		out += "\n";

	}

	return out;
}
