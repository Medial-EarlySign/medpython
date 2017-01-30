#include "MedPlot.h"
#include <fstream>
#include <ctime>
#include <iostream>

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
	snprintf(res, sizeof(res) ,"%2.4f", num);
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

void createHtmlGraph(string outPath, vector<map<float, float>> data, string title, string xName, string yName, vector<string> seriesNames, int refreshTime)
{


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

	string rep = "";

	for (size_t i = 0; i < data.size(); ++i)
	{
		map<float, float> dmap = data[i];	

		rep += "var series" + to_string(i) + " = {\n mode: 'lines+markers',\n x: [";
		for (auto it = dmap.begin(); it != dmap.end(); ++it) {
			rep += float2Str(it->first) + ", ";
		}
		if (rep[rep.size() - 2] == ',') {
			rep = rep.substr(0, rep.size() - 2);
		}
		rep += "], \n";

		rep += "y: [";
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
		snprintf(buf, 100 ,"setTimeout(function() { window.location.reload(1); }, %d);", refreshTime);
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
	rep += title;
	rep += "', \n xaxis: { title : '";
	rep += xName;
	rep += "'}, \n yaxis: { title: '";
	rep += yName;
	rep += "' }, \n height: 800, \n    width: 1200 \n }; ";

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
