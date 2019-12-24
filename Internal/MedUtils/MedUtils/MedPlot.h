#ifndef __MED_PLOT_H__
#define __MED_PLOT_H__

#include <map>
#include <vector>
#include <string>

/** @file
* A Library to plot graphs in HTML using plotly.js \n
* Example Code: \n
* vector<map<float, float>> data(2); //plot 2 series of two lines \n
* \n
* //create data for lines: \n
* int numPnt = 1000; \n
* float m1 = 2;\n
* float n1 = 3;\n
* float m2 = -4;\n
* float n2 = 9;\n
* for (int i = 0; i < numPnt; ++i)\n
* {\n
* float x = (i - numPnt / 2) / float(100.0);\n
* float y = m1 * (i - numPnt / 2) / 100 + n1;\n
* data[0][x] = y;\n
*\n
* y = m2 * (i - numPnt / 2) / 100 + n2;\n
* data[1][x] = y;\n
* }\n
* //end creation of data, now  plot:\n
* \n
* vector<string> seriesNames = {"line_1", "line_2"};\n
* createHtmlGraph("test.html", data, "Graph Title", "x", "y", seriesNames);\n
*/

using namespace std;

extern string BaseResourcePath;
/// <summary>
/// prety print float number
/// </summary>
string float2Str(float num);

/// <summary>
///makes histogram for vector of numbers and stores it in map object 
/// </summary>
map<float, float> BuildHist(vector<float> featNums);

/// <summary>
/// proccess data to plot x,y. x is vector of signals and combineFeat will be used if this vector size > 1 to transform each row of signals into single number
/// that will be considred as X. aggFunction will be used to select which value of Y to return for each transformed X value - it could by mean, median, max, min, prctile..
/// </summary>
map<float, float> BuildAggeration(const vector<vector<float>> &vec_x, const vector<float> &y,
	float(*aggFunction)(const vector<float> &),
	float(*combineFeat)(const vector<float>&) = NULL);

/// <summary>
/// proccess data to plot x,y,z. x is vector of signals and combineFeat will be used if this vector size > 1 to transform each row of signals into single number
/// that will be considred as X. aggFunction will be used to select which value of Y to return for each transformed X value - it could by mean, median, max, min, prctile..
/// </summary>
void Build3Data(const vector<float> &x1, const vector<float> &x2,
	const vector<float> &y,
	float(*aggFunction)(const vector<float> &), vector<vector<float>> &data, int min_filter_cnt = 10);

/// <summary>
/// Will create Html Graph string - you will decide where to save it to disk. 
/// @param outPath the location to save the html file (recommend ending file ext with .html)
/// @param data is vector of series to plot with coresponding names in vector seriesNames. each element in the vector is series to plot represented by map<float, float> object
/// the plot will print the iteration on the keys with their corresponding values. the map object is used to store vector of tuples (x,y) to plot in each series
/// @param title graph title
/// @param xName x Axis name
/// @param yName y Axis name
/// @param seriesNames same size vector to data with coresponding labels to each data
/// @param refreshTime Time in milliseconds for the file to be refreshed by the browser (default 0, taken as do not refresh)
/// @param chart type Can be: "scatter", "bar", "pie"
/// </summary>
void createHtmlGraph(const string &outPath, const vector<map<float, float>> &data, const string &title = "", const string &xName = "", const string &yName = "",
	const vector<string> &seriesNames = vector<string>(), int refreshTime = 0, 
	const string &chart_type = "scatter", const string &mode = "lines", const string &template_str = "");

/// <summary>
/// Will create Html Graph string - you will decide where to save it to disk. 
/// @param outPath the location to save the html file (recommend ending file ext with .html)
/// @param data is vector of series to plot with coresponding names in vector seriesNames. each element in the vector is series to plot represented by map<float, float> object
/// the plot will print the iteration on the keys with their corresponding values. the map object is used to store vector of tuples (x,y) to plot in each series
/// @param title graph title
/// @param xName x Axis name
/// @param yName y Axis name
/// @param seriesNames same size vector to data with coresponding labels to each data
/// @param refreshTime Time in milliseconds for the file to be refreshed by the browser (default 0, taken as do not refresh)
/// @param chart type Can be: "scatter", "bar", "pie"
/// </summary>
void createScatterHtmlGraph(const string &outPath, const vector<vector<pair<float, float>>> &data, const string &title = "", 
	const string &xName = "", const string &yName = "", const vector<string> &seriesNames = vector<string>(), 
	int refreshTime = 0, const string &chart_type = "scatter", const string &mode = "markers", const string &template_str = "");

/// <summary>
/// Plot of 3D graph data
/// @param outPath The output file (recomanded html)
/// @param vec3d the 3d vector, first dim is vector of all series each sereis in
/// diffrent color. second dim is of size 3 for x,y,z axis data and thirds dim is the data in each axis
/// @param seriesNames same size vector to data with coresponding labels to each data 
/// @param heatmap - if true will print heatmap, else 3d graph
/// @param title - graph title
/// @param xName - the x axis name
/// @param yName - the y axis name
/// @param zName - the z axis name
/// </summary>
void createHtml3D(const string &outPath, const vector<vector<vector<float>>> &vec3d, const vector<string> &seriesNames,
	bool heatmap = true, const string &title = "", const string &xName = "x", const string &yName = "y", const string &zName = "z");

/// <summary>
/// returns a csv string content of all features with header name for each feature to save in csv format
/// </summary>
/// <returns>
/// returns the csv string content of all features with header name
/// </returns>
string createCsvFile(const map<float, float> &data);

/// <summary>
/// returns a csv string content of all features with header name for each feature to save in csv format
/// </summary>
/// <returns>
/// returns the csv string content of all features with header name
/// </returns>
string createCsvFile(const vector<vector<float>> &data, const vector<string> &headers);

/// <summary>
/// Down sampling the number of points in the graph to points_count if has more
/// points in the data. the interpulation is linear.
/// </summary>
void down_sample_graph(map<float, float> &points, int points_count = 10000);
extern vector<bool> empty_bool_arr;
/// <summary>
/// calculates true_rate, false_rate, ppv based on labels(y) and predictions scores(preds)
/// indexes is used for filtering samples
/// </summary>
/// <returns>
/// updates true_rate, false_rate, ppv
/// </returns>
void get_ROC_working_points(const vector<float> &preds, const vector<float> &y, const vector<float> &weights,
	vector<float> &pred_threshold, vector<float> &true_rate, vector<float> &false_rate, vector<float> &ppv, vector<float> &pr,
	const vector<bool> &indexes = empty_bool_arr);
/// <summary>
/// plot AUC Graph for all scores and each score has diffrent label size
/// </summary>
void plotAUC(const vector<vector<float>> &all_preds, const vector<vector<float>> &y, const vector<vector<float>> &weights, const vector<string> &modelNames,
	string baseOut, bool print_y = true);
/// <summary>
/// Plot AUC Graph for all scores for same labels
/// </summary>
void plotAUC(const vector<vector<float>> &all_preds, const vector<float> &y, const vector<string> &modelNames,
	string baseOut, const vector<bool> &indexes = empty_bool_arr, const vector<float> *weights = NULL);
#endif