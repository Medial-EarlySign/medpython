#ifndef __MED_PLOT_H__
#define __MED_PLOT_H__

#include <map>
#include <vector>
#include <string>

/* Example Code:

vector<map<float, float>> data(2); //plot 2 series of two lines

//create data for lines:
int numPnt = 1000;
float m1 = 2;
float n1 = 3;
float m2 = -4;
float n2 = 9;
for (int i = 0; i < numPnt; ++i)
{
float x = (i - numPnt / 2) / float(100.0);
float y = m1 * (i - numPnt / 2) / 100 + n1;
data[0][x] = y;

y = m2 * (i - numPnt / 2) / 100 + n2;
data[1][x] = y;
}
//end creation of data, now  plot:

vector<string> seriesNames = {"line_1", "line_2"};
createHtmlGraph("test.html", data, "Graph Title", "x", "y", seriesNames);

*/

using namespace std;

extern string BaseResourcePath;
//prety print for float number
string float2Str(float num);

//makes histogram for vector of numbers and stores it in map object 
map<float, float> BuildHist(vector<float> featNums);

//proccess data to plot x,y. x is vector of signals and combineFeat will be used if this vector size > 1 to transform each row of signals into single number
// that will be considred as X. aggFunction will be used to select which value of Y to return for each transformed X value - it could by mean, median, max, min, prctile..
map<float, float> BuildAggeration(const vector<vector<float>> &vec_x, const vector<float> &y,
	float(*aggFunction)(const vector<float> &),
	float(*combineFeat)(const vector<float>&) = NULL);

void Build3Data(const vector<float> &x1, const vector<float> &x2,
	const vector<float> &y,
	float(*aggFunction)(const vector<float> &), vector<vector<float>> &data, int min_filter_cnt = 10);

//Will create Html Graph string - you will decide where to save it to disk. outPath - the location to save the html file (recommend ending file ext with .html)
// data - is vector of series to plot with coresponding names in vector seriesNames. each element in the vector is series to plot represented by map<float, float> object
// the plot will print the iteration on the keys with their corresponding values. the map object is used to store vector of tuples (x,y) to plot in each series
// title - graph title, xName-  xAxis name, yName - y axis name
//refreshTime - time in milliseconds for the file to be refreshed by the browser (default 0, taken as do not refresh)
void createHtmlGraph(string outPath, vector<map<float, float>> data, string title = "", string xName = "", string yName = "", vector<string> seriesNames = vector<string>(), int refreshTime = 0);

void createHtml3D(string outPath, const vector<vector<float>> &vec3d, bool heatmap = true, string title = "", string xName = "x", string yName = "y", string zName = "z");

//returns a csv string content of the map data that represents X->Y function. will return Csv text with header X,Y and the map values
string createCsvFile(map<float, float> data);

//returns a csv string content of all features with header name for each feature to save in csv format
string createCsvFile(vector<vector<float>> &data, const vector<string> &headers);

#endif