// This is the main DLL file.
#ifdef _WIN32
#pragma warning(disable : 4996)
#endif

#include "SigSum_AlgoMarker.h"
// This is the main DLL file.

#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <iterator>
#include <fstream>
#include <cstring>
#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>

#ifdef __linux__
#include <wordexp.h>
#define AM_LOG_CONTROL_ENV_VAR "${AM_LOG}"
#define AM_LOG_FILE_ENV_VAR "${AM_LOG_FILE}"
#define AM_LOG__DEFAULT_FILE "/var/log/signalssum/signalssum.log"
#elif _WIN32
#include "windows.h"
#define AM_LOG_CONTROL_ENV_VAR "%AM_LOG%"
#define AM_LOG_FILE_ENV_VAR "%AM_LOG_FILE%"
#define AM_LOG__DEFAULT_FILE "%TEMP%\\AlgoMarkers_Log.txt"
#endif

// #include <climits>
// using namespace std;

string expandEnvVars(const string &str)
{
	string ret = "";
#ifdef __linux__
	wordexp_t p;
	char **w;
	wordexp(str.c_str(), &p, 0);
	w = p.we_wordv;
	for (size_t i = 0; i < p.we_wordc; i++)
		ret += w[i];
	wordfree(&p);
#elif _WIN32
	DWORD max_str_len = 4 * 1024;
	auto buf = new char[max_str_len];
	DWORD req_len = ExpandEnvironmentStrings(str.c_str(), buf, max_str_len);
	if (req_len > max_str_len)
	{
		delete buf;
		buf = new char[req_len];
		req_len = ExpandEnvironmentStrings(str.c_str(), buf, req_len);
	}
	if (req_len > 0)
		ret = buf;
	delete buf;
#endif
	return ret;
}

#define LOCAL_SECTION LOG_APP
#define LOCAL_LEVEL LOG_DEF_LEVEL
//-----------------------------------------------------------------------------------
void AMMessages::get_messages(int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	if (need_to_update_args)
	{
		args.clear();
		for (auto &s : args_strs)
			args.push_back((char *)s.c_str());
		need_to_update_args = 0;
	}

	*n_msgs = get_n_msgs();
	if (*n_msgs > 0)
	{
		*msgs_codes = &codes[0];
		*msgs_args = &args[0];
	}
	else
	{
		*msgs_codes = NULL;
		*msgs_args = NULL;
	}
}

//-----------------------------------------------------------------------------------
void AMMessages::insert_message(int code, const char *arg_ch)
{
	string arg;
	if (arg_ch == NULL)
	{
		arg_ch = "";
	}
	else
	{
		arg = string(arg_ch);
	}

	codes.push_back(code);
	args_strs.push_back(arg);
	need_to_update_args = 1;

	////string arg;
	// codes.push_back(code);
	////args_strs.push_back(arg);
	// if(arg_ch == NULL)
	//{
	//	args.push_back(NULL);
	//	return;
	// }
	// const size_t strSize = strlen(arg_ch) + 1;
	// char *cstr = new char[strSize];
	// strcpy_s(cstr, strSize, arg_ch);
	// args.push_back(cstr);
}

////-----------------------------------------------------------------------------------
// if does not exist returns -1.
int AMResponses::get_response_index_by_point(int _pid, long long _timestamp)
{
	// pair<int, long long> p(_pid, _timestamp);

	if (_pid != 1) // point2response_idx.find(p) == point2response_idx.end())
		return LogData::EndFunction(AM_FAIL_RC, __func__);

	return 0;
}
//
////-----------------------------------------------------------------------------------
//// if does not exist returns NULL
// AMResponse *AMResponses::get_response_by_point(int _pid, long long _timestamp)
//{
//	pair<int, long long> p(_pid, _timestamp);
//
//	if (point2response_idx.find(p) == point2response_idx.end())
//		return NULL;
//
//	return &responses[point2response_idx[p]];
// }

//-----------------------------------------------------------------------------------
void AMResponses::get_score_types(int *n_score_types, char ***_score_types)
{
	*n_score_types = (int)score_types.size();
	if (n_score_types == 0)
		*_score_types = NULL;
	else
		*_score_types = &score_types[0];
}

//-----------------------------------------------------------------------------------
int AMResponses::get_score(int _pid, long long _timestamp, char *_score_type, float *out_score)
{
	pair<int, long long> p(_pid, _timestamp);

	if (_pid != 1)
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	// int pidx = point2response_idx[p];

	return get_score_by_type(0, _score_type, out_score);
}

//-----------------------------------------------------------------------------------
int AMResponses::get_score_by_type(int index, char *_score_type, float *out_score)
{
	string s = string(_score_type);

	if (index < 0 || index >= get_n_responses())
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	// if (stype2idx.find(s) == stype2idx.end())
	//	return LogData::EndFunction(AM_FAIL_RC, __func__);
	// int sidx = stype2idx[s];
	char *dummy_type;

	if (responses[index].get_score(0, out_score, &dummy_type) != AM_OK_RC)
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	if (strcmp(dummy_type, _score_type) != 0)
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
void AMResponses::insert_score_types(char **_score_type, int n_score_types)
{
	for (int i = 0; i < n_score_types; i++)
	{
		string s = string(_score_type[i]);
		score_types_str.push_back(s);
		// stype2idx[s] = (int)score_types.size() - 1;
	}

	for (int i = 0; i < n_score_types; i++)
		score_types.push_back((char *)score_types_str[i].c_str());
}

//-----------------------------------------------------------------------------------
AMResponse *AMResponses::create_point_response(int _pid, long long _timestamp)
{
	pair<int, long long> p(_pid, _timestamp);

	AMResponse response;

	response.set_patient_id(_pid);
	response.set_timestamp(_timestamp);
	response.init_scores((int)score_types.size());

	responses.push_back(response);

	// point2response_idx[p] = (int)responses.size() - 1;

	return &responses.back();
}

//-----------------------------------------------------------------------------------
int AlgoMarker::IsScoreTypeSupported(const char *_stype)
{
	string stype = string(_stype);

	for (auto &s : supported_score_types)
		if (stype == s)
			return 1;
	return 0;
}

//-----------------------------------------------------------------------------------
AlgoMarker *AlgoMarker::make_algomarker(AlgoMarkerType am_type)
{

	if (am_type == AM_TYPE_MEDIAL_INFRA)
		return new MedialInfraAlgoMarker;

	return NULL;
}

#define SWITCH_CHANGE_PARAMETER "ChangeParameter"

#define SWITCH_ARRAY_EMPTY_ARRAY "array_with_empty_array"
#define SWITCH_ARRAY_EMPTY_STRING "array_with_empty_string"
#define SWITCH_ARRAY_WITH_NULL "array_with_null_array"

#define SWITCH_EMPTY_ARRAY "empty_array"
#define SWITCH_EMPTY_STRING "empty_string"
#define SWITCH_NULL "null"

#define SWITCH_SET_SCORE "SET_SCORE"
#define SWITCH_MSG "MSG"
#define SWITCH_RMSG "RMSG"
#define SWITCH_SMSG "SMSG"
#define SWITCH_MSG2 "911"
#define SWITCH_RMSG2 "9112"
#define SWITCH_SMSG2 "9113"

bool LogData::TurnOnLogs = false;
vector<string> LogData::Arguments;
const char *LogData::FunctionName;
ofstream *LogData::logFileStream = nullptr;
bool LogData::DoesEnterFunction;

void LogData::WriteToLog(string line)
{
	// printf("LogData::TurnOnLogs = %d\n", (int)LogData::TurnOnLogs);
	if (!LogData::TurnOnLogs)
		return;
	if (logFileStream == nullptr)
	{
		string filename = expandEnvVars(AM_LOG_FILE_ENV_VAR);
		if (filename == "")
			filename = expandEnvVars(AM_LOG__DEFAULT_FILE);
		// printf("filename = '%s'\n", filename.c_str());
		// printf("AM_LOG_FILE_ENV_VAR = '%s'\n", expandEnvVars(AM_LOG_FILE_ENV_VAR).c_str());
		// printf("AM_LOG__DEFAULT_FILE = '%s'\n", expandEnvVars(AM_LOG__DEFAULT_FILE).c_str());
		logFileStream = new ofstream();
		logFileStream->open(filename, ios::app);
	}
	(*logFileStream) << line + "\r\n";
}

void LogData::closeLogFile()
{
	if (logFileStream != nullptr)
	{
		logFileStream->flush();
		logFileStream->close();
		logFileStream = nullptr;
	}
}

void LogData::FlushLogData()
{
	if (!LogData::TurnOnLogs)
		return;

	string line = "";

	if (LogData::Arguments.size() > 0)
	{
		line += "(";

		for (size_t argumrntIndex = 0; argumrntIndex < LogData::Arguments.size(); argumrntIndex++)
		{
			line += LogData::Arguments[argumrntIndex];
			if (argumrntIndex != LogData::Arguments.size() - 1)
			{
				line += "; ";
			}
		}
		line += ")";
	}
	LogData::Arguments.clear();

	if (DoesEnterFunction)
	{
		line = "Enter Function " + string(LogData::FunctionName) + ". " + line;
	}
	else
	{
		line = "End Function " + string(LogData::FunctionName) + ". " + line + "\r\n";
	}

	WriteToLog(line);
}

void LogData::StartFunction(const char *function)
{
	if (!TurnOnLogs)
		return;
	LogData::DoesEnterFunction = true;
	LogData::FunctionName = function;
	FlushLogData();
}

void LogData::EndFunction(const char *function)
{
	if (!TurnOnLogs)
		return;
	LogData::DoesEnterFunction = false;
	LogData::FunctionName = function;
	FlushLogData();
}

int LogData::StartFunction(int returnValue, const char *function)
{
	if (!LogData::TurnOnLogs)
		return returnValue;
	LogData::StartFunction(function);
	return returnValue;
}

int LogData::EndFunction(int returnValue, const char *function)
{
	if (!LogData::TurnOnLogs)
		return returnValue;
	LogData::AddArgument("return", returnValue);
	LogData::EndFunction(function);
	return returnValue;
}

void LogData::AddArgument(string argumentName, string argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::Arguments.push_back("\"" + argumentName + "\":" + argument);
}
void LogData::AddArgument(string argumentName, char *argument)
{
	LogData::AddArgument(argumentName, string(argument));
}
void LogData::AddArgument(string argumentName, int *argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, "[" + to_string(argument[0]) + "]");
}
void LogData::AddArgument(string argumentName, float *argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, "[" + to_string(argument[0]) + "]");
}
void LogData::AddArgument(string argumentName, long long *argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, "[" + to_string(argument[0]) + "]");
}
void LogData::AddArgument(string argumentName, int argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, to_string(argument));
}
void LogData::AddArgument(string argumentName, unsigned long long argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, to_string(argument));
}
void LogData::AddArgument(string argumentName, float argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, to_string(argument));
}
void LogData::AddArgument(string argumentName, long long argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, to_string(argument));
}
void LogData::AddArgument(string argumentName, char **argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, argument[0]);
}
void LogData::AddArgument(string argumentName, AlgoMarker *argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument);
}
void LogData::AddArgument(string argumentName, AlgoMarker **argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument[0]);
}
void LogData::AddArgument(string argumentName, AMRequest *argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument);
}
void LogData::AddArgument(string argumentName, AMRequest **argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument[0]);
}
void LogData::AddArgument(string argumentName, AMResponses *argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument);
}
void LogData::AddArgument(string argumentName, AMResponses **argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument[0]);
}
void LogData::AddArgument(string argumentName, AMResponse *argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument);
}
void LogData::AddArgument(string argumentName, AMResponse **argument)
{
	if (!LogData::TurnOnLogs)
		return;
	LogData::AddArgument(argumentName, (unsigned long long)argument[0]);
}
void LogData::AddArgument(string argumentName, int **argument, int arraySize)
{
	if (!LogData::TurnOnLogs)
		return;
	string argumentsString = "[ [";
	for (size_t argumrntIndex = 0; argumrntIndex < arraySize; argumrntIndex++)
	{
		argumentsString += to_string(argument[0][argumrntIndex]);
		if (argumrntIndex != arraySize - 1)
		{
			argumentsString += ", ";
		}
	}
	argumentsString += "] ]";

	LogData::AddArgument(argumentName, argumentsString);
}
void LogData::AddArgument(string argumentName, char ***argument, int arraySize)
{
	if (!LogData::TurnOnLogs)
		return;
	string argumentsString = "[ [";
	for (size_t argumrntIndex = 0; argumrntIndex < arraySize; argumrntIndex++)
	{
		argumentsString += string(argument[0][argumrntIndex]);
		if (argumrntIndex != arraySize - 1)
		{
			argumentsString += ", ";
		}
	}
	argumentsString += "] ]";
	LogData::AddArgument(argumentName, argumentsString);
}
void LogData::AddArgument(string argumentName, char **argument, int arraySize)
{
	if (!LogData::TurnOnLogs)
		return;
	string argumentsString = "[";
	for (size_t argumrntIndex = 0; argumrntIndex < arraySize; argumrntIndex++)
	{
		argumentsString += string(argument[argumrntIndex]);
		if (argumrntIndex != arraySize - 1)
		{
			argumentsString += ", ";
		}
	}
	argumentsString += "]";
	LogData::AddArgument(argumentName, argumentsString);
}
void LogData::AddArgument(string argumentName, long long *argument, int arraySize)
{
	if (!LogData::TurnOnLogs)
		return;
	string argumentsString = "[";
	for (size_t argumrntIndex = 0; argumrntIndex < arraySize; argumrntIndex++)
	{
		argumentsString += to_string(argument[argumrntIndex]);
		if (argumrntIndex != arraySize - 1)
		{
			argumentsString += ", ";
		}
	}
	argumentsString += "]";
	LogData::AddArgument(argumentName, argumentsString);
}
void LogData::AddArgument(string argumentName, int *argument, int arraySize)
{
	if (!LogData::TurnOnLogs)
		return;
	string argumentsString = "[";
	for (size_t argumrntIndex = 0; argumrntIndex < arraySize; argumrntIndex++)
	{
		argumentsString += to_string(argument[argumrntIndex]);
		if (argumrntIndex != arraySize - 1)
		{
			argumentsString += ", ";
		}
	}
	argumentsString += "]";
	LogData::AddArgument(argumentName, argumentsString);
}
void LogData::AddArgument(string argumentName, float *argument, int arraySize)
{
	if (!LogData::TurnOnLogs)
		return;
	string argumentsString = "[";
	for (size_t argumrntIndex = 0; argumrntIndex < arraySize; argumrntIndex++)
	{
		argumentsString += to_string(argument[argumrntIndex]);
		if (argumrntIndex != arraySize - 1)
		{
			argumentsString += ", ";
		}
	}
	argumentsString += "]";
	LogData::AddArgument(argumentName, argumentsString);
}

void ChangeParameter(vector<string> config, const char *function, int *parameter)
{
	if (strcmp(config[0].c_str(), function) != 0)
		return;
	if (strcmp(config[1].c_str(), SWITCH_CHANGE_PARAMETER) != 0)
		return;

	if (strcmp(config[2].c_str(), "return") != 0)
		return;

	stringstream geek(config[3]);
	int x = 0;
	geek >> x;
	*parameter = x;
}
int SpecialReturn(vector<string> config, const char *function, int returnValue = AM_OK_RC)
{
	ChangeParameter(config, function, &returnValue);
	LogData::EndFunction(returnValue, function);
	return returnValue;
}
int SpecialReturn(AlgoMarker *ref, const char *function, int returnValue = AM_OK_RC)
{
	if (ref == NULL)
		return returnValue;
	return SpecialReturn(ref->verificationConfig, function, returnValue);
}
int SpecialReturn(AMResponse *ref, const char *function, int returnValue = AM_OK_RC)
{
	if (ref == NULL)
		return returnValue;
	return SpecialReturn(ref->verificationConfig, function, returnValue);
}
int SpecialReturn(AMRequest *ref, const char *function, int returnValue = AM_OK_RC)
{
	if (ref == NULL)
		return returnValue;
	return SpecialReturn(ref->verificationConfig, function, returnValue);
}
int SpecialReturn(AMResponses *ref, const char *function, int returnValue = AM_OK_RC)
{
	if (ref == NULL)
		return returnValue;
	return SpecialReturn(ref->verificationConfig, function, returnValue);
}

void ChangeParameter(vector<string> config, const char *function, char **parameter)
{
	if (config.size() == 0)
	{
		cout << "Error: ChangeParameter got config size = 0\n";
		return;
	}

	if (strcmp(config[0].c_str(), function) != 0)
		return;

	if (config.size() <= 1)
	{
		cout << "Error: ChangeParameter got config size <=1\n";
		return;
	}

	if (strcmp(config[1].c_str(), SWITCH_CHANGE_PARAMETER) != 0)
		return;

	if (config.size() <= 2)
	{
		cout << "Error: ChangeParameter got config size <=2\n";
		return;
	}

	if (strcmp(config[2].c_str(), "name") == 0 ||
		strcmp(config[2].c_str(), "msgs_args") == 0 ||
		strcmp(config[2].c_str(), "_score_type") == 0 ||
		strcmp(config[2].c_str(), "requestId") == 0)
	{
		if (config.size() <= 4)
		{
			cout << "Error: ChangeParameter got config size <=4\n";
			return;
		}

		if (strcmp(config[3].c_str(), SWITCH_NULL) == 0 ||
			strcmp(config[3].c_str(), SWITCH_ARRAY_WITH_NULL) == 0)
		{
			*parameter = NULL;
		}
		else if (strcmp(config[3].c_str(), SWITCH_EMPTY_ARRAY) == 0 ||
				 strcmp(config[3].c_str(), SWITCH_ARRAY_EMPTY_ARRAY) == 0)
		{
			*parameter = new char[0];
		}
		else if (strcmp(config[3].c_str(), SWITCH_EMPTY_STRING) == 0 ||
				 strcmp(config[3].c_str(), SWITCH_ARRAY_EMPTY_STRING) == 0)
		{
			*parameter = "";
		}
		else
		{
			const size_t strSize = config[3].length() + 1;
			char *cstr = new char[strSize];
			strncpy(cstr, config[3].c_str(), strSize);
			cstr[strSize - 1] = 0;
			*parameter = cstr;
		}
	}
}

void ChangeParameter(vector<string> config, const char *function, char ***parameter)
{
	if (strcmp(config[0].c_str(), function) != 0)
		return;
	if (strcmp(config[1].c_str(), SWITCH_CHANGE_PARAMETER) != 0)
		return;

	if (strcmp(config[2].c_str(), "msgs_args") == 0)
	{
		if (strcmp(config[3].c_str(), SWITCH_NULL) == 0)
		{
			*parameter = NULL;
		}
		else if (strcmp(config[3].c_str(), SWITCH_EMPTY_ARRAY) == 0)
		{
			*parameter = new char *[0];
		}
		else
		{
			ChangeParameter(config, function, *parameter);
		}
	}
}

bool replace(std::string &str, const std::string &from, const std::string &to)
{
	size_t start_pos = str.find(from);
	if (start_pos == std::string::npos)
		return false;
	str.replace(start_pos, from.length(), to);
	return true;
}

//===========================================================================================================
//===========================================================================================================
// MedialInfraAlgoMarker Implementations ::
// Follows is an implementation of an AlgoMarker , which basically means filling in the:
// Load , Unload, ClearData, AddData and Calculate APIs. ( + private internal functions)
// This specific implementation uses medial internal infrastructure for holding data, models, and getting
// predictions.
//===========================================================================================================
//===========================================================================================================
//-----------------------------------------------------------------------------------
// Load() - reading a config file and initializing repository and model
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Load(const char *config_f)
{
	return SpecialReturn(this, __func__);
}

//------------------------------------------------------------------------------------------------
// UnLoad() - clears all data, repository and model, making object ready to be deleted and freed
//------------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Unload()
{
	return SpecialReturn(this, __func__);
}

//-----------------------------------------------------------------------------------
// ClearData() - clearing current data inserted inside.
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::ClearData()
{
	for (int i = 0; i < data.size(); i++)
	{
		data[i]->clear();
	}
	data.clear();
	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// AddData() - adding data for a signal with values and timestamps
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddData(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, float *Values)
{
	// At the moment MedialInfraAlgoMarker only loads timestamps given as ints.
	// This may change in the future as needed.
	AMData *newData = new AMData();
	newData->patient_id = patient_id;
	newData->signalName = string(signalName);
	newData->TimeStamps_len = TimeStamps_len;
	newData->TimeStamps = new long long[TimeStamps_len];
	for (size_t i = 0; i < TimeStamps_len; i++)
	{
		newData->TimeStamps[i] = TimeStamps[i];
	}

	newData->Values = new float[Values_len];
	for (size_t i = 0; i < Values_len; i++)
	{
		newData->Values[i] = Values[i];
	}
	newData->Values_len = Values_len;

	data.push_back(newData);

	// int sizeOfData = data.size();
	// data.resize(sizeOfData + 1);
	// data[sizeOfData] = newData;

	return AM_OK_RC;
}

//-----------------------------------------------------------------------------------
// AddDataStr() - adding data for a signal with values and timestamps
//-----------------------------------------------------------------------------------
int MedialInfraAlgoMarker::AddDataStr(int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, char **Values)
{
	// At the moment MedialInfraAlgoMarker only loads timestamps given as ints.
	// This may change in the future as needed.
	AMDataStr *newData = new AMDataStr();
	newData->patient_id = patient_id;
	newData->signalName = string(signalName);
	newData->TimeStamps_len = TimeStamps_len;
	newData->TimeStamps = new long long[TimeStamps_len];
	for (size_t i = 0; i < TimeStamps_len; i++)
	{
		newData->TimeStamps[i] = TimeStamps[i];
	}
	newData->Values = new char*[Values_len];
	// To send to AddData
	float *values = new float[Values_len];
	for (size_t i = 0; i < Values_len; i++)
	{
		values[i] = 0.0;
	}
	for (size_t i = 0; i < Values_len; i++)
	{
		newData->Values[i] = Values[i];
		std::string stringValue(Values[i]);
		for (char &c : stringValue)
		{
			values[i] = values[i] + (int)c;
		}
	}
	newData->Values_len = Values_len;

	AddData(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, values);
	// data.push_back(newData);

	// int sizeOfData = data.size();
	// data.resize(sizeOfData + 1);
	// data[sizeOfData] = newData;

	return AM_OK_RC;
}

//------------------------------------------------------------------------------------------
// Calculate() - after data loading : get a request, get predictions, and pack as responses
//------------------------------------------------------------------------------------------
int MedialInfraAlgoMarker::Calculate(AMRequest *request, AMResponses *responses)
{

	responses->set_request_id(request->get_request_id());
	for (int i = 0; i < request->get_n_score_types(); i++)
	{
		char *stype = request->get_score_type(i);
		responses->insert_score_types(&stype, 1);
	}

	AMMessages *shared_msgs = responses->get_shared_messages();
	// shared_msgs->insert_message(390, __func__);

	if (request->get_n_score_types() != 1)
	{
		shared_msgs->insert_message(AM_GENERAL_FATAL, "Verification Error: More than one score Types");
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
	if (request == NULL)
	{
		string msg = "(" + to_string(AM_MSG_NULL_REQUEST) + " ) NULL request in Calculate()";
		shared_msgs->insert_message(AM_GENERAL_FATAL, msg.c_str());
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}

	for (int i = 0; i < this->verificationConfig.size(); i++)
	{
		request->verificationConfig.push_back(this->verificationConfig[i]);
		responses->verificationConfig.push_back(this->verificationConfig[i]);
	}

	// We now have to prepare samples for the requested points
	// again - we only deal with int times in this class, so we convert the long long stamps to int
	/*int n_points =*/request->get_n_points();

	int pid = 1;
	if (data.size() > 0)
	{
		data.front()->patient_id;
	}
	// string msg = "Shared message text";
	// shared_msgs->insert_message(390, msg.c_str());

	int realScoreIndex = -1;

	for (int timestampIndex = 0; timestampIndex < request->get_n_points(); timestampIndex++)
	{
		AMResponse *response = responses->create_point_response(pid, request->get_timestamp(timestampIndex));
		for (int i = 0; i < this->verificationConfig.size(); i++)
		{
			response->verificationConfig.push_back(this->verificationConfig[i]);
		}

		// Set Scores
		int numOfScores = 0;
		response->need_to_update_scoreTypes = 1;
		string realScoreType = string(request->get_score_type(0));

		string setScorePrefix = string(SWITCH_SET_SCORE) + "-";
		string realScoreTypeSwitch = setScorePrefix + realScoreType;
		string realScoreTypeSwitch2 = "NO_SCORE";
		bool ramoveRealScoreType = false;
		int dataIndex = 0;
		for (dataIndex = 0; dataIndex < data.size(); dataIndex++) // Before real score
		{

			// When compiling to linux - use strcasecmp instead of stricmp
			if (strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch.c_str()) == 0 || strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch2.c_str()) == 0)
			{
				ramoveRealScoreType = true;
				break;
			}
		}

		for (dataIndex = 0; dataIndex < data.size(); dataIndex++) // Before real score
		{
			if (data[dataIndex]->signalName.find(SWITCH_SET_SCORE) == string::npos ||
				strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch.c_str()) == 0 ||
				strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch2.c_str()) == 0)
			{
				break;
			}
			string scoreType = string(data[dataIndex]->signalName);
			replace(scoreType, setScorePrefix, "");
			map<string, int>::iterator it = response->scoresIndexs.find(scoreType);
			if (it == response->scoresIndexs.end())
			{
				response->scoresIndexs.insert(pair<string, int>(scoreType, numOfScores++));
			}
		}
		while (dataIndex < data.size())
		{
			if (data[dataIndex]->signalName.find(SWITCH_SET_SCORE) != string::npos)
			{
				break;
			}
			dataIndex++;
		}
		if (!ramoveRealScoreType)
		{
			for (size_t scoreIndex = 0; scoreIndex < request->get_n_score_types(); scoreIndex++)
			{
				// if (!strcmp(request->get_score_type(scoreIndex), "Score") == 0 && !strcmp(request->get_score_type(scoreIndex), "TsScore") == 0)
				//	return LogData::EndFunction(AM_FAIL_RC, __func__);

				realScoreIndex = numOfScores;
				response->scoresIndexs.insert(pair<string, int>(request->get_score_type((int)scoreIndex), numOfScores++));
			}
		}
		for (; dataIndex < data.size(); dataIndex++) // Before real score
		{
			if (data[dataIndex]->signalName.find(SWITCH_SET_SCORE) == string::npos ||
				strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch.c_str()) == 0 ||
				strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch2.c_str()) == 0)
			{
				continue;
			}
			string scoreType = string(data[dataIndex]->signalName);
			replace(scoreType, setScorePrefix, "");
			map<string, int>::iterator it = response->scoresIndexs.find(scoreType);
			if (it == response->scoresIndexs.end())
			{
				response->scoresIndexs.insert(pair<string, int>(scoreType, numOfScores++));
			}
		}
		response->init_scores(numOfScores);

		float score = 0;
		long long ts = 0;
		AMMessages *messages = NULL;

		if (realScoreIndex != -1)
		{
			messages = response->get_score_msgs(realScoreIndex);
		}

		for (dataIndex = 0; dataIndex < data.size(); dataIndex++)
		{
			AMData *singleData = data[dataIndex];
			string signalName = singleData->signalName;

			AMMessages *responseMessages = response->get_msgs();
			if (singleData->signalName.find("SET_SCORE") != string::npos &&
				(strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch.c_str()) != 0 ||
				 strcasecmp(data[dataIndex]->signalName.c_str(), realScoreTypeSwitch2.c_str()) != 0))
			{
				string scoreType = string(data[dataIndex]->signalName);
				replace(scoreType, setScorePrefix, "");
				map<string, int>::iterator it = response->scoresIndexs.find(scoreType);
				if (it == response->scoresIndexs.end())
				{
					shared_msgs->insert_message(AM_GENERAL_FATAL, "Verification Error: Cannot find scoreType in map");
					return LogData::EndFunction(AM_FAIL_RC, __func__);
				}

				int scoreindex = it->second;
				AMMessages *messages = response->get_score_msgs(scoreindex);
				float fakeScore = 0;
				if (singleData->Values_len > 0)
				{
					fakeScore = singleData->Values[0];
				}

				if (singleData->Values_len > 1)
				{
					int msgCode = static_cast<int>(singleData->Values[1]);
					messages->insert_message(msgCode, "Fake Message Arg");
				}

				response->set_score(it->second, fakeScore, (char *)it->first.c_str());
			}
			else if (singleData->signalName.find(SWITCH_MSG) != string::npos || singleData->signalName.find(SWITCH_MSG2) != string::npos)
			{
				int msgCode = 999;
				if (singleData->Values_len > 0)
				{
					msgCode = static_cast<int>(singleData->Values[0]);
				}

				string msgArg = signalName;

				for (int index = 1; index < singleData->Values_len; index++)
				{
					stringstream stream;
					stream << fixed << setprecision(4) << singleData->Values[index];
					if (index == 1)
					{
						msgArg = stream.str();
					}
					else
					{
						msgArg = msgArg + "|" + stream.str();
					}
				}

				const char *msgCArg;
				if (msgArg.size() == 0)
				{
					msgCArg = NULL;
				}
				else
				{
					msgCArg = msgArg.c_str();
				}

				if (singleData->signalName.find(SWITCH_SMSG) != string::npos || singleData->signalName.find(SWITCH_SMSG2) != string::npos)
				{
					if (msgCArg == NULL)
					{
						shared_msgs->insert_message(msgCode, NULL);
					}
					else
					{
						shared_msgs->insert_message(msgCode, ("Shared:" + string(msgCArg)).c_str());
					}
				}
				else if (singleData->signalName.find(SWITCH_RMSG) != string::npos || singleData->signalName.find(SWITCH_RMSG2) != string::npos)
				{
					if (msgCArg == NULL)
					{
						responseMessages->insert_message(msgCode, NULL);
					}
					else
					{
						responseMessages->insert_message(msgCode, ("Response:" + string(msgCArg)).c_str());
					}
				}
				else if (realScoreIndex != -1)
				{
					messages->insert_message(msgCode, msgCArg);
				}
			}
			else if (realScoreIndex != -1)
			{
				for (size_t index = 0; index < singleData->Values_len; index++)
				{
					score = score + (*singleData).Values[index];
					if (strcmp(get_name(), "SignalsSum2") == 0 || strcmp(get_name(), "XSignalsSum2") == 0)
					{
						score = score + 1.0f;
					}
				}

				for (size_t index = 0; index < singleData->TimeStamps_len; index++)
				{
					ts = ts + (*singleData).TimeStamps[index] - 20000000;
					if (strcmp(get_name(), "SignalsSum2") == 0 || strcmp(get_name(), "XSignalsSum2") == 0)
					{
						ts = ts + 1;
					}
				}
			}
		}

		if (realScoreIndex != -1)
		{
			for (map<string, int>::iterator iter = response->scoresIndexs.begin(); iter != response->scoresIndexs.end(); ++iter)
			{
				float realScore = score;
				if (strcmp(iter->first.c_str(), "TsScore") == 0)
				{
					realScore = static_cast<float>(ts);
				}
				else if (strcmp(iter->first.c_str(), "S_SUM2") == 0)
				{
					realScore += data.size();
				}
				response->set_score(iter->second, realScore, (char *)iter->first.c_str());
			}
		}
	}

	return AM_OK_RC;
}

template <typename Out>
void split(const std::string &s, const char delim, Out result)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
	{
		*(result++) = item;
	}
}

//-----------------------------------------------------------------------------------
// private internals for class MedialInfraAlgoMarker
//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------

//===========================================================================================================
//===========================================================================================================
//===========================================================================================================
// A P I   I M P L E M E N T A T I O N S
//===========================================================================================================
//===========================================================================================================
//===========================================================================================================

//-----------------------------------------------------------------------------------------------------------
// create a new AlgoMarker of type am_type and init its name
//-----------------------------------------------------------------------------------------------------------
int AM_API_Create(int am_type, AlgoMarker **new_am)
{
	LogData::AddArgument("am_type", am_type);
	LogData::AddArgument("new_am", new_am);
	LogData::StartFunction(__func__);

	try
	{
		*new_am = AlgoMarker::make_algomarker((AlgoMarkerType)am_type);

		if (new_am == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		LogData::AddArgument("am_type", am_type);
		LogData::AddArgument("new_am", new_am);
		return LogData::EndFunction(AM_OK_RC, __func__);
		return AM_OK_RC;
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// loading AlgoMarker and making it ready to get Requests
//-----------------------------------------------------------------------------------------------------------
int AM_API_Load(AlgoMarker *pAlgoMarker, const char *config_fname)
{
	vector<string> elems;
	split(config_fname, '|', std::back_inserter(elems));
	pAlgoMarker->set_name(elems[0].c_str());

	for (auto index = 1; index < 6; index++)
	{
		if (index < elems.size())
		{
			pAlgoMarker->verificationConfig.push_back(elems[index]);
		}
		else
		{
			pAlgoMarker->verificationConfig.push_back("");
		}
	}

	size_t start_pos = string(pAlgoMarker->get_name()).find("_LOG");
	if (start_pos != std::string::npos)
	{
		LogData::TurnOnLogs = true;
	}
	if (expandEnvVars(AM_LOG_CONTROL_ENV_VAR) == "1")
		LogData::TurnOnLogs = true;
	if (expandEnvVars(AM_LOG_CONTROL_ENV_VAR) == "0")
		LogData::TurnOnLogs = false;
	// printf("AM_LOG_CONTROL_ENV_VAR = %s\n", expandEnvVars(AM_LOG_CONTROL_ENV_VAR).c_str());
	// printf("LogData::TurnOnLogs = %d\n", (int)LogData::TurnOnLogs);

	LogData::AddArgument("pAlgoMarker", pAlgoMarker);
	LogData::AddArgument("config_fname", config_fname);
	LogData::StartFunction(__func__);

	try
	{
		LogData::AddArgument("config_fname", config_fname);
		LogData::AddArgument("pAlgoMarker", pAlgoMarker);

		if (pAlgoMarker == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		return SpecialReturn(pAlgoMarker, __func__);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// clearing data from AlgoMarker (recommended at the start and/or end of each query session
//-----------------------------------------------------------------------------------------------------------
int AM_API_ClearData(AlgoMarker *pAlgoMarker)
{
	LogData::AddArgument("pAlgoMarker", pAlgoMarker);
	LogData::StartFunction(__func__);

	try
	{
		LogData::AddArgument("pAlgoMarker", pAlgoMarker);
		return SpecialReturn(pAlgoMarker, __func__, pAlgoMarker->ClearData());
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// adding data to an AlgoMarker
// this API allows adding a specific signal, with matching arrays of times and values
//-----------------------------------------------------------------------------------------------------------
int AM_API_AddData(AlgoMarker *pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, float *Values)
{
	LogData::AddArgument("patient_id", patient_id);
	LogData::AddArgument("signalName", signalName);
	LogData::AddArgument("TimeStamps_len", TimeStamps_len);
	LogData::AddArgument("TimeStamps", TimeStamps, TimeStamps_len);
	LogData::AddArgument("Values_len", Values_len);
	LogData::AddArgument("Values", Values, Values_len);
	LogData::StartFunction(__func__);

	try
	{
		if (pAlgoMarker == NULL)
		{
			return LogData::EndFunction(AM_FAIL_RC, __func__);
		}

		auto result = pAlgoMarker->AddData(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
		LogData::AddArgument("patient_id", patient_id);
		LogData::AddArgument("signalName", signalName);
		LogData::AddArgument("TimeStamps_len", TimeStamps_len);
		LogData::AddArgument("TimeStamps", TimeStamps, TimeStamps_len);
		LogData::AddArgument("Values_len", Values_len);
		LogData::AddArgument("Values", Values, Values_len);
		return SpecialReturn(pAlgoMarker, __func__, result);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// adding string data to an AlgoMarker
// this API allows adding a specific signal, with matching arrays of times and string values
//-----------------------------------------------------------------------------------------------------------
int AM_API_AddDataStr(AlgoMarker *pAlgoMarker, int patient_id, const char *signalName, int TimeStamps_len, long long *TimeStamps, int Values_len, char **Values)
{
	LogData::AddArgument("patient_id", patient_id);
	LogData::AddArgument("signalName", signalName);
	LogData::AddArgument("TimeStamps_len", TimeStamps_len);
	LogData::AddArgument("TimeStamps", TimeStamps, TimeStamps_len);
	LogData::AddArgument("Values_len", Values_len);
	LogData::AddArgument("Values", Values, Values_len);
	LogData::StartFunction(__func__);

	try
	{
		if (pAlgoMarker == NULL)
		{
			return LogData::EndFunction(AM_FAIL_RC, __func__);
		}

		auto result = pAlgoMarker->AddDataStr(patient_id, signalName, TimeStamps_len, TimeStamps, Values_len, Values);
		LogData::AddArgument("patient_id", patient_id);
		LogData::AddArgument("signalName", signalName);
		LogData::AddArgument("TimeStamps_len", TimeStamps_len);
		LogData::AddArgument("TimeStamps", TimeStamps, TimeStamps_len);
		LogData::AddArgument("Values_len", Values_len);
		LogData::AddArgument("Values", Values, Values_len);
		return SpecialReturn(pAlgoMarker, __func__, result);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Prepare a Request
// Null RC means failure
// pids and timestamps here are the timepoints to give predictions at
//-----------------------------------------------------------------------------------------------------------
int AM_API_CreateRequest(char *requestId, char **_score_types, int n_score_types, int *patient_ids, long long *time_stamps, int n_points, AMRequest **new_req)
{
	LogData::AddArgument("requestId", requestId);
	LogData::AddArgument("_score_types", _score_types, n_score_types);
	LogData::AddArgument("patient_ids", patient_ids, 1);
	LogData::AddArgument("time_stamps", time_stamps, n_points);
	LogData::AddArgument("new_req", new_req);
	LogData::StartFunction(__func__);

	try
	{
		(*new_req) = new AMRequest;

		if ((*new_req) == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		(*new_req)->set_request_id(requestId);
		(*new_req)->insert_score_types(_score_types, n_score_types);
		for (int i = 0; i < n_points; i++)
			(*new_req)->insert_point(patient_ids[i], time_stamps[i]);

		LogData::AddArgument("requestId", requestId);
		LogData::AddArgument("_score_types", _score_types, n_score_types);
		LogData::AddArgument("patient_ids", patient_ids, 1);
		LogData::AddArgument("time_stamps", time_stamps, n_points);
		LogData::AddArgument("new_req", new_req);
		return LogData::EndFunction(AM_OK_RC, __func__);
		;
	}
	catch (...)
	{
		(*new_req) = NULL;
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Get scores for a ready request
//-----------------------------------------------------------------------------------------------------------
int AM_API_Calculate(AlgoMarker *pAlgoMarker, AMRequest *request, AMResponses *responses)
{
	LogData::AddArgument("pAlgoMarker", pAlgoMarker);
	LogData::AddArgument("request", request);
	LogData::AddArgument("responses", responses);
	LogData::StartFunction(__func__);

	try
	{
		if (pAlgoMarker == NULL || request == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		const int result = pAlgoMarker->Calculate(request, responses);
		// ChangeParameter(pAlgoMarker->verificationConfig, __func__, &responses);
		LogData::AddArgument("pAlgoMarker", pAlgoMarker);
		LogData::AddArgument("request", request);
		LogData::AddArgument("responses", responses);
		return SpecialReturn(responses, __func__, result);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Create a new empty responses object to be later used
//-----------------------------------------------------------------------------------------------------------
int AM_API_CreateResponses(AMResponses **new_responses)
{
	LogData::AddArgument("new_responses", new_responses);
	LogData::StartFunction(__func__);

	try
	{
		(*new_responses) = new AMResponses;

		if ((*new_responses) == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		LogData::AddArgument("new_responses", new_responses);
		return LogData::EndFunction(AM_OK_RC, __func__);
	}
	catch (...)
	{
		(*new_responses) = NULL;
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of AlgoMarker - free all memory
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeAlgoMarker(AlgoMarker *pAlgoMarker)
{
	LogData::AddArgument("pAlgoMarker", pAlgoMarker);
	LogData::StartFunction(__func__);

	try
	{
		if (pAlgoMarker == NULL)
			return;

		pAlgoMarker->Unload();
		auto p = (MedialInfraAlgoMarker *)pAlgoMarker;

		delete p;
	}
	catch (...)
	{
	}

	LogData::AddArgument("pAlgoMarker", pAlgoMarker);
	LogData::EndFunction(__func__);

	LogData::closeLogFile();
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of AMRequest - free all memory
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeRequest(AMRequest *pRequest)
{
	LogData::AddArgument("pRequest", pRequest);
	LogData::StartFunction(__func__);

	try
	{
		if (pRequest == NULL)
			return;
		pRequest->clear();
		delete pRequest;
	}
	catch (...)
	{
	}

	LogData::AddArgument("pRequest", pRequest);
	LogData::EndFunction(__func__);

	if (LogData::TurnOnLogs)
	{
		LogData::WriteToLog("\r\n--------------- :) ##### (: ##### ;) ##### (; ##### :-) ##### (-: ##### ;-) ##### (+: ##### :) ##### (-: ##### :) ##### (-; ##### 8) ##### (: ---------------\r\n\r\n");
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Dispose of responses - free all memory
//-----------------------------------------------------------------------------------------------------------
void AM_API_DisposeResponses(AMResponses *responses)
{
	LogData::AddArgument("responses", responses);
	LogData::StartFunction(__func__);

	try
	{
		if (responses == NULL)
			return;
		delete responses;
	}
	catch (...)
	{
	}

	LogData::AddArgument("responses", responses);
	LogData::EndFunction(__func__);
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get number of responses (= no. of pid,time result points)
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponsesNum(AMResponses *responses)
{
	LogData::AddArgument("responses", responses);
	LogData::StartFunction(__func__);

	try
	{
		if (responses == NULL)
			return LogData::EndFunction(0, __func__);

		const int result = responses->get_n_responses();

		LogData::AddArgument("responses", responses);
		return SpecialReturn(responses, __func__, result);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get shared msgs. Not a copy - direct pointers, so do not free.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetSharedMessages(AMResponses *responses, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	LogData::AddArgument("n_msgs", n_msgs);
	LogData::AddArgument("msgs_codes", msgs_codes, 0);
	LogData::AddArgument("msgs_args", msgs_args, 0);
	LogData::StartFunction(__func__);

	try
	{
		if (responses == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		AMMessages *shared_m = responses->get_shared_messages();
		shared_m->get_messages(n_msgs, msgs_codes, msgs_args);

		ChangeParameter(responses->verificationConfig, __func__, msgs_args);
		LogData::AddArgument("n_msgs", n_msgs);
		LogData::AddArgument("msgs_codes", msgs_codes, *n_msgs);
		LogData::AddArgument("msgs_args", msgs_args, *n_msgs);
		return SpecialReturn(responses, __func__);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get an index of a specific pid,time response, or -1 if it doesn't exist
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseIndex(AMResponses *responses, int _pid, long long _timestamp)
{
	LogData::AddArgument("responses", responses);
	LogData::AddArgument("_pid", _pid);
	LogData::AddArgument("_timestamp", _timestamp);
	LogData::StartFunction(__func__);

	try
	{
		const int result = responses->get_response_index_by_point(_pid, _timestamp);

		LogData::AddArgument("responses", responses);
		LogData::AddArgument("_pid", _pid);
		LogData::AddArgument("_timestamp", _timestamp);
		return SpecialReturn(responses, __func__, result);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get scores for a scpefic response given its index.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseAtIndex(AMResponses *responses, int res_index, AMResponse **res)
{
	LogData::AddArgument("responses", responses);
	LogData::AddArgument("res_index", res_index);
	LogData::AddArgument("res", res);
	LogData::StartFunction(__func__);

	try
	{
		*res = NULL;
		if (responses == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		if (res_index < 0 || res_index >= responses->get_n_responses())
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		*res = responses->get_response(res_index);

		LogData::AddArgument("responses", responses);
		LogData::AddArgument("res_index", res_index);
		LogData::AddArgument("res", res);
		return SpecialReturn(responses, __func__);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get number of scores in a response (could contain several score types)
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoresNum(AMResponse *response, int *n_scores)
{
	LogData::AddArgument("response", response);
	LogData::AddArgument("n_scores", n_scores);
	LogData::StartFunction(__func__);

	try
	{
		if (response == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		*n_scores = response->get_n_scores();

		ChangeParameter(response->verificationConfig, __func__, n_scores);
		LogData::AddArgument("response", response);
		LogData::AddArgument("n_scores", n_scores);
		return SpecialReturn(response, __func__, AM_OK_RC);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// given a score index , return all we need about it : pid , timestamp, score and score type
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseScoreByIndex(AMResponse *response, int score_index, float *_score, char **_score_type)
{
	LogData::AddArgument("response", response);
	LogData::AddArgument("score_index", score_index);
	LogData::AddArgument("_score", _score);
	LogData::AddArgument("_score_type", "[]");
	LogData::StartFunction(__func__);

	try
	{
		if (response == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		if (score_index < 0 || score_index >= response->get_n_scores())
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		/**pid = response->get_patient_id();
		 *timestamp = response->get_timestamp();*/
		if (response->get_score(score_index, _score, _score_type) != AM_OK_RC)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		ChangeParameter(response->verificationConfig, __func__, _score_type);
		LogData::AddArgument("response", response);
		LogData::AddArgument("score_index", score_index);
		LogData::AddArgument("_score", _score);
		LogData::AddArgument("_score_type", _score_type);
		return SpecialReturn(response, __func__);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get all messages for a specific response given its index
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponseMessages(AMResponse *response, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	LogData::AddArgument("response", response);
	LogData::AddArgument("n_msgs", n_msgs);
	LogData::AddArgument("msgs_codes", msgs_codes, *n_msgs);
	LogData::AddArgument("msgs_args", msgs_args, *n_msgs);
	LogData::StartFunction(__func__);

	try
	{
		if (response == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		response->get_msgs()->get_messages(n_msgs, msgs_codes, msgs_args);
		LogData::AddArgument("response", response);
		LogData::AddArgument("n_msgs", n_msgs);
		LogData::AddArgument("msgs_codes", msgs_codes, *n_msgs);
		LogData::AddArgument("msgs_args", msgs_args, *n_msgs);
		return SpecialReturn(response, __func__);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get all messages for a specific response given its index
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetScoreMessages(AMResponse *response, int score_index, int *n_msgs, int **msgs_codes, char ***msgs_args)
{
	LogData::AddArgument("response", response);
	LogData::AddArgument("score_index", score_index);
	LogData::AddArgument("n_msgs", n_msgs);
	LogData::AddArgument("msgs_codes", msgs_codes, 0);
	LogData::AddArgument("msgs_args", msgs_args, 0);
	LogData::StartFunction(__func__);

	try
	{
		if (response == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		if (score_index < 0 || score_index >= response->get_n_scores())
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		response->get_score_msgs(score_index)->get_messages(n_msgs, msgs_codes, msgs_args);
		ChangeParameter(response->verificationConfig, __func__, msgs_args);
		LogData::AddArgument("response", response);
		LogData::AddArgument("score_index", score_index);
		LogData::AddArgument("n_msgs", n_msgs);
		LogData::AddArgument("msgs_codes", msgs_codes, n_msgs[0]);
		LogData::AddArgument("msgs_args", msgs_args, n_msgs[0]);
		return SpecialReturn(response, __func__);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get pid and timestamp of a response
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponsePoint(AMResponse *response, int *pid, long long *timestamp)
{
	LogData::AddArgument("response", response);
	LogData::AddArgument("pid", pid);
	LogData::AddArgument("timestamp", timestamp);
	LogData::StartFunction(__func__);

	try
	{
		if (response == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		*pid = response->get_patient_id();
		*timestamp = response->get_timestamp();
		LogData::AddArgument("response", response);
		LogData::AddArgument("pid", pid);
		LogData::AddArgument("timestamp", timestamp);
		return SpecialReturn(response, __func__);
		;
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get request id . Direct pointer so do not free.
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetResponsesRequestId(AMResponses *responses, char **requestId)
{
	LogData::AddArgument("responses", responses);
	LogData::AddArgument("requestId", "");
	LogData::StartFunction(__func__);

	try
	{
		if (responses == NULL)
			return LogData::EndFunction(AM_FAIL_RC, __func__);

		*requestId = responses->get_request_id();
		ChangeParameter(responses->verificationConfig, __func__, requestId);

		LogData::AddArgument("responses", responses);
		LogData::AddArgument("requestId", requestId);
		return SpecialReturn(responses, __func__);
	}
	catch (...)
	{
		LogData::AddArgument("responses", responses);
		LogData::AddArgument("requestId", requestId);
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// get the nameof an algo marker
//-----------------------------------------------------------------------------------------------------------
int AM_API_GetName(AlgoMarker *pAlgoMarker, char **name)
{
	LogData::AddArgument("pAlgoMarker", pAlgoMarker);
	LogData::AddArgument("name", "[]");
	LogData::StartFunction(__func__);
	try
	{
		*name = NULL;
		if (pAlgoMarker == NULL)
		{
			LogData::AddArgument("name", name);
			return LogData::EndFunction(AM_FAIL_RC, __func__);
		}

		*name = pAlgoMarker->get_name();
		ChangeParameter(pAlgoMarker->verificationConfig, __func__, name);
		LogData::AddArgument("pAlgoMarker", pAlgoMarker);
		LogData::AddArgument("name", name);
		return SpecialReturn(pAlgoMarker, __func__);
	}
	catch (...)
	{
		return LogData::EndFunction(AM_FAIL_RC, __func__);
	}
}
//-----------------------------------------------------------------------------------------------------------

int AM_API_GetResponseScoreByType(AMResponses *responses, int res_index, char *_score_type, float *out_score)
{
	throw runtime_error("AM_API_GetResponseScoreByType Not implemented in this AM");
	return AM_FAIL_RC;
}
