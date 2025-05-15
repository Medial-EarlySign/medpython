#pragma once

#define AM_UNDEFINED_VALUE -9999.99

//=====================================================================
// External Error Codes
//=====================================================================
#define AM_ELIGIBILITY_ERROR						300
#define AM_RESPONSES_ELIGIBILITY_ERROR				301
#define AM_DATA_BAD_FORMAT_FATAL					330
#define AM_DATA_UNKNOWN_SIGNAL						331
#define AM_DATA_BAD_STRUCTURE						332
#define AM_DATA_GENERAL_ERROR						333
#define AM_DATA_BAD_FORMAT_NON_FATAL				334

#define AM_THRESHOLD_ERROR_NON_FATAL				392 //350
#define AM_GENERAL_FATAL							392

//=====================================================================
// Internal Error codes used for AlgoMarkers 
//=====================================================================

// OK RESPONSE
#define AM_OK_RC									0

// General FAIL RC
#define AM_FAIL_RC									-1

// SPECIFIC ERR CODES

// Create() Errors
#define AM_ERROR_CREATE_FAILED						1001

// Load() Errors
#define AM_ERROR_LOAD_NO_CONFIG_FILE				1101
#define AM_ERROR_LOAD_BAD_CONFIG_FILE				1102
#define AM_ERROR_LOAD_NON_MATCHING_TYPE				1103
#define AM_ERROR_LOAD_READ_REP_ERR					1104
#define AM_ERROR_LOAD_READ_MODEL_ERR				1105
#define AM_ERROR_LOAD_BAD_NAME						1106
#define AM_ERROR_LOAD_BAD_TESTERS_FILE				1107
#define AM_ERROR_LOAD_MISSING_REQ_SIGS				1108
#define AM_ERROR_MUST_BE_LOADED						1109
#define AM_ERROR_UNKNOWN_LOAD_TYPE					1110
#define AM_ERROR_READING_DICT_FILE					1111
#define AM_ERROR_PARSING_JSON_DICT					1112

// AddData Errors
#define AM_ERROR_ADD_DATA_FAILED					1201 //Only in old API - AddData, not AddDataByType
#define AM_ERROR_DATA_JSON_PARSE					1202
#define AM_ERROR_DATA_UNKNOWN_ADD_DATA_TYPE			1203 //Can't happen
// Messages Codes For Responses - in old API Calculate
#define AM_MSG_NULL_REQUEST							101
#define AM_MSG_BAD_PREDICTION_POINT					102
#define AM_MSG_BAD_SCORE_TYPE						103
#define AM_MSG_RAW_SCORES_ERROR						104



