#pragma once
//=====================================================================
// External Error Codes
//=====================================================================
#define AM_ELIGIBILITY_ERROR						300
#define AM_INSUFFICIENT_DATA						310
#define AM_INVALID_DATA								311
#define AM_OUTLIER_NON_FATAL						320
#define AM_OUTLIER_FATAL							321
#define AM_GENERAL_INFO								390
#define AM_GENERAL_NON_FATAL						391
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

// AddData Errors
#define AM_ERROR_ADD_DATA_FAILED					1201

// Messages Codes For Responses
#define AM_MSG_NULL_REQUEST							101
#define AM_MSG_BAD_PREDICTION_POINT					102
#define AM_MSG_BAD_SCORE_TYPE						103
#define AM_MSG_RAW_SCORES_ERROR						104



