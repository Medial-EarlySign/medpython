#include "AlgoMarkerInternal.h"

//-------------------------------------------------------------------------------------------------------------------------
InputTester *InputTester::make_input_tester(InputTesterType it_type)
{
	if (it_type == INPUT_TESTER_TYPE_SIMPLE)
		return new InputTesterSimple;

	return NULL;
}

//-------------------------------------------------------------------------------------------------------------------------
void InputTesterSimple::input_from_string(const string &in_str)
{
	sf.init_from_string(in_str);
}

//-------------------------------------------------------------------------------------------------------------------------
// 1: good to go 0: did not pass -1: could not test
int InputTesterSimple::test_if_ok(MedRepository &rep, int pid, long long timestamp, int &nvals, int &noutliers)
{
	MedSample s;

	s.id = pid;
	s.time = (int)timestamp;

	return sf.test_filter(s, rep, nvals, noutliers);
}
