//
// General utils
//

#ifndef __UTILS_H__
#define __UTILS_H__

#include <cstdlib>
#include <vector>
#include <random>
using namespace std;

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

//#define FRAND30 ((double)((unsigned int)(my_rand()<<15)+(unsigned int)(my_rand()))/(double)(1<<30))
//#define FRAND30 ((double)((unsigned int)(rand()<<15)+(unsigned int)(rand()))/(double)(1<<30))
//#define FRAND30_R(s) ((double)((unsigned int)(rand_r(s)<<15)+(unsigned int)(rand_r(s)))/(double)(1<<30))
//#define IRAND(S) ((int)((double)(S) * FRAND30))
//#define IRAND_R(S,s) ((int)((double)(S) * FRAND30_R(s)))

#define SQUARE(x) ((x)*(x))
#define POWER_2(x) ((x)*(x))
#define POWER_3(x) ((x)*(x)*(x))

double Get_AUC(vector<float> &scores, vector<char> &Y);
double Get_AUC(vector<float> &scores, vector<char> &Y, vector<double> &spe, vector<double> &sen);

class QRFglobalRNG
{
public:
	static std::minstd_rand::result_type rand() { return getInstance()._rng(); };
	static unsigned int rand30() { return ((getInstance()._rng() << 15) ^ getInstance()._rng()) & 0x4fffffff; };
	static void srand(std::minstd_rand::result_type val) { getInstance()._rng.seed(val); };
	static const std::minstd_rand::result_type max() { return getInstance()._rng.max(); };

private:
	std::minstd_rand _rng;

	static QRFglobalRNG& getInstance()
	{
		static QRFglobalRNG instance; // instansiated on first call
		return instance;
	}

	QRFglobalRNG() : _rng(20150715) {}; // constructor

	QRFglobalRNG(QRFglobalRNG const&) { fprintf(stderr, "Error: copying is forbidden for the globalRNG object\n"); exit(-1); }; // no copy constructor
	void operator=(QRFglobalRNG const&) { fprintf(stderr, "Error: assignment is forbidden for the globalRNG object\n"); exit(-1); }; // no assignment operator

};

#endif