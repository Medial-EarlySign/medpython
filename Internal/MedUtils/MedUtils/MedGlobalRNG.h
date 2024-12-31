#ifndef __GLOBAL_RNG_H__
#define __GLOBAL_RNG_H__

#include <random>

class globalRNG
{
public:
  static std::minstd_rand::result_type rand() { return getInstance()._rng(); };
  static unsigned int rand30() { return ((getInstance()._rng() << 15) ^ getInstance()._rng()) & 0x4fffffff; };
  static void srand(std::minstd_rand::result_type val) { getInstance()._rng.seed(val); };
  static const std::minstd_rand::result_type max() { return getInstance()._rng.max(); };
  static std::mt19937 &get_engine() { return getInstance()._rng_new; };

private:
  std::minstd_rand _rng;
  std::mt19937 _rng_new;

  static globalRNG &getInstance()
  {
    static globalRNG instance; // instansiated on first call
    return instance;
  }

  globalRNG() : _rng(20150715), _rng_new(20150715) {}; // constructor

  globalRNG(globalRNG const &)
  {
    fprintf(stderr, "Error: copying is forbidden for the globalRNG object\n");
    exit(-1);
  }; // no copy constructor
  void operator=(globalRNG const &)
  {
    fprintf(stderr, "Error: assignment is forbidden for the globalRNG object\n");
    exit(-1);
  }; // no assignment operator
};

// inline float rand_1();
// inline int rand_N(int N);
// inline float rand_range(float from, float to);

// float random number in 0...1
inline float rand_1() { return (float)globalRNG::rand() / (float)(globalRNG::max() + 1.0); }

// int random number in 0..N-1
inline int rand_N(int N) { return ((int)((float)N * rand_1() - (float)1e-8)); }
inline size_t rand_N_i64(size_t N) { return ((size_t)((double)N * rand_1() - (double)1e-20)); }

// float number in from...to range
inline float rand_range(float from, float to) { return from + rand_1() * (to - from); }

#endif
