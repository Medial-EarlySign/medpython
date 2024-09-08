#ifndef __GLOBAL_RNG_H__
#define __GLOBAL_RNG_H__ 
#pragma once

#include <random>

class globalRNG
{
 public:
  static std::minstd_rand::result_type rand() {return getInstance()._rng();};
  static unsigned int rand30() {return ((getInstance()._rng() << 15) ^ getInstance()._rng()) & 0x3fffffff;};
  static void srand(std::minstd_rand::result_type val) {getInstance()._rng.seed(val);};
  static const std::minstd_rand::result_type max() {return getInstance()._rng.max();};

 private:
  std::minstd_rand _rng;

  static globalRNG& getInstance()
  {
    static globalRNG instance; // instansiated on first call
    return instance;
  }

  globalRNG() : _rng(20150715) {}; // constructor

  globalRNG(globalRNG const&) {fprintf(stderr, "Error: copying is forbidden for the globalRNG object\n"); exit(-1);}; // no copy constructor
  void operator=(globalRNG const&) {fprintf(stderr, "Error: assignment is forbidden for the globalRNG object\n"); exit(-1);}; // no assignment operator

};

#endif
