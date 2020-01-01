#!/bin/bash
CXX=/usr/local/bin/c++; CC=/usr/local/bin/gcc ./configure --with-boost-libdir=/server/Linux/alon/packages/boost_1_63_0/stage/lib --enable-parallelization
#change cinfig.status BOOST_CPPFLAGS to /server/Linux/alon/packages/boost_1_63_0
make -j 8