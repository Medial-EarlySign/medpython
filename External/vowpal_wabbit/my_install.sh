#!/bin/bash
./configure --with-boost-libdir=/server/Linux/alon/packages/boost_1_63_0/stage/lib
#change cinfig.status BOOST_CPPFLAGS to /server/Linux/alon/packages/boost_1_63_0
make -j 8