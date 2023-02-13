#include "glibc_extension.h"

int feenableexcept(int e) throw () {
	return 0;
	//gcc -Wall -Werror -shared -static -static-libstdc++ -o libalpine.so lib_ext.cpp
}