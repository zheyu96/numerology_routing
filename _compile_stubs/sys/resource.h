// Windows syntax-check stub for <sys/resource.h> (real header exists only on
// POSIX).  Only the pieces main.cpp touches are declared.
#pragma once

#define RUSAGE_SELF 0

struct rusage {
    long ru_maxrss;
};

inline int getrusage(int, struct rusage* ru) {
    if(ru) ru->ru_maxrss = 0;
    return -1;
}
