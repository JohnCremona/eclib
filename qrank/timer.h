#ifdef LiDIA_INTS
#include "LiDIA/timer.h"
#else
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#endif

void init_time();
void start_time();
void stop_time();
void show_time();

#ifndef TIME_CONICS
#define TIME_CONICS 0
#endif

void init_conic_timer();
void start_conic_timer();
void stop_conic_timer();
void show_conic_timer();
