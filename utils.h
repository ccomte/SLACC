#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <unistd.h>   /* for getopt */
#include <ctype.h>    /* for isprint */
#include <stdint.h>   /* unsigned integer type uint64_t */



static inline double uniform () {
  return((double) (random() % (RAND_MAX-2) + 1) / RAND_MAX);
}

static inline int bernoulli (double p) {
  return (uniform() < p);
}

static inline double exponential (double lambda) {
  return -1. / lambda * log(uniform());
}



#endif
