#ifndef UTIL_H_

#define UTIL_H_

double second();

extern double elapsed_time;

#define T1 {elapsed_time=second();}
#define T2(a) { \
  elapsed_time=second()-elapsed_time;  \
  printf("Step %s: elapsed time: %f secs.\n",(a),elapsed_time);  \
  fflush(stdout); \
}

#endif /* UTIL_H_ */
