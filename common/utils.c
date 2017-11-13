/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "utils.h"
#include <ctype.h>
#include <stdlib.h>
#include <sys/time.h>

#define ARRAYSIZE 1024  //line length can be quite large if there is lots of layers

void read_paraml(FILE *file, long *par) {
  char line[ARRAYSIZE];
  fgets(line, sizeof(line), file);
  fgets(line, sizeof(line), file);
  long tmpi;
  if( sscanf(line, "%ld", &tmpi) == 1 ) *par = tmpi;
}

void read_parami(FILE *file, int *par) {
  char line[ARRAYSIZE];
  fgets(line, sizeof(line), file);
  fgets(line, sizeof(line), file);
  int tmpi;
  if( sscanf(line, "%d", &tmpi) == 1 ) *par = tmpi;
}

void read_paramd(FILE *file, real_t *par) {
  char line[ARRAYSIZE];
  fgets(line, sizeof(line), file);
  fgets(line, sizeof(line), file);
//TODO should we read in double precicion if needed, instead?
  float tmpf;
  if( sscanf(line, "%f", &tmpf) == 1 ) *par = tmpf;
}

void read_arrayl(FILE *file, long *arr, int nmax) {
  char line[ARRAYSIZE];
  char *tail = line;
  int i = 0;
  fgets(line, sizeof(line), file);
  fgets(line, sizeof(line), file);
  while((*tail) && (i < nmax)){
    while (isspace (*tail)) tail++;
    arr[i++] = strtol(tail,&tail,10);
  }
}

void read_arrayi(FILE *file, int *arr, int nmax) {
  char line[ARRAYSIZE];
  char *tail = line;
  int i = 0;
  fgets(line, sizeof(line), file);
  fgets(line, sizeof(line), file);
  while((*tail) && (i < nmax)){
    while (isspace (*tail)) tail++;
    arr[i++] = (int)strtol(tail,&tail,10);
  }
}

void read_arrayd(FILE *file, real_t *arr, int nmax) {
  char line[ARRAYSIZE];
  char *tail = line;
  int i = 0;
  fgets(line, sizeof(line), file);
  fgets(line, sizeof(line), file);
  while((*tail) && (i < nmax)){
    while (isspace (*tail)) tail++;
    arr[i++] = strtod(tail,&tail);
  }
}


double cWtime(void)
{
  struct timeval tp;
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec;
}
