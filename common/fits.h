/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef FITS_H
#define FITS_H

#include "moao_defs.h"

#ifdef __cplusplus
extern "C"{
#endif
int fitsExists(char *fileName);
void concatenateFileName(char *fileName, char *filePath, char *typeMat, char *suffix);
void getFitsDims(char *fileName,long *naxes);
void readFits(char *fileName,double *data);
void writeFits(char *file_name, int nx, int ny, real_t* data);

#ifdef __cplusplus
}
#endif


#endif //FITS_H
