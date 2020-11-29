/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef FITS_H
#define FITS_H


int fitsExists(const char *fileName);
void concatenateFileName(char *fileName, char *filePath, char *typeMat, char *suffix);
void getFitsDims(const char *fileName,long *naxes);
template<typename T>
void readFits(const char *fileName,T *data);
template<typename T>
void writeFits(const char *file_name, int nx, int ny, T* data);



#endif //FITS_H
