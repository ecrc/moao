/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef MATCOV_H
#define MATCOV_H

#include "tomo_struct/tomo_struct.hpp"

namespace MOAO_COMMON{
template<typename T>
void matcov_tile(
  T* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  Tomo_struct<T> *tomo, int part);
}

#endif //MATCOV_H
