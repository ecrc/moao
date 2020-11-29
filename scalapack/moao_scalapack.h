#include "myscalapack.h"
#include "moao_defs.h"

#define HANDLE_ERROR(_result, _func)               \
  if(0 != _result) {                  \
    fprintf(stderr,"An error occured (%d), when calling %s at line: %d in file %s \n\n", _result, _func, __LINE__ -1, __FILE__);   \
    return _result;                                         \
  }
  
int reconstructor(int nmeas, int nmeasts, real_t *Cmm, int *descCmm, real_t *Ctm, int *descCtm, long maxgal, long offset);
int compute_Cee_Cvv(int nmeas, int nmeasts, int nact, real_t *Cmm, int *descCmm, real_t *Cpp, int *descCpp, real_t *Cpm, int *descCpm, real_t *R, int *descR, real_t *Dx, int *descDx, real_t *Cee, int *descCee, real_t *Cvv, int *descCvv, real_t *Tmp, int *descTmp, double* flops);
