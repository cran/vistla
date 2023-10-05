#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_thread_num() 0
 #define omp_get_max_threads() 1
 #define omp_set_num_threads(x)
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h> 

#include <stdint.h>
#include <stdbool.h>

typedef uint32_t u32;

//Hash table

#include "ht.h"

//Heap (priority queue)

#include "heap.h"

//Common stuff

#include "shared.h"

//Input conversion

#include "convert.h"

//Algorithm

#include "vistla.h"

//Registration

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_CallMethodDef R_CallDef[]={
 CALLDEF(C_heapTest,3),
 CALLDEF(C_heapTiedTest,2),
 CALLDEF(C_convertTest,2),
 CALLDEF(C_vistla,8),
 {NULL,NULL,0}
};

void attribute_visible R_init_vistla(DllInfo *dll){
 R_registerRoutines(dll,NULL,R_CallDef,NULL,NULL);
 R_useDynamicSymbols(dll,FALSE);
 R_forceSymbols(dll,TRUE);
}
