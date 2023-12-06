#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "ks.h"

/* .C calls */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
   CALLDEF(pSmirnov2x, 3),
   CALLDEF(pKS2, 2),
   {NULL, NULL, 0}
};

void attribute_visible R_init_twang(DllInfo *dll)
{
   R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
   R_forceSymbols(dll, TRUE);
}
