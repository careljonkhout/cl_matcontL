#include "fem.h"
/* --------------------------------------------------- */
/* Automatically generated by mwrap                    */
/* --------------------------------------------------- */

#include <mex.h>
#include <stdio.h>
#include <string.h>


#ifndef ulong
#  define ulong unsigned long
#endif
#ifndef uint
#  define uint  unsigned int
#endif
#ifndef uchar
#  define uchar unsigned char
#endif


/*
 * Support routines for copying data into and out of the MEX stubs
 */

void* mxWrapGetP(const mxArray* a, const char* fmt)
{
    void* p = 0;
    if (mxGetClassID(a) == mxDOUBLE_CLASS && 
        mxGetM(a)*mxGetN(a) == 1 && *mxGetPr(a) == 0)
        return p;
    if (mxIsChar(a)) {
        char pbuf[128];
        mxGetString(a, pbuf, sizeof(pbuf));
        sscanf(pbuf, fmt, &p);
    }
    if (p == 0)
        mexErrMsgTxt("Invalid pointer");
    return p;
}

mxArray* mxWrapCreateP(void* p, const char* fmt)
{
    if (p == 0) {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    } else {
        char pbuf[128];
        sprintf(pbuf, fmt, p);
        return mxCreateString(pbuf);
    }
}

mxArray* mxWrapStrncpy(const char* s)
{
    if (s) {
        return mxCreateString(s);
    } else {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    }
}

double mxWrapGetScalar(const mxArray* a)
{
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS || mxGetM(a)*mxGetN(a) != 1)
        mexErrMsgTxt("Invalid scalar argument");
    return *mxGetPr(a);
}

char* mxWrapGetString(const mxArray* a)
{
    char* s;
    int slen;
    if (!a || (!mxIsChar(a) && mxGetM(a)*mxGetN(a) > 0))
        mexErrMsgTxt("Invalid string argument");
    slen = mxGetM(a)*mxGetN(a) + 1;
    s = (char*) mxMalloc(slen);
    if (mxGetM(a)*mxGetN(a) == 0)
        *s = 0;
    else
        mxGetString(a, s, slen);
    return s;
}


#define mxWrapGetArrayDef(func, T) \
T* func(const mxArray* a) \
{ \
    T* array; \
    int arraylen; \
    int i; \
    T* p; \
    double* q; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) \
        mexErrMsgTxt("Invalid array argument"); \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    q = mxGetPr(a); \
    for (i = 0; i < arraylen; ++i) \
        *p++ = (T) (*q++); \
    return array; \
}


#define mxWrapCopyDef(func, T) \
void func(mxArray* a, const T* q, int n) \
{ \
    int i; \
    double* p = mxGetPr(a); \
    for (i = 0; i < n; ++i) \
        *p++ = *q++; \
}


#define mxWrapReturnDef(func, T) \
mxArray* func(const T* q, int m, int n) \
{ \
    int i; \
    double* p; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxREAL); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxREAL); \
        p = mxGetPr(a); \
        for (i = 0; i < m*n; ++i) \
            *p++ = *q++; \
        return a; \
    } \
}


#define mxWrapGetScalarZDef(func, T, ZT, setz) \
void func(T* z, const mxArray* a) \
{ \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    setz(z, (ZT) *pr, (pi ? (ZT) *pi : (ZT) 0)); \
}


#define mxWrapGetArrayZDef(func, T, ZT, setz) \
T* func(const mxArray* a) \
{ \
    T* array; \
    int arraylen; \
    int i; \
    T* p; \
    double* qr; \
    double* qi; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) \
        mexErrMsgTxt("Invalid array argument"); \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    qr = mxGetPr(a); \
    qi = mxGetPi(a); \
    for (i = 0; i < arraylen; ++i) { \
        ZT val_qr = *qr++; \
        ZT val_qi = (qi ? (ZT) *qi++ : (ZT) 0); \
        setz(p, val_qr, val_qi); \
        ++p; \
    } \
    return array; \
}


#define mxWrapCopyZDef(func, T, real, imag) \
void func(mxArray* a, const T* q, int n) \
{ \
    int i; \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    for (i = 0; i < n; ++i) { \
        *pr++ = real(*q); \
        *pi++ = imag(*q); \
        ++q; \
    } \
}


#define mxWrapReturnZDef(func, T, real, imag) \
mxArray* func(const T* q, int m, int n) \
{ \
    int i; \
    double* pr; \
    double* pi; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxCOMPLEX); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxCOMPLEX); \
        pr = mxGetPr(a); \
        pi = mxGetPi(a); \
        for (i = 0; i < m*n; ++i) { \
            *pr++ = real(*q); \
            *pi++ = imag(*q); \
            ++q; \
        } \
        return a; \
    } \
}

/* Array copier definitions */
mxWrapGetArrayDef(mxWrapGetArray_bool, bool)
mxWrapCopyDef    (mxWrapCopy_bool,     bool)
mxWrapReturnDef  (mxWrapReturn_bool,   bool)
mxWrapGetArrayDef(mxWrapGetArray_char, char)
mxWrapCopyDef    (mxWrapCopy_char,     char)
mxWrapReturnDef  (mxWrapReturn_char,   char)
mxWrapGetArrayDef(mxWrapGetArray_double, double)
mxWrapCopyDef    (mxWrapCopy_double,     double)
mxWrapReturnDef  (mxWrapReturn_double,   double)
mxWrapGetArrayDef(mxWrapGetArray_float, float)
mxWrapCopyDef    (mxWrapCopy_float,     float)
mxWrapReturnDef  (mxWrapReturn_float,   float)
mxWrapGetArrayDef(mxWrapGetArray_int, int)
mxWrapCopyDef    (mxWrapCopy_int,     int)
mxWrapReturnDef  (mxWrapReturn_int,   int)
mxWrapGetArrayDef(mxWrapGetArray_long, long)
mxWrapCopyDef    (mxWrapCopy_long,     long)
mxWrapReturnDef  (mxWrapReturn_long,   long)
mxWrapGetArrayDef(mxWrapGetArray_size_t, size_t)
mxWrapCopyDef    (mxWrapCopy_size_t,     size_t)
mxWrapReturnDef  (mxWrapReturn_size_t,   size_t)
mxWrapGetArrayDef(mxWrapGetArray_uchar, uchar)
mxWrapCopyDef    (mxWrapCopy_uchar,     uchar)
mxWrapReturnDef  (mxWrapReturn_uchar,   uchar)
mxWrapGetArrayDef(mxWrapGetArray_uint, uint)
mxWrapCopyDef    (mxWrapCopy_uint,     uint)
mxWrapReturnDef  (mxWrapReturn_uint,   uint)
mxWrapGetArrayDef(mxWrapGetArray_ulong, ulong)
mxWrapCopyDef    (mxWrapCopy_ulong,     ulong)
mxWrapReturnDef  (mxWrapReturn_ulong,   ulong)

/* ---- fem.mw: 14 ----
 * block2d4(double x1, double y1, double x2, double y2, int nx, int ny, output int[4, nix] ix, output double[2, nnodex] nodex);
 */
const char* stubids1_ = "block2d4(i double, i double, i double, i double, i int, i int, o int[xx], o double[xx])";

void mexStub1(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    double      in0_;   /* x1         */
    double      in1_;   /* y1         */
    double      in2_;   /* x2         */
    double      in3_;   /* y2         */
    int         in4_;   /* nx         */
    int         in5_;   /* ny         */
    int*        out0_;  /* ix         */
    double*     out1_;  /* nodex      */
    int         dim6_;  /* 4          */
    int         dim7_;  /* nix        */
    int         dim8_;  /* 2          */
    int         dim9_;  /* nnodex     */

    dim6_ = (int) mxWrapGetScalar(prhs[6]);
    dim7_ = (int) mxWrapGetScalar(prhs[7]);
    dim8_ = (int) mxWrapGetScalar(prhs[8]);
    dim9_ = (int) mxWrapGetScalar(prhs[9]);

    in0_ = (double) mxWrapGetScalar(prhs[0]);
    in1_ = (double) mxWrapGetScalar(prhs[1]);
    in2_ = (double) mxWrapGetScalar(prhs[2]);
    in3_ = (double) mxWrapGetScalar(prhs[3]);
    in4_ = (int) mxWrapGetScalar(prhs[4]);
    in5_ = (int) mxWrapGetScalar(prhs[5]);
    out0_ = (int*) mxMalloc(dim6_*dim7_*sizeof(int));
    out1_ = (double*) mxMalloc(dim8_*dim9_*sizeof(double));
    block2d4(in0_, in1_, in2_, in3_, in4_, in5_, out0_, out1_);
    plhs[0] = mxCreateDoubleMatrix(dim6_, dim7_, mxREAL);
    mxWrapCopy_int(plhs[0], out0_, dim6_*dim7_);
    plhs[1] = mxCreateDoubleMatrix(dim8_, dim9_, mxREAL);
    mxWrapCopy_double(plhs[1], out1_, dim8_*dim9_);
    mxFree(out0_);
    mxFree(out1_);
}

/* ---- fem.mw: 22 ----
 * block2d9(double x1, double y1, double x2, double y2, int nx, int ny, output int[9, nix] ix, output double[2, nnodex] nodex);
 */
const char* stubids2_ = "block2d9(i double, i double, i double, i double, i int, i int, o int[xx], o double[xx])";

void mexStub2(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    double      in0_;   /* x1         */
    double      in1_;   /* y1         */
    double      in2_;   /* x2         */
    double      in3_;   /* y2         */
    int         in4_;   /* nx         */
    int         in5_;   /* ny         */
    int*        out0_;  /* ix         */
    double*     out1_;  /* nodex      */
    int         dim6_;  /* 9          */
    int         dim7_;  /* nix        */
    int         dim8_;  /* 2          */
    int         dim9_;  /* nnodex     */

    dim6_ = (int) mxWrapGetScalar(prhs[6]);
    dim7_ = (int) mxWrapGetScalar(prhs[7]);
    dim8_ = (int) mxWrapGetScalar(prhs[8]);
    dim9_ = (int) mxWrapGetScalar(prhs[9]);

    in0_ = (double) mxWrapGetScalar(prhs[0]);
    in1_ = (double) mxWrapGetScalar(prhs[1]);
    in2_ = (double) mxWrapGetScalar(prhs[2]);
    in3_ = (double) mxWrapGetScalar(prhs[3]);
    in4_ = (int) mxWrapGetScalar(prhs[4]);
    in5_ = (int) mxWrapGetScalar(prhs[5]);
    out0_ = (int*) mxMalloc(dim6_*dim7_*sizeof(int));
    out1_ = (double*) mxMalloc(dim8_*dim9_*sizeof(double));
    block2d9(in0_, in1_, in2_, in3_, in4_, in5_, out0_, out1_);
    plhs[0] = mxCreateDoubleMatrix(dim6_, dim7_, mxREAL);
    mxWrapCopy_int(plhs[0], out0_, dim6_*dim7_);
    plhs[1] = mxCreateDoubleMatrix(dim8_, dim9_, mxREAL);
    mxWrapCopy_double(plhs[1], out1_, dim8_*dim9_);
    mxFree(out0_);
    mxFree(out1_);
}

/* ---- fem.mw: 46 ----
 * element_neohook(double[2, nshape] enodex, double[2, nshape] enodeu, double E, double nu, output double[nshape2, nshape2] Kelt, output double[nshape2] felt, int nshape);
 */
const char* stubids3_ = "element_neohook(i double[xx], i double[xx], i double, i double, o double[xx], o double[x], i int)";

void mexStub3(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    double*     in0_;   /* enodex     */
    double*     in1_;   /* enodeu     */
    double      in2_;   /* E          */
    double      in3_;   /* nu         */
    int         in4_;   /* nshape     */
    double*     out0_;  /* Kelt       */
    double*     out1_;  /* felt       */
    int         dim5_;  /* 2          */
    int         dim6_;  /* nshape     */
    int         dim7_;  /* 2          */
    int         dim8_;  /* nshape     */
    int         dim9_;  /* nshape2    */
    int         dim10_;  /* nshape2    */
    int         dim11_;  /* nshape2    */

    dim5_ = (int) mxWrapGetScalar(prhs[5]);
    dim6_ = (int) mxWrapGetScalar(prhs[6]);
    dim7_ = (int) mxWrapGetScalar(prhs[7]);
    dim8_ = (int) mxWrapGetScalar(prhs[8]);
    dim9_ = (int) mxWrapGetScalar(prhs[9]);
    dim10_ = (int) mxWrapGetScalar(prhs[10]);
    dim11_ = (int) mxWrapGetScalar(prhs[11]);

    if (mxGetM(prhs[0]) != dim5_ ||
        mxGetN(prhs[0]) != dim6_)
        mexErrMsgTxt("Bad argument size: enodex");

    if (mxGetM(prhs[1]) != dim7_ ||
        mxGetN(prhs[1]) != dim8_)
        mexErrMsgTxt("Bad argument size: enodeu");

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0)
        in0_ = mxGetPr(prhs[0]);
    else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0)
        in1_ = mxGetPr(prhs[1]);
    else
        in1_ = NULL;
    in2_ = (double) mxWrapGetScalar(prhs[2]);
    in3_ = (double) mxWrapGetScalar(prhs[3]);
    in4_ = (int) mxWrapGetScalar(prhs[4]);
    out0_ = (double*) mxMalloc(dim9_*dim10_*sizeof(double));
    out1_ = (double*) mxMalloc(dim11_*sizeof(double));
    element_neohook(in0_, in1_, in2_, in3_, out0_, out1_, in4_);
    plhs[0] = mxCreateDoubleMatrix(dim9_, dim10_, mxREAL);
    mxWrapCopy_double(plhs[0], out0_, dim9_*dim10_);
    plhs[1] = mxCreateDoubleMatrix(dim11_, 1, mxREAL);
    mxWrapCopy_double(plhs[1], out1_, dim11_);
    mxFree(out0_);
    mxFree(out1_);
}

/* ----
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    char id[512];
    if (nrhs == 0) {
        mexPrintf("Mex function installed\n");
        return;
    }

    if (mxGetString(prhs[0], id, sizeof(id)) != 0)
        mexErrMsgTxt("Identifier should be a string");
    else if (strcmp(id, stubids1_) == 0)
        mexStub1(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids2_) == 0)
        mexStub2(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids3_) == 0)
        mexStub3(nlhs,plhs, nrhs-1,prhs+1);
    else
        mexErrMsgTxt("Unknown identifier");
}

