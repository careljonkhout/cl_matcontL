#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <math.h>
#include <mex.h>
#include "user_data.h"

#if JACOBIAN_STORAGE == DENSE
#define JAC(i,j) (SM_ELEMENT_D(jacobian, (i), (j)))
#endif

#if JACOBIAN_STORAGE == BANDED
#define JAC(i,j) (SM_ELEMENT_B(jacobian, (i), (j)))
#endif

int jacobian_dydt(
               realtype t, N_Vector y_vec, N_Vector fy, SUNMatrix jacobian,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  
  double* y          = N_VGetArrayPointer(y_vec);
  double* parameters = ((UserData) user_data)->parameters;

  #if JACOBIAN_STORAGE == DENSE
  if (SUNMatGetID(jacobian) != SUNMATRIX_DENSE) {
    mexErrMsgIdAndTxt("cvode:wrong_matrix_type",
            "expected type %d = SUNMATRIX_DENSE, got type %d",
            SUNMATRIX_DENSE, SUNMatGetID(jacobian));
  }
  #endif
  
  #if JACOBIAN_STORAGE == BANDED
  if (SUNMatGetID(jacobian) != SUNMATRIX_BAND) {
    mexErrMsgIdAndTxt("cvode:wrong_matrix_type",
            "expected type %d = SUNMATRIX_BAND, got type %d",
            SUNMATRIX_DENSE, SUNMatGetID(jacobian));
  }
  #endif

  JAC(0,0) =  parameters[1]-parameters[0]*exp(-parameters[5]*(1.0/y[1]-1.0))-(parameters[1]*9.8E+1+4.802E+3)/parameters[1];
  JAC(1,0) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[1]-1.0));
  JAC(2,0) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(0,1) =  -parameters[0]*parameters[5]*y[0]*1.0/(y[1]*y[1])*exp(-parameters[5]*(1.0/y[1]-1.0));
  JAC(1,1) =  -parameters[3]+parameters[2]-(parameters[2]*9.8E+1+4.802E+3)/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[0]*1.0/(y[1]*y[1])*exp(-parameters[5]*(1.0/y[1]-1.0));
  JAC(3,1) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(0,2) =  4.802E+3/parameters[1];
  JAC(2,2) =  -parameters[0]*exp(-parameters[5]*(1.0/y[3]-1.0))-4.802E+3/parameters[1];
  JAC(3,2) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[3]-1.0));
  JAC(4,2) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(1,3) =  4.802E+3/parameters[2];
  JAC(2,3) =  -parameters[0]*parameters[5]*y[2]*1.0/(y[3]*y[3])*exp(-parameters[5]*(1.0/y[3]-1.0));
  JAC(3,3) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[2]*1.0/(y[3]*y[3])*exp(-parameters[5]*(1.0/y[3]-1.0));
  JAC(5,3) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(2,4) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(4,4) =  -parameters[0]*exp(-parameters[5]*(1.0/y[5]-1.0))-4.802E+3/parameters[1];
  JAC(5,4) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[5]-1.0));
  JAC(6,4) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(3,5) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(4,5) =  -parameters[0]*parameters[5]*y[4]*1.0/(y[5]*y[5])*exp(-parameters[5]*(1.0/y[5]-1.0));
  JAC(5,5) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[4]*1.0/(y[5]*y[5])*exp(-parameters[5]*(1.0/y[5]-1.0));
  JAC(7,5) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(4,6) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(6,6) =  -parameters[0]*exp(-parameters[5]*(1.0/y[7]-1.0))-4.802E+3/parameters[1];
  JAC(7,6) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[7]-1.0));
  JAC(8,6) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(5,7) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(6,7) =  -parameters[0]*parameters[5]*y[6]*1.0/(y[7]*y[7])*exp(-parameters[5]*(1.0/y[7]-1.0));
  JAC(7,7) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[6]*1.0/(y[7]*y[7])*exp(-parameters[5]*(1.0/y[7]-1.0));
  JAC(9,7) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(6,8) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(8,8) =  -parameters[0]*exp(-parameters[5]*(1.0/y[9]-1.0))-4.802E+3/parameters[1];
  JAC(9,8) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[9]-1.0));
  JAC(10,8) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(7,9) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(8,9) =  -parameters[0]*parameters[5]*y[8]*1.0/(y[9]*y[9])*exp(-parameters[5]*(1.0/y[9]-1.0));
  JAC(9,9) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[8]*1.0/(y[9]*y[9])*exp(-parameters[5]*(1.0/y[9]-1.0));
  JAC(11,9) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(8,10) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(10,10) =  -parameters[0]*exp(-parameters[5]*(1.0/y[11]-1.0))-4.802E+3/parameters[1];
  JAC(11,10) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[11]-1.0));
  JAC(12,10) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(9,11) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(10,11) =  -parameters[0]*parameters[5]*y[10]*1.0/(y[11]*y[11])*exp(-parameters[5]*(1.0/y[11]-1.0));
  JAC(11,11) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[10]*1.0/(y[11]*y[11])*exp(-parameters[5]*(1.0/y[11]-1.0));
  JAC(13,11) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(10,12) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(12,12) =  -parameters[0]*exp(-parameters[5]*(1.0/y[13]-1.0))-4.802E+3/parameters[1];
  JAC(13,12) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[13]-1.0));
  JAC(14,12) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(11,13) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(12,13) =  -parameters[0]*parameters[5]*y[12]*1.0/(y[13]*y[13])*exp(-parameters[5]*(1.0/y[13]-1.0));
  JAC(13,13) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[12]*1.0/(y[13]*y[13])*exp(-parameters[5]*(1.0/y[13]-1.0));
  JAC(15,13) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(12,14) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(14,14) =  -parameters[0]*exp(-parameters[5]*(1.0/y[15]-1.0))-4.802E+3/parameters[1];
  JAC(15,14) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[15]-1.0));
  JAC(16,14) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(13,15) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(14,15) =  -parameters[0]*parameters[5]*y[14]*1.0/(y[15]*y[15])*exp(-parameters[5]*(1.0/y[15]-1.0));
  JAC(15,15) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[14]*1.0/(y[15]*y[15])*exp(-parameters[5]*(1.0/y[15]-1.0));
  JAC(17,15) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(14,16) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(16,16) =  -parameters[0]*exp(-parameters[5]*(1.0/y[17]-1.0))-4.802E+3/parameters[1];
  JAC(17,16) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[17]-1.0));
  JAC(18,16) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(15,17) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(16,17) =  -parameters[0]*parameters[5]*y[16]*1.0/(y[17]*y[17])*exp(-parameters[5]*(1.0/y[17]-1.0));
  JAC(17,17) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[16]*1.0/(y[17]*y[17])*exp(-parameters[5]*(1.0/y[17]-1.0));
  JAC(19,17) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(16,18) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(18,18) =  -parameters[0]*exp(-parameters[5]*(1.0/y[19]-1.0))-4.802E+3/parameters[1];
  JAC(19,18) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[19]-1.0));
  JAC(20,18) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(17,19) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(18,19) =  -parameters[0]*parameters[5]*y[18]*1.0/(y[19]*y[19])*exp(-parameters[5]*(1.0/y[19]-1.0));
  JAC(19,19) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[18]*1.0/(y[19]*y[19])*exp(-parameters[5]*(1.0/y[19]-1.0));
  JAC(21,19) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(18,20) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(20,20) =  -parameters[0]*exp(-parameters[5]*(1.0/y[21]-1.0))-4.802E+3/parameters[1];
  JAC(21,20) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[21]-1.0));
  JAC(22,20) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(19,21) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(20,21) =  -parameters[0]*parameters[5]*y[20]*1.0/(y[21]*y[21])*exp(-parameters[5]*(1.0/y[21]-1.0));
  JAC(21,21) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[20]*1.0/(y[21]*y[21])*exp(-parameters[5]*(1.0/y[21]-1.0));
  JAC(23,21) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(20,22) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(22,22) =  -parameters[0]*exp(-parameters[5]*(1.0/y[23]-1.0))-4.802E+3/parameters[1];
  JAC(23,22) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[23]-1.0));
  JAC(24,22) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(21,23) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(22,23) =  -parameters[0]*parameters[5]*y[22]*1.0/(y[23]*y[23])*exp(-parameters[5]*(1.0/y[23]-1.0));
  JAC(23,23) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[22]*1.0/(y[23]*y[23])*exp(-parameters[5]*(1.0/y[23]-1.0));
  JAC(25,23) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(22,24) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(24,24) =  -parameters[0]*exp(-parameters[5]*(1.0/y[25]-1.0))-4.802E+3/parameters[1];
  JAC(25,24) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[25]-1.0));
  JAC(26,24) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(23,25) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(24,25) =  -parameters[0]*parameters[5]*y[24]*1.0/(y[25]*y[25])*exp(-parameters[5]*(1.0/y[25]-1.0));
  JAC(25,25) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[24]*1.0/(y[25]*y[25])*exp(-parameters[5]*(1.0/y[25]-1.0));
  JAC(27,25) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(24,26) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(26,26) =  -parameters[0]*exp(-parameters[5]*(1.0/y[27]-1.0))-4.802E+3/parameters[1];
  JAC(27,26) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[27]-1.0));
  JAC(28,26) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(25,27) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(26,27) =  -parameters[0]*parameters[5]*y[26]*1.0/(y[27]*y[27])*exp(-parameters[5]*(1.0/y[27]-1.0));
  JAC(27,27) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[26]*1.0/(y[27]*y[27])*exp(-parameters[5]*(1.0/y[27]-1.0));
  JAC(29,27) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(26,28) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(28,28) =  -parameters[0]*exp(-parameters[5]*(1.0/y[29]-1.0))-4.802E+3/parameters[1];
  JAC(29,28) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[29]-1.0));
  JAC(30,28) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(27,29) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(28,29) =  -parameters[0]*parameters[5]*y[28]*1.0/(y[29]*y[29])*exp(-parameters[5]*(1.0/y[29]-1.0));
  JAC(29,29) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[28]*1.0/(y[29]*y[29])*exp(-parameters[5]*(1.0/y[29]-1.0));
  JAC(31,29) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(28,30) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(30,30) =  -parameters[0]*exp(-parameters[5]*(1.0/y[31]-1.0))-4.802E+3/parameters[1];
  JAC(31,30) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[31]-1.0));
  JAC(32,30) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(29,31) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(30,31) =  -parameters[0]*parameters[5]*y[30]*1.0/(y[31]*y[31])*exp(-parameters[5]*(1.0/y[31]-1.0));
  JAC(31,31) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[30]*1.0/(y[31]*y[31])*exp(-parameters[5]*(1.0/y[31]-1.0));
  JAC(33,31) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(30,32) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(32,32) =  -parameters[0]*exp(-parameters[5]*(1.0/y[33]-1.0))-4.802E+3/parameters[1];
  JAC(33,32) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[33]-1.0));
  JAC(34,32) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(31,33) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(32,33) =  -parameters[0]*parameters[5]*y[32]*1.0/(y[33]*y[33])*exp(-parameters[5]*(1.0/y[33]-1.0));
  JAC(33,33) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[32]*1.0/(y[33]*y[33])*exp(-parameters[5]*(1.0/y[33]-1.0));
  JAC(35,33) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(32,34) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(34,34) =  -parameters[0]*exp(-parameters[5]*(1.0/y[35]-1.0))-4.802E+3/parameters[1];
  JAC(35,34) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[35]-1.0));
  JAC(36,34) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(33,35) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(34,35) =  -parameters[0]*parameters[5]*y[34]*1.0/(y[35]*y[35])*exp(-parameters[5]*(1.0/y[35]-1.0));
  JAC(35,35) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[34]*1.0/(y[35]*y[35])*exp(-parameters[5]*(1.0/y[35]-1.0));
  JAC(37,35) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(34,36) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(36,36) =  -parameters[0]*exp(-parameters[5]*(1.0/y[37]-1.0))-4.802E+3/parameters[1];
  JAC(37,36) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[37]-1.0));
  JAC(38,36) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(35,37) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(36,37) =  -parameters[0]*parameters[5]*y[36]*1.0/(y[37]*y[37])*exp(-parameters[5]*(1.0/y[37]-1.0));
  JAC(37,37) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[36]*1.0/(y[37]*y[37])*exp(-parameters[5]*(1.0/y[37]-1.0));
  JAC(39,37) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(36,38) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(38,38) =  -parameters[0]*exp(-parameters[5]*(1.0/y[39]-1.0))-4.802E+3/parameters[1];
  JAC(39,38) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[39]-1.0));
  JAC(40,38) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(37,39) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(38,39) =  -parameters[0]*parameters[5]*y[38]*1.0/(y[39]*y[39])*exp(-parameters[5]*(1.0/y[39]-1.0));
  JAC(39,39) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[38]*1.0/(y[39]*y[39])*exp(-parameters[5]*(1.0/y[39]-1.0));
  JAC(41,39) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(38,40) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(40,40) =  -parameters[0]*exp(-parameters[5]*(1.0/y[41]-1.0))-4.802E+3/parameters[1];
  JAC(41,40) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[41]-1.0));
  JAC(42,40) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(39,41) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(40,41) =  -parameters[0]*parameters[5]*y[40]*1.0/(y[41]*y[41])*exp(-parameters[5]*(1.0/y[41]-1.0));
  JAC(41,41) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[40]*1.0/(y[41]*y[41])*exp(-parameters[5]*(1.0/y[41]-1.0));
  JAC(43,41) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(40,42) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(42,42) =  -parameters[0]*exp(-parameters[5]*(1.0/y[43]-1.0))-4.802E+3/parameters[1];
  JAC(43,42) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[43]-1.0));
  JAC(44,42) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(41,43) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(42,43) =  -parameters[0]*parameters[5]*y[42]*1.0/(y[43]*y[43])*exp(-parameters[5]*(1.0/y[43]-1.0));
  JAC(43,43) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[42]*1.0/(y[43]*y[43])*exp(-parameters[5]*(1.0/y[43]-1.0));
  JAC(45,43) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(42,44) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(44,44) =  -parameters[0]*exp(-parameters[5]*(1.0/y[45]-1.0))-4.802E+3/parameters[1];
  JAC(45,44) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[45]-1.0));
  JAC(46,44) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(43,45) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(44,45) =  -parameters[0]*parameters[5]*y[44]*1.0/(y[45]*y[45])*exp(-parameters[5]*(1.0/y[45]-1.0));
  JAC(45,45) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[44]*1.0/(y[45]*y[45])*exp(-parameters[5]*(1.0/y[45]-1.0));
  JAC(47,45) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(44,46) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(46,46) =  -parameters[0]*exp(-parameters[5]*(1.0/y[47]-1.0))-4.802E+3/parameters[1];
  JAC(47,46) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[47]-1.0));
  JAC(48,46) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(45,47) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(46,47) =  -parameters[0]*parameters[5]*y[46]*1.0/(y[47]*y[47])*exp(-parameters[5]*(1.0/y[47]-1.0));
  JAC(47,47) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[46]*1.0/(y[47]*y[47])*exp(-parameters[5]*(1.0/y[47]-1.0));
  JAC(49,47) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(46,48) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(48,48) =  -parameters[0]*exp(-parameters[5]*(1.0/y[49]-1.0))-4.802E+3/parameters[1];
  JAC(49,48) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[49]-1.0));
  JAC(50,48) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(47,49) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(48,49) =  -parameters[0]*parameters[5]*y[48]*1.0/(y[49]*y[49])*exp(-parameters[5]*(1.0/y[49]-1.0));
  JAC(49,49) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[48]*1.0/(y[49]*y[49])*exp(-parameters[5]*(1.0/y[49]-1.0));
  JAC(51,49) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(48,50) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(50,50) =  -parameters[0]*exp(-parameters[5]*(1.0/y[51]-1.0))-4.802E+3/parameters[1];
  JAC(51,50) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[51]-1.0));
  JAC(52,50) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(49,51) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(50,51) =  -parameters[0]*parameters[5]*y[50]*1.0/(y[51]*y[51])*exp(-parameters[5]*(1.0/y[51]-1.0));
  JAC(51,51) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[50]*1.0/(y[51]*y[51])*exp(-parameters[5]*(1.0/y[51]-1.0));
  JAC(53,51) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(50,52) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(52,52) =  -parameters[0]*exp(-parameters[5]*(1.0/y[53]-1.0))-4.802E+3/parameters[1];
  JAC(53,52) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[53]-1.0));
  JAC(54,52) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(51,53) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(52,53) =  -parameters[0]*parameters[5]*y[52]*1.0/(y[53]*y[53])*exp(-parameters[5]*(1.0/y[53]-1.0));
  JAC(53,53) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[52]*1.0/(y[53]*y[53])*exp(-parameters[5]*(1.0/y[53]-1.0));
  JAC(55,53) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(52,54) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(54,54) =  -parameters[0]*exp(-parameters[5]*(1.0/y[55]-1.0))-4.802E+3/parameters[1];
  JAC(55,54) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[55]-1.0));
  JAC(56,54) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(53,55) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(54,55) =  -parameters[0]*parameters[5]*y[54]*1.0/(y[55]*y[55])*exp(-parameters[5]*(1.0/y[55]-1.0));
  JAC(55,55) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[54]*1.0/(y[55]*y[55])*exp(-parameters[5]*(1.0/y[55]-1.0));
  JAC(57,55) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(54,56) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(56,56) =  -parameters[0]*exp(-parameters[5]*(1.0/y[57]-1.0))-4.802E+3/parameters[1];
  JAC(57,56) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[57]-1.0));
  JAC(58,56) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(55,57) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(56,57) =  -parameters[0]*parameters[5]*y[56]*1.0/(y[57]*y[57])*exp(-parameters[5]*(1.0/y[57]-1.0));
  JAC(57,57) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[56]*1.0/(y[57]*y[57])*exp(-parameters[5]*(1.0/y[57]-1.0));
  JAC(59,57) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(56,58) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(58,58) =  -parameters[0]*exp(-parameters[5]*(1.0/y[59]-1.0))-4.802E+3/parameters[1];
  JAC(59,58) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[59]-1.0));
  JAC(60,58) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(57,59) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(58,59) =  -parameters[0]*parameters[5]*y[58]*1.0/(y[59]*y[59])*exp(-parameters[5]*(1.0/y[59]-1.0));
  JAC(59,59) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[58]*1.0/(y[59]*y[59])*exp(-parameters[5]*(1.0/y[59]-1.0));
  JAC(61,59) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(58,60) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(60,60) =  -parameters[0]*exp(-parameters[5]*(1.0/y[61]-1.0))-4.802E+3/parameters[1];
  JAC(61,60) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[61]-1.0));
  JAC(62,60) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(59,61) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(60,61) =  -parameters[0]*parameters[5]*y[60]*1.0/(y[61]*y[61])*exp(-parameters[5]*(1.0/y[61]-1.0));
  JAC(61,61) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[60]*1.0/(y[61]*y[61])*exp(-parameters[5]*(1.0/y[61]-1.0));
  JAC(63,61) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(60,62) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(62,62) =  -parameters[0]*exp(-parameters[5]*(1.0/y[63]-1.0))-4.802E+3/parameters[1];
  JAC(63,62) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[63]-1.0));
  JAC(64,62) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(61,63) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(62,63) =  -parameters[0]*parameters[5]*y[62]*1.0/(y[63]*y[63])*exp(-parameters[5]*(1.0/y[63]-1.0));
  JAC(63,63) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[62]*1.0/(y[63]*y[63])*exp(-parameters[5]*(1.0/y[63]-1.0));
  JAC(65,63) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(62,64) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(64,64) =  -parameters[0]*exp(-parameters[5]*(1.0/y[65]-1.0))-4.802E+3/parameters[1];
  JAC(65,64) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[65]-1.0));
  JAC(66,64) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(63,65) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(64,65) =  -parameters[0]*parameters[5]*y[64]*1.0/(y[65]*y[65])*exp(-parameters[5]*(1.0/y[65]-1.0));
  JAC(65,65) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[64]*1.0/(y[65]*y[65])*exp(-parameters[5]*(1.0/y[65]-1.0));
  JAC(67,65) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(64,66) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(66,66) =  -parameters[0]*exp(-parameters[5]*(1.0/y[67]-1.0))-4.802E+3/parameters[1];
  JAC(67,66) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[67]-1.0));
  JAC(68,66) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(65,67) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(66,67) =  -parameters[0]*parameters[5]*y[66]*1.0/(y[67]*y[67])*exp(-parameters[5]*(1.0/y[67]-1.0));
  JAC(67,67) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[66]*1.0/(y[67]*y[67])*exp(-parameters[5]*(1.0/y[67]-1.0));
  JAC(69,67) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(66,68) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(68,68) =  -parameters[0]*exp(-parameters[5]*(1.0/y[69]-1.0))-4.802E+3/parameters[1];
  JAC(69,68) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[69]-1.0));
  JAC(70,68) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(67,69) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(68,69) =  -parameters[0]*parameters[5]*y[68]*1.0/(y[69]*y[69])*exp(-parameters[5]*(1.0/y[69]-1.0));
  JAC(69,69) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[68]*1.0/(y[69]*y[69])*exp(-parameters[5]*(1.0/y[69]-1.0));
  JAC(71,69) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(68,70) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(70,70) =  -parameters[0]*exp(-parameters[5]*(1.0/y[71]-1.0))-4.802E+3/parameters[1];
  JAC(71,70) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[71]-1.0));
  JAC(72,70) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(69,71) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(70,71) =  -parameters[0]*parameters[5]*y[70]*1.0/(y[71]*y[71])*exp(-parameters[5]*(1.0/y[71]-1.0));
  JAC(71,71) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[70]*1.0/(y[71]*y[71])*exp(-parameters[5]*(1.0/y[71]-1.0));
  JAC(73,71) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(70,72) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(72,72) =  -parameters[0]*exp(-parameters[5]*(1.0/y[73]-1.0))-4.802E+3/parameters[1];
  JAC(73,72) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[73]-1.0));
  JAC(74,72) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(71,73) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(72,73) =  -parameters[0]*parameters[5]*y[72]*1.0/(y[73]*y[73])*exp(-parameters[5]*(1.0/y[73]-1.0));
  JAC(73,73) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[72]*1.0/(y[73]*y[73])*exp(-parameters[5]*(1.0/y[73]-1.0));
  JAC(75,73) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(72,74) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(74,74) =  -parameters[0]*exp(-parameters[5]*(1.0/y[75]-1.0))-4.802E+3/parameters[1];
  JAC(75,74) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[75]-1.0));
  JAC(76,74) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(73,75) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(74,75) =  -parameters[0]*parameters[5]*y[74]*1.0/(y[75]*y[75])*exp(-parameters[5]*(1.0/y[75]-1.0));
  JAC(75,75) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[74]*1.0/(y[75]*y[75])*exp(-parameters[5]*(1.0/y[75]-1.0));
  JAC(77,75) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(74,76) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(76,76) =  -parameters[0]*exp(-parameters[5]*(1.0/y[77]-1.0))-4.802E+3/parameters[1];
  JAC(77,76) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[77]-1.0));
  JAC(78,76) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(75,77) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(76,77) =  -parameters[0]*parameters[5]*y[76]*1.0/(y[77]*y[77])*exp(-parameters[5]*(1.0/y[77]-1.0));
  JAC(77,77) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[76]*1.0/(y[77]*y[77])*exp(-parameters[5]*(1.0/y[77]-1.0));
  JAC(79,77) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(76,78) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(78,78) =  -parameters[0]*exp(-parameters[5]*(1.0/y[79]-1.0))-4.802E+3/parameters[1];
  JAC(79,78) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[79]-1.0));
  JAC(80,78) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(77,79) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(78,79) =  -parameters[0]*parameters[5]*y[78]*1.0/(y[79]*y[79])*exp(-parameters[5]*(1.0/y[79]-1.0));
  JAC(79,79) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[78]*1.0/(y[79]*y[79])*exp(-parameters[5]*(1.0/y[79]-1.0));
  JAC(81,79) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(78,80) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(80,80) =  -parameters[0]*exp(-parameters[5]*(1.0/y[81]-1.0))-4.802E+3/parameters[1];
  JAC(81,80) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[81]-1.0));
  JAC(82,80) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(79,81) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(80,81) =  -parameters[0]*parameters[5]*y[80]*1.0/(y[81]*y[81])*exp(-parameters[5]*(1.0/y[81]-1.0));
  JAC(81,81) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[80]*1.0/(y[81]*y[81])*exp(-parameters[5]*(1.0/y[81]-1.0));
  JAC(83,81) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(80,82) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(82,82) =  -parameters[0]*exp(-parameters[5]*(1.0/y[83]-1.0))-4.802E+3/parameters[1];
  JAC(83,82) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[83]-1.0));
  JAC(84,82) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(81,83) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(82,83) =  -parameters[0]*parameters[5]*y[82]*1.0/(y[83]*y[83])*exp(-parameters[5]*(1.0/y[83]-1.0));
  JAC(83,83) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[82]*1.0/(y[83]*y[83])*exp(-parameters[5]*(1.0/y[83]-1.0));
  JAC(85,83) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(82,84) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(84,84) =  -parameters[0]*exp(-parameters[5]*(1.0/y[85]-1.0))-4.802E+3/parameters[1];
  JAC(85,84) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[85]-1.0));
  JAC(86,84) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(83,85) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(84,85) =  -parameters[0]*parameters[5]*y[84]*1.0/(y[85]*y[85])*exp(-parameters[5]*(1.0/y[85]-1.0));
  JAC(85,85) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[84]*1.0/(y[85]*y[85])*exp(-parameters[5]*(1.0/y[85]-1.0));
  JAC(87,85) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(84,86) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(86,86) =  -parameters[0]*exp(-parameters[5]*(1.0/y[87]-1.0))-4.802E+3/parameters[1];
  JAC(87,86) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[87]-1.0));
  JAC(88,86) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(85,87) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(86,87) =  -parameters[0]*parameters[5]*y[86]*1.0/(y[87]*y[87])*exp(-parameters[5]*(1.0/y[87]-1.0));
  JAC(87,87) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[86]*1.0/(y[87]*y[87])*exp(-parameters[5]*(1.0/y[87]-1.0));
  JAC(89,87) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(86,88) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(88,88) =  -parameters[0]*exp(-parameters[5]*(1.0/y[89]-1.0))-4.802E+3/parameters[1];
  JAC(89,88) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[89]-1.0));
  JAC(90,88) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(87,89) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(88,89) =  -parameters[0]*parameters[5]*y[88]*1.0/(y[89]*y[89])*exp(-parameters[5]*(1.0/y[89]-1.0));
  JAC(89,89) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[88]*1.0/(y[89]*y[89])*exp(-parameters[5]*(1.0/y[89]-1.0));
  JAC(91,89) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(88,90) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(90,90) =  -parameters[0]*exp(-parameters[5]*(1.0/y[91]-1.0))-4.802E+3/parameters[1];
  JAC(91,90) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[91]-1.0));
  JAC(92,90) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(89,91) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(90,91) =  -parameters[0]*parameters[5]*y[90]*1.0/(y[91]*y[91])*exp(-parameters[5]*(1.0/y[91]-1.0));
  JAC(91,91) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[90]*1.0/(y[91]*y[91])*exp(-parameters[5]*(1.0/y[91]-1.0));
  JAC(93,91) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(90,92) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(92,92) =  -parameters[0]*exp(-parameters[5]*(1.0/y[93]-1.0))-4.802E+3/parameters[1];
  JAC(93,92) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[93]-1.0));
  JAC(94,92) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(91,93) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(92,93) =  -parameters[0]*parameters[5]*y[92]*1.0/(y[93]*y[93])*exp(-parameters[5]*(1.0/y[93]-1.0));
  JAC(93,93) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[92]*1.0/(y[93]*y[93])*exp(-parameters[5]*(1.0/y[93]-1.0));
  JAC(95,93) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(92,94) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(94,94) =  -parameters[0]*exp(-parameters[5]*(1.0/y[95]-1.0))-4.802E+3/parameters[1];
  JAC(95,94) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[95]-1.0));
  JAC(96,94) =  2.401E+3/parameters[1]+4.9E+1/2.0;
  JAC(93,95) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(94,95) =  -parameters[0]*parameters[5]*y[94]*1.0/(y[95]*y[95])*exp(-parameters[5]*(1.0/y[95]-1.0));
  JAC(95,95) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[94]*1.0/(y[95]*y[95])*exp(-parameters[5]*(1.0/y[95]-1.0));
  JAC(97,95) =  2.401E+3/parameters[2]+4.9E+1/2.0;
  JAC(94,96) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(96,96) =  -parameters[0]*exp(-parameters[5]*(1.0/y[97]-1.0))-4.802E+3/parameters[1];
  JAC(97,96) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[97]-1.0));
  JAC(98,96) =  4.802E+3/parameters[1];
  JAC(95,97) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(96,97) =  -parameters[0]*parameters[5]*y[96]*1.0/(y[97]*y[97])*exp(-parameters[5]*(1.0/y[97]-1.0));
  JAC(97,97) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[96]*1.0/(y[97]*y[97])*exp(-parameters[5]*(1.0/y[97]-1.0));
  JAC(99,97) =  4.802E+3/parameters[2];
  JAC(96,98) =  2.401E+3/parameters[1]-4.9E+1/2.0;
  JAC(98,98) =  -parameters[0]*exp(-parameters[5]*(1.0/y[99]-1.0))-4.802E+3/parameters[1];
  JAC(99,98) =  parameters[6]*parameters[0]*exp(-parameters[5]*(1.0/y[99]-1.0));
  JAC(97,99) =  2.401E+3/parameters[2]-4.9E+1/2.0;
  JAC(98,99) =  -parameters[0]*parameters[5]*y[98]*1.0/(y[99]*y[99])*exp(-parameters[5]*(1.0/y[99]-1.0));
  JAC(99,99) =  -parameters[3]-4.802E+3/parameters[2]+parameters[6]*parameters[0]*parameters[5]*y[98]*1.0/(y[99]*y[99])*exp(-parameters[5]*(1.0/y[99]-1.0));
  
  return 0;
}