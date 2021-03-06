#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <math.h>
#include <mex.h>
#include "user_data.h"

#define X(i) u[2*(i)]
#define Y(i) u[2*(i)+1]

#define SX(i) s[2*(i)]
#define SY(i) s[2*(i)+1]

#define DSX_DT(i) dsdt[2*(i)]
#define DSY_DT(i) dsdt[2*(i)+1]

#define L  parameters[0]
#define A  parameters[1]
#define B  parameters[2]
#define DX parameters[3]
#define DY parameters[4]


int d_sensitivity_dt(int Ns, realtype t,
                N_Vector y_vector, 
                N_Vector ydot,
                int iS,
                N_Vector s_vector, 
                N_Vector dsdt_vector,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2) {
   
  //double* y        = N_VGetArrayPointer(   y_vector);
  double* u          = N_VGetArrayPointer(   y_vector);
  double* s          = N_VGetArrayPointer(   s_vector);
  double* dsdt       = N_VGetArrayPointer(dsdt_vector);
  double* parameters = ((UserData) user_data)->parameters;
    
  double cx = DX * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
  double cy = DY * (N_MESH_POINTS + 1) * (N_MESH_POINTS+1) / (L*L);
  
  double jac_xx = - 2 * cx - (B+1) + 2 * X(0) * Y(0);
  double jac_xy = X(0) * X(0);
  DSX_DT(0)     = jac_xx * SX(0) + cx * SX(1) + jac_xy * SY(0);
  
	double jac_yy = -2 * cy - X(0) * X(0);
  double jac_yx = B - 2 * X(0) * Y(0);
	DSY_DT(0)     = jac_yy * SY(0) + cy * SY(1) + jac_yx * SX(0);
  
  for (int i = 1; i < N_MESH_POINTS - 1; i++) {
    jac_xx    = - 2 * cx     - (B+1) + 2 * X(i) * Y(i);
    jac_xy    = X(i) * X(i);
    DSX_DT(i) = cx * SX(i-1) + jac_xx * SX(i) + cx * SX(i+1) + jac_xy * SY(i);
    
    jac_yy    = - 2 * cy         - X(i) * X(i);
    jac_yx    = B - 2 * X(i) * Y(i);
    DSY_DT(i) = cy * SY(i-1) + jac_yy * SY(i) + cy * SY(i+1) + jac_yx * SX(i);
  }
  
  int i = N_MESH_POINTS - 1;
  
  jac_xx    = - 2 * cx     - (B+1) + 2 * X(i) * Y(i);
  jac_xy    = X(i) * X(i);
  DSX_DT(i) = cx * SX(i-1) + jac_xx * SX(i) + jac_xy * SY(i);
  
  jac_yy    = - 2 * cy         - X(i) * X(i);
  jac_yx    = B - 2 * X(i) * Y(i);
  DSY_DT(i) = cy * SY(i-1) + jac_yy * SY(i) + jac_yx * SX(i);
 
  
  double ds[100];

  /*
  ds[0] =  -s[0]*(B-y[0]*y[1]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[1]*(y[0]*y[0])+DX*1.0/(L*L)*s[2]*2.601E+3;
  ds[1] =  -s[1]*(DY*1.0/(L*L)*5.202E+3+y[0]*y[0])+s[0]*(B-y[0]*y[1]*2.0)+DY*1.0/(L*L)*s[3]*2.601E+3;
  ds[2] =  -s[2]*(B - y[2]*y[3]*2.0 + DX*1.0/(L*L)*5.202E+3 + 1.0)   +   s[3]*(y[2]*y[2])+DX*1.0/(L*L)*s[0]*2.601E+3 + DX*1.0/(L*L)*s[4]*2.601E+3;
  ds[3] =  -s[3]*(DY*1.0/(L*L)*5.202E+3+y[2]*y[2])+s[2]*(B-y[2]*y[3]*2.0)+DY*1.0/(L*L)*s[1]*2.601E+3+DY*1.0/(L*L)*s[5]*2.601E+3;
  ds[4] =  -s[4]*(B-y[4]*y[5]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[5]*(y[4]*y[4])+DX*1.0/(L*L)*s[2]*2.601E+3+DX*1.0/(L*L)*s[6]*2.601E+3;
  ds[5] =  -s[5]*(DY*1.0/(L*L)*5.202E+3+y[4]*y[4])+s[4]*(B-y[4]*y[5]*2.0)+DY*1.0/(L*L)*s[3]*2.601E+3+DY*1.0/(L*L)*s[7]*2.601E+3;
  ds[6] =  -s[6]*(B-y[6]*y[7]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[7]*(y[6]*y[6])+DX*1.0/(L*L)*s[4]*2.601E+3+DX*1.0/(L*L)*s[8]*2.601E+3;
  ds[7] =  -s[7]*(DY*1.0/(L*L)*5.202E+3+y[6]*y[6])+s[6]*(B-y[6]*y[7]*2.0)+DY*1.0/(L*L)*s[5]*2.601E+3+DY*1.0/(L*L)*s[9]*2.601E+3;
  ds[8] =  -s[8]*(B-y[8]*y[9]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[9]*(y[8]*y[8])+DX*1.0/(L*L)*s[6]*2.601E+3+DX*1.0/(L*L)*s[10]*2.601E+3;
  ds[9] =  -s[9]*(DY*1.0/(L*L)*5.202E+3+y[8]*y[8])+s[8]*(B-y[8]*y[9]*2.0)+DY*1.0/(L*L)*s[7]*2.601E+3+DY*1.0/(L*L)*s[11]*2.601E+3;
  ds[10] =  -s[10]*(B-y[10]*y[11]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[11]*(y[10]*y[10])+DX*1.0/(L*L)*s[8]*2.601E+3+DX*1.0/(L*L)*s[12]*2.601E+3;
  ds[11] =  -s[11]*(DY*1.0/(L*L)*5.202E+3+y[10]*y[10])+s[10]*(B-y[10]*y[11]*2.0)+DY*1.0/(L*L)*s[9]*2.601E+3+DY*1.0/(L*L)*s[13]*2.601E+3;
  ds[12] =  -s[12]*(B-y[12]*y[13]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[13]*(y[12]*y[12])+DX*1.0/(L*L)*s[10]*2.601E+3+DX*1.0/(L*L)*s[14]*2.601E+3;
  ds[13] =  -s[13]*(DY*1.0/(L*L)*5.202E+3+y[12]*y[12])+s[12]*(B-y[12]*y[13]*2.0)+DY*1.0/(L*L)*s[11]*2.601E+3+DY*1.0/(L*L)*s[15]*2.601E+3;
  ds[14] =  -s[14]*(B-y[14]*y[15]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[15]*(y[14]*y[14])+DX*1.0/(L*L)*s[12]*2.601E+3+DX*1.0/(L*L)*s[16]*2.601E+3;
  ds[15] =  -s[15]*(DY*1.0/(L*L)*5.202E+3+y[14]*y[14])+s[14]*(B-y[14]*y[15]*2.0)+DY*1.0/(L*L)*s[13]*2.601E+3+DY*1.0/(L*L)*s[17]*2.601E+3;
  ds[16] =  -s[16]*(B-y[16]*y[17]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[17]*(y[16]*y[16])+DX*1.0/(L*L)*s[14]*2.601E+3+DX*1.0/(L*L)*s[18]*2.601E+3;
  ds[17] =  -s[17]*(DY*1.0/(L*L)*5.202E+3+y[16]*y[16])+s[16]*(B-y[16]*y[17]*2.0)+DY*1.0/(L*L)*s[15]*2.601E+3+DY*1.0/(L*L)*s[19]*2.601E+3;
  ds[18] =  -s[18]*(B-y[18]*y[19]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[19]*(y[18]*y[18])+DX*1.0/(L*L)*s[16]*2.601E+3+DX*1.0/(L*L)*s[20]*2.601E+3;
  ds[19] =  -s[19]*(DY*1.0/(L*L)*5.202E+3+y[18]*y[18])+s[18]*(B-y[18]*y[19]*2.0)+DY*1.0/(L*L)*s[17]*2.601E+3+DY*1.0/(L*L)*s[21]*2.601E+3;
  ds[20] =  -s[20]*(B-y[20]*y[21]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[21]*(y[20]*y[20])+DX*1.0/(L*L)*s[18]*2.601E+3+DX*1.0/(L*L)*s[22]*2.601E+3;
  ds[21] =  -s[21]*(DY*1.0/(L*L)*5.202E+3+y[20]*y[20])+s[20]*(B-y[20]*y[21]*2.0)+DY*1.0/(L*L)*s[19]*2.601E+3+DY*1.0/(L*L)*s[23]*2.601E+3;
  ds[22] =  -s[22]*(B-y[22]*y[23]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[23]*(y[22]*y[22])+DX*1.0/(L*L)*s[20]*2.601E+3+DX*1.0/(L*L)*s[24]*2.601E+3;
  ds[23] =  -s[23]*(DY*1.0/(L*L)*5.202E+3+y[22]*y[22])+s[22]*(B-y[22]*y[23]*2.0)+DY*1.0/(L*L)*s[21]*2.601E+3+DY*1.0/(L*L)*s[25]*2.601E+3;
  ds[24] =  -s[24]*(B-y[24]*y[25]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[25]*(y[24]*y[24])+DX*1.0/(L*L)*s[22]*2.601E+3+DX*1.0/(L*L)*s[26]*2.601E+3;
  ds[25] =  -s[25]*(DY*1.0/(L*L)*5.202E+3+y[24]*y[24])+s[24]*(B-y[24]*y[25]*2.0)+DY*1.0/(L*L)*s[23]*2.601E+3+DY*1.0/(L*L)*s[27]*2.601E+3;
  ds[26] =  -s[26]*(B-y[26]*y[27]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[27]*(y[26]*y[26])+DX*1.0/(L*L)*s[24]*2.601E+3+DX*1.0/(L*L)*s[28]*2.601E+3;
  ds[27] =  -s[27]*(DY*1.0/(L*L)*5.202E+3+y[26]*y[26])+s[26]*(B-y[26]*y[27]*2.0)+DY*1.0/(L*L)*s[25]*2.601E+3+DY*1.0/(L*L)*s[29]*2.601E+3;
  ds[28] =  -s[28]*(B-y[28]*y[29]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[29]*(y[28]*y[28])+DX*1.0/(L*L)*s[26]*2.601E+3+DX*1.0/(L*L)*s[30]*2.601E+3;
  ds[29] =  -s[29]*(DY*1.0/(L*L)*5.202E+3+y[28]*y[28])+s[28]*(B-y[28]*y[29]*2.0)+DY*1.0/(L*L)*s[27]*2.601E+3+DY*1.0/(L*L)*s[31]*2.601E+3;
  ds[30] =  -s[30]*(B-y[30]*y[31]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[31]*(y[30]*y[30])+DX*1.0/(L*L)*s[28]*2.601E+3+DX*1.0/(L*L)*s[32]*2.601E+3;
  ds[31] =  -s[31]*(DY*1.0/(L*L)*5.202E+3+y[30]*y[30])+s[30]*(B-y[30]*y[31]*2.0)+DY*1.0/(L*L)*s[29]*2.601E+3+DY*1.0/(L*L)*s[33]*2.601E+3;
  ds[32] =  -s[32]*(B-y[32]*y[33]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[33]*(y[32]*y[32])+DX*1.0/(L*L)*s[30]*2.601E+3+DX*1.0/(L*L)*s[34]*2.601E+3;
  ds[33] =  -s[33]*(DY*1.0/(L*L)*5.202E+3+y[32]*y[32])+s[32]*(B-y[32]*y[33]*2.0)+DY*1.0/(L*L)*s[31]*2.601E+3+DY*1.0/(L*L)*s[35]*2.601E+3;
  ds[34] =  -s[34]*(B-y[34]*y[35]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[35]*(y[34]*y[34])+DX*1.0/(L*L)*s[32]*2.601E+3+DX*1.0/(L*L)*s[36]*2.601E+3;
  ds[35] =  -s[35]*(DY*1.0/(L*L)*5.202E+3+y[34]*y[34])+s[34]*(B-y[34]*y[35]*2.0)+DY*1.0/(L*L)*s[33]*2.601E+3+DY*1.0/(L*L)*s[37]*2.601E+3;
  ds[36] =  -s[36]*(B-y[36]*y[37]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[37]*(y[36]*y[36])+DX*1.0/(L*L)*s[34]*2.601E+3+DX*1.0/(L*L)*s[38]*2.601E+3;
  ds[37] =  -s[37]*(DY*1.0/(L*L)*5.202E+3+y[36]*y[36])+s[36]*(B-y[36]*y[37]*2.0)+DY*1.0/(L*L)*s[35]*2.601E+3+DY*1.0/(L*L)*s[39]*2.601E+3;
  ds[38] =  -s[38]*(B-y[38]*y[39]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[39]*(y[38]*y[38])+DX*1.0/(L*L)*s[36]*2.601E+3+DX*1.0/(L*L)*s[40]*2.601E+3;
  ds[39] =  -s[39]*(DY*1.0/(L*L)*5.202E+3+y[38]*y[38])+s[38]*(B-y[38]*y[39]*2.0)+DY*1.0/(L*L)*s[37]*2.601E+3+DY*1.0/(L*L)*s[41]*2.601E+3;
  ds[40] =  -s[40]*(B-y[40]*y[41]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[41]*(y[40]*y[40])+DX*1.0/(L*L)*s[38]*2.601E+3+DX*1.0/(L*L)*s[42]*2.601E+3;
  ds[41] =  -s[41]*(DY*1.0/(L*L)*5.202E+3+y[40]*y[40])+s[40]*(B-y[40]*y[41]*2.0)+DY*1.0/(L*L)*s[39]*2.601E+3+DY*1.0/(L*L)*s[43]*2.601E+3;
  ds[42] =  -s[42]*(B-y[42]*y[43]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[43]*(y[42]*y[42])+DX*1.0/(L*L)*s[40]*2.601E+3+DX*1.0/(L*L)*s[44]*2.601E+3;
  ds[43] =  -s[43]*(DY*1.0/(L*L)*5.202E+3+y[42]*y[42])+s[42]*(B-y[42]*y[43]*2.0)+DY*1.0/(L*L)*s[41]*2.601E+3+DY*1.0/(L*L)*s[45]*2.601E+3;
  ds[44] =  -s[44]*(B-y[44]*y[45]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[45]*(y[44]*y[44])+DX*1.0/(L*L)*s[42]*2.601E+3+DX*1.0/(L*L)*s[46]*2.601E+3;
  ds[45] =  -s[45]*(DY*1.0/(L*L)*5.202E+3+y[44]*y[44])+s[44]*(B-y[44]*y[45]*2.0)+DY*1.0/(L*L)*s[43]*2.601E+3+DY*1.0/(L*L)*s[47]*2.601E+3;
  ds[46] =  -s[46]*(B-y[46]*y[47]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[47]*(y[46]*y[46])+DX*1.0/(L*L)*s[44]*2.601E+3+DX*1.0/(L*L)*s[48]*2.601E+3;
  ds[47] =  -s[47]*(DY*1.0/(L*L)*5.202E+3+y[46]*y[46])+s[46]*(B-y[46]*y[47]*2.0)+DY*1.0/(L*L)*s[45]*2.601E+3+DY*1.0/(L*L)*s[49]*2.601E+3;
  ds[48] =  -s[48]*(B-y[48]*y[49]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[49]*(y[48]*y[48])+DX*1.0/(L*L)*s[46]*2.601E+3+DX*1.0/(L*L)*s[50]*2.601E+3;
  ds[49] =  -s[49]*(DY*1.0/(L*L)*5.202E+3+y[48]*y[48])+s[48]*(B-y[48]*y[49]*2.0)+DY*1.0/(L*L)*s[47]*2.601E+3+DY*1.0/(L*L)*s[51]*2.601E+3;
  ds[50] =  -s[50]*(B-y[50]*y[51]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[51]*(y[50]*y[50])+DX*1.0/(L*L)*s[48]*2.601E+3+DX*1.0/(L*L)*s[52]*2.601E+3;
  ds[51] =  -s[51]*(DY*1.0/(L*L)*5.202E+3+y[50]*y[50])+s[50]*(B-y[50]*y[51]*2.0)+DY*1.0/(L*L)*s[49]*2.601E+3+DY*1.0/(L*L)*s[53]*2.601E+3;
  ds[52] =  -s[52]*(B-y[52]*y[53]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[53]*(y[52]*y[52])+DX*1.0/(L*L)*s[50]*2.601E+3+DX*1.0/(L*L)*s[54]*2.601E+3;
  ds[53] =  -s[53]*(DY*1.0/(L*L)*5.202E+3+y[52]*y[52])+s[52]*(B-y[52]*y[53]*2.0)+DY*1.0/(L*L)*s[51]*2.601E+3+DY*1.0/(L*L)*s[55]*2.601E+3;
  ds[54] =  -s[54]*(B-y[54]*y[55]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[55]*(y[54]*y[54])+DX*1.0/(L*L)*s[52]*2.601E+3+DX*1.0/(L*L)*s[56]*2.601E+3;
  ds[55] =  -s[55]*(DY*1.0/(L*L)*5.202E+3+y[54]*y[54])+s[54]*(B-y[54]*y[55]*2.0)+DY*1.0/(L*L)*s[53]*2.601E+3+DY*1.0/(L*L)*s[57]*2.601E+3;
  ds[56] =  -s[56]*(B-y[56]*y[57]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[57]*(y[56]*y[56])+DX*1.0/(L*L)*s[54]*2.601E+3+DX*1.0/(L*L)*s[58]*2.601E+3;
  ds[57] =  -s[57]*(DY*1.0/(L*L)*5.202E+3+y[56]*y[56])+s[56]*(B-y[56]*y[57]*2.0)+DY*1.0/(L*L)*s[55]*2.601E+3+DY*1.0/(L*L)*s[59]*2.601E+3;
  ds[58] =  -s[58]*(B-y[58]*y[59]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[59]*(y[58]*y[58])+DX*1.0/(L*L)*s[56]*2.601E+3+DX*1.0/(L*L)*s[60]*2.601E+3;
  ds[59] =  -s[59]*(DY*1.0/(L*L)*5.202E+3+y[58]*y[58])+s[58]*(B-y[58]*y[59]*2.0)+DY*1.0/(L*L)*s[57]*2.601E+3+DY*1.0/(L*L)*s[61]*2.601E+3;
  ds[60] =  -s[60]*(B-y[60]*y[61]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[61]*(y[60]*y[60])+DX*1.0/(L*L)*s[58]*2.601E+3+DX*1.0/(L*L)*s[62]*2.601E+3;
  ds[61] =  -s[61]*(DY*1.0/(L*L)*5.202E+3+y[60]*y[60])+s[60]*(B-y[60]*y[61]*2.0)+DY*1.0/(L*L)*s[59]*2.601E+3+DY*1.0/(L*L)*s[63]*2.601E+3;
  ds[62] =  -s[62]*(B-y[62]*y[63]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[63]*(y[62]*y[62])+DX*1.0/(L*L)*s[60]*2.601E+3+DX*1.0/(L*L)*s[64]*2.601E+3;
  ds[63] =  -s[63]*(DY*1.0/(L*L)*5.202E+3+y[62]*y[62])+s[62]*(B-y[62]*y[63]*2.0)+DY*1.0/(L*L)*s[61]*2.601E+3+DY*1.0/(L*L)*s[65]*2.601E+3;
  ds[64] =  -s[64]*(B-y[64]*y[65]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[65]*(y[64]*y[64])+DX*1.0/(L*L)*s[62]*2.601E+3+DX*1.0/(L*L)*s[66]*2.601E+3;
  ds[65] =  -s[65]*(DY*1.0/(L*L)*5.202E+3+y[64]*y[64])+s[64]*(B-y[64]*y[65]*2.0)+DY*1.0/(L*L)*s[63]*2.601E+3+DY*1.0/(L*L)*s[67]*2.601E+3;
  ds[66] =  -s[66]*(B-y[66]*y[67]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[67]*(y[66]*y[66])+DX*1.0/(L*L)*s[64]*2.601E+3+DX*1.0/(L*L)*s[68]*2.601E+3;
  ds[67] =  -s[67]*(DY*1.0/(L*L)*5.202E+3+y[66]*y[66])+s[66]*(B-y[66]*y[67]*2.0)+DY*1.0/(L*L)*s[65]*2.601E+3+DY*1.0/(L*L)*s[69]*2.601E+3;
  ds[68] =  -s[68]*(B-y[68]*y[69]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[69]*(y[68]*y[68])+DX*1.0/(L*L)*s[66]*2.601E+3+DX*1.0/(L*L)*s[70]*2.601E+3;
  ds[69] =  -s[69]*(DY*1.0/(L*L)*5.202E+3+y[68]*y[68])+s[68]*(B-y[68]*y[69]*2.0)+DY*1.0/(L*L)*s[67]*2.601E+3+DY*1.0/(L*L)*s[71]*2.601E+3;
  ds[70] =  -s[70]*(B-y[70]*y[71]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[71]*(y[70]*y[70])+DX*1.0/(L*L)*s[68]*2.601E+3+DX*1.0/(L*L)*s[72]*2.601E+3;
  ds[71] =  -s[71]*(DY*1.0/(L*L)*5.202E+3+y[70]*y[70])+s[70]*(B-y[70]*y[71]*2.0)+DY*1.0/(L*L)*s[69]*2.601E+3+DY*1.0/(L*L)*s[73]*2.601E+3;
  ds[72] =  -s[72]*(B-y[72]*y[73]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[73]*(y[72]*y[72])+DX*1.0/(L*L)*s[70]*2.601E+3+DX*1.0/(L*L)*s[74]*2.601E+3;
  ds[73] =  -s[73]*(DY*1.0/(L*L)*5.202E+3+y[72]*y[72])+s[72]*(B-y[72]*y[73]*2.0)+DY*1.0/(L*L)*s[71]*2.601E+3+DY*1.0/(L*L)*s[75]*2.601E+3;
  ds[74] =  -s[74]*(B-y[74]*y[75]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[75]*(y[74]*y[74])+DX*1.0/(L*L)*s[72]*2.601E+3+DX*1.0/(L*L)*s[76]*2.601E+3;
  ds[75] =  -s[75]*(DY*1.0/(L*L)*5.202E+3+y[74]*y[74])+s[74]*(B-y[74]*y[75]*2.0)+DY*1.0/(L*L)*s[73]*2.601E+3+DY*1.0/(L*L)*s[77]*2.601E+3;
  ds[76] =  -s[76]*(B-y[76]*y[77]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[77]*(y[76]*y[76])+DX*1.0/(L*L)*s[74]*2.601E+3+DX*1.0/(L*L)*s[78]*2.601E+3;
  ds[77] =  -s[77]*(DY*1.0/(L*L)*5.202E+3+y[76]*y[76])+s[76]*(B-y[76]*y[77]*2.0)+DY*1.0/(L*L)*s[75]*2.601E+3+DY*1.0/(L*L)*s[79]*2.601E+3;
  ds[78] =  -s[78]*(B-y[78]*y[79]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[79]*(y[78]*y[78])+DX*1.0/(L*L)*s[76]*2.601E+3+DX*1.0/(L*L)*s[80]*2.601E+3;
  ds[79] =  -s[79]*(DY*1.0/(L*L)*5.202E+3+y[78]*y[78])+s[78]*(B-y[78]*y[79]*2.0)+DY*1.0/(L*L)*s[77]*2.601E+3+DY*1.0/(L*L)*s[81]*2.601E+3;
  ds[80] =  -s[80]*(B-y[80]*y[81]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[81]*(y[80]*y[80])+DX*1.0/(L*L)*s[78]*2.601E+3+DX*1.0/(L*L)*s[82]*2.601E+3;
  ds[81] =  -s[81]*(DY*1.0/(L*L)*5.202E+3+y[80]*y[80])+s[80]*(B-y[80]*y[81]*2.0)+DY*1.0/(L*L)*s[79]*2.601E+3+DY*1.0/(L*L)*s[83]*2.601E+3;
  ds[82] =  -s[82]*(B-y[82]*y[83]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[83]*(y[82]*y[82])+DX*1.0/(L*L)*s[80]*2.601E+3+DX*1.0/(L*L)*s[84]*2.601E+3;
  ds[83] =  -s[83]*(DY*1.0/(L*L)*5.202E+3+y[82]*y[82])+s[82]*(B-y[82]*y[83]*2.0)+DY*1.0/(L*L)*s[81]*2.601E+3+DY*1.0/(L*L)*s[85]*2.601E+3;
  ds[84] =  -s[84]*(B-y[84]*y[85]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[85]*(y[84]*y[84])+DX*1.0/(L*L)*s[82]*2.601E+3+DX*1.0/(L*L)*s[86]*2.601E+3;
  ds[85] =  -s[85]*(DY*1.0/(L*L)*5.202E+3+y[84]*y[84])+s[84]*(B-y[84]*y[85]*2.0)+DY*1.0/(L*L)*s[83]*2.601E+3+DY*1.0/(L*L)*s[87]*2.601E+3;
  ds[86] =  -s[86]*(B-y[86]*y[87]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[87]*(y[86]*y[86])+DX*1.0/(L*L)*s[84]*2.601E+3+DX*1.0/(L*L)*s[88]*2.601E+3;
  ds[87] =  -s[87]*(DY*1.0/(L*L)*5.202E+3+y[86]*y[86])+s[86]*(B-y[86]*y[87]*2.0)+DY*1.0/(L*L)*s[85]*2.601E+3+DY*1.0/(L*L)*s[89]*2.601E+3;
  ds[88] =  -s[88]*(B-y[88]*y[89]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[89]*(y[88]*y[88])+DX*1.0/(L*L)*s[86]*2.601E+3+DX*1.0/(L*L)*s[90]*2.601E+3;
  ds[89] =  -s[89]*(DY*1.0/(L*L)*5.202E+3+y[88]*y[88])+s[88]*(B-y[88]*y[89]*2.0)+DY*1.0/(L*L)*s[87]*2.601E+3+DY*1.0/(L*L)*s[91]*2.601E+3;
  ds[90] =  -s[90]*(B-y[90]*y[91]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[91]*(y[90]*y[90])+DX*1.0/(L*L)*s[88]*2.601E+3+DX*1.0/(L*L)*s[92]*2.601E+3;
  ds[91] =  -s[91]*(DY*1.0/(L*L)*5.202E+3+y[90]*y[90])+s[90]*(B-y[90]*y[91]*2.0)+DY*1.0/(L*L)*s[89]*2.601E+3+DY*1.0/(L*L)*s[93]*2.601E+3;
  ds[92] =  -s[92]*(B-y[92]*y[93]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[93]*(y[92]*y[92])+DX*1.0/(L*L)*s[90]*2.601E+3+DX*1.0/(L*L)*s[94]*2.601E+3;
  ds[93] =  -s[93]*(DY*1.0/(L*L)*5.202E+3+y[92]*y[92])+s[92]*(B-y[92]*y[93]*2.0)+DY*1.0/(L*L)*s[91]*2.601E+3+DY*1.0/(L*L)*s[95]*2.601E+3;
  ds[94] =  -s[94]*(B-y[94]*y[95]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[95]*(y[94]*y[94])+DX*1.0/(L*L)*s[92]*2.601E+3+DX*1.0/(L*L)*s[96]*2.601E+3;
  ds[95] =  -s[95]*(DY*1.0/(L*L)*5.202E+3+y[94]*y[94])+s[94]*(B-y[94]*y[95]*2.0)+DY*1.0/(L*L)*s[93]*2.601E+3+DY*1.0/(L*L)*s[97]*2.601E+3;
  ds[96] =  -s[96]*(B-y[96]*y[97]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[97]*(y[96]*y[96])+DX*1.0/(L*L)*s[94]*2.601E+3+DX*1.0/(L*L)*s[98]*2.601E+3;
  ds[97] =  -s[97]*(DY*1.0/(L*L)*5.202E+3+y[96]*y[96])+s[96]*(B-y[96]*y[97]*2.0)+DY*1.0/(L*L)*s[95]*2.601E+3+DY*1.0/(L*L)*s[99]*2.601E+3;
  ds[98] =  -s[98]*(B-y[98]*y[99]*2.0+DX*1.0/(L*L)*5.202E+3+1.0)+s[99]*(y[98]*y[98])+DX*1.0/(L*L)*s[96]*2.601E+3;
  ds[99] =  -s[99]*(DY*1.0/(L*L)*5.202E+3+y[98]*y[98])+s[98]*(B-y[98]*y[99]*2.0)+DY*1.0/(L*L)*s[97]*2.601E+3;
  
  for ( i = 0; i < 100; i++ ) {
    if (fabs(dsdt[i]-ds[i]) > 1e-12) {
      mexErrMsgIdAndTxt("dsdt:wrong", "dsdt_%d is not correct, discrepancy: %.8e\n", i, fabs(dsdt[i]-ds[i]));
      
    }
  }
  */
  return 0;
}