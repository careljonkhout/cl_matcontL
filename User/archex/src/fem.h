#ifndef FEM_H
#define FEM_H

void element_neohook(const double* nodex, const double* nodeu, 
                     double E, double nu, double* KT, double* f,
                     int nshape);

void block2d(double x1, double y1, double x2, double y2, int nx, int ny,
             int* ix, double* nodex);

void block2d4(double x1, double y1, double x2, double y2, int nx, int ny,
              int* ix, double* nodex);

void block2d9(double x1, double y1, double x2, double y2, int nx, int ny,
              int* ix, double* nodex);

#endif /* FEM_H */
