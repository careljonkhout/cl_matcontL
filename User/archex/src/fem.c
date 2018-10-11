#include "fem.h"
#include "gaussquad.h"

#define MAXNSHAPE 9


void neohook(double* sigma, double* D, const double* F, double E, double nu)
{
    /* <generator matexpr>

    input E, nu;
    input F(2,2);
    output D(3,3);
    output sigma(3);

    b = F*F';
    I = eye(2);
    J = F(1,1)*F(2,2) - F(2,1)*F(1,2);

    lambda = E*nu/(1+nu)/(1-2*nu);
    mu     = E/2/(1+nu);

    sigmat = mu/J*(b-I) + lambda*(J-1)*I;
    sigma  = [sigmat(1,1); sigmat(2,2); sigmat(1,2)];

    a = mu/J - lambda*(J-1);
    b = lambda*(2*J-1);
    D = [2*a+b,     b, 0;
             b, 2*a+b, 0;
             0,     0, a];
    */
}


void remap_gradients(double* FF, double* dN, double* J, int nshape)
{
    int i;
    double invF[2*2];
    *J = (FF[0]*FF[3]-FF[1]*FF[2]);

    invF[0] =  FF[3]/ *J; invF[2] = -FF[2]/ *J;
    invF[1] = -FF[1]/ *J, invF[3] =  FF[0]/ *J;

    for (i = 0; i < nshape; ++i) {
        double* dNi = dN + i;

        /* <generator matexpr>
        inout dNi[nshape](1,2);
        input invF(2,2);
        dNi = dNi*invF;
        */
    }
}


void get_shapes(const double* nodex, const double *XX, double* xx,
                double* N, double* dN, double* J, int nshape)
{
    double X = XX[0];
    double Y = XX[1];
    double FF[2*2];

    if (nshape == 4) {
        /* <generator matexpr>
        input X, Y;
        input nodex(2,4);

        output N(4);
        output dN(4,2);
        output xx(2);
        output FF(2,2);

        N1x = (1-X)/2;  N1y = (1-Y)/2;
        N2x = (1+X)/2;  N2y = (1+Y)/2;    

        N  = [N1x*N1y; N2x*N1y; N2x*N2y; N1x*N2y];
        dN = deriv(N, [X,Y]);

        xx = nodex*N;
        FF = nodex*dN;
        */
    } else if (nshape == 9) {
        /* <generator matexpr>
        input X, Y;
        input nodex(2,9);

        output N(9);
        output dN(9,2);
        output xx(2);
        output FF(2,2);

        N1x = (X-1)*X/2;    N1y = (Y-1)*Y/2; 
        N2x = (1-X)*(1+X);  N2y = (1-Y)*(1+Y); 
        N3x = (X+1)*X/2;    N3y = (Y+1)*Y/2;

        N  = [N1x*N1y; N2x*N1y; N3x*N1y;
              N3x*N2y; N3x*N3y; N2x*N3y;
              N1x*N3y; N1x*N2y; N2x*N2y];
        dN = deriv(N, [X,Y]);

        xx = nodex*N;
        FF = nodex*dN;
        */
    }
    remap_gradients(FF, dN, J, nshape);
}


void evalF(const double* nodex, const double* nodeu, const double* XX,
           double* F, int nshape)
{
    double xx[2];
    double N [MAXNSHAPE];
    double dN[MAXNSHAPE*2];
    double J;
    int i;

    get_shapes(nodex, XX, xx, N, dN, &J, nshape);

    /* <generator matexpr>
    output F(2,2);
    F = eye(2);
    */

    for (i = 0; i < nshape; ++i) {
        const double* nodeui = nodeu + 2*i;
        double*       dNi    = dN + i;

        /* <generator matexpr>
        input nodeui(2,1);
        input dNi[nshape](1,2);
        inout F(2,2);
        F += nodeui * dNi;
        */
    }
}


void add_Kf(const double* nodex, const double* nodeu, 
            const double* XX, double wt, 
            const double* sigma, const double* DT,
            double* K, double* f, int nshape)
{
    int i, j;
    double nodexu[MAXNSHAPE*2];
    double J;
    double xx[2];
    double N [MAXNSHAPE];
    double dN[MAXNSHAPE*2];

    for (i = 0; i < nshape*2; ++i)
        nodexu[i] = nodex[i] + nodeu[i];

    get_shapes(nodexu, XX, xx, N, dN, &J, nshape);

    for (j = 0; j < nshape; ++j) {
        double* dNj = dN + j;
        double* fj  = f + 2*j;

        /* <generator matexpr>
        input J, wt;
        input dNj[nshape](1,2);
        input sigma(3);
        inout fj(2);

        Bj = [dNj(1), 0;
                   0, dNj(2);
              dNj(2), dNj(1)];

        fj += -(Bj'*sigma)*(J*wt);
        */

        for (i = 0; i < nshape; ++i) {
            double* dNi = dN + i;
            double* Kij = K + (4*nshape*j + 2*i);
            int nshape2 = nshape*2;

            /* <generator matexpr>
            input J, wt;
            input dNi[nshape](1,2);
            input dNj[nshape](1,2);
            input DT(3,3);
            input sigma(3);
            inout Kij[nshape2](2,2);

            Bi = [dNi(1), 0;
                       0, dNi(2);
                  dNi(2), dNi(1)];

            Bj = [dNj(1), 0;
                       0, dNj(2);
                  dNj(2), dNj(1)];

            sig  = [sigma(1), sigma(3); 
                    sigma(3), sigma(2)];

            KMij = Bi'*DT*Bj;
            Gij  = dNi*sig*dNj';
            Kij += (KMij + Gij*eye(2))*(J*wt);
            */
        }
    }
}


void clear_Kf(double* Kelt, double* felt, int nshape)
{
    int i;
    int nshape2 = nshape*2;
    for (i = 0; i < nshape2*nshape2; ++i)
        Kelt[i] = 0;
    for (i = 0; i < nshape2; ++i)
        felt[i] = 0;
}


void element_neohook(const double* nodex, const double* nodeu, 
                     double E, double nu, double* K, double* f,
                     int nshape)
{
    int ix, iy;
    int ngauss = (nshape == 4) ? 2 : 3;
    clear_Kf(K, f, nshape);

    for (ix = 0; ix < ngauss; ++ix) {
        for (iy = 0; iy < ngauss; ++iy) {
            double XX[2];
            double F[4];
            double D[3*3];
            double sigma[3];
            double wt;

            XX[0] = gauss_point(ix, ngauss);
            XX[1] = gauss_point(iy, ngauss);
            wt = gauss_weight(ix,ngauss)*gauss_weight(iy,ngauss);

            evalF(nodex, nodeu, XX, F, nshape);
            neohook(sigma, D, F, E, nu);
            add_Kf(nodex, nodeu, XX, wt, sigma, D, K, f, nshape);
        }
    }
}


void block2d4(double x1, double y1, double x2, double y2, int nx, int ny,
              int* ix, double* nodex)
{
    int jx, jy, k;
    int mx = nx+1;
    int my = ny+1;
    int ixt[2*2];

    ixt[0] = 0;
    ixt[1] = my;
    ixt[2] = my+1;
    ixt[3] = 1;

    for (jx = 0; jx < mx; ++jx) {
        double x = x1+jx*(x2-x1)/(mx-1);
        for (jy = 0; jy < my; ++jy) {
            double y = y1+jy*(y2-y1)/(my-1);
            *nodex++ = x;
            *nodex++ = y;
        }
    }

    for (jx = 0; jx < nx; ++jx)
        for (jy = 0; jy < ny; ++jy)
            for (k = 0; k < 4; ++k)
                *ix++ = (jx*my+jy+1) + ixt[k];
}


void block2d9(double x1, double y1, double x2, double y2, int nx, int ny,
              int* ix, double* nodex)
{
    int jx, jy, k;
    int mx = 2*nx+1;
    int my = 2*ny+1;
    int ixt[3*3];

    ixt[0] = 0;      ixt[1] = my;     ixt[2] = 2*my;
    ixt[3] = 2*my+1; ixt[4] = 2*my+2; ixt[5] = my+2;
    ixt[6] = 2;      ixt[7] = 1;      ixt[8] = my+1;

    for (jx = 0; jx < mx; ++jx) {
        double x = x1+jx*(x2-x1)/(mx-1);
        for (jy = 0; jy < my; ++jy) {
            double y = y1+jy*(y2-y1)/(my-1);
            *nodex++ = x;
            *nodex++ = y;
        }
    }

    for (jx = 0; jx < nx; ++jx)
        for (jy = 0; jy < ny; ++jy)
            for (k = 0; k < 9; ++k)
                *ix++ = (2*jx*my+2*jy+1) + ixt[k];
}
