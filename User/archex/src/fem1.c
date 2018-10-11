#line 1 "fem.c"
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
    /* <generated matexpr> */ {
#line 11 "fem.c"
    double tmp2_ = E;
#line 11 "fem.c"
    double tmp4_ = nu;
#line 12 "fem.c"
    double tmp6_ = F[0*2+0];
    double tmp7_ = F[0*2+1];
    double tmp8_ = F[1*2+0];
    double tmp9_ = F[1*2+1];
#line 16 "fem.c"
    double tmp11_ = tmp6_ * tmp6_;
    double tmp12_ = tmp8_ * tmp8_;
    double tmp13_ = tmp11_ + tmp12_;
    double tmp14_ = tmp6_ * tmp7_;
    double tmp15_ = tmp8_ * tmp9_;
    double tmp16_ = tmp14_ + tmp15_;
    double tmp17_ = tmp7_ * tmp7_;
    double tmp18_ = tmp9_ * tmp9_;
    double tmp19_ = tmp17_ + tmp18_;
#line 17 "fem.c"
    double tmp21_ = 1;
    double tmp22_ = 0;
#line 18 "fem.c"
    double tmp24_ = tmp6_ * tmp9_;
    double tmp25_ = tmp7_ * tmp8_;
    double tmp26_ = tmp24_ - tmp25_;
#line 20 "fem.c"
    double tmp28_ = tmp2_ * tmp4_;
    double tmp29_ = 1 + tmp4_;
    double tmp30_ = tmp28_ / tmp29_;
    double tmp31_ = 2.0;
    double tmp32_ = 2.0 * tmp4_;
    double tmp33_ = 1 - tmp32_;
    double tmp34_ = tmp30_ / tmp33_;
#line 21 "fem.c"
    double tmp36_ = tmp2_ / 2.0;
    double tmp37_ = tmp36_ / tmp29_;
#line 23 "fem.c"
    double tmp39_ = tmp37_ / tmp26_;
    double tmp40_ = tmp13_ - 1;
    double tmp41_ = tmp19_ - 1;
    double tmp42_ = tmp39_ * tmp40_;
    double tmp43_ = tmp39_ * tmp16_;
    double tmp44_ = tmp39_ * tmp41_;
    double tmp45_ = tmp26_ - 1;
    double tmp46_ = tmp34_ * tmp45_;
    double tmp47_ = tmp42_ + tmp46_;
    double tmp48_ = tmp44_ + tmp46_;
#line 24 "fem.c"
#line 26 "fem.c"
    double tmp51_ = tmp39_ - tmp46_;
#line 27 "fem.c"
    double tmp53_ = 2.0 * tmp26_;
    double tmp54_ = tmp53_ - 1;
    double tmp55_ = tmp34_ * tmp54_;
#line 30 "fem.c"
    double tmp57_ = 2.0 * tmp51_;
    double tmp58_ = tmp57_ + tmp55_;
#line 13 "fem.c"
    D[0*3+0] = tmp58_;
    D[0*3+1] = tmp55_;
    D[0*3+2] = 0;
    D[1*3+0] = tmp55_;
    D[1*3+1] = tmp58_;
    D[1*3+2] = 0;
    D[2*3+0] = 0;
    D[2*3+1] = 0;
    D[2*3+2] = tmp51_;
#line 14 "fem.c"
    sigma[0*3+0] = tmp47_;
    sigma[0*3+1] = tmp48_;
    sigma[0*3+2] = tmp43_;
    } /* </generated matexpr> */
#line 32 "fem.c"
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
        /* <generated matexpr> */ {
#line 48 "fem.c"
        double tmp2_ = dNi[0*nshape+0];
        double tmp3_ = dNi[1*nshape+0];
#line 49 "fem.c"
        double tmp5_ = invF[0*2+0];
        double tmp6_ = invF[0*2+1];
        double tmp7_ = invF[1*2+0];
        double tmp8_ = invF[1*2+1];
#line 50 "fem.c"
        double tmp10_ = tmp2_ * tmp5_;
        double tmp11_ = tmp3_ * tmp6_;
        double tmp12_ = tmp10_ + tmp11_;
        double tmp13_ = tmp2_ * tmp7_;
        double tmp14_ = tmp3_ * tmp8_;
        double tmp15_ = tmp13_ + tmp14_;
#line 48 "fem.c"
        dNi[0*nshape+0] = tmp12_;
        dNi[1*nshape+0] = tmp15_;
        } /* </generated matexpr> */
#line 52 "fem.c"
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
        /* <generated matexpr> */ {
#line 65 "fem.c"
        double tmp2_ = X;
#line 65 "fem.c"
        double tmp4_ = Y;
#line 66 "fem.c"
        double tmp6_ = nodex[0*2+0];
        double tmp7_ = nodex[0*2+1];
        double tmp8_ = nodex[1*2+0];
        double tmp9_ = nodex[1*2+1];
        double tmp10_ = nodex[2*2+0];
        double tmp11_ = nodex[2*2+1];
        double tmp12_ = nodex[3*2+0];
        double tmp13_ = nodex[3*2+1];
#line 73 "fem.c"
        double tmp15_ = 1.0;
        double tmp16_ = 1.0 - tmp2_;
        double tmp17_ = 2.0;
        double tmp18_ = tmp16_ / 2.0;
#line 73 "fem.c"
        double tmp20_ = 1.0 - tmp4_;
        double tmp21_ = tmp20_ / 2.0;
#line 74 "fem.c"
        double tmp23_ = 1.0 + tmp2_;
        double tmp24_ = tmp23_ / 2.0;
#line 74 "fem.c"
        double tmp26_ = 1.0 + tmp4_;
        double tmp27_ = tmp26_ / 2.0;
#line 76 "fem.c"
        double tmp29_ = tmp18_ * tmp21_;
        double tmp30_ = tmp24_ * tmp21_;
        double tmp31_ = tmp24_ * tmp27_;
        double tmp32_ = tmp18_ * tmp27_;
#line 77 "fem.c"
        double tmp34_ = 0;
        double tmp35_ = 0 - 1.0;
        double tmp36_ = 2.0 * tmp35_;
        double tmp37_ = 2.0 * 2.0;
        double tmp38_ = tmp36_ / tmp37_;
        double tmp39_ = tmp21_ * tmp38_;
        double tmp40_ = 2.0 / tmp37_;
        double tmp41_ = tmp21_ * tmp40_;
        double tmp42_ = tmp27_ * tmp40_;
        double tmp43_ = tmp27_ * tmp38_;
        double tmp44_ = tmp18_ * tmp38_;
        double tmp45_ = tmp24_ * tmp38_;
        double tmp46_ = tmp24_ * tmp40_;
        double tmp47_ = tmp18_ * tmp40_;
#line 79 "fem.c"
        double tmp49_ = tmp6_ * tmp29_;
        double tmp50_ = tmp8_ * tmp30_;
        double tmp51_ = tmp49_ + tmp50_;
        double tmp52_ = tmp10_ * tmp31_;
        double tmp53_ = tmp51_ + tmp52_;
        double tmp54_ = tmp12_ * tmp32_;
        double tmp55_ = tmp53_ + tmp54_;
        double tmp56_ = tmp7_ * tmp29_;
        double tmp57_ = tmp9_ * tmp30_;
        double tmp58_ = tmp56_ + tmp57_;
        double tmp59_ = tmp11_ * tmp31_;
        double tmp60_ = tmp58_ + tmp59_;
        double tmp61_ = tmp13_ * tmp32_;
        double tmp62_ = tmp60_ + tmp61_;
#line 80 "fem.c"
        double tmp64_ = tmp6_ * tmp39_;
        double tmp65_ = tmp8_ * tmp41_;
        double tmp66_ = tmp64_ + tmp65_;
        double tmp67_ = tmp10_ * tmp42_;
        double tmp68_ = tmp66_ + tmp67_;
        double tmp69_ = tmp12_ * tmp43_;
        double tmp70_ = tmp68_ + tmp69_;
        double tmp71_ = tmp7_ * tmp39_;
        double tmp72_ = tmp9_ * tmp41_;
        double tmp73_ = tmp71_ + tmp72_;
        double tmp74_ = tmp11_ * tmp42_;
        double tmp75_ = tmp73_ + tmp74_;
        double tmp76_ = tmp13_ * tmp43_;
        double tmp77_ = tmp75_ + tmp76_;
        double tmp78_ = tmp6_ * tmp44_;
        double tmp79_ = tmp8_ * tmp45_;
        double tmp80_ = tmp78_ + tmp79_;
        double tmp81_ = tmp10_ * tmp46_;
        double tmp82_ = tmp80_ + tmp81_;
        double tmp83_ = tmp12_ * tmp47_;
        double tmp84_ = tmp82_ + tmp83_;
        double tmp85_ = tmp7_ * tmp44_;
        double tmp86_ = tmp9_ * tmp45_;
        double tmp87_ = tmp85_ + tmp86_;
        double tmp88_ = tmp11_ * tmp46_;
        double tmp89_ = tmp87_ + tmp88_;
        double tmp90_ = tmp13_ * tmp47_;
        double tmp91_ = tmp89_ + tmp90_;
#line 68 "fem.c"
        N[0*4+0] = tmp29_;
        N[0*4+1] = tmp30_;
        N[0*4+2] = tmp31_;
        N[0*4+3] = tmp32_;
#line 69 "fem.c"
        dN[0*4+0] = tmp39_;
        dN[0*4+1] = tmp41_;
        dN[0*4+2] = tmp42_;
        dN[0*4+3] = tmp43_;
        dN[1*4+0] = tmp44_;
        dN[1*4+1] = tmp45_;
        dN[1*4+2] = tmp46_;
        dN[1*4+3] = tmp47_;
#line 70 "fem.c"
        xx[0*2+0] = tmp55_;
        xx[0*2+1] = tmp62_;
#line 71 "fem.c"
        FF[0*2+0] = tmp70_;
        FF[0*2+1] = tmp77_;
        FF[1*2+0] = tmp84_;
        FF[1*2+1] = tmp91_;
        } /* </generated matexpr> */
#line 82 "fem.c"
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
        /* <generated matexpr> */ {
#line 84 "fem.c"
        double tmp2_ = X;
#line 84 "fem.c"
        double tmp4_ = Y;
#line 85 "fem.c"
        double tmp6_ = nodex[0*2+0];
        double tmp7_ = nodex[0*2+1];
        double tmp8_ = nodex[1*2+0];
        double tmp9_ = nodex[1*2+1];
        double tmp10_ = nodex[2*2+0];
        double tmp11_ = nodex[2*2+1];
        double tmp12_ = nodex[3*2+0];
        double tmp13_ = nodex[3*2+1];
        double tmp14_ = nodex[4*2+0];
        double tmp15_ = nodex[4*2+1];
        double tmp16_ = nodex[5*2+0];
        double tmp17_ = nodex[5*2+1];
        double tmp18_ = nodex[6*2+0];
        double tmp19_ = nodex[6*2+1];
        double tmp20_ = nodex[7*2+0];
        double tmp21_ = nodex[7*2+1];
        double tmp22_ = nodex[8*2+0];
        double tmp23_ = nodex[8*2+1];
#line 92 "fem.c"
        double tmp25_ = 1.0;
        double tmp26_ = tmp2_ - 1.0;
        double tmp27_ = tmp26_ * tmp2_;
        double tmp28_ = 2.0;
        double tmp29_ = tmp27_ / 2.0;
#line 92 "fem.c"
        double tmp31_ = tmp4_ - 1.0;
        double tmp32_ = tmp31_ * tmp4_;
        double tmp33_ = tmp32_ / 2.0;
#line 93 "fem.c"
        double tmp35_ = 1.0 - tmp2_;
        double tmp36_ = 1.0 + tmp2_;
        double tmp37_ = tmp35_ * tmp36_;
#line 93 "fem.c"
        double tmp39_ = 1.0 - tmp4_;
        double tmp40_ = 1.0 + tmp4_;
        double tmp41_ = tmp39_ * tmp40_;
#line 94 "fem.c"
        double tmp43_ = tmp2_ + 1.0;
        double tmp44_ = tmp43_ * tmp2_;
        double tmp45_ = tmp44_ / 2.0;
#line 94 "fem.c"
        double tmp47_ = tmp4_ + 1.0;
        double tmp48_ = tmp47_ * tmp4_;
        double tmp49_ = tmp48_ / 2.0;
#line 98 "fem.c"
        double tmp51_ = tmp29_ * tmp33_;
        double tmp52_ = tmp37_ * tmp33_;
        double tmp53_ = tmp45_ * tmp33_;
        double tmp54_ = tmp45_ * tmp41_;
        double tmp55_ = tmp45_ * tmp49_;
        double tmp56_ = tmp37_ * tmp49_;
        double tmp57_ = tmp29_ * tmp49_;
        double tmp58_ = tmp29_ * tmp41_;
        double tmp59_ = tmp37_ * tmp41_;
#line 99 "fem.c"
        double tmp61_ = 0;
        double tmp62_ = tmp26_ + tmp2_;
        double tmp63_ = 2.0 * tmp62_;
        double tmp64_ = 2.0 * 2.0;
        double tmp65_ = tmp63_ / tmp64_;
        double tmp66_ = tmp33_ * tmp65_;
        double tmp67_ = 0 - 1.0;
        double tmp68_ = tmp36_ * tmp67_;
        double tmp69_ = tmp35_ + tmp68_;
        double tmp70_ = tmp33_ * tmp69_;
        double tmp71_ = tmp43_ + tmp2_;
        double tmp72_ = 2.0 * tmp71_;
        double tmp73_ = tmp72_ / tmp64_;
        double tmp74_ = tmp33_ * tmp73_;
        double tmp75_ = tmp41_ * tmp73_;
        double tmp76_ = tmp49_ * tmp73_;
        double tmp77_ = tmp49_ * tmp69_;
        double tmp78_ = tmp49_ * tmp65_;
        double tmp79_ = tmp41_ * tmp65_;
        double tmp80_ = tmp41_ * tmp69_;
        double tmp81_ = tmp31_ + tmp4_;
        double tmp82_ = 2.0 * tmp81_;
        double tmp83_ = tmp82_ / tmp64_;
        double tmp84_ = tmp29_ * tmp83_;
        double tmp85_ = tmp37_ * tmp83_;
        double tmp86_ = tmp45_ * tmp83_;
        double tmp87_ = tmp40_ * tmp67_;
        double tmp88_ = tmp39_ + tmp87_;
        double tmp89_ = tmp45_ * tmp88_;
        double tmp90_ = tmp47_ + tmp4_;
        double tmp91_ = 2.0 * tmp90_;
        double tmp92_ = tmp91_ / tmp64_;
        double tmp93_ = tmp45_ * tmp92_;
        double tmp94_ = tmp37_ * tmp92_;
        double tmp95_ = tmp29_ * tmp92_;
        double tmp96_ = tmp29_ * tmp88_;
        double tmp97_ = tmp37_ * tmp88_;
#line 101 "fem.c"
        double tmp99_ = tmp6_ * tmp51_;
        double tmp100_ = tmp8_ * tmp52_;
        double tmp101_ = tmp99_ + tmp100_;
        double tmp102_ = tmp10_ * tmp53_;
        double tmp103_ = tmp101_ + tmp102_;
        double tmp104_ = tmp12_ * tmp54_;
        double tmp105_ = tmp103_ + tmp104_;
        double tmp106_ = tmp14_ * tmp55_;
        double tmp107_ = tmp105_ + tmp106_;
        double tmp108_ = tmp16_ * tmp56_;
        double tmp109_ = tmp107_ + tmp108_;
        double tmp110_ = tmp18_ * tmp57_;
        double tmp111_ = tmp109_ + tmp110_;
        double tmp112_ = tmp20_ * tmp58_;
        double tmp113_ = tmp111_ + tmp112_;
        double tmp114_ = tmp22_ * tmp59_;
        double tmp115_ = tmp113_ + tmp114_;
        double tmp116_ = tmp7_ * tmp51_;
        double tmp117_ = tmp9_ * tmp52_;
        double tmp118_ = tmp116_ + tmp117_;
        double tmp119_ = tmp11_ * tmp53_;
        double tmp120_ = tmp118_ + tmp119_;
        double tmp121_ = tmp13_ * tmp54_;
        double tmp122_ = tmp120_ + tmp121_;
        double tmp123_ = tmp15_ * tmp55_;
        double tmp124_ = tmp122_ + tmp123_;
        double tmp125_ = tmp17_ * tmp56_;
        double tmp126_ = tmp124_ + tmp125_;
        double tmp127_ = tmp19_ * tmp57_;
        double tmp128_ = tmp126_ + tmp127_;
        double tmp129_ = tmp21_ * tmp58_;
        double tmp130_ = tmp128_ + tmp129_;
        double tmp131_ = tmp23_ * tmp59_;
        double tmp132_ = tmp130_ + tmp131_;
#line 102 "fem.c"
        double tmp134_ = tmp6_ * tmp66_;
        double tmp135_ = tmp8_ * tmp70_;
        double tmp136_ = tmp134_ + tmp135_;
        double tmp137_ = tmp10_ * tmp74_;
        double tmp138_ = tmp136_ + tmp137_;
        double tmp139_ = tmp12_ * tmp75_;
        double tmp140_ = tmp138_ + tmp139_;
        double tmp141_ = tmp14_ * tmp76_;
        double tmp142_ = tmp140_ + tmp141_;
        double tmp143_ = tmp16_ * tmp77_;
        double tmp144_ = tmp142_ + tmp143_;
        double tmp145_ = tmp18_ * tmp78_;
        double tmp146_ = tmp144_ + tmp145_;
        double tmp147_ = tmp20_ * tmp79_;
        double tmp148_ = tmp146_ + tmp147_;
        double tmp149_ = tmp22_ * tmp80_;
        double tmp150_ = tmp148_ + tmp149_;
        double tmp151_ = tmp7_ * tmp66_;
        double tmp152_ = tmp9_ * tmp70_;
        double tmp153_ = tmp151_ + tmp152_;
        double tmp154_ = tmp11_ * tmp74_;
        double tmp155_ = tmp153_ + tmp154_;
        double tmp156_ = tmp13_ * tmp75_;
        double tmp157_ = tmp155_ + tmp156_;
        double tmp158_ = tmp15_ * tmp76_;
        double tmp159_ = tmp157_ + tmp158_;
        double tmp160_ = tmp17_ * tmp77_;
        double tmp161_ = tmp159_ + tmp160_;
        double tmp162_ = tmp19_ * tmp78_;
        double tmp163_ = tmp161_ + tmp162_;
        double tmp164_ = tmp21_ * tmp79_;
        double tmp165_ = tmp163_ + tmp164_;
        double tmp166_ = tmp23_ * tmp80_;
        double tmp167_ = tmp165_ + tmp166_;
        double tmp168_ = tmp6_ * tmp84_;
        double tmp169_ = tmp8_ * tmp85_;
        double tmp170_ = tmp168_ + tmp169_;
        double tmp171_ = tmp10_ * tmp86_;
        double tmp172_ = tmp170_ + tmp171_;
        double tmp173_ = tmp12_ * tmp89_;
        double tmp174_ = tmp172_ + tmp173_;
        double tmp175_ = tmp14_ * tmp93_;
        double tmp176_ = tmp174_ + tmp175_;
        double tmp177_ = tmp16_ * tmp94_;
        double tmp178_ = tmp176_ + tmp177_;
        double tmp179_ = tmp18_ * tmp95_;
        double tmp180_ = tmp178_ + tmp179_;
        double tmp181_ = tmp20_ * tmp96_;
        double tmp182_ = tmp180_ + tmp181_;
        double tmp183_ = tmp22_ * tmp97_;
        double tmp184_ = tmp182_ + tmp183_;
        double tmp185_ = tmp7_ * tmp84_;
        double tmp186_ = tmp9_ * tmp85_;
        double tmp187_ = tmp185_ + tmp186_;
        double tmp188_ = tmp11_ * tmp86_;
        double tmp189_ = tmp187_ + tmp188_;
        double tmp190_ = tmp13_ * tmp89_;
        double tmp191_ = tmp189_ + tmp190_;
        double tmp192_ = tmp15_ * tmp93_;
        double tmp193_ = tmp191_ + tmp192_;
        double tmp194_ = tmp17_ * tmp94_;
        double tmp195_ = tmp193_ + tmp194_;
        double tmp196_ = tmp19_ * tmp95_;
        double tmp197_ = tmp195_ + tmp196_;
        double tmp198_ = tmp21_ * tmp96_;
        double tmp199_ = tmp197_ + tmp198_;
        double tmp200_ = tmp23_ * tmp97_;
        double tmp201_ = tmp199_ + tmp200_;
#line 87 "fem.c"
        N[0*9+0] = tmp51_;
        N[0*9+1] = tmp52_;
        N[0*9+2] = tmp53_;
        N[0*9+3] = tmp54_;
        N[0*9+4] = tmp55_;
        N[0*9+5] = tmp56_;
        N[0*9+6] = tmp57_;
        N[0*9+7] = tmp58_;
        N[0*9+8] = tmp59_;
#line 88 "fem.c"
        dN[0*9+0] = tmp66_;
        dN[0*9+1] = tmp70_;
        dN[0*9+2] = tmp74_;
        dN[0*9+3] = tmp75_;
        dN[0*9+4] = tmp76_;
        dN[0*9+5] = tmp77_;
        dN[0*9+6] = tmp78_;
        dN[0*9+7] = tmp79_;
        dN[0*9+8] = tmp80_;
        dN[1*9+0] = tmp84_;
        dN[1*9+1] = tmp85_;
        dN[1*9+2] = tmp86_;
        dN[1*9+3] = tmp89_;
        dN[1*9+4] = tmp93_;
        dN[1*9+5] = tmp94_;
        dN[1*9+6] = tmp95_;
        dN[1*9+7] = tmp96_;
        dN[1*9+8] = tmp97_;
#line 89 "fem.c"
        xx[0*2+0] = tmp115_;
        xx[0*2+1] = tmp132_;
#line 90 "fem.c"
        FF[0*2+0] = tmp150_;
        FF[0*2+1] = tmp167_;
        FF[1*2+0] = tmp184_;
        FF[1*2+1] = tmp201_;
        } /* </generated matexpr> */
#line 104 "fem.c"
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
    /* <generated matexpr> */ {
#line 122 "fem.c"
    double tmp2_ = 1;
    double tmp3_ = 0;
#line 121 "fem.c"
    F[0*2+0] = 1;
    F[0*2+1] = 0;
    F[1*2+0] = 0;
    F[1*2+1] = 1;
    } /* </generated matexpr> */
#line 124 "fem.c"

    for (i = 0; i < nshape; ++i) {
        const double* nodeui = nodeu + 2*i;
        double*       dNi    = dN + i;

        /* <generator matexpr>
        input nodeui(2,1);
        input dNi[nshape](1,2);
        inout F(2,2);
        F += nodeui * dNi;
        */
        /* <generated matexpr> */ {
#line 130 "fem.c"
        double tmp2_ = nodeui[0*2+0];
        double tmp3_ = nodeui[0*2+1];
#line 131 "fem.c"
        double tmp5_ = dNi[0*nshape+0];
        double tmp6_ = dNi[1*nshape+0];
#line 132 "fem.c"
        double tmp8_ = F[0*2+0];
        double tmp9_ = F[0*2+1];
        double tmp10_ = F[1*2+0];
        double tmp11_ = F[1*2+1];
#line 133 "fem.c"
        double tmp13_ = tmp2_ * tmp5_;
        double tmp14_ = tmp3_ * tmp5_;
        double tmp15_ = tmp2_ * tmp6_;
        double tmp16_ = tmp3_ * tmp6_;
        double tmp17_ = tmp8_ + tmp13_;
        double tmp18_ = tmp9_ + tmp14_;
        double tmp19_ = tmp10_ + tmp15_;
        double tmp20_ = tmp11_ + tmp16_;
#line 132 "fem.c"
        F[0*2+0] = tmp17_;
        F[0*2+1] = tmp18_;
        F[1*2+0] = tmp19_;
        F[1*2+1] = tmp20_;
        } /* </generated matexpr> */
#line 135 "fem.c"
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
        /* <generated matexpr> */ {
#line 161 "fem.c"
        double tmp2_ = J;
#line 161 "fem.c"
        double tmp4_ = wt;
#line 162 "fem.c"
        double tmp6_ = dNj[0*nshape+0];
        double tmp7_ = dNj[1*nshape+0];
#line 163 "fem.c"
        double tmp9_ = sigma[0*3+0];
        double tmp10_ = sigma[0*3+1];
        double tmp11_ = sigma[0*3+2];
#line 164 "fem.c"
        double tmp13_ = fj[0*2+0];
        double tmp14_ = fj[0*2+1];
#line 168 "fem.c"
#line 170 "fem.c"
        double tmp17_ = tmp6_ * tmp9_;
        double tmp18_ = tmp7_ * tmp11_;
        double tmp19_ = tmp17_ + tmp18_;
        double tmp20_ = tmp7_ * tmp10_;
        double tmp21_ = tmp6_ * tmp11_;
        double tmp22_ = tmp20_ + tmp21_;
        double tmp23_ = -tmp19_;
        double tmp24_ = -tmp22_;
        double tmp25_ = tmp2_ * tmp4_;
        double tmp26_ = tmp23_ * tmp25_;
        double tmp27_ = tmp24_ * tmp25_;
        double tmp28_ = tmp13_ + tmp26_;
        double tmp29_ = tmp14_ + tmp27_;
#line 164 "fem.c"
        fj[0*2+0] = tmp28_;
        fj[0*2+1] = tmp29_;
        } /* </generated matexpr> */
#line 172 "fem.c"

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
            /* <generated matexpr> */ {
#line 179 "fem.c"
            double tmp2_ = J;
#line 179 "fem.c"
            double tmp4_ = wt;
#line 180 "fem.c"
            double tmp6_ = dNi[0*nshape+0];
            double tmp7_ = dNi[1*nshape+0];
#line 181 "fem.c"
            double tmp9_ = dNj[0*nshape+0];
            double tmp10_ = dNj[1*nshape+0];
#line 182 "fem.c"
            double tmp12_ = DT[0*3+0];
            double tmp13_ = DT[0*3+1];
            double tmp14_ = DT[0*3+2];
            double tmp15_ = DT[1*3+0];
            double tmp16_ = DT[1*3+1];
            double tmp17_ = DT[1*3+2];
            double tmp18_ = DT[2*3+0];
            double tmp19_ = DT[2*3+1];
            double tmp20_ = DT[2*3+2];
#line 183 "fem.c"
            double tmp22_ = sigma[0*3+0];
            double tmp23_ = sigma[0*3+1];
            double tmp24_ = sigma[0*3+2];
#line 184 "fem.c"
            double tmp26_ = Kij[0*nshape2+0];
            double tmp27_ = Kij[0*nshape2+1];
            double tmp28_ = Kij[1*nshape2+0];
            double tmp29_ = Kij[1*nshape2+1];
#line 188 "fem.c"
#line 192 "fem.c"
#line 195 "fem.c"
#line 197 "fem.c"
            double tmp34_ = tmp6_ * tmp12_;
            double tmp35_ = tmp7_ * tmp14_;
            double tmp36_ = tmp34_ + tmp35_;
            double tmp37_ = tmp7_ * tmp13_;
            double tmp38_ = tmp6_ * tmp14_;
            double tmp39_ = tmp37_ + tmp38_;
            double tmp40_ = tmp6_ * tmp15_;
            double tmp41_ = tmp7_ * tmp17_;
            double tmp42_ = tmp40_ + tmp41_;
            double tmp43_ = tmp7_ * tmp16_;
            double tmp44_ = tmp6_ * tmp17_;
            double tmp45_ = tmp43_ + tmp44_;
            double tmp46_ = tmp6_ * tmp18_;
            double tmp47_ = tmp7_ * tmp20_;
            double tmp48_ = tmp46_ + tmp47_;
            double tmp49_ = tmp7_ * tmp19_;
            double tmp50_ = tmp6_ * tmp20_;
            double tmp51_ = tmp49_ + tmp50_;
            double tmp52_ = tmp36_ * tmp9_;
            double tmp53_ = tmp48_ * tmp10_;
            double tmp54_ = tmp52_ + tmp53_;
            double tmp55_ = tmp39_ * tmp9_;
            double tmp56_ = tmp51_ * tmp10_;
            double tmp57_ = tmp55_ + tmp56_;
            double tmp58_ = tmp42_ * tmp10_;
            double tmp59_ = tmp48_ * tmp9_;
            double tmp60_ = tmp58_ + tmp59_;
            double tmp61_ = tmp45_ * tmp10_;
            double tmp62_ = tmp51_ * tmp9_;
            double tmp63_ = tmp61_ + tmp62_;
#line 198 "fem.c"
            double tmp65_ = tmp6_ * tmp22_;
            double tmp66_ = tmp7_ * tmp24_;
            double tmp67_ = tmp65_ + tmp66_;
            double tmp68_ = tmp6_ * tmp24_;
            double tmp69_ = tmp7_ * tmp23_;
            double tmp70_ = tmp68_ + tmp69_;
            double tmp71_ = tmp67_ * tmp9_;
            double tmp72_ = tmp70_ * tmp10_;
            double tmp73_ = tmp71_ + tmp72_;
#line 199 "fem.c"
            double tmp75_ = tmp54_ + tmp73_;
            double tmp76_ = tmp63_ + tmp73_;
            double tmp77_ = tmp2_ * tmp4_;
            double tmp78_ = tmp75_ * tmp77_;
            double tmp79_ = tmp57_ * tmp77_;
            double tmp80_ = tmp60_ * tmp77_;
            double tmp81_ = tmp76_ * tmp77_;
            double tmp82_ = tmp26_ + tmp78_;
            double tmp83_ = tmp27_ + tmp79_;
            double tmp84_ = tmp28_ + tmp80_;
            double tmp85_ = tmp29_ + tmp81_;
#line 184 "fem.c"
            Kij[0*nshape2+0] = tmp82_;
            Kij[0*nshape2+1] = tmp83_;
            Kij[1*nshape2+0] = tmp84_;
            Kij[1*nshape2+1] = tmp85_;
            } /* </generated matexpr> */
#line 201 "fem.c"
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
