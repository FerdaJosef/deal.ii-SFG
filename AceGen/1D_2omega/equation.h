#pragma once

#include <vector>
#include <cmath>

#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

using namespace dealii;

/*************************************************************
* AceGen    9.105 Windows (19 Oct 25)                        *
*           Co. J. Korelc  2020           12 Apr 26 22:54:17 *
**************************************************************
User     : Limited evaluation version
Notebook : DoubleDitch2D
Evaluation time                 : 1 s     Mode  : Optimal
Number of formulae              : 36      Method: Automatic
Subroutine                      : RandomEquation size: 537
Total size of Mathematica  code : 537 subexpressions
Total size of C code            : 1872 bytes */

/******************* S U B R O U T I N E *********************/
template <int dim, int n>
inline void equation(
    std::vector<double> &v,
    const Vector<double> &U,
    const Vector<double> &U0,
    const std::vector<Tensor<1,dim>> &GradU,
    Vector<double> &dPsiDu,
    std::vector<Tensor<1,dim>> &dPsidGradU,
    FullMatrix<double> &dPsiDu2,
    std::vector<std::vector<Tensor<1,dim>>> &dPsidUdGradU,
    std::vector<std::vector<Tensor<2,dim>>> &dPsidGradU2,
    double (*dt))
{
int i01;
v[7]=U[0];
v[38]=1e0-v[7];
v[41]=(v[38]*v[38]);
v[36]=2e0*v[7];
v[8]=U[1];
v[43]=1e0-v[8];
v[46]=(v[43]*v[43]);
v[37]=2e0*v[8];
v[9]=U0[0];
v[10]=U0[1];
v[11]=GradU[0][0];
v[12]=GradU[1][0];
v[50]=0.125e1/(*dt);
v[25]=v[43]-v[7];
v[39]=1e0-v[25];
v[40]=(v[39]*v[39]);
v[35]=-2e0*v[25];
v[70]=-0.6000000000000005e-2*v[35]*v[39];
v[48]=3e0+v[35];
v[31]=(v[25]*v[25]);
v[69]=-0.6000000000000005e-2*v[31];
v[88]=0.8e-1*v[35]-0.6000000000000005e-2*v[40]+0.4e-1*v[48]+v[69]+2e0*v[70];
v[89]=0.3e-1*v[50]+v[88];
v[78]=100e0*v[36]*v[37]-0.15000000000000013e-2*v[50]+v[88];
v[87]=0.4e-1*v[31]+v[35]*(-0.30000000000000027e-2*v[40]+0.2e-1*v[48])+v[39]*v[69];
v[55]=0.15000000000000013e-2*(-1e0+v[10]+v[25]+v[9]);
v[60]=-0.15000000000000013e-2*(v[11]+v[12]);
v[28]=(v[7]*v[7]);
v[80]=0.126e0*v[28];
v[29]=(v[8]*v[8]);
v[76]=0.126e0*v[29];
dPsiDu[0]=v[36]*(100e0*v[29]+0.63e-1*v[41])-v[38]*v[80]+v[87]+v[50]*(v[55]+0.315e-1*(v[7]-v[9]));
dPsiDu[1]=v[37]*(100e0*v[28]+0.63e-1*v[46])-v[43]*v[76]+v[50]*(v[55]+0.315e-1*(-v[10]+v[8]))+v[87];
dPsidGradU[0][0]=0.315e-1*v[11]+v[60];
dPsidGradU[1][0]=0.315e-1*v[12]+v[60];
dPsiDu2[0][0]=200e0*v[29]-0.252e0*v[36]*v[38]+0.126e0*v[41]+v[80]+v[89];
dPsiDu2[0][1]=v[78];
dPsiDu2[1][0]=v[78];
dPsiDu2[1][1]=200e0*v[28]-0.252e0*v[37]*v[43]+0.126e0*v[46]+v[76]+v[89];
dPsidUdGradU[0][0][0]=0e0;
dPsidUdGradU[0][1][0]=0e0;
dPsidUdGradU[1][0][0]=0e0;
dPsidUdGradU[1][1][0]=0e0;
dPsidGradU2[0][0][0][0]=0.3e-1;
dPsidGradU2[0][1][0][0]=-0.15000000000000013e-2;
dPsidGradU2[1][0][0][0]=-0.15000000000000013e-2;
dPsidGradU2[1][1][0][0]=0.3e-1;
};
