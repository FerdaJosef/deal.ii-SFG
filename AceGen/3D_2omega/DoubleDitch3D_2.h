#pragma once

#include <vector>
#include <cmath>

#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

using namespace dealii;

/*************************************************************
* AceGen    9.105 Windows (19 Oct 25)                        *
*           Co. J. Korelc  2020           9 Apr 26 22:08:32  *
**************************************************************
User     : Limited evaluation version
Notebook : DoubleDitch2D
Evaluation time                 : 1 s     Mode  : Optimal
Number of formulae              : 42      Method: Automatic
Subroutine                      : RandomEquation size: 995
Total size of Mathematica  code : 995 subexpressions
Total size of C code            : 3444 bytes */

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
v[5]=GradU[0][0];
v[6]=GradU[0][1];
v[7]=GradU[0][2];
v[8]=GradU[1][0];
v[9]=GradU[1][1];
v[10]=GradU[1][2];
v[11]=U[0];
v[48]=1e0-v[11];
v[51]=(v[48]*v[48]);
v[46]=2e0*v[11];
v[12]=U[1];
v[53]=1e0-v[12];
v[56]=(v[53]*v[53]);
v[47]=2e0*v[12];
v[13]=U0[0];
v[14]=U0[1];
v[60]=0.125e1/(*dt);
v[33]=-v[11]+v[53];
v[49]=1e0-v[33];
v[50]=(v[49]*v[49]);
v[45]=-2e0*v[33];
v[86]=-0.6000000000000005e-2*v[45]*v[49];
v[58]=3e0+v[45];
v[41]=(v[33]*v[33]);
v[85]=-0.6000000000000005e-2*v[41];
v[104]=0.8e-1*v[45]-0.6000000000000005e-2*v[50]+0.4e-1*v[58]+v[85]+2e0*v[86];
v[105]=v[104]+0.3e-1*v[60];
v[94]=v[104]+100e0*v[46]*v[47]-0.15000000000000013e-2*v[60];
v[103]=0.4e-1*v[41]+v[45]*(-0.30000000000000027e-2*v[50]+0.2e-1*v[58])+v[49]*v[85];
v[65]=0.15000000000000013e-2*(-1e0+v[13]+v[14]+v[33]);
v[72]=-0.15000000000000013e-2*(v[5]+v[8]);
v[74]=-0.15000000000000013e-2*(v[6]+v[9]);
v[76]=-0.15000000000000013e-2*(v[10]+v[7]);
v[38]=(v[11]*v[11]);
v[96]=0.126e0*v[38];
v[39]=(v[12]*v[12]);
v[92]=0.126e0*v[39];
dPsiDu[0]=v[103]+v[46]*(100e0*v[39]+0.63e-1*v[51])+v[60]*(0.315e-1*(v[11]-v[13])+v[65])
 -v[48]*v[96];
dPsiDu[1]=v[103]+v[47]*(100e0*v[38]+0.63e-1*v[56])+v[60]*(0.315e-1*(v[12]-v[14])+v[65])
 -v[53]*v[92];
dPsidGradU[0][0]=0.315e-1*v[5]+v[72];
dPsidGradU[0][1]=0.315e-1*v[6]+v[74];
dPsidGradU[0][2]=0.315e-1*v[7]+v[76];
dPsidGradU[1][0]=v[72]+0.315e-1*v[8];
dPsidGradU[1][1]=v[74]+0.315e-1*v[9];
dPsidGradU[1][2]=0.315e-1*v[10]+v[76];
dPsiDu2[0][0]=v[105]+200e0*v[39]-0.252e0*v[46]*v[48]+0.126e0*v[51]+v[96];
dPsiDu2[0][1]=v[94];
dPsiDu2[1][0]=v[94];
dPsiDu2[1][1]=v[105]+200e0*v[38]-0.252e0*v[47]*v[53]+0.126e0*v[56]+v[92];
dPsidUdGradU[0][0][0]=0e0;
dPsidUdGradU[0][0][1]=0e0;
dPsidUdGradU[0][0][2]=0e0;
dPsidUdGradU[0][1][0]=0e0;
dPsidUdGradU[0][1][1]=0e0;
dPsidUdGradU[0][1][2]=0e0;
dPsidUdGradU[1][0][0]=0e0;
dPsidUdGradU[1][0][1]=0e0;
dPsidUdGradU[1][0][2]=0e0;
dPsidUdGradU[1][1][0]=0e0;
dPsidUdGradU[1][1][1]=0e0;
dPsidUdGradU[1][1][2]=0e0;
dPsidGradU2[0][0][0][0]=0.3e-1;
dPsidGradU2[0][0][0][1]=0e0;
dPsidGradU2[0][0][0][2]=0e0;
dPsidGradU2[0][0][1][0]=0e0;
dPsidGradU2[0][0][1][1]=0.3e-1;
dPsidGradU2[0][0][1][2]=0e0;
dPsidGradU2[0][0][2][0]=0e0;
dPsidGradU2[0][0][2][1]=0e0;
dPsidGradU2[0][0][2][2]=0.3e-1;
dPsidGradU2[0][1][0][0]=-0.15000000000000013e-2;
dPsidGradU2[0][1][0][1]=0e0;
dPsidGradU2[0][1][0][2]=0e0;
dPsidGradU2[0][1][1][0]=0e0;
dPsidGradU2[0][1][1][1]=-0.15000000000000013e-2;
dPsidGradU2[0][1][1][2]=0e0;
dPsidGradU2[0][1][2][0]=0e0;
dPsidGradU2[0][1][2][1]=0e0;
dPsidGradU2[0][1][2][2]=-0.15000000000000013e-2;
dPsidGradU2[1][0][0][0]=-0.15000000000000013e-2;
dPsidGradU2[1][0][0][1]=0e0;
dPsidGradU2[1][0][0][2]=0e0;
dPsidGradU2[1][0][1][0]=0e0;
dPsidGradU2[1][0][1][1]=-0.15000000000000013e-2;
dPsidGradU2[1][0][1][2]=0e0;
dPsidGradU2[1][0][2][0]=0e0;
dPsidGradU2[1][0][2][1]=0e0;
dPsidGradU2[1][0][2][2]=-0.15000000000000013e-2;
dPsidGradU2[1][1][0][0]=0.3e-1;
dPsidGradU2[1][1][0][1]=0e0;
dPsidGradU2[1][1][0][2]=0e0;
dPsidGradU2[1][1][1][0]=0e0;
dPsidGradU2[1][1][1][1]=0.3e-1;
dPsidGradU2[1][1][1][2]=0e0;
dPsidGradU2[1][1][2][0]=0e0;
dPsidGradU2[1][1][2][1]=0e0;
dPsidGradU2[1][1][2][2]=0.3e-1;
};
