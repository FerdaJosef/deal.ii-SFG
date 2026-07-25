#pragma once

#include <vector>
#include <cmath>

#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

using namespace dealii;

/*************************************************************
* AceGen    9.201 Linux (9 Dec 25)                           *
*           Co. J. Korelc  2020           11 Jun 26 18:43:37 *
**************************************************************
User     : Limited evaluation version
Notebook : DoubleDitch4vars
Evaluation time                 : 1 s     Mode  : Optimal
Number of formulae              : 72      Method: Automatic
Subroutine                      : RandomEquation size: 1882
Total size of Mathematica  code : 1882 subexpressions
Total size of C code            : 6714 bytes */

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
int i01;int i02;int i03;
v[9]=GradU[0][0];
v[10]=GradU[0][1];
v[11]=GradU[1][0];
v[12]=GradU[1][1];
v[13]=GradU[2][0];
v[14]=GradU[2][1];
v[15]=GradU[3][0];
v[16]=GradU[3][1];
v[17]=U[0];
v[63]=1e0-v[17];
v[66]=(v[63]*v[63]);
v[59]=2e0*v[17];
v[18]=U[1];
v[68]=1e0-v[18];
v[71]=(v[68]*v[68]);
v[60]=2e0*v[18];
v[165]=100e0*v[60];
v[19]=U[2];
v[73]=1e0-v[19];
v[74]=(v[73]*v[73]);
v[61]=2e0*v[19];
v[166]=100e0*v[61];
v[20]=U[3];
v[76]=1e0-v[20];
v[77]=(v[76]*v[76]);
v[62]=2e0*v[20];
v[21]=U0[0];
v[22]=U0[1];
v[23]=U0[2];
v[24]=U0[3];
v[81]=1e0/(2e0*(*dt));
v[45]=-v[17]-v[18]-v[19]+v[76];
v[64]=1e0-v[45];
v[65]=(v[64]*v[64]);
v[58]=-2e0*v[45];
v[116]=-0.1199999999999999e-1*v[58]*v[64];
v[79]=3e0+v[58];
v[114]=0.8e-1*v[58]+0.4e-1*v[79];
v[54]=(v[45]*v[45]);
v[115]=-0.1199999999999999e-1*v[54];
v[163]=v[115]+2e0*v[116]-0.1199999999999999e-1*v[65];
v[164]=v[114]+v[163]-0.29999999999999977e-2*v[81];
v[155]=0.4e-1*v[54]+v[115]*v[64]+v[58]*(-0.5999999999999995e-2*v[65]+0.2e-1*v[79]);
v[86]=0.29999999999999977e-2*(-1e0+v[21]+v[22]+v[23]+v[24]+v[45]);
v[98]=-0.29999999999999977e-2*(v[11]+v[13]+v[15]+v[9]);
v[100]=-0.29999999999999977e-2*(v[10]+v[12]+v[14]+v[16]);
v[49]=(v[17]*v[17]);
v[161]=100e0*v[49];
v[146]=0.132e0*v[49];
v[139]=v[163]+200e0*v[49];
v[50]=(v[18]*v[18]);
v[158]=100e0*v[50];
v[147]=200e0*v[50];
v[140]=0.132e0*v[50];
v[132]=v[139]+v[147];
v[51]=(v[19]*v[19]);
v[160]=100e0*v[51];
v[156]=v[158]+v[160];
v[141]=200e0*v[51];
v[162]=v[114]+v[141]+0.3e-1*v[81];
v[133]=0.132e0*v[51];
v[52]=(v[20]*v[20]);
v[157]=100e0*v[52];
v[159]=v[157]+v[161];
v[167]=v[162]+200e0*v[52];
v[124]=0.132e0*v[52];
dPsiDu[0]=v[155]-v[146]*v[63]+v[59]*(v[156]+v[157]+0.6599999999999999e-1*v[66])+v[81]*
 (0.32999999999999996e-1*(v[17]-v[21])+v[86]);
dPsiDu[1]=v[155]-v[140]*v[68]+v[60]*(v[159]+v[160]+0.6599999999999999e-1*v[71])+v[81]*
 (0.32999999999999996e-1*(v[18]-v[22])+v[86]);
dPsiDu[2]=v[155]-v[133]*v[73]+v[61]*(v[158]+v[159]+0.6599999999999999e-1*v[74])+v[81]*
 (0.32999999999999996e-1*(v[19]-v[23])+v[86]);
dPsiDu[3]=v[155]-v[124]*v[76]+v[62]*(v[156]+v[161]+0.6599999999999999e-1*v[77])+v[81]*
 (0.32999999999999996e-1*(v[20]-v[24])+v[86]);
dPsidGradU[0][0]=0.32999999999999996e-1*v[9]+v[98];
dPsidGradU[0][1]=0.32999999999999996e-1*v[10]+v[100];
dPsidGradU[1][0]=0.32999999999999996e-1*v[11]+v[98];
dPsidGradU[1][1]=v[100]+0.32999999999999996e-1*v[12];
dPsidGradU[2][0]=0.32999999999999996e-1*v[13]+v[98];
dPsidGradU[2][1]=v[100]+0.32999999999999996e-1*v[14];
dPsidGradU[3][0]=0.32999999999999996e-1*v[15]+v[98];
dPsidGradU[3][1]=v[100]+0.32999999999999996e-1*v[16];
dPsiDu2[0][0]=v[146]+v[147]+v[163]+v[167]-0.264e0*v[59]*v[63]+0.132e0*v[66];
dPsiDu2[0][1]=v[164]+v[165]*v[59];
dPsiDu2[0][2]=v[164]+v[166]*v[59];
dPsiDu2[0][3]=v[164]+100e0*v[59]*v[62];
dPsiDu2[1][1]=v[139]+v[140]+v[167]-0.264e0*v[60]*v[68]+0.132e0*v[71];
dPsiDu2[1][2]=v[164]+1.e0*v[166]*v[60];
dPsiDu2[1][3]=v[164]+v[165]*v[62];
dPsiDu2[2][2]=v[132]+v[133]-v[141]+v[167]-0.264e0*v[61]*v[73]+0.132e0*v[74];
dPsiDu2[2][3]=v[164]+v[166]*v[62];
dPsiDu2[3][3]=v[124]+v[132]+v[162]-0.264e0*v[62]*v[76]+0.132e0*v[77];
for(i02=1;i02<4;i02++){
 for(i03=0;i03<i02;i03++){
    dPsiDu2[i02][i03]=dPsiDu2[i03][i02];}};
dPsidUdGradU[0][0][0]=0e0;
dPsidUdGradU[0][0][1]=0e0;
dPsidUdGradU[0][1][0]=0e0;
dPsidUdGradU[0][1][1]=0e0;
dPsidUdGradU[0][2][0]=0e0;
dPsidUdGradU[0][2][1]=0e0;
dPsidUdGradU[0][3][0]=0e0;
dPsidUdGradU[0][3][1]=0e0;
dPsidUdGradU[1][0][0]=0e0;
dPsidUdGradU[1][0][1]=0e0;
dPsidUdGradU[1][1][0]=0e0;
dPsidUdGradU[1][1][1]=0e0;
dPsidUdGradU[1][2][0]=0e0;
dPsidUdGradU[1][2][1]=0e0;
dPsidUdGradU[1][3][0]=0e0;
dPsidUdGradU[1][3][1]=0e0;
dPsidUdGradU[2][0][0]=0e0;
dPsidUdGradU[2][0][1]=0e0;
dPsidUdGradU[2][1][0]=0e0;
dPsidUdGradU[2][1][1]=0e0;
dPsidUdGradU[2][2][0]=0e0;
dPsidUdGradU[2][2][1]=0e0;
dPsidUdGradU[2][3][0]=0e0;
dPsidUdGradU[2][3][1]=0e0;
dPsidUdGradU[3][0][0]=0e0;
dPsidUdGradU[3][0][1]=0e0;
dPsidUdGradU[3][1][0]=0e0;
dPsidUdGradU[3][1][1]=0e0;
dPsidUdGradU[3][2][0]=0e0;
dPsidUdGradU[3][2][1]=0e0;
dPsidUdGradU[3][3][0]=0e0;
dPsidUdGradU[3][3][1]=0e0;
dPsidGradU2[0][0][0][0]=0.3e-1;
dPsidGradU2[0][0][0][1]=0e0;
dPsidGradU2[0][0][1][0]=0e0;
dPsidGradU2[0][0][1][1]=0.3e-1;
dPsidGradU2[0][1][0][0]=-0.29999999999999977e-2;
dPsidGradU2[0][1][0][1]=0e0;
dPsidGradU2[0][1][1][0]=0e0;
dPsidGradU2[0][1][1][1]=-0.29999999999999977e-2;
dPsidGradU2[0][2][0][0]=-0.29999999999999977e-2;
dPsidGradU2[0][2][0][1]=0e0;
dPsidGradU2[0][2][1][0]=0e0;
dPsidGradU2[0][2][1][1]=-0.29999999999999977e-2;
dPsidGradU2[0][3][0][0]=-0.29999999999999977e-2;
dPsidGradU2[0][3][0][1]=0e0;
dPsidGradU2[0][3][1][0]=0e0;
dPsidGradU2[0][3][1][1]=-0.29999999999999977e-2;
dPsidGradU2[1][0][0][0]=-0.29999999999999977e-2;
dPsidGradU2[1][0][0][1]=0e0;
dPsidGradU2[1][0][1][0]=0e0;
dPsidGradU2[1][0][1][1]=-0.29999999999999977e-2;
dPsidGradU2[1][1][0][0]=0.3e-1;
dPsidGradU2[1][1][0][1]=0e0;
dPsidGradU2[1][1][1][0]=0e0;
dPsidGradU2[1][1][1][1]=0.3e-1;
dPsidGradU2[1][2][0][0]=-0.29999999999999977e-2;
dPsidGradU2[1][2][0][1]=0e0;
dPsidGradU2[1][2][1][0]=0e0;
dPsidGradU2[1][2][1][1]=-0.29999999999999977e-2;
dPsidGradU2[1][3][0][0]=-0.29999999999999977e-2;
dPsidGradU2[1][3][0][1]=0e0;
dPsidGradU2[1][3][1][0]=0e0;
dPsidGradU2[1][3][1][1]=-0.29999999999999977e-2;
dPsidGradU2[2][0][0][0]=-0.29999999999999977e-2;
dPsidGradU2[2][0][0][1]=0e0;
dPsidGradU2[2][0][1][0]=0e0;
dPsidGradU2[2][0][1][1]=-0.29999999999999977e-2;
dPsidGradU2[2][1][0][0]=-0.29999999999999977e-2;
dPsidGradU2[2][1][0][1]=0e0;
dPsidGradU2[2][1][1][0]=0e0;
dPsidGradU2[2][1][1][1]=-0.29999999999999977e-2;
dPsidGradU2[2][2][0][0]=0.3e-1;
dPsidGradU2[2][2][0][1]=0e0;
dPsidGradU2[2][2][1][0]=0e0;
dPsidGradU2[2][2][1][1]=0.3e-1;
dPsidGradU2[2][3][0][0]=-0.29999999999999977e-2;
dPsidGradU2[2][3][0][1]=0e0;
dPsidGradU2[2][3][1][0]=0e0;
dPsidGradU2[2][3][1][1]=-0.29999999999999977e-2;
dPsidGradU2[3][0][0][0]=-0.29999999999999977e-2;
dPsidGradU2[3][0][0][1]=0e0;
dPsidGradU2[3][0][1][0]=0e0;
dPsidGradU2[3][0][1][1]=-0.29999999999999977e-2;
dPsidGradU2[3][1][0][0]=-0.29999999999999977e-2;
dPsidGradU2[3][1][0][1]=0e0;
dPsidGradU2[3][1][1][0]=0e0;
dPsidGradU2[3][1][1][1]=-0.29999999999999977e-2;
dPsidGradU2[3][2][0][0]=-0.29999999999999977e-2;
dPsidGradU2[3][2][0][1]=0e0;
dPsidGradU2[3][2][1][0]=0e0;
dPsidGradU2[3][2][1][1]=-0.29999999999999977e-2;
dPsidGradU2[3][3][0][0]=0.3e-1;
dPsidGradU2[3][3][0][1]=0e0;
dPsidGradU2[3][3][1][0]=0e0;
dPsidGradU2[3][3][1][1]=0.3e-1;
};
