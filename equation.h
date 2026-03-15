/*************************************************************
* AceGen    9.105 Windows (19 Oct 25)                        *
*           Co. J. Korelc  2020           29 Jan 26 10:35:07 *
**************************************************************
User     : Limited evaluation version
Notebook : RandomEquation
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 5       Method: Automatic
Subroutine                      : RandomEquation size: 550
Total size of Mathematica  code : 550 subexpressions
Total size of C code            : 1935 bytes */

/******************* S U B R O U T I N E *********************/

using namespace dealii;

void equation(std::vector<double> &v,
    Vector<double> &U,
    Vector<double> &U0,
    std::vector<Tensor<1,3>> &GradU,
    Vector<double> &dPsiDu, 
    std::vector<Tensor<1,3>> &dPsidGradU,
    FullMatrix<double> &dPsiDu2,
    std::vector<std::vector<Tensor<1,3>>> &dPsidUdGradU, 
    std::vector<std::vector<Tensor<2,3>>> &dPsidGradU2,
    double (*dt))
{
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
v[87]=-0.6000000000000005e-2*v[45]*v[49];
v[41]=(v[33]*v[33]);
v[86]=-0.6000000000000005e-2*v[41];
v[58]=3e0-2e0*v[41];
v[105]=-0.4e0*v[41]-0.6000000000000005e-2*v[50]+0.4e-1*v[58]+v[86]+2e0*v[87];
v[106]=v[105]+0.3e-1*v[60];
v[95]=v[105]+100e0*v[46]*v[47]-0.15000000000000013e-2*v[60];
v[104]=v[45]*(-0.4e-1*v[41]-0.30000000000000027e-2*v[50]+0.2e-1*v[58])+v[49]*v[86];
v[65]=0.15000000000000013e-2*(-1e0+v[13]+v[14]+v[33]);
v[72]=-0.15000000000000013e-2*(v[5]+v[8]);
v[74]=-0.15000000000000013e-2*(v[6]+v[9]);
v[76]=-0.15000000000000013e-2*(v[10]+v[7]);
v[38]=(v[11]*v[11]);
v[97]=0.126e0*v[38];
v[39]=(v[12]*v[12]);
v[93]=0.126e0*v[39];
dPsiDu[0]=v[104]+v[46]*(100e0*v[39]+0.63e-1*v[51])+v[60]*(0.315e-1*(v[11]-v[13])+v[65])
 -v[48]*v[97];
dPsiDu[1]=v[104]+v[47]*(100e0*v[38]+0.63e-1*v[56])+v[60]*(0.315e-1*(v[12]-v[14])+v[65])
 -v[53]*v[93];
dPsidGradU[0][0]=0.315e-1*v[5]+v[72];
dPsidGradU[0][1]=0.315e-1*v[6]+v[74];
dPsidGradU[0][2]=0.315e-1*v[7]+v[76];
dPsidGradU[1][0]=v[72]+0.315e-1*v[8];
dPsidGradU[1][1]=v[74]+0.315e-1*v[9];
dPsidGradU[1][2]=0.315e-1*v[10]+v[76];
dPsiDu2[0][0]=v[106]+200e0*v[39]-0.252e0*v[46]*v[48]+0.126e0*v[51]+v[97];
dPsiDu2[0][1]=v[95];
dPsiDu2[1][0]=v[95];
dPsiDu2[1][1]=v[106]+200e0*v[38]-0.252e0*v[47]*v[53]+0.126e0*v[56]+v[93];
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
