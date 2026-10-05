/*************************************************************
* AceGen    9.105 Windows (19 Oct 25)                        *
*           Co. J. Korelc  2020           9 Apr 26 22:07:57  *
**************************************************************
User     : Limited evaluation version
Notebook : DoubleDitch2D
Evaluation time                 : 1 s     Mode  : Optimal
Number of formulae              : 39      Method: Automatic
Subroutine                      : RandomEquation size: 730
Total size of Mathematica  code : 730 subexpressions
Total size of C code            : 2544 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void RandomEquation(double v[207],double U[2],double U0[2],double 
     (*dt),double GradU[2][2],double dPsiDu[2],double dPsidGradU[2][2]
     ,double dPsiDu2[2][2],double dPsidUdGradU[2][2][2]
     ,double dPsidGradU2[2][2][2][2])
{
int i01;
v[9]=U[0];
v[43]=1e0-v[9];
v[46]=(v[43]*v[43]);
v[41]=2e0*v[9];
v[10]=U[1];
v[48]=1e0-v[10];
v[51]=(v[48]*v[48]);
v[42]=2e0*v[10];
v[11]=U0[0];
v[12]=U0[1];
v[13]=GradU[0][0];
v[14]=GradU[0][1];
v[15]=GradU[1][0];
v[16]=GradU[1][1];
v[55]=0.125e1/(*dt);
v[29]=v[48]-v[9];
v[44]=1e0-v[29];
v[45]=(v[44]*v[44]);
v[40]=-2e0*v[29];
v[78]=-0.6000000000000005e-2*v[40]*v[44];
v[53]=3e0+v[40];
v[36]=(v[29]*v[29]);
v[77]=-0.6000000000000005e-2*v[36];
v[96]=0.8e-1*v[40]-0.6000000000000005e-2*v[45]+0.4e-1*v[53]+v[77]+2e0*v[78];
v[97]=0.3e-1*v[55]+v[96];
v[86]=100e0*v[41]*v[42]-0.15000000000000013e-2*v[55]+v[96];
v[95]=0.4e-1*v[36]+v[40]*(-0.30000000000000027e-2*v[45]+0.2e-1*v[53])+v[44]*v[77];
v[60]=0.15000000000000013e-2*(-1e0+v[11]+v[12]+v[29]);
v[66]=-0.15000000000000013e-2*(v[13]+v[15]);
v[68]=-0.15000000000000013e-2*(v[14]+v[16]);
v[33]=(v[9]*v[9]);
v[88]=0.126e0*v[33];
v[34]=(v[10]*v[10]);
v[84]=0.126e0*v[34];
dPsiDu[0]=v[41]*(100e0*v[34]+0.63e-1*v[46])-v[43]*v[88]+v[55]*(v[60]+0.315e-1*(-v[11]+v[9]))+v[95];
dPsiDu[1]=v[42]*(100e0*v[33]+0.63e-1*v[51])+v[55]*(0.315e-1*(v[10]-v[12])+v[60])-v[48]*v[84]+v[95];
dPsidGradU[0][0]=0.315e-1*v[13]+v[66];
dPsidGradU[0][1]=0.315e-1*v[14]+v[68];
dPsidGradU[1][0]=0.315e-1*v[15]+v[66];
dPsidGradU[1][1]=0.315e-1*v[16]+v[68];
dPsiDu2[0][0]=200e0*v[34]-0.252e0*v[41]*v[43]+0.126e0*v[46]+v[88]+v[97];
dPsiDu2[0][1]=v[86];
dPsiDu2[1][0]=v[86];
dPsiDu2[1][1]=200e0*v[33]-0.252e0*v[42]*v[48]+0.126e0*v[51]+v[84]+v[97];
dPsidUdGradU[0][0][0]=0e0;
dPsidUdGradU[0][0][1]=0e0;
dPsidUdGradU[0][1][0]=0e0;
dPsidUdGradU[0][1][1]=0e0;
dPsidUdGradU[1][0][0]=0e0;
dPsidUdGradU[1][0][1]=0e0;
dPsidUdGradU[1][1][0]=0e0;
dPsidUdGradU[1][1][1]=0e0;
dPsidGradU2[0][0][0][0]=0.3e-1;
dPsidGradU2[0][0][0][1]=0e0;
dPsidGradU2[0][0][1][0]=0e0;
dPsidGradU2[0][0][1][1]=0.3e-1;
dPsidGradU2[0][1][0][0]=-0.15000000000000013e-2;
dPsidGradU2[0][1][0][1]=0e0;
dPsidGradU2[0][1][1][0]=0e0;
dPsidGradU2[0][1][1][1]=-0.15000000000000013e-2;
dPsidGradU2[1][0][0][0]=-0.15000000000000013e-2;
dPsidGradU2[1][0][0][1]=0e0;
dPsidGradU2[1][0][1][0]=0e0;
dPsidGradU2[1][0][1][1]=-0.15000000000000013e-2;
dPsidGradU2[1][1][0][0]=0.3e-1;
dPsidGradU2[1][1][0][1]=0e0;
dPsidGradU2[1][1][1][0]=0e0;
dPsidGradU2[1][1][1][1]=0.3e-1;
};
