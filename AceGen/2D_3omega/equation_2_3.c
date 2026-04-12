/*************************************************************
* AceGen    9.105 Windows (19 Oct 25)                        *
*           Co. J. Korelc  2020           10 Apr 26 17:08:11 *
**************************************************************
User     : Limited evaluation version
Notebook : DoubleDitch3
Evaluation time                 : 2 s     Mode  : Optimal
Number of formulae              : 59      Method: Automatic
Subroutine                      : RandomEquation size: 1282
Total size of Mathematica  code : 1282 subexpressions
Total size of C code            : 4353 bytes */

/******************* S U B R O U T I N E *********************/

void equation(double v[242],double U[3],double U0[3],double 
     (*dt),double GradU[3][2],double dPsiDu[3],double dPsidGradU[3][2]
     ,double dPsiDu2[3][3],double dPsidUdGradU[3][3][2]
     ,double dPsidGradU2[3][3][2][2])
{
int i01;
v[7]=GradU[0][0];
v[8]=GradU[0][1];
v[9]=GradU[1][0];
v[10]=GradU[1][1];
v[11]=GradU[2][0];
v[12]=GradU[2][1];
v[13]=U[0];
v[53]=1e0-v[13];
v[56]=(v[53]*v[53]);
v[50]=2e0*v[13];
v[14]=U[1];
v[58]=1e0-v[14];
v[61]=(v[58]*v[58]);
v[51]=2e0*v[14];
v[125]=100e0*v[51];
v[15]=U[2];
v[63]=1e0-v[15];
v[64]=(v[63]*v[63]);
v[52]=2e0*v[15];
v[16]=U0[0];
v[17]=U0[1];
v[18]=U0[2];
v[68]=0.125e1/(*dt);
v[37]=-v[13]-v[14]+v[63];
v[54]=1e0-v[37];
v[55]=(v[54]*v[54]);
v[49]=-2e0*v[37];
v[97]=-0.6000000000000005e-2*v[49]*v[54];
v[66]=3e0+v[49];
v[95]=0.8e-1*v[49]+0.4e-1*v[66];
v[45]=(v[37]*v[37]);
v[96]=-0.6000000000000005e-2*v[45];
v[126]=-0.6000000000000005e-2*v[55]+v[96]+2e0*v[97];
v[124]=v[126]-0.15000000000000013e-2*v[68]+v[95];
v[113]=v[124]+v[125]*v[50];
v[107]=v[124]+v[125]*v[52];
v[106]=v[124]+100e0*v[50]*v[52];
v[127]=0.4e-1*v[45]+v[49]*(-0.30000000000000027e-2*v[55]+0.2e-1*v[66])+v[54]*v[96];
v[73]=0.15000000000000013e-2*(-1e0+v[16]+v[17]+v[18]+v[37]);
v[82]=-0.15000000000000013e-2*(v[11]+v[7]+v[9]);
v[84]=-0.15000000000000013e-2*(v[10]+v[12]+v[8]);
v[41]=(v[13]*v[13]);
v[130]=100e0*v[41];
v[115]=0.126e0*v[41];
v[42]=(v[14]*v[14]);
v[128]=100e0*v[42];
v[116]=200e0*v[42];
v[131]=v[116]+0.3e-1*v[68]+v[95];
v[132]=v[126]+v[131]+200e0*v[41];
v[111]=0.126e0*v[42];
v[43]=(v[15]*v[15]);
v[129]=100e0*v[43];
v[117]=200e0*v[43];
v[104]=0.126e0*v[43];
dPsiDu[0]=v[127]-v[115]*v[53]+v[50]*(v[128]+v[129]+0.63e-1*v[56])+v[68]*(0.315e-1*(v[13]-v[16])
 +v[73]);
dPsiDu[1]=v[127]-v[111]*v[58]+v[51]*(v[129]+v[130]+0.63e-1*v[61])+v[68]*(0.315e-1*(v[14]-v[17])
 +v[73]);
dPsiDu[2]=v[127]-v[104]*v[63]+v[52]*(v[128]+v[130]+0.63e-1*v[64])+v[68]*(0.315e-1*(v[15]-v[18])
 +v[73]);
dPsidGradU[0][0]=0.315e-1*v[7]+v[82];
dPsidGradU[0][1]=0.315e-1*v[8]+v[84];
dPsidGradU[1][0]=v[82]+0.315e-1*v[9];
dPsidGradU[1][1]=0.315e-1*v[10]+v[84];
dPsidGradU[2][0]=0.315e-1*v[11]+v[82];
dPsidGradU[2][1]=0.315e-1*v[12]+v[84];
dPsiDu2[0][0]=v[115]+v[117]+v[126]+v[131]-0.252e0*v[50]*v[53]+0.126e0*v[56];
dPsiDu2[0][1]=v[113];
dPsiDu2[0][2]=v[106];
dPsiDu2[1][0]=v[113];
dPsiDu2[1][1]=v[111]-v[116]+v[117]+v[132]-0.252e0*v[51]*v[58]+0.126e0*v[61];
dPsiDu2[1][2]=v[107];
dPsiDu2[2][0]=v[106];
dPsiDu2[2][1]=v[107];
dPsiDu2[2][2]=v[104]+v[132]-0.252e0*v[52]*v[63]+0.126e0*v[64];
dPsidUdGradU[0][0][0]=0e0;
dPsidUdGradU[0][0][1]=0e0;
dPsidUdGradU[0][1][0]=0e0;
dPsidUdGradU[0][1][1]=0e0;
dPsidUdGradU[0][2][0]=0e0;
dPsidUdGradU[0][2][1]=0e0;
dPsidUdGradU[1][0][0]=0e0;
dPsidUdGradU[1][0][1]=0e0;
dPsidUdGradU[1][1][0]=0e0;
dPsidUdGradU[1][1][1]=0e0;
dPsidUdGradU[1][2][0]=0e0;
dPsidUdGradU[1][2][1]=0e0;
dPsidUdGradU[2][0][0]=0e0;
dPsidUdGradU[2][0][1]=0e0;
dPsidUdGradU[2][1][0]=0e0;
dPsidUdGradU[2][1][1]=0e0;
dPsidUdGradU[2][2][0]=0e0;
dPsidUdGradU[2][2][1]=0e0;
dPsidGradU2[0][0][0][0]=0.3e-1;
dPsidGradU2[0][0][0][1]=0e0;
dPsidGradU2[0][0][1][0]=0e0;
dPsidGradU2[0][0][1][1]=0.3e-1;
dPsidGradU2[0][1][0][0]=-0.15000000000000013e-2;
dPsidGradU2[0][1][0][1]=0e0;
dPsidGradU2[0][1][1][0]=0e0;
dPsidGradU2[0][1][1][1]=-0.15000000000000013e-2;
dPsidGradU2[0][2][0][0]=-0.15000000000000013e-2;
dPsidGradU2[0][2][0][1]=0e0;
dPsidGradU2[0][2][1][0]=0e0;
dPsidGradU2[0][2][1][1]=-0.15000000000000013e-2;
dPsidGradU2[1][0][0][0]=-0.15000000000000013e-2;
dPsidGradU2[1][0][0][1]=0e0;
dPsidGradU2[1][0][1][0]=0e0;
dPsidGradU2[1][0][1][1]=-0.15000000000000013e-2;
dPsidGradU2[1][1][0][0]=0.3e-1;
dPsidGradU2[1][1][0][1]=0e0;
dPsidGradU2[1][1][1][0]=0e0;
dPsidGradU2[1][1][1][1]=0.3e-1;
dPsidGradU2[1][2][0][0]=-0.15000000000000013e-2;
dPsidGradU2[1][2][0][1]=0e0;
dPsidGradU2[1][2][1][0]=0e0;
dPsidGradU2[1][2][1][1]=-0.15000000000000013e-2;
dPsidGradU2[2][0][0][0]=-0.15000000000000013e-2;
dPsidGradU2[2][0][0][1]=0e0;
dPsidGradU2[2][0][1][0]=0e0;
dPsidGradU2[2][0][1][1]=-0.15000000000000013e-2;
dPsidGradU2[2][1][0][0]=-0.15000000000000013e-2;
dPsidGradU2[2][1][0][1]=0e0;
dPsidGradU2[2][1][1][0]=0e0;
dPsidGradU2[2][1][1][1]=-0.15000000000000013e-2;
dPsidGradU2[2][2][0][0]=0.3e-1;
dPsidGradU2[2][2][0][1]=0e0;
dPsidGradU2[2][2][1][0]=0e0;
dPsidGradU2[2][2][1][1]=0.3e-1;
};
