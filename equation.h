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
    std::vector<Tensor<1,3>> &GradU,
    Vector<double> &dPsiDu, 
    std::vector<Tensor<1,3>> &dPsidGradU,
    FullMatrix<double> &dPsiDu2,
    std::vector<std::vector<Tensor<1,3>>> &dPsidUdGradU, 
    std::vector<std::vector<Tensor<2,3>>> &dPsidGradU2)
{
v[3]=GradU[0][0];
v[4]=GradU[0][1];
v[5]=GradU[0][2];
v[6]=GradU[1][0];
v[87]=(v[6]*v[6]);
v[7]=GradU[1][1];
v[8]=GradU[1][2];
v[9]=U[0];
v[53]=(v[3]*v[3]);
v[36]=2e0*v[3];
v[60]=(v[4]*v[4]);
v[37]=2e0*v[4];
v[68]=(v[5]*v[5]);
v[38]=2e0*v[5];
v[61]=2e0*(v[53]+v[60]+v[68]);
v[39]=2e0*v[6];
v[40]=2e0*v[7];
v[41]=2e0*v[8];
v[26]=(v[7]*v[7])+v[87]+(v[8]*v[8]);
v[43]=std::pow(v[26],4);
v[84]=5e0*v[43];
v[85]=4e0*v[41]*v[84];
v[88]=4e0*v[40]*v[84];
v[76]=4e0*std::pow(v[26],5);
v[21]=v[3]*v[6]+v[4]*v[7]+v[5]*v[8];
v[86]=3e0*(v[21]*v[21]);
v[52]=v[5]*v[86];
v[78]=v[37]*v[52]+v[7]*v[85];
v[75]=v[36]*v[52]+v[6]*v[85];
v[67]=v[40]*v[52];
v[59]=v[39]*v[52];
v[51]=v[4]*v[86];
v[74]=v[36]*v[51]+v[6]*v[88];
v[71]=v[41]*v[51];
v[58]=v[39]*v[51];
v[50]=v[3]*v[86];
v[70]=v[41]*v[50];
v[64]=v[40]*v[50];
v[49]=v[8]*v[86];
v[63]=v[37]*v[38]+v[40]*v[49];
v[56]=v[36]*v[38]+v[39]*v[49];
v[48]=v[7]*v[86];
v[55]=v[36]*v[37]+v[39]*v[48];
v[65]=2e0*std::pow(v[21],3);
v[72]=v[41]*v[52]+v[65];
v[66]=v[40]*v[51]+v[65];
v[57]=v[39]*v[50]+v[65];
dPsiDu[0]=0.15e1*(v[9]*v[9]);
dPsiDu[1]=U[1]/4e0;
dPsidGradU[0][0]=v[3]*v[61]+v[6]*v[65];
dPsidGradU[0][1]=v[4]*v[61]+v[65]*v[7];
dPsidGradU[0][2]=v[5]*v[61]+v[65]*v[8];
dPsidGradU[1][0]=v[3]*v[65]+v[6]*v[76];
dPsidGradU[1][1]=v[4]*v[65]+v[7]*v[76];
dPsidGradU[1][2]=v[5]*v[65]+v[76]*v[8];
dPsiDu2[0][0]=3e0*v[9];
dPsiDu2[0][1]=0e0;
dPsiDu2[1][0]=0e0;
dPsiDu2[1][1]=0.25e0;
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
dPsidGradU2[0][0][0][0]=4e0*v[53]+v[61]+v[39]*v[6]*v[86];
dPsidGradU2[0][0][0][1]=v[55];
dPsidGradU2[0][0][0][2]=v[56];
dPsidGradU2[0][0][1][0]=v[55];
dPsidGradU2[0][0][1][1]=v[40]*v[48]+4e0*v[60]+v[61];
dPsidGradU2[0][0][1][2]=v[63];
dPsidGradU2[0][0][2][0]=v[56];
dPsidGradU2[0][0][2][1]=v[63];
dPsidGradU2[0][0][2][2]=v[41]*v[49]+v[61]+4e0*v[68];
dPsidGradU2[0][1][0][0]=v[57];
dPsidGradU2[0][1][0][1]=v[58];
dPsidGradU2[0][1][0][2]=v[59];
dPsidGradU2[0][1][1][0]=v[64];
dPsidGradU2[0][1][1][1]=v[66];
dPsidGradU2[0][1][1][2]=v[67];
dPsidGradU2[0][1][2][0]=v[70];
dPsidGradU2[0][1][2][1]=v[71];
dPsidGradU2[0][1][2][2]=v[72];
dPsidGradU2[1][0][0][0]=v[57];
dPsidGradU2[1][0][0][1]=v[64];
dPsidGradU2[1][0][0][2]=v[70];
dPsidGradU2[1][0][1][0]=v[58];
dPsidGradU2[1][0][1][1]=v[66];
dPsidGradU2[1][0][1][2]=v[71];
dPsidGradU2[1][0][2][0]=v[59];
dPsidGradU2[1][0][2][1]=v[67];
dPsidGradU2[1][0][2][2]=v[72];
dPsidGradU2[1][1][0][0]=v[36]*v[50]+v[76]+40e0*v[43]*v[87];
dPsidGradU2[1][1][0][1]=v[74];
dPsidGradU2[1][1][0][2]=v[75];
dPsidGradU2[1][1][1][0]=v[74];
dPsidGradU2[1][1][1][1]=v[37]*v[51]+v[76]+v[7]*v[88];
dPsidGradU2[1][1][1][2]=v[78];
dPsidGradU2[1][1][2][0]=v[75];
dPsidGradU2[1][1][2][1]=v[78];
dPsidGradU2[1][1][2][2]=v[38]*v[52]+v[76]+v[8]*v[85];
};
