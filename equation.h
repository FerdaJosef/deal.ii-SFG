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
v[81]=(v[6]*v[6]);
v[7]=GradU[1][1];
v[8]=GradU[1][2];
v[50]=(v[3]*v[3]);
v[33]=2e0*v[3];
v[57]=(v[4]*v[4]);
v[34]=2e0*v[4];
v[65]=(v[5]*v[5]);
v[35]=2e0*v[5];
v[58]=2e0*(v[50]+v[57]+v[65]);
v[36]=2e0*v[6];
v[37]=2e0*v[7];
v[38]=2e0*v[8];
v[24]=(v[7]*v[7])+v[81]+(v[8]*v[8]);
v[40]=std::pow(v[24],4);
v[78]=5e0*v[40];
v[79]=6e0*v[38]*v[78];
v[82]=6e0*v[37]*v[78];
v[73]=6e0*std::pow(v[24],5);
v[19]=v[3]*v[6]+v[4]*v[7]+v[5]*v[8];
v[80]=3e0*(v[19]*v[19]);
v[49]=v[5]*v[80];
v[75]=v[34]*v[49]+v[7]*v[79];
v[72]=v[33]*v[49]+v[6]*v[79];
v[64]=v[37]*v[49];
v[56]=v[36]*v[49];
v[48]=v[4]*v[80];
v[71]=v[33]*v[48]+v[6]*v[82];
v[68]=v[38]*v[48];
v[55]=v[36]*v[48];
v[47]=v[3]*v[80];
v[67]=v[38]*v[47];
v[61]=v[37]*v[47];
v[46]=v[8]*v[80];
v[60]=v[34]*v[35]+v[37]*v[46];
v[53]=v[33]*v[35]+v[36]*v[46];
v[45]=v[7]*v[80];
v[52]=v[33]*v[34]+v[36]*v[45];
v[62]=2e0*std::pow(v[19],3);
v[69]=v[38]*v[49]+v[62];
v[63]=v[37]*v[48]+v[62];
v[54]=v[36]*v[47]+v[62];
dPsiDu[0]=U[0];
dPsiDu[1]=U[1];
dPsidGradU[0][0]=v[3]*v[58]+v[6]*v[62];
dPsidGradU[0][1]=v[4]*v[58]+v[62]*v[7];
dPsidGradU[0][2]=v[5]*v[58]+v[62]*v[8];
dPsidGradU[1][0]=v[3]*v[62]+v[6]*v[73];
dPsidGradU[1][1]=v[4]*v[62]+v[7]*v[73];
dPsidGradU[1][2]=v[5]*v[62]+v[73]*v[8];
dPsiDu2[0][0]=1e0;
dPsiDu2[0][1]=0e0;
dPsiDu2[1][0]=0e0;
dPsiDu2[1][1]=1e0;
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
dPsidGradU2[0][0][0][0]=4e0*v[50]+v[58]+v[36]*v[6]*v[80];
dPsidGradU2[0][0][0][1]=v[52];
dPsidGradU2[0][0][0][2]=v[53];
dPsidGradU2[0][0][1][0]=v[52];
dPsidGradU2[0][0][1][1]=v[37]*v[45]+4e0*v[57]+v[58];
dPsidGradU2[0][0][1][2]=v[60];
dPsidGradU2[0][0][2][0]=v[53];
dPsidGradU2[0][0][2][1]=v[60];
dPsidGradU2[0][0][2][2]=v[38]*v[46]+v[58]+4e0*v[65];
dPsidGradU2[0][1][0][0]=v[54];
dPsidGradU2[0][1][0][1]=v[55];
dPsidGradU2[0][1][0][2]=v[56];
dPsidGradU2[0][1][1][0]=v[61];
dPsidGradU2[0][1][1][1]=v[63];
dPsidGradU2[0][1][1][2]=v[64];
dPsidGradU2[0][1][2][0]=v[67];
dPsidGradU2[0][1][2][1]=v[68];
dPsidGradU2[0][1][2][2]=v[69];
dPsidGradU2[1][0][0][0]=v[54];
dPsidGradU2[1][0][0][1]=v[61];
dPsidGradU2[1][0][0][2]=v[67];
dPsidGradU2[1][0][1][0]=v[55];
dPsidGradU2[1][0][1][1]=v[63];
dPsidGradU2[1][0][1][2]=v[68];
dPsidGradU2[1][0][2][0]=v[56];
dPsidGradU2[1][0][2][1]=v[64];
dPsidGradU2[1][0][2][2]=v[69];
dPsidGradU2[1][1][0][0]=v[33]*v[47]+v[73]+60e0*v[40]*v[81];
dPsidGradU2[1][1][0][1]=v[71];
dPsidGradU2[1][1][0][2]=v[72];
dPsidGradU2[1][1][1][0]=v[71];
dPsidGradU2[1][1][1][1]=v[34]*v[48]+v[73]+v[7]*v[82];
dPsidGradU2[1][1][1][2]=v[75];
dPsidGradU2[1][1][2][0]=v[72];
dPsidGradU2[1][1][2][1]=v[75];
dPsidGradU2[1][1][2][2]=v[35]*v[49]+v[73]+v[79]*v[8];
};
