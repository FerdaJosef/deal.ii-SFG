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
v[9]=U[0];
v[10]=U[1];
dPsiDu[0]=v[10]+std::pow(v[9],3);
dPsiDu[1]=std::pow(v[10],3)+v[9];
dPsidGradU[0][0]=GradU[0][0];
dPsidGradU[0][1]=GradU[0][1];
dPsidGradU[0][2]=GradU[0][2];
dPsidGradU[1][0]=GradU[1][0];
dPsidGradU[1][1]=GradU[1][1];
dPsidGradU[1][2]=GradU[1][2];
dPsiDu2[0][0]=3e0*(v[9]*v[9]);
dPsiDu2[0][1]=1e0;
dPsiDu2[1][0]=1e0;
dPsiDu2[1][1]=3e0*(v[10]*v[10]);
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
dPsidGradU2[0][0][0][0]=1e0;
dPsidGradU2[0][0][0][1]=0e0;
dPsidGradU2[0][0][0][2]=0e0;
dPsidGradU2[0][0][1][0]=0e0;
dPsidGradU2[0][0][1][1]=1e0;
dPsidGradU2[0][0][1][2]=0e0;
dPsidGradU2[0][0][2][0]=0e0;
dPsidGradU2[0][0][2][1]=0e0;
dPsidGradU2[0][0][2][2]=1e0;
dPsidGradU2[0][1][0][0]=0e0;
dPsidGradU2[0][1][0][1]=0e0;
dPsidGradU2[0][1][0][2]=0e0;
dPsidGradU2[0][1][1][0]=0e0;
dPsidGradU2[0][1][1][1]=0e0;
dPsidGradU2[0][1][1][2]=0e0;
dPsidGradU2[0][1][2][0]=0e0;
dPsidGradU2[0][1][2][1]=0e0;
dPsidGradU2[0][1][2][2]=0e0;
dPsidGradU2[1][0][0][0]=0e0;
dPsidGradU2[1][0][0][1]=0e0;
dPsidGradU2[1][0][0][2]=0e0;
dPsidGradU2[1][0][1][0]=0e0;
dPsidGradU2[1][0][1][1]=0e0;
dPsidGradU2[1][0][1][2]=0e0;
dPsidGradU2[1][0][2][0]=0e0;
dPsidGradU2[1][0][2][1]=0e0;
dPsidGradU2[1][0][2][2]=0e0;
dPsidGradU2[1][1][0][0]=1e0;
dPsidGradU2[1][1][0][1]=0e0;
dPsidGradU2[1][1][0][2]=0e0;
dPsidGradU2[1][1][1][0]=0e0;
dPsidGradU2[1][1][1][1]=1e0;
dPsidGradU2[1][1][1][2]=0e0;
dPsidGradU2[1][1][2][0]=0e0;
dPsidGradU2[1][1][2][1]=0e0;
dPsidGradU2[1][1][2][2]=1e0;
};
