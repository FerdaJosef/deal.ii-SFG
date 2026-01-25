/*************************************************************
* AceGen    9.105 Windows (19 Oct 25)                        *
*           Co. J. Korelc  2020           14 Jan 26 22:01:02 *
**************************************************************
User     : Limited evaluation version
Notebook : RandomEquation
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 6       Method: Automatic
Subroutine                      : RandomEquation size: 569
Total size of Mathematica  code : 569 subexpressions
Total size of C code            : 1966 bytes */

/******************* S U B R O U T I N E *********************/
void RandomEquation(double v[125],double U[2],double GradU[2][3]
     ,double dPsiDu[2],double dPsidGradU[2][3],double dPsiDu2[2][2]
     ,double dPsidGradU2[2][3][2][3],double dPsidUdGradU[2][2][3])
{
v[10]=U[0]+U[1];
dPsiDu[0]=v[10];
dPsiDu[1]=v[10];
dPsidGradU[0][0]=-GradU[0][0];
dPsidGradU[0][1]=-GradU[0][1];
dPsidGradU[0][2]=-GradU[0][2];
dPsidGradU[1][0]=-GradU[1][0];
dPsidGradU[1][1]=-GradU[1][1];
dPsidGradU[1][2]=-GradU[1][2];
dPsiDu2[0][0]=1e0;
dPsiDu2[0][1]=1e0;
dPsiDu2[1][0]=1e0;
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
dPsidGradU2[0][0][0][0]=-1e0;
dPsidGradU2[0][0][0][1]=0e0;
dPsidGradU2[0][0][0][2]=0e0;
dPsidGradU2[0][0][1][0]=0e0;
dPsidGradU2[0][0][1][1]=-1e0;
dPsidGradU2[0][0][1][2]=0e0;
dPsidGradU2[0][0][2][0]=0e0;
dPsidGradU2[0][0][2][1]=0e0;
dPsidGradU2[0][0][2][2]=-1e0;
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
dPsidGradU2[1][1][0][0]=-1e0;
dPsidGradU2[1][1][0][1]=0e0;
dPsidGradU2[1][1][0][2]=0e0;
dPsidGradU2[1][1][1][0]=0e0;
dPsidGradU2[1][1][1][1]=-1e0;
dPsidGradU2[1][1][1][2]=0e0;
dPsidGradU2[1][1][2][0]=0e0;
dPsidGradU2[1][1][2][1]=0e0;
dPsidGradU2[1][1][2][2]=-1e0;
};
