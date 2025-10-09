/*************************************************************
* AceGen    9.102 Windows (8 Aug 25)                         *
*           Co. J. Korelc  2020           2 Sep 25 18:51:36  *
**************************************************************
User     : Limited evaluation version
Notebook : simple_code_for_Marc
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 6       Method: Automatic
Subroutine                      : laplace size: 161
Total size of Mathematica  code : 161 subexpressions
Total size of C code            : 721 bytes */

/******************* S U B R O U T I N E *********************/
unsigned int dim=3;
void laplace(std::vector<double> &v,double (*Uin),double (*Uold)
     ,dealii::Tensor<1, 3> &GradU,double (*dPsiDu),dealii::Tensor<1,3> &dPsidGradU,double (*dPsiDu2)
     ,dealii::Tensor<2,3> &dPsidGradU2,dealii::Tensor<1,3>& dPsidUdGradU,double (*deltat))
{
int i01;
v[6]=(*deltat);
(*dPsiDu)=((*Uin)-(*Uold))/v[6];
dPsidGradU[0]=0.1e0*GradU[0];
dPsidGradU[1]=0.1e0*GradU[1];
dPsidGradU[2]=0.1e0*GradU[2];
(*dPsiDu2)=1e0/v[6];
dPsidUdGradU[0]=0e0;
dPsidUdGradU[1]=0e0;
dPsidUdGradU[2]=0e0;
dPsidGradU2[0][0]=0.1e0;
dPsidGradU2[0][1]=0e0;
dPsidGradU2[0][2]=0e0;
dPsidGradU2[1][0]=0e0;
dPsidGradU2[1][1]=0.1e0;
dPsidGradU2[1][2]=0e0;
dPsidGradU2[2][0]=0e0;
dPsidGradU2[2][1]=0e0;
dPsidGradU2[2][2]=0.1e0;
};
