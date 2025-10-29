/*************************************************************
* AceGen    9.102 Windows (8 Aug 25)                         *
*           Co. J. Korelc  2020           23 Oct 25 12:34:24 *
**************************************************************
User     : Limited evaluation version
Notebook : simple_code_for_Marc
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 7       Method: Automatic
Subroutine                      : laplace-nonlinear size: 177
Total size of Mathematica  code : 177 subexpressions
Total size of C code            : 790 bytes */

/******************* S U B R O U T I N E *********************/

void laplace_nonlinear(std::vector<double> &v,double (*Uin),double (*Uold)
     ,dealii::Tensor<1, 3> &GradU, double (*dPsiDu), dealii::Tensor<1,3> &dPsidGradU, double (*dPsiDu2)
     ,dealii::Tensor<2,3> &dPsidGradU2,dealii::Tensor<1,3>& dPsidUdGradU,double (*deltat))
{
v[1]=(*Uin);
v[6]=(*deltat);
(*dPsiDu)=-1e0+(v[1]*v[1])+(-(*Uold)+v[1])/v[6];
dPsidGradU[0]=0.1e0*GradU[0];
dPsidGradU[1]=0.1e0*GradU[1];
dPsidGradU[2]=0.1e0*GradU[2];
(*dPsiDu2)=2e0*v[1]+1e0/v[6];
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