/*************************************************************
* AceGen    9.105 Windows (19 Oct 25)                        *
*           Co. J. Korelc  2020           4 Dec 25 10:01:16  *
**************************************************************
User     : Limited evaluation version
Notebook : AllenCahn
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 11      Method: Automatic
Subroutine                      : AllenCahn size: 222
Total size of Mathematica  code : 222 subexpressions
Total size of C code            : 902 bytes */

/******************* S U B R O U T I N E *********************/
void laplace_nonlinear(std::vector<double> &v,double (*Uin),double (*Uold)
     ,dealii::Tensor<1, 3> &GradU, double (*dPsiDu), dealii::Tensor<1,3> &dPsidGradU, double (*dPsiDu2)
     ,dealii::Tensor<2,3> &dPsidGradU2,dealii::Tensor<1,3>& dPsidUdGradU,double (*deltat))
{
int i01;
v[1]=(*Uin);
v[16]=1e0-v[1];
v[17]=(v[16]*v[16]);
v[13]=2e0*v[1];
v[6]=(*deltat);
v[26]=0.19999999999999998e2*v[13];
v[11]=(v[1]*v[1]);
v[31]=0.24e3*v[11];
(*dPsiDu)=-0.2e0*(-6e0*v[11]+3e0*v[13])+6e0*v[17]*v[26]-v[16]*v[31]+(-(*Uold)+v[1])/v[6];
dPsidGradU[0]=0.6000000000000001e0*GradU[0];
dPsidGradU[1]=0.6000000000000001e0*GradU[1];
dPsidGradU[2]=0.6000000000000001e0*GradU[2];
(*dPsiDu2)=-0.2e0*(6e0-6e0*v[13])+0.24e3*v[17]-24e0*v[16]*v[26]+v[31]+1e0/v[6];
dPsidUdGradU[0]=0e0;
dPsidUdGradU[1]=0e0;
dPsidUdGradU[2]=0e0;
dPsidGradU2[0][0]=0.6000000000000001e0;
dPsidGradU2[0][1]=0e0;
dPsidGradU2[0][2]=0e0;
dPsidGradU2[1][0]=0e0;
dPsidGradU2[1][1]=0.6000000000000001e0;
dPsidGradU2[1][2]=0e0;
dPsidGradU2[2][0]=0e0;
dPsidGradU2[2][1]=0e0;
dPsidGradU2[2][2]=0.6000000000000001e0;
}
