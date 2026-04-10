#include "AceGen/equation.h"

template <int dim, int n>
void wrapper(std::vector<double> &v,
              Vector<double> &U,
              Vector<double> &U0,
              std::vector<Tensor<1,dim>> &GradU,
              Vector<double> &dPsiDu,
              std::vector<Tensor<1,dim>> &dPsidGradU,
              FullMatrix<double> &dPsiDu2,
              std::vector<std::vector<Tensor<1,dim>>> &dPsidUdGradU,
              std::vector<std::vector<Tensor<2,dim>>> &dPsidGradU2,
              double *dt)
{
    double U_raw[n];
    double U0_raw[n];
    double GradU_raw[n][dim];

    double dPsiDu_raw[n];
    double dPsidGradU_raw[n][dim];
    double dPsiDu2_raw[n][n];
    double dPsidUdGradU_raw[n][n][dim];
    double dPsidGradU2_raw[n][n][dim][dim];

    for (unsigned int i = 0; i < n; ++i)
    {
        U_raw[i] = U[i];
        U0_raw[i] = U0[i];

        for (unsigned int d = 0; d < dim; ++d)
            GradU_raw[i][d] = GradU[i][d];

        for (unsigned int d = dim; d < dim; ++d)
            GradU_raw[i][d] = 0.0;
    }

    equation(v.data(),
                   U_raw,
                   U0_raw,
                   dt,
                   GradU_raw,
                   dPsiDu_raw,
                   dPsidGradU_raw,
                   dPsiDu2_raw,
                   dPsidUdGradU_raw,
                   dPsidGradU2_raw);

    for (unsigned int i = 0; i < n; ++i)
    {
        dPsiDu[i] = dPsiDu_raw[i];

        for (unsigned int d = 0; d < dim; ++d)
            dPsidGradU[i][d] = dPsidGradU_raw[i][d];

        for (unsigned int j = 0; j < n; ++j)
        {
            dPsiDu2(i,j) = dPsiDu2_raw[i][j];

            for (unsigned int d = 0; d < dim; ++d)
            {
                dPsidUdGradU[i][j][d] =
                    dPsidUdGradU_raw[i][j][d];

                for (unsigned int e = 0; e < dim; ++e)
                    dPsidGradU2[i][j][d][e] =
                        dPsidGradU2_raw[i][j][d][e];
            }
        }
    }
};