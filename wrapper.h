#define WRAPPER_H

#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <vector>

using namespace dealii;

template <int dim, int n>
inline void equation_wrapper(
    std::vector<double> &v,
    const Vector<double> &U,
    const Vector<double> &U0,
    const std::vector<Tensor<1,dim>> &GradU,
    Vector<double> &dPsiDu,
    std::vector<Tensor<1,dim>> &dPsidGradU,
    FullMatrix<double> &dPsiDu2,
    std::vector<std::vector<Tensor<1,dim>>> &dPsidUdGradU,
    std::vector<std::vector<Tensor<2,dim>>> &dPsidGradU2,
    double *dt)
{
    // ===== RAW ARRAYS (STACK, FAST) =====

    Assert(U.size() == n, ExcInternalError());
    Assert(dPsiDu.size() == n, ExcInternalError());

    double U_raw[n];
    double U0_raw[n];
    double GradU_raw[n][dim] = {};

    double dPsiDu_raw[n];
    double dPsidGradU_raw[n][dim];
    double dPsiDu2_raw[n][n];
    double dPsidUdGradU_raw[n][n][dim];
    double dPsidGradU2_raw[n][n][dim][dim];

    // ===== INPUT COPY (VERY CHEAP) =====
    for (unsigned int i = 0; i < n; ++i)
    {
        U_raw[i]  = U[i];
        U0_raw[i] = U0[i];

        for (unsigned int d = 0; d < dim; ++d)
            GradU_raw[i][d] = GradU[i][d];
    }

    // ===== CALL ACEGEN =====
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

    // ===== OUTPUT COPY =====
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
                dPsidUdGradU[i][j][d] = dPsidUdGradU_raw[i][j][d];

                for (unsigned int e = 0; e < dim; ++e)
                    dPsidGradU2[i][j][d][e] =
                        dPsidGradU2_raw[i][j][d][e];
            }
        }
    }
}