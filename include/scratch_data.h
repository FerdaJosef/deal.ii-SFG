#include "model.h"

template <int dim, int n>
Step3<dim, n>::AssemblyScratchData::AssemblyScratchData(const FiniteElement<dim> &fe)
  :
  fe_values(fe,
            QGauss<dim>(fe.degree + 1),
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values)
{
  const unsigned int n_q_points =
      fe_values.get_quadrature().size();

  values_newton.resize(n_q_points, Vector<double>(n));
  values_old.resize(n_q_points, Vector<double>(n));

  gradients_newton.resize(n_q_points,
                          std::vector<Tensor<1,dim>>(n));

  dPsiDu.reinit(n);
  dPsiDu2.reinit(n, n);

  dPsidGradU.resize(n);

  dPsidUdGradU.resize(n, std::vector<Tensor<1,dim>>(n));

  dPsidGradU2.resize(n, std::vector<Tensor<2,dim>>(n));

  acegen_scratch.resize(256);
}

template <int dim, int n>
Step3<dim, n>::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : fe_values(scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
    , values_newton(scratch_data.values_newton)
    , values_old(scratch_data.values_old)
    , gradients_newton(scratch_data.gradients_newton)
    , acegen_scratch(scratch_data.acegen_scratch)
    , dPsiDu(scratch_data.dPsiDu)
    , dPsidGradU(scratch_data.dPsidGradU)
    , dPsiDu2(scratch_data.dPsiDu2)
    , dPsidUdGradU(scratch_data.dPsidUdGradU)
    , dPsidGradU2(scratch_data.dPsidGradU2)
  {}
