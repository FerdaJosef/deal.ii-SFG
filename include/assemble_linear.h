#include "model.h"

template <int dim, int n>
void Step3<dim, n>::assemble_system()
{
  system_matrix = 0;
  system_rhs = 0;

  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &Step3<dim, n>::local_assemble_system,
                  &Step3<dim, n>::copy_local_to_global,
                  AssemblyScratchData(fe),
                  AssemblyCopyData());
}

template <int dim, int n>
void Step3<dim, n>::local_assemble_system(
        const typename DoFHandler<dim>::active_cell_iterator &cell,
        AssemblyScratchData                                  &scratch_data,
        AssemblyCopyData                                     &copy_data)
{
  const QGauss<dim> quadrature_formula(fe.degree + 1);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  //n_q_points    = quadrature_formula.size();

  copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  copy_data.cell_rhs.reinit(dofs_per_cell);

  copy_data.local_dof_indices.resize(dofs_per_cell);
 
  scratch_data.fe_values.reinit(cell);

  //======= ACEGEN=======

  scratch_data.fe_values.reinit(cell);

  copy_data.cell_matrix = 0;
  copy_data.cell_rhs    = 0;

  scratch_data.fe_values.get_function_gradients(solution, scratch_data.gradients_newton);
  scratch_data.fe_values.get_function_values(oldsolution, scratch_data.values_old);
  scratch_data.fe_values.get_function_values(solution, scratch_data.values_newton);

  const unsigned int cell_id = cell->active_cell_index();

  //right_hand_side(fe_values.get_quadrature_points(), rhs_values);

  for (const unsigned int q_index : scratch_data.fe_values.quadrature_point_indices())
    {
    //const auto &x_q = fe_values.quadrature_point(q_index);

    equation<dim, n>(scratch_data.acegen_scratch,
            scratch_data.values_newton[q_index],
            scratch_data.values_old[q_index],
        scratch_data.gradients_newton[q_index],
        scratch_data.dPsiDu,
        scratch_data.dPsidGradU,
        scratch_data.dPsiDu2,
        scratch_data.dPsidUdGradU,
        scratch_data.dPsidGradU2, // 2, 2, 3, 3
        &this->delta_t 
      );

      const auto &sd = scratch_data;
      for (const unsigned int i : scratch_data.fe_values.dof_indices())
        {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        
        // number of degrees of freedom is equql to number "geometric dof" = support points
        // function fe.system_to_component_index(i) returns a pair:
        // first is the corresponding component of vector system
        // second is the index of support point
        // shape value known about it through index i
        // ted tim vzdycky zredukuju dimenzi objektu a prevedu to na skalarni pripad

          copy_data.cell_rhs(i) -= (sd.fe_values.shape_value(i, q_index) * 
                          sd.dPsiDu[component_i]

                      + sd.fe_values.shape_grad(i, q_index) * sd.dPsidGradU[component_i]
                          + sd.fe_values.shape_value(i,q_index) * rhs_values[cell_id][q_index][component_i]) *
                          sd.fe_values.JxW(q_index);
          }
      for (const unsigned int i : scratch_data.fe_values.dof_indices())
      {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        
          for (const unsigned int j : scratch_data.fe_values.dof_indices())
          {
            const unsigned int component_j = fe.system_to_component_index(j).first;

            copy_data.cell_matrix(i, j) +=
                (sd.fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    sd.dPsidGradU2[component_i][component_j] *           // dPsi/d2(grad u)
                sd.fe_values.shape_grad(j, q_index)    // grad phi_j(x_q)
                + 
                sd.fe_values.shape_value(i, q_index) * // phi_i(x_q)
                    sd.dPsiDu2(component_i, component_j) *           // dPsi/d2( u)
                sd.fe_values.shape_value(j, q_index)     // phi_j(x_q)
                +
                sd.fe_values.shape_value(i, q_index) * //  phi_i(x_q)
                    sd.dPsidUdGradU[component_i][component_j] *           // dPsi/d(grad u)du
                sd.fe_values.shape_grad(j, q_index)     // grad phi_j(x_q)
                +
                sd.fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    sd.dPsidUdGradU[component_i][component_j] *           // dPsi/d(grad u)du
                sd.fe_values.shape_value(j, q_index)     // phi_j(x_q)
                ) *
                sd.fe_values.JxW(q_index);           // dx
        }
      }
    }
  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim, int n>
void Step3<dim, n>::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(
      copy_data.cell_matrix,
      copy_data.cell_rhs,
      copy_data.local_dof_indices,
      system_matrix,
      system_rhs);
  }