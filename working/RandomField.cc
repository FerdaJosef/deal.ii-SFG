#include "RandomField.h"
#include <deal.II/dofs/dof_handler.h>
#include <cmath>

template <int dim, int n>
void RandomField<dim, n>::reinit(unsigned int n_cells, unsigned int n_q_points)
{
    // Alokace vnějšího vektoru pro buňky
    data.resize(n_cells);
    for (auto &cell_data : data)
    {
        // Alokace vnitřního vektoru pro kvadraturní body
        cell_data.resize(n_q_points);
    }
}

template <int dim, int n>
void RandomField<dim, n>::generate(const DoFHandler<dim> &dof_handler,
                                   const unsigned int n_q_points,
                                   const double delta_t,
                                   const double amplitude)
{
    // Faktor pro škálování bílého šumu v čase
    const double scale = amplitude / std::sqrt(delta_t);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        const unsigned int cell_id = cell->active_cell_index();

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
            for (unsigned int v = 0; v < n; ++v)
            {
                // Přístup k datům a naplnění šumem
                data[cell_id][q][v] = scale * dist(gen);
            }
        }
    }
}

// Explicitní instanciace pro tvé kombinace
template class RandomField<1, 2>;