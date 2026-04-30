#include "RandomField.h"
#include <cmath>

using namespace dealii;

template <int dim, int n>
void RandomField<dim, n>::reinit(unsigned int n_cells, unsigned int n_q_points)
{
    data.resize(n_cells);
    for (auto &cell_data : data)
    {
        cell_data.resize(n_q_points);
    }
}

template <int dim, int n>
void RandomField<dim, n>::generate(const unsigned int cell_count,
                                   const unsigned int n_q_points,
                                   const double delta_t,
                                   const double amplitude)
{
    const double scale = amplitude / std::sqrt(delta_t);

    for (unsigned int c = 0; c < cell_count; ++c)
    {
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
            for (unsigned int v = 0; v < n; ++v)
            {
                data[c][q][v] = scale * dist(gen);
            }
        }
    }
}

// Explicitní instanciace
template class RandomField<1, 2>;