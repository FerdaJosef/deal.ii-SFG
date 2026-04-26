template <int dim, int n>
class RandomField {
public:
    RandomField() = default;

    void reinit(unsigned int n_cells, unsigned int n_q_points);

    void generate(const double delta_t, const double amplitude, const unsigned int n_active_cells);

    const Tensor<1,n>& get_value(unsigned int cell_id, unsigned int q_point) const {
        return data[cell_id][q_point];
    }

private:
    std::vector<std::vector<Tensor<1,n>>> data;
    std::mt19937 gen{std::random_device{}()};
    std::normal_distribution<double> dist{0.0, 1.0};
};