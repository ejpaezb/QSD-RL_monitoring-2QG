#include <vector>
#include <string>

namespace experiment {

    struct StateConfig {
        std::string spin_state_type;
        std::string motion_state_type;

        std::vector<double> n_bar;
        double cutoff_probability;

        std::vector<int> spin;
        std::vector<int> motion;

        int cutoff_mode = 0;
        std::vector<int> auto_cutoffs;
        std::vector<int> phonon_cutoffs;

        std::vector<int> dynamic_degrees;

        bool dynamic_cutoff = false;
        double cutoff_epsilon;
        int cutoff_pad_size;

        bool phonon_moving_basis = false;
        double shift_accuracy;
        bool use_x_basis = false;
    };
}