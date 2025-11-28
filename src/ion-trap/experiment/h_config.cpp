#include <cmath>
#include <vector>
#include <numeric>

namespace experiment {

    struct Ion {
        int idx;
        bool interacts;
        std::vector<double> crosstalk_scale_factors;
    };

    class HConfig {
    public:
        HConfig() {};

        int mode;
        int n_ions;
        int n_gate_ions;
        int n_interacting_ions;
        int expansion_order;
        int spectator;

        double stark_scale_factor;
        double crosstalk_scale_factor;

        bool kerr_status;
        bool carrier_status;
        bool crosstalk_status;
        bool stark_status;
        bool constant_light_shift;
        bool monitoring_status;

        std::vector<int> gate_ions_idx;
        std::vector<Ion> interacting_ions;        

        std::vector<int> get_gate_ions() { return gate_ions; }

        void set_gate_ions(std::vector<int> _gate_ions) {
            gate_ions = _gate_ions; 
            n_gate_ions = gate_ions.size();

            // Update the list of neighbours
            std::vector<int>().swap(gate_ions_idx);
            std::vector<Ion>().swap(interacting_ions);
            for (int i = 0; i < n_ions; i++) { // n_ions

                Ion ion;
                ion.idx = i;
                ion.interacts = false;
                ion.crosstalk_scale_factors.resize(n_gate_ions);

                for (int j = 0; j < n_gate_ions; j++) { // gate_ions

                    int distance = std::fabs(i - gate_ions[j]);
                    if (distance == 0) {
                        gate_ions_idx.push_back(interacting_ions.size());

                        ion.interacts = true;
                        ion.crosstalk_scale_factors[j] = 1.0;
                    } else if (distance > 0 && distance < 2) { // 1st neighbourest neighbour
                        if (crosstalk_status) {
                            ion.interacts = true;
                            ion.crosstalk_scale_factors[j] = crosstalk_scale_factor;
                        }
                    } else {
                        ion.crosstalk_scale_factors[j] = 0.0;
                    }
                }

                if (ion.interacts) 
                    interacting_ions.push_back(ion);
            }
            n_interacting_ions = interacting_ions.size();
            if (monitoring_status) {n_interacting_ions += 1;}
        }

    private:
        std::vector<int> gate_ions;        
    };
}