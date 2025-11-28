//
// Created by Shakib Vedaie on 2019-09-20.
//

#ifndef CVCQED_CV_CQED_H
#define CVCQED_CV_CQED_H

#include <qsd/ACG.h>
#include <qsd/Traject.h>
#include <qsd/State.h>
#include <qsd/Operator.h>
#include <qsd/AtomOp.h>
#include <qsd/FieldOp.h>
#include <qsd/SpinOp.h>
#include <qsd/Complex.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/foreach.hpp>

#include <cnpy/cnpy.h>
#include <qsd/qsd_experiment.h>

namespace experiment {

    /*
     * CV_CQED simulator
     */
class cv_cqed : public qsd::qsd_experiment {
    public:
        cv_cqed();
        cv_cqed(double _dt, int _numdts, int _numsteps, double _delta, std::vector<double> &_pulse);
        cv_cqed(const std::string &filename);

        enum Hamiltonian {
            TWO_ION_FULL, TWO_ION_SIMPLE, FIVE_ION, FIVE_ION_KERR, SMSPDC
        };
        enum Lindblads {
            FULL, NONE
        };
        enum Processor {
            FIDELITY, INFIDELITY, PARITY
        };
        struct States {
            std::vector<int> spin;
            std::vector<int> motion;
        };

        void set_flist(std::vector<string> _flist);

        std::vector<string> get_flist();
        std::vector<qsd::Operator> get_outlist();

        void set_processor(std::string type);

        void set_processor(Processor type);

        double processor(qsd::Result result);


        void parse_config(const std::string &config_path);

        qsd::Operator outlist_processor(const std::string &_operator);

        void set_stark_status(bool _stark_status);

        void set_delta(double _delta);

        void set_delta(int _delta_idx);

        double get_delta();

        void set_pulse(std::vector<double> &_pulse);

        void set_pulse(int _pulse_idx);

        std::vector<double> get_delta_list();

        void set_time(double _dt, int _numdts, int _numsteps);

        void set_time(double t);

        void set_nu_list(std::vector<double> _nu_list);

        void set_eta_list(std::vector<std::vector<double>> _eta_list);

        void set_ion_pair(std::vector<int> _ion_pair);

        void set_phonon_cutoffs(std::vector<int> _phonon_cutoffs, bool _dynamic_cutoff = false, double _epsilon = 1e-4,
                                int _pad_size = 2);

        bool get_dynamic_cutoff();

        std::vector<int> get_dynamic_degrees();

        double get_cutoff_epsilon();

        int get_cutoff_pad_size();

        void set_H(std::string hamiltonian);

        void set_H(Hamiltonian hamiltonian);

        qsd::Operator get_H();

        void set_L(std::string lindblads);

        void set_L(Lindblads lindblads);

        qsd::Operator *get_L();

        int get_nL();

        void set_state(States states);

        qsd::State get_state();


        void set_outlist(std::vector<qsd::Operator> _outlist);

        void set_phi_res(double _phi_res);

        double dt;      // basic time step
        int numdts;     // time interval between outputs = numdts * dt
        int numsteps;   // total integration time = numsteps * numdts * dt

        Complex I{0.0, 1.0};

        /*
         * Operators
         * */
        qsd::SigmaX sx_1{0};
        qsd::SigmaX sx_2{1};

        qsd::SigmaY sy_1{0};
        qsd::SigmaY sy_2{1};

        qsd::SigmaZ sz_1{0};
        qsd::SigmaZ sz_2{1};

        qsd::SigmaPlus sp_1{0};
        qsd::SigmaPlus sp_2{1};

        qsd::AnnihilationOperator a_1{2};
        qsd::AnnihilationOperator a_2{3};
        qsd::AnnihilationOperator a_3{4};
        qsd::AnnihilationOperator a_4{5};
        qsd::AnnihilationOperator a_5{6};

        qsd::NumberOperator N_1{2};
        qsd::NumberOperator N_2{3};
        qsd::NumberOperator N_3{4};
        qsd::NumberOperator N_4{5};
        qsd::NumberOperator N_5{6};

        qsd::IdentityOperator id_s1{0};
        qsd::IdentityOperator id_s2{1};
        qsd::IdentityOperator id_m1{2};
        qsd::IdentityOperator id_m2{3};
        qsd::IdentityOperator id_m3{4};
        qsd::IdentityOperator id_m4{5};
        qsd::IdentityOperator id_m5{6};

        qsd::Operator id_2 = id_s1 * id_s2 * id_m1 * id_m2;
        qsd::Operator id_5 = id_s1 * id_s2 * id_m1 * id_m2 * id_m3 * id_m4 * id_m5;

        qsd::Operator sm_1 = sp_1.hc();
        qsd::Operator sm_2 = sp_2.hc();

        qsd::Operator p00_1 = sm_1 * sp_1;
        qsd::Operator p00_2 = sm_2 * sp_2;

        qsd::Operator p11_1 = sp_1 * sm_1;
        qsd::Operator p11_2 = sp_2 * sm_2;

        qsd::Operator p01_1 = sm_1;
        qsd::Operator p01_2 = sm_2;

        qsd::Operator p10_1 = sp_1;
        qsd::Operator p10_2 = sp_2;

        // // Operator rho_ideal_1 = 0.5 * (p00_1 * p00_2 + p11_1 * p11_2 + I * (-p01_1 * p01_2 + p10_1 * p10_2));
        // // Operator rho_ideal_2 = 0.5 * (p00_1 * p00_2 + p11_1 * p11_2 + I * (p01_1 * p01_2 - p10_1 * p10_2));

        qsd::Operator rho_ideal_1 = 0.5 * (p00_1 * p00_2 + p11_1 * p11_2 + I * (-sm_1 * sm_2 + sp_1 * sp_2));
        qsd::Operator rho_ideal_2 =
                0.5 * (p00_1 * p00_2 + p11_1 * p11_2 + I * (sm_1 * sm_2 - sp_1 * sp_2)); // \psi_{ideal} of Manning

        qsd::Operator ad_1 = a_1.hc();
        qsd::Operator ad_2 = a_2.hc();
        qsd::Operator ad_3 = a_3.hc();
        qsd::Operator ad_4 = a_4.hc();
        qsd::Operator ad_5 = a_5.hc();

        qsd::Operator rho_cubic_1 = 0.16666 * (a_1 * a_1 * a_1 - ad_1 * ad_1 * ad_1);
        qsd::Operator rho_cubic_2 = 0.16666 * (a_1 * a_1 * a_1 + ad_1 * ad_1 * ad_1);

        qsd::Operator x_1 = ad_1 + a_1;
        qsd::Operator x_2 = ad_2 + a_2;
        qsd::Operator x_3 = ad_3 + a_3;
        qsd::Operator x_4 = ad_4 + a_4;
        qsd::Operator x_5 = ad_5 + a_5;
        qsd::Operator x2 = x_3 * x_3;

        qsd::Operator p_1 = ad_1 - a_1;
        qsd::Operator p_2 = ad_2 - a_2;
        qsd::Operator p_3 = ad_3 - a_3;
        qsd::Operator p_4 = ad_4 - a_4;
        qsd::Operator p_5 = ad_5 - a_5;
        qsd::Operator p2 = p_3 * p_3;

    private:
        std::string name;
        std::string description;

        qsd::Operator H;
        qsd::State psi0;
        qsd::Operator *L;
        int nL = 0;
        std::vector<string> flist;      // Output files
        std::vector<qsd::Operator> outlist;  // Operators to output

        Processor processor_type;
        std::vector<int> ion_pair;

        const double stark_scale_factor = 0.0004938763;  // The scale factor used to compute the Stark shift for each pulse segment
        std::vector<double> stark_shift;
        bool stark_status;

        std::vector<std::vector<double>> pulse_list;
        std::vector<std::vector<double>> eta_list;      // Lamb-Dicke parameter for the ions coupling to the motional modes
        std::vector<double> delta_list;
        std::vector<double> nu_list;                    // Motional mode frequencies

        std::vector<double> pulse;
        int pulse_steps;

        double delta;
        double t_gate;

        const double gamma_heating = 5e-4;          // Heating rate [1/us] (300 quanta/second)
        const double gamma_Rayleigh = 2 * 7.5e-10;  // Rayleigh scattering [1/us] (2 * 7.5e-4 1/s)
        const double gamma_Raman = 2 * 15e-6;       // Raman scattering [1/us] (2 * 15 1/s)
        const double amp_fluct = 2 * 1e-3;

        std::vector<int> motional_states;
        std::vector<int> phonon_cutoffs;
        std::vector<int> dynamic_degrees;
        bool dynamic_cutoff = false;
        double cutoff_epsilon;
        int cutoff_pad_size;

        double phi_res;

        // Full Hamiltonian time-dependent coefficients for 2 ions
        double H0_coeff(double t);

        Complex H1_coeff(double t);

        Complex H2_coeff(double t);

        Complex H3_coeff(double t);

        Complex H4_coeff(double t);

        double H5_coeff(double t);

        double H6_coeff(double t);

        Complex H7_coeff(double t);

        Complex H8_coeff(double t);

        Complex H9_coeff(double t);

        Complex H10_coeff(double t);

        Complex H11_coeff(double t);

        Complex H12_coeff(double t);

        Complex H13_coeff(double t);

        Complex H14_coeff(double t);

        // Simple Hamiltonian time-dependent coefficients for 2 ions
        Complex H1_coeff_simple(double t);

        Complex H2_coeff_simple(double t);

        Complex H3_coeff_simple(double t);

        Complex H4_coeff_simple(double t);

        // Simple Hamiltonian time-dependent coefficients for 5 ions
        double H0_coeff_5(double t);

        Complex H1_coeff_5(double t);

        Complex H2_coeff_5(double t);

        Complex H3_coeff_5(double t);

        Complex H4_coeff_5(double t);

        Complex H5_coeff_5(double t);

        Complex H6_coeff_5(double t);

        Complex H7_coeff_5(double t);

        Complex H8_coeff_5(double t);

        Complex H9_coeff_5(double t);

        Complex H10_coeff_5(double t);

        Complex omega(double t);

        Complex omega_d(double t);

        Complex g_t(double t);


        // Lindblad operators time-dependent coefficients
        double gamma_motion_dephasing_coeff(double t);

        double gamma_spin_dephasing_coeff(double t);

        double gamma_spin_flip_coeff(double t);

        void set_stark_shift();

        // Utility methods
        void load_1d_list(std::string filename, std::vector<double> &_list);

        void load_2d_list(std::string filename, std::vector<std::vector<double>> &_list);
    };
}
#endif //IONTRAP_ION_TRAP_H
