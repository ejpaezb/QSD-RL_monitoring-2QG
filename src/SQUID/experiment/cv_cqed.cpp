//
// Created by Shakib Vedaie on 2019-09-20.
//

#include "cv_cqed.h"

namespace experiment {

    /*
     * continuous variable CQED simulator.
     * */

    cv_cqed::cv_cqed() {}

    cv_cqed::cv_cqed(double _dt, int _numdts, int _numsteps, double _delta, std::vector<double> &_pulse) {
        dt = _dt;
        numdts = _numdts;
        numsteps = _numsteps;

        t_gate = _dt * (double) (_numdts * _numsteps);

        delta = _delta;

        pulse.resize(_pulse.size());
        for (int i = 0; i < _pulse.size(); i++) {
            pulse[i] = _pulse[i];
        }
        pulse_steps = _pulse.size();

        set_stark_shift();
    }

    cv_cqed::cv_cqed(const std::string &filename) { parse_config(filename); }

    qsd::Operator cv_cqed::get_H() { return H; }

    qsd::Operator* cv_cqed::get_L() { return L; }

    int cv_cqed::get_nL() { return nL; }

    qsd::State cv_cqed::get_state() { return psi0; }

    std::vector<string> cv_cqed::get_flist() { return flist; }

    std::vector<qsd::Operator> cv_cqed::get_outlist() { return outlist; }

    void cv_cqed::parse_config(const std::string &config_path)
    {
        // Create empty property tree object
        boost::property_tree::ptree tree;

        // Parse the INFO into the property tree.
        boost::property_tree::read_info(config_path, tree);

        // Read the experiment configurations
        name = tree.get<std::string>("experiment.name");
        description = tree.get<std::string>("experiment.description");

        load_1d_list(tree.get<std::string>("experiment.delta_list_path"), delta_list);
        load_2d_list(tree.get<std::string>("experiment.pulse_list_path"), pulse_list);
        load_1d_list(tree.get<std::string>("experiment.nu_list_path"), nu_list);
        load_2d_list(tree.get<std::string>("experiment.eta_list_path"), eta_list);

        set_nu_list(nu_list);
        set_eta_list(eta_list);

        std::vector<int> _ion_pair;
        for (auto& item : tree.get_child("experiment.ion_pair")) {
            _ion_pair.emplace_back(item.second.get_value<int>());
        }
        set_ion_pair(_ion_pair);

        std::vector<int> _phonon_cutoffs;
        for (auto& item : tree.get_child("experiment.phonon_cutoffs")) {
            _phonon_cutoffs.emplace_back(item.second.get_value<int>());
        }

        bool _dynamic_cutoff = tree.get<bool>("experiment.dynamic_cutoff");
        double _cutoff_epsilon = tree.get<double>("experiment.cutoff_epsilon");
        int _cutoff_pad_size = tree.get<int>("experiment.cutoff_pad_size");

        set_phonon_cutoffs(_phonon_cutoffs, _dynamic_cutoff, _cutoff_epsilon, _cutoff_pad_size);

        std::string _hamiltonian = tree.get<std::string>("experiment.hamiltonian");
        set_H(_hamiltonian);

        std::string _lindblads = tree.get<std::string>("experiment.lindblads");
        set_L(_lindblads);

        States _states;
        for (auto& item : tree.get_child("experiment.state.spin.PURE")) {
            _states.spin.emplace_back(item.second.get_value<int>());
        }
        for (auto& item : tree.get_child("experiment.state.motion.PURE")) {
            _states.motion.emplace_back(item.second.get_value<int>());
        }
        set_state(_states);

        set_stark_status(tree.get<bool>("experiment.stark_status"));

        std::vector<qsd::Operator> _outlist;
        for (auto& item : tree.get_child("experiment.outlist")) {
            _outlist.emplace_back(outlist_processor(item.second.get_value<std::string>()));
        }
        set_outlist(_outlist);

        std::vector<std::string> _flist;
        for (auto& item : tree.get_child("experiment.flist")) {
            _flist.emplace_back(item.second.get_value<std::string>());
        }
        set_flist(_flist);

        double _dt = tree.get<double>("experiment.time.dt");
        int _numdts = tree.get<int>("experiment.time.numdts");
        int _numsteps = tree.get<int>("experiment.time.numsteps");

        set_time(_dt, _numdts, _numsteps);

        std::string _processor = (tree.get<std::string>("experiment.processor"));
        set_processor(_processor);

        if (_processor == "PARITY")
            set_phi_res(tree.get<double>("experiment.PARITY.phi_res"));
    }

    qsd::Operator cv_cqed::outlist_processor(const std::string &_operator)
    {
        /*
        * Operators
        * */

        if (_operator == "sx_1")
            return sx_1;
        else if (_operator == "sx_2")
            return sx_2;

        else if (_operator == "sy_1")
            return sy_1;
        else if (_operator == "sy_2")
            return sy_2;

        else if (_operator == "sz_1")
            return sz_1;
        else if (_operator == "sz_2")
            return sz_2;

        else if (_operator == "sp_1")
            return sp_1;
        else if (_operator == "sp_2")
            return sp_2;

        else if (_operator == "a_1")
            return a_1;
        else if (_operator == "a_2")
            return a_2;
        else if (_operator == "a_3")
            return a_3;
        else if (_operator == "a_4")
            return a_4;
        else if (_operator == "a_5")
            return a_5;

        else if (_operator == "N_1")
            return N_1;
        else if (_operator == "N_2")
            return N_2;
        else if (_operator == "N_3")
            return N_3;
        else if (_operator == "N_4")
            return N_4;
        else if (_operator == "N_5")
            return N_5;

        else if (_operator == "id_s1")
            return id_s1;
        else if (_operator == "id_s2")
            return id_s2;
        else if (_operator == "id_m1")
            return id_m1;
        else if (_operator == "id_m2")
            return id_m2;
        else if (_operator == "id_m3")
            return id_m3;
        else if (_operator == "id_m4")
            return id_m4;
        else if (_operator == "id_m5")
            return id_m5;

        else if (_operator == "id_2")
            return id_2;
        else if (_operator == "id_5")
            return id_5;

        else if (_operator == "sm_1")
            return sm_1;
        else if (_operator == "sm_2")
            return sm_2;

        else if (_operator == "p00_1")
            return p00_1;
        else if (_operator == "p00_2")
            return p00_2;

        else if (_operator == "p11_1")
            return p11_1;
        else if (_operator == "p11_2")
            return p11_2;

        else if (_operator == "rho_ideal_1")
            return rho_ideal_1;
        else if (_operator == "rho_ideal_2")
            return rho_ideal_2;
        else if (_operator == "rho_cubic_1")
            return rho_cubic_1;
        else if (_operator == "rho_cubic_2")
            return rho_cubic_2;

        else if (_operator == "ad_1")
            return ad_1;
        else if (_operator == "ad_2")
            return ad_2;
        else if (_operator == "ad_3")
            return ad_3;
        else if (_operator == "ad_4")
            return ad_4;
        else if (_operator == "ad_5")
            return ad_5;

        else if (_operator == "x_1")
            return x_1;
        else if (_operator == "x_2")
            return x_2;
        else if (_operator == "x_3")
            return x_3;
        else if (_operator == "x_4")
            return x_4;
        else if (_operator == "x_5")
            return x_5;

        else if (_operator == "p_1")
            return p_1;
        else if (_operator == "p_2")
            return p_2;
        else if (_operator == "p_3")
            return p_3;
        else if (_operator == "p_4")
            return p_4;
        else if (_operator == "p_5")
            return p_5;

        else if (_operator == "x2")
            return x2;
        else if (_operator == "p2")
            return p2;

        std::cout << "Error: Unrecognized operator in outlist.\n";
        return qsd::NullOperator();
    }

    void cv_cqed::set_stark_status(bool _stark_status)
    {
        stark_status = _stark_status;
    }

    void cv_cqed::set_delta(double _delta) {

        delta = _delta;
    }

    void cv_cqed::set_delta(int _delta_idx) {

        set_delta(delta_list[_delta_idx]);
    }

    double cv_cqed::get_delta() { return delta; }

    void cv_cqed::set_pulse(std::vector<double> &_pulse) {

        pulse.resize(_pulse.size());
        for (int i = 0; i < _pulse.size(); i++) {
            pulse[i] = _pulse[i];
        }
        pulse_steps = _pulse.size();

        set_stark_shift();
    }

    void cv_cqed::set_pulse(int _pulse_idx) {

        set_pulse(pulse_list[_pulse_idx]);
    }

    std::vector<double> cv_cqed::get_delta_list() { return delta_list; }

    void cv_cqed::set_time(double _dt, int _numdts, int _numsteps) {
        dt = _dt;
        numdts = _numdts;
        numsteps = _numsteps;

        t_gate = _dt * (double) (_numdts * _numsteps);
    }

    void cv_cqed::set_time(double t) {
        dt = 0.1;
        numdts = 1;
        numsteps = (int) (t / 0.1);

        t_gate = 0.1 * (double) (1 * (int) (t / 0.1));
    }

    void cv_cqed::set_nu_list(std::vector<double> _nu_list) { nu_list = _nu_list; } // The input is in Mrad/s
    void cv_cqed::set_eta_list(std::vector<std::vector<double>> _eta_list) { eta_list = _eta_list; }

    void cv_cqed::set_ion_pair(std::vector<int> _ion_pair) { ion_pair = _ion_pair; }

    void cv_cqed::set_phonon_cutoffs(std::vector<int> _phonon_cutoffs, bool _dynamic_cutoff, double _epsilon,
                                      int _pad_size) {
        phonon_cutoffs = _phonon_cutoffs;
        dynamic_cutoff = _dynamic_cutoff;
        cutoff_epsilon = _epsilon;
        cutoff_pad_size = _pad_size;
    }
    bool cv_cqed::get_dynamic_cutoff() { return dynamic_cutoff; };
    std::vector<int> cv_cqed::get_dynamic_degrees() { return dynamic_degrees; };
    double cv_cqed::get_cutoff_epsilon() { return cutoff_epsilon; };
    int cv_cqed::get_cutoff_pad_size() { return cutoff_pad_size; };

// Full Hamiltonian time-dependent coefficients
    double cv_cqed::H0_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H1_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[0] * t) - I * sin(nu_list[0] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H2_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[0] * t) + I * sin(nu_list[0] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H3_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[1] * t) - I * sin(nu_list[1] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H4_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[1] * t) + I * sin(nu_list[1] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    double cv_cqed::H5_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * 2 * cos((delta + stark_shift[idx]) * t);
    }

    double cv_cqed::H6_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H7_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(2 * nu_list[0] * t) - I * sin(2 * nu_list[0] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H8_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(2 * nu_list[0] * t) + I * sin(2 * nu_list[0] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H9_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(2 * nu_list[1] * t) - I * sin(2 * nu_list[1] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H10_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(2 * nu_list[1] * t) + I * sin(2 * nu_list[1] * t)) * 2 * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H11_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos((nu_list[0] + nu_list[1]) * t) - I * sin((nu_list[0] + nu_list[1]) * t)) * 2 *
               cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H12_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos((-nu_list[0] + nu_list[1]) * t) + I * sin((-nu_list[0] + nu_list[1]) * t)) * 2 *
               cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H13_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos((nu_list[0] - nu_list[1]) * t) + I * sin((nu_list[0] - nu_list[1]) * t)) * 2 *
               cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H14_coeff(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos((nu_list[0] + nu_list[1]) * t) + I * sin((nu_list[0] + nu_list[1]) * t)) * 2 *
               cos((delta + stark_shift[idx]) * t);
    }

    // Simple Hamiltonian time-dependent coefficients
    Complex cv_cqed::H1_coeff_simple(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[0] * t) - I * sin(nu_list[0] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H2_coeff_simple(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[0] * t) + I * sin(nu_list[0] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H3_coeff_simple(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[1] * t) - I * sin(nu_list[1] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H4_coeff_simple(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[1] * t) + I * sin(nu_list[1] * t)) * cos((delta + stark_shift[idx]) * t);
    }


    // Time-dependent coefficients for a 5 ion-crystal
    double cv_cqed::H0_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H1_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[0] * t) + I * sin(nu_list[0] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H2_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[1] * t) + I * sin(nu_list[1] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H3_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[2] * t) + I * sin(nu_list[2] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H4_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[3] * t) + I * sin(nu_list[3] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H5_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[4] * t) + I * sin(nu_list[4] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H6_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[0] * t) - I * sin(nu_list[0] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H7_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[1] * t) - I * sin(nu_list[1] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H8_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[2] * t) - I * sin(nu_list[2] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H9_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[3] * t) - I * sin(nu_list[3] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::H10_coeff_5(double t) {
        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return pulse[idx] * (cos(nu_list[4] * t) - I * sin(nu_list[4] * t)) * cos((delta + stark_shift[idx]) * t);
    }

    Complex cv_cqed::omega(double t){
        return cos(nu_list[1] * t) - I * sin(nu_list[1] * t);
    }

    Complex cv_cqed::omega_d(double t){
        return cos(nu_list[1] * t) + I * sin(nu_list[1] * t);
    }

    Complex cv_cqed::g_t(double t) {
        return cos(3 * nu_list[1] * t) - I * sin(3 * nu_list[1] * t);
    }

    // Lindblad operators time-dependent coefficients
    double cv_cqed::gamma_motion_dephasing_coeff(double t) {
        return 0.0;
    }

    double cv_cqed::gamma_spin_dephasing_coeff(double t) {

        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps) {
            idx -= 1;
        }

        return sqrt(fabs(pulse[idx]) * gamma_Rayleigh / 4.0);
    }

    double cv_cqed::gamma_spin_flip_coeff(double t) {

        int idx = (int) (t * (pulse_steps / t_gate));
        if (idx == pulse_steps){
            idx -= 1;
        }

        return sqrt(fabs(pulse[idx]) * gamma_Raman);
    }

    void cv_cqed::set_stark_shift() {
        /*
         * The input pulse is in units of Mrad/s
         * */

        stark_shift.resize(pulse_steps);
        if (stark_status) {
            for (int i = 0; i < pulse_steps; i++) {
                stark_shift[i] = stark_scale_factor * pulse[i];
            }
        } else {
            for (int i = 0; i < pulse_steps; i++) {
                stark_shift[i] = 0.0;
            }
        }
    }

    void cv_cqed::set_H(std::string hamiltonian) {
        if (hamiltonian == "TWO_ION_SIMPLE")
            set_H(Hamiltonian::TWO_ION_SIMPLE);
        else if (hamiltonian == "TWO_ION_FULL")
            set_H(Hamiltonian::TWO_ION_FULL);
        else if (hamiltonian == "FIVE_ION")
            set_H(Hamiltonian::FIVE_ION);
        else if (hamiltonian == "FIVE_ION_KERR")
            set_H(Hamiltonian::FIVE_ION_KERR);
        else if (hamiltonian == "SMSPDC")
            set_H(Hamiltonian::SMSPDC);
    }

    void cv_cqed::set_H(Hamiltonian hamiltonian) {

        switch (hamiltonian) {
            case TWO_ION_SIMPLE: {
                /*
                 * Time dependencies
                 * */

                // Two-ion simple Hamiltonian (Manning thesis)
                qsd::ComplexFunction H1_t_simple = std::bind(&cv_cqed::H1_coeff_simple, this, std::placeholders::_1);
                qsd::ComplexFunction H2_t_simple = std::bind(&cv_cqed::H2_coeff_simple, this, std::placeholders::_1);
                qsd::ComplexFunction H3_t_simple = std::bind(&cv_cqed::H3_coeff_simple, this, std::placeholders::_1);
                qsd::ComplexFunction H4_t_simple = std::bind(&cv_cqed::H4_coeff_simple, this, std::placeholders::_1);

                // Two-ion simple Hamiltonian (Manning thesis)
                H = a_1 * (eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * sx_2) * H1_t_simple     // H1
                    + ad_1 * (eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * sx_2) * H2_t_simple  // H2
                    + a_2 * (eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * sx_2) * H3_t_simple   // H3
                    + ad_2 * (eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * sx_2) * H4_t_simple; // H4

                break;
            }

            case TWO_ION_FULL: {
                /*
                * Time dependencies
                * */

                //  Two-ion full Hamiltonian
                qsd::RealFunction H0_t = std::bind(&cv_cqed::H0_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H1_t = std::bind(&cv_cqed::H1_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H2_t = std::bind(&cv_cqed::H2_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H3_t = std::bind(&cv_cqed::H3_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H4_t = std::bind(&cv_cqed::H4_coeff, this, std::placeholders::_1);
                qsd::RealFunction H5_t = std::bind(&cv_cqed::H5_coeff, this, std::placeholders::_1);
                qsd::RealFunction H6_t = std::bind(&cv_cqed::H6_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H7_t = std::bind(&cv_cqed::H7_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H8_t = std::bind(&cv_cqed::H8_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H9_t = std::bind(&cv_cqed::H9_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H10_t = std::bind(&cv_cqed::H10_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H11_t = std::bind(&cv_cqed::H11_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H12_t = std::bind(&cv_cqed::H12_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H13_t = std::bind(&cv_cqed::H13_coeff, this, std::placeholders::_1);
                qsd::ComplexFunction H14_t = std::bind(&cv_cqed::H14_coeff, this, std::placeholders::_1);

                // Two-ion full Hamiltonian
                H = 0.5 * (sx_1 + sx_2) * H0_t // H0
                    - 0.5 * a_1 * (eta_list[0][ion_pair[0]] * sy_1 + eta_list[0][ion_pair[1]] * sy_2) * H1_t  // H1
                    - 0.5 * ad_1 * (eta_list[0][ion_pair[0]] * sy_1 + eta_list[0][ion_pair[1]] * sy_2) * H2_t // H2
                    - 0.5 * a_2 * (eta_list[1][ion_pair[0]] * sy_1 + eta_list[1][ion_pair[1]] * sy_2) * H3_t  // H3
                    - 0.5 * ad_2 * (eta_list[1][ion_pair[0]] * sy_1 + eta_list[1][ion_pair[1]] * sy_2) * H4_t // H4
                    - 0.25 * (2 * N_1 + id_2) * (eta_list[0][ion_pair[0]] * eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * eta_list[0][ion_pair[1]] * sx_2) * H5_t  // H5
                    - 0.25 * (2 * N_2 + id_2) * (eta_list[1][ion_pair[0]] * eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * eta_list[1][ion_pair[1]] * sx_2) * H6_t  // H6
                    - 0.25 * a_1 * a_1 * (eta_list[0][ion_pair[0]] * eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * eta_list[0][ion_pair[1]] * sx_2) * H7_t       // H7
                    - 0.25 * ad_1 * ad_1 * (eta_list[0][ion_pair[0]] * eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * eta_list[0][ion_pair[1]] * sx_2) * H8_t     // H8
                    - 0.25 * a_2 * a_2 * (eta_list[1][ion_pair[0]] * eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * eta_list[1][ion_pair[1]] * sx_2) * H9_t       // H9
                    - 0.25 * ad_2 * ad_2 * (eta_list[1][ion_pair[0]] * eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * eta_list[1][ion_pair[1]] * sx_2) * H10_t    // H10
                    - 0.25 * (2 * eta_list[0][ion_pair[0]] * eta_list[1][ion_pair[0]] * sx_1 + 2 * eta_list[0][ion_pair[1]] * eta_list[1][ion_pair[1]] * sx_2) * a_1 * a_2 * H11_t      // H11
                    - 0.25 * (2 * eta_list[0][ion_pair[0]] * eta_list[1][ion_pair[0]] * sx_1 + 2 * eta_list[0][ion_pair[1]] * eta_list[1][ion_pair[1]] * sx_2) * a_1 * ad_2 * H12_t     // H12
                    - 0.25 * (2 * eta_list[0][ion_pair[0]] * eta_list[1][ion_pair[0]] * sx_1 + 2 * eta_list[0][ion_pair[1]] * eta_list[1][ion_pair[1]] * sx_2) * ad_1 * a_2 * H13_t     // H13
                    - 0.25 * (2 * eta_list[0][ion_pair[0]] * eta_list[1][ion_pair[0]] * sx_1 + 2 * eta_list[0][ion_pair[1]] * eta_list[1][ion_pair[1]] * sx_2) * ad_1 * ad_2 * H14_t;   // H14

                break;
            }
            case FIVE_ION: {
                /*
                 * Time dependencies
                 * */

                // Five-ion simple Hamiltonian (Overleaf: https://www.overleaf.com/2163481312xkrbdkmybfks)
                qsd::ComplexFunction H0_t5 = std::bind(&cv_cqed::H0_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H1_t5 = std::bind(&cv_cqed::H1_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H2_t5 = std::bind(&cv_cqed::H2_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H3_t5 = std::bind(&cv_cqed::H3_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H4_t5 = std::bind(&cv_cqed::H4_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H5_t5 = std::bind(&cv_cqed::H5_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H6_t5 = std::bind(&cv_cqed::H6_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H7_t5 = std::bind(&cv_cqed::H7_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H8_t5 = std::bind(&cv_cqed::H8_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H9_t5 = std::bind(&cv_cqed::H9_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H10_t5 = std::bind(&cv_cqed::H10_coeff_5, this, std::placeholders::_1);

                // Five-ion simple Hamiltonian (Manning thesis)
                H = (1.0 * sy_1 + 1.0 * sy_2) * H0_t5 // H0
                    + ad_1 * (eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * sx_2) * H1_t5  // H1
                    + ad_2 * (eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * sx_2) * H2_t5  // H2
                    + ad_3 * (eta_list[2][ion_pair[0]] * sx_1 + eta_list[2][ion_pair[1]] * sx_2) * H3_t5  // H3
                    + ad_4 * (eta_list[3][ion_pair[0]] * sx_1 + eta_list[3][ion_pair[1]] * sx_2) * H4_t5  // H4
                    + ad_5 * (eta_list[4][ion_pair[0]] * sx_1 + eta_list[4][ion_pair[1]] * sx_2) * H5_t5  // H5
                    + a_1 * (eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * sx_2) * H6_t5   // H6
                    + a_2 * (eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * sx_2) * H7_t5   // H7
                    + a_3 * (eta_list[2][ion_pair[0]] * sx_1 + eta_list[2][ion_pair[1]] * sx_2) * H8_t5   // H8
                    + a_4 * (eta_list[3][ion_pair[0]] * sx_1 + eta_list[3][ion_pair[1]] * sx_2) * H9_t5   // H9
                    + a_5 * (eta_list[4][ion_pair[0]] * sx_1 + eta_list[4][ion_pair[1]] * sx_2) * H10_t5; // H10

                break;
            }

            case FIVE_ION_KERR: {
                /*
                 * Time dependencies
                 * */

                // Five-ion simple Hamiltonian (Overleaf: https://www.overleaf.com/2163481312xkrbdkmybfks)
                qsd::ComplexFunction H0_t5 = std::bind(&cv_cqed::H0_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H1_t5 = std::bind(&cv_cqed::H1_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H2_t5 = std::bind(&cv_cqed::H2_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H3_t5 = std::bind(&cv_cqed::H3_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H4_t5 = std::bind(&cv_cqed::H4_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H5_t5 = std::bind(&cv_cqed::H5_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H6_t5 = std::bind(&cv_cqed::H6_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H7_t5 = std::bind(&cv_cqed::H7_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H8_t5 = std::bind(&cv_cqed::H8_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H9_t5 = std::bind(&cv_cqed::H9_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H10_t5 = std::bind(&cv_cqed::H10_coeff_5, this, std::placeholders::_1);

                // Five-ion simple Hamiltonian (Manning thesis)
                H = (1.0 * sy_1 + 1.0 * sy_2) * H0_t5 // H0
                    + ad_1 * (eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * sx_2) * H1_t5  // H1
                    + ad_2 * (eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * sx_2) * H2_t5  // H2
                    + ad_3 * (eta_list[2][ion_pair[0]] * sx_1 + eta_list[2][ion_pair[1]] * sx_2) * H3_t5  // H3
                    + ad_4 * (eta_list[3][ion_pair[0]] * sx_1 + eta_list[3][ion_pair[1]] * sx_2) * H4_t5  // H4
                    + ad_5 * (eta_list[4][ion_pair[0]] * sx_1 + eta_list[4][ion_pair[1]] * sx_2) * H5_t5  // H5
                    + a_1 * (eta_list[0][ion_pair[0]] * sx_1 + eta_list[0][ion_pair[1]] * sx_2) * H6_t5   // H6
                    + a_2 * (eta_list[1][ion_pair[0]] * sx_1 + eta_list[1][ion_pair[1]] * sx_2) * H7_t5   // H7
                    + a_3 * (eta_list[2][ion_pair[0]] * sx_1 + eta_list[2][ion_pair[1]] * sx_2) * H8_t5   // H8
                    + a_4 * (eta_list[3][ion_pair[0]] * sx_1 + eta_list[3][ion_pair[1]] * sx_2) * H9_t5   // H9
                    + a_5 * (eta_list[4][ion_pair[0]] * sx_1 + eta_list[4][ion_pair[1]] * sx_2) * H10_t5  // H10
                    // + N_1 * (N_2 + N_3 + N_4 + N_5)
                    + TWOPI * 3e-4 * N_2 * (N_3 + N_4 + N_5)
                    + TWOPI * 3e-4 * N_3 * (N_4 + N_5)
                    + TWOPI * 3e-4 * N_4 * (N_5);

                break;
            }

            case SMSPDC: {
                // Single mode parametric spontaneous down conversion
                qsd::ComplexFunction G_t = std::bind(&cv_cqed::g_t, this, std::placeholders::_1);
                qsd::ComplexFunction w = std::bind(&cv_cqed::omega, this, std::placeholders::_1);
                qsd::ComplexFunction w_d = std::bind(&cv_cqed::omega_d, this, std::placeholders::_1);


                H = (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * G_t * 1;//(ad_1 * w_d + a_1 * w) * G_t * 0.0001 // first term

//                    + (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * (1.0) * 0.1
//                    + (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * (ad_1 * w_d + a_1 * w) * (-1) * (G_t) * 0.1;

            }
        }
    }

    void cv_cqed::set_state(States states) {
        /*
        * The initial state
        * */

        std::vector<int>().swap(dynamic_degrees);

        std::vector<qsd::State> psilist;
        for (int i = 0; i < nu_list.size() + 2; i++) {
            if (i < 2) psilist.emplace_back(2, qsd::SPIN);
            else {
                dynamic_degrees.emplace_back(i);
                psilist.emplace_back(phonon_cutoffs[i - 2], states.motion[i - 2], qsd::FIELD);
            }
        }

        psi0 = qsd::State(psilist.size(), psilist.data());

        if (states.spin[0] == 1)
            psi0 *= sp_1;
        if (states.spin[1] == 1)
            psi0 *= sp_2;
    }

    void cv_cqed::set_L(std::string lindblads)
    {
        if (lindblads == "NONE") {
            set_L(Lindblads::NONE);
        } else if (lindblads == "FULL") {
            set_L(Lindblads::FULL);
        }
    }

    void cv_cqed::set_L(Lindblads lindblads) {

        /*
        * The Lindblad operators
        *
        * Example:
        *  const int nOfLindblads = 1;
        *  Operator L[nOfLindblads] = {sz_1};
        * */

        /*
         * Section 4.4, Christopher J. Ballance thesis
         */

        switch (lindblads) {
            case FULL: {

                nL = 38;
                L = new qsd::Operator[nL];

                /*
                * Time dependencies
                * */

                qsd::RealFunction L_md_t = std::bind(&cv_cqed::gamma_motion_dephasing_coeff, this, std::placeholders::_1);
                qsd::RealFunction L_sd_t = std::bind(&cv_cqed::gamma_spin_dephasing_coeff, this, std::placeholders::_1); // gamma_Rayleigh
                qsd::RealFunction L_sf_t = std::bind(&cv_cqed::gamma_spin_flip_coeff, this, std::placeholders::_1); // gamma_Raman

                //Power laser fluctuations
                qsd::ComplexFunction H0_t5 = std::bind(&cv_cqed::H0_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H1_t5 = std::bind(&cv_cqed::H1_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H2_t5 = std::bind(&cv_cqed::H2_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H3_t5 = std::bind(&cv_cqed::H3_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H4_t5 = std::bind(&cv_cqed::H4_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H5_t5 = std::bind(&cv_cqed::H5_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H6_t5 = std::bind(&cv_cqed::H6_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H7_t5 = std::bind(&cv_cqed::H7_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H8_t5 = std::bind(&cv_cqed::H8_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H9_t5 = std::bind(&cv_cqed::H9_coeff_5, this, std::placeholders::_1);
                qsd::ComplexFunction H10_t5 = std::bind(&cv_cqed::H10_coeff_5, this, std::placeholders::_1);

                /*
                * Lindblad operators
                **/

                int nL_cnt = 0;

                // Motional Heating
                L[nL_cnt++] = sqrt(gamma_heating) * a_1;  // Heating of the 1st motional mode
                L[nL_cnt++] = sqrt(gamma_heating) * ad_1; // ...

                L[nL_cnt++] = sqrt(gamma_heating) * a_2;  // Heating of the 2nd motional mode
                L[nL_cnt++] = sqrt(gamma_heating) * ad_2; // ...

                L[nL_cnt++] = sqrt(gamma_heating) * a_3;  // Heating of the 3nd motional mode
                L[nL_cnt++] = sqrt(gamma_heating) * ad_3; // ...

                L[nL_cnt++] = sqrt(gamma_heating) * a_4;  // Heating of the 4nd motional mode
                L[nL_cnt++] = sqrt(gamma_heating) * ad_4; // ...

                L[nL_cnt++] = sqrt(gamma_heating) * a_5;  // Heating of the 5nd motional mode
                L[nL_cnt++] = sqrt(gamma_heating) * ad_5; // ...

                // Motional Dephasing
                // L[] = L_md_t * (ad_1 * a_1); // Dephasing of the 1st motional mode
                // L[] = L_md_t * (ad_2 * a_2); // Dephasing of the 2nd motional mode

                // Spin Dephasing
                L[nL_cnt++] = sz_1 * L_sd_t;  // Dephasing of the 1st ion
                L[nL_cnt++] = sz_2 * L_sd_t;  // Dephasing of the 2nd ion

                // Spin Flip
                L[nL_cnt++] =  sm_1 * L_sf_t;  // Spin Flip of the 1st ion
                L[nL_cnt++] =  sp_1 * L_sf_t;  // ...

                L[nL_cnt++] = sm_2 * L_sf_t;  // Spin Flip of the 2st ion
                L[nL_cnt++] = sp_2 * L_sf_t;  // ...

                L[nL_cnt++] = amp_fluct * sy_1 * H0_t5;
                L[nL_cnt++] = amp_fluct * sy_2 * H0_t5;
                L[nL_cnt++] = amp_fluct * ad_1 * sx_1 * H1_t5;
                L[nL_cnt++] = amp_fluct * ad_1 * sx_2 * H1_t5;
                L[nL_cnt++] = amp_fluct * ad_2 * sx_1 * H2_t5;
                L[nL_cnt++] = amp_fluct * ad_2 * sx_2 * H2_t5;
                L[nL_cnt++] = amp_fluct * ad_3 * sx_1 * H3_t5;
                L[nL_cnt++] = amp_fluct * ad_3 * sx_2 * H3_t5;
                L[nL_cnt++] = amp_fluct * ad_4 * sx_1 * H4_t5;
                L[nL_cnt++] = amp_fluct * ad_4 * sx_2 * H4_t5;
                L[nL_cnt++] = amp_fluct * ad_5 * sx_1 * H5_t5;
                L[nL_cnt++] = amp_fluct * ad_5 * sx_2 * H5_t5;

                L[nL_cnt++] = amp_fluct * a_1 * sx_1 * H6_t5;
                L[nL_cnt++] = amp_fluct * a_1 * sx_2 * H6_t5;
                L[nL_cnt++] = amp_fluct * a_2 * sx_1 * H7_t5;
                L[nL_cnt++] = amp_fluct * a_2 * sx_2 * H7_t5;
                L[nL_cnt++] = amp_fluct * a_3 * sx_1 * H8_t5;
                L[nL_cnt++] = amp_fluct * a_3 * sx_2 * H8_t5;
                L[nL_cnt++] = amp_fluct * a_4 * sx_1 * H9_t5;
                L[nL_cnt++] = amp_fluct * a_4 * sx_2 * H9_t5;
                L[nL_cnt++] = amp_fluct * a_5 * sx_1 * H10_t5;
                L[nL_cnt++] = amp_fluct * a_5 * sx_2 * H10_t5;

                break;
            }

            case NONE: {
                nL = 0;
                break;
            }
        }
    }

    void cv_cqed::set_flist(std::vector<string> _flist) { flist = _flist; }
    void cv_cqed::set_outlist(std::vector<qsd::Operator> _outlist) { outlist = _outlist; }

    void cv_cqed::set_processor(std::string type)
    {
        if (type == "FIDELITY")
            set_processor(Processor::FIDELITY);
        else if (type == "INFIDELITY")
            set_processor(Processor::INFIDELITY);
        else if (type == "PARITY")
            set_processor(Processor::PARITY);
    }

    void cv_cqed::set_processor(Processor type) { processor_type = type; }

    double cv_cqed::processor(qsd::Result result) {
        switch (processor_type) {
            case FIDELITY: {
                return std::max(result.data[0].back()[0], result.data[1].back()[0]);
            }

            case INFIDELITY: {
                return fabs(1.0 - std::max(result.data[0].back()[0], result.data[1].back()[0]));
            }

            case PARITY: {

                /*
                 * \rho_mn,pq = <mn| \rho |pq>
                 * Parity = (\rho_00,00 + \rho_11,11) - (\rho_01,01 + \rho_10,10)
                 */

                /*
                 * The main formula for calculation of fidelity
                 * fid = 0.5 * (_rho_00_00.real() + _rho_11_11.real() - 2 * abs(_rho_00_11) * sin(-arg(_rho_00_11))); // sin(arg(_rho_00_11))) gives the other state
                 *
                 * The formula for calculation of fidelity based on the contrast of parity curve
                 * fid = 0.5 * (_rho_00_00.real() + _rho_11_11.real() + 2 * abs(_rho_00_11));
                 *
                 * The formula for calculation of fidelity based on the contrast of parity curve that is implemented in the code
                 * fid = 0.5 * (_rho_00_00.real() + _rho_11_11.real() + two_abs_rho_00_11);
                 */

                // Open the parity file
                std::vector<ofstream> parity_files;

                for (int i = 0; i < result.state.size(); i++) {
                    if (result.state.size() == 1) {
                        parity_files.emplace_back(std::ofstream{"parity.txt"});
                    } else {
                        parity_files.emplace_back(std::ofstream{"parity_traj_" + std::to_string(i + 1) + ".txt"});
                    }
                    parity_files[i] << "phi, parity" << std::endl;
                }

                std::vector<double> fidelities;

                for (int i = 0; i < result.state.size(); i++) {

                    ///////////////////////////////////
                    qsd::State psi_1 = result.state[i];
                    qsd::State psi_2 = result.state[i];

                    // \rho_00,00
                    psi_2 = psi_1;

                    psi_2 *= p00_1;
                    psi_2 *= p00_2;

                    Complex _rho_00_00 = psi_1 * psi_2;

                    // \rho_11,11
                    psi_2 = psi_1;

                    psi_2 *= p11_1;
                    psi_2 *= p11_2;

                    Complex _rho_11_11 = psi_1 * psi_2;
                    ///////////////////////////////////

                    // Global PI/2 rotation (Analyzer pulse) - Manning thesis page 110
                    double phi = 0.0; // Analyzer pulse phase

                    std::vector<double> parity;
                    while (phi < TWOPI) {

                        qsd::State psi_1 = result.state[i];
                        qsd::State psi_2 = result.state[i];

                        qsd::Operator rg_1 = cos(PIHALF / 2.0) * id_s1
                                        - I * sin(PIHALF / 2.0) * (cos(phi) * sx_1 + sin(phi) * sy_1);
                        qsd::Operator rg_2 = cos(PIHALF / 2.0) * id_s2
                                        - I * sin(PIHALF / 2.0) * (cos(phi) * sx_2 + sin(phi) * sy_2);

                        psi_1 *= rg_1;
                        psi_1 *= rg_2;

                        // \rho_00,00
                        psi_2 = psi_1;

                        psi_2 *= p00_1;
                        psi_2 *= p00_2;

                        Complex _rho_00_00 = psi_1 * psi_2;

                        // \rho_11,11
                        psi_2 = psi_1;

                        psi_2 *= p11_1;
                        psi_2 *= p11_2;

                        Complex _rho_11_11 = psi_1 * psi_2;

                        // \rho_00,11
                        psi_2 = psi_1;

                        psi_2 *= p01_1;
                        psi_2 *= p01_2;

                        Complex _rho_00_11 = psi_1 * psi_2;

                        // \rho_11,00
                        psi_2 = psi_1;

                        psi_2 *= p10_1;
                        psi_2 *= p10_2;

                        Complex _rho_11_00 = psi_1 * psi_2;

                        // \rho_01,01
                        psi_2 = psi_1;

                        psi_2 *= p00_1;
                        psi_2 *= p11_2;

                        Complex _rho_01_01 = psi_1 * psi_2;

                        // \rho_10,10
                        psi_2 = psi_1;

                        psi_2 *= p11_1;
                        psi_2 *= p00_2;

                        Complex _rho_10_10 = psi_1 * psi_2;

                        double _parity = real(_rho_00_00 + _rho_11_11) - real(_rho_01_01 + _rho_10_10);

                        parity_files[i] << phi << ", " << _parity << std::endl;
                        parity.push_back(_parity);

                        phi += phi_res;
                    }

                    auto max_parity = max_element(std::begin(parity), std::end(parity));
                    auto min_parity = min_element(std::begin(parity), std::end(parity));
                    double two_abs_rho_00_11 = (*max_parity - *min_parity) / 2;

                    fidelities.push_back(0.5 * (_rho_00_00.real() + _rho_11_11.real() + two_abs_rho_00_11));
                }

                return accumulate(fidelities.begin(), fidelities.end(), 0.0) / fidelities.size();
            }

            default: {
                return 0.0;
            }
        }
    }

    void cv_cqed::set_phi_res(double _phi_res)
    {
        phi_res = _phi_res;
    }

    // Utility methods
    void cv_cqed::load_1d_list(std::string filename, std::vector<double>& _list)
    {
        // Clear the delta_list
        std::vector<double>().swap(_list);

        // load the data file into a new array
        cnpy::NpyArray npy_data = cnpy::npy_load(filename);
        auto data_ptr = npy_data.data<double>();

        for(auto i=0; i< npy_data.num_vals; i++){
            _list.push_back(data_ptr[i]);
        }
    }

    void cv_cqed::load_2d_list(std::string filename, std::vector<std::vector<double>>& _list)
    {
        // Clear the pulse_list
        std::vector<std::vector<double>>().swap(_list);

        // load the data file into a new array
        cnpy::NpyArray npy_data = cnpy::npy_load(filename);
        auto data_ptr = npy_data.data<double>();

        int nrows = npy_data.shape[0];
        int ncols = npy_data.shape[1];

        _list.reserve(nrows);
        for(int row = 0; row < nrows; row++) {
            _list.emplace_back(ncols);
            for(int col = 0; col < ncols; col++) {
                _list[row][col] = data_ptr[ncols * row + col];
            }
        }
    }

}