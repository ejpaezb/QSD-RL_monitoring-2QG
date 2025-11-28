
//
// Created by Shakib Vedaie on 2019-09-20.
//

#ifndef IONTRAP_ION_TRAP_H
#define IONTRAP_ION_TRAP_H

#include <map>
#include <cmath>
#include <random>
#include <complex>
#include <iomanip>

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
#include <boost/math/constants/constants.hpp>

#include <qsd/qsd_experiment.h>

#include "h_config.cpp"
#include "l_config.cpp"
#include "state_config.cpp"

#include "utilities.h"

namespace experiment {
/*
 * Ion-trap simulator.
 * */

class ion_trap : public qsd::qsd_experiment {
public:
    ion_trap();
    ion_trap(double _dt, int _numdts, int _numsteps, double _delta, std::vector<double> &_pulse);
    ion_trap(const std::string &filename);
	ion_trap(ion_trap &_ion_trap);

    struct PulseCfg {

        int profile_am;
		int addressing_am;
		int steps_am;

	  	int profile_fm;
	  	int addressing_fm;
	  	int steps_fm;

	  	bool SYM_AM;
        bool SYM_FM;

	  	double scale_am;
	  	double scale_fm;
    };

	struct ThermalState {

	  std::vector<std::pair<double, std::vector<int>>> states;

	  double cutoff_probability;
	  std::vector<double> n_bar;
	  std::vector<int> max_n;
	  int n_states;
	  int n_states_cutoff;
	};

    HConfig h_cfg;
    qsd::Operator H;
	qsd::Operator H_fluc;
    qsd::Operator H_fluc_1;

    qsd::State psi0;
    StateConfig state_cfg;
    
    std::vector<qsd::Operator> L;

    std::vector<std::string> flist;      // Output files

	std::map<std::string, qsd::Operator> outlist;  // Operators to output
	std::vector<std::string> outlist_string;

    double t_gate;

    std::string name;
    std::string description;

    LConfig l_cfg;

    int processor_type;
    int simulation_mode;

    std::vector<std::vector<double>> pulse_list;
    std::vector<std::vector<double>>  shelving_pulse_list;
    std::vector<std::vector<double>> eta_list;      // Lamb-Dicke parameter for the ions coupling to the motional modes
    std::vector<double> delta_list;
    std::vector<double> nu_list;                    // Motional mode frequencies [Mrad/s]
    std::vector<std::vector<double>> shaving_phase_list;

    double nu_error; // Detuning error in the motional-mode frequencies [Mrad/s]

	PulseCfg pulse_cfg;
	std::vector<std::vector<double>> g_data;

	std::vector<std::vector<double>> pulse;
    std::vector<std::vector<double>> shelving_pulse;
    int num_pulses;
    int shelving_flag;
    std::vector<int> shelving_idx;
    std::vector<double> shelving_time;
    double energy_drift;
    std::vector<std::vector<double>> EIT_pulse;
    std::vector<std::vector<double>> shaving_phase;

    std::vector<double> Delta_pump;

    std::vector<double> pulse_sampled;

	std::vector<std::vector<double>> delta;
    std::string pulse_sequence;
    std::string prep_sequence;
    int prep_state;

	std::vector<double> rabi_freq_ratio;

    double phi_res;
    double Rabi_s;
    double Delta_L;

	void initialize_expt(unsigned int seed);

    void parse_config(const std::string &config_path);

    qsd::Operator outlist_processor(const std::string &_operator);

    void set_delta(std::vector<std::vector<double>> _delta);
    void set_delta(std::vector<double> _delta);
    void set_delta(double _delta);
    void set_delta(int _delta_idx);
	void set_external_delta(boost::property_tree::ptree external_delta_cfg, std::string select);
    void set_shelving_phase(std::vector<double> &_phas);
    void set_shelving_phase(int _phase_idx);

    std::vector<double> get_delta();
	double get_delta_new(double t);
    double get_shelving_phase(double t);
    double get_energy_drift(int idx, double t);
    double get_coloured_noise(double t, int idx);

    void set_pulse(std::vector<double> &_pulse);
    void set_pulse(std::vector<std::vector<double>> &_pulse);
    void set_pulse(int _pulse_idx);
    void set_EIT_pulse(std::vector<double> _pulse);
    void set_Delta_pump(double det);
    void set_shelving_pulse(std::vector<double> &_pulse);
    void set_shelving_pulse(int _pulse_idx);

	void set_pulse_external(boost::property_tree::ptree external_pulse_cfg, std::string select);

    void set_time(double _dt, int _numdts, int _numsteps);
    void set_time(double t);
    
    void set_phonon_cutoffs();

    void initialize_H();
    void initialize_L();

    void set_state();
    void set_spin_state(std::vector<qsd::State> spin_state);

    void set_outlist(std::vector<std::string> _outlist_string);

    double processor(qsd::State qsd_state);

	qsd::Expectation processor_fidelity(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);
	qsd::Expectation processor_infidelity(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);
	std::vector<double> processor_parity(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);
	qsd::TrajectoryResult processor_average_trajectory(std::vector<qsd::TrajectoryResult> qsd_result, std::vector<double> weights=std::vector<double>(), bool save_per_traj=false, std::ofstream* log_file=NULL);
	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> processor_concatenated_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> processor_reduced_concatenated_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> processor_concatenated_multilevel_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result,int atom, bool save_per_traj=false, std::ofstream* log_file=NULL);
    std::vector<std::vector<std::vector<qsd::Expectation>>> processor_average_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);
    std::vector<std::vector<std::vector<qsd::Expectation>>> processor_average_multilevel_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, int atom, bool save_per_traj=false, std::ofstream* log_file=NULL);
    std::vector<std::vector<std::vector<qsd::Expectation>>> processor_reduced_average_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);
    qsd::Expectation processor_even_population(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj=false, std::ofstream* log_file=NULL);

	ion_trap::ThermalState prepare_thermal_state();

	double get_g(int idx, double t);
	double get_g_fluc(int idx, double t);

    double get_omega(int idx, double t);
    double get_EIT_pulse(int idx);
    double get_Delta_pump(int idx);
    double get_stark_shift(int idx, double t);

    double get_shelving_pulse(int idx, double t, int flag);
    int get_shelving_idx(int idx, double t);

    bool findParity(int x);

    int get_n_ions();
    int get_n_spins();

	double h_cos(int idx, double t);

    double dt;      // basic time step
    int numdts;     // time interval between outputs = numdts * dt
    int numsteps;   // total integration time = numsteps * numdts * dt

    std::complex<double> IM {0.0, 1.0}; // for operations with complex/real functions
    std::complex<double> mIM{0, -1};
    const qsd::ImaginaryUnit Im = qsd::Im; // for explicit operations with operators (not states)
    const qsd::ImaginaryUnit mIm = qsd::mIm;

    void initialize_operators();

    /*
     * Operators
     * */

    std::vector<qsd::SigmaX> sx;
    std::vector<qsd::SigmaY> sy;
    std::vector<qsd::SigmaZ> sz;

    std::vector<qsd::SigmaPlus> sp;
    std::vector<qsd::Operator> sm;

    std::vector<qsd::XOperator> x;
    std::vector<qsd::POperator> p;

    std::vector<qsd::AnnihilationOperator> a;
    std::vector<qsd::Operator> ad;

    std::vector<qsd::NumberOperator> N;

    qsd::Operator id;
//    Operators for EIT cooling dynamics.
    std::vector<qsd::TransitionOperator> S_01;
    std::vector<qsd::TransitionOperator> S_13;
    std::vector<qsd::TransitionOperator> S_03;
    std::vector<qsd::TransitionOperator> S_23;
    std::vector<qsd::TransitionOperator> S_00;
    std::vector<qsd::TransitionOperator> S_11;
    std::vector<qsd::TransitionOperator> S_22;
    std::vector<qsd::TransitionOperator> S_33;

        //    Operators for shelving/deshelving fidelities
    std::vector<qsd::TransitionOperator> S0_D0; //Equivalent to a Sigma^+
    std::vector<qsd::TransitionOperator> D0_S0; //Equivalent to a Sigma^-
    std::vector<qsd::TransitionOperator> S0_D1;
    std::vector<qsd::TransitionOperator> S1_D0;
    std::vector<qsd::TransitionOperator> S1_D1;
    std::vector<qsd::TransitionOperator> D1_S1;
    std::vector<qsd::TransitionOperator> S0_S1;
    std::vector<qsd::TransitionOperator> S1_S0;
    std::vector<qsd::TransitionOperator> Y_S0_D0; //Equivalent to a Sigma^y
    std::vector<qsd::TransitionOperator> X_S0_D0; //Equivalent to a Sigma^x
    std::vector<qsd::TransitionOperator> X_S1_D1;
    std::vector<qsd::TransitionOperator> Y_S1_D1;
    std::vector<qsd::TransitionOperator> Y_S0_S1;
    std::vector<qsd::TransitionOperator> pulse_seq_0;
    std::vector<qsd::TransitionOperator> pulse_seq_1;
    std::vector<qsd::TransitionOperator> S0;
    std::vector<qsd::TransitionOperator> S1;
    std::vector<qsd::TransitionOperator> D0;
    std::vector<qsd::TransitionOperator> D1;
    std::vector<qsd::TransitionOperator> X_30_0;
    std::vector<qsd::TransitionOperator> X_60_0;
    std::vector<qsd::TransitionOperator> Y_30_0;
    std::vector<qsd::TransitionOperator> Y_60_0;
    std::vector<qsd::TransitionOperator> X_30_1;
    std::vector<qsd::TransitionOperator> X_60_1;
    std::vector<qsd::TransitionOperator> Y_30_1;
    std::vector<qsd::TransitionOperator> Y_60_1;

    std::vector<qsd::IdentityOperator> id_s;
    std::vector<qsd::IdentityOperator> id_m;

    std::vector<std::vector<qsd::Operator>> projectors_z;
    std::vector<std::vector<qsd::Operator>> projectors_x;
    std::vector<qsd::TransitionOperator> proj_z0_basis;
    std::vector<qsd::TransitionOperator> proj_z1_basis;
    std::vector<qsd::Operator> joint_x;
    std::vector<qsd::TransitionOperator> proj_z3_basis;
    std::vector<qsd::TransitionOperator> sigma_X;
    std::vector<qsd::TransitionOperator> sigma_Y;
    std::vector<qsd::Operator> joint_p;

    std::vector<qsd::Operator> rho_ideals;

private:
    };

}
#endif //IONTRAP_ION_TRAP_H
