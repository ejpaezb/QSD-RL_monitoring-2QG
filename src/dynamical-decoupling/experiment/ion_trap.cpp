//
// Created by Shakib Vedaie on 2019-09-20.
//

#include "ion_trap.h"
//#define pi = 3.1416;

namespace experiment {
/*
 * Ion-trap simulator.
 * */

ion_trap::ion_trap() {}

ion_trap::ion_trap(const std::string &filename) { parse_config(filename); }

ion_trap::ion_trap(ion_trap &_ion_trap) {

	// Read the experiment information
	this->name = _ion_trap.name;
	this->description = _ion_trap.description;
    this->simulation_mode = _ion_trap.simulation_mode;

	// Read the detuning list and the pulse list
	this->delta_list = _ion_trap.delta_list;
    this->pulse_list = _ion_trap.pulse_list;
    this->shelving_pulse_list = _ion_trap.shelving_pulse_list;

    this->pulse_sequence = _ion_trap.pulse_sequence;
    this->prep_sequence = _ion_trap.prep_sequence;
    this->prep_state = _ion_trap.prep_state;

	// Read the trap parameters
	this->nu_list = _ion_trap.nu_list;
	this->eta_list = _ion_trap.eta_list;

	// Read the Hamiltonian
	this->h_cfg = _ion_trap.h_cfg;

	// Read the Lindblads
	this->l_cfg = _ion_trap.l_cfg;
//    this->mode_cfg = _ion_trap.mode_cfg;

	// Read the state
	this->state_cfg = _ion_trap.state_cfg;

	// Read the pulse config
	// AM & FM
	this->pulse_cfg = _ion_trap.pulse_cfg;

	// Initialize the operators
	this->initialize_operators();

	// Read the outlist
	this->outlist_string = _ion_trap.outlist_string;
	this->set_outlist(this->outlist_string);

	// Read the flist
	this->flist = _ion_trap.flist;

	// Read the time
	this->dt = _ion_trap.dt;
	this->numdts = _ion_trap.numdts;
	this->numsteps = _ion_trap.t_gate;

	this->t_gate = _ion_trap.t_gate;
    this->shelving_flag = _ion_trap.shelving_flag;

	// Initialize the Hamiltonian and the Lindblads
	this->initialize_H();
	this->initialize_L();

	// Set delta and pulse

    this->delta = _ion_trap.delta;
    this->shaving_phase_list = _ion_trap.shaving_phase_list;

    this->shaving_phase = _ion_trap.shaving_phase;

    this->set_phonon_cutoffs();

	this->pulse = _ion_trap.pulse;
    this->EIT_pulse = _ion_trap.EIT_pulse;
    this->shelving_pulse = _ion_trap.shelving_pulse;
    this->num_pulses = _ion_trap.num_pulses;
    this->shelving_idx = _ion_trap.shelving_idx;
    this->shelving_time = _ion_trap.shelving_time;
    this->energy_drift = _ion_trap.energy_drift;
    this->Delta_pump = _ion_trap.Delta_pump;
	this->pulse_sampled = _ion_trap.pulse_sampled;
    this-> Rabi_s = _ion_trap.Rabi_s;
    this-> Delta_L = _ion_trap.Delta_L;

	this->rabi_freq_ratio = _ion_trap.rabi_freq_ratio;

	this->g_data = _ion_trap.g_data;


}

void ion_trap::initialize_expt(unsigned int seed) {
	// Not implemented!
}

void ion_trap::parse_config(const std::string &config_path){

	// Create empty property tree object
	boost::property_tree::ptree trap_cfg;
	boost::property_tree::ptree solutions_cfg;
	boost::property_tree::ptree experiment_cfg;
	boost::property_tree::ptree interaction_cfg;
    boost::property_tree::ptree mode_cfg;

    // Parse the INFO into the property tree.
    boost::property_tree::read_info(config_path + "trap.cfg", trap_cfg);
	boost::property_tree::read_info(config_path + "solutions.cfg", solutions_cfg);
	boost::property_tree::read_info(config_path + "experiment.cfg", experiment_cfg);
	boost::property_tree::read_info(config_path + "interaction.cfg", interaction_cfg);

    simulation_mode = experiment_cfg.get<int>("experiment.simulation_mode");

    if (simulation_mode==8){
        boost::property_tree::read_info(config_path + "mode-8.cfg", mode_cfg);
        load_2d_list(mode_cfg.get<std::string>("mode-8.phase_list_path"), shaving_phase_list);
        load_2d_list(mode_cfg.get<std::string>("mode-8.pulse_list_path"), shelving_pulse_list);
        pulse_sequence = mode_cfg.get<std::string>("mode-8.pulse_sequence");
        prep_state = mode_cfg.get<int>("mode-8.prep_state");

    } else if (simulation_mode==7){
        boost::property_tree::read_info(config_path + "mode-7.cfg", mode_cfg);
    }

    // Read the experiment information
    name = experiment_cfg.get<std::string>("experiment.name");
    description = experiment_cfg.get<std::string>("experiment.description");
//    Delta_pump = mode_cfg.get<double>("mode-7.Delta_pump");

    // Read the detuning list and the pulse list
    load_1d_list(solutions_cfg.get<std::string>("solutions.internal.delta_list_path"), delta_list);
    load_2d_list(solutions_cfg.get<std::string>("solutions.internal.pulse_list_path"), pulse_list);

	// Read the trap parameters
	int nu_list_mode = trap_cfg.get<int>("trap.nu_list.mode");
	int nu_list_unit = trap_cfg.get<int>("trap.nu_list.unit");

	switch (nu_list_mode) {

		case 0: { // Manual

			std::vector<double> nu_list_manual;

			std::stringstream string_stream(trap_cfg.get<std::string>("trap.nu_list.nu_list_manual"));
			while(string_stream.good())
			{
				std::string substr;
				std::getline(string_stream, substr, ',');
				if (substr != "")
					nu_list_manual.push_back(std::stod(substr));
			}

			nu_list = nu_list_manual;

			break;
		}

		case 1: { // File

			load_1d_list(trap_cfg.get<std::string>("trap.nu_list.nu_list_path"), nu_list);
			break;
		}

		default: {
			break;
		}
	}

	for (int k = 0; k < nu_list.size(); k++) {
		if (nu_list_unit == 0) { // MHz
			using namespace boost::math::double_constants;
			nu_list[k] = two_pi * nu_list[k];
		}
	}


	int eta_list_mode = trap_cfg.get<int>("trap.eta_list.mode");
	switch (eta_list_mode) {

		case 0: { // Manual

			std::vector<std::vector<double>> eta_list_manual;

			for (auto& item : trap_cfg.get_child("trap.eta_list.eta_list_manual")) {

				std::vector<double> eta_list_manual_k;

				std::stringstream string_stream(item.second.get_value<std::string>());
				while(string_stream.good())
				{
					std::string substr;
					std::getline(string_stream, substr, ',');
					if (substr != "")
						eta_list_manual_k.push_back(std::stod(substr));
				}

				eta_list_manual.push_back(eta_list_manual_k);
			}

			eta_list = eta_list_manual;

			break;
		}

		case 1: { // File

			load_2d_list(trap_cfg.get<std::string>("trap.eta_list.eta_list_path"), eta_list);
			break;
		}

		default: {
			break;
		}
	}

	// Apply a Gaussian noise to the motional-mode frequencies [Mrad/s]
	std::random_device rd{};
	std::mt19937 gen{rd()};

	double rnd_std = trap_cfg.get<double>("trap.nu_std");
	for (int i = 0; i < nu_list.size(); i++) {

		double rnd_mean = nu_list[i];
		std::normal_distribution<> d{rnd_mean, rnd_std};

		double rnd_nu = d(gen);
		nu_list[i] = rnd_nu;

		// // std::cout << nu_list[i] << std::endl;
	}

    // Read and apply the detuning error in the motional-mode frequencies [Mrad/s]
    nu_error = trap_cfg.get<double>("trap.nu_error");
    for (int i = 0; i < nu_list.size(); i++) {
		nu_list[i] += nu_error;
	}

	// Read and apply mode status
	std::vector<int> mode_status;
	parse_tree(trap_cfg, "trap.mode_status", mode_status);

	for (int k = 0; k < mode_status.size(); k++) {
		if (mode_status[k] == 0) {
			for (int i = 0; i < eta_list[k].size();i++) {
				eta_list[k][i] = 0.0;
			}
		}
	}

	// Apply a scale to the eta_list
	double eta_scale = trap_cfg.get<double>("trap.eta_scale");
	for (int k = 0; k < eta_list.size(); k++) {
		for (int j = 0; j < eta_list[j].size(); j++) {
			eta_list[k][j] *= eta_scale;
		}
	}

    // Read the Hamiltonian
    h_cfg.stark_status = interaction_cfg.get<bool>("interaction.hamiltonian.stark_status");
    h_cfg.stark_scale_factor = interaction_cfg.get<double>("interaction.hamiltonian.stark_scale_factor");
    h_cfg.constant_light_shift = interaction_cfg.get<bool>("interaction.hamiltonian.constant_light_shift");

    h_cfg.mode = interaction_cfg.get<int>("interaction.hamiltonian.mode");

    h_cfg.expansion_order = interaction_cfg.get<int>("interaction.hamiltonian.MS.expansion_order");

    h_cfg.kerr_status = interaction_cfg.get<bool>("interaction.hamiltonian.MS.kerr");
    h_cfg.carrier_status = interaction_cfg.get<bool>("interaction.hamiltonian.MS.carrier");
    h_cfg.crosstalk_status = interaction_cfg.get<bool>("interaction.hamiltonian.MS.crosstalk");
    h_cfg.monitoring_status = interaction_cfg.get<bool>("interaction.hamiltonian.monitoring_status");
    h_cfg.spectator = interaction_cfg.get<int>("interaction.hamiltonian.spectator_ion");

    h_cfg.crosstalk_scale_factor = interaction_cfg.get<double>("interaction.hamiltonian.MS.crosstalk_scale_factor");

    h_cfg.n_ions = interaction_cfg.get<int>("interaction.hamiltonian.n_ions");

    std::vector<int> _gate_ions;
    parse_tree(interaction_cfg, "interaction.hamiltonian.gate_ions", _gate_ions);
    h_cfg.set_gate_ions(_gate_ions);

    // Read the Lindblads
    l_cfg.mode = interaction_cfg.get<int>("interaction.lindblads.mode");

    l_cfg.motional_heating=interaction_cfg.get<bool>("interaction.lindblads.motional_heating");
    l_cfg.motional_dephasing=interaction_cfg.get<bool>("interaction.lindblads.motional_dephasing");
    l_cfg.spin_dephasing_rayleigh=interaction_cfg.get<bool>("interaction.lindblads.spin_dephasing_rayleigh");
    l_cfg.spin_flip_raman=interaction_cfg.get<bool>("interaction.lindblads.spin_flip_raman");
    l_cfg.spontaneous_decay_01=interaction_cfg.get<bool>("interaction.lindblads.spontaneous_decay_01");
    l_cfg.spontaneous_decay_12=interaction_cfg.get<bool>("interaction.lindblads.spontaneous_decay_12");
    l_cfg.laser_dephasing_S=interaction_cfg.get<bool>("interaction.lindblads.laser_dephasing_S");
    l_cfg.laser_power_fluc_2Atom=interaction_cfg.get<bool>("interaction.lindblads.laser_power_fluc_2Atom");
    l_cfg.pink_noise=interaction_cfg.get<bool>("interaction.lindblads.pink_noise");

    l_cfg.laser_dephasing=interaction_cfg.get<bool>("interaction.lindblads.laser_dephasing");
    l_cfg.laser_power_fluctuation=interaction_cfg.get<bool>("interaction.lindblads.laser_power_fluctuation");
    l_cfg.position_operator=interaction_cfg.get<bool>("interaction.lindblads.position_operator");

    parse_tree(interaction_cfg, "interaction.lindblads.gamma_heating", l_cfg.gamma_heating);
    parse_tree(interaction_cfg, "interaction.lindblads.motional_coherence_time", l_cfg.motional_coherence_time);

    l_cfg.laser_coherence_time = interaction_cfg.get<double>("interaction.lindblads.laser_coherence_time");

    l_cfg.gamma_Rayleigh = interaction_cfg.get<double>("interaction.lindblads.gamma_Rayleigh");
    l_cfg.gamma_Raman = interaction_cfg.get<double>("interaction.lindblads.gamma_Raman");
    l_cfg.gamma_01 = interaction_cfg.get<double>("interaction.lindblads.gamma_01");
    l_cfg.gamma_12 = interaction_cfg.get<double>("interaction.lindblads.gamma_12");
    l_cfg.pink_cumm_1 = interaction_cfg.get<double>("interaction.lindblads.pink_cumm_1");

    l_cfg.amp_fluct = interaction_cfg.get<double>("interaction.lindblads.amp_fluct");

    l_cfg.position_meas_coupling = interaction_cfg.get<double>("interaction.lindblads.position_meas_coupling");

    // Read the state
    state_cfg.cutoff_mode = interaction_cfg.get<int>("interaction.state.motion.cutoff.mode");

    parse_tree(interaction_cfg, "interaction.state.motion.cutoff.auto_cutoffs", state_cfg.auto_cutoffs);
    parse_tree(interaction_cfg, "interaction.state.motion.cutoff.phonon_cutoffs", state_cfg.phonon_cutoffs);

    state_cfg.phonon_moving_basis = interaction_cfg.get<bool>("interaction.state.motion.moving_basis.status");
    state_cfg.shift_accuracy = interaction_cfg.get<double>("interaction.state.motion.moving_basis.shift_accuracy");
    state_cfg.use_x_basis = interaction_cfg.get<bool>("interaction.state.motion.moving_basis.use_x_basis");

    state_cfg.dynamic_cutoff = interaction_cfg.get<bool>("interaction.state.motion.dynamic_cutoff.status");
    state_cfg.cutoff_epsilon = interaction_cfg.get<double>("interaction.state.motion.dynamic_cutoff.cutoff_epsilon");
    state_cfg.cutoff_pad_size = interaction_cfg.get<int>("interaction.state.motion.dynamic_cutoff.cutoff_pad_size");

    state_cfg.spin_state_type = interaction_cfg.get<std::string>("interaction.state.spin.type");
    parse_tree(interaction_cfg, "interaction.state.spin.PURE", state_cfg.spin);

    state_cfg.motion_state_type = interaction_cfg.get<std::string>("interaction.state.motion.type");
    parse_tree(interaction_cfg, "interaction.state.motion.PURE", state_cfg.motion);
	parse_tree(interaction_cfg, "interaction.state.motion.THERMAL.n_bar", state_cfg.n_bar);
    state_cfg.cutoff_probability = interaction_cfg.get<double>("interaction.state.motion.THERMAL.cutoff_probability");


    // Read the pulse config
	// AM
	pulse_cfg.profile_am = interaction_cfg.get<int>("interaction.pulse.AM.profile");
	pulse_cfg.addressing_am = interaction_cfg.get<int>("interaction.pulse.AM.addressing");
	pulse_cfg.steps_am = interaction_cfg.get<int>("interaction.pulse.AM.steps");
	pulse_cfg.SYM_AM = interaction_cfg.get<bool>("interaction.pulse.AM.SYM_AM");

	// FM
	pulse_cfg.profile_fm = interaction_cfg.get<int>("interaction.pulse.FM.profile");
	pulse_cfg.addressing_fm = interaction_cfg.get<int>("interaction.pulse.FM.addressing");
	pulse_cfg.steps_fm = interaction_cfg.get<int>("interaction.pulse.FM.steps");
	pulse_cfg.SYM_FM = interaction_cfg.get<bool>("interaction.pulse.FM.SYM_FM");

	// Read the GData
	if (pulse_cfg.profile_am == 4) {
		load_2d_list(solutions_cfg.get<std::string>("solutions.g.g_list_path"), g_data);
	}

	// Read the sampled \Omega(t)
	if (pulse_cfg.profile_am == 5) {
		load_1d_list(solutions_cfg.get<std::string>("solutions.sampled.omega_path"), pulse_sampled);
	}

	// Read the rabi_freq_ratio
	parse_tree(interaction_cfg, "interaction.pulse.rabi_freq_ratio", rabi_freq_ratio);

    // Initialize the operators
    initialize_operators();

    // Read the outlist
	parse_tree(experiment_cfg, "experiment.outlist", outlist_string);
    set_outlist(outlist_string);


    // Read the flist
    parse_tree(experiment_cfg, "experiment.flist", flist);


    // Read the time
    double _dt = interaction_cfg.get<double>("interaction.time.dt");
    int _numdts = interaction_cfg.get<int>("interaction.time.numdts");
    int _numsteps = interaction_cfg.get<int>("interaction.time.numsteps");

    set_time(_dt, _numdts, _numsteps);


    // Read the processor
    // // processor_type = experiment_cfg.get<int>("experiment.processor");
    // // phi_res = experiment_cfg.get<double>("experiment.PARITY.phi_res");

    // Initialize the Hamiltonian and the Lindblads
    initialize_H();

    initialize_L();
}

qsd::Operator ion_trap::outlist_processor(const std::string &_operator) {
    /*
    * Operators
    * */

	std::vector<std::string> operator_string;

	std::stringstream string_stream(_operator);
	while(string_stream.good())
	{
		std::string substr;
		std::getline(string_stream, substr, '_');
		if (substr != "")
			operator_string.push_back(substr);
	}

    if (operator_string[0] == "sig"){
        if (operator_string[1] == "00")
            return S_00[0];
        else if (operator_string[1] == "11")
            return S_11[0];
        else if (operator_string[1] == "22")
            return S_22[0];
    }

    if (operator_string[0] == "S0") {
        if (operator_string[1] == "0")
            return S0[0];
        else if (operator_string[1] == "1")
            return S0[1];
    }

    if (operator_string[0] == "S1") {
        if (operator_string[1] == "0")
            return S1[0];
        else if (operator_string[1] == "1")
            return S1[1];
    }
    if (operator_string[0] == "D0") {
        if (operator_string[1] == "0")
            return D0[0];
        else if (operator_string[1] == "1")
            return D0[1];
    }
    if (operator_string[0] == "D1") {
        if (operator_string[1] == "0")
            return D1[0];
        else if (operator_string[1] == "1")
            return D1[1];
    }

    if (operator_string[0] == "sx")
	{
		if (operator_string[1] == "1")
			return sx[0];
		else if (operator_string[1] == "2")
			return sx[1];
	}

	else if (operator_string[0] == "sy")
	{
		if (operator_string[1] == "1")
			return sy[0];
		else if (operator_string[1] == "2")
			return sy[1];
	}

	else if (operator_string[0] == "sz")
	{std::cout<<'here'<<std::endl;
		if (operator_string[1] == "1")
			return sz[0];
		else if (operator_string[1] == "2")
			return sz[1];
	}

	else if (operator_string[0] == "sp")
	{
		if (operator_string[1] == "1")
			return sp[0];
		else if (operator_string[1] == "2")
			return sp[1];
	}

	else if (operator_string[0] == "a")
	{
		if (operator_string[1] == "1")
			return a[0];
		else if (operator_string[1] == "2")
			return a[1];
		else if (operator_string[1] == "3")
			return a[2];
		else if (operator_string[1] == "4")
			return a[3];
		else if (operator_string[1] == "5")
			return a[4];
		else if (operator_string[1] == "6")
			return a[5];
		else if (operator_string[1] == "7")
			return a[6];
        else if (operator_string[1] == "8")
            return a[7];
	}

	else if (operator_string[0] == "N")
	{
		if (operator_string[1] == "1")
			return N[0];
		else if (operator_string[1] == "2")
			return N[1];
		else if (operator_string[1] == "3")
			return N[2];
		else if (operator_string[1] == "4")
			return N[3];
		else if (operator_string[1] == "5")
			return N[4];
		else if (operator_string[1] == "6")
			return N[5];
		else if (operator_string[1] == "7")
			return N[6];
        else if (operator_string[1] == "8")
            return a[7];
	}

	else if (operator_string[0] == "id")
	{
		if (operator_string[1] == "s1")
			return id_s[0];
		else if (operator_string[1] == "s2")
			return id_s[1];

		else if (operator_string[1] == "m1")
			return id_m[0];
		else if (operator_string[1] == "m2")
			return id_m[1];
		else if (operator_string[1] == "m3")
			return id_m[2];
		else if (operator_string[1] == "m4")
			return id_m[3];
		else if (operator_string[1] == "m5")
			return id_m[4];
		else if (operator_string[1] == "m6")
			return id_m[5];
		else if (operator_string[1] == "m7")
			return id_m[6];

		else if (operator_string.size() == 1)
			return id;
	}

	else if (operator_string[0] == "sm")
	{
		if (operator_string[1] == "1")
			return sm[0];
		else if (operator_string[1] == "2")
			return sm[1];
	}

	else if (operator_string[0] == "p00")
	{
		if (operator_string[1] == "1")
			return projectors_z[0][0];
		else if (operator_string[1] == "2")
			return projectors_z[1][0];
	}

	else if (operator_string[0] == "p11")
	{
		if (operator_string[1] == "1")
			return projectors_z[0][3];
		else if (operator_string[1] == "2")
			return projectors_z[1][3];
	}

	else if (operator_string[0] == "rho" && operator_string[1] == "ideal")
	{
		if (operator_string[2] == "1")
			return rho_ideals[0];
		else if (operator_string[2] == "2")
			return rho_ideals[1];
	}

	else if (operator_string[0] == "ad")
	{
		if (operator_string[1] == "1")
			return ad[0];
		else if (operator_string[1] == "2")
			return ad[1];
		else if (operator_string[1] == "3")
			return ad[2];
		else if (operator_string[1] == "4")
			return ad[3];
		else if (operator_string[1] == "5")
			return ad[4];
		else if (operator_string[1] == "6")
			return ad[5];
		else if (operator_string[1] == "7")
			return ad[6];
	}

	else if (operator_string[0] == "x")
	{
		if (operator_string[1] == "1")
			return x[0];
		else if (operator_string[1] == "2")
			return x[1];
		else if (operator_string[1] == "3")
			return x[2];
		else if (operator_string[1] == "4")
			return x[3];
		else if (operator_string[1] == "5")
			return x[4];
		else if (operator_string[1] == "6")
			return x[5];
		else if (operator_string[1] == "7")
			return x[6];
	}

	else if (operator_string[0] == "p")
	{
		if (operator_string[1] == "1")
			return p[0];
		else if (operator_string[1] == "2")
			return p[1];
		else if (operator_string[1] == "3")
			return p[2];
		else if (operator_string[1] == "4")
			return p[3];
		else if (operator_string[1] == "5")
			return p[4];
		else if (operator_string[1] == "6")
			return p[5];
		else if (operator_string[1] == "7")
			return p[6];
	}

	else if (operator_string[0] == "joint")
	{
		// joint operator between spins
		// joint_(spin_basis)_(spin_1_projector)_(spin_2_projector)
		if (operator_string.size() == 4) {

			std::string spin_basis = operator_string[1];

			int spin_1_projector = std::stoi(operator_string[2]);
			int spin_2_projector = std::stoi(operator_string[3]);

			if (spin_basis == "x")
				return projectors_x[0][spin_1_projector] * projectors_x[1][spin_2_projector];
			else if (spin_basis == "z")
				return projectors_z[0][spin_1_projector] * projectors_z[1][spin_2_projector];
		}

		// joint operator between spin and motion
		// joint_(spin_basis)_(spin_1_projector)_(spin_2_projector)_(x or p)_(mode_idx)
		if (operator_string.size() == 6) {

			std::string spin_basis = operator_string[1];

			int spin_1_projector = std::stoi(operator_string[2]);
			int spin_2_projector = std::stoi(operator_string[3]);


			std::string x_or_p = operator_string[4];
			int mode_idx = std::stoi(operator_string[5]);
            if (spin_1_projector==1)
                spin_1_projector = 3;
            if (spin_2_projector==1)
                spin_2_projector = 3;

			if (x_or_p == "x") {

				if (spin_basis == "x")
					return projectors_x[0][spin_1_projector] * projectors_x[1][spin_2_projector] * x[mode_idx];
				else if (spin_basis == "z")
					return projectors_z[0][spin_1_projector] * projectors_z[1][spin_2_projector] * x[mode_idx];

			} else if (x_or_p == "p") {

				if (spin_basis == "x")
					return projectors_x[0][spin_1_projector] * projectors_x[1][spin_2_projector] * p[mode_idx];
				else if (spin_basis == "z")
					return projectors_z[0][spin_1_projector] * projectors_z[1][spin_2_projector] * p[mode_idx];
			}
		}
	}

    else if (operator_string[0] == "xx")
        return projectors_x[0][0] * projectors_x[1][0];

    std::cout << "Error: Unrecognized operator in outlist.\n";
    return qsd::NullOperator();
}

void ion_trap::set_delta(std::vector<double> _delta) {

	std::vector<std::vector<double>>().swap(delta);

    delta.resize(h_cfg.n_gate_ions);
    for (int i = 0; i < h_cfg.n_gate_ions; i++) {
        delta[i] = _delta;
    }

    set_phonon_cutoffs();
}

void ion_trap::set_delta(std::vector<std::vector<double>> _delta) {

	std::vector<std::vector<double>>().swap(delta);

    delta.resize(h_cfg.n_gate_ions);
    for (int i = 0; i < h_cfg.n_gate_ions; i++) {
        delta[i] = _delta[i];
    }

    set_phonon_cutoffs();
}

void ion_trap::set_delta(double _delta) {

	std::vector<std::vector<double>>().swap(delta);

    delta.resize(h_cfg.n_gate_ions);
    for (int i = 0; i < h_cfg.n_gate_ions; i++) {
        delta[i].resize(pulse_cfg.steps_fm);
        for (int j = 0; j < pulse_cfg.steps_fm; j++) {
            delta[i][j] = _delta;
        }
    }
    set_phonon_cutoffs();
}

void ion_trap::set_delta(int _delta_idx) {
    set_delta(delta_list[_delta_idx]);
}


void ion_trap::set_external_delta(boost::property_tree::ptree external_delta_cfg, std::string select) {
	/*
	* Load the external delta
	* */

	int delta_mode = external_delta_cfg.get<int>("delta_mode");

	switch (delta_mode) {

		case 0: { // Absolute

			std::vector<double> external_delta;

			char delimiter = external_delta_cfg.get<char>("delimiter");

			std::stringstream string_stream(external_delta_cfg.get<std::string>(select + ".external_delta"));
			while(string_stream.good())
			{
				std::string substr;
				std::getline(string_stream, substr, delimiter);
				if (substr != "")
					external_delta.push_back(std::stod(substr));
			}

			/*for (auto& item : tree.get_child("experiment.simulation.mode-0.external_delta")) {
				external_delta = item.second.get_value<double>();
			}*/

			set_delta(external_delta);

			break;
		}

		case 1: { // Relative

			std::vector<int> reference_modes;
			std::vector<double> detuning_ratios;

			std::stringstream string_stream_modes(external_delta_cfg.get<std::string>(select + ".reference_modes"));
			while(string_stream_modes.good())
			{
				std::string substr;
				std::getline(string_stream_modes, substr, ':');
				if (substr != "")
					reference_modes.push_back(std::stoi(substr));
			}

			std::stringstream string_stream_ratios(external_delta_cfg.get<std::string>(select + ".detuning_ratios"));
			while(string_stream_ratios.good())
			{
				std::string substr;
				std::getline(string_stream_ratios, substr, ':');
				if (substr != "")
					detuning_ratios.push_back(std::stod(substr));
			}

			double detuning = 0;
			if (std::fabs(detuning_ratios[1] - 1.0) < 1e-6) {
				detuning = nu_list[reference_modes[1] - 1] + (nu_list[reference_modes[0] - 1] - nu_list[reference_modes[1] - 1]) / (1.0 - detuning_ratios[0]);
			} else {
				detuning = nu_list[reference_modes[0] - 1] + (nu_list[reference_modes[1] - 1] - nu_list[reference_modes[0] - 1]) / (1.0 - detuning_ratios[1]);
			}

			set_delta(detuning);

			break;
		}

		default: {
			break;
		}

	}
}

std::vector<double> ion_trap::get_delta() {

    switch (pulse_cfg.addressing_fm) {
        case 0: { // Coupled
            return delta[0];
        }

        case 1: { // Individual
            return {0.0};
        }

        default: {
            std::cout << "Error: Unknown addressing FM mode.\n";
            return {0.0};
        }
    }
}

double ion_trap::get_energy_drift(int idx, double t) {
    return energy_drift;
}


void ion_trap::set_pulse(std::vector<double> &_pulse) {

	std::vector<std::vector<double>>().swap(pulse);

    pulse.resize(h_cfg.n_gate_ions);
    for (int i = 0; i < h_cfg.n_gate_ions; i++) {
        pulse[i] = _pulse;
    }

}

void ion_trap::set_pulse(std::vector<std::vector<double>> &_pulse) {

	std::vector<std::vector<double>>().swap(pulse);

    pulse.resize(h_cfg.n_gate_ions);
    for (int i = 0; i < h_cfg.n_gate_ions; i++) {
        pulse[i] = _pulse[i];
    }
}

void ion_trap::set_pulse(int _pulse_idx){
    set_pulse(pulse_list[_pulse_idx]);
}

void ion_trap::set_EIT_pulse(std::vector<double> _pulse){
    EIT_pulse.resize(2);
    for (int i=0;i<2;i++) {
        EIT_pulse[i]=_pulse;
    }
}

void ion_trap::set_shelving_pulse(std::vector<double> &_pulse) {
    std::vector<std::vector<double>>().swap(shelving_pulse);
    shelving_pulse.resize(h_cfg.n_gate_ions);
    for (int i = 0; i < h_cfg.n_gate_ions; i++){
        shelving_pulse[i] = _pulse;
    }
}

void ion_trap::set_shelving_pulse(int _pulse_idx) {
    set_shelving_pulse(shelving_pulse_list[_pulse_idx]);
}

void ion_trap::set_shelving_phase(std::vector<double> &_phas) {
    std::vector<std::vector<double>>().swap(shaving_phase);
    shaving_phase.resize(h_cfg.n_gate_ions);
    for (int i = 0; i < h_cfg.n_gate_ions; i++) {
        shaving_phase[i].resize(_phas.size());
        shaving_phase[i] = _phas;
    }
}

void ion_trap::set_shelving_phase(int _phase_idx) {
    set_shelving_phase(shaving_phase_list[_phase_idx]);
}

void ion_trap::set_Delta_pump(double det){
    Delta_pump.resize(1);
    Delta_pump[0] = det;
}

void ion_trap::set_pulse_external(boost::property_tree::ptree external_pulse_cfg, std::string select) {
	/*
	* Load the external pulse
	* */
	std::vector<double> external_pulse;

	char delimiter = external_pulse_cfg.get<char>("delimiter");
	double pulse_scale = external_pulse_cfg.get<double>(select + ".pulse_scale");

	std::stringstream string_stream(external_pulse_cfg.get<std::string>(select + ".external_pulse"));
	while(string_stream.good())
	{
		std::string substr;
		std::getline(string_stream, substr, delimiter);
		if (substr != "")
			external_pulse.push_back(std::stod(substr));
	}

	/*for (auto& item : tree.get_child("experiment.simulation.mode-0.external_pulse")) {
		external_pulse.emplace_back(item.second.get_value<double>());
	}*/

	set_pulse(external_pulse);
	pulse_cfg.scale_am = pulse_scale;
}

void ion_trap::set_time(double _dt, int _numdts, int _numsteps) {
    dt = _dt;
    numdts = _numdts;
    numsteps = _numsteps;

    t_gate = _dt * (double) (_numdts * _numsteps);
}

void ion_trap::set_time(double t) {
    dt = 0.1;
    numdts = 1;
    numsteps = (int) (t / 0.1);

    t_gate = 0.1 * (double) (1 * (int) (t / 0.1));
}

void ion_trap::set_phonon_cutoffs()
{
//    std::cout<< simulation_mode << std::endl;

    if (state_cfg.cutoff_mode == 1) { // Autonomous assignment
        std::vector<std::pair<int, double>> modes;
        for (int i = 0; i < h_cfg.n_ions; i++) {
            modes.emplace_back(std::make_pair(i, nu_list[i]));
        }

        std::sort(modes.begin(), modes.end(), [&](const std::pair<int, double> &left, const std::pair<int, double> &right) {
            return std::abs(left.second - delta[0][0]) < std::abs(right.second - delta[0][0]);
        });

        state_cfg.phonon_cutoffs.resize(h_cfg.n_ions);
        for (int i = 0; i < h_cfg.n_ions; i++) {
            if (i <= 1) {
                state_cfg.phonon_cutoffs[modes[i].first] = state_cfg.auto_cutoffs[0];
            } else if (i >= 2 && i <= 3) {
                state_cfg.phonon_cutoffs[modes[i].first] = state_cfg.auto_cutoffs[1];
            } else {
                state_cfg.phonon_cutoffs[modes[i].first] = state_cfg.auto_cutoffs[2];
            }
        }
    }

    set_state();
}

double ion_trap::get_stark_shift(int idx, double t) {
    /*
    * The input pulse is in units of Mrad/s
    * */
    if (h_cfg.stark_status)
        if (h_cfg.constant_light_shift)
            return h_cfg.stark_scale_factor * 1;
        else
            return h_cfg.stark_scale_factor * get_omega(idx, t);

    else
        return 0.0;
}

bool ion_trap::findParity(int x)
{
    int y = x ^ (x >> 1);
    y = y ^ (y >> 2);
    y = y ^ (y >> 4);
    y = y ^ (y >> 8);
    y = y ^ (y >> 16);

    // Rightmost bit of y holds the parity value
    // if (y&1) is 1 then parity is odd else even
    if (y & 1)
        return 1;
    return 0;
}

void ion_trap::initialize_operators() {
    double _pi = 3.141592653589793238;

    for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
        sx.emplace_back(qsd::SigmaX{i});
        sy.emplace_back(qsd::SigmaY{i});
        sz.emplace_back(qsd::SigmaZ{i});
        sp.emplace_back(qsd::SigmaPlus{i});
        id_s.emplace_back(qsd::IdentityOperator{i});
    }

    for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
        sx[i].resetStackPointer();
        sy[i].resetStackPointer();
        sz[i].resetStackPointer();
        sp[i].resetStackPointer();
        id_s[i].resetStackPointer();
    }

    for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
        sm.emplace_back(sp[i].hc());
    }

    /*
    qsd::Operator pp_1 = 0.5 * (p00_1 + p11_1 + p01_1 + p10_1);
    qsd::Operator pp_2 = 0.5 * (p00_2 + p11_2 + p01_2 + p10_2);

    qsd::Operator mm_1 = 0.5 * (p00_1 + p11_1 - p01_1 - p10_1);
    qsd::Operator mm_2 = 0.5 * (p00_2 + p11_2 - p01_2 - p10_2);
    */

    /*
     * Projectors
     *
     * Notice:
     * 1) In "projectors_z[i][j]", "i" is the spin index and "j" is the projector index.
     *         j = 0: |0><0|
     *         j = 1: |0><1|
     *         j = 2: |1><0|
     *         j = 3: |1><1|
     *
     * 2) In "projectors_x[i][j]", "i" is the spin index and "j" is the projector index.
     *         j = 0: |+><+|
     *         j = 1: |+><-|
     *         j = 2: |-><+|
     *         j = 3: |-><-|
     */

    std::vector<std::vector<qsd::Operator>>().swap(projectors_z);
    for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
        std::vector<qsd::Operator> proj_z;

        proj_z.push_back(sm[i] * sp[i]); // |0><0|
        proj_z.push_back(sm[i]); // |0><1|
        proj_z.push_back(sp[i]); // |1><0|
        proj_z.push_back(sp[i] * sm[i]); // |1><1|
        projectors_z.push_back(proj_z);
    }

     for (int i = 0; i < h_cfg.n_interacting_ions; i++) {

        S_01.emplace_back(qsd::TransitionOperator{1,0,i}); // |1><0|
        S_12.emplace_back(qsd::TransitionOperator{2,1,i}); // |2><1|
        S_02.emplace_back(qsd::TransitionOperator{2,0,i}); // |2><1|
        S_00.emplace_back(qsd::TransitionOperator{0,0,i}); // |0><0|
        S_11.emplace_back(qsd::TransitionOperator{1,1,i}); // |1><1|
        S_22.emplace_back(qsd::TransitionOperator{2,2,i}); // |2><2|

        S0_D0.emplace_back(qsd::TransitionOperator{2,0,i}); //
        S0_D1.emplace_back(qsd::TransitionOperator{3,0,i}); //
        S1_D0.emplace_back(qsd::TransitionOperator{2,1,i}); //
        S1_D1.emplace_back(qsd::TransitionOperator{3,1,i}); //
        S0_S1.emplace_back(qsd::TransitionOperator{1,0,i}); //
        S1_S0.emplace_back(qsd::TransitionOperator{0,1,i}); //
        D0_S0.emplace_back(qsd::TransitionOperator{0,2,i}); //
        D1_S1.emplace_back(qsd::TransitionOperator{1,3,i}); //

        S0.emplace_back(qsd::TransitionOperator{0,0,i}); //
        S1.emplace_back(qsd::TransitionOperator{1,1,i}); //
        D0.emplace_back(qsd::TransitionOperator{2,2,i}); //
        D1.emplace_back(qsd::TransitionOperator{3,3,i}); //

     }

    for (int i = 0; i < h_cfg.n_interacting_ions; i++){
         S_01[i].resetStackPointer();
         S_12[i].resetStackPointer();
         S_02[i].resetStackPointer();
         S_00[i].resetStackPointer();
         S_11[i].resetStackPointer();
         S_22[i].resetStackPointer();
         S0_D0[i].resetStackPointer();
         D0_S0[i].resetStackPointer();
         S0_D1[i].resetStackPointer();
         S1_D0[i].resetStackPointer();
         S1_D1[i].resetStackPointer();
         D1_S1[i].resetStackPointer();
         S0_S1[i].resetStackPointer();
         S0[i].resetStackPointer();
         S1[i].resetStackPointer();
         D0[i].resetStackPointer();
         D1[i].resetStackPointer();
    }


    proj_z0_basis.emplace_back(S0[0]);
    proj_z0_basis.emplace_back(S0_D0[0]);
    proj_z0_basis.emplace_back(D0_S0[0]);
    proj_z0_basis.emplace_back(D0[0]);

    proj_z1_basis.emplace_back(S1[0]);
    proj_z1_basis.emplace_back(S1_D1[0]);
    proj_z1_basis.emplace_back(D1_S1[0]);
    proj_z1_basis.emplace_back(D1[0]);

    proj_z3_basis.emplace_back(S1[0]);
    proj_z3_basis.emplace_back(S0_S1[0]);
    proj_z3_basis.emplace_back(S1_S0[0]);
    proj_z3_basis.emplace_back(D1[0]);

    for (int i = 0; i < 4; ++i) {
            proj_z0_basis[i].resetStackPointer();
            proj_z1_basis[i].resetStackPointer();
    }


    std::stringstream string_stream(pulse_sequence);
//    std::cout<<pulse_sequence<<std::endl;
    int len = 0;
    std::string substr;

    while(string_stream.good())
    {
        qsd::TransitionOperator Aux0;
        qsd::TransitionOperator Aux1;
        qsd::TransitionOperator Aux2;
        qsd::TransitionOperator Aux3;

        std::getline(string_stream, substr, '_');
        if (substr == "x") {
//            std::cout<<"x"<<std::endl;

            Aux0 = S0_D0[0];
            Aux1 = D0_S0[0];
            Aux1 += Aux0;

            pulse_seq_0.emplace_back(Aux1); //

            Aux0 = S1_D1[0];
            Aux1 = D1_S1[0];
            Aux1 += Aux0;
            pulse_seq_1.emplace_back(Aux1); //
        } else if (substr == "mx") {
//            std::cout<<"mx"<<std::endl;

            Aux0 = S0_D0[0];
            Aux1 = D0_S0[0];
            Aux1 += Aux0;
            Aux1 *= -1;

            pulse_seq_0.emplace_back(Aux1); //

            Aux0 = S1_D1[0];
            Aux1 = D1_S1[0];
            Aux1 += Aux0;
            Aux1 *= -1;
            pulse_seq_1.emplace_back(Aux1); //
        } else if (substr == "y") {
//            std::cout<<"y"<<std::endl;

            Aux0 = S0_D0[0];
            Aux0 *= Im;
            Aux1 = D0_S0[0];
            Aux1 *= mIm;
            Aux1 += Aux0;

            pulse_seq_0.emplace_back(Aux1); //

            Aux0 = S1_D1[0];
            Aux0 *= Im;
            Aux1 = D1_S1[0];
            Aux1 *= mIm;
            Aux1 += Aux0;
            pulse_seq_1.emplace_back(Aux1); //

         } else if (substr == "my") {
//            std::cout<<"y"<<std::endl;

            Aux0 = S0_D0[0];
            Aux0 *= Im;
            Aux1 = D0_S0[0];
            Aux1 *= mIm;
            Aux1 += Aux0;
            Aux1 *= -1;

            pulse_seq_0.emplace_back(Aux1); //

            Aux0 = S1_D1[0];
            Aux0 *= Im;
            Aux1 = D1_S1[0];
            Aux1 *= mIm;
            Aux1 += Aux0;
            Aux1 *= -1;

            pulse_seq_1.emplace_back(Aux1); //

        } else if (substr == "z") {
//            std::cout<<"ZZ"<<std::endl;
            Aux0 = S0[0];
            Aux0 -= D0[0];
            Aux1 = S1[0];
            Aux1 -= D1[0];
            pulse_seq_0.emplace_back(Aux0);
            pulse_seq_1.emplace_back(Aux1);
        } else if (substr == "60") {
            double phi = _pi/3;
            Aux0 = S0_D0[0];
            Aux1 = D0_S0[0];
            Aux1 += Aux0; // \sigma^x
            Aux1 *= cos(phi);

            Aux2 = S0_D0[0];
            Aux2 *= Im;
            Aux3 = D0_S0[0];
            Aux3 *= mIm;
            Aux3 += Aux2;
            Aux3 *= sin(phi);
//            std::cout<<sin(phi)<<std::endl;

            Aux3 += Aux1;
            pulse_seq_0.emplace_back(Aux3);

            Aux0 = S1_D1[0];
            Aux1 = D1_S1[0];
            Aux1 += Aux0; // \sigma^x
            Aux1 *= cos(phi);

            Aux2 = S1_D1[0];
            Aux2 *= Im;
            Aux3 = D1_S1[0];
            Aux3 *= mIm;
            Aux3 += Aux2;
            Aux3 *= sin(phi);

            Aux3 += Aux1;

            pulse_seq_1.emplace_back(Aux3);
        } else if (substr == "30") {
//            std::cout<<"30"<<std::endl;
            double phi = _pi/6;
            Aux0 = S0_D0[0];
            Aux1 = D0_S0[0];
            Aux1 += Aux0; // \sigma^x
            Aux1 *= cos(phi);

            Aux2 = S0_D0[0];
            Aux2 *= Im;
            Aux3 = D0_S0[0];
            Aux3 *= mIm;
            Aux3 += Aux2;
            Aux3 *= sin(phi);

            Aux3 += Aux1;
            pulse_seq_0.emplace_back(Aux3);

            Aux0 = S1_D1[0];
            Aux1 = D1_S1[0];
            Aux1 += Aux0; // \sigma^x
            Aux1 *= cos(phi);

            Aux2 = S1_D1[0];
            Aux2 *= Im;
            Aux3 = D1_S1[0];
            Aux3 *= mIm;
            Aux3 += Aux2;
            Aux3 *= sin(phi);

            Aux3 += Aux1;

            pulse_seq_1.emplace_back(Aux3);
        } else if (substr == "p"){
            Aux0 = S0_S1[0];
            Aux0 += S0_S1[0].hc();
            pulse_seq_0.emplace_back(Aux0);

            Aux1 = S0_S1[1];
            Aux1 += S0_S1[1].hc();
            pulse_seq_1.emplace_back(Aux1);
        }
        len += 1;
    }
    num_pulses = len;


    qsd::TransitionOperator Aux0;
    qsd::TransitionOperator Aux1;
    Aux0 = S0_D0[0];
    Aux0 += D0_S0[0];
    sigma_X.emplace_back(Aux0);
    Aux0 = S1_D1[0];
    Aux0 += D1_S1[0];
    sigma_X.emplace_back(Aux0);

    Aux0 = S0_D0[0];
    Aux0 *= Im;
    Aux1 = D0_S0[0];
    Aux1 *= mIm;
    Aux1 += Aux0;
    sigma_Y.emplace_back(Aux1);
    Aux0 = S1_D1[0];
    Aux0 *= Im;
    Aux1 = D1_S1[0];
    Aux1 *= mIm;
    Aux1 += Aux0;
    sigma_Y.emplace_back(Aux1);

    std::vector<std::vector<qsd::Operator>>().swap(projectors_x);
    for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
        std::vector<qsd::Operator> proj_x;

        proj_x.push_back(0.5 * (projectors_z[i][0] + projectors_z[i][1] + projectors_z[i][2] + projectors_z[i][3])); // |+><+|
        proj_x.push_back(0.5 * (projectors_z[i][0] - projectors_z[i][1] + projectors_z[i][2] - projectors_z[i][3])); // |+><-|
        proj_x.push_back(0.5 * (projectors_z[i][0] + projectors_z[i][1] - projectors_z[i][2] - projectors_z[i][3])); // |-><+|
        proj_x.push_back(0.5 * (projectors_z[i][0] - projectors_z[i][1] - projectors_z[i][2] + projectors_z[i][3])); // |-><-|

        projectors_x.push_back(proj_x);
//        S_01[i].resetStackPointer();
    }

    /*
     * \rho ideal
     *
     * rho_ideal_1 = 0.5 * (p00_1 * p00_2 + p11_1 * p11_2 + IM * (-p01_1 * p01_2 + p10_1 * p10_2))
     * rho_ideal_2 = 0.5 * (p00_1 * p00_2 + p11_1 * p11_2 + IM * (p01_1 * p01_2 - p10_1 * p10_2)) // \psi_{ideal} of Manning
     */
    std::vector<qsd::Operator>().swap(rho_ideals);

	if (h_cfg.n_gate_ions > 1) {
		rho_ideals.push_back(0.5 * (projectors_z[h_cfg.gate_ions_idx[0]][0] * projectors_z[h_cfg.gate_ions_idx[1]][0] + projectors_z[h_cfg.gate_ions_idx[0]][3] * projectors_z[h_cfg.gate_ions_idx[1]][3]
                                    + IM * (-projectors_z[h_cfg.gate_ions_idx[0]][1] * projectors_z[h_cfg.gate_ions_idx[1]][1] + projectors_z[h_cfg.gate_ions_idx[0]][2] * projectors_z[h_cfg.gate_ions_idx[1]][2])));

		rho_ideals.push_back(0.5 * (projectors_z[h_cfg.gate_ions_idx[0]][0] * projectors_z[h_cfg.gate_ions_idx[1]][0] + projectors_z[h_cfg.gate_ions_idx[0]][3] * projectors_z[h_cfg.gate_ions_idx[1]][3]
                                    + IM * (projectors_z[h_cfg.gate_ions_idx[0]][1] * projectors_z[h_cfg.gate_ions_idx[1]][1] - projectors_z[h_cfg.gate_ions_idx[0]][2] * projectors_z[h_cfg.gate_ions_idx[1]][2]))); // \psi_{ideal} of Manning
	} else {
		rho_ideals.push_back(id_s[0]);
		rho_ideals.push_back(id_s[0]);
	}

    for (int i = 0; i < h_cfg.n_ions; i++) {
        x.emplace_back(qsd::XOperator(h_cfg.n_interacting_ions + i));
        p.emplace_back(qsd::POperator(h_cfg.n_interacting_ions + i));

        a.emplace_back(qsd::AnnihilationOperator{h_cfg.n_interacting_ions + i});

        N.emplace_back(qsd::NumberOperator{h_cfg.n_interacting_ions + i});

        id_m.emplace_back(qsd::IdentityOperator{h_cfg.n_interacting_ions + i});
    }

    for (int i = 0; i < h_cfg.n_ions; i++) {
        x[i].resetStackPointer();
        p[i].resetStackPointer();

        a[i].resetStackPointer();

        N[i].resetStackPointer();

        id_m[i].resetStackPointer();
    }

    for (int i = 0; i < h_cfg.n_ions; i++) {
        ad.emplace_back(a[i].hc());
    }

    // The identity operator
    id = id_s[0];
    for (int i = 1; i < h_cfg.n_interacting_ions; i++)
        id = id * id_s[i];

    for (int i = 0; i < h_cfg.n_ions; i++)
        id = id * id_m[i];
}

double ion_trap::get_g(int idx, double t) {

	if (pulse_cfg.profile_am != 4) { // Not the custom g(t) mode

		return get_omega(idx, t) * h_cos(idx, t);

	} else { // Custom g(t) function

		double _pulse = 0.0;
		for (int i = 0; i < h_cfg.n_gate_ions; i++) {
			for (int j = 0; j < g_data[0].size(); j++) { // n_list = g_data[0], a_list = g_data[1]

				double _amp = h_cfg.interacting_ions[idx].crosstalk_scale_factors[i] * pulse_cfg.scale_am * g_data[1][j];
				_amp *= rabi_freq_ratio[h_cfg.interacting_ions[idx].idx];

				if (h_cfg.stark_status)
					_pulse += _amp * sin((2.0 * PI * g_data[0][j] / t_gate + h_cfg.stark_scale_factor * _amp) * t);
				else
					_pulse += _amp * sin((2.0 * PI * g_data[0][j] / t_gate) * t);
			}
		}
		return _pulse;
	}
}

double ion_trap::get_g_fluc(int idx, double t) {

	if (pulse_cfg.profile_am != 4) { // Not the custom g(t) mode

		return h_cos(idx, t);

	} else { // Custom g(t) function
		// Todo: implement the effect of Stark Shift

		double _pulse = 0.0;

		for (int i = 0; i < g_data[0].size(); i++) { // n_list = g_data[0], a_list = g_data[1]
			_pulse += sin((2.0 * PI * g_data[0][i] / t_gate) * t);
		}

		return _pulse;
	}
}

double ion_trap::get_EIT_pulse(int idx){
    return EIT_pulse[0][idx];
}

double ion_trap::get_shelving_pulse(int idx, double t, int _flag=0){
    int steps_AM = (int) shelving_pulse[0].size();

    if (_flag==0) {

        double segments = std::accumulate(shelving_time.begin(), shelving_time.end(), 0.0);
        std::vector<double> track(num_pulses);
//        std::cout<<steps_AM<<std::endl;


        for (int i = 0; i <= num_pulses; ++i) {
            track[i] = 0.0;
            for (int j = 0; j <= i; j++) {
                track[i] += shelving_time[j];
            }
        }

        int idx_AM;

        int j = 0;
        if (idx == 0) {
            if (t <= t_gate * track[0] / segments) {
//            std::cout<<idx<<std::endl;
                return shelving_pulse[0][idx] * pulse_cfg.scale_am;
            } else {
                return 0.0;
            }
        } else {
            if (t > t_gate * track[idx - 1] / segments & t < t_gate * track[idx] / segments) {
//            if (idx==5){
//                std::cout<<idx<<std::endl;
//            }
                return shelving_pulse[0][idx] * pulse_cfg.scale_am;

            } else {
                return 0.0;
            }
        }

        idx_AM = (int) (t * (shelving_time[idx] / t_gate));

        if (idx_AM == steps_AM) {
            idx_AM -= 1;
        }
    } else {

//        int steps_AM = pulse_cfg.steps_am * (pulse_cfg.SYM_AM ? 2.0 : 1.0);

        int idx_AM = 0;
        double dt_AM = 0.0;

        double tt = 0.0;
        if (pulse_cfg.profile_am == 0) { // piecewise constant
            idx_AM = (int) (t * (steps_AM / t_gate));

            if (idx_AM == steps_AM) {
                idx_AM -= 1;
            }

            return shelving_pulse[0][idx_AM]* pulse_cfg.scale_am;
        } else {
            return 0.0;
        }
    }
}

double ion_trap::get_Delta_pump(int idx) {
    return Delta_pump[0];
}

double ion_trap::get_coloured_noise(double t, int idx) {
    return Delta_pump[0];
}

double ion_trap::get_omega(int idx, double t) {

    int steps_AM = pulse_cfg.steps_am * (pulse_cfg.SYM_AM ? 2.0 : 1.0);

    int idx_AM = 0;
    double dt_AM = 0.0;

    double tt = 0.0;
    if (pulse_cfg.profile_am == 0) { // piecewise constant
        idx_AM = (int) (t * (steps_AM / t_gate));

		if (idx_AM == steps_AM) {
			idx_AM -= 1;
		}

		if (pulse_cfg.SYM_AM) {
			if (idx_AM >= pulse_cfg.steps_am) {
				idx_AM = pulse_cfg.steps_am - idx_AM % pulse_cfg.steps_am - 1;
			}
		}

		double _pulse = 0.0;
		for (int i = 0; i < h_cfg.n_gate_ions; i++) {
			try {
				_pulse += h_cfg.interacting_ions[idx].crosstalk_scale_factors[i] * pulse[i].at(idx_AM);
			}
			catch (std::out_of_range) {
				// // std::cout << "std::out_of_range" << std::endl;
				_pulse += h_cfg.interacting_ions[idx].crosstalk_scale_factors[i] * pulse[i].at(idx_AM - 1);
			}
		}
		return rabi_freq_ratio[h_cfg.interacting_ions[idx].idx] * pulse_cfg.scale_am * _pulse;

    } else if (pulse_cfg.profile_am == 1) { // error function
		// Todo: scale and rabi ratio are not implemented!

        dt_AM = t_gate / pulse_cfg.steps_am;
        tt = t - dt_AM / 2;

        idx_AM = (int) (tt * (steps_AM / t_gate));

		if (idx_AM == steps_AM) {
			idx_AM -= 1;
		}

		if (pulse_cfg.SYM_AM) {
			if (idx_AM >= pulse_cfg.steps_am) {
				idx_AM = pulse_cfg.steps_am - idx_AM % pulse_cfg.steps_am - 1;
			}
		}

		double t_idx_AM_begin = idx_AM * dt_AM;
		double t_idx_AM_end = (idx_AM + 1.0) * dt_AM;

		if (!pulse_cfg.SYM_AM) {
			if (idx_AM == steps_AM - 1)
				return pulse[idx][idx_AM];
			else
				return (pulse[idx][idx_AM] + pulse[idx][idx_AM + 1]) / 2.0 + (pulse[idx][idx_AM + 1] - pulse[idx][idx_AM]) / 2 * erf(5.0 / dt_AM * (tt - (t_idx_AM_begin + t_idx_AM_end) / 2.0));
		} else {
			return pulse[idx][idx_AM];
		}

    } else if (pulse_cfg.profile_am == 2) { // Fourier series

		double _pulse = pulse[idx][0];

		int n = 1;
		for (int i = 1; i < pulse_cfg.steps_am; i++) {

			if (i % 2 == 0) { // even index (sin term)
				_pulse += pulse[idx][i] * sin((2.0 * PI) / t_gate * n * t);
				n += 1;
			} else { // odd index (cos term)
				_pulse += pulse[idx][i] * cos((2.0 * PI) / t_gate * n * t);
			}
		}

		return rabi_freq_ratio[h_cfg.interacting_ions[idx].idx] * pulse_cfg.scale_am * _pulse;

	} else if (pulse_cfg.profile_am == 3) { // Custom Fourier series

		double _pulse = pulse[idx][0];

		// cos(2.0 * PI * f - phi)
		_pulse += pulse[idx][1] * cos((2.0 * PI) * pulse[idx][2] * t - pulse[idx][3]);
		_pulse += pulse[idx][4] * cos((2.0 * PI) * pulse[idx][5] * t - pulse[idx][6]);
		_pulse += pulse[idx][7] * cos((2.0 * PI) * pulse[idx][8] * t - pulse[idx][9]);

		// cos(2.0 * PI * f)
		// // _pulse += pulse[idx][1] * cos((2.0 * PI) * pulse[idx][2] * t);
		// // _pulse += pulse[idx][3] * cos((2.0 * PI) * pulse[idx][4] * t);
		// // _pulse += pulse[idx][5] * cos((2.0 * PI) * pulse[idx][6] * t);

		// cos term
		// _pulse += pulse[idx][1] * cos((2.0 * PI) * pulse[idx][3] * t);
		// sin term
		// _pulse += pulse[idx][2] * sin((2.0 * PI) * pulse[idx][3] * t);

		// cos term
		// _pulse += pulse[idx][4] * cos((2.0 * PI) * pulse[idx][6] * t);
		// sin term
		// _pulse += pulse[idx][5] * sin((2.0 * PI) * pulse[idx][6] * t);

		// cos term
		// _pulse += pulse[idx][7] * cos((2.0 * PI) * pulse[idx][9] * t);
		// sin term
		// _pulse += pulse[idx][8] * sin((2.0 * PI) * pulse[idx][9] * t);

		return rabi_freq_ratio[h_cfg.interacting_ions[idx].idx] * pulse_cfg.scale_am * _pulse;

	} else if (pulse_cfg.profile_am == 4) { // Custom g(t) function (Used only by Rabi dependent jump operators)

		double _pulse_rms = 0.0;
		for (int i = 0; i < h_cfg.n_gate_ions; i++) {

			for (int j = 0; j < g_data[0].size(); j++) { // n_list = g_data[0], a_list = g_data[1]
				_pulse_rms += h_cfg.interacting_ions[idx].crosstalk_scale_factors[i] * (pulse_cfg.scale_am * g_data[1][j]) * (pulse_cfg.scale_am * g_data[1][j]) / 2.0;
			}

		}

		return rabi_freq_ratio[h_cfg.interacting_ions[idx].idx] * _pulse_rms;
	} else if (pulse_cfg.profile_am == 5) { // Sampled \Omega(t)

		steps_AM = pulse_sampled.size();
		idx_AM = (int) (t * (steps_AM / t_gate));

		if (idx_AM == steps_AM) {
			idx_AM -= 1;
		}

		//	if (pulse_cfg.SYM_AM) {
		//		if (idx_AM >= pulse_cfg.steps_am) {
		//			idx_AM = pulse_cfg.steps_am - idx_AM % pulse_cfg.steps_am - 1;
		//		}
		//	}

		double _pulse = 0.0;
		for (int i = 0; i < h_cfg.n_gate_ions; i++) {
			try {
				_pulse += h_cfg.interacting_ions[idx].crosstalk_scale_factors[i] * pulse_sampled.at(idx_AM);
			}
			catch (std::out_of_range) {
				// // std::cout << "std::out_of_range" << std::endl;
				_pulse += h_cfg.interacting_ions[idx].crosstalk_scale_factors[i] * pulse_sampled.at(idx_AM - 1);
			}

		}
		return rabi_freq_ratio[h_cfg.interacting_ions[idx].idx] * pulse_cfg.scale_am * _pulse;

	} else if (pulse_cfg.profile_am == 6) { // Fourier series (cosine)

		double _pulse = pulse[idx][0];

		int n = 1;
		for (int i = 1; i < pulse_cfg.steps_am; i++) {

			if (i % 2 != 0) { // odd index amplitudes, even phases
				_pulse += pulse[idx][i] * cos((2.0 * PI) / t_gate * n * t - pulse[idx][i + 1]);
				n += 1;
			}
		}

		return rabi_freq_ratio[h_cfg.interacting_ions[idx].idx] * pulse_cfg.scale_am * _pulse;

	} else if (pulse_cfg.profile_am == 7) { // EIT cooling
		return pulse[idx][0];
	} else {
        std::cout<<"Unknown pulse shape"<<std::endl;
        return 0.0;
    }
}

double ion_trap::get_shelving_phase(double t){
    int steps_FM = (int) shaving_phase[0].size();

    int idx_FM = 0;
    double dt_FM = 0.0;
    idx_FM = (int) (t * (steps_FM / t_gate));
//    std::cout<<idx_FM<<std::endl;
    return shaving_phase[0][idx_FM];
}

double ion_trap::h_cos(int idx, double t) {
//    std::cout<<t<<std::endl;

	std::vector<double> crosstalk_scale_factors = h_cfg.interacting_ions[idx].crosstalk_scale_factors;
	std::vector<double>::iterator max_element = std::max_element(crosstalk_scale_factors.begin(), crosstalk_scale_factors.end());
	int arg_max = std::distance(crosstalk_scale_factors.begin(), max_element);

    if (crosstalk_scale_factors[arg_max] >= 1.0)
        // // return cos((delta[arg_max][idx_FM] + get_stark_shift(i, t)) * t);
        return cos((get_delta_new(t) + get_stark_shift(idx, t)) * t);
    else
        // // return cos((delta[0][idx_FM] + get_stark_shift(i, t)) * t);
        return cos((get_delta_new(t) + get_stark_shift(idx, t)) * t);

// use if shaving_phase imbalance is used
//	if (crosstalk_scale_factors[arg_max] >= 1.0){
//        double t_0 = std::fmod(t , 15.0);
//        double t_edge = 0.15;
//
////        std::cout<< t <<"  " << tt <<std::endl;
//		if (t_0 < t_edge){
////            std::cout << t << std::endl;
//            return cos((get_delta_new(t) + get_stark_shift(idx, t) + 0.05) * t + 0.1);}
//        else if(std::abs(t_0-5) < t_edge){
////            std::cout << t << std::endl;
//            return cos((get_delta_new(t) + get_stark_shift(idx, t) - 0.05) * t);
//        }
//        else
//            return cos((get_delta_new(t) + get_stark_shift(idx, t)) * t);}
//	else{
//		// // return cos((delta[0][idx_FM] + get_stark_shift(i, t)) * t);
//        if (t < 50.0){
//            std::cout<<"here"<<std::endl;
//            return cos((get_delta_new(t) + get_stark_shift(idx, t) + 0.0) * t);}
//        else
//            return cos((get_delta_new(t) + get_stark_shift(idx, t)) * t);}
}

double ion_trap::get_delta_new(double t) {
	int steps_FM = pulse_cfg.steps_fm * (pulse_cfg.SYM_FM ? 2.0 : 1.0);

	int idx_FM = 0;
	double dt_FM = 0.0;

	double tt = 0.0;
	if (pulse_cfg.profile_fm == 0) { // piecewise constant
		idx_FM = (int) (t * (steps_FM / t_gate));
	} else if (pulse_cfg.profile_fm == 1) { // error function
		dt_FM = t_gate / pulse_cfg.steps_fm;
		tt = t - dt_FM / 2;

		idx_FM = (int) (tt * (steps_FM / t_gate));
	} else if (pulse_cfg.profile_fm == 2) { // Fourier series

		double _delta = delta[0][0];

		int n = 1;
		for (int i = 1; i < pulse_cfg.steps_fm; i++) {

			if (i % 2 == 0) { // even index (sin term)
				_delta += delta[0][i] * sin((2.0 * PI) / t_gate * n * t);
				n += 1;
			} else { // odd index (cos term)
				_delta += delta[0][i] * cos((2.0 * PI) / t_gate * n * t);
			}
		}

		return _delta;
	}

	if (idx_FM == steps_FM) {
		idx_FM -= 1;
	}

	if (pulse_cfg.SYM_FM) {
		if (idx_FM >= pulse_cfg.steps_fm) {
			idx_FM = pulse_cfg.steps_fm - idx_FM % pulse_cfg.steps_fm - 1;
		}
	}

	if (pulse_cfg.profile_fm == 0) { // piecewise constant
//		double _pulse = 0.0;
//		for (int i = 0; i < h_cfg.n_gate_ions; i++) {
//			_pulse += h_cfg.interacting_ions[idx].crosstalk_scale_factors[i] * pulse[i][idx_FM];
//		}
//		return pulse_steps.scale * _pulse;
//        std::cout<<delta[0][idx_FM]<<std::endl;

		return delta[0][idx_FM];

    } else if (pulse_cfg.profile_fm == 1) { // error function
		double t_idx_FM_begin = idx_FM * dt_FM;
		double t_idx_FM_end = (idx_FM + 1.0) * dt_FM;

		if (!pulse_cfg.SYM_FM) {
			if (idx_FM == steps_FM - 1)
				return delta[0][idx_FM];
			else
				return (delta[0][idx_FM] + delta[0][idx_FM + 1]) / 2.0 + (delta[0][idx_FM + 1] - delta[0][idx_FM]) / 2 * erf(5.0 / dt_FM * (tt - (t_idx_FM_begin + t_idx_FM_end) / 2.0));
		} else {
			return delta[0][idx_FM];
		}
	} else {
        std::cout<<"Unknown frequency modulation function"<<std::endl;
        return 0.0;
    }
}

void ion_trap::initialize_H() {

    switch (h_cfg.mode) {
        case 0: // Raman transition
		// Todo: Add the new g(t) implementation
		{
            /*
            * Time dependencies
            * */

			// Hamiltonian generator
			std::vector<qsd::ComplexFunction> exp_iomegat_p(h_cfg.n_ions);
			std::vector<qsd::ComplexFunction> exp_iomegat_m(h_cfg.n_ions);

			std::vector<qsd::ComplexFunction> exp_imut_m(h_cfg.n_ions);

			std::vector<qsd::RealFunction> omega(h_cfg.n_interacting_ions);
			// // std::vector<qsd::RealFunction> h_cos(h_cfg.n_interacting_ions);

			for (int k = 0; k < h_cfg.n_ions; k++) {
				exp_iomegat_p[k] = [&, k](double t) {return cos(nu_list[k] * t) + IM * sin(nu_list[k] * t);};
				exp_iomegat_m[k] = [&, k](double t) {return cos(nu_list[k] * t) + mIM * sin(nu_list[k] * t);};

				exp_imut_m[k] = [&](double t) {return cos(get_delta_new(t) * t) + mIM * sin(get_delta_new(t) * t);};
			}

			for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
				omega[i] = [&, i](double t) {
				  return get_omega(i, t);
				};
			}

			// Calculating \beta_j = \sum^N_{k=1} \eta_{jk} (a_k \exp{-i\omega_k t} + a^\dagger_k \exp{+i\omega_k t})
			std::vector<qsd::Operator> beta(h_cfg.n_interacting_ions);
			for (int j = 0; j < h_cfg.n_interacting_ions; j++) {
				for (int k = 0; k < h_cfg.n_ions; k++) {
					if (k == 0) {
						beta[j] = eta_list[k][h_cfg.interacting_ions[j].idx] * (a[k] * exp_iomegat_m[k] + ad[k] * exp_iomegat_p[k]);
					} else {
						beta[j] += eta_list[k][h_cfg.interacting_ions[j].idx] * (a[k] * exp_iomegat_m[k] + ad[k] * exp_iomegat_p[k]);
					}
				}
			}

			for (int j = 0; j < h_cfg.n_interacting_ions; j++) {

				qsd::Operator _H;

				// Carrier transition
				if (h_cfg.carrier_status) {
					_H = sp[j];
				}

				std::complex<double> I_n = IM;
				int n_factorial = 1;
				qsd::Operator beta_n = beta[j];
				for (int n = 1; n <= h_cfg.expansion_order; n++) {

					if (!h_cfg.carrier_status) {
						if (n == 1) {
							_H = (I_n/double(n_factorial))*beta_n*sp[j];
						} else {
							_H += (I_n/double(n_factorial))*beta_n*sp[j];
						}
					} else {
						_H += (I_n/double(n_factorial))*beta_n*sp[j];
					}

					I_n *= IM;
					n_factorial *= (n + 1);
					beta_n = beta_n * beta[j];
				}

				_H *= 0.5;
				_H *= omega[j];
				_H *= exp_imut_m[j];

				if (j == 0) {
					H = _H;
				} else {
					H += _H;
				}
			}

			H += H.hc();

			break;
        }
        case 1: // MS gate
		{
            /*
             * Time dependencies
             * */

            // Hamiltonian generator
            std::vector<qsd::ComplexFunction> exp_iomegat_p(h_cfg.n_ions);
            std::vector<qsd::ComplexFunction> exp_iomegat_m(h_cfg.n_ions);
            std::vector<qsd::ComplexFunction> exp_iomegat_p_s(h_cfg.n_ions);
            std::vector<qsd::ComplexFunction> exp_iomegat_m_s(h_cfg.n_ions);
            std::vector<qsd::ComplexFunction> exp_iomegat_c_p(h_cfg.n_ions);
            std::vector<qsd::ComplexFunction> exp_iomegat_c_m(h_cfg.n_ions);
            Rabi_s = 0.2;
            Delta_L = 10;

            int MS_ions;

            if (h_cfg.monitoring_status){MS_ions = h_cfg.n_interacting_ions-1;} //interacting ions is greater by one to ensure proper initialization of operators (sx,sy, etc.) for spectator ion
            else{MS_ions =h_cfg.n_interacting_ions;}

            std::vector<qsd::RealFunction> g(MS_ions);
            std::vector<qsd::RealFunction> g_fluc(MS_ions);
            // // std::vector<qsd::RealFunction> h_cos(h_cfg.n_interacting_ions);

            for (int k = 0; k < h_cfg.n_ions; k++) {
                exp_iomegat_p[k] = [&, k](double t) {return cos(nu_list[k] * t) + IM * sin(nu_list[k] * t);};
                exp_iomegat_m[k] = [&, k](double t) {return cos(nu_list[k] * t) + mIM * sin(nu_list[k] * t);};
            }

            for (int k = 0; k < h_cfg.n_ions; k++) {
                exp_iomegat_p_s[k] = [&, k](double t) {return cos((nu_list[k]-TWOPI*Delta_L) * t) + IM * sin((nu_list[k] - TWOPI * Delta_L) * t);};
                exp_iomegat_m_s[k] = [&, k](double t) {return cos((nu_list[k]-TWOPI*Delta_L) * t) + mIM * sin((nu_list[k] - TWOPI * Delta_L) * t);};
            }

            for (int k = 0; k < h_cfg.n_ions; k++) {
                exp_iomegat_c_p[k] = [&, k](double t) {return cos((TWOPI*Delta_L) * t) + IM * sin((TWOPI * Delta_L) * t);};
                exp_iomegat_c_m[k] = [&, k](double t) {return cos((TWOPI*Delta_L) * t) + mIM * sin((TWOPI * Delta_L) * t);};
            }

            for (int i = 0; i < MS_ions; i++) {

                g[i] = [&, i](double t) {
                    return get_g(i, t);
                };

				g_fluc[i] = [&, i](double t) {
				  return get_g_fluc(i, t);
				};
            }

            // Calculating \beta_j = \sum^N_{k=1} \eta_{jk} (a_k \exp{-i\omega_k t} + a^\dagger_k \exp{+i\omega_k t})
            std::vector<qsd::Operator> beta(MS_ions);
            for (int j = 0; j < MS_ions; j++) {
                for (int k = 0; k < h_cfg.n_ions; k++) {
                    if (k == 0) {
                        beta[j] = eta_list[k][h_cfg.interacting_ions[j].idx] * (a[k] * exp_iomegat_m[k] + ad[k] * exp_iomegat_p[k]);
                    } else {
                        beta[j] += eta_list[k][h_cfg.interacting_ions[j].idx] * (a[k] * exp_iomegat_m[k] + ad[k] * exp_iomegat_p[k]);
                    }
                }
            }

            std::vector<qsd::Operator> beta_s(1);
            for (int k = 0; k < h_cfg.n_ions; k++) {
                if (k == 0) {
                    beta_s[0] =  IM *eta_list[k][h_cfg.spectator] * (a[k] * exp_iomegat_m_s[k] + ad[k] * exp_iomegat_p_s[k]);
                } else {
                    beta_s[0] +=  IM *eta_list[k][h_cfg.spectator] * (a[k] * exp_iomegat_m_s[k] + ad[k] * exp_iomegat_p_s[k]);
                }
            }

            for (int j = 0; j < MS_ions; j++) {

                qsd::Operator _H;
				qsd::Operator _H_fluc;

                // Carrier transition
                if (h_cfg.carrier_status) {
                    _H = sy[j];
                }

//                std::complex<double> I_n{0,1};
                std::complex<double> I_n = IM;
                int n_factorial = 1;
                qsd::Operator beta_n = beta[j];
                for (int n = 1; n <= h_cfg.expansion_order; n++) {

                    if (!h_cfg.carrier_status) {
						if (n == 1) {
							_H = (I_n / double(n_factorial)) * beta_n * (n % 2 ? (mIm * sx[j]) : sy[j]); // n is odd ? (if odd) : (if even)
//                            std::cout<< n <<std::endl;
						} else {
							_H += (I_n / double(n_factorial)) * beta_n * (n % 2 ? (mIm * sx[j]) : sy[j]); // n is odd ? (if odd) : (if even)
						}
                    } else {
                        _H += (I_n / double(n_factorial)) * beta_n * (n % 2 ? (mIm * sx[j]) : sy[j]); // n is odd ? (if odd) : (if even)
                    }

                    I_n *= IM;
                    n_factorial *= (n + 1);
                    beta_n = beta_n * beta[j];
                }

				_H_fluc = _H;
				_H_fluc *= g_fluc[j];

                _H *= g[j];
				// // _H *= h_cos[j];

                if (j == 0) {
                    H = _H;
					H_fluc = _H_fluc;
                } else {
                    H += _H;
					H_fluc += _H_fluc;
                }
            }

            // Kerr
            if (h_cfg.kerr_status) {
                // H += N_1 * (N_2 + N_3 + N_4 + N_5)
                H += TWOPI * 50e-4 * N[1] * (N[2]);
//                H += TWOPI * 3e-4 * N[2] * (N[2]);
//                H += TWOPI * 3e-4 * N[1] * (N[2] + N[3] + N[4]);
//                H += TWOPI * 3e-4 * N[2] * (N[3] + N[4]);
//                H += TWOPI * 3e-4 * N[3] * (N[4]);
            }

            // Adding the spectator ions for continuous monitoring
            if (h_cfg.monitoring_status){
//                std::cout<<"sd";
                H += TWOPI * 0.5 * Rabi_s * ( sp[h_cfg.spectator] * exp_iomegat_c_m[h_cfg.spectator] + sm[h_cfg.spectator] * exp_iomegat_c_p[h_cfg.spectator] );
                H += TWOPI * 0.5 * Rabi_s * ( sp[h_cfg.spectator] * beta_s[0] + sm[h_cfg.spectator] * beta_s[0].hc() );
            }

            break;
        }
		case 2: // Debug
		{
			/*
			* Time dependencies
			* */

			// Omega(t)
			std::vector<qsd::RealFunction> omega(h_cfg.n_interacting_ions);
			for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
				omega[i] = [&, i](double t) {
				  return get_omega(i, t);
				};
			}

			// blue sideband transition
			H = sp[0];
			H *= omega[0];

			break;
		}
        case 3: // Cooling dynamics
        {
            /*
            * Time dependencies
            * */

            // Omega(t)
            std::vector<qsd::RealFunction> omega(2);
            std::vector<qsd::RealFunction> detuning_12(h_cfg.n_ions);
            qsd::RealFunction det;

//            int processor = mode_cfg.get<int>("mode-8.processor");


            detuning_12[0] = [&](double t) {
                return (get_delta_new(0));
            };

            detuning_12[1] = [&](double t) {
                return (get_delta_new(t) - get_Delta_pump(0));
            };

            det = [&](double t) {
//                return (-get_Delta_pump(0) + get_delta_new(t));
                return (get_delta_new(t));
            };

            for (int k = 0; k < 2; k++){
                omega[k] = [&, k](double t) {
                    return get_EIT_pulse(k);
                };
            }

            for (int j = 0; j < h_cfg.n_interacting_ions; j++) {

                qsd::Operator _H;

                _H = omega[1]*S_02[j] + omega[0]*S_12[j];
                _H += eta_list[0][j] * Im* (a[j] + ad[j]) * S_02[j]*omega[1] + eta_list[0][j] * Im*(a[j] + ad[j]) * S_12[j]*omega[0];
                _H -= det*S_22[j] + detuning_12[1]* S_11[j];
                if (j == 0) {
                     H = _H;
                } else {
                    H += _H;
                }
            }

            H += H.hc();
//            H += (ad[0]*a[0])*nu_list[0] + id_m[0]*nu_list[0]*0.5 + 0.0*(1.0*nu_list[1] + 0.5)*ad[1]*a[1]; incorrect scaling of vacuum fluctuations to Mrad/s

            break;
        }
        case 4: // Shelving/deshelving fidelities
        {
            int omega_seg;

            /*
            * Time dependencies
            * */
//            std::cout<<"12"<<std::endl;

            // Omega(t)
            if (shelving_flag==1){omega_seg = 2;} else {omega_seg = num_pulses;}

            std::vector<qsd::RealFunction> omega(omega_seg);
            std::vector<qsd::ComplexFunction> exp_imut_m(h_cfg.n_ions);
            std::vector<qsd::RealFunction> energy_drift_SD(h_cfg.n_ions);
            std::vector<qsd::RealFunction> costheta(1);
            std::vector<qsd::RealFunction> sintheta(1);

            costheta[0] = [&](double t) {
                return cos(3.141592653589793238 * get_shelving_phase(t));
            };

            sintheta[0] = [&](double t) {
                return sin(3.141592653589793238 * get_shelving_phase(t));
            };

            for (int i = 0; i < h_cfg.n_ions; ++i) {
                energy_drift_SD[i] = [&](double t) {
                    return TWOPI*get_energy_drift(i,t);
                };
            }

            for (int k = 0; k < omega_seg; k++) {
                omega[k] = [&, k](double t) {
                    return get_shelving_pulse(k,t, shelving_flag);
                };
            }

            qsd::Operator _H;
            qsd::Operator _H_fluc, _H_fluc_1, _Hmot;

            if (shelving_flag==1){
                H = 0.5*( (sigma_X[0]+sigma_X[1])*costheta[0] + (sigma_Y[0]+sigma_Y[1])*sintheta[0] ) * omega[0];
            } else {
                for (int j = 0; j < num_pulses; j++) {

                    _H = 0.5*(pulse_seq_0[j] + pulse_seq_1[j]) * omega[j];
    //                _H_fluc = 0.5*(pulse_seq_0[j]) * omega[j];
    //                _H_fluc_1 = 0.5*(pulse_seq_1[j]) * omega[j];

                    if ( j == 0){
                        H = _H;
//                        H += _Hmot + _Hmot.hc();
    //                    H_fluc = _H_fluc;
    //                    H_fluc_1 = _H_fluc_1;
                    } else {
                        H += _H;
//                        H += _Hmot + _Hmot.hc();
    //                    H_fluc += _H_fluc;
    //                    H_fluc_1 += _H_fluc_1;
                    }
                }
    //            std::cout<<energy_drift<<std::endl;

//                H += (D0[0] + D1[0] - (S0[0] + S1[0]))*(energy_drift_SD[0]);
                _Hmot = IM * (a[0] + ad[0] ) * (S0_D0[0] + S1_D1[0]) ;
                _Hmot *= 0.5*0.2* 0.10134;
                H +=0.14*TWOPI*( N[0] ) + _Hmot + _Hmot.hc();
    //            H += H.hc();
            }
            break;
        }
    }
}

void ion_trap::set_state() {
    /*
    * The initial state
    * */
    std::vector<int>().swap(state_cfg.dynamic_degrees);

    std::vector<qsd::State> psilist;

    if ( (simulation_mode < 7) ) {
//        std::cout << "sim mode:" << simulation_mode<<std::endl;

        for (int i = 0; i < nu_list.size() +h_cfg.n_interacting_ions; i++) {
            if (i < h_cfg.n_interacting_ions) {
                psilist.emplace_back(2, qsd::SPIN);
            } else {
                state_cfg.dynamic_degrees.emplace_back(i);
                psilist.emplace_back(state_cfg.phonon_cutoffs[i - h_cfg.n_interacting_ions],
                                     state_cfg.motion[i - h_cfg.n_interacting_ions], qsd::FIELD);
            }
        }
            psi0 = qsd::State(psilist.size(), psilist.data());

        for (int i = 0; i < h_cfg.n_gate_ions; i++) {
            if (state_cfg.spin[i] == 1) {
                psi0 *= sp[h_cfg.gate_ions_idx[i]];
            }
        }
    } else if (simulation_mode == 7){

//        std::cout << "here..."<<simulation_mode <<std::endl;
        for (int i = 0; i < 2*h_cfg.n_interacting_ions; i++) {   //    for (int i = 0; i < nu_list.size() + h_cfg.n_interacting_ions; i++) {
            if (i < h_cfg.n_interacting_ions) {
                psilist.emplace_back(3,qsd::ATOM);

            } else {
                state_cfg.dynamic_degrees.emplace_back(i);
                psilist.emplace_back(state_cfg.phonon_cutoffs[i - h_cfg.n_interacting_ions],
                                     state_cfg.motion[i - h_cfg.n_interacting_ions], qsd::FIELD);
            }
        }
        psi0 = qsd::State(psilist.size(), psilist.data());
//        psi0 *= sig_01;
    } else if(simulation_mode == 8) {
        double _phi =  3.141592653589793238/4.0;
        std::complex<double> imag{0,1};
        std::complex<double> angle = cos(_phi) + imag*sin(_phi);

        qsd::TransitionOperator Aux0{0,0,0};
        qsd::TransitionOperator Aux1{1,1,0};
        qsd::TransitionOperator Aux2{1,0,0};
        qsd::TransitionOperator Aux3{0,1,0};
        qsd::TransitionOperator Aux_p{1,1,0};

        for (int i = 0; i < nu_list.size() +h_cfg.n_interacting_ions; i++) {
            if (i < h_cfg.n_interacting_ions) {
                psilist.emplace_back(4,qsd::ATOM);
            } else {
                state_cfg.dynamic_degrees.emplace_back(i);
                psilist.emplace_back(state_cfg.phonon_cutoffs[i - h_cfg.n_interacting_ions],
                                     state_cfg.motion[i - h_cfg.n_interacting_ions], qsd::FIELD);
            }
        }

        psi0 = qsd::State(psilist.size(), psilist.data());
        psi0 *= S0_S1[1];

        if (prep_state==1) {
            Aux0 *= 1;
            Aux1 *= -1;
            Aux1 += Aux0;
            Aux1 *= sqrt(0.5);

            Aux3 += Aux2;
//            Aux0 += S0_S1[0].hc();
//            Aux0 = S0_S1[0];
            Aux3 *= sqrt(0.5);
            Aux3 += Aux1;
            Aux_p *= angle;
            Aux_p += Aux0;

//            std::cout << prep_sequence << std::endl;
            std::stringstream string_stream(prep_sequence);
            std::string substr;

            while (string_stream.good()) {
                std::getline(string_stream, substr, '_');
                if (substr == "H"){
//                    psi0 *= Aux0;
                    psi0 *= Aux3;
                    psi0 *= Aux_p;
                }
            }
        }
    }
}

void ion_trap::set_spin_state(std::vector<qsd::State> spin_state){
    /*
    * The initial state
    * */

    std::vector<int>().swap(state_cfg.dynamic_degrees);
//
    std::vector<qsd::State> psilist;

    for (int i = 0; i < nu_list.size() + h_cfg.n_interacting_ions; i++) {
        if (i < h_cfg.n_interacting_ions) psilist.push_back(spin_state[i]);
        else {
            state_cfg.dynamic_degrees.emplace_back(i);
            psilist.emplace_back(state_cfg.phonon_cutoffs[i - h_cfg.n_interacting_ions], state_cfg.motion[i - h_cfg.n_interacting_ions], qsd::FIELD);
        }
    }
//    psilist.emplace_back(3,qsd::ATOM);

    psi0 = qsd::State(psilist.size(), psilist.data());
}

void ion_trap::initialize_L() {

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

    switch (l_cfg.mode) {
        case 1: { // Full

            /*
            * Time dependencies (Lindblad operators time-dependent coefficients)
            * */

            // Motional dephasing
            /*
            qsd::RealFunction L_md_t = [&](double t) {
                return 0.0;
            };
            */

			qsd::RealFunction gamma_amp = [&](double t) {
				// Todo: Double check for laser spill-over effect
			  return sqrt(l_cfg.amp_fluct);
			};

            // Spin dephasing (gamma_Rayleigh)
            std::vector<qsd::RealFunction> L_sd_t(h_cfg.n_interacting_ions);
            for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
                L_sd_t[i] = [&, i](double t) {
                    return sqrt(std::abs(get_omega(i, t)) * l_cfg.gamma_Rayleigh / 4.0);
                };
            }

            // Spin flip (gamma_Raman)
            std::vector<qsd::RealFunction> L_sf_t(h_cfg.n_interacting_ions);
            for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
                L_sf_t[i] = [&, i](double t) {
                    return sqrt(std::abs(get_omega(i, t)) * l_cfg.gamma_Raman);
                };
            }

            std::vector<qsd::RealFunction> pink_noise(h_cfg.n_interacting_ions);
            for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
                pink_noise[i] = [&, i](double t) {
                    return sqrt(std::exp( -t* l_cfg.pink_cumm_1) );
                };
            }
            // Laser power fluctuations
			// // std::vector<qsd::RealFunction> L_sf_t(h_cfg.n_interacting_ions);

            // // std::vector<qsd::ComplexFunction> H0_t(h_cfg.n_interacting_ions);
            // // std::vector<std::vector<qsd::ComplexFunction>> H_t_1(h_cfg.n_interacting_ions);
            // // std::vector<std::vector<qsd::ComplexFunction>> H_t_2(h_cfg.n_interacting_ions);

            /*
            * Lindblad operators
            **/

            // Motional Heating
            if (l_cfg.motional_heating) {
                for (int i = 0; i < h_cfg.n_ions; i++) {
                    L.emplace_back(sqrt(l_cfg.gamma_heating[i]+0*1e-6) * a[i]);  // Heating of the ith motional mode
                    L.emplace_back(sqrt(l_cfg.gamma_heating[i]) * ad[i]); // ...
                }
            }

            // Motional Dephasing
            if (l_cfg.motional_dephasing) {
                for (int i = 0; i < h_cfg.n_ions; i++) {
                    L.emplace_back(sqrt(2.0 / l_cfg.motional_coherence_time[i]) * (ad[i] * a[i])); // Dephasing of the ith motional mode
                }
            }

            // Spin Dephasing (Rayleigh)
            if (l_cfg.spin_dephasing_rayleigh) {
                for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
                    L.emplace_back(sz[i] * L_sd_t[i]);  // Dephasing of the ith ion
                }
            }

            // Spin Flip (Raman)
            if (l_cfg.spin_flip_raman) {
                for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
                    L.emplace_back(sm[i] * L_sf_t[i]);  // Spin Flip of the ith ion
                    L.emplace_back(sp[i] * L_sf_t[i]);  // ...
                }
            }

            // Laser dephasing (Laser coherence time)
            if (l_cfg.laser_dephasing) {
                for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
                    L.emplace_back(sqrt(1.0 / l_cfg.laser_coherence_time) * sz[i]);  // Dephasing of the ith ion
                }
            }

            //Position operator measurement
            if (l_cfg.position_operator) {
                for (int i = 0; i < h_cfg.n_ions; i++) {
                    L.emplace_back(sqrt(l_cfg.position_meas_coupling) * eta_list[h_cfg.spectator][h_cfg.spectator] * x[i]);  // Measurement of X of the ith ion
                }
            }

            // Laser power fluctuations
            if (l_cfg.laser_power_fluctuation) {

				/*
                for (int i = 0; i < h_cfg.n_interacting_ions; i ++) {
                    L.emplace_back(sqrt(l_cfg.amp_fluct) * sy[i] * H0_t[i]);
                }
                */

				/*
                for (int k = 0; k < h_cfg.n_ions; k++) {
                    for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
                        L.emplace_back(sqrt(l_cfg.amp_fluct) * ad[k] * sx[i] * H_t_1[i][k]);
                        L.emplace_back(sqrt(l_cfg.amp_fluct) * a[k] * sx[i] * H_t_2[i][k]);
                    }
                }
                */

				// // L.emplace_back(sqrt(l_cfg.amp_fluct) * H); // Todo: this implementation does not work with cross-Kerr coupling enabled
            	L.emplace_back(gamma_amp * H_fluc); // Todo: this implementation does not work with cross-Kerr coupling enabled
			}

// ****************************************** EIT dynamics Lindbladians ******************************************

            // Spontaneous decay from upper to ground state
            if (l_cfg.spontaneous_decay_01) {
                for (int i = 0; i < 1; i++) {
                    L.emplace_back(sqrt(l_cfg.gamma_01) * S_02[i].hc()); // Spontaneous decay from upper to ground state
                }
            }

//            Spontaneous decay from upper to lower state
            if (l_cfg.spontaneous_decay_12) {
                for (int i = 0; i < 1; i++) {
                    L.emplace_back(sqrt(l_cfg.gamma_12) * S_12[i].hc()); // Spontaneous decay from upper to lower state
                }
            }

// ****************************************** Shelving/deshelving dynamics Lindbladians ******************************************

            // Laser dephasing (Laser coherence time)
            if (l_cfg.laser_dephasing_S) {
                L.emplace_back(sqrt(1.0 / l_cfg.laser_coherence_time) * (D0[0] - S0[0] + (D1[0] - S1[0])) );  // Dephasing of the ith ion
//                L.emplace_back(sqrt(1.0 / l_cfg.laser_coherence_time) * (S1[0] - D1[0]));
            }

            if (l_cfg.laser_power_fluc_2Atom){
                L.emplace_back(gamma_amp * (H_fluc + H_fluc_1));
//                L.emplace_back(gamma_amp * (H_fluc_1));
            }

            if (l_cfg.pink_noise){
                L.emplace_back((D0[0] - S0[0] + (D1[0] - S1[0]) ) * pink_noise[0] );//                L.emplace_back(gamma_amp * (H_fluc_1));
            }

            // blue sideband
			// Omega(t)
			/*
			std::vector<qsd::RealFunction> omega(h_cfg.n_interacting_ions);
			for (int i = 0; i < h_cfg.n_interacting_ions; i++) {
				omega[i] = [&, i](double t) {
				  return get_omega(i, t);
				};
			}
			*/

            break;
        }

        case 0: { // None
            break;
        }
    }
}

void ion_trap::set_outlist(std::vector<std::string> _outlist_string)
{
	std::map<std::string, qsd::Operator> _outlist;

	for (auto key : _outlist_string) {
		outlist[key] = outlist_processor(key);
	}
}

qsd::Expectation ion_trap::processor_fidelity(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {

	if (qsd_result.size() > 0) {

		std::vector<double> expec_0;
		std::vector<double> expec_1;

		int ntraj = qsd_result.size();

		// save_per_traj
		if (save_per_traj) {
			for (int i = 0; i < ntraj; i++) {

				double expec_0_mean = qsd_result[i].observables["rho_ideal_1"].back().mean.real();
				double expec_1_mean = qsd_result[i].observables["rho_ideal_2"].back().mean.real();

				if (expec_0_mean > expec_1_mean) {
					(*log_file) << std::setprecision(8) << expec_0_mean << std::endl;
				} else {
					(*log_file) << std::setprecision(8) << expec_1_mean << std::endl;
				}
			}
		}

		for (auto &traj_result : qsd_result) {
			expec_0.push_back(traj_result.observables["rho_ideal_1"].back().mean.real());
			expec_1.push_back(traj_result.observables["rho_ideal_2"].back().mean.real());
		}

		// compute the mean for expec_0 and expec_1
		double expec_0_mean = std::accumulate(expec_0.begin(), expec_0.end(), 0.0) / ntraj;
		double expec_1_mean = std::accumulate(expec_1.begin(), expec_1.end(), 0.0) / ntraj;

		// compute the variance for expec_0 and expec_1
		double expec_0_var = 0.0;
		double expec_1_var = 0.0;

		if (ntraj > 1) {
			expec_0_var = variance(expec_0);
			expec_1_var = variance(expec_1);
		}

		// return std::max(qsd_result.data[0].back()[0], qsd_result.data[1].back()[0]);
		if (expec_0_mean > expec_1_mean) {
			return qsd::Expectation(expec_0_mean, 0.0, expec_0_var, 0.0);
		} else {
			return qsd::Expectation(expec_1_mean, 0.0, expec_1_var, 0.0);
		}

	} else {
		return qsd::Expectation(0.0, 0.0, 0.0, 0.0);
	}
}

qsd::Expectation ion_trap::processor_infidelity(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {

	if (qsd_result.size() > 0) {

		std::vector<double> expec_0;
		std::vector<double> expec_1;

		int ntraj = qsd_result.size();

		for (auto &traj_result : qsd_result) {
			expec_0.push_back(1.0 - traj_result.observables["rho_ideal_1"].back().mean.real());
			expec_1.push_back(1.0 - traj_result.observables["rho_ideal_2"].back().mean.real());
		}

		// compute the mean for expec_0 and expec_1
		double expec_0_mean = std::accumulate(expec_0.begin(), expec_0.end(), 0.0) / ntraj;
		double expec_1_mean = std::accumulate(expec_1.begin(), expec_1.end(), 0.0) / ntraj;

		// compute the variance for expec_0 and expec_1
		double expec_0_var = 0.0;
		double expec_1_var = 0.0;

		if (ntraj > 1) {
			expec_0_var = standard_deviation(expec_0);
			expec_1_var = standard_deviation(expec_1);
		}

		if (expec_0_mean < expec_1_mean) {
			return qsd::Expectation(expec_0_mean, 0.0, expec_0_var, 0.0);
		} else {
			return qsd::Expectation(expec_1_mean, 0.0, expec_1_var, 0.0);
		}

	} else {
		return qsd::Expectation(1.0, 0.0, 0.0, 0.0);
	}

	/*
	if (!qsd_result.timeout && !qsd_result.error)
		return (1.0 - std::max(qsd_result.data[0].back()[0], qsd_result.data[1].back()[0]));
	else
		return 1.0;
	*/
}

std::vector<double> ion_trap::processor_parity(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {

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
	std::vector<std::ofstream> parity_files;

	for (int i = 0; i < qsd_result.size(); i++) {
		if (qsd_result.size() == 1) {
			parity_files.emplace_back(std::ofstream{"parity.txt"});
		} else {
			parity_files.emplace_back(std::ofstream{"parity_traj_" + std::to_string(i + 1) + ".txt"});
		}
		parity_files[i] << "phi, parity" << std::endl;
	}

	std::vector<double> fidelities;
	std::vector<std::complex<double>> _rho;

	double parity_contrast = 0.0;
	for (int i = 0; i < qsd_result.size(); i++) {

		///////////////////////////////////
		qsd::State psi_1 = qsd_result[i].state;
		qsd::State psi_2 = qsd_result[i].state;

		// Create more projector_z operators for more ions (those involved in the gate) and do the parity equation for more qubits
		// important: the small phi rotation has to be also done on those ions, meaning more rg_n

		// // for (int j = 0; j < pow(2, h_cfg.n_interacting_ions); j++){
		for (int j = 0; j < pow(2, h_cfg.n_gate_ions); j++){
			psi_2 = psi_1;

			if (findParity(j)==0){
				std::cout<<"Even Parity\n";
				// cout<< j <<endl;
				std::complex<double> aux;


				int x = 1;
				for (int k = 0; k < h_cfg.n_gate_ions; k++){
				// // for (int k = 0; k < h_cfg.n_interacting_ions; k++){
					// cout<< x <<endl;
					bool m = j & x;
					psi_2 *= projectors_z[h_cfg.gate_ions_idx[k]][3 * m];
					x <<= 1;
				}
				aux = psi_1 * psi_2;
				_rho.push_back(aux);
				std::cout << aux << std::endl;
			}
		}
		///////////////////////////////////

		// Global PI/2 rotation (Analyzer pulse) - Manning thesis page 110
		double phi = 0.0; // Analyzer pulse shaving_phase

		std::vector<double> parity;
		while (phi < TWOPI) {

			qsd::State psi_1 = qsd_result[i].state;
			qsd::State psi_2 = qsd_result[i].state;
			std::vector<qsd::Operator> rg;

			// // for (int i=0; i < h_cfg.n_interacting_ions; i++){
			for (int i=0; i < h_cfg.n_gate_ions; i++) {
				// cout << cos(phi) <<endl;
				qsd::Operator aux = cos(PIHALF / 2.0) * id_s[h_cfg.gate_ions_idx[i]] + mIM * sin(PIHALF / 2.0) * (cos(phi) * sx[h_cfg.gate_ions_idx[i]] + sin(phi) * sy[h_cfg.gate_ions_idx[i]]);
				rg.push_back(aux); // |0><0|
			}

			// // for (int i=0; i < h_cfg.n_interacting_ions; i++){
			for (int i=0; i < h_cfg.n_gate_ions; i++){
				psi_1 *= rg[i];
			}

			std::complex<double> _rho;
			double _parity = 0;
			std::complex<double> aux1;

			// // for (int w = 0; w < pow(2, h_cfg.n_interacting_ions); w++){
			for (int w = 0; w < pow(2, h_cfg.n_gate_ions); w++){
				psi_2 = psi_1;

				int x = 1;
				// // for (int j = 0; j < h_cfg.n_interacting_ions; j++){
				for (int j = 0; j < h_cfg.n_gate_ions; j++){

					bool k = w & x;
					psi_2 *= projectors_z[h_cfg.gate_ions_idx[j]][3 * k];
					x <<= 1;
					// cout<< 3 * k <<endl;
				}
				aux1 = psi_1 * psi_2;
				// cout<< aux1 <<endl;

				if (findParity(w) == 1){
					_parity += real(aux1);
				} else {
					_parity -= real(aux1);
				}
				// cout<< _parity <<endl;
			}

			parity_files[i] << phi << ", " << _parity << std::endl;
			parity.push_back(_parity);

			phi += phi_res;
		}

		auto max_parity = max_element(std::begin(parity), std::end(parity));
		auto min_parity = min_element(std::begin(parity), std::end(parity));
		double two_abs_rho_00_11 = (*max_parity - *min_parity) / 2;
		std::cout<< two_abs_rho_00_11  << std::endl;
		parity_contrast = two_abs_rho_00_11;

		std::complex<double> _rhoSum;
		_rhoSum = accumulate(_rho.begin(), _rho.end(), _rhoSum);
		// // fidelities.push_back(_rhoSum.real()/pow(2, h_cfg.n_interacting_ions - 1) + 0.5 * two_abs_rho_00_11);
		fidelities.push_back(_rhoSum.real()/pow(2, h_cfg.n_gate_ions - 1) + 0.5 * two_abs_rho_00_11);
		std::cout<< _rhoSum << std::endl;
	}

	// // return accumulate(fidelities.begin(), fidelities.end(), 0.0) / fidelities.size();

	// return the parity contrast
	return {parity_contrast, 0.0};
}

std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> ion_trap::processor_concatenated_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {
// [i][j][dt][ntraj]

/* computed matrix
|00><00| |00><01| |01><00| |01><01|
|00><10| |00><11| |01><10| |01><01|
|10><00| |10><01| |11><00| |11><01|
|10><10| |10><11| |11><10| |11><11|
*/

/* density_matrix
|00><00| |00><01| |00><10| |00><11|
|01><00| |01><01| |01><10| |01><11|
|10><00| |10><01| |10><10| |10><11|
|11><00| |11><01| |11><10| |11><11|
*/


	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> density_matrix_evol(4);
	if (qsd_result.size() > 0) {

		// // int ntraj = qsd_result.size();
		for (auto& traj_result: qsd_result) {
            std::vector<qsd::State> state_v = traj_result.state_v;

            for (int k = 0; k < state_v.size(); k++) {
//                density_matrix_expec[i].resize(4);
                for (int i = 0; i < 4; ++i) {
                    density_matrix_evol[i].resize(4);

                    for (int j = 0; j < 4; j++) {
                        density_matrix_evol[i][j].resize(state_v.size());
                        qsd::State psi_tmp_v, psi_v;

                        psi_v = state_v[k];
                        psi_tmp_v = psi_v;

                        psi_tmp_v *= projectors_z[h_cfg.gate_ions_idx[0]][i] * projectors_z[h_cfg.gate_ions_idx[1]][j];
                        density_matrix_evol[i][j][k].push_back(psi_v * psi_tmp_v);
                    }
                }
            }
        }

        std::vector<std::vector<std::complex<double>>> aux; //swap entries to match Density Matrix
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if ((i == 0 or i == 2) and j > 1) {
                    aux = density_matrix_evol[i + 1][j - 2];
                    density_matrix_evol[i + 1][j - 2] = density_matrix_evol[i][j];
                    density_matrix_evol[i][j] = aux;
//                            std::cout << "here" << std::endl;
                }
            }
        }
	}

	return density_matrix_evol;
}


std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> ion_trap::processor_reduced_concatenated_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {
// [i][j][dt][ntraj]
/* computed matrix
|00><00| |00><01| |01><00| |01><01|
|00><10| |00><11| |01><10| |01><01|
|10><00| |10><01| |11><00| |11><01|
|10><10| |10><11| |11><10| |11><11|
*/

/* density_matrix
|0><0| |0><1|
|1><0| |1><1|
*/


    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> density_matrix_evol(2);
    if (qsd_result.size() > 0) {

//		  int ntraj = qsd_result.size();
        for (auto& traj_result: qsd_result) {
            std::vector<qsd::State> state_v = traj_result.state_v;
//            std::cout<< "here2: " << state_v.size() << std::endl;

            for (int k = 0; k < state_v.size(); k++) {
                for (int i = 0; i < 2; ++i) {
                    density_matrix_evol[i].resize(2);

                    for (int j = 0; j < 2; j++) {
                        density_matrix_evol[i][j].resize(state_v.size());
                        qsd::State psi_tmp_v, psi_v;

                        psi_v = state_v[k];
                        psi_tmp_v = psi_v;
                        int index = i*1+j*2;

                        psi_tmp_v *= projectors_z[h_cfg.gate_ions_idx[0]][index];
//                        std::cout<<index<<std::endl;
                        density_matrix_evol[i][j][k].push_back(psi_v * psi_tmp_v);
                    }
                }
            }
        }

        std::vector<std::vector<std::complex<double>>> aux; //swap entries to match Density Matrix
        aux = density_matrix_evol[0][1];
        density_matrix_evol[0][1] = density_matrix_evol[1][0];
        density_matrix_evol[1][0] = aux;

    }


    return density_matrix_evol;
}

std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> ion_trap::processor_concatenated_multilevel_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, int atom, bool save_per_traj, std::ofstream* log_file) {
// [row][column][dt][ntraj]

/* computed matrix
|00><00| |00><01| |01><00| |01><01|
|00><10| |00><11| |01><10| |01><01|
|10><00| |10><01| |11><00| |11><01|
|10><10| |10><11| |11><10| |11><11|
*/

/* density_matrix
|00><00| |00><01| |00><10| |00><11|
|01><00| |01><01| |01><10| |01><11|
|10><00| |10><01| |10><10| |10><11|
|11><00| |11><01| |11><10| |11><11|
*/

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> density_matrix_evol(2);
    if (qsd_result.size() > 0) {

        // // int ntraj = qsd_result.size();
        for (auto& traj_result: qsd_result) {
            std::vector<qsd::State> state_v = traj_result.state_v;

            for (int k = 0; k < state_v.size(); k++) {
                for (int i = 0; i < 2; ++i) {
                    density_matrix_evol[i].resize(2);

                    for (int j = 0; j < 2; j++) {
                        density_matrix_evol[i][j].resize(state_v.size()); //resizing for time length
                        qsd::State psi_tmp_v, psi_v;

                        psi_v = state_v[k];
                        psi_tmp_v = psi_v;
                        int ind = i*1 + j*2;

                        if (atom == 0){
//                            psi_tmp_v *= proj_z0_basis[ind];
                            psi_tmp_v *= proj_z0_basis[ind];
                        } else if (atom == 1){
                            psi_tmp_v *= proj_z1_basis[ind];
                        } else {
//                            psi_tmp_v *= proj_z3_basis[ind];
                            std::cout<<"Atom number not implemented. Using Identity as projector"<<std::endl;
                            psi_tmp_v *= id;
                        }
//                        psi_tmp_v *= -1*S0_D0[0] + 1*S0_D0[0].hc();
//                        psi_tmp_v *= proj_z_basis[0][0] * proj_z_basis[1][0];

                        density_matrix_evol[i][j][k].push_back(psi_v * psi_tmp_v);
                    }
                }
            }
        }
    }

    return density_matrix_evol;
}


std::vector<std::vector<std::vector<qsd::Expectation>>> ion_trap::processor_average_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {

    /* density_matrix
    |00><00| |00><01| |01><00| |01><01|
    |00><10| |00><11| |01><10| |01><01|
    |10><00| |10><01| |11><00| |11><01|
    |10><10| |10><11| |11><10| |11><11|
    */

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> density_matrix_expec = processor_concatenated_density_matrix(qsd_result);
    std::vector<std::vector<std::vector<qsd::Expectation>>> density_matrix(4);
    int kmax = density_matrix_expec[0][0].size();
//    std::cout << kmax<< std::endl;


    if (qsd_result.size() > 0) {

        int ntraj = qsd_result.size();
//        std::cout<<ntraj<<std::endl;

        for (int i = 0; i < 4; i++) {
            density_matrix[i].resize(4);
            for (int j = 0; j < 4; ++j) {
                density_matrix[i][j].resize(kmax,qsd::Expectation(0.0, 0.0, 0.0, 0.0));
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < kmax; ++k) {

                    std::complex<double> expec_0_mean = std::accumulate(density_matrix_expec[i][j][k].begin(), density_matrix_expec[i][j][k].end(), std::complex<double>(0.0, 0.0))/double(ntraj);
                    std::complex<double> expec_0_var = std::complex<double>(0.0, 0.0);

//                    if (ntraj > 1) {
//                        expec_0_var = variance(density_matrix_expec[i][j][k]);
//                    }
                    density_matrix[i][j][k] = qsd::Expectation(expec_0_mean.real(), expec_0_mean.imag(), expec_0_var.real(), expec_0_var.imag());
                }
            }
        }

    }

    return density_matrix; //i & j: expectation value for measurement operators on target ions. k: time stamp evolution.
}

    std::vector<std::vector<std::vector<qsd::Expectation>>> ion_trap::processor_average_multilevel_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, int atom, bool save_per_traj, std::ofstream* log_file) {

        /* density_matrix
        |00><00| |00><01| |01><00| |01><01|
        |00><10| |00><11| |01><10| |01><01|
        |10><00| |10><01| |11><00| |11><01|
        |10><10| |10><11| |11><10| |11><11|
        */

        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> density_matrix_expec = processor_concatenated_multilevel_density_matrix(qsd_result, atom);
        std::vector<std::vector<std::vector<qsd::Expectation>>> density_matrix(2);
        int kmax = density_matrix_expec[0][0].size();
//    std::cout << kmax<< std::endl;

        if (qsd_result.size() > 0) {

            int ntraj = qsd_result.size();
//        std::cout<<ntraj<<std::endl;

            for (int i = 0; i < 2; i++) {
                density_matrix[i].resize(2);
                for (int j = 0; j < 2; ++j) {
                    density_matrix[i][j].resize(kmax,qsd::Expectation(0.0, 0.0, 0.0, 0.0));
                }
            }

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < kmax; ++k) {

                        std::complex<double> expec_0_mean = std::accumulate(density_matrix_expec[i][j][k].begin(), density_matrix_expec[i][j][k].end(), std::complex<double>(0.0, 0.0))/double(ntraj);
                        std::complex<double> expec_0_var = std::complex<double>(0.0, 0.0);

//                    if (ntraj > 1) {
//                        expec_0_var = variance(density_matrix_expec[i][j][k]);
//                    }
                        density_matrix[i][j][k] = qsd::Expectation(expec_0_mean.real(), expec_0_mean.imag(), expec_0_var.real(), expec_0_var.imag());
                    }
                }
            }

            // Print the density_matrix
            /*
                std::cout << "Density Matrix \n";
            std::cout << density_matrix[0][0][99].mean << density_matrix[0][1][99].mean << density_matrix[1][0][99].mean << density_matrix[1][1][99].mean << std::endl;
            std::cout << density_matrix[0][2][99].mean << density_matrix[0][3][99].mean << density_matrix[1][2][99].mean << density_matrix[1][3][99].mean << std::endl;
            std::cout << density_matrix[2][0][99].mean << density_matrix[2][1][99].mean << density_matrix[3][0][99].mean << density_matrix[3][1][99].mean << std::endl;
            std::cout << density_matrix[2][2][99].mean << density_matrix[2][3][99].mean << density_matrix[3][2][99].mean << density_matrix[3][3][99].mean << std::endl;
            */
        }
        return density_matrix; //i & j: expectation value for measurement operators on target ions. k: time stamp evolution.
    }

std::vector<std::vector<std::vector<qsd::Expectation>>> ion_trap::processor_reduced_average_density_matrix(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {

    /* density_matrix
    |00><00| |00><01| |01><00| |01><01|
    |00><10| |00><11| |01><10| |01><01|
    |10><00| |10><01| |11><00| |11><01|
    |10><10| |10><11| |11><10| |11><11|
    */

    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> density_matrix_expec = processor_reduced_concatenated_density_matrix(qsd_result);
    std::vector<std::vector<std::vector<qsd::Expectation>>> density_matrix(2);
    int kmax = density_matrix_expec[0][0].size();
    int dim = 2;
//    std::cout << kmax<< std::endl;

    if (qsd_result.size() > 0) {

        int ntraj = qsd_result.size();
//        std::cout<<"trajectory: " <<ntraj<<std::endl;

        for (int i = 0; i < dim; i++) {
            density_matrix[i].resize(dim);
            for (int j = 0; j < dim; ++j) {
                density_matrix[i][j].resize(kmax,qsd::Expectation(0.0, 0.0, 0.0, 0.0));
            }
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < kmax; ++k) {

                    std::complex<double> expec_0_mean = std::accumulate(density_matrix_expec[i][j][k].begin(), density_matrix_expec[i][j][k].end(), std::complex<double>(0.0, 0.0))/double(ntraj);
                    std::complex<double> expec_0_var = std::complex<double>(0.0, 0.0);

//                    if (ntraj > 1) {
//                        expec_0_var = variance(density_matrix_expec[i][j][k]);
//                    }
                    density_matrix[i][j][k] = qsd::Expectation(expec_0_mean.real(), expec_0_mean.imag(), expec_0_var.real(), expec_0_var.imag());
                }
            }
        }

        // Print the density_matrix
        /*
            std::cout << "Density Matrix \n";
        std::cout << density_matrix[0][0][99].mean << density_matrix[0][1][99].mean << density_matrix[1][0][99].mean << density_matrix[1][1][99].mean << std::endl;
        std::cout << density_matrix[0][2][99].mean << density_matrix[0][3][99].mean << density_matrix[1][2][99].mean << density_matrix[1][3][99].mean << std::endl;
        std::cout << density_matrix[2][0][99].mean << density_matrix[2][1][99].mean << density_matrix[3][0][99].mean << density_matrix[3][1][99].mean << std::endl;
        std::cout << density_matrix[2][2][99].mean << density_matrix[2][3][99].mean << density_matrix[3][2][99].mean << density_matrix[3][3][99].mean << std::endl;
        */
    }

    return density_matrix; //i & j: expectation value for measurement operators on target ions. k: time stamp evolution.
}

qsd::Expectation ion_trap::processor_even_population(std::vector<qsd::TrajectoryResult> qsd_result, bool save_per_traj, std::ofstream* log_file) {

    qsd::Expectation even_pop(0.0, 0.0, 0.0, 0.0);

    if (qsd_result.size() > 0) {

        int ntraj = qsd_result.size();

        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> density_matrix_expec = processor_concatenated_density_matrix(qsd_result);
        int kmax = density_matrix_expec[0][0].size();
        // save_per_traj
        if (save_per_traj) {
            for (int i = 0; i < ntraj; i++) {

                double _even_pop_00 = density_matrix_expec[0][0][kmax][i].real();
                double _even_pop_11 = density_matrix_expec[3][3][kmax][i].real();
                double _even_pop = _even_pop_00 + _even_pop_11;

                (*log_file) << std::setprecision(8) << _even_pop << std::endl;
            }
        }

        // Compute the even and odd populations
        std::complex<double> even_pop_00_mean = std::accumulate(density_matrix_expec[0][0][kmax].begin(), density_matrix_expec[0][0][kmax].end(), std::complex<double>(0.0, 0.0)) / double(ntraj);
        std::complex<double> even_pop_00_var = std::complex<double>(0.0, 0.0);

        std::complex<double> even_pop_11_mean = std::accumulate(density_matrix_expec[3][3][kmax].begin(), density_matrix_expec[3][3][kmax].end(), std::complex<double>(0.0, 0.0)) / double(ntraj);
        std::complex<double> even_pop_11_var = std::complex<double>(0.0, 0.0);

        if (ntraj > 1) {
            even_pop_00_var = variance(density_matrix_expec[0][0][kmax]);
            even_pop_11_var = variance(density_matrix_expec[3][3][kmax]);
        }

        even_pop = qsd::Expectation(even_pop_00_mean.real() + even_pop_11_mean.real(), 0.0, even_pop_00_var.real() + even_pop_11_var.real(), 0.0);
    }

    return even_pop;
}

qsd::TrajectoryResult ion_trap::processor_average_trajectory(std::vector<qsd::TrajectoryResult> qsd_result, std::vector<double> weights, bool save_per_traj, std::ofstream* log_file)
{
	qsd::TrajectoryResult avg_traj_result;
	if (qsd_result.size() > 0) {

		int ntraj = qsd_result.size();
		int numsteps = qsd_result[0].t.size();

		avg_traj_result.error = false;
		avg_traj_result.timeout = false;

		avg_traj_result.t = qsd_result[0].t;

		// compute the avg_traj_result.observables
		std::map<std::string, std::vector<qsd::Expectation>> observables;
		for (auto pair : qsd_result[0].observables) {
			observables[pair.first] = std::vector<qsd::Expectation>();
		}

		for (int n = 0; n < numsteps; n++) {

			std::map<std::string, std::vector<std::complex<double>>> expec_list_mean;
			std::map<std::string, std::vector<std::complex<double>>> expec_list_var;

			// // for (auto &traj_result : qsd_result) {
			for (int i = 0; i < qsd_result.size(); i++) {

				// // for (auto pair : traj_result.observables) {
				for (auto pair : qsd_result[i].observables) {

					if (weights.size() > 0) {
						expec_list_mean[pair.first].push_back(weights[i] * pair.second[n].mean);
						expec_list_var[pair.first].push_back(weights[i] *  pair.second[n].var); // is this correct?
					} else {
						expec_list_mean[pair.first].push_back(pair.second[n].mean);
						expec_list_var[pair.first].push_back(pair.second[n].var);
					}
				}
			}

			// compute the mean and var for expec_list
			for (auto pair : expec_list_mean) {

				std::complex<double> mean = std::accumulate(pair.second.begin(), pair.second.end(), std::complex<double>{0.0, 0.0});
				// // std::complex<double> var = std::accumulate(expec_list_var[pair.first].begin(), expec_list_var[pair.first].end(), std::complex<double>{0.0, 0.0});
				std::complex<double> var = experiment::variance(pair.second); // todo: how do we compute variance of weighted variables?

				if (weights.size() == 0) {
					mean /= double(ntraj);
					// // var /= double(ntraj); // how do we compute variance of weighted variables?
				}

				// // if (ntraj > 1) {
				// // expec_list_var = experiment::variance(pair.second);
				// // }

				// store the results in the "observables" variable
				observables[pair.first].push_back(qsd::Expectation(mean.real(), mean.imag(), var.real(), var.imag()));
			}
		}

		avg_traj_result.observables = observables;

	} else {
		std::cout << "qsd_result.size() == 0" << std::endl;
	}
	return avg_traj_result;
}

double ion_trap::processor(qsd::State qsd_state) {

    std::vector<double> results;

	for (auto pair : outlist) { // compute expectation values for output...

		std::string key = pair.first;
		qsd::Operator& observable = pair.second;

		qsd::State tmp_qsd_state = qsd_state;
		tmp_qsd_state *= observable;
		std::complex<double> result = qsd_state * tmp_qsd_state;
		results.push_back(result.real());
	}

    switch (processor_type) {

        case 0: { // Fidelity
            return std::max(results[0], results[1]);
        }

        case 1: { // Infidelity
            return (1.0 - std::max(results[0], results[1]));
        }
        default:
            std::cout<<"Unknown processor type"<<std::endl;
            return 0.0;
    }
}

ion_trap::ThermalState ion_trap::prepare_thermal_state() {

	ThermalState thermal_state;

	double cutoff_probability = state_cfg.cutoff_probability;

	std::vector<int> max_n;
	for (auto n_bar : state_cfg.n_bar) {
		max_n.push_back((int) (log(cutoff_probability * (n_bar + 1.0)) / log(n_bar / (n_bar + 1))));
	}
	int global_max_n = *std::max_element(max_n.begin(), max_n.end());
	int N_ions = state_cfg.motion.size();

	//// Experimental (permutation generator) ////
	std::vector<std::vector<int>> permutations;

	int length = N_ions;
	std::vector<int> list;
	for (int i = 0; i <= global_max_n; i++) list.push_back(i);

	permutations_with_repetition(list, std::vector<int>(), length, permutations);  //Note: this function works on all cases and not just the case above

	/*for (int k = 0; k < permutations.size(); k++) {
		for (auto item: permutations[k]) {
			std::cout << item;
		}
		cout << std::endl;
	}*/

	//// Experimental (filter permutations) ////
	std::vector<std::vector<double>> occupation_prob; // // (global_max_n + 1);

	for (int i = 0; i < N_ions; i++) {

		std::vector<double> _occupation_prob;
		// // for (int n = 0; n <= max_n[i]; n++) {
		for (int n = 0; n <= global_max_n; n++) {
			_occupation_prob.push_back(std::pow(state_cfg.n_bar[i], n) / std::pow(state_cfg.n_bar[i] + 1.0, n + 1.0));
		}
		occupation_prob.push_back(_occupation_prob);
	}

	for (auto permutation : permutations) {

		double _prob = occupation_prob[0][permutation[0]];
		for (int i = 1; i < permutation.size(); i++) {
			_prob *= occupation_prob[i][permutation[i]];
		}

		if (_prob >= cutoff_probability) {
			thermal_state.states.push_back(std::make_pair(_prob, permutation));
		}
	}

	//// Experimental (Normalize the thermal state) ////
	double total_prob = 0.0;
	for (auto permutation : thermal_state.states) {
		total_prob += permutation.first;
	}
	double prob_scale_factor = 1.0 / total_prob;

	for (auto &permutation : thermal_state.states) {
		permutation.first *= prob_scale_factor;
	}

	thermal_state.cutoff_probability = cutoff_probability;
	thermal_state.n_bar = state_cfg.n_bar;
	thermal_state.max_n = max_n;
	thermal_state.n_states = permutations.size();
	thermal_state.n_states_cutoff = thermal_state.states.size();

	return thermal_state;
}
}