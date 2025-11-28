/*
 * Copyright (c) 2019, Seyed Shakib Vedaie & Eduardo J. Paez
 */

#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include <boost/asio/post.hpp>
#include <boost/filesystem.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include <pybind11/stl.h>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

#include <qsd/Traject.h>
#include <qsd/qsd_simulation.h>

#include "ion_trap.h"
//#include "fidelity/fidelity.h"

boost::property_tree::ptree mode_cfg;
//boost::property_tree::ptree trap_cfg;
boost::property_tree::ptree solutions_cfg;
boost::property_tree::ptree experiment_cfg;
boost::property_tree::ptree interaction_cfg;

pybind11::object IonTrapEnv;
pybind11::object RLAgent;

pybind11::object ion_trap_env;
pybind11::object rl_agent;

int last_idx_AM = 0;
std::vector<double> rabi_rate;
std::vector<std::vector<double>> state_debug;


void reinforcement_learning_callback_debug(experiment::ion_trap *IonTrap, qsd::TrajectoryResult &traj_result, int n, double t)
{
	std::vector<double> state;

	/*
    for (int i = 2; i <= 7; i++) // {x_1, p_1, x_2, p_2, ...}
		state.push_back(OutputExp[0].ExpR[n][i]);
	*/

	state.push_back(traj_result.observables["joint_x_0_0_x_0"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_1_x_1"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_0_x_2"][n].mean.real());

	state.push_back(traj_result.observables["joint_x_0_0_p_0"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_1_p_1"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_0_p_2"][n].mean.real());

    int idx_AM = 0;
	double _idx = t * (IonTrap->pulse_cfg.steps_am / IonTrap->t_gate);
	if (std::fabs(_idx - int(_idx + 1)) < 1e-6) {
		idx_AM = int(_idx) + 1;
	} else {
		idx_AM = int(_idx);
	}

	if ((std::fabs(t) < 1e-6) || (idx_AM > last_idx_AM)) {
        state_debug.push_back(state);
//        std::cout << "idx_AM: " << idx_AM << std::endl;
//        for (int i = 0; i < state.size(); i++) {
//            std::cout << "state[" << i << "] = " << state[i] << std::endl;
//        }
		last_idx_AM = idx_AM;
	}
}

void reinforcement_learning_debug(const std::string& config_path, std::ofstream* log_file)
{
	std::random_device rd;
	bool auto_qsd_simulation_seed = mode_cfg.get<bool>("mode-6.auto_qsd_simulation_seed");

	experiment::ion_trap IonTrap(config_path);

	// Set the pulse scale to 1.0
	IonTrap.pulse_cfg.scale_am = 1.0;

	qsd_simulation<experiment::ion_trap> QSDSim(config_path);
	if (auto_qsd_simulation_seed) {

		unsigned int seed = rd();
		QSDSim.set_seed(seed);
		std::cout << "\nQSDSim.set_seed(" << seed << ")" << std::endl;
	}

	bool use_external_delta = solutions_cfg.get<bool>("solutions.use_external_delta");
	if (use_external_delta) {

		std::string select = solutions_cfg.get<std::string>("solutions.external.select");
		boost::property_tree::ptree external_delta_cfg = solutions_cfg.get_child("solutions.external");
		IonTrap.set_external_delta(external_delta_cfg, select);

	} else {

		IonTrap.set_delta(solutions_cfg.get<int>("solutions.internal.delta_idx"));
	}

	bool use_external_pulse = solutions_cfg.get<bool>("solutions.use_external_pulse");
	if (use_external_pulse) {

		std::string select = solutions_cfg.get<std::string>("solutions.external.select");
		boost::property_tree::ptree external_pulse_cfg = solutions_cfg.get_child("solutions.external");
		IonTrap.set_pulse_external(external_pulse_cfg, select);

	} else {

		IonTrap.pulse_cfg.scale_am = solutions_cfg.get<double>("solutions.internal.pulse_scale");
		IonTrap.set_pulse(solutions_cfg.get<int>("solutions.internal.pulse_idx"));
	}

	std::vector<qsd::TrajectoryResult> qsd_result;
	qsd::Expectation result(0.0, 0.0, 0.0, 0.0);

	experiment::StateConfig state_cfg = IonTrap.state_cfg;
	if (state_cfg.motion_state_type=="PURE") {
		auto begin = std::chrono::high_resolution_clock::now();

		if (IonTrap.state_cfg.use_x_basis) {
			// Not implemented
			// ...
		}
		else {
			qsd_result = QSDSim.run(IonTrap, std::bind(reinforcement_learning_callback_debug, &IonTrap, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
			result = IonTrap.processor_fidelity(qsd_result);
		}

		(*log_file) << "Result: " << std::setprecision(8) << result.mean.real() << std::endl;

		auto end = std::chrono::high_resolution_clock::now();
		(*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s)\n";
	}
	else if (state_cfg.motion_state_type=="THERMAL") {
		// Not implemented
		// ...
	}

	// Prepare the state
	std::vector<double> state;
	if (qsd_result.size() > 0) {

		/*
		for (int i = 2; i <= 7; i++) // {x_1, p_1, x_2, p_2, ...}
			state.push_back(qsd_result.data[i].back()[0]);
		 */

		state.push_back(qsd_result[0].observables["joint_x_0_0_x_0"].back().mean.real());
		state.push_back(qsd_result[0].observables["joint_x_0_1_x_1"].back().mean.real());
		state.push_back(qsd_result[0].observables["joint_x_0_0_x_2"].back().mean.real());

		state.push_back(qsd_result[0].observables["joint_x_0_0_p_0"].back().mean.real());
		state.push_back(qsd_result[0].observables["joint_x_0_1_p_1"].back().mean.real());
		state.push_back(qsd_result[0].observables["joint_x_0_0_p_2"].back().mean.real());
	}
    state_debug.push_back(state);

	std::cout << "idx_AM: last" << std::endl;

    for (int i = 0; i < state.size(); i++) {
        for (int j = 0; j <= (last_idx_AM + 1); j++) {
            std::cout << state_debug[j][i] << ", ";
//            std::cout << "state[" << i << "] = " << state[i] << std::endl;
        }
        std::cout << "Operator: " << i << std::endl;
    }

//	for (int i = 0; i < state.size(); i++) {
//		std::cout << "state[" << i << "] = " << state[i] << std::endl;
//	}
}

void reinforcement_learning_callback(experiment::ion_trap *IonTrap, qsd::TrajectoryResult &traj_result, int n, double t)
{
	std::vector<double> state;
	/*
	for (int i = 2; i <= 7; i++) // {x_1, p_1, x_2, p_2, ...}
		state.push_back(OutputExp[0].ExpR[n][i]);
	 */

	state.push_back(traj_result.observables["joint_x_0_0_x_0"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_1_x_1"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_0_x_2"][n].mean.real());

	state.push_back(traj_result.observables["joint_x_0_0_p_0"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_1_p_1"][n].mean.real());
	state.push_back(traj_result.observables["joint_x_0_0_p_2"][n].mean.real());

	int idx_AM = 0;
	double _idx = t * (IonTrap->pulse_cfg.steps_am / IonTrap->t_gate);
	if (std::fabs(_idx - int(_idx + 1)) < 1e-6) {
		idx_AM = int(_idx) + 1;
	} else {
		idx_AM = int(_idx);
	}

	if (std::fabs(t) < 1e-6) {
//        std::cout << std::setprecision(8) << t << std::endl;

		std::vector<double> action = rl_agent.attr("take_action")().cast<std::vector<double>>();

		rabi_rate.push_back(action[0]);

		IonTrap->set_pulse(rabi_rate);
	}

	if (idx_AM > last_idx_AM) {
//        std::cout << std::setprecision(8) << t << std::endl;

        // Prepare the state
		state.push_back((double) idx_AM);
//        std::cout<<state[0]<<" -- "<<state[1]<<" -- "<<state[2]<<std::endl;

		pybind11::object step_act = rl_agent.attr("observe_and_learn")(state);
		std::vector<double> action = rl_agent.attr("take_action")().cast<std::vector<double>>();

		rabi_rate.push_back(action[0]);
        IonTrap->set_pulse(rabi_rate);

		last_idx_AM = idx_AM;
	}
}

void reinforcement_learning(const std::string& config_path, std::ofstream* log_file)
{
	int mode = mode_cfg.get<int>("mode-6.mode");

	std::string py_dir = mode_cfg.get<std::string>("mode-6.py_dir");
	std::string save_dir = mode_cfg.get<std::string>("mode-6.save_dir");

	double max_xp = mode_cfg.get<double>("mode-6.max_xp");
	int max_step = interaction_cfg.get<int>("interaction.pulse.AM.steps"); // experiment_cfg.get<int>("experiment.pulse.FM.steps");
	double max_rabi = mode_cfg.get<double>("mode-6.max_rabi");

	int max_epi = mode_cfg.get<int>("mode-6.max_epi");
	double termination_condition = mode_cfg.get<double>("mode-6.termination_condition");

	int avg_window = mode_cfg.get<int>("mode-6.avg_window");

	bool auto_qsd_simulation_seed = mode_cfg.get<bool>("mode-6.auto_qsd_simulation_seed");

	pybind11::scoped_interpreter guard{}; // start the interpreter and keep it alive

	pybind11::module sys = pybind11::module::import("sys");
	pybind11::module os = pybind11::module::import("os");

	sys.attr("path").attr("append")(os.attr("path").attr("dirname")(os.attr("getcwd")()).cast<std::string>()+py_dir);

	IonTrapEnv = pybind11::module::import("ion_trap_env").attr("IonTrapEnv");
	ion_trap_env = IonTrapEnv(max_xp, max_step, max_rabi);

	bool training;
	if (mode == 0) { // training
		training = true;
	} else { // test or debug
		training = false;
	}

	std::string rl_agent_select = mode_cfg.get<std::string>("mode-6.rl_agent.select");

	if (rl_agent_select == "ddpg") { // ddpg

		double critic_lr = mode_cfg.get<double>("mode-6.rl_agent.ddpg.critic_lr");
		double actor_lr = mode_cfg.get<double>("mode-6.rl_agent.ddpg.actor_lr");
		double gamma = mode_cfg.get<double>("mode-6.rl_agent.ddpg.gamma");
		double tau = mode_cfg.get<double>("mode-6.rl_agent.ddpg.tau");
		double ou_noise_sd = mode_cfg.get<double>("mode-6.rl_agent.ddpg.ou_noise_sd");
		int buffer_capacity = mode_cfg.get<int>("mode-6.rl_agent.ddpg.buffer_capacity");
		int batch_size = mode_cfg.get<int>("mode-6.rl_agent.ddpg.batch_size");

		RLAgent = pybind11::module::import("ddpg").attr("DDPG");
		rl_agent = RLAgent(ion_trap_env, training, critic_lr, actor_lr, gamma, tau, ou_noise_sd, buffer_capacity, batch_size);

	} else if (rl_agent_select == "ppo") { // ppo

		double critic_lr = mode_cfg.get<double>("mode-6.rl_agent.ppo.critic_lr");
		double actor_lr = mode_cfg.get<double>("mode-6.rl_agent.ppo.actor_lr");
		double epsilon = mode_cfg.get<double>("mode-6.rl_agent.ppo.epsilon");
		double entropy_weight = mode_cfg.get<double>("mode-6.rl_agent.ppo.entropy_weight");
		double gamma = mode_cfg.get<double>("mode-6.rl_agent.ppo.gamma");
		double tau = mode_cfg.get<double>("mode-6.rl_agent.ppo.tau");
		int epoch = mode_cfg.get<int>("mode-6.rl_agent.ppo.epoch");
		int batch_size = mode_cfg.get<int>("mode-6.rl_agent.ppo.batch_size");
		int rollout_len = mode_cfg.get<int>("mode-6.rl_agent.ppo.rollout_len");

		RLAgent = pybind11::module::import("ppo").attr("PPO");
		rl_agent = RLAgent(ion_trap_env, training, critic_lr, actor_lr, epsilon, entropy_weight, gamma, tau, epoch, batch_size, rollout_len);
	}

	if (mode == 1) { // test
		rl_agent.attr("load_weights")(save_dir);
	}

	// Open the log files
	std::ofstream rewards_file;
	std::ofstream fidelities_file;

	rewards_file.open(save_dir + "rewards.txt");
	fidelities_file.open(save_dir + "fidelities.txt");
	//

	std::random_device rd;
	std::vector<double> ep_reward_list; // reward history of each episode
	std::vector<double> ep_fidelity_list; // fidelity history of each episode
	bool termination_condition_status = false;
	for (int episode = 0; episode < max_epi; episode++) {

		last_idx_AM = 0;
		std::vector<double>().swap(rabi_rate);
		rl_agent.attr("reset")();

		experiment::ion_trap IonTrap(config_path);

		// Set the pulse scale to 1.0
		IonTrap.pulse_cfg.scale_am = 1.0;

		qsd_simulation<experiment::ion_trap> QSDSim(config_path);
		if (auto_qsd_simulation_seed) {

			unsigned int seed = rd();
			QSDSim.set_seed(seed);
			std::cout << "\nQSDSim.set_seed(" << seed << ")" << std::endl;
		}

		bool use_external_delta = solutions_cfg.get<bool>("solutions.use_external_delta");
		if (use_external_delta) {

			std::string select = solutions_cfg.get<std::string>("solutions.external.select");
			boost::property_tree::ptree external_delta_cfg = solutions_cfg.get_child("solutions.external");
			IonTrap.set_external_delta(external_delta_cfg, select);

		} else {
			IonTrap.set_delta(solutions_cfg.get<int>("solutions.internal.delta_idx"));
		}

//		bool use_external_pulse = trap_solution_cfg.get<bool>("solutions.use_external_pulse");
//		if (use_external_pulse) {
//
//			boost::property_tree::ptree external_pulse_cfg = trap_solution_cfg.get_child("solutions.external");
//			IonTrap.set_pulse_external(external_pulse_cfg);
//
//		} else {
//
//			IonTrap.pulse_cfg.scale_am = trap_solution_cfg.get<double>("solutions.internal.pulse_scale");
//			IonTrap.set_pulse(trap_solution_cfg.get<int>("solutions.internal.pulse_idx"));
//		}

		std::vector<qsd::TrajectoryResult> qsd_result;
		qsd::Expectation result(0.0, 0.0, 0.0, 0.0);

		experiment::StateConfig state_cfg = IonTrap.state_cfg;
		if (state_cfg.motion_state_type=="PURE") {
			auto begin = std::chrono::high_resolution_clock::now();

			if (IonTrap.state_cfg.use_x_basis) {
				// Not implemented
				// ...
			}
			else {
				qsd_result = QSDSim.run(IonTrap, std::bind(reinforcement_learning_callback, &IonTrap, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));

                result = IonTrap.processor_fidelity(qsd_result);
			}

			(*log_file) << "Result: " << std::setprecision(8) << result.mean.real() << std::endl;

			auto end = std::chrono::high_resolution_clock::now();
			(*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s)\n";
		}
		else if (state_cfg.motion_state_type=="THERMAL") {
			// Not implemented
			// ...
		}

		// Prepare the state
		std::vector<double> state;
		if (qsd_result.size() > 0) {
			/*
			for (int i = 2; i <= 7; i++) // {x_1, p_1, x_2, p_2, ...}
				state.push_back(qsd_result.data[i].back()[0]);
			 */

			state.push_back(qsd_result[0].observables["joint_x_0_0_x_0"].back().mean.real());
			state.push_back(qsd_result[0].observables["joint_x_0_1_x_1"].back().mean.real());
			state.push_back(qsd_result[0].observables["joint_x_0_0_x_2"].back().mean.real());

			state.push_back(qsd_result[0].observables["joint_x_0_0_p_0"].back().mean.real());
			state.push_back(qsd_result[0].observables["joint_x_0_1_p_1"].back().mean.real());
			state.push_back(qsd_result[0].observables["joint_x_0_0_p_2"].back().mean.real());
		}
		state.push_back((double) (last_idx_AM + 1));

		rl_agent.attr("observe_and_learn")(state, result.mean.real());

		ep_reward_list.push_back(rl_agent.attr("episodic_reward").cast<double>());
		ep_fidelity_list.push_back(result.mean.real());

		std::cout << "Episode * " << episode << " * Reward is ==> " << ep_reward_list.back() << " (fidelity = " << ep_fidelity_list.back() << ")" << std::endl;
		std::cout << "Pulse: ";
		for (auto step : IonTrap.pulse[0]) {
			std::cout << step << ", ";
		}
		std::cout << std::endl;

		if ((episode > 0) && ((episode % (avg_window - 1)) == 0)) {
			double avg_ep_reward = std::accumulate(ep_reward_list.begin() + (int(episode / (avg_window - 1)) - 1) * (avg_window - 1), ep_reward_list.end(), 0.0) / avg_window;
			double avg_ep_fidelity = std::accumulate(ep_fidelity_list.begin() + (int(episode / (avg_window - 1)) - 1) * (avg_window - 1), ep_fidelity_list.end(), 0.0) / avg_window;

			if (qsd_result.size() > 0) {
				rewards_file << avg_ep_reward << std::endl;
				fidelities_file << avg_ep_fidelity << std::endl;

				// // rewards_file << ep_reward_list.back() << std::endl;
				// // fidelities_file << ep_fidelity_list.back() << std::endl;
			}
		}

		if (training && (ep_fidelity_list.back() >= termination_condition)) {
			std::cout << "termination_condition" << std::endl;
			termination_condition_status = true;
			break;
		}
	}

	if (training) {
		rl_agent.attr("save_weights")(save_dir, termination_condition_status);
	}
}

int main(int argc, char *argv[])
{
    // Read the config file path
    std::string config_path;
    if (argc == 1) {
        std::cout << "Error: Config file path is missing.\n";
        std::cout << "IonTrap-RL config_file_path\n";
        return 0;
    } else {
        config_path = std::string(argv[1]);
    }

    // Read the simulation configurations
	boost::property_tree::read_info(config_path + "experiment.cfg", experiment_cfg);
	boost::property_tree::read_info(config_path + "interaction.cfg", interaction_cfg);

    // Open the log file
	std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");
	std::string log_filename = experiment_cfg.get<std::string>("experiment.log_filename");

	boost::filesystem::path dir(log_dir);
	if(!(boost::filesystem::exists(dir))){
		if (boost::filesystem::create_directory(dir))
			std::cout << "log_dir successfully created!" << std::endl;
	}

	std::ofstream log_file;
	log_file.open(log_dir + log_filename);

    int simulation_mode = experiment_cfg.get<int>("experiment.simulation_mode");

    switch (simulation_mode) {
        case 6: {
            /*
             * Reinforcement Learning
             * */

			boost::property_tree::read_info(config_path + "mode-6.cfg", mode_cfg);

			// Read the trap solution config
			boost::property_tree::read_info(config_path + "solutions.cfg", solutions_cfg);

			int reinforcement_learning_mode = mode_cfg.get<int>("mode-6.mode");
			if (reinforcement_learning_mode == 2) { // debug
				reinforcement_learning_debug(config_path, &log_file);
			} else {
				reinforcement_learning(config_path, &log_file);
			}

            break;
        }

        default: {
            std::cout << "Error: Unknown simulation mode.\n";
			break;
        }
    }

    return 0;
}