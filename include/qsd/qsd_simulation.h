//
// Created by Shakib Vedaie on 2019-02-28.
//

#ifndef IONTRAP_QSD_SIMULATION_H
#define IONTRAP_QSD_SIMULATION_H

#include "ACG.h"
#include "CmplxRan.h"
#include "Traject.h"

#include <boost/foreach.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include <thread>
#include <random>
#include <numeric>
#include <cstdlib>

/*
 * QSD simulation.
 * */

template <class experiment>
class qsd_simulation {
public:
    qsd_simulation();
    qsd_simulation(double _timeout, int _n_processors, int _ntraj, int _accuracy_mode, double _accuracy_max_acc,
            int _accuracy_max_gen, unsigned int _seed, bool _save_on_disk, bool _save_per_traj);
    qsd_simulation(const std::string& config_path);

    void initialize(double _timeout, int _n_processors, int _ntraj, int _accuracy_mode, double _accuracy_max_acc,
                    int _accuracy_max_gen, unsigned int _seed, bool _save_on_disk, bool _save_per_traj);
    void initialize(const std::string& config_path);

    void set_accuracy(double _accuracy);

    int get_accuracy_mode();
    double get_accuracy_max_acc();
    int get_accuracy_max_gen();
    int get_ntraj();

	void set_seed(unsigned int _seed);

    qsd::State concatenate_qsd_states(experiment &expt, std::vector<qsd::TrajectoryResult> qsd_results);

    void parse_config(const std::string &config_path);
	std::vector<qsd::TrajectoryResult> run(experiment &expt, std::function<void(qsd::TrajectoryResult&, int, double)> callback=std::function<void(qsd::TrajectoryResult&, int, double)>());

private:
	std::string log_dir;

    double dt;      // basic time step
    int numdts;     // time interval between outputs = numdts * dt
    int numsteps;   // total integration time = numsteps * numdts * dt

    double timeout; // Simulation timeout in seconds
    int n_processors;

    int ntraj;
	unsigned int seed;

    int accuracy_mode;
    double accuracy;
    double accuracy_max_acc;
    int accuracy_max_gen;

	bool log_progress;
	bool log_errors;

    bool save_on_disk;
    bool save_per_traj;

	bool clone_expt;
    int unravelling_mode;
};

template <class experiment>
qsd_simulation<experiment>::qsd_simulation() {}

template <class experiment>
qsd_simulation<experiment>::qsd_simulation(double _timeout, int _n_processors, int _ntraj, int _accuracy_mode,
        double _accuracy_max_acc, int _accuracy_max_gen, unsigned int _seed, bool _save_on_disk, bool _save_per_traj)
{
    timeout = _timeout;
    n_processors = _n_processors;

    ntraj = _ntraj;
    seed = _seed;

    accuracy_mode = _accuracy_mode;
    accuracy_max_acc = _accuracy_max_acc;
    accuracy_max_gen = _accuracy_max_gen;

    if (accuracy_mode == 0)
        accuracy = accuracy_max_acc;

    save_on_disk = _save_on_disk;
    save_per_traj = _save_per_traj;
}

template <class experiment>
qsd_simulation<experiment>::qsd_simulation(const std::string& config_path)
{
    parse_config(config_path);
}

template <class experiment>
void qsd_simulation<experiment>::initialize(const std::string& config_path)
{
    parse_config(config_path);
}

template <class experiment>
void qsd_simulation<experiment>::initialize(double _timeout, int _n_processors, int _ntraj, int _accuracy_mode,
        double _accuracy_max_acc, int _accuracy_max_gen, unsigned int _seed, bool _save_on_disk, bool _save_per_traj)
{
    timeout = _timeout;
    n_processors = _n_processors;

    ntraj = _ntraj;
    seed = _seed;

    accuracy_mode = _accuracy_mode;
    accuracy_max_acc = _accuracy_max_acc;
    accuracy_max_gen = _accuracy_max_gen;

    if (accuracy_mode == 0)
        accuracy = accuracy_max_acc;

    save_on_disk = _save_on_disk;
    save_per_traj = _save_per_traj;
}

template <class experiment>
void qsd_simulation<experiment>::set_accuracy(double _accuracy)
{
    accuracy = _accuracy;
}

template <class experiment>
int qsd_simulation<experiment>::get_ntraj()
{
    return ntraj;
}

template <class experiment>
void qsd_simulation<experiment>::set_seed(unsigned int _seed)
{
	seed = _seed;
}

template <class experiment>
int qsd_simulation<experiment>::get_accuracy_mode() { return accuracy_mode; };

template <class experiment>
double qsd_simulation<experiment>::get_accuracy_max_acc() { return accuracy_max_acc; };

template <class experiment>
int qsd_simulation<experiment>::get_accuracy_max_gen() { return accuracy_max_gen; };

template <class experiment>
void qsd_simulation<experiment>::parse_config(const std::string &config_path)
{
    // Create empty property tree object
	boost::property_tree::ptree experiment_cfg;
    // // boost::property_tree::ptree qsd_simulation_cfg;

    // Parse the INFO into the property tree.
	boost::property_tree::read_info(config_path + "experiment.cfg", experiment_cfg);
    // // boost::property_tree::read_info(config_path + "qsd_simulation.cfg", qsd_simulation_cfg);

	// Read the log_dir
	log_dir = experiment_cfg.get<std::string>("experiment.log_dir");

    // Read the simulation configurations
    timeout  = experiment_cfg.get<double>("experiment.timeout");
    n_processors = experiment_cfg.get<int>("experiment.n_processors");

    ntraj = experiment_cfg.get<int>("experiment.ntraj");
    seed = experiment_cfg.get<unsigned int>("experiment.seed");

    accuracy_mode = experiment_cfg.get<int>("experiment.accuracy.mode");
    accuracy_max_acc = experiment_cfg.get<double>("experiment.accuracy.max_acc");
    accuracy_max_gen = experiment_cfg.get<int>("experiment.accuracy.max_gen");

    if (accuracy_mode == 0)
        accuracy = accuracy_max_acc;

	log_progress = experiment_cfg.get<bool>("experiment.log.progress");
	log_errors = experiment_cfg.get<bool>("experiment.log.errors");

    save_on_disk = experiment_cfg.get<bool>("experiment.log.save_on_disk");
    save_per_traj = experiment_cfg.get<bool>("experiment.log.save_per_traj");

	clone_expt = experiment_cfg.get<bool>("experiment.clone_expt");

    unravelling_mode = experiment_cfg.get<bool>("experiment.unravelling_mode");
}

template <class experiment>
std::vector<qsd::TrajectoryResult> qsd_simulation<experiment>::run(experiment &expt, std::function<void(qsd::TrajectoryResult&, int, double)> callback)
{
    /*
     *  Check if the input state is mixed. For pure states a single run of the experiment is enough
     *
     */

    // The random number generator
    ACG gen(seed, 55);

    // Stepsize and integration time
    double dt = expt.dt;            // basic time step
    int numdts = expt.numdts;       // time interval between outputs = numdts * dt
    int numsteps = expt.numsteps;   // total integration time = numsteps * numdts * dt

    // // std::cout << n_processors << " concurrent threads are available.\n";

    // // std::vector<double> results(n_processors);
    std::vector<std::vector<qsd::TrajectoryResult>> qsd_results_workers;
    auto worker = [&] (int worker_id, int n_traj) {

		// Prepare the trajectory config file
        qsd::TrajectoryConfig cfg;

		cfg.log_dir = log_dir;

		cfg.log_progress = log_progress;
        cfg.save_on_disk = save_on_disk;
        cfg.save_per_traj = save_per_traj;

        cfg.dynamic_degrees = expt.state_cfg.dynamic_degrees;

        cfg.dynamic_cutoff = expt.state_cfg.dynamic_cutoff;
        cfg.cutoff_epsilon = expt.state_cfg.cutoff_epsilon;
        cfg.cutoff_pad_size = expt.state_cfg.cutoff_pad_size;

        cfg.phonon_moving_basis = expt.state_cfg.phonon_moving_basis;
        cfg.shift_accuracy = expt.state_cfg.shift_accuracy;

		// The random number generator for "qsd::Trajectory"
		ACG _gen(gen.asLong(), 55);
	  	std::mt19937 rnd_gen(gen.asLong());

		ComplexNormal rand1(&_gen);

	  	ComplexNormal *rand1_pt = &rand1;
	  	std::mt19937 *rnd_gen_pt = &rnd_gen;

		if (expt.L.size() == 0) {
			rand1_pt = nullptr;
			rnd_gen_pt = nullptr;
		}

		// Clone the expt
		if (clone_expt) {
			experiment expt_clone(expt);

			// Initialize the copy of experiment "expt" with a random seed
			unsigned int expt_clone_seed = gen.asLong();
			std::cout << "seed: " << expt_clone_seed << std::endl; // ......
			expt_clone.initialize_expt(expt_clone_seed);

            // deterministic part: adaptive stepsize 4th/5th order Runge Kutta
			// stochastic part: fixed stepsize Euler

            switch (unravelling_mode) {
                case 0:{
                    // AdaptiveStep (Quantum State Diffusion)

        			qsd::AdaptiveStep theStepper(expt_clone.psi0, expt_clone.H, expt_clone.L, accuracy);
                    qsd::Trajectory theTraject(expt_clone.psi0, dt, theStepper, rand1_pt, rnd_gen_pt);
                    std::vector<qsd::TrajectoryResult> qsd_results_worker = theTraject.sumExp_vector(expt_clone.outlist, expt_clone.flist, numdts, numsteps, timeout, n_traj, worker_id, cfg, callback);
                    qsd_results_workers[worker_id] = qsd_results_worker;
                    break;
                }
                case 1: {
                    // AdaptiveJump
                    qsd::AdaptiveJump theStepper(expt_clone.psi0, expt_clone.H, expt_clone.L, accuracy);
                    qsd::Trajectory theTraject(expt_clone.psi0, dt, theStepper, rand1_pt, rnd_gen_pt);
                    std::vector<qsd::TrajectoryResult> qsd_results_worker = theTraject.sumExp_vector(expt_clone.outlist, expt_clone.flist, numdts, numsteps, timeout, n_traj, worker_id, cfg, callback);
                    qsd_results_workers[worker_id] = qsd_results_worker;
                    break;
                }
                default: {
                    std::cerr<<"Unknown unravelling mode!"<<std::endl;
                    break;
                }
            }

		} else {

			// deterministic part: adaptive stepsize 4th/5th order Runge Kutta
			// stochastic part: fixed stepsize Euler

			// AdaptiveStep (Quantum State Diffusion)
			qsd::AdaptiveStep theStepper(expt.psi0, expt.H, expt.L, accuracy);


            // AdaptiveJump (Quantum jump)
//            qsd::AdaptiveJump theStepper(expt.psi0, expt.H, expt.L, accuracy);

			qsd::Trajectory theTraject(expt.psi0, dt, theStepper, rand1_pt, rnd_gen_pt);

			std::vector<qsd::TrajectoryResult> qsd_results_worker = theTraject.sumExp_vector(expt.outlist, expt.flist, numdts, numsteps, timeout, n_traj, worker_id, cfg, callback);
			qsd_results_workers[worker_id] = qsd_results_worker;
		}

		std::cout << "qsd worker " << worker_id << " done!" << std::endl;
    };

	if (n_processors > 1) {

		//// Experimental ////
		qsd_results_workers.resize(ntraj);

		// Launch the pool with "n_processors" threads.
		boost::asio::thread_pool pool(n_processors);

		// Submit a function to the pool.
		for (int worker_id = 0; worker_id < ntraj; worker_id++) {
			boost::asio::post(pool, std::bind(worker, worker_id, 1));
		}

		// Wait for all tasks in the pool to complete.
		pool.join();
		//////////////////////

	} else {
		qsd_results_workers.resize(1);
		worker(0, ntraj);
	}

    // // return std::accumulate(results.begin(), results.end(), 0.0) / results.size();

    // filter all the QSD results
    std::vector<qsd::TrajectoryResult> result;

	int corrupted_cnt = 0;
	for (auto worker_results : qsd_results_workers) {
		for (auto traj_result : worker_results) {
			if (!traj_result.error && !traj_result.timeout) {
				result.push_back(traj_result);
			} else {
				corrupted_cnt++;
			}
		}
	}

	if (log_errors) {
		std::cout << corrupted_cnt << " of the trajectories were corrupted!" << std::endl;
	}

    // result.t = qsd_results[0].t;

	/*
    result.data.resize(qsd_results[0].data.size());
    for (int i = 0; i < qsd_results[0].data.size(); i++) {
        for (int j = 0; j < qsd_results[0].data[i].size(); j++) {
            result.data[i].push_back({0.0, 0.0, 0.0, 0.0}); // Re(<...>), Im(<...>), Re(<...^2>)
        }
    }
	*/

	/*
    for (qsd::Result &qsd_res : qsd_results) {

        if (qsd_res.error)
            result.error = true;

        if (qsd_res.timeout)
            result.timeout = true;

        // Concatenate the data
        for (auto &qsd_res_data : qsd_res.data) {
			result.data.push_back(qsd_res_data);
		}


//		for (int i = 0; i < qsd_res.data.size(); i++) {
//            for (int j = 0; j < qsd_res.data[i].size(); j++) {
//                result.data[i][j] = {result.data[i][j][0] + (qsd_res.state.size() * qsd_res.data[i][j][0]) / ntraj,
//                                     result.data[i][j][1] + (qsd_res.state.size() * qsd_res.data[i][j][1]) / ntraj,
//                                     0.0, 0.0};
//            }
//        }


        // Concatenate the qsd::States
        for (qsd::State &qsd_res_state : qsd_res.state) {
            result.state.push_back(qsd_res_state);
        }
    }
	*/

    return result;
}

template <class experiment>
qsd::State qsd_simulation<experiment>::concatenate_qsd_states(experiment &expt, std::vector<qsd::TrajectoryResult> qsd_results) {

    // Make the cutoffs of the state match the physical size in memory and
    // move the center of coordinates to a new position in shaving_phase space.
    for (int i = 0; i < qsd_results.size(); i++) {
        // for (int j = 0; j < qsd_results[i].state.size(); j++) {
		qsd_results[i].state.fullSize();
		for (int k : expt.state_cfg.dynamic_degrees) {
			qsd_results[i].state.moveToCoords({0.0, 0.0}, k, expt.state_cfg.shift_accuracy);
		}
        // }
    }

    // Put the results together
    qsd::State psi_f;
    for (int i = 0; i < qsd_results.size(); i++) {

		// for (int j = 0; j < qsd_results[i].state.size(); j++) {

		if (i == 0) {
			psi_f = qsd_results[i].state;
		} else {
			psi_f = psi_f + qsd_results[i].state;
		}

        // }
    }
    psi_f.normalize();

    return psi_f;
}

#endif //IONTRAP_QSD_SIMULATION_H
