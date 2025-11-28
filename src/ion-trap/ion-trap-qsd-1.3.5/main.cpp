/*
 * Copyright (c) 2019, Seyed Shakib Vedaie & Eduardo J. Paez
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>

#include <boost/asio/post.hpp>
#include <boost/filesystem.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include <qsd/qsd_simulation.h>
#include <Eigen/LU>
//#include <Eigen/Dense>
//#include <Eigen/SparseCore>
#include <unsupported/Eigen/MatrixFunctions>


#include "ion_trap.h"
#include "utilities.h"
//#include "fidelity/fidelity.h"
using namespace Eigen;

boost::property_tree::ptree mode_cfg;
boost::property_tree::ptree trap_cfg;
boost::property_tree::ptree solutions_cfg;
boost::property_tree::ptree experiment_cfg;
boost::property_tree::ptree interaction_cfg;

typedef std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> time_point;

void print_time(time_point begin, time_point end, std::ofstream* log_file)
{
	int time_s = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
	int time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

	if (time_s > 0) {
		(*log_file) << " (took: " << time_s << " s)\n";
	} else {
		(*log_file) << " (took: " << time_ms << " ms)\n";
	}
}

void save_trajectory_to_disk(qsd::TrajectoryResult avg_traj_result) {

	std::map<std::string, FILE*> fp_traj;
	std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");

	for (auto pair : avg_traj_result.observables) {

		std::string filename = pair.first + "_avg_traj";
		fp_traj[filename] = fopen((log_dir + filename).c_str(), "w");

		fprintf(fp_traj[filename], "Average Trajectory\n");
		fflush(fp_traj[filename]);
	}

	for (int n = 0; n < avg_traj_result.t.size(); n++) {
		for (auto pair : avg_traj_result.observables) {

			std::string key = pair.first;
			std::vector<qsd::Expectation> &expec = pair.second;

			std::string filename = key + "_avg_traj";

			fprintf(fp_traj[filename],
					"%lG %lG %lG %lG %lG\n",
					avg_traj_result.t[n],
					expec[n].mean.real(),
					expec[n].mean.imag(),
					expec[n].var.real(),
					expec[n].var.imag());

			fflush(fp_traj[filename]);
		}
	}

	for (auto pair : fp_traj) fclose(pair.second);
}

void save_trajectory_to_disk_freq(qsd::TrajectoryResult avg_traj_result, int freq) {

    std::map<std::string, FILE*> fp_traj;
    std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");

    for (auto pair : avg_traj_result.observables) {

        std::string filename = pair.first + "_avg_traj_" + std::to_string(freq);
        fp_traj[filename] = fopen((log_dir + filename).c_str(), "w");

        fprintf(fp_traj[filename], "Average Trajectory\n");
        fflush(fp_traj[filename]);
    }

    for (int n = 0; n < avg_traj_result.t.size(); n++) {
        for (auto pair : avg_traj_result.observables) {

            std::string key = pair.first;
            std::vector<qsd::Expectation> &expec = pair.second;

            std::string filename = key + "_avg_traj_" + std::to_string(freq);

            fprintf(fp_traj[filename],
                    "%lG %lG %lG %lG %lG\n",
                    avg_traj_result.t[n],
                    expec[n].mean.real(),
                    expec[n].mean.imag(),
                    expec[n].var.real(),
                    expec[n].var.imag());

            fflush(fp_traj[filename]);
        }
    }

    for (auto pair : fp_traj) fclose(pair.second);
}
double fidelity_x_basis(experiment::ion_trap& IonTrap, qsd_simulation<experiment::ion_trap>& QSDSim)
{
    std::vector<std::vector<qsd::State>> states;
    std::vector<qsd::State> x_basis;

    // Defining the basis states
    qsd::State b_0(2, qsd::SPIN);
    qsd::State b_1(2, qsd::SPIN);

    qsd::SigmaPlus sigma_plus;
    b_1 *= sigma_plus;

    x_basis.push_back(b_0 + b_1);
    x_basis.push_back(b_0 - b_1);
    x_basis[0].normalize();
    x_basis[1].normalize();

    // Preparing psi0
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            states.push_back({x_basis[i], x_basis[j]});
        }
    }

    std::vector<qsd::TrajectoryResult> qsd_results;
    for (int i = 0; i < states.size(); i++) {
        IonTrap.set_spin_state(states[i]);

        std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

		for (auto traj_result : qsd_result) {
			qsd_results.push_back(traj_result);
		}
    }

    qsd::State psi_f = QSDSim.concatenate_qsd_states(IonTrap, qsd_results);
    double result = IonTrap.processor(psi_f);
    // // (*log_file) << "Result: " << std::setprecision(8) << result_2 << std::endl;

    return result;
}

void set_shelving_dynamics(experiment::ion_trap &IonTrap){
    IonTrap.pulse_cfg.scale_am = mode_cfg.get<double>("mode-8.pulse_scale");
    IonTrap.pulse_sequence = mode_cfg.get<std::string>("mode-8.pulse_sequence");
    IonTrap.prep_sequence = mode_cfg.get<std::string>("mode-8.prep_sequence");
    IonTrap.energy_drift = mode_cfg.get<double>("mode-8.energy_drift");
    IonTrap.shelving_flag = mode_cfg.get<int>("mode-8.flag");

    int len = 0;
    std::stringstream string_stream(IonTrap.pulse_sequence);
    std::string substr;
    while(string_stream.good())
    {
        std::getline(string_stream, substr, '_');
        len += 1;
    }

    IonTrap.num_pulses = len;
    std::cout<<IonTrap.num_pulses<< " Composite pulses" <<std::endl;
    double total_t = 0.0;

    std::stringstream string_stream_1(mode_cfg.get<std::string>("mode-8.pulses_duration"));
    while(string_stream_1.good())
    {
        std::string substr_1;
        std::getline(string_stream_1, substr_1, ',');
        if (substr_1 != ""){
            IonTrap.shelving_time.push_back(std::stod(substr_1));}
    }
}

void set_ion_trap_solution(experiment::ion_trap &IonTrap)
{
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

	if (IonTrap.pulse_cfg.profile_am == 4) {
		IonTrap.pulse_cfg.scale_am = solutions_cfg.get<double>("solutions.g.pulse_scale");
	}

	if (IonTrap.pulse_cfg.profile_am == 5) {
		IonTrap.pulse_cfg.scale_am = solutions_cfg.get<double>("solutions.sampled.pulse_scale");
	}

    IonTrap.Delta_L = solutions_cfg.get<double>("solutions.internal.Delta_L");
    IonTrap.Rabi_s = solutions_cfg.get<double>("solutions.internal.Rabi_s");

}

qsd::Expectation switch_fidelity_mode(int fidelity_mode, experiment::ion_trap &IonTrap, qsd_simulation<experiment::ion_trap> &QSDSim, bool save_per_traj=false, std::ofstream* log_file=NULL)
{
	qsd::Expectation result(0.0, 0.0, 0.0, 0.0);

	switch (fidelity_mode) {

		case 0: { // Numerical fidelity

			if (IonTrap.state_cfg.use_x_basis) {
				result = qsd::Expectation(fidelity_x_basis(IonTrap, QSDSim), 0.0, 0.0, 0.0);
			} else {
				std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
				result = IonTrap.processor_fidelity(qsd_result, save_per_traj, log_file);
			}

			break;
		}
//
//		case 1: { // Manning fidelity
//
//			// Analytic info
//			experiment::Fidelity _fidelity;
//			// // _fidelity.initialize();
//
//			result = qsd::Expectation(_fidelity.fidelity_1(IonTrap), 0.0, 0.0, 0.0);
//
//			break;
//		}
//
//		case 2: { // Magnus expansion
//
//			// Analytic info
//			experiment::Fidelity _fidelity;
//			_fidelity.initialize();
//
//			result = qsd::Expectation(_fidelity.fidelity_2(IonTrap), 0.0, 0.0, 0.0);
//
//			break;
//		}
//
//		case 3: { // Manning fidelity (est.)
//
//			// Analytic info
//			experiment::Fidelity _fidelity;
//			// // _fidelity.initialize();
//
//			result = qsd::Expectation(_fidelity.fidelity_1_estimate(IonTrap), 0.0, 0.0, 0.0);
//
//			break;
//		}

		case 5: { // Even population

			std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
			result  = IonTrap.processor_even_population(qsd_result);

			break;
		}

		default: {
			break;
		}
	}

	return result;
}

void run_through_processor_fidelity(int fidelity_mode, experiment::ion_trap &IonTrap, qsd_simulation<experiment::ion_trap> &QSDSim, qsd::Expectation &result, time_point &begin, time_point &end, bool save_per_traj=false, std::ofstream* log_file=NULL) {

	experiment::StateConfig state_cfg = IonTrap.state_cfg;
	if (state_cfg.motion_state_type == "PURE") {

		begin = std::chrono::steady_clock::now();
		result = switch_fidelity_mode(fidelity_mode, IonTrap, QSDSim, save_per_traj, log_file);
		end = std::chrono::steady_clock::now();

	} else if (state_cfg.motion_state_type == "THERMAL") {

		////////////////////////////////////////////////////
		// // (*log_file) << "Thermal state mode:" << std::endl;
		// // (*log_file) << "Total number of ions: " << N_ions << std::endl;
		// // (*log_file) << "Occupation probability cutoff: " << cutoff_probability << std::endl;

		// // (*log_file) << "Average phonon number per mode: ";
		// // for (auto n_bar : state_cfg.n_bar) {
		// // (*log_file) << n_bar << ", ";
		// // }
		// // (*log_file) << std::endl;

		// // (*log_file) << "Maximum phonon excitations allowed per mode: ";
		// // for (auto max_n_k : max_n) {
		// // (*log_file) << max_n_k << ", ";
		// // }
		// // (*log_file) << std::endl;

		// // (*log_file) << "Total number of states: " << permutations.size() << std::endl;
		// // (*log_file) << "Total number of states after applying cutoff: " << thermal_state.size() << std::endl;
		////////////////////////////////////////////////////

		// ...
		experiment::ion_trap::ThermalState thermal_state = IonTrap.prepare_thermal_state();

		qsd::Expectation avg_fidelity(0.0,0.0,0.0,0.0);
		auto begin = std::chrono::steady_clock::now();

		for (auto permutation : thermal_state.states) {

			state_cfg.motion = permutation.second;
			IonTrap.state_cfg = state_cfg;
			IonTrap.set_state();

			// // (*log_file) << "\nState = ";
			// // for (auto item : permutation.second) (*log_file) << item;
			// // (*log_file) << " Probability = " << permutation.first << std::endl;

			auto begin_i = std::chrono::steady_clock::now();

			qsd::Expectation result = switch_fidelity_mode(0, IonTrap, QSDSim);

			// // (*log_file) << "Result: " << '{' << result.mean.real() << ", " << std::sqrt(result.var.real()) << '}';

			auto end_i = std::chrono::steady_clock::now();
			// // (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count() << "s)\n";

			avg_fidelity.mean += permutation.first * result.mean.real();
		}
		auto end = std::chrono::steady_clock::now();
		// // (*log_file) << "\nAverage fidelity: " << avg_fidelity;
		// // (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s)\n";
	}
}

void single_run(const std::string& config_path, std::ofstream* log_file)
{
    experiment::ion_trap IonTrap(config_path);
    qsd_simulation<experiment::ion_trap> QSDSim(config_path);

	set_ion_trap_solution(IonTrap);

	int processor = mode_cfg.get<int>("mode-0.processor");
	bool save_per_traj = mode_cfg.get<bool>("mode-0.save_per_traj");
    int fidelity_mode = mode_cfg.get<int>("mode-0.fidelity.mode");

	if ((processor == 0) || (processor == 1)) // fidelity or infidelity
	{
		time_point begin;
		time_point end;
//        std::ofstream sim_log_file;
        std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");
//        sim_log_file.open(log_dir + "simulation_log.txt");
		qsd::Expectation result(0.0,0.0,0.0,0.0);

		if (save_per_traj) {

			std::ofstream log_file_processor;
//			std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");
			log_file_processor.open(log_dir + "fidelity.txt");

			run_through_processor_fidelity(fidelity_mode, IonTrap, QSDSim, result, begin, end, save_per_traj, &log_file_processor);

			log_file_processor.close();
		} else {
			run_through_processor_fidelity(fidelity_mode, IonTrap, QSDSim, result, begin, end);
            std::cout << "Result: " << std::setprecision(8) << '{' << result.mean.real() << ", " << std::sqrt(result.var.real()) << '}' << " (" << QSDSim.get_ntraj() << " trajectories)" << std::endl;
		}

		(*log_file) << "Result: " << std::setprecision(8) << '{' << result.mean.real() << ", " << std::sqrt(result.var.real()) << '}' << " (" << QSDSim.get_ntraj() << " trajectories)" << std::endl;
		print_time(begin, end, log_file);
	}

	if (processor == 2) // Parity
	{
		std::cout << "Not implemented!" << std::endl;
	}

    if (processor == 3) // Density matrix
    {
        experiment::StateConfig state_cfg = IonTrap.state_cfg;
        if (state_cfg.motion_state_type == "PURE") {
            time_point begin = std::chrono::steady_clock::now();

            std::vector<std::vector<std::vector<qsd::Expectation>>> density_matrix(4);
            std::vector<std::vector<std::vector<qsd::Expectation>>> reduced_density_matrix(2);

            std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
            density_matrix = IonTrap.processor_average_density_matrix(qsd_result);
            reduced_density_matrix = IonTrap.processor_reduced_average_density_matrix(qsd_result);

            int kmax = density_matrix[0][0].size();
//            std::cout << "here" <<kmax << std::endl;

            std::vector<std::complex<double>> B(kmax);
            std::vector<Eigen::MatrixXcd> SA(kmax);


            for (int i = 0; i < kmax; ++i) {
                SA[i].resize(4,4);
            }

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    for (int k = 0; k < kmax; ++k) {
                        SA[k](i,j) = density_matrix[i][j][k].mean;
                    }
                }
            }

            std::ofstream log_file_vneumann;
            std::ofstream DM;

            std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");
            log_file_vneumann.open(log_dir + "vonNeumann.txt");
            DM.open(log_dir + "DM.txt");

            B[0] = 0.0;
            log_file_vneumann << "von Neumann entropy" << std::endl << 0 << " " << B[0].real() <<" "<< B[0].imag() << std::endl;
            DM << "Density matrix" << std::endl;

            for (int i = 1; i < kmax; ++i){
//                B[i] =-( SA[i] * (SA[i].array().log()).matrix() ).trace(); //von Neuman entropy
//                B[i] =-( SA[i] * SA[i].log() ).trace(); //von Neuman entropy
//                log_file_vneumann << i << " " << B[i].real() <<" "<< B[i].imag() << std::endl;
                DM << i << std::endl << SA[i] << std::endl;
            }

            log_file_vneumann.close();
            DM.close();
            std::cout << "Computed von Neumann entropy with traced-out motional modes \n";

//            std::cout << "here" <<kmax << std::endl;

//            std::vector<std::complex<double>> red_B(kmax);
//            std::vector<Eigen::MatrixXcd> red_SA(kmax);
//            for (int i = 0; i < kmax; ++i) {
//                red_SA[i].resize(2,2);
//            }
//
//            for (int i = 0; i < 2; ++i) {
//                for (int j = 0; j < 2; ++j) {
//                    for (int k = 0; k < kmax; ++k) {
//                        red_SA[k](i,j) = reduced_density_matrix[i][j][k].mean;
//                    }
//                }
//            }
//
//            std::ofstream log_file_red_holevo;
//            std::ofstream red_DM;
//
//            log_file_red_holevo.open(log_dir + "red_holevo.txt");
//            red_DM.open(log_dir + "red_DM.txt");
//
//            std::cout << "Computed von Neumann entropy for reduced density matrix\n";
//            red_B[0] = 0.0;
//            log_file_red_holevo << "von Neumann Information" << std::endl << 0 << " " << red_B[0].real() <<" "<< red_B[0].imag() << std::endl;
//            red_DM << "Density matrix" << std::endl;
//
//            for (int i = 1; i < kmax; ++i) {
//                red_B[i] =-(red_SA[i] * (red_SA[i]).log() ).trace(); //von Neuman entropy
//                log_file_red_holevo << i << " " << red_B[i].real() <<" "<< red_B[i].imag() << std::endl;
//                red_DM << i << std::endl << red_SA[i] << std::endl;
//            }
//
//            log_file_red_holevo.close();
//            red_DM.close();
//            std::cout << "Log(A) \n" << B << "\n\n";
//            std::cout << "Trace = \n" << B.trace() << "\n\n";

            if (fidelity_mode==5){
                std::vector<qsd::Expectation> even_pop(kmax,qsd::Expectation(0.0, 0.0, 0.0, 0.0));         // Compute the even and odd populations
                std::vector<qsd::Expectation> odd_pop(kmax, qsd::Expectation(0.0, 0.0, 0.0, 0.0));
                for (int i = 0; i < kmax; ++i) {

                    even_pop[i] = qsd::Expectation(
                            1.0*density_matrix[0][0][i].mean.real() + 1.0*density_matrix[3][3][i].mean.real(), 0.0,
                            density_matrix[0][0][i].var.real() + density_matrix[3][3][i].var.real(), 0.0);
                    odd_pop[i] = qsd::Expectation(
                            density_matrix[0][3][i].mean.real() + density_matrix[3][0][i].mean.real(), 0.0,
                            density_matrix[0][3][i].var.real() + density_matrix[3][0][i].var.real(), 0.0);
//                std::cout << even_pop[i].mean.real() << " ; "<< even_pop[i].mean.imag() << std::endl;
                }
                (*log_file) << "Result: even_pop " << std::setprecision(8) << '{' << even_pop[kmax-1].mean.real() << ", " << std::sqrt(even_pop[kmax-1].var.real()) << '}' << " (" << QSDSim.get_ntraj() << " trajectories)" << std::endl;
            }
            else if (fidelity_mode==0){
                qsd::Expectation result(0.0, 0.0, 0.0, 0.0);
                result = IonTrap.processor_fidelity(qsd_result, save_per_traj, log_file);
                (*log_file) << "Result: " << std::setprecision(8) << '{' << result.mean.real() << ", " << std::sqrt(result.var.real()) << '}' << " (" << QSDSim.get_ntraj() << " trajectories)" << std::endl;
            }

            time_point end = std::chrono::steady_clock::now();
			print_time(begin, end, log_file);

		} else if (state_cfg.motion_state_type == "THERMAL") {
			std::cout << "Not implemented!" << std::endl;
		}
	}

	if (processor == 4) // Average trajectory
	{
		experiment::StateConfig state_cfg = IonTrap.state_cfg;
		if (state_cfg.motion_state_type == "PURE") {
			time_point begin = std::chrono::steady_clock::now();

			std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
			qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

			// store the results on disk
			save_trajectory_to_disk(avg_traj_result);

			time_point end = std::chrono::steady_clock::now();
			print_time(begin, end, log_file);

		} else if (state_cfg.motion_state_type == "THERMAL") {

			// ...
			time_point begin = std::chrono::steady_clock::now();
			experiment::ion_trap::ThermalState thermal_state = IonTrap.prepare_thermal_state();

			std::vector<double> weights;
			std::vector<qsd::TrajectoryResult> avg_traj_results;

			for (auto permutation : thermal_state.states) {

				state_cfg.motion = permutation.second;
				IonTrap.state_cfg = state_cfg;
				IonTrap.set_state();

				(*log_file) << "\nState = ";
				for (auto item : permutation.second) (*log_file) << item;
				(*log_file) << " Probability = " << permutation.first;

				time_point begin_i = std::chrono::steady_clock::now();

				std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

				// compute the average trajectory
				qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

				weights.push_back(permutation.first);
				avg_traj_results.push_back(avg_traj_result);

				time_point end_i = std::chrono::steady_clock::now();
				(*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count() << "s)\n";
			}

			// compute the average trajectory
			qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(avg_traj_results, weights);

			// store the results on disk
			save_trajectory_to_disk(avg_traj_result);

			time_point end = std::chrono::steady_clock::now();
			print_time(begin, end, log_file);
		}
	}

	if (processor == 5) // Even population
	{
		experiment::StateConfig state_cfg = IonTrap.state_cfg;
		if (state_cfg.motion_state_type == "PURE") {
			time_point begin = std::chrono::steady_clock::now();

			std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

			qsd::Expectation even_pop(0.0, 0.0, 0.0, 0.0);
			if (save_per_traj) {

				std::ofstream log_file_processor;
				std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");
				log_file_processor.open(log_dir + "even_pop.txt");

				even_pop  = IonTrap.processor_even_population(qsd_result, save_per_traj, &log_file_processor);

				log_file_processor.close();
			} else {
				even_pop  = IonTrap.processor_even_population(qsd_result);
			}

			(*log_file) << "Result: even_pop " << std::setprecision(8) << '{' << even_pop.mean.real() << ", " << std::sqrt(even_pop.var.real()) << '}' << " (" << QSDSim.get_ntraj() << " trajectories)" << std::endl;

			time_point end = std::chrono::steady_clock::now();
			print_time(begin, end, log_file);

		} else if (state_cfg.motion_state_type == "THERMAL") {
			std::cout << "Not implemented!" << std::endl;
		}
	}

	if (processor == 7) // EIT
	{
        std::cout<<"here3"<<std::endl;


        experiment::StateConfig state_cfg = IonTrap.state_cfg;
        if (state_cfg.motion_state_type == "PURE") {
            time_point begin = std::chrono::steady_clock::now();

            std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
            qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

            // store the results on disk
            save_trajectory_to_disk(avg_traj_result);

            time_point end = std::chrono::steady_clock::now();
            print_time(begin, end, log_file);

        } else if (state_cfg.motion_state_type == "THERMAL") {

            // ...
            time_point begin = std::chrono::steady_clock::now();
            experiment::ion_trap::ThermalState thermal_state = IonTrap.prepare_thermal_state();

            std::vector<double> weights;
            std::vector<qsd::TrajectoryResult> avg_traj_results;

            for (auto permutation : thermal_state.states) {

                state_cfg.motion = permutation.second;
                IonTrap.state_cfg = state_cfg;
                IonTrap.set_state();

                (*log_file) << "\nState = ";
                for (auto item : permutation.second) (*log_file) << item;
                (*log_file) << " Probability = " << permutation.first;

                time_point begin_i = std::chrono::steady_clock::now();

                std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

                // compute the average trajectory
                qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

                weights.push_back(permutation.first);
                avg_traj_results.push_back(avg_traj_result);

                time_point end_i = std::chrono::steady_clock::now();
                (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count() << "s)\n";
            }

            // compute the average trajectory
            qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(avg_traj_results, weights);

            // store the results on disk
            save_trajectory_to_disk(avg_traj_result);

            time_point end = std::chrono::steady_clock::now();
            print_time(begin, end, log_file);
        }
	}
}

void fidelity_vs_detuning(const std::string& config_path, std::ofstream* log_file)
{
    int type = mode_cfg.get<int>("simulation.mode-1.type");

    std::vector<int> list;
    switch (type) {
        case 0: {
            // Read the range
            std::vector<int> range;
            for (auto& item : mode_cfg.get_child("simulation.mode-1.range")) {
                range.emplace_back(item.second.get_value<int>());
            }

            for (int i = range[0]; i <= range[1]; i++) {
                list.push_back(i);
            }
            break;
        }

        case 1: {
            // Read the list
            std::stringstream string_stream(mode_cfg.get<std::string>("simulation.mode-1.list"));
            while(string_stream.good())
            {
                std::string substr;
                std::getline(string_stream, substr, ',');
                if (substr != "")
                    list.push_back(std::stoi(substr));
            }
            break;
        }

        default:
            std::cout << "Error: Unrecognized type for simulation mode 1.\n";
    }

    int n_processors = mode_cfg.get<int>("simulation.mode-1.n_processors");
	int mode = mode_cfg.get<int>("simulation.mode-1.mode");
	bool store_on_the_fly = mode_cfg.get<bool>("simulation.mode-1.store_on_the_fly");

    struct fidelity_result {
        int idx;
        qsd::Expectation fidelity;
        double time;

        fidelity_result(int _idx, qsd::Expectation _fidelity, double _time) : idx(_idx), fidelity(_fidelity), time(_time) {}
    };
    // // std::vector<std::vector<fidelity_result>> results(n_processors);
	std::vector<fidelity_result> results;

    // // auto worker = [&] (int worker_id, std::vector<int> idx_list) {
	auto worker = [&] (std::vector<int> idx_list) {

        for (int idx : idx_list) {

            try
            {
                // Initialize the IonTrap and the QSDSim
                experiment::ion_trap IonTrap(config_path);
                qsd_simulation<experiment::ion_trap> QSDSim(config_path);

                // Set delta and pulse
                IonTrap.set_delta(idx);
                IonTrap.set_pulse(idx);
				IonTrap.pulse_cfg.scale_am = mode_cfg.get<double>("simulation.mode-1.pulse_scale");

                auto begin_i = std::chrono::steady_clock::now();

                qsd::Expectation fidelity(0.0, 0.0, 0.0, 0.0);

                //// Experimental ////
				// // double fidelity = 0.0;
				switch (mode) {

					case 0: { // Numerical fidelity

						if (IonTrap.state_cfg.use_x_basis) {
							fidelity = qsd::Expectation(fidelity_x_basis(IonTrap, QSDSim), 0.0, 0.0, 0.0);
						} else {
							std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
							fidelity = IonTrap.processor_fidelity(qsd_result);
						}

						break;
					}

					default: {
						break;
					}
				}
                //////////////////////

                auto end_i = std::chrono::steady_clock::now();

				results.emplace_back(idx, fidelity, std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count());

				if (store_on_the_fly) {
					(*log_file) << "idx: " << idx << " " << '{' << fidelity.mean.real() << ", " << std::sqrt(fidelity.var.real()) << '}' << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count() << "s)" << std::endl;
				}
            }
            catch (const char* msg)
            {
				results.emplace_back(idx, qsd::Expectation(0.0, 0.0, 0.0, 0.0), 0.0);

				if (store_on_the_fly) {
					(*log_file) << "idx: " << idx << " " << 0.0 << std::endl;
				}

                std::cout << msg << std::endl;
            }
        }
    };

    auto begin = std::chrono::steady_clock::now();

    // Shuffle the list
	auto rng = std::default_random_engine {};
	std::shuffle(std::begin(list), std::end(list), rng);

    if (n_processors > 1) {

    	//// Experimental ////
		// Launch the pool with "n_processors" threads.
		boost::asio::thread_pool pool(n_processors);

		// Submit a function to the pool.
		for  (auto idx : list) {
			std::vector<int> idx_list = {idx};
			boost::asio::post(pool, std::bind(worker, idx_list));
		}

		// Wait for all tasks in the pool to complete.
		pool.join();
    	//////////////////////

    } else {
		worker(list);
    }

    // Sort the list
	std::sort(results.begin(), results.end(), [] (const fidelity_result& lhs, const fidelity_result& rhs) {
		return lhs.idx < rhs.idx;
	});

    // Store the results to the log_file
	if (store_on_the_fly) {
		(*log_file) << std::endl;
	}

	for (fidelity_result result : results) {
            (*log_file) << "idx: " << result.idx << " " << '{' << result.fidelity.mean.real() << ", " << std::sqrt(result.fidelity.var.real()) << '}' << " (took: " << result.time << "s)" << std::endl;
    }

    auto end = std::chrono::steady_clock::now();
    (*log_file) << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s\n";
}

void robustness_test(const std::string& config_path, std::ofstream* log_file)
{
	int mode = mode_cfg.get<int>("mode-3.mode");
    int type = mode_cfg.get<int>("mode-3.type");
    int n_processors = mode_cfg.get<int>("mode-3.n_processors");
	bool store_on_the_fly = mode_cfg.get<bool>("mode-3.store_on_the_fly");

    switch (type) {
        case 0: { // detuning error
            /*
            * Check against a DC offset
            * */

			// Read the steps
			int steps = mode_cfg.get<int>("mode-3.detuning_error.steps");

			// Read the offset range
			std::vector<double> range;
			for (auto& item : mode_cfg.get_child("mode-3.detuning_error.range")) {
				range.emplace_back(item.second.get_value<double>());
			}

			// Create a list of detunings
			double delta = (range[1] - range[0]) / (double) steps;

			std::vector<double> offsets;
			for (int i = 0; i <= steps; i++) {
				offsets.emplace_back(range[0] + i * delta);
			}

			struct robustness_result {
			  double offset;
			  qsd::Expectation fidelity;
			  double contrast;
			  double odd_pop;
			  double time;

			  robustness_result(double _offset, qsd::Expectation _fidelity, double _contrast, double _odd_pop, double _time) :
			  offset(_offset), fidelity(_fidelity), contrast(_contrast), odd_pop(_odd_pop), time(_time) {}
			};
			std::vector<robustness_result> results;

			auto worker = [&] (std::vector<double> offsets_list) {

			  // // for (int i = 0; i < offsets.size(); i++) {
			  for (double offset : offsets_list) {

				  try
				  {
					  // Initialize the IonTrap and the QSDSim
					  experiment::ion_trap IonTrap(config_path);
					  qsd_simulation<experiment::ion_trap> QSDSim(config_path);

					  // Set delta and pulse
					  set_ion_trap_solution(IonTrap);

					  // Get the central detuning
					  /*
					  std::vector<double> detuning(IonTrap.get_delta());

					  std::vector<double> _detuning;
					  for (auto val : detuning) {
						  _detuning.emplace_back(val + offset);
					  }

					  IonTrap.set_delta(_detuning);
					  */

					  //// Experimental ////
					  // Get the vibrational frequencies
					  std::vector<double> nu_list(IonTrap.nu_list);

					  // Read and apply the detuning error in the motional-mode frequencies [Mrad/s]
					  std::vector<double> _nu_list;
					  for (auto val : nu_list) {
						  _nu_list.emplace_back(val + offset);
					  }

					  IonTrap.nu_list = _nu_list;
					  //////////////////////

					  auto begin_i = std::chrono::steady_clock::now();

					  qsd::Expectation result_fidelity(0.0, 0.0, 0.0, 0.0);
					  double result_parity = 0.0;
					  double result_density_matrix = 0.0;

					  if (IonTrap.state_cfg.use_x_basis) {
						  result_fidelity = qsd::Expectation(fidelity_x_basis(IonTrap, QSDSim), 0.0, 0.0, 0.0);
					  } else {
						  std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

						  // Experimental //
						  std::ofstream log_file_processor;
						  std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");
						  log_file_processor.open(log_dir + "fidelity_" + std::to_string(offset) + ".txt");

						  bool save_per_traj = true;

						  result_fidelity = IonTrap.processor_fidelity(qsd_result, save_per_traj, &log_file_processor);

						  log_file_processor.close();
						  // Experimental //

						  // // result_parity = IonTrap.processor_parity(qsd_result)[0];
						  // // result_density_matrix = IonTrap.processor_density_matrix(qsd_result)[0];
					  }

					  auto end_i = std::chrono::steady_clock::now();

					  if (store_on_the_fly) {
						  (*log_file) << "Offset: " << offset << " fidelity: " << '{' << result_fidelity.mean.real() << ", " << std::sqrt(result_fidelity.var.real()) << '}' << " contrast: " << result_parity << " odd pop: " << result_density_matrix;
						  (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i-begin_i).count() << "s)" << std::endl;
					  }

					  results.emplace_back(offset, result_fidelity, result_parity, 0.0, std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count());
				  }
				  catch (const char* msg)
				  {
					  results.emplace_back(offset, qsd::Expectation(0.0, 0.0, 0.0, 0.0), 0.0, 0.0, 0.0);

					  if (store_on_the_fly) {
						  std::cout << msg << std::endl;
					  }
				  }
			  }

			};

			auto begin = std::chrono::steady_clock::now();

			/////
			if (n_processors > 1) {

				//// Experimental ////
				// Launch the pool with "n_processors" threads.
				boost::asio::thread_pool pool(n_processors);

				// Submit a function to the pool.
				for (int i = 0; i < offsets.size(); i++) {
					std::vector<double> offsets_list = {offsets[i]};
					boost::asio::post(pool, std::bind(worker, offsets_list));
				}

				// Wait for all tasks in the pool to complete.
				pool.join();
				//////////////////////

			} else {
				worker(offsets);
			}
			/////

			// Sort the list
			std::sort(results.begin(), results.end(), [] (const robustness_result& lhs, const robustness_result& rhs) {
			  return lhs.offset < rhs.offset;
			});

			// Store the results to the log_file
			if (store_on_the_fly) {
				(*log_file) << std::endl;
			}

			for (robustness_result worker_results : results) {
//				for (robustness_result result : worker_results)
				(*log_file) << std::setprecision(8) << "Offset: " << worker_results.offset << " fidelity: " << '{' << worker_results.fidelity.mean.real() << ", " << std::sqrt(worker_results.fidelity.var.real()) << '}' << " contrast: " << worker_results.contrast << " odd pop: " << worker_results.odd_pop << " (took: " << worker_results.time << "s)" << std::endl;
			}

			auto end = std::chrono::steady_clock::now();
			(*log_file) << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;

            break;
        }

        case 1: { // pulse scale
            /*
            * Check against a scale in the pulse
            * */

			// Read the scale list
			std::vector<double> scales;

			int test_type = mode_cfg.get<int>("mode-3.test.type");
			switch (test_type) {

				case 0: { // range

					// Read the steps
					int steps = mode_cfg.get<int>("mode-3.test.steps");

					// Read the range
					std::vector<double> range;
					for (auto& item : mode_cfg.get_child("mode-3.test.range")) {
						range.emplace_back(item.second.get_value<double>());
					}

					// Create a list of pulse scales
					double delta = (range[1] - range[0]) / (double) steps;

					for (int i = 0; i <= steps; i++) {
						scales.emplace_back(range[0] + i * delta);
					}

					break;
				}

				case 1: { // list

					std::stringstream string_stream(mode_cfg.get<std::string>("mode-3.test.list"));
					while(string_stream.good())
					{
						std::string substr;
						std::getline(string_stream, substr, ',');
						if (substr != "")
							scales.push_back(std::stod(substr));
					}

					break;
				}

				default:
				std::cout << "Error: Unrecognized type for simulation mode 3.\n";
			}

            struct robustness_result {
                double scale;
                qsd::Expectation fidelity;
                double time;

                robustness_result(double _scale, qsd::Expectation _fidelity, double _time) : scale(_scale), fidelity(_fidelity), time(_time) {}
            };
            std::vector<std::vector<robustness_result>> results(n_processors);

            auto worker = [&] (int worker_id, std::vector<double> scale_list) {

                for (double scale : scale_list) {
                    try
                    {
						time_point begin_i;
						time_point end_i;

                        // Initialize the IonTrap and the QSDSim
                        experiment::ion_trap IonTrap(config_path);
                        qsd_simulation<experiment::ion_trap> QSDSim(config_path);

                        // Set delta and pulse
						set_ion_trap_solution(IonTrap);
                        IonTrap.pulse_cfg.scale_am = scale;

						qsd::Expectation fidelity(0.0,0.0,0.0,0.0);
						run_through_processor_fidelity(mode, IonTrap, QSDSim, fidelity, begin_i, end_i);

						if (store_on_the_fly) {
							(*log_file) << "Scale: " << scale << " " << '{' << fidelity.mean.real() << ", " << std::sqrt(fidelity.var.real()) << '}';
							(*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count() << "s)" << std::endl;
						}

                        results[worker_id].emplace_back(scale, fidelity, std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count());
                    }
                    catch (const char* msg)
                    {
                        results[worker_id].emplace_back(scale, qsd::Expectation(0.0, 0.0, 0.0, 0.0), 0.0);

                        std::cout << msg << std::endl;
                    }
                }
            };

            auto begin = std::chrono::steady_clock::now();

            if (n_processors > 1) {
                std::vector<std::thread> threads(n_processors);
                const int batch_size = scales.size() / n_processors;

                for (int k = 0; k < threads.size() - 1; k++) {
                    std::cout << batch_size << " task(s) added\n";
                    std::vector<double> scale_list = {scales.begin() + k * batch_size, scales.begin() + (k + 1) * batch_size };
                    threads[k] = std::thread(worker, k, scale_list);
                }
                std::cout << scales.size() - int((n_processors - 1) * batch_size) << " task(s) added\n";
                std::vector<double> scale_list = {scales.begin() + (n_processors - 1) * batch_size, scales.end()};
                threads.back() = std::thread(worker, threads.size() - 1, scale_list);

                for(auto&& i : threads) {
                    i.join();
                }
            } else {
                worker(0, scales);
            }

            // Store the results to the log_file
			if (store_on_the_fly) {
				(*log_file) << std::endl;
			}

			for (std::vector<robustness_result> worker_results : results) {
                for (robustness_result result : worker_results)
                    (*log_file) << "scale: " << result.scale << " " << '{' << result.fidelity.mean.real() << ", " << std::sqrt(result.fidelity.var.real()) << '}' << " (took: " << result.time << "s)" << std::endl;
            }

            auto end = std::chrono::steady_clock::now();
            (*log_file) << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s\n";

            break;
        }

		case 2: { // time error

			// Read dt and numdts
			double dt = mode_cfg.get<double>("mode-3.time_error.dt");
			int numdts = mode_cfg.get<int>("mode-3.time_error.numdts");

			// Read the numsteps range
			std::vector<int> range;
			for (auto& item : mode_cfg.get_child("mode-3.time_error.numsteps")) {
				range.emplace_back(item.second.get_value<int>());
			}

			// Create a list of numsteps
			std::vector<int> numsteps_list;
			for (int i = range[0]; i <= range[1]; i++) {
				numsteps_list.push_back(i);
			}

			struct robustness_result {
			  double t_gate;
			  qsd::Expectation fidelity;
			  double contrast;
			  double odd_pop;
			  double time;

			  robustness_result(double _t_gate, qsd::Expectation _fidelity, double _contrast, double _odd_pop, double _time) :
					  t_gate(_t_gate), fidelity(_fidelity), contrast(_contrast), odd_pop(_odd_pop), time(_time) {}
			};
			std::vector<std::vector<robustness_result>> robustness_results_workers;

			auto worker = [&] (int worker_id, std::vector<int> _numsteps_list) {

			  for (int numsteps : _numsteps_list) {

				  try
				  {
					  // Initialize the IonTrap and the QSDSim
					  experiment::ion_trap IonTrap(config_path);
					  qsd_simulation<experiment::ion_trap> QSDSim(config_path);

					  // Set delta and pulse
					  set_ion_trap_solution(IonTrap);

					  // Set the gate time
					  IonTrap.set_time(dt, numdts, numsteps);
					  double t_gate = IonTrap.t_gate;

					  auto begin_i = std::chrono::steady_clock::now();

					  qsd::Expectation result_fidelity(0.0, 0.0, 0.0, 0.0);
					  double result_parity = 0.0;
					  double result_density_matrix = 0.0;

					  if (IonTrap.state_cfg.use_x_basis) {
						  result_fidelity = qsd::Expectation(fidelity_x_basis(IonTrap, QSDSim), 0.0, 0.0, 0.0);
					  } else {
						  std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

						  result_fidelity = IonTrap.processor_fidelity(qsd_result);
						  // // result_parity = IonTrap.processor_parity(qsd_result)[0];
						  // // result_density_matrix = IonTrap.processor_density_matrix(qsd_result)[0];
					  }

					  auto end_i = std::chrono::steady_clock::now();

					  if (store_on_the_fly) {
						  (*log_file) << "t_gate: " << t_gate << " fidelity: " << '{' << result_fidelity.mean.real() << ", " << std::sqrt(result_fidelity.var.real()) << '}' << " contrast: " << result_parity << " odd pop: " << result_density_matrix;
						  (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i-begin_i).count() << "s)" << std::endl;
					  }

					  robustness_results_workers[worker_id].emplace_back(t_gate, result_fidelity, result_parity, result_density_matrix, std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count());
				  }
				  catch (const char* msg)
				  {
					  robustness_results_workers[worker_id].emplace_back(0.0, qsd::Expectation(0.0, 0.0, 0.0, 0.0), 0.0, 0.0, 0.0);

					  std::cout << msg << std::endl;
				  }
			  }

			};

			auto begin = std::chrono::steady_clock::now();

			if (n_processors > 1) {

				robustness_results_workers.resize(numsteps_list.size());

				// Launch the pool with "n_processors" threads.
				boost::asio::thread_pool pool(n_processors);

				// Submit a function to the pool.
				for (int worker_id = 0; worker_id < numsteps_list.size(); worker_id++) {
					std::vector<int> _numsteps_list = {numsteps_list[worker_id]};
					boost::asio::post(pool, std::bind(worker, worker_id, _numsteps_list));
				}

				// Wait for all tasks in the pool to complete.
				pool.join();

			} else {
				robustness_results_workers.resize(1);
				worker(0, numsteps_list);
			}

			// Store the results to the log_file
			if (store_on_the_fly) {
				(*log_file) << std::endl;
			}

			for (auto worker_result : robustness_results_workers) {
				for (auto result : worker_result) {
						(*log_file) << "t_gate: " << result.t_gate << " fidelity: " << '{' << result.fidelity.mean.real() << ", " << std::sqrt(result.fidelity.var.real()) << '}' << " contrast: " << result.contrast << " odd pop: " << result.odd_pop << " (took: " << result.time << "s)" << std::endl;
				}
			}

			auto end = std::chrono::steady_clock::now();
			(*log_file) << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;

			break;
		}

		default: {
			break;
		}

    }
}

void playground(const std::string& config_path, std::ofstream* log_file)
{
	experiment::ion_trap IonTrap(config_path);
	qsd_simulation<experiment::ion_trap> QSDSim(config_path);

	int mode = mode_cfg.get<int>("simulation.mode-5.mode");

	// set pulse and delta
	set_ion_trap_solution(IonTrap);

	// Test the "get_delta_new()" implementation
	/*
	double t = 0.0;
	double dt = 2.0;

	while (t <= 200.0) {

		std::cout << IonTrap.get_delta_new(t) << std::endl;
		t += dt;
	}
	*/

	// Test the "get_omega()" implementation
	/*
	double t = 0.0;
	double dt = 1.0;

	while (t <= 100.0) {

		std::cout << IonTrap.get_omega(0, t) << std::endl;
		t += dt;
	}
	 */

}


int main(int argc, char *argv[])
{
    // Read the config file path
    std::string config_path;
    if (argc == 1) {
        std::cout << "Error: Config file path is missing.\n";
        std::cout << "IonTrap-QSD config_file_path\n";
        return 0;
    } else {
        config_path = std::string(argv[1]);
    }

    // Read the simulation configurations
	boost::property_tree::read_info(config_path + "trap.cfg", trap_cfg);
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
//    std::cout<<simulation_mode<<std::endl;

    switch (simulation_mode) {
        case 0: {
            /*
             * Single run
             * */

			boost::property_tree::read_info(config_path + "mode-0.cfg", mode_cfg);

			// Read the trap solution config
			boost::property_tree::read_info(config_path + "solutions.cfg", solutions_cfg);

            single_run(config_path, &log_file);

            break;
        }

        case 1: {
            /*
            * Fidelity vs. detuning
            * */

			boost::property_tree::read_info(config_path + "mode-1.cfg", mode_cfg);
            fidelity_vs_detuning(config_path, &log_file);


            break;
        }

        case 3: {
            /*
            * Robustness check
            * */

			boost::property_tree::read_info(config_path + "mode-3.cfg", mode_cfg);

			// Read the trap solution config
			boost::property_tree::read_info(config_path + "solutions.cfg", solutions_cfg);

            robustness_test(config_path, &log_file);

            break;
        }

		case 5: {
			/*
			 * Playground
			 */

			boost::property_tree::read_info(config_path + "mode-5.cfg", mode_cfg);
			playground(config_path, &log_file);

			break;
		}

        default: {
            std::cout << "Error: Unknown simulation mode.\n";
        }
    }

    log_file.close();

    return 0;
}