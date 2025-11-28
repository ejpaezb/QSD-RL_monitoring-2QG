/*
 * Copyright (c) 2023, Eduardo J. Paez
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

#include "qsd/qsd_simulation.h"
#include <Eigen/LU>
//#include <Eigen/Dense>
//#include <Eigen/SparseCore>
#include <unsupported/Eigen/MatrixFunctions>

#include "ion_trap.h"
#include "utilities.h"
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

void set_cooling_dynamics(experiment::ion_trap &IonTrap){
    IonTrap.set_delta(solutions_cfg.get<int>("solutions.internal.delta_idx"));
    double pump = mode_cfg.get<double>("mode-7.pump");
    double scan = mode_cfg.get<double>("mode-7.scan");
    double det = mode_cfg.get<double>("mode-7.Delta_pump");
    IonTrap.set_EIT_pulse({pump,scan});
    IonTrap.set_Delta_pump(det);
}

void playground(const std::string& config_path, std::ofstream* log_file)
{
    experiment::ion_trap IonTrap(config_path);
    qsd_simulation<experiment::ion_trap> QSDSim(config_path);

    int mode = mode_cfg.get<int>("simulation.mode-5.mode");

    // set pulse and delta
//    set_ion_trap_solution(IonTrap);

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

void DOPPLER(const std::string& config_path, std::ofstream* log_file)
{
    int processor = mode_cfg.get<int>("mode-7.processor");
    int n_processors = mode_cfg.get<int>("mode-7.n_processors");

    bool save_per_traj = mode_cfg.get<bool>("mode-7.save_per_traj");
    int type = mode_cfg.get<int>("mode-7.type");

    experiment::ion_trap IonTrap(config_path);
    qsd_simulation<experiment::ion_trap> QSDSim(config_path);
    set_cooling_dynamics(IonTrap);
    experiment::StateConfig state_cfg = IonTrap.state_cfg;
    std::vector<int> list;

    switch (type) {
        case 0: {
            // Read the range
            std::vector<int> range;
            for (auto& item : mode_cfg.get_child("mode-7.range")) {
                range.emplace_back(item.second.get_value<int>());
            }

            for (int i = range[0]; i <= range[1]; i++) {
                list.push_back(i);
            }

            break;
        }

        case 1: {
            // Read the list
            std::stringstream string_stream(mode_cfg.get<std::string>("mode-7.list"));
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
    auto worker = [&] (std::vector<int> idx_list) {

        for (int idx: idx_list) {

            // Set delta and pulse
            IonTrap.set_delta(idx);
            IonTrap.pulse_cfg.scale_am = mode_cfg.get<double>("mode-7.pulse_scale");

            if (state_cfg.motion_state_type == "PURE") {
                time_point begin = std::chrono::steady_clock::now();

                std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
                qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

                // store the results on disk
                save_trajectory_to_disk_freq(avg_traj_result, idx);

                time_point end = std::chrono::steady_clock::now();
                print_time(begin, end, log_file);

            } else if (state_cfg.motion_state_type == "THERMAL") {

                // ...
                time_point begin = std::chrono::steady_clock::now();
                experiment::ion_trap::ThermalState thermal_state = IonTrap.prepare_thermal_state();

                std::vector<double> weights;
                std::vector<qsd::TrajectoryResult> avg_traj_results;

                for (auto permutation: thermal_state.states) {

                    state_cfg.motion = permutation.second;
                    IonTrap.state_cfg = state_cfg;
                    IonTrap.set_state();

                    (*log_file) << "\nState = ";
                    for (auto item: permutation.second) (*log_file) << item;
                    (*log_file) << " Probability = " << permutation.first;

                    time_point begin_i = std::chrono::steady_clock::now();

                    std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

                    // compute the average trajectory
                    qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

                    weights.push_back(permutation.first);
                    avg_traj_results.push_back(avg_traj_result);

                    time_point end_i = std::chrono::steady_clock::now();
                    (*log_file) << " (took: "
                                << std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count()
                                << "s)\n";
                }
            }
        }
    };

    auto begin = std::chrono::steady_clock::now();

    // Shuffle the list
//    auto rng = std::default_random_engine {};
//    std::shuffle(std::begin(list), std::end(list), rng);

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

    auto end = std::chrono::steady_clock::now();
}


void EIT(const std::string& config_path, std::ofstream* log_file)
{
    int processor = mode_cfg.get<int>("mode-7.processor");
    bool save_per_traj = mode_cfg.get<bool>("mode-7.save_per_traj");
    int type = mode_cfg.get<int>("mode-7.type");

    experiment::ion_trap IonTrap(config_path);
    qsd_simulation<experiment::ion_trap> QSDSim(config_path);
    set_cooling_dynamics(IonTrap);
    experiment::StateConfig state_cfg = IonTrap.state_cfg;
    std::vector<int> list;

    switch (type) {
        case 0: {
            // Read the range
            std::vector<int> range;
            for (auto& item : mode_cfg.get_child("mode-7.range")) {
                range.emplace_back(item.second.get_value<int>());
            }

            for (int i = range[0]; i <= range[1]; i++) {
                list.push_back(i);
            }

            break;
        }

        case 1: {
            // Read the list
            std::stringstream string_stream(mode_cfg.get<std::string>("mode-7.list"));
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

    for (int idx:list){

        // Set delta and pulse
        IonTrap.set_delta(idx);
        IonTrap.pulse_cfg.scale_am = mode_cfg.get<double>("mode-7.pulse_scale");

        if (state_cfg.motion_state_type == "PURE") {
            time_point begin = std::chrono::steady_clock::now();

            std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);
            qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

            // store the results on disk
            save_trajectory_to_disk_freq(avg_traj_result, idx);

            time_point end = std::chrono::steady_clock::now();
            print_time(begin, end, log_file);

        } else if (state_cfg.motion_state_type == "THERMAL") {

            // ...
            time_point begin = std::chrono::steady_clock::now();
            experiment::ion_trap::ThermalState thermal_state = IonTrap.prepare_thermal_state();

            std::vector<double> weights;
            std::vector<qsd::TrajectoryResult> avg_traj_results;

            for (auto permutation: thermal_state.states) {

                state_cfg.motion = permutation.second;
                IonTrap.state_cfg = state_cfg;
                IonTrap.set_state();

                (*log_file) << "\nState = ";
                for (auto item: permutation.second) (*log_file) << item;
                (*log_file) << " Probability = " << permutation.first;

                time_point begin_i = std::chrono::steady_clock::now();

                std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

                // compute the average trajectory
                qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

                weights.push_back(permutation.first);
                avg_traj_results.push_back(avg_traj_result);

                time_point end_i = std::chrono::steady_clock::now();
                (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count()
                            << "s)\n";
            }
        }
    }
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
        case 5: {
            /*
             * Playground
             */

            boost::property_tree::read_info(config_path + "mode-5.cfg", mode_cfg);
            playground(config_path, &log_file);

            break;
        }
        case 7: {
            /*
             * Cooling dynamics
             */
            boost::property_tree::read_info(config_path + "mode-7.cfg", mode_cfg);
            boost::property_tree::read_info(config_path + "solutions.cfg", solutions_cfg);
            int processor = mode_cfg.get<int>("mode-7.processor");

            switch (processor){
                case 0: {
                    EIT(config_path, &log_file);
                    break;
                }
                case 1: {
                    DOPPLER(config_path, &log_file);
                    break;
                }

                default: {
                    std::cout << "Error: Unknown processor for cooling dynamics.\n";
                }
            }
            // Read the trap solution config

            break;
        }

        default: {
            std::cout << "Error: Unknown simulation mode.\n";
        }
    }

    log_file.close();

    return 0;
}