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

    IonTrap.Delta_L = solutions_cfg.get<double>("solutions.internal.Delta_L");;
    IonTrap.Rabi_s = solutions_cfg.get<double>("solutions.internal.Rabi_s");;

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

void shelving(const std::string& config_path, std::ofstream* log_file)
{
    int processor = mode_cfg.get<int>("mode-8.processor");
    bool save_per_traj = mode_cfg.get<bool>("mode-8.save_per_traj");
    int type = mode_cfg.get<int>("mode-8.type");

    experiment::ion_trap IonTrap(config_path);
    set_shelving_dynamics(IonTrap);

    qsd_simulation<experiment::ion_trap> QSDSim(config_path);

    experiment::StateConfig state_cfg = IonTrap.state_cfg;
    std::vector<int> list;

    switch (type) {
        case 0: {
            // Read the range
            std::vector<int> range;
            for (auto& item : mode_cfg.get_child("mode-8.range")) {
                range.emplace_back(item.second.get_value<int>());
            }

            for (int i = range[0]; i <= range[1]; i++) {
                list.push_back(i);
            }

            break;
        }

        case 1: {
            // Read the list
            std::stringstream string_stream(mode_cfg.get<std::string>("mode-8.list"));
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

        IonTrap.set_shelving_phase(idx);
        IonTrap.set_delta(idx);

//        for (int i = 0; i < 30; ++i) {
//            std::cout<<IonTrap.shaving_phase_list[0][i]<<std::endl;
//        }

        IonTrap.pulse_cfg.scale_am = mode_cfg.get<double>("mode-8.pulse_scale");
        IonTrap.set_shelving_pulse(idx);

        if (state_cfg.motion_state_type == "PURE") {
            time_point begin = std::chrono::steady_clock::now();

            std::vector<qsd::TrajectoryResult> qsd_result = QSDSim.run(IonTrap);

            qsd::TrajectoryResult avg_traj_result = IonTrap.processor_average_trajectory(qsd_result);

            // store the results on disk
            save_trajectory_to_disk_freq(avg_traj_result, idx);
            if (processor==3){
                std::vector<std::vector<std::vector<qsd::Expectation>>> density_matrix(4);

                int transitions = 2;

                for (int j = 0; j < transitions; ++j) {
                    density_matrix = IonTrap.processor_average_multilevel_density_matrix(qsd_result, j);

                    int kmax = density_matrix[0][0].size();
                    int c= 2;

                    std::vector<Eigen::MatrixXcd> SA(kmax);
                    for (int i = 0; i < kmax; ++i) {
                        SA[i].resize(c,c);
                    }

                    for (int i = 0; i < c; ++i) {
                        for (int j = 0; j < c; ++j) {
                            for (int k = 0; k < kmax; ++k) {
                                SA[k](i,j) = density_matrix[i][j][k].mean;
                            }
                        }
                    }

                    std::cout << "Computed average density matrix with traced-out motional modes. Ion 0, Transition: " + std::to_string(j) + "\n";
                    std::ofstream DM;

                    std::string log_dir = experiment_cfg.get<std::string>("experiment.log_dir");
                    DM.open(log_dir + "DM" + std::to_string(j) + ".txt");

                    for (int i = 0; i < kmax; ++i) {
                        DM << i << std::endl << SA[i] << std::endl;
                    }
                    DM.close();
                }
            }

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

        case 8:{
            /*
            * Shelving fidelities
            */

            boost::property_tree::read_info(config_path + "mode-8.cfg", mode_cfg);

            // Read the trap solution config
            boost::property_tree::read_info(config_path + "solutions.cfg", solutions_cfg);

            shelving(config_path, &log_file);

            break;

        }

        default: {
            std::cout << "Error: Unknown simulation mode.\n";
        }
    }

    log_file.close();

    return 0;
}