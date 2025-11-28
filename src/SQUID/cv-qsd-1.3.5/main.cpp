/*
 * Copyright (c) 2019, Seyed Shakib Vedaie & Eduardo J. Paez
 */

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/foreach.hpp>

#include <qsd/qsd_simulation.h>
#include "cv_cqed.h"

boost::property_tree::ptree tree;

void single_run(const std::string& config_path, ofstream* log_file)
{
    experiment::cv_cqed CvQed(config_path);
    qsd_simulation<experiment::cv_cqed> QSDSim(config_path);

    bool use_external_delta = tree.get<bool>("experiment.simulation.mode-0.use_external_delta");
    if (use_external_delta) {
        /*
        * Load the external delta
        * */
        double external_delta;
        for (auto& item : tree.get_child("experiment.simulation.mode-0.external_delta")) {
            external_delta = item.second.get_value<double>();
        }

        CvQed.set_delta(external_delta);
    } else {

        CvQed.set_delta(tree.get<int>("experiment.simulation.mode-0.delta_idx"));
    }

    bool use_external_pulse = tree.get<bool>("experiment.simulation.mode-0.use_external_pulse");
    if (use_external_pulse) {
        /*
        * Load the external pulse
        * */
        std::vector<double> external_pulse;
        for (auto& item : tree.get_child("experiment.simulation.mode-0.external_pulse")) {
            external_pulse.emplace_back(item.second.get_value<double>());
        }

        CvQed.set_pulse(external_pulse);
    } else {

        CvQed.set_pulse(tree.get<int>("experiment.simulation.mode-0.pulse_idx"));
    }

    auto begin = std::chrono::high_resolution_clock::now();

    (*log_file) << "Result: " << QSDSim.run(CvQed);

    auto end = std::chrono::high_resolution_clock::now();
    (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s)\n";
}

void fidelity_vs_detuning(const std::string& config_path, ofstream* log_file)
{
    experiment::cv_cqed CvQed(config_path);
    qsd_simulation<experiment::cv_cqed> QSDSim(config_path);

    // Read the range
    std::vector<int> range;
    for (auto& item : tree.get_child("experiment.simulation.mode-1.range")) {
        range.emplace_back(item.second.get_value<int>());
    }

    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<double> fidelity_list;
    for (int i = range[0]; i <= range[1]; i++) {
        try
        {
            // Set delta and pulse
            CvQed.set_delta(i);
            CvQed.set_pulse(i);

            auto begin_i = std::chrono::high_resolution_clock::now();

            double fidelity = QSDSim.run(CvQed);
            fidelity_list.push_back(fidelity);

            auto end_i = std::chrono::high_resolution_clock::now();
            //(*log_file) << "Progress: " << i + 1 << '/' << CvQed.get_delta_list().size() << " (took: " <<
            //    std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count() << "s)" << std::endl;

            (*log_file) << "idx: " << i << " " << fidelity << " (took: " <<
                        std::chrono::duration_cast<std::chrono::seconds>(end_i - begin_i).count() << "s)" << std::endl;
        }
        catch (const char* msg)
        {
            fidelity_list.push_back(0.0);
            (*log_file) << "idx: " << i << " " << 0.0 << std::endl;

            cout << msg << std::endl;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    (*log_file) << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s\n";
}

void robustness_test(const std::string& config_path, ofstream* log_file)
{
    experiment::cv_cqed CvQed(config_path);
    qsd_simulation<experiment::cv_cqed> QSDSim(config_path);

    int robustness_test_type = tree.get<int>("experiment.simulation.mode-3.robustness_test_type");
    switch (robustness_test_type) {
        case 0: {
            /*
            * Check against a DC offset
            * */
            std::string type = tree.get<std::string>("experiment.simulation.mode-3.DC.type");

            if (type == "DETUNING") {

                CvQed.set_delta(tree.get<int>("experiment.simulation.mode-3.DC.delta_idx"));

                bool use_external_pulse = tree.get<bool>("experiment.simulation.mode-3.DC.use_external_pulse");
                if (use_external_pulse) {
                    /*
                    * Load the external pulse
                    * */
                    std::vector<double> external_pulse;
                    for (auto& item : tree.get_child("experiment.simulation.mode-3.DC.external_pulse")) {
                        external_pulse.emplace_back(item.second.get_value<double>());
                    }

                    CvQed.set_pulse(external_pulse);
                } else {

                    CvQed.set_pulse(tree.get<int>("experiment.simulation.mode-3.DC.pulse_idx"));
                }

                // Get the central detuning
                double detuning = CvQed.get_delta();

                // Read the steps
                int steps = tree.get<int>("experiment.simulation.mode-3.DC.steps");

                // Read the offset range
                std::vector<double> range;
                for (auto& item : tree.get_child("experiment.simulation.mode-3.DC.range")) {
                    range.emplace_back(item.second.get_value<double>());
                }

                // Create a list of detunings
                double delta = (range[1] - range[0]) / (double) steps;

                std::vector<double> offsets;
                for (int i = 0; i <= steps; i++) {
                    offsets.emplace_back(range[0] + i * delta);
                }

                for (int i = 0; i < offsets.size(); i++) {

                    CvQed.set_delta(detuning + offsets[i]);

                    auto begin = std::chrono::high_resolution_clock::now();

                    (*log_file) << "Offset: " << offsets[i] << " " << QSDSim.run(CvQed);

                    auto end = std::chrono::high_resolution_clock::now();
                    (*log_file) << " (took: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s)\n";
                }

            } else if (type == "PULSE") {

            } else if (type == "PHASE") {

            }

            break;
        }

    }
}

int main(int argc, char *argv[])
{
    // Read the config file path
    std::string config_path;
    if (argc == 1) {
        std::cout << "Error: Config file path is missing.\n";
        std::cout << "CvQed-QSD config_file_path\n";
        return 0;
    } else {
        config_path = std::string(argv[1]);
    }

    // Read the simulation configurations
    boost::property_tree::read_info(config_path, tree);

    // Open the log file
    ofstream log_file;
    log_file.open(tree.get<std::string>("experiment.log_file"));

    int simulation_mode = tree.get<int>("experiment.simulation.mode");

    switch (simulation_mode) {
        case 0: {
            /*
             * Single run
             * */

            single_run(config_path, &log_file);

            break;
        }

        case 1: {
            /*
            * Fidelity vs. detuning
            * */

            fidelity_vs_detuning(config_path, &log_file);

            break;
        }

        case 3: {
            /*
            * Robustness check
            * */

            robustness_test(config_path, &log_file);

            break;
        }

        default: {
            std::cout << "Error: Unknown simulation mode.\n";
        }
    }

    return 0;
}