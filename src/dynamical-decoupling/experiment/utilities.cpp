#include "utilities.h"

namespace experiment {

    // Utility methods
    void load_1d_list(std::string filename, std::vector<double> &_list) {
        // Clear the delta_list
        std::vector<double>().swap(_list);

        // load the data file into a new array
        cnpy::NpyArray npy_data = cnpy::npy_load(filename);
        auto data_ptr = npy_data.data<double>();

        for (auto i = 0; i < npy_data.num_vals; i++) {
            _list.push_back(data_ptr[i]);
        }
    }

    void load_2d_list(std::string filename, std::vector<std::vector<double>> &_list) {
        // Clear the pulse_list
        std::vector<std::vector<double>>().swap(_list);

        // load the data file into a new array
        cnpy::NpyArray npy_data = cnpy::npy_load(filename);
        auto data_ptr = npy_data.data<double>();

        int nrows = npy_data.shape[0];
        int ncols = npy_data.shape[1];

//        std::cout<<filename<<"   "<<nrows<<" rows "<<ncols<<std::endl;

        _list.reserve(nrows);
        for (int row = 0; row < nrows; row++) {
            _list.emplace_back(ncols);
            for (int col = 0; col < ncols; col++) {
                _list[row][col] = data_ptr[ncols * row + col];
//                std::cout<<data_ptr[ncols * row + col]<<std::endl;
            }
        }
    }

    void parse_tree(boost::property_tree::ptree& tree, std::string path, std::vector<int>& list) {
        // Clear the list
        std::vector<int>().swap(list);

        std::stringstream stream(tree.get<std::string>(path));
        while (stream.good()) {
            std::string substr;
            std::getline(stream, substr, ',');
            if (substr != "")
                list.push_back(std::stoi(substr));
        }
    }

    void parse_tree(boost::property_tree::ptree& tree, std::string path, std::vector<double>& list) {
        // Clear the list
        std::vector<double>().swap(list);

        std::stringstream stream(tree.get<std::string>(path));
        while (stream.good()) {
            std::string substr;
            std::getline(stream, substr, ',');
            if (substr != "")
                list.push_back(std::stod(substr));
        }
    }

    void parse_tree(boost::property_tree::ptree& tree, std::string path, std::vector<std::string>& list) {
        // Clear the list
        std::vector<std::string>().swap(list);

        std::stringstream stream(tree.get<std::string>(path));
        while (stream.good()) {
            std::string substr;
            std::getline(stream, substr, ',');
            if (substr != "")
                list.push_back(substr);
        }
    }

    double variance(std::vector<double> &samples)
    {
        int size = samples.size();

        double variance = 0;
        double t = samples[0];
        for (int i = 1; i < size; i++)
        {
            t += samples[i];
            double diff = ((i + 1) * samples[i]) - t;
            variance += (diff * diff) / ((i + 1.0) * i);
        }

        return variance / (size - 1);
    }

	double variance(std::vector<std::complex<double>> &samples)
	{
		int size = samples.size();

		double variance = 0;
		std::complex<double> t = samples[0];
		for (int i = 1; i < size; i++)
		{
			t += samples[i];
			std::complex<double> diff = ((i + 1.0) * samples[i]) - t;
			variance += std::real((diff * std::conj(diff)) / ((i + 1.0) * i));
		}

		return variance / (size - 1);
	}

	double standard_deviation(std::vector<double> &samples) {
		return std::sqrt(variance(samples));
	}

	double standard_deviation(std::vector<std::complex<double>> &samples) {
		return std::sqrt(variance(samples));
	}

	void permutations_with_repetition(std::vector<int> str, std::vector<int> prefix, const int lenght, std::vector<std::vector<int>> &permutations)
	{
		if (lenght == 1) {
			for (int j = 0; j < str.size(); j++) {
				std::vector<int> _permutation(prefix);
				_permutation.push_back(str[j]);
				permutations.push_back(_permutation);
			}
		} else {
			for (int i = 0; i < str.size(); i++) {
				std::vector<int> _permutation(prefix);
				_permutation.push_back(str[i]);
				permutations_with_repetition(str, _permutation, lenght - 1, permutations);
			}
		}
	}
}