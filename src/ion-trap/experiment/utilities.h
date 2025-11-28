#include <cmath>
#include <vector>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include <cnpy/cnpy.h>
#include <complex>

namespace experiment {

    // Utility methods
    void load_1d_list(std::string filename, std::vector<double> &_list);
    void load_2d_list(std::string filename, std::vector<std::vector<double>> &_list);

    void parse_tree(boost::property_tree::ptree& tree, std::string path, std::vector<int>& list);
    void parse_tree(boost::property_tree::ptree& tree, std::string path, std::vector<double>& list);
    void parse_tree(boost::property_tree::ptree& tree, std::string path, std::vector<std::string>& list);

	double variance(std::vector<double> &samples);
	double variance(std::vector<std::complex<double>> &samples);

	double standard_deviation(std::vector<double> &samples);
	double standard_deviation(std::vector<std::complex<double>> &samples);

	void permutations_with_repetition(std::vector<int> str, std::vector<int> prefix, const int lenght, std::vector<std::vector<int>> &permutations);
}