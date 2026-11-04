#include "system.h"

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

using Config = map<string, string>;

Config load_config(const string& filename)
{
    ifstream file(filename);
    if (!file)
    {
        throw runtime_error("Cannot open " + filename);
    }

    Config config;
    string line;
    while (getline(file, line))
    {
        const size_t comment = line.find('#');
        if (comment != string::npos)
        {
            line.erase(comment);
        }

        istringstream input(line);
        string key;
        string equal;
        string value;

        if (!(input >> key))
        {
            continue;
        }

        if (!(input >> equal >> value) || equal != "=")
        {
            throw runtime_error("Invalid config line: " + line);
        }

        config[key] = value;
    }

    return config;
}

int main()
{
    const clock_t start = clock();

    try
    {
        const Config config = load_config("input/config.txt");

        const std::filesystem::path output_dir = "output";
        std::filesystem::create_directories(output_dir);
        std::filesystem::current_path(output_dir);
        const double pressure_end = stod(config.at("pressure_end"));
        const double packing_fraction_end =
            stod(config.at("packing_fraction_end"));
        const int mc_sweeps = stoi(config.at("mc_sweeps"));

        System simulation(config);
        simulation.set_config_random();
        simulation.mc_npt_annealing(
            pressure_end, packing_fraction_end, mc_sweeps);
    }
    catch (const exception& error)
    {
        cerr << "Simulation error: " << error.what() << '\n';
        return 1;
    }

    cout << "NPT simulation used "
        << static_cast<double>(clock() - start) / CLOCKS_PER_SEC / 60.0
        << " minutes.\n";
    return 0;
}
