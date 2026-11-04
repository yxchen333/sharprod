#pragma once

#include "particle.h"

#include <array>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

class MonteCarloEngine;

class System
{
public:
    enum class BoundaryMode
    {
        CLOSED,
        PERIODIC
    };

    // Change only this alias to switch the simulated particle model.
    //using ParticleType = SharpRodParticle;
    using ParticleType = EllipseParticle;

    System();
    explicit System(const std::map<std::string, std::string>& config);
    ~System();

    void set_config_random();
    void cell_list_on();

    void mc_npt(double pressure, int sweeps);
    void mc_npt_annealing(double pressure_end,
        double packing_fraction_end, int sweeps_per_stage);

    bool overlap_unit_check(int pid) const;
    bool overlap_sys_check() const;

    void set_center(int pid, const vec& new_center);
    void set_theta(int pid, double new_theta);
    void scale_particle_sizes(double factor);
    int choose_particle() const;

    void output(std::ofstream& structure, int time_step) const;

private:
    friend class MonteCarloEngine;

    struct ChainPointer
    {
        int pre = -1;
        int next = -1;
        int cx = 0;
        int cy = 0;
    };

    bool overlap_unit_cell(int pid) const;
    bool overlap_sys_cell() const;

    bool update_cell(int pid);
    bool rebuild_cell_list();
    int collect_neighbor_cells(int cx, int cy,
        std::array<std::pair<int, int>, 9>& cells) const;
    std::pair<int, int> locate_cell(const vec& center) const;
    int periodic_cell_index(int index) const;
    double wrap_coordinate(double value) const;
    void update_packing_fraction();

    int N = 0;
    std::vector<ParticleType> particles;

    double aspect_ratio = 1.0;
    double box_length = 0.0;
    double P = 0.0;
    double beta = 1.0;
    double V = 0.0; // Effective NPT area; the reference box stays fixed.
    double packing_fraction = 0.0;

    BoundaryMode boundary_mode = BoundaryMode::PERIODIC;

    int cell_number = 0;
    double cell_width = 0.0;
    std::vector<std::vector<int>> head_list;
    std::vector<ChainPointer> particle_list;

    std::unique_ptr<MonteCarloEngine> mc_engine;
};
