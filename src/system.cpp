#include "system.h"

#include "hard_overlap.h"
#include "mc_engine.h"
#include "random.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
int required_int(const std::map<std::string, std::string>& config,
    const std::string& key)
{
    const auto item = config.find(key);
    if (item == config.end()) throw std::runtime_error("Missing parameter: " + key);
    return std::stoi(item->second);
}

double required_double(const std::map<std::string, std::string>& config,
    const std::string& key)
{
    const auto item = config.find(key);
    if (item == config.end()) throw std::runtime_error("Missing parameter: " + key);
    return std::stod(item->second);
}

}

System::System()
    : mc_engine(std::make_unique<MonteCarloEngine>(*this))
{
}

System::System(const std::map<std::string, std::string>& config)
    : N(required_int(config, "N")),
      aspect_ratio(required_double(config, "aspect_ratio")),
      mc_engine(std::make_unique<MonteCarloEngine>(*this))
{
    if (N <= 0) throw std::runtime_error("N must be positive.");
    if (aspect_ratio < 1.0) throw std::runtime_error("aspect_ratio must be at least 1.");

    const auto boundary = config.find("boundary_mode");
    if (boundary == config.end() || boundary->second == "periodic")
        boundary_mode = BoundaryMode::PERIODIC;
    else if (boundary->second == "closed")
        boundary_mode = BoundaryMode::CLOSED;
    else
        throw std::runtime_error("boundary_mode must be periodic or closed.");

    const double initial_length = 1.0;
    const double initial_width = 1.0 / aspect_ratio;
    particles.assign(N, ParticleType(vec(), 0.0, initial_length, initial_width));

    const double initial_packing_fraction = 0.3;
    box_length = std::sqrt(N * particles.front().get_area()
        / initial_packing_fraction);
    V = box_length * box_length;
    packing_fraction = initial_packing_fraction;

    // At eta <= 1 the longest particle is sqrt(1/0.3) < 2. A periodic
    // reference box shorter than 2 could therefore overlap its own image.
    if (boundary_mode == BoundaryMode::PERIODIC && box_length <= 2.0)
    {
        throw std::runtime_error(
            "Periodic box is too small for the maximum allowed particle length.");
    }
}

System::~System() = default;

void System::set_config_random()
{
    // Initial particle-particle overlaps are intentionally allowed. At the
    // initial packing fraction 0.3, the following NVT stage can disperse the
    // random configuration before pressure annealing begins.
    double available_length = box_length;
    double coordinate_min = -0.5 * box_length;
    if (boundary_mode == BoundaryMode::CLOSED)
    {
        // Avoid trapping a particle across a hard wall. The half-length is a
        // safe orientation-independent wall margin for both particle models.
        const double margin = particles.front().get_bounding_radius();
        available_length -= 2.0 * margin;
        coordinate_min += margin;
        if (available_length <= 0.0)
            throw std::runtime_error("Closed box is too small for one particle.");
    }

    for (int pid = 0; pid < N; ++pid)
    {
        const vec center(
            coordinate_min + available_length * dis(gen),
            coordinate_min + available_length * dis(gen));
        particles[pid].translate(center);
        particles[pid].rotate(2.0 * constants::pi * dis(gen));
    }

    cell_list_on();
}

void System::cell_list_on()
{
    // Keep an equal tiling of the fixed reference box. floor() makes the
    // actual cell width at least 2 rather than leaving a narrow final cell.
    cell_number = std::max(1, static_cast<int>(std::floor(box_length / 2.0)));

    cell_width = box_length / cell_number;
    rebuild_cell_list();
}

std::pair<int, int> System::locate_cell(const vec& center) const
{
    int cx = static_cast<int>(std::floor((center.x + 0.5 * box_length) / cell_width));
    int cy = static_cast<int>(std::floor((center.y + 0.5 * box_length) / cell_width));
    cx = std::clamp(cx, 0, cell_number - 1);
    cy = std::clamp(cy, 0, cell_number - 1);
    return { cx, cy };
}

bool System::rebuild_cell_list()
{
    head_list.assign(cell_number, std::vector<int>(cell_number, -1));
    particle_list.assign(N, ChainPointer());
    for (int pid = 0; pid < N; ++pid)
    {
        const auto [cx, cy] = locate_cell(particles[pid].get_center());
        const int old_head = head_list[cx][cy];
        particle_list[pid].cx = cx;
        particle_list[pid].cy = cy;
        particle_list[pid].pre = -1;
        particle_list[pid].next = old_head;
        if (old_head != -1) particle_list[old_head].pre = pid;
        head_list[cx][cy] = pid;
    }
    return true;
}

bool System::update_cell(int pid)
{
    const auto [new_x, new_y] = locate_cell(particles[pid].get_center());
    ChainPointer& node = particle_list[pid];
    if (node.cx == new_x && node.cy == new_y) return false;

    if (node.pre == -1)
        head_list[node.cx][node.cy] = node.next;
    else
        particle_list[node.pre].next = node.next;
    if (node.next != -1) particle_list[node.next].pre = node.pre;

    const int old_head = head_list[new_x][new_y];
    node.pre = -1;
    node.next = old_head;
    node.cx = new_x;
    node.cy = new_y;
    if (old_head != -1) particle_list[old_head].pre = pid;
    head_list[new_x][new_y] = pid;
    return true;
}

int System::periodic_cell_index(int index) const
{
    if (index < 0) return index + cell_number;
    if (index >= cell_number) return index - cell_number;
    return index;
}

int System::collect_neighbor_cells(int cx, int cy,
    std::array<std::pair<int, int>, 9>& cells) const
{
    int count = 0;
    for (int dx = -1; dx <= 1; ++dx)
    {
        int nx = cx + dx;
        if (boundary_mode == BoundaryMode::PERIODIC)
            nx = periodic_cell_index(nx);
        else if (nx < 0 || nx >= cell_number)
            continue;

        for (int dy = -1; dy <= 1; ++dy)
        {
            int ny = cy + dy;
            if (boundary_mode == BoundaryMode::PERIODIC)
                ny = periodic_cell_index(ny);
            else if (ny < 0 || ny >= cell_number)
                continue;

            const std::pair<int, int> cell(nx, ny);
            bool duplicate = false;
            for (int index = 0; index < count; ++index)
            {
                if (cells[index] == cell)
                {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) cells[count++] = cell;
        }
    }
    return count;
}

bool System::overlap_unit_cell(int pid) const
{
    if (boundary_mode == BoundaryMode::CLOSED
        && hard_overlap::boundary_overlap(particles[pid], box_length))
    {
        return true;
    }

    const int cx = particle_list[pid].cx;
    const int cy = particle_list[pid].cy;
    std::array<std::pair<int, int>, 9> cells;
    const int cell_count = collect_neighbor_cells(cx, cy, cells);
    for (int index = 0; index < cell_count; ++index)
    {
        const auto [nx, ny] = cells[index];
        int other = head_list[nx][ny];
        while (other != -1)
        {
            if (other != pid)
            {
                const bool does_overlap = (boundary_mode == BoundaryMode::PERIODIC)
                    ? hard_overlap::overlap_periodic(
                        particles[pid], particles[other], box_length)
                    : hard_overlap::overlap(particles[pid], particles[other]);
                if (does_overlap) return true;
            }
            other = particle_list[other].next;
        }
    }
    return false;
}

bool System::overlap_sys_cell() const
{
    for (int pid = 0; pid < N; ++pid)
    {
        if (boundary_mode == BoundaryMode::CLOSED
            && hard_overlap::boundary_overlap(particles[pid], box_length))
        {
            return true;
        }

        const int cx = particle_list[pid].cx;
        const int cy = particle_list[pid].cy;
        std::array<std::pair<int, int>, 9> cells;
        const int cell_count = collect_neighbor_cells(cx, cy, cells);
        for (int index = 0; index < cell_count; ++index)
        {
            const auto [nx, ny] = cells[index];
            int other = head_list[nx][ny];
            while (other != -1)
            {
                // Each unordered particle pair is tested exactly once.
                if (other > pid)
                {
                    const bool does_overlap = (boundary_mode == BoundaryMode::PERIODIC)
                        ? hard_overlap::overlap_periodic(
                            particles[pid], particles[other], box_length)
                        : hard_overlap::overlap(particles[pid], particles[other]);
                    if (does_overlap) return true;
                }
                other = particle_list[other].next;
            }
        }
    }
    return false;
}

bool System::overlap_unit_check(int pid) const
{
    return overlap_unit_cell(pid);
}

bool System::overlap_sys_check() const
{
    return overlap_sys_cell();
}

double System::wrap_coordinate(double value) const
{
    double wrapped = std::fmod(value + 0.5 * box_length, box_length);
    if (wrapped < 0.0) wrapped += box_length;
    return wrapped - 0.5 * box_length;
}

void System::set_center(int pid, const vec& new_center)
{
    vec center = new_center;
    if (boundary_mode == BoundaryMode::PERIODIC)
    {
        center.x = wrap_coordinate(center.x);
        center.y = wrap_coordinate(center.y);
    }
    particles[pid].translate(center);
    update_cell(pid);
}

void System::set_theta(int pid, double new_theta)
{
    new_theta = std::fmod(new_theta, 2.0 * constants::pi);
    if (new_theta < 0.0) new_theta += 2.0 * constants::pi;
    particles[pid].rotate(new_theta);
}

void System::scale_particle_sizes(double factor)
{
    for (ParticleType& particle : particles) particle.scale_size(factor);
}

int System::choose_particle() const
{
    return std::min(N - 1, static_cast<int>(dis(gen) * N));
}

void System::update_packing_fraction()
{
    packing_fraction = N * particles.front().get_area()
        / (box_length * box_length);
}

void System::mc_npt(double pressure, int sweeps)
{
    P = pressure;
    mc_engine->run_npt(sweeps, std::numeric_limits<double>::infinity());
}

void System::mc_npt_annealing(double pressure_end,
    double packing_fraction_end, int sweeps_per_stage)
{
    if (pressure_end <= 0.0)
        throw std::runtime_error("pressure_end must be positive.");
    if (packing_fraction_end <= 0.0 || packing_fraction_end >= 1.0)
        throw std::runtime_error(
            "packing_fraction_end must be between 0 and 1.");
    if (sweeps_per_stage <= 0)
        throw std::runtime_error("mc_sweeps must be positive.");

    double ln_pressure = 0.5;
    P = std::exp(ln_pressure);
    std::ofstream structure_dump("structure.dump");
    std::ofstream pressure_dump("P_out.dump");
    std::ofstream pressure_packing("p-pkf.dump");
    if (!structure_dump || !pressure_dump || !pressure_packing)
        throw std::runtime_error("Could not open simulation output files.");

    int structure_step = 0;
    int pressure_step = 0;
    mc_engine->run_nvt(sweeps_per_stage);

    double packing_checkpoint = std::min(
        packing_fraction + 0.01, packing_fraction_end);

    while (packing_fraction < packing_fraction_end && P <= pressure_end)
    {
        const bool completed = mc_engine->run_npt(
            sweeps_per_stage, packing_checkpoint);
        if (!completed)
        {
            mc_engine->run_nvt(sweeps_per_stage);
            pressure_packing << P << " " << packing_fraction << "\n";
            output(structure_dump, structure_step++);
            packing_checkpoint = std::min(
                packing_checkpoint + 0.01, packing_fraction_end);
        }
        else
        {
            std::cout << "Packing fraction = " << packing_fraction
                << ", P = " << P << "\n";

            // Save the equilibrated state whenever the next pressure level
            // is about to be applied. P_out.dump is a LAMMPS snapshot file.
            output(pressure_dump, pressure_step++);
            ln_pressure += 0.1;
            P = std::exp(ln_pressure);
        }
    }

    if (packing_fraction >= packing_fraction_end)
    {
        std::cout << "Target packing fraction reached: "
            << packing_fraction << "\n";
    }
    else
    {
        std::cout << "Pressure safety limit reached before target packing "
            << "fraction: P_end = " << pressure_end
            << ", packing fraction = " << packing_fraction << "\n";
    }
    output(structure_dump, structure_step);
}

void System::output(std::ofstream& structure, int time_step) const
{
    structure << "ITEM: TIMESTEP\n" << time_step << "\n";
    structure << "ITEM: NUMBER OF ATOMS\n" << N << "\n";
    const char* flag = (boundary_mode == BoundaryMode::PERIODIC) ? "pp" : "ff";
    structure << "ITEM: BOX BOUNDS " << flag << " " << flag << " ff\n";
    structure << -0.5 * box_length << " " << 0.5 * box_length << "\n";
    structure << -0.5 * box_length << " " << 0.5 * box_length << "\n";
    structure << "-0.5 0.5\n";
    structure << "ITEM: ATOMS id type x y z theta radius "
        << "quati quatj quatk quatw valid\n";
    for (int pid = 0; pid < N; ++pid)
    {
        const ParticleType& particle = particles[pid];
        const double half_theta = 0.5 * particle.get_theta();
        structure << pid << " 1 " << particle.get_center().x << " "
            << particle.get_center().y << " 0 " << particle.get_theta() << " "
            << particle.get_length() << " 0 0 " << std::sin(half_theta) << " "
            << std::cos(half_theta) << " "
            << (overlap_unit_check(pid) ? 0 : 1) << "\n";
    }
}
