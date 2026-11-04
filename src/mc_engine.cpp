#include "mc_engine.h"

#include "random.h"
#include "system.h"

#include <algorithm>
#include <cmath>

MonteCarloEngine::MonteCarloEngine(System& system)
    : sys(system), delta_theta(constants::pi / 20.0)
{
}

bool MonteCarloEngine::attempt_translate(int pid)
{
    const vec old_center = sys.particles[pid].get_center();
    const vec displacement(
        2.0 * (dis(gen) - 0.5) * delta_x,
        2.0 * (dis(gen) - 0.5) * delta_x);

    sys.set_center(pid, old_center + displacement);
    ++trial_x;
    if (!sys.overlap_unit_check(pid))
    {
        ++accept_x;
        return true;
    }

    // Restoring through System also restores the particle's cell membership.
    sys.set_center(pid, old_center);
    return false;
}

bool MonteCarloEngine::attempt_rotate(int pid)
{
    const double old_theta = sys.particles[pid].get_theta();
    const double new_theta = old_theta
        + 2.0 * (dis(gen) - 0.5) * delta_theta;

    sys.set_theta(pid, new_theta);
    ++trial_theta;
    if (!sys.overlap_unit_check(pid))
    {
        ++accept_theta;
        return true;
    }

    sys.set_theta(pid, old_theta);
    return false;
}

bool MonteCarloEngine::attempt_area()
{
    const double delta_ln = delta_ln_area * (dis(gen) - 0.5);
    const double scale = std::exp(-0.5 * delta_ln);
    const double trial_packing = sys.packing_fraction * scale * scale;
    ++area_trials;

    // This also guarantees a maximum particle length below sqrt(1/0.3) < 2,
    // which is the fixed cell-list cutoff used by System.
    if (trial_packing >= 1.0 - 1.0e-12) return false;

    sys.scale_particle_sizes(scale);
    const bool compression = delta_ln < 0.0;
    const bool invalid = compression && sys.overlap_sys_check();

    const double new_area = sys.V * std::exp(delta_ln);
    const double exponent = -sys.beta * sys.P * (new_area - sys.V)
        + static_cast<double>(sys.N + 1) * delta_ln;
    const double acceptance = std::exp(std::min(0.0, exponent));

    if (!invalid && dis(gen) < acceptance)
    {
        sys.V = new_area;
        sys.update_packing_fraction();
        ++area_accepts;
        return true;
    }

    sys.scale_particle_sizes(1.0 / scale);
    return false;
}

void MonteCarloEngine::adapt_particle_steps()
{
    const auto adapt = [](double& delta, int accepted, int attempted,
        double upper_limit)
    {
        if (attempted == 0) return;
        const double rate = static_cast<double>(accepted) / attempted;
        if (rate < 0.20) delta *= 0.5;
        else if (rate > 0.50) delta = std::min(2.0 * delta, upper_limit);
    };

    adapt(delta_x, accept_x, trial_x, 0.25 * sys.box_length);
    adapt(delta_theta, accept_theta, trial_theta, 0.5 * constants::pi);
    accept_x = trial_x = 0;
    accept_theta = trial_theta = 0;
}

bool MonteCarloEngine::run_nvt(int sweeps)
{
    if (sys.N <= 0 || sweeps <= 0) return true;
    change_step = std::max(1, 100 * sys.N);

    for (long long step = 0; step < static_cast<long long>(sweeps) * sys.N; ++step)
    {
        const int pid = sys.choose_particle();
        if (dis(gen) < 0.5) attempt_translate(pid);
        else attempt_rotate(pid);

        ++step_counter;
        if (step_counter % change_step == 0) adapt_particle_steps();
    }
    return true;
}

bool MonteCarloEngine::run_npt(int sweeps, double packing_checkpoint)
{
    if (sys.N <= 0 || sweeps <= 0) return true;
    change_step = std::max(1, 100 * sys.N);

    for (long long step = 0; step < static_cast<long long>(sweeps) * sys.N; ++step)
    {
        const int pid = sys.choose_particle();
        if (dis(gen) < 0.5) attempt_translate(pid);
        else attempt_rotate(pid);

        if (dis(gen) < 1.0 / static_cast<double>(sys.N + 1))
        {
            const bool accepted = attempt_area();
            if (accepted && sys.packing_fraction >= packing_checkpoint)
                return false;

            if (area_trials >= 100)
            {
                const double rate = static_cast<double>(area_accepts) / area_trials;
                if (rate < 0.20) delta_ln_area *= 0.5;
                else if (rate > 0.50) delta_ln_area *= 2.0;
                area_trials = 0;
                area_accepts = 0;
            }
        }

        ++step_counter;
        if (step_counter % change_step == 0) adapt_particle_steps();
    }
    return true;
}
