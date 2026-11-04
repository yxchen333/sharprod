#pragma once

class System;

class MonteCarloEngine
{
public:
    explicit MonteCarloEngine(System& system);

    bool run_nvt(int sweeps);
    // Returns false when an accepted area move reaches packing_checkpoint.
    bool run_npt(int sweeps, double packing_checkpoint);

private:
    bool attempt_translate(int pid);
    bool attempt_rotate(int pid);
    bool attempt_area();
    void adapt_particle_steps();

    System& sys;

    long long step_counter = 0;
    int change_step = 1;
    int accept_x = 0;
    int trial_x = 0;
    int accept_theta = 0;
    int trial_theta = 0;
    int area_accepts = 0;
    int area_trials = 0;

    double delta_x = 0.1;
    double delta_theta = 0.0;
    double delta_ln_area = 0.001;
};
