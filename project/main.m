%% MAIN.M
% Single entry point for the PEMFC membrane degradation simulator.
% Run this file from the project root (no manual edits required) to
% execute the master closed-loop lifetime simulation and, optionally, the
% supporting mechanism-isolation / sensitivity / diagnostic analyses.
%
% USAGE:
%   >> main
%
% This will:
%   1. Run the closed-loop lifetime_simulation.m for p.default_hours hours
%      (edit the value below to change simulation length; 5000 h is the
%      project's target duration, but note this takes considerable wall-
%      clock time - see VALIDATION_REPORT.md for benchmarked runtimes and
%      the ode15s/ode45 solver note in lifetime_simulation.m).
%   2. Save the resulting history struct to sim_history.mat.
%   3. Generate the core diagnostic figures (voltage, conductivity, HF,
%      ECSA, Pt, crossover, H2O2, side-chain concentration) via
%      generate_presentation_data.m.
%
% The secondary analysis scripts (species_evolution.m, time_scale_analysis.m,
% mechanism_ranking.m, parameter_sensitivity.m, interaction_analysis.m,
% heatmap_analysis.m) can be run independently; see README.md for their
% individual calling conventions and current verification status.

clear; clc;
addpath(pwd);

fprintf('=====================================================\n');
fprintf(' PEMFC Membrane Degradation Simulator - main.m\n');
fprintf('=====================================================\n\n');

%% Configuration
simulation_hours = 120;   % EDIT: set to 5000 for the full target lifetime
                          % (see VALIDATION_REPORT.md for runtime guidance)

%% 1. Run master closed-loop lifetime simulation
fprintf('[1/2] Running lifetime_simulation.m for %d hours...\n', simulation_hours);
history = lifetime_simulation(simulation_hours);
save('sim_history.mat', 'history');
fprintf('      Saved sim_history.mat\n\n');

%% 2. Generate presentation/diagnostic figures
fprintf('[2/2] Generating figures via generate_presentation_data.m...\n');
try
    generate_presentation_data();
    fprintf('      Figures generated in presentation_exports/.\n');
catch err
    fprintf('      NOTE: generate_presentation_data.m raised: %s\n', err.message);
    fprintf('      (Core simulation results in sim_history.mat are unaffected.)\n');
end

fprintf('\n=====================================================\n');
fprintf(' Done. See sim_history.mat for full results.\n');
fprintf(' See README.md and VALIDATION_REPORT.md for details.\n');
fprintf('=====================================================\n');
