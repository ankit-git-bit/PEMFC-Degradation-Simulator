# PEMFC Membrane Degradation Simulator

A physics-based, modular MATLAB/Octave simulator for PEM fuel cell membrane
chemical degradation, coupling five mechanisms into a closed feedback loop:
peroxide generation (from H2/O2 crossover) -> Fenton chemistry -> Pt
dissolution/redeposition -> membrane side-chain attack -> backbone
unzipping -> HF release -> IEC/conductivity/thickness loss -> increased
hydrogen crossover -> increased peroxide generation (closes the loop).

**Read `VALIDATION_REPORT.md` before using simulation outputs for anything
beyond qualitative trend illustration.** It documents every bug found and
fixed, every parameter that is a literature value vs. an engineering
assumption, and the current verification status of every file.

## Quick start

```matlab
main
```

Runs a 120-hour demonstration simulation (edit `simulation_hours` in
`main.m` to run longer, e.g. 5000 h for the full target lifetime -- see
the runtime note below) and saves `sim_history.mat`.

## Execution order / module map

```
parameters.m              -- all physical/initial-condition constants (p)
kinetic_parameters.m      -- all rate constants (kinetics)
operating_conditions.m    -- T, P, RH, current density, etc. (op)
mechanism_switch.m        -- on/off flags for each mechanism

reaction_rates.m          -- R1-R23 Arrhenius rates + QSSA radical solve
  |-- peroxide_generation.m   (crossover -> H2O2 source term)
  |     `-- hydrogen_crossover.m
  `-- pt_dissolution.m        (Pt kinetics, single source of truth)

ode_system.m               -- 21-state chemistry ODE (maps rates -> dydt)

lifetime_simulation.m      -- MASTER DRIVER. For each simulated hour:
  |   1. integrate ode_system.m (fast chemistry)
  |   2. degradation_update.m   (SOLE authority for IEC, sigma, L_mem,
  |                               crossover_H2, ECSA)
  |   3. bidirectional_coupling.m (feeds crossover back into op.H2/op.O2
  |                                 for the next hour -- closes the loop)
  `-- 4. performance_update.m   (-> voltage, power, efficiency)

membrane_properties.m      -- diagnostics only (water uptake, porosity,
                               mechanical/thermal estimates); reuses
                               degradation_update.m's authoritative values,
                               does not recompute them independently

main.m                     -- single entry point, calls the above
```

## Known-working vs. not-fully-verified components

**Fully debugged, tested, and verified in this engagement** (see
VALIDATION_REPORT.md for the full bug list): `parameters.m`,
`kinetic_parameters.m`, `reaction_rates.m`, `peroxide_generation.m`,
`hydrogen_crossover.m`, `pt_dissolution.m`, `ode_system.m`,
`degradation_update.m`, `bidirectional_coupling.m`, `performance_update.m`,
`membrane_properties.m`, `lifetime_simulation.m`, `main.m`.

**Present in the project but NOT independently re-verified in this pass**
(they use their own self-contained, simplified local ODE implementations,
separate from the master pipeline above -- a design already present before
this engagement): `species_evolution.m`, `time_scale_analysis.m`,
`mechanism_ranking.m`, `parameter_sensitivity.m`, `interaction_analysis.m`,
`heatmap_analysis.m`, `generate_presentation_data.m`,
`get_audited_parameters.m`. `interaction_analysis.m` uses MATLAB's
`digraph`/`plot` graph objects and `mechanism_ranking.m` uses MATLAB's
`table` -- both unavailable in GNU Octave, so these two could not be
executed at all in the Octave-based test environment used for this
engagement and should be checked directly in MATLAB. See
VALIDATION_REPORT.md, Section 5, for specifics and recommended fixes.

## Runtime note (important)

The reaction network (Fenton radical chemistry + backbone unzipping) is
numerically stiff: hydroxyl-radical timescales are nanosecond-scale while
the simulation spans thousands of hours. `reaction_rates.m` uses a
quasi-steady-state approximation (QSSA) for OH/OOH/H to make this
tractable (see code comments there for the numerical-stability rationale).

`lifetime_simulation.m` tries `ode15s` (fast, MATLAB's variable-order BDF
solver) first and falls back to `ode45` if it fails to converge. In the
GNU Octave environment used to develop and test this project, `ode15s`
(backed by SUNDIALS/IDA in Octave) consistently failed to converge on this
system, and the `ode45` fallback was used for all validation runs
reported here (roughly 3-9 seconds of wall-clock time per simulated hour
in Octave). **In MATLAB proper, `ode15s` is expected to succeed and should
be substantially faster** -- this was not verified against a real MATLAB
license in this engagement and should be checked by the user. A full 5000 h
run took longer than the compute budget available in this engagement; a
120 h demonstration dataset (`checkpoint.mat`, `history_export.csv`,
`figures/`) is included instead. See VALIDATION_REPORT.md for the observed
120 h trends and guidance on extrapolating/extending the run.

## Requirements

MATLAB (R2018b+ recommended, for `ode15s`, `digraph`, `table`) or GNU
Octave 8+ with its bundled `ode45`/`ode15s` solvers. No toolboxes beyond
base MATLAB/Octave are required for the master pipeline.

## Literature basis

The Fenton/radical/iron-redox rate constants (R1-R23 pre-exponential
factors and activation energies in `kinetic_parameters.m`, Section 1) are
taken from Fruhwirt, P. et al., "Holistic approach to chemical degradation
of Nafion membranes in fuel cells: modelling and predictions," *Phys.
Chem. Chem. Phys.* 22 (2020) 5647-5666, and were verified for correct
L/(mol*s) -> m^3/(mol*s) unit conversion (cross-checked against the
self-contained audit in `get_audited_parameters.m`). Pt dissolution
functional form follows Darling & Meyers (2003) and Rinaldo et al. (2011);
specific rate constants for Pt, coupling/percolation, crossover, and
performance parameters are engineering order-of-magnitude estimates, not
individually fit to a published dataset for this project -- every such
parameter is flagged `TODO/ASSUMPTION` with its rationale directly in the
source code (`parameters.m` and `kinetic_parameters.m`). See
VALIDATION_REPORT.md and the accompanying LaTeX report for the complete
reference list and an explicit inventory of literature-sourced vs.
assumed parameters.
