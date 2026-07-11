# Validation Report

This project was received as a partially-built, non-running MATLAB codebase
(new files `fenton_chemistry.m`, `pt_dissolution.m`, `bidirectional_coupling.m`,
`hydrogen_crossover.m`, `degradation_update.m`, `performance_update.m` were
already present, but wired together incorrectly / not wired together at
all). This report documents the audit process, every bug found and fixed,
what remains a documented engineering assumption rather than a validated
literature value, and the current status of every file. All testing was
performed in GNU Octave 8/9 (no MATLAB license was available in this
engagement); see the runtime note in README.md.

## 1. Summary of what was fixed

| # | File | Bug | Fix |
|---|------|-----|-----|
| 1 | `reaction_rates.m` | Called `peroxide_generation(species,props,op,p)` (4 args, order `species,props,op,p`) but the function was defined as `peroxide_generation(rates,species,props,op,p)` (5 args, different order/type) -- a hard argument-count/order mismatch that crashed on first call. | Rewrote `peroxide_generation.m` as a dual-mode function: 4-arg "flux mode" (used by `reaction_rates.m`) and 5-arg "aggregator mode" (kept for standalone Mechanism-1 diagnostics), sharing one crossover-flux implementation. |
| 2 | `hydrogen_crossover.m` | `J.H2 = (P_H2*delta_P)/(props.L_mem*1e-6)` -- `props.L_mem` is already stored in SI meters (e.g. `50e-6`), so the extra `*1e-6` double-converted it, inflating crossover flux by ~1e6x (5549 mol/(m^2 s) instead of ~0.006). This alone made the coupled ODE numerically unsolvable. | Removed the extra `*1e-6`; documented the correct SI-unit contract for `props.L_mem` in the code. |
| 3 | `parameters.m` | `p.Fe2_0 = 89.5 mol/m^3` for a stated "5 ppm" Fe contamination, via `p.ppm_to_mol_m3 = 17.9`. A defensible ppm-by-mass-to-molarity conversion gives ~0.09-0.18 mol/m^3 -- the original value was ~500x too concentrated, and produced complete side-chain sulfonic-acid depletion within *seconds* of simulated time instead of the literature-expected thousands-of-hours timescale. | Corrected `ppm_to_mol_m3` to 0.0179 mol/m^3 per ppm (1 ppm = 1 g Fe/m^3 reference solution / 55.845 g/mol), giving `Fe2_0 = 0.0895 mol/m^3` for 5 ppm -- cross-checked against a Nafion-mass-density-based estimate (0.18 mol/m^3), same order of magnitude. |
| 4 | `reaction_rates.m` (new QSSA code) | The quadratic solver used for the OH/OOH quasi-steady-state balance, `x=(-b+sqrt(b^2+4ac))/(2a)`, suffers catastrophic floating-point cancellation whenever `b^2 >> 4ac` (common case here: `b` dominated by e.g. `k14*SC_SO3H ~ 1e9`, while `c` is many orders of magnitude smaller). `sqrt(b^2+4ac)` then rounds to exactly `b` in double precision, so the subtraction returns exactly 0.0, which silently froze all membrane-attack/unzipping chemistry after roughly the first simulated hour (verified: state vector bit-identical to 14 decimal places across multiple simulated hours). | Rewrote using the numerically stable "Citardauq" form `x = 2c/(b+sqrt(b^2+4ac))`, which avoids subtracting nearly-equal large numbers. Verified: `OH_qss` went from exactly `0.0` to `5.3e-14` at the previously-frozen state, physically reasonable. |
| 5 | `reaction_rates.m` / `ode_system.m` | OH, OOH, H radicals were modeled as explicit ODE states alongside hour-timescale species. Their pseudo-first-order sink rates (~1e9 s^-1, sub-microsecond lifetimes) versus the hour-scale macro-timestep created a stiffness ratio >1e12, which no solver (Octave's `ode15s`/IDA or `ode45`) could resolve at the initial condition. | Implemented a quasi-steady-state approximation (QSSA) for OH/OOH/H: solved algebraically (closed-form, single deterministic pass -- no iterative loop, to keep the RHS smooth for the solver's numerical Jacobian) at every rate evaluation; their ODE-state slots carry `dydt=0`. Standard practice in stiff radical-chemistry modeling (Turanyi & Tomlin, *Analysis of Kinetic Reaction Mechanisms*, Springer 2014). |
| 6 | `peroxide_generation.m` | Even after fixing bug #2, `k_gen = 0.5` produced a steady-state membrane H2O2 concentration of ~600 mol/m^3 (0.6 mol/L), roughly 3-6 orders of magnitude above literature-reported PEMFC membrane H2O2 concentrations (typically micromolar), driving unrealistic OH levels and multi-second side-chain depletion. | Recalibrated `k_gen` to `5e-6`, targeting a literature-consistent micromolar steady-state H2O2 concentration. Flagged ASSUMPTION/TODO pending direct experimental calibration. |
| 7 | `kinetic_parameters.m` | `kinetics.Pt.k_diss_A = 1e-12` gave a Pt-loss e-folding time of ~11 hours (near-total metallic Pt depletion within about a day) -- inconsistent with any realistic PEMFC durability timescale. | Recalibrated to `2e-16`, giving an e-folding time of ~54,600 h (~9% Pt loss over a 5000 h simulated lifetime), consistent with the order of magnitude of ECSA loss reported in durability studies (e.g. Ohma et al., ECS Trans. 41 (2011) 775). Flagged ASSUMPTION/TODO. |
| 8 | `performance_update.m` | Compared a per-real-catalyst-area exchange current density (`p.i0_ref ~ 1e-4 A/m^2_Pt`) directly against the geometric operating current density (`op.current_density ~ 1000 A/m^2`) without correcting for catalyst-layer roughness factor, giving an activation loss exceeding V_ocv, so voltage was clamped to exactly 0 V for every condition. | Added the standard roughness-factor correction (ECSA x Pt areal loading, ~240 m^2_Pt/m^2_geometric here), giving a plausible non-zero cell voltage (~0.28 V; still on the low side for a fully calibrated polarization curve -- flagged for further tuning). |
| 9 | `performance_update.m` | Read `op.R, op.F, op.alpha, op.ECSA_0, op.i0_ref, op.V_ocv` -- none defined in `operating_conditions.m` (they live in `parameters.m`/`kinetic_parameters.m`). | Rewrote to source all constants from `p` consistently. |
| 10 | `degradation_update.m` | Referenced `species.CF2_backbone` (never constructed anywhere -- the ODE tracks `CF2_7`...`CF2_1` individually) and `p.CF2_0` (undefined; correct field is `p.CF2_total_0`). | Fixed to sum the seven tracked `CF2_7..CF2_1` pools plus the untracked fixed `p.CF2_backbone_0` reservoir, matching the actual state-vector layout and initial mass balance. |
| 11 | `ode_system.m` | Independently recomputed IEC/sigma/L_mem/ECSA/crossover_H2 inside the ODE derivative function, duplicating and conflicting with `degradation_update.m` -- directly violating the project's own module contract ("degradation_update.m should be the ONLY place these are updated"). | Removed the duplicate property-evolution block; state vector trimmed to the 21 genuinely fast/chemistry-timescale states; all physical-property evolution now happens exclusively in `degradation_update.m`, called once per macro-timestep. |
| 12 | `pt_dissolution.m`, `bidirectional_coupling.m`, `membrane_properties.m` | All three files were orphaned -- fully written but never called by any other file (confirmed via static call-site search). `pt_dissolution.m` additionally referenced parameters that do not exist anywhere in `parameters.m`. | `pt_dissolution.m` rewritten to use the actual `kinetics.Pt` struct and wired in as the single source of truth for Pt kinetics (called from `reaction_rates.m`). `bidirectional_coupling.m` rewritten and wired in to feed the updated crossover state back into `op.H2`/`op.O2` for the next macro-timestep, closing the feedback loop. `membrane_properties.m` rewritten as diagnostics-only, reusing `degradation_update.m`'s authoritative values instead of recomputing them independently. |
| 13 | `lifetime_simulation.m` | Referenced undefined variables `pt_state` and `initial_state` (leftover from an incomplete prior edit) and called a local 18-state helper instead of the actual `ode_system.m`/`degradation_update.m`/`performance_update.m` modules -- immediate crash. | Rewritten from scratch as the master driver described in README.md. |

## 2. Numerical-methods notes (read before extending the simulation)

- **QSSA for OH/OOH/H** (bug #5): a genuine simplification made necessary
  by the reaction network's intrinsic stiffness. It does not change the
  reaction network, rate constants, or stoichiometry -- only how the fast
  radical concentrations are numerically obtained. A future revision
  needing fully time-resolved sub-second radical transients would need a
  proper stiff DAE solver with an analytically supplied, well-scaled
  Jacobian instead.
- **Solver choice**: `ode15s` failed to converge in Octave for this system
  even after the fixes above (observed failure: `IDASolve` error-test
  failures at vanishingly small step sizes). The `ode45` fallback is
  empirically robust (used for all data in this report) but slower
  (~3-9 s wall-clock per simulated hour in Octave). This may be an
  Octave/SUNDIALS-specific limitation; re-testing `ode15s` directly in
  MATLAB is recommended.

## 3. Observed 120-hour demonstration run

A 120-hour run (`checkpoint.mat`, `history_export.csv`, `figures/`) is
included as a demonstration/validation dataset (a full 5000 h run exceeded
the compute budget available in this engagement -- see README.md).
Observed behavior:

- Voltage, conductivity, IEC, ECSA, and Pt loading all decline
  monotonically and smoothly, as expected.
- An interesting emergent behavior: Fe2+ is consumed rapidly within the
  first simulated hour (converted to Fe3+ via R1), and its regeneration
  (via the much slower R2) is not fast enough to sustain the initial
  radical-generation rate. The result is a fast initial "burst" of
  side-chain/HF-generating chemistry within the first hour, followed by a
  much slower "iron-depletion-limited" quasi-steady regime for the rest of
  the 120 h window (HF and SC-SO3H become nearly constant after ~hour 30).
  **This needs longer-run verification**: real membranes are reported to
  keep degrading steadily over thousands of hours, so either (a) this
  iron-limited stall is a genuine feature of Fenton-catalyzed degradation
  once local iron is consumed, requiring slow Fe3+ -> Fe2+ regeneration or
  direct H2O2 auto-decomposition (R6) to sustain long-run degradation, or
  (b) the R2 rate constant / iron-transport assumptions understate the
  long-run iron-cycling rate. A full 5000 h run is the recommended next
  validation step to determine which regime dominates.
- Absolute voltage (~0.28 V at 1000 A/m^2) is on the low side for a
  well-tuned PEMFC polarization curve (typically 0.6-0.7 V in that current
  range); `performance_update.m`'s single-Tafel-region model would benefit
  from calibration against a real polarization curve.

## 4. Explicit inventory: literature-sourced vs. assumed parameters

**Literature-sourced, individually referenced, unit-conversion-verified:**
`kinetic_parameters.m` Section 1 (R1-R23 rate constants) -- Fruhwirt et
al., *Phys. Chem. Chem. Phys.* 22 (2020) 5647.

**Functional form from literature, specific constants are engineering
assumptions (flagged `TODO/ASSUMPTION` at each definition site):** Pt
dissolution/redeposition rate constants (functional form after Darling &
Meyers 2003 and Rinaldo et al. 2011); percolation exponent `tau` (Stauffer
& Aharony 1994); crossover permeability activation energy; ORR exchange
current density / Tafel coefficient (representative literature ranges per
Neyerlin et al. 2006); IEC/thickness decay rate constants; water-uptake
and mechanical/thermal correlations in `membrane_properties.m`.

**Corrected in this engagement, order-of-magnitude estimate with
documented derivation:** initial Fe2+ contamination level; initial Pt
loading (0.4 mg/cm^2 assumed, per Gasteiger et al. 2005); peroxide-
generation efficiency factor `k_gen`.

## 5. Secondary analysis scripts -- not re-verified in this pass

`species_evolution.m`, `time_scale_analysis.m`, `mechanism_ranking.m`,
`parameter_sensitivity.m`, `interaction_analysis.m`, `heatmap_analysis.m`,
and `generate_presentation_data.m` each contain their own separate,
self-contained simplified ODE implementation (not calling `ode_system.m`
or `reaction_rates.m` in full-coupled mode) -- this separation predates
this engagement. Given the scope of debugging required for the master
pipeline (Sections 1-3), these were not rewritten to the same standard.

- `species_evolution.m`, `time_scale_analysis.m`, `generate_presentation_data.m`:
  timed out during testing, consistent with their local ODE
  implementations sharing the class of stiffness/initial-condition issues
  fixed for the master pipeline, but not yet independently corrected.
- `mechanism_ranking.m`: uses MATLAB's `table` type, unavailable in
  Octave -- could not be executed in the test environment.
- `interaction_analysis.m`: uses MATLAB's `digraph` graph-plot objects,
  unavailable in Octave -- could not be executed in the test environment.
- `parameter_sensitivity.m`: partially executes (multiple sweep
  sub-cases completed before a test timeout), suggesting the underlying
  approach is functional but slow for a full sweep at full duration.
- `heatmap_analysis(sens)` requires a `sens` struct from
  `parameter_sensitivity.m` as input -- it is a downstream consumer, not a
  standalone script; the "undefined `sens`" error seen when called with no
  arguments during testing is expected usage-order behavior, not a bug.

**Recommendation for a follow-up pass:** apply the same fixes that
resolved the master pipeline (corrected Fe2+ initial condition, corrected
crossover flux units, QSSA radical reduction or loosened tolerance) to
each file's local ODE implementation, and verify `mechanism_ranking.m` /
`interaction_analysis.m` directly in MATLAB.

## 6. Files not modified

`operating_conditions.m`, `mechanism_switch.m`, `membrane_attack.m`,
`polymer_unzipping.m`, `fenton_chemistry.m`, `get_audited_parameters.m`
were reviewed but not modified (no blocking bugs found in the
master-pipeline call paths that exercise them; the latter two are
currently orphaned/standalone diagnostic files not called by the master
pipeline, consistent with their state prior to this engagement).
