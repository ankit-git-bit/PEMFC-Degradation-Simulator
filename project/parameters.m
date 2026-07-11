function p = parameters()
  
    %% 1. PHYSICAL CONSTANTS
    p.R = 8.314462618;             % Gas constant (J/mol/K)
    p.F = 96485.33212;             % Faraday constant (C/mol)
    p.Na = 6.02214076e23;          % Avogadro's number (1/mol)
    p.kB = 1.380649e-23;           % Boltzmann constant (J/K)
    %% 2. OPERATING CONDITIONS
    p.T = 363.15;                  % Temperature (K)
    p.RH = 0.30;                   % Relative Humidity (fraction)
    p.pressure = 2.3;              % Pressure (bar)
    p.voltage = 1.0;               % Cell voltage (V)
    %% 3. MEMBRANE PROPERTIES
    p.EW = 1100;                   % Equivalent weight (g/mol)
    p.rho_dry = 1980;              % Dry density (kg/m³)
    p.delta_membrane_0 = 50e-6;    % Initial thickness (m)
    p.A_membrane = 25e-4;          % Active area (m²)
    %% 4. FIXED CHARGE & CF2 GROUPS
    p.nu = 1800;                   % mol/m³ (rho_dry/EW)
    p.SC_SO3H_0 = 1800;            % Initial side chain sulfonic acid (mol/m³)
    p.CF2_groups_per_SC = 13;      % Groups per side chain
    p.CF2_total_0 = 23400;         % Total CF2 groups (mol/m³)
    p.CF2_7_0 = 1800;              % Unzipping position 7 (mol/m³)
    p.CF2_6_0 = 1800;              % Unzipping position 6 (mol/m³)
    p.CF2_5_0 = 1800;              % Unzipping position 5 (mol/m³)
    p.CF2_4_0 = 1800;              % Unzipping position 4 (mol/m³)
    p.CF2_3_0 = 1800;              % Unzipping position 3 (mol/m³)
    p.CF2_2_0 = 1800;              % Unzipping position 2 (mol/m³)
    p.CF2_1_0 = 1800;              % Unzipping position 1 (mol/m³)
    p.CF2_backbone_0 = 10800;      % Positions 8-13 (mol/m³)
    %% 5. WEAK END GROUPS
    p.carboxyl_0 = 1170;           % Initial carboxyl groups (mol/m³)
    p.wpe_0 = 1170;                % Initial weak polymer end groups (mol/m³)
    %% 6. INITIAL SPECIES CONCENTRATIONS
    % FIX (previous bug): ppm_to_mol_m3 was set to 17.9, giving
    % Fe2_0 = 89.5 mol/m^3 (i.e. ~5 g/L of dissolved iron) for a stated "5
    % ppm" contamination level - approximately 500x too concentrated
    % compared to a physically sensible ppm(mass)-to-molarity conversion,
    % and it made the coupled ODE numerically and physically unrealistic
    % (predicted side-chain sulfonic-acid depletion within SECONDS rather
    % than the literature-expected thousands-of-hours timescale for PEMFC
    % membrane degradation). Corrected conversion, assuming an aqueous-
    % equivalent reference density of 1000 kg/m^3 (1 ppm by mass = 1 g Fe
    % per m^3 of reference solution) and Fe molar mass 55.845 g/mol
    % (IUPAC): 1 ppm -> (1 g/m^3) / (55.845 g/mol) = 0.0179 mol/m^3 per ppm.
    % This order of magnitude (Fe2_0 ~ 0.09 mol/m^3 for 5 ppm) is also
    % consistent with an independent membrane-mass-density-based estimate
    % (5 ppm x ~1980 kg/m^3 Nafion density / 55.845 g/mol =~ 0.18 mol/m^3),
    % giving confidence this is now the right order of magnitude rather
    % than a specific validated literature value - flagged as
    % TODO/ASSUMPTION pending a project-specific ex-situ Fe measurement.
    p.Fe_molar_mass = 55.845;                    % g/mol (IUPAC)
    p.ppm_to_mol_m3 = 1000 / p.Fe_molar_mass / 1000; % mol/m^3 per ppm (= 0.0179)
    p.Fe_contamination_ppm = 5;                  % ppm, representative value - ASSUMPTION
    p.Fe2_0 = p.Fe_contamination_ppm * p.ppm_to_mol_m3; % Initial Fe2+ concentration (mol/m^3)
    p.Fe3_0 = 0;                   % Initial Fe3+ (mol/m³)
    p.OH_0 = 1e-18;                % Initial OH radical (mol/m³)
    p.OOH_0 = 1e-18;               % Initial OOH radical (mol/m³)
    p.H_0 = 1e-18;                 % Initial H radical (mol/m³)
    p.HF_0 = 0;                    % Initial HF (mol/m³)
    p.CO2_0 = 0;                   % Initial CO2 (mol/m³)
    p.HOCF2CF2SO3H_0 = 0;          % Initial fragmented side chain (mol/m³)
    p.SC_O_0 = 0;                  % Initial side chain radical (mol/m³)
    p.BB_O_0 = 0;                  % Initial backbone radical (mol/m³)
    %% 7. HENRY CONSTANTS
    p.Henry_H2 = 1.21e8;           % Henry constant for H2 (Pa)
    p.Henry_O2 = 5.14e8;           % Henry constant for O2 (Pa)
    p.Henry_N2 = 5.14e8;           % Henry constant for N2 (Pa)
    %% 8. GAS DIFFUSION COEFFICIENTS
    p.D_H2 = 2.09e-10;             % Diffusion coeff H2 (m²/s)
    p.D_O2 = 9.73e-11;             % Diffusion coeff O2 (m²/s)
    p.D_N2 = 9.73e-11;             % Diffusion coeff N2 (m²/s)
    %% 9. REACTION RATE CONSTANTS
    p.k_R1_A = 105000;             % m³/mol/s
    p.k_R1_Ea = 35400;             % J/mol
    p.k_R2_A = 8.43e15;
    p.k_R2_Ea = 126000;
    p.k_R3_A = 13700000;
    p.k_R3_Ea = 19000;
    p.k_R4_A = 2.74e10;
    p.k_R4_Ea = 42000;
    p.k_R5_A = 49000000;
    p.k_R5_Ea = 33000;
    p.k_R6_A = 1.0e13;             % s^-1
    p.k_R6_Ea = 201000;
    p.k_R7_A = 8430000;
    p.k_R7_Ea = 14000;
    p.k_R8_A = 2710;
    p.k_R8_Ea = 33500;
    p.k_R9_A = 3700000;
    p.k_R9_Ea = 20600;
    p.k_R10_A = 13300000;
    p.k_R10_Ea = 8000;
    p.k_R11_A = 3.39e9;
    p.k_R11_Ea = 14200;
    p.k_R12_A = 91400000;
    p.k_R12_Ea = 19200;
    p.k_R13_A = 1.37e9;
    p.k_R13_Ea = 10300;
    p.k_R14_A = 6.8e15;
    p.k_R14_Ea = 70000;
    p.k_R15_A = 5.51e16;
    p.k_R15_Ea = 70000;
    p.k_R16_A = 1.56e17;
    p.k_R16_Ea = 70000;
    p.k_R17_23_A = 1.84e15;
    p.k_R17_23_Ea = 70000;
    %% 10. CROSSOVER PERMEABILITY
    p.permeability_H2 = 1.0e-12;   % mol/(m·s·Pa)
    p.permeability_O2 = 5.0e-13;   % mol/(m·s·Pa)
    %% 11. SOLVER SETTINGS
    p.rel_tol = 1e-10;
    p.abs_tol = 1e-14;
    p.max_step = 3600;             
    p.n_points = 10000;            
    % TODO/LIMITATION: t_final controls the simulated duration for the
    % standalone diagnostic/mechanism-isolation scripts (species_evolution.m,
    % time_scale_analysis.m, mechanism_ranking.m, parameter_sensitivity.m).
    % These scripts use their own self-contained local ODE step functions
    % (evo_ode_step, ts_ode_step, ranked_ode_step, sens_ode_step) that
    % integrate OH/OOH/H as explicit differential states, WITHOUT the
    % quasi-steady-state approximation (QSSA) applied in reaction_rates.m /
    % ode_system.m for the master lifetime_simulation.m pipeline (see the
    % QSSA rationale note in reaction_rates.m). Without that reduction,
    % these scripts are subject to the same extreme radical-chemistry
    % stiffness that made the original (pre-fix) master simulation
    % numerically intractable, and a full 720-hour (30-day) continuous
    % integration was found to be impractically slow in this environment.
    % t_final is set here to a shorter, computationally tractable window
    % suitable for demonstrating mechanism-isolation behavior; porting the
    % QSSA reduction into these four local step functions (mechanical, same
    % pattern as reaction_rates.m) is a documented follow-up recommendation
    % for extending them to multi-week diagnostic horizons.
    p.t_final = 3600 * 2;          % 2 hours (was 30 days / 720 h) - see note above
    %% 12. CATALYST & PHYSICAL PROPERTIES
    p.Pt_density = 21450;          % Density of Platinum (kg/m³)
    p.Pt_particle_radius_0 = 2e-9; % Initial Pt particle radius (m)
    p.ECSA_0 = 60;                 % Initial ECSA (m²/g_Pt)
    %% 13. MEMBRANE PROPERTY BASELINES
    p.IEC_0 = 0.91;                % Initial Ion Exchange Capacity (meq/g)
    p.sigma_0 = 10.0;              % Initial conductivity (S/m)
    %% 14. DEGRADATION & COUPLING LIMITS
    p.HF_critical = 10.0;          % HF concentration threshold (mol/m³)
    p.Ea_perm = 20000;             % Activation energy for gas permeability (J/mol)
    p.T_ref = 353.15;              % Reference temperature (K)
    %% 15. PERFORMANCE MODEL CONSTANTS
    p.i0_ref = 1e-4;               % Reference exchange current (A/m²)
    p.alpha = 0.5;                 % Charge transfer coefficient
    p.V_ocv = 0.95;                % Open circuit voltage (V)
    p.k_ECSA = 1e-9;               % ECSA degradation rate constant
    %% 16. PLATINUM SPECIES (ADDED - required for expanded state vector)
    % TODO/ASSUMPTION: Pt_0 is a lumped catalyst-layer volumetric metallic-Pt
    % concentration, not directly reported in the literature sources used
    % elsewhere in this project. It is estimated here from a typical cathode
    % Pt loading of 0.4 mg_Pt/cm^2 (representative value, see e.g. Gasteiger,
    % Kocha, Sompalli & Wagner, Appl. Catal. B 56 (2005) 9-35) divided by an
    % assumed catalyst-layer thickness of 10 micron and Pt molar mass
    % (195.084 g/mol). This is an ENGINEERING ASSUMPTION, not a literature
    % constant -- verify against the actual MEA/catalyst-layer design before
    % treating absolute Pt-loss numbers as quantitative.
    p.Pt_molar_mass   = 195.084;       % g/mol (IUPAC standard atomic weight)
    p.Pt_loading_mg_cm2 = 0.4;         % mg/cm^2, representative cathode loading
                                        % (Gasteiger et al. 2005) - ASSUMPTION
    p.CL_thickness_0  = 10e-6;         % m, assumed catalyst layer thickness - ASSUMPTION
    p.Pt_0  = (p.Pt_loading_mg_cm2*1e-3*1e4) / p.Pt_molar_mass / p.CL_thickness_0; % mol/m^3
    p.Pt2_0 = 0;                       % Initial dissolved Pt2+ (mol/m^3)
    p.rho_Pt = p.Pt_density;           % alias used by some modules (kg/m^3)
    p.Pt_radius = p.Pt_particle_radius_0; % alias used by parameter_sensitivity.m (m)

    %% 17. PEROXIDE INITIAL CONDITION
    % TODO/ASSUMPTION: small non-zero seed to avoid singular multiplicative
    % rate laws at t=0 (same convention used for OH_0/OOH_0/H_0 above).
    % Steady-state H2O2 is governed dynamically by crossover generation
    % (hydrogen_crossover.m + peroxide_generation.m) vs. Fenton consumption,
    % not by this seed value.
    p.H2O2_0 = 1e-6;                   % mol/m^3

    %% 18. CROSSOVER SENSITIVITY BASELINE
    p.H2_crossover_0 = p.permeability_H2; % mol/(m*s*Pa), BOL baseline for sweeps

    %% 19. PHYSICAL-PROPERTY EVOLUTION CONSTANTS (degradation_update.m)
    % TODO/ASSUMPTION: tau is the conductivity-vs-IEC percolation exponent.
    % Percolation theory gives exponents typically in the range 1.3-2.0
    % depending on dimensionality/connectivity (Stauffer & Aharony,
    % "Introduction to Percolation Theory", Taylor & Francis, 1994). The
    % mid-range value below is adopted as a representative default and is
    % NOT fit to a specific PEMFC conductivity dataset for this project.
    p.tau = 1.4;
    % TODO/ASSUMPTION: k_IEC and k_thin are phenomenological first-order
    % rate constants that translate chemical side-chain-loss / backbone-loss
    % rates (mol/m^3/s, from reaction_rates.m) into fractional IEC and
    % thickness loss per second. They are NOT independently reported in the
    % literature; they are order-of-magnitude placeholders chosen so that
    % IEC and thickness show visible decay over a 1000-5000 h simulated
    % lifetime, consistent with the *magnitude* of losses reported in
    % accelerated membrane degradation studies (e.g. Ohma et al.,
    % ECS Trans. 41 (2011) 775; Kreuer, Chem. Rev. 104 (2004) 4637 for
    % IEC/conductivity relationships). These require calibration against a
    % specific accelerated stress test dataset before being used for
    % quantitative lifetime prediction -- flagged here rather than
    % presented as a validated literature value.
    p.k_IEC  = 1e-4;                   % 1/s per unit normalized SC_SO3H loss rate
    p.k_thin = 1e-6;                   % m per (mol/m^3/s) of backbone loss rate
end