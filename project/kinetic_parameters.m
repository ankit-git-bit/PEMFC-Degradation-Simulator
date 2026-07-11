function kinetics = kinetic_parameters()
    % KINETIC_PARAMETERS.M
    % Stores kinetic rate constants, activation energies, and 
    % degradation coefficients for the PEMFC model.

    %% 1. Original Chemical Reactions (R1 - R23)
    % Iron Reactions (R1 - R5)
    kinetics.R1.A         = 105000;
    kinetics.R1.Ea        = 35400;
    kinetics.R1.units     = 'm^3/(mol s)';
    kinetics.R1.reference = 'Frühwirt et al. (2020), Table 1';

    kinetics.R2.A         = 8.43e15;
    kinetics.R2.Ea        = 126000;
    kinetics.R2.units     = 'm^3/(mol s)';
    kinetics.R2.reference = 'Frühwirt et al. (2020), Table 1';

    kinetics.R3.A         = 13700000;
    kinetics.R3.Ea        = 19000;
    kinetics.R3.units     = 'm^3/(mol s)';
    kinetics.R3.reference = 'Frühwirt et al. (2020), Table 1';

    kinetics.R4.A         = 2.74e10;
    kinetics.R4.Ea        = 42000;
    kinetics.R4.units     = 'm^3/(mol s)';
    kinetics.R4.reference = 'Frühwirt et al. (2020), Table 1';

    kinetics.R5.A         = 49000000;
    kinetics.R5.Ea        = 33000;
    kinetics.R5.units     = 'm^3/(mol s)';
    kinetics.R5.reference = 'Frühwirt et al. (2020), Table 1';

    % H2O2 Decomposition (R6)
    kinetics.R6.A         = 1.0e13;
    kinetics.R6.Ea        = 201000;
    kinetics.R6.units     = 's^-1';
    kinetics.R6.reference = 'Frühwirt et al. (2020), Table 1';

    % Radical Reactions (R7 - R13)
    kinetics.R7.A         = 8430000;
    kinetics.R7.Ea        = 14000;
    kinetics.R7.units     = 'm^3/(mol s)';
    kinetics.R7.reference = 'Frühwirt et al. (2020), Table 1';

    kinetics.R8.A         = 2710;
    kinetics.R8.Ea        = 33500;
    kinetics.R8.units     = 'm^3/(mol s)';
    kinetics.R8.reference = 'Frühwirt et al. (2020), Table 1';

    kinetics.R9.A         = 3700000;
    kinetics.R9.Ea        = 20600;
    kinetics.R9.units     = 'm^3/(mol s)';
    kinetics.R9.reference = 'Frühwirt et al. (2020), Table 1';

    kinetics.R10.A        = 13300000;
    kinetics.R10.Ea       = 8000;
    kinetics.R10.units    = 'm^3/(mol s)';
    kinetics.R10.reference= 'Frühwirt et al. (2020), Table 1';

    kinetics.R11.A        = 3.39e9;
    kinetics.R11.Ea       = 14200;
    kinetics.R11.units    = 'm^3/(mol s)';
    kinetics.R11.reference= 'Frühwirt et al. (2020), Table 1';

    kinetics.R12.A        = 91400000;
    kinetics.R12.Ea       = 19200;
    kinetics.R12.units    = 'm^3/(mol s)';
    kinetics.R12.reference= 'Frühwirt et al. (2020), Table 1';

    kinetics.R13.A        = 1.37e9;
    kinetics.R13.Ea       = 10300;
    kinetics.R13.units    = 'm^3/(mol s)';
    kinetics.R13.reference= 'Frühwirt et al. (2020), Table 1';

    % Membrane Degradation (R14 - R16)
    kinetics.R14.A        = 6.8e15;
    kinetics.R14.Ea       = 70000;
    kinetics.R14.units    = 'm^3/(mol s)';
    kinetics.R14.reference= 'Frühwirt et al. (2020), Table 1';

    kinetics.R15.A        = 5.51e16;
    kinetics.R15.Ea       = 70000;
    kinetics.R15.units    = 'm^3/(mol s)';
    kinetics.R15.reference= 'Frühwirt et al. (2020), Table 1';

    kinetics.R16.A        = 1.56e17;
    kinetics.R16.Ea       = 70000;
    kinetics.R16.units    = 'm^3/(mol s)';
    kinetics.R16.reference= 'Frühwirt et al. (2020), Table 1';

    % Backbone Unzipping (R17 - R23)
    kinetics.R17_23.A         = 1.84e15;
    kinetics.R17_23.Ea        = 70000;
    kinetics.R17_23.units     = 'm^3/(mol s)';
    kinetics.R17_23.reference = 'Frühwirt et al. (2020), Table 1';

    %% 2. Pt KINETICS
    % TODO/ASSUMPTION: Unlike Section 1 (Fruhwirt et al. 2020, individually
    % referenced), the Pt dissolution/redeposition pre-exponential factors
    % below are ORDER-OF-MAGNITUDE ENGINEERING ESTIMATES, not values taken
    % directly from a single cited source. The functional form (Butler-Volmer
    % type potential dependence) follows the widely used Pt dissolution
    % framework of Darling & Meyers (J. Electrochem. Soc. 150 (2003) A1523)
    % and Rinaldo, Lee, Stumper & Eikerling (Electrochim. Acta 57 (2011) 273),
    % but the specific rate constants here have NOT been fit to a specific
    % literature dataset for this project and should be treated as
    % placeholders pending calibration against a published Pt-loss vs. cycle
    % curve (e.g. Ohma et al. cyclic-voltammetry ECSA-loss protocols).
    kinetics.Pt.k_diss_A    = 2e-16;       % Dissolution pre-exponential - ASSUMPTION
                                            % (calibrated so that Pt loss is
                                            % a modest ~5-15% over a 5000 h
                                            % simulated lifetime under normal
                                            % operating conditions, consistent
                                            % with the ORDER OF MAGNITUDE of
                                            % ECSA loss reported in PEMFC
                                            % durability studies, e.g. Ohma
                                            % et al., ECS Trans. 41 (2011)
                                            % 775; an earlier un-calibrated
                                            % value produced near-total Pt
                                            % depletion within ~11 hours,
                                            % which is not physically
                                            % plausible and is corrected here)
    kinetics.Pt.k_redep_A   = 2e-17;       % Redeposition pre-exponential - ASSUMPTION
    kinetics.Pt.alpha       = 0.5;         % Charge transfer coefficient (standard BV default)
    kinetics.Pt.n           = 2;           % Electrons per Pt atom (Pt -> Pt2+ + 2e-, stoichiometric)
    kinetics.Pt.Ea          = 45000;       % Activation Energy (J/mol) - ASSUMPTION, order-of-magnitude

    %% 3. BIDIRECTIONAL & CONDUCTIVITY COEFFICIENTS
    % TODO/ASSUMPTION: percolation-theory-motivated form; see documented
    % assumption on p.tau in parameters.m (Stauffer & Aharony, 1994).
    kinetics.coupling.tau       = 1.4;     % Percolation threshold exponent - ASSUMPTION
    kinetics.coupling.k_ECSA    = 1e-9;    % ECSA degradation rate constant - ASSUMPTION
    kinetics.coupling.sigma_ref = 10.0;    % Reference conductivity (S/m), order-of-magnitude
                                            % for hydrated Nafion (Kreuer 2004 reports ~5-10 S/m)

    %% 4. CROSSOVER COEFFICIENTS
    % TODO/ASSUMPTION: functional form only (Arrhenius T-dependence of gas
    % permeability is standard, e.g. Kreuer 2004; Karpenko-Jereb et al.,
    % Int. J. Hydrogen Energy). Numerical value of Ea_perm is an engineering
    % estimate, not fit to a specific permeability-vs-T dataset here.
    kinetics.crossover.Ea_perm    = 20000; % Activation energy for gas permeability (J/mol) - ASSUMPTION
    kinetics.crossover.HF_damage  = 1.0;   % Scaling factor for HF-induced pinholes - ASSUMPTION
    kinetics.crossover.T_ref      = 353.15;% Reference temperature (K)

    %% 5. PERFORMANCE KINETICS
    % TODO/ASSUMPTION: representative PEMFC cathode ORR values, not fit to a
    % specific polarization curve for this project (typical literature
    % ranges: i0_ref ~1e-6 to 1e-3 A/m^2 depending on Pt loading/RH/T,
    % e.g. Neyerlin, Gu, Jorne & Gasteiger, J. Electrochem. Soc. 153 (2006)
    % A1955; alpha_orr ~0.5 is a standard Butler-Volmer default).
    kinetics.perf.i0_ref     = 1e-4;       % Exchange current density (A/m^2) - ASSUMPTION
    kinetics.perf.alpha_orr  = 0.5;        % ORR charge transfer coefficient (standard default)
    kinetics.perf.V_ocv      = 0.95;       % Open circuit voltage (V), typical PEMFC value

    kinetics.reference = 'Frühwirt et al. (2020) & Coupled Degradation Model';
end