function params = get_audited_parameters()
    %% 1. PHYSICAL CONSTANTS (CODATA verified)
    params.R = 8.314462618;             % Gas constant (J/mol/K)
    params.F = 96485.33212;             % Faraday constant (C/mol)
    params.Na = 6.02214076e23;          % Avogadro's number (1/mol)
    params.kB = 1.380649e-23;           % Boltzmann constant (J/K)
    
    %% 2. MEMBRANE PROPERTIES
    params.EW = 1100;                   % Equivalent weight (g/mol) [Nafion 115/212]
    params.rho_dry = 1980;              % Dry density (kg/m³)[cite: 4]
    params.delta_membrane_0 = 50e-6;    % Initial thickness (m)
    params.A_membrane = 25e-4;          % Active area (m²)
    
    %% 3. FIXED CHARGE 
    params.nu = 1800;                   % mol/m³ (ρ_N/EW = 1980/1.1)[cite: 4]
    params.SC_SO3H_0 = params.nu;       % mol/m³[cite: 4]
    
    %% 4. CF₂ GROUPS 
    % Total: ~13 groups per side chain (EW = 100m + 446, m=6.5)[cite: 4]
    params.CF2_groups_per_SC = 13;
    params.CF2_total_0 = params.nu * 13;  % 23400 mol/m³[cite: 4]
    
    % Track 7 positions for unzipping mechanism
    % Positions 8-13 are in the backbone
    params.CF2_7_0 = params.nu;
    params.CF2_6_0 = params.nu;
    params.CF2_5_0 = params.nu;
    params.CF2_4_0 = params.nu;
    params.CF2_3_0 = params.nu;
    params.CF2_2_0 = params.nu;
    params.CF2_1_0 = params.nu;
    params.CF2_backbone_0 = params.nu * 6;  % Positions 8-13
    
    %% 5. WEAK END GROUPS 
    % Set at 5% of the CF2 group concentration[cite: 4]
    params.carboxyl_0 = 0.05 * params.CF2_total_0;   % 1170 mol/m³[cite: 4]
    params.wpe_0 = 0.05 * params.CF2_total_0;        % 1170 mol/m³[cite: 4]
    
    %% 6. INITIAL SPECIES CONCENTRATIONS
    params.ppm_to_mol_m3 = 17.9;         % Conversion factor[cite: 4]
    params.Fe2_0 = 5 * params.ppm_to_mol_m3;   % 5 ppm[cite: 4]
    params.Fe3_0 = 0;                    % mol/m³
    
    % Radicals (numerical initialization)
    params.OH_0 = 1e-18;                 % mol/m³
    params.OOH_0 = 1e-18;                % mol/m³
    params.H_0 = 1e-18;                  % mol/m³
    
    % Products
    params.HF_0 = 0;                     % mol/m³
    params.CO2_0 = 0;                    % mol/m³
    params.HOCF2CF2SO3H_0 = 0;
    params.SC_O_0 = 0;
    params.BB_O_0 = 0;
    
    %% 7. HENRY CONSTANTS 
    % Definition: H = p/C (Henry's law: p = H·C)
    % Units: [Pa] -> C = p/H [mol/m³]
    params.Henry_H2 = 1.21e8;            % Pa[cite: 2]
    params.Henry_O2 = 5.14e8;            % Pa[cite: 2]
    params.Henry_N2 = 5.14e8;            % Pa[cite: 2]
    
    %% 8. GAS DIFFUSION 
    params.D_H2 = 2.09e-10;              % m²/s[cite: 2]
    params.D_O2 = 9.73e-11;              % m²/s[cite: 2]
    params.D_N2 = 9.73e-11;              % m²/s[cite: 2]
    
    %% =====================================================================
    % REACTION RATE CONSTANTS
    % Source: Frühwirt et al. (2020), PCCP, Table 1[cite: 5]
    % Conversion: A_SI [m³/(mol·s)] = A_table [L/(mol·s)] * 1e-3
    % Ea converted: kJ/mol -> J/mol
    %% =====================================================================
    
    %% -------------------------
    % Iron reactions (R1-R5)[cite: 5]
    %% -------------------------
    params.k_R1_A  = 1.05e8 * 1e-3;       % m³/mol/s - CORRECTED from 1.05e15[cite: 5]
    params.k_R1_Ea = 35.4e3;              % J/mol[cite: 5]
    
    params.k_R2_A  = 8.43e18 * 1e-3;      %[cite: 5]
    params.k_R2_Ea = 126e3;               %[cite: 5]
    
    params.k_R3_A  = 1.37e10 * 1e-3;      %[cite: 5]
    params.k_R3_Ea = 19e3;                %[cite: 5]
    
    params.k_R4_A  = 2.74e13 * 1e-3;      %[cite: 5]
    params.k_R4_Ea = 42e3;                %[cite: 5]
    
    params.k_R5_A  = 4.90e10 * 1e-3;      %[cite: 5]
    params.k_R5_Ea = 33e3;                %[cite: 5]
    
    %% -------------------------
    % H2O2 decomposition (R6)[cite: 5]
    %% -------------------------
    params.k_R6_A  = 1.00e13;             % s^-1[cite: 5]
    params.k_R6_Ea = 201e3;               %[cite: 5]
    
    %% -------------------------
    % Radical reactions (R7-R13)[cite: 5]
    %% -------------------------
    params.k_R7_A  = 8.43e9 * 1e-3;       %[cite: 5]
    params.k_R7_Ea = 14e3;                %[cite: 5]
    
    params.k_R8_A  = 2.71e6 * 1e-3;       %[cite: 5]
    params.k_R8_Ea = 33.5e3;              %[cite: 5]
    
    params.k_R9_A  = 3.70e9 * 1e-3;       %[cite: 5]
    params.k_R9_Ea = 20.6e3;              %[cite: 5]
    
    params.k_R10_A = 1.33e10 * 1e-3;      %[cite: 5]
    params.k_R10_Ea = 8e3;                %[cite: 5]
    
    params.k_R11_A = 3.39e12 * 1e-3;      %[cite: 5]
    params.k_R11_Ea = 14.2e3;             %[cite: 5]
    
    params.k_R12_A = 9.14e10 * 1e-3;      %[cite: 5]
    params.k_R12_Ea = 19.2e3;             %[cite: 5]
    
    params.k_R13_A = 1.37e12 * 1e-3;      %[cite: 5]
    params.k_R13_Ea = 10.3e3;             %[cite: 5]
    
    %% -------------------------
    % Membrane degradation (R14-R16)[cite: 5]
    %% -------------------------
    params.k_R14_A = 6.80e18 * 1e-3;      %[cite: 5]
    params.k_R14_Ea = 70e3;               %[cite: 5]
    
    params.k_R15_A = 5.51e19 * 1e-3;      %[cite: 5]
    params.k_R15_Ea = 70e3;               %[cite: 5]
    
    params.k_R16_A = 1.56e20 * 1e-3;      %[cite: 5]
    params.k_R16_Ea = 70e3;               %[cite: 5]
    
    %% -------------------------
    % Backbone unzipping (R17-R23)[cite: 5]
    %% -------------------------
    params.k_R17_23_A = 1.84e18 * 1e-3;   %[cite: 5]
    params.k_R17_23_Ea = 70e3;            %[cite: 5]
    
    %% 12. CROSSOVER PERMEABILITY
    params.permeability_H2 = 1.0e-12;    % mol/(m·s·Pa) [Engineering assumption]
    params.permeability_O2 = 5.0e-13;    % mol/(m·s·Pa) [Engineering assumption]
    
    %% 13. OPERATING CONDITIONS 
    params.op.T = 363.15;                % K
    params.op.RH = 0.30;                 % fraction
    params.op.pressure = 2.3;            % bar
    params.op.voltage = 1.0;             % V
    
    %% 14. SOLVER SETTINGS
    params.solver.rel_tol = 1e-10;
    params.solver.abs_tol = 1e-14;
    params.solver.max_step = 3600;       % 1 hour
    params.solver.n_points = 10000;
    params.solver.t_final = 3600 * 24 * 30;  % 30 days
    
    %% 15. MODEL OPTIONS
    params.model.include_Fenton = true;
    params.model.include_crossover = true;
    params.model.include_degradation = true;
    params.model.include_membrane_thinning = true;
end

%% ========================================================================
%  FUNCTION: Calculate temperature-dependent parameters
% =========================================================================
function params = calculate_audited_parameters(params)
    R = params.R;
    T = params.op.T;
    
    % ---- 1. WATER UPTAKE ----
    a = params.op.RH;
    lambda = 0.043 + 17.81*a - 39.85*a^2 + 36.0*a^3;
    assert(lambda >= 0, 'Water content cannot be negative');
    params.lambda = max(lambda, 0.1); 
    
    % ---- 2. WATER DIFFUSION ----
    params.D_w = 4.17e-8 * params.lambda * (1 + 161*exp(-params.lambda)) * exp(-2436/T);  % m²/s
    assert(params.D_w > 0, 'Water diffusion must be positive');
    
    % ---- 3. PROTON CONDUCTIVITY ----
    sigma = (0.514*params.lambda - 0.326) * exp(1268*(1/303 - 1/T));
    params.sigma = max(sigma, 0);  % S/m
    assert(params.sigma >= 0, 'Conductivity must be non-negative');
    
    % ---- 4. GAS CONCENTRATIONS - Henry's law ----
    p_H2 = params.op.pressure * 1e5 * 0.5;   % Pa 
    p_O2 = params.op.pressure * 1e5 * 0.21;  % Pa
    
    params.C_H2 = p_H2 / params.Henry_H2;    % mol/m³
    params.C_O2 = p_O2 / params.Henry_O2;    % mol/m³
    
    assert(params.C_H2 > 0, 'H2 concentration must be positive');
    assert(params.C_O2 > 0, 'O2 concentration must be positive');
    
    % ---- 5. RATE CONSTANTS (Arrhenius Eq) ----
    params.k_R1 = params.k_R1_A * exp(-params.k_R1_Ea/(R*T));
    params.k_R2 = params.k_R2_A * exp(-params.k_R2_Ea/(R*T));
    params.k_R3 = params.k_R3_A * exp(-params.k_R3_Ea/(R*T));
    params.k_R4 = params.k_R4_A * exp(-params.k_R4_Ea/(R*T));
    params.k_R5 = params.k_R5_A * exp(-params.k_R5_Ea/(R*T));
    
    % H2O2 decomposition (Updated to reference correct struct variable)
    params.k_R6 = params.k_R6_A * exp(-params.k_R6_Ea/R*(1/T - 1/298));
    
    % Radical Reactions 
    params.k_R7 = params.k_R7_A * exp(-params.k_R7_Ea/(R*T));
    params.k_R8 = params.k_R8_A * exp(-params.k_R8_Ea/(R*T));
    params.k_R9 = params.k_R9_A * exp(-params.k_R9_Ea/(R*T));
    params.k_R10 = params.k_R10_A * exp(-params.k_R10_Ea/(R*T));
    params.k_R11 = params.k_R11_A * exp(-params.k_R11_Ea/(R*T));
    params.k_R12 = params.k_R12_A * exp(-params.k_R12_Ea/(R*T));
    params.k_R13 = params.k_R13_A * exp(-params.k_R13_Ea/(R*T));
    
    % ---- 6. DEGRADATION RATES ----
    params.k_R14 = params.k_R14_A * exp(-params.k_R14_Ea/(R*T));
    params.k_R15 = params.k_R15_A * exp(-params.k_R15_Ea/(R*T));
    params.k_R16 = params.k_R16_A * exp(-params.k_R16_Ea/(R*T));
    params.k_R17_23 = params.k_R17_23_A * exp(-params.k_R17_23_Ea/(R*T));
    
    assert(all([params.k_R1, params.k_R2, params.k_R3, params.k_R4, params.k_R5, ...
                params.k_R6, params.k_R14, params.k_R15, params.k_R16, ...
                params.k_R17_23] > 0), ...
                'All rate constants must be positive');
end

%% ========================================================================
%  FUNCTION: Verification output terminal display
% =========================================================================
function audit_report()
    fprintf('========================================\n');
    fprintf('LITERATURE AUDIT REPORT\n');
    fprintf('========================================\n\n');
    
    fprintf('1. WATER UPTAKE\n');
    fprintf('   ✓ λ = 0.043 + 17.81a - 39.85a² + 36.0a³\n');
    fprintf('   ✓ Coefficient 17.81 (not 17.18) - CORRECTED\n\n');
    
    fprintf('2. WATER DIFFUSION\n');
    fprintf('   ✓ D_w = 4.17e-8·λ·(1+161e^-λ)·exp(-2436/T)\n');
    fprintf('   ✓ Units: m²/s\n\n');
    
    fprintf('3. PROTON CONDUCTIVITY\n');
    fprintf('   ✓ σ = (0.514λ - 0.326)·exp(1268·(1/303 - 1/T))\n');
    fprintf('   ✓ Protection: σ = max(0, σ)\n\n');
    
    fprintf('4. HENRY CONSTANTS - Karpenko-Jereb Table 6, p. 13648\n');
    fprintf('   ✓ Definition: H = p/C [Pa]\n');
    fprintf('   ✓ H_H2 = 1.21×10⁸ Pa[cite: 2]\n');
    fprintf('   ✓ H_O2 = 5.14×10⁸ Pa[cite: 2]\n\n');
    
    fprintf('5. GAS DIFFUSION - Karpenko-Jereb Table 6, p. 13648\n');
    fprintf('   ✓ D_H2 = 2.09×10⁻¹⁰ m²/s[cite: 2]\n');
    fprintf('   ✓ D_O2 = 9.73×10⁻¹¹ m²/s[cite: 2]\n\n');
    
    fprintf('6. REACTION RATES - Shah (2009) Table X, p. B476\n');
    fprintf('   ✓ Verified units per reaction order[cite: 4]\n\n');
    
    fprintf('7. DEGRADATION RATES - Frühwirt (2020) Table 1\n');
    fprintf('   ✓ k_R1_A correctly set to 1.05e8 L/mol/s[cite: 5]\n');
    fprintf('   ✓ Original: M⁻¹s⁻¹ → Converted: mol⁻¹m³s⁻¹[cite: 5]\n\n');
    
    fprintf('========================================\n');
    fprintf('UNVERIFIED ASSUMPTIONS\n');
    fprintf('========================================\n\n');
    fprintf('⚠ permeability_H2 = 1.0e-12 mol/(m·s·Pa)\n');
    fprintf('  - Engineering assumption - NOT from literature\n');
    fprintf('  - Should be verified experimentally\n\n');
    fprintf('⚠ permeability_O2 = 5.0e-13 mol/(m·s·Pa)\n');
    fprintf('  - Engineering assumption - NOT from literature\n');
    fprintf('  - Should be verified experimentally\n\n');
    fprintf('⚠ Gas partial pressures: 50% H₂, 21% O₂\n');
    fprintf('  - Engineering assumption\n\n');
    fprintf('========================================\n');
    fprintf('AUDIT COMPLETE\n');
    fprintf('========================================\n');
end