function d_mech2 = membrane_attack(rates, species, kinetics, op)
    % MEMBRANE_ATTACK.M
    % Mechanism 2: Initial membrane attack and side-chain cleavage[cite: 5].
    % Uses ONLY reaction rates R14, R15, and R16.
    % Outputs the partial derivatives for membrane degradation species
    % and calculates IEC loss due to side-chain cleavage.
    %
    % Inputs:
    %   rates   - Structure containing reaction rates R14, R15, R16
    %   species - Structure containing species concentrations (SC_SO3H, SC_O, BB_O)
    %   kinetics- Structure containing kinetic parameters
    %   op      - Operating conditions structure
    %
    % Outputs:
    %   d_mech2 - Structure containing derivatives for:
    %             SC, SCO, BBO, HF, CO2, HOCF2CF2SO3H, and IEC

    %% 1. Membrane Degradation (Side-chain and Backbone Cleavage)[cite: 5]
    % Primary species derivatives
    d_mech2.SC   = -rates.R14;          % d[SC-SO3H]/dt - Side chain loss
    d_mech2.SCO  =  rates.R14 - rates.R15;  % d[SC-O]/dt - Side chain radical intermediate
    d_mech2.BBO  =  rates.R15 - rates.R16;  % d[BB-O]/dt - Backbone radical intermediate
    
    %% 2. Degradation Products Generation
    % HF generation from side-chain and backbone cleavage
    % R15 produces 6 HF molecules per reaction
    % R16 produces 3 HF molecules per reaction
    d_mech2.HF   = 6*rates.R15 + 3*rates.R16;
    
    % CO2 generation from backbone cleavage
    % R15 produces 3 CO2 molecules per reaction
    % R16 produces 1 CO2 per reaction (for initial backbone attack)
    d_mech2.CO2  = 3*rates.R15 + rates.R16;
    
    % HOCF2CF2SO3H (side-chain product) generation
    d_mech2.HOCF2CF2SO3H = rates.R14;
    
    %% 3. Ion Exchange Capacity (IEC) Loss
    % IEC loss is directly proportional to side-chain cleavage (R14)
    % Each SC-SO3H group lost reduces IEC by one equivalent
    % Conversion factor: 1 mol SC loss = 1 equivalent IEC loss
    % Scale by initial IEC for relative loss calculation
    
    % Absolute IEC loss rate (equivalents/m³·s)
    d_mech2.IEC_loss_abs = rates.R14;
    
    % Relative IEC loss rate (fraction of initial IEC per second)
    if isfield(species, 'IEC_0') && species.IEC_0 > 0
        d_mech2.IEC_loss_rel = rates.R14 / species.IEC_0;
    else
        % Default initial IEC if not provided (typical Nafion value)
        IEC_0_default = 1.0; % equivalents/kg
        d_mech2.IEC_loss_rel = rates.R14 / IEC_0_default;
    end
    
    % IEC decay rate constant (if using empirical degradation)
    if isfield(kinetics, 'membrane') && isfield(kinetics.membrane, 'k_IEC')
        d_mech2.IEC_empirical = -kinetics.membrane.k_IEC * species.SC_SO3H;
    else
        d_mech2.IEC_empirical = -rates.R14; % Rate-limited IEC loss
    end
    
    %% 4. Membrane Thickness Loss
    % Membrane thinning due to backbone and side-chain loss
    % Each mole of backbone loss results in a specific volume reduction
    % Molar volume of Nafion ~ 0.5 L/mol (approximate)
    MOLAR_VOLUME_NAFION = 0.0005; % m³/mol
    
    % Volume loss rate (m³/m³·s)
    d_mech2.volume_loss = MOLAR_VOLUME_NAFION * (rates.R15 + rates.R16);
    
    % Thickness loss rate (m/s) assuming uniform area
    if isfield(species, 'L_mem') && species.L_mem > 0
        d_mech2.thickness_loss = d_mech2.volume_loss * species.L_mem;
    else
        % Default membrane thickness
        L_mem_default = 1e-4; % 100 μm
        d_mech2.thickness_loss = d_mech2.volume_loss * L_mem_default;
    end
    
    %% 5. Diagnostic and Cumulative Metrics
    % Total side-chain degradation rate
    d_mech2.total_SC_loss = rates.R14;
    
    % Total backbone degradation rate
    d_mech2.total_BB_loss = rates.R15 + rates.R16;
    
    % Total membrane degradation rate (sum of all attack processes)
    d_mech2.total_membrane_loss = rates.R14 + rates.R15 + rates.R16;
    
    % Degradation fraction per time step
    if isfield(species, 'SC_0') && species.SC_0 > 0
        d_mech2.SC_loss_fraction = rates.R14 / species.SC_0;
    end
    
    %% 6. Coupling with Other Degradation Mechanisms
    % OH radical consumption in membrane attack
    d_mech2.OH_consumption = rates.R14 + 3*rates.R15 + rates.R16;
    
    % Radical attack efficiency (fraction of OH used for membrane attack)
    if isfield(rates, 'OH_total') && rates.OH_total > 0
        d_mech2.attack_efficiency = d_mech2.OH_consumption / rates.OH_total;
    else
        d_mech2.attack_efficiency = 0.5; % Default assumption
    end
    
    %% 7. IEC Loss Due to Side-Chain Cleavage
    % Calculate cumulative IEC loss from SC degradation
    % This is the main output for coupling with property evolution
    
    % IEC loss rate (equivalents/m³·s) - absolute value
    d_mech2.IEC_loss_rate = rates.R14;
    
    % IEC loss as percentage of initial IEC per second
    if isfield(species, 'IEC_0') && species.IEC_0 > 0
        d_mech2.IEC_loss_percent = (rates.R14 / species.IEC_0) * 100;
    else
        d_mech2.IEC_loss_percent = rates.R14 * 100; % Assuming IEC_0 = 1
    end
    
    % Normalized IEC derivative (dIEC/dt / IEC)
    if isfield(species, 'IEC') && species.IEC > 0
        d_mech2.IEC_normalized = -rates.R14 / species.IEC;
    end
    
    %% 8. Membrane Porosity Change
    % Porosity increases as membrane degrades
    % Loss of ionomer volume creates pores
    if isfield(species, 'porosity') && species.porosity >= 0
        d_mech2.porosity_change = d_mech2.volume_loss / (1 - species.porosity);
    else
        d_mech2.porosity_change = d_mech2.volume_loss;
    end
end