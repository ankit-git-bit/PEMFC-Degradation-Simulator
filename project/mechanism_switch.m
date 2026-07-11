function switches = mechanism_switch()
    % MECHANISM_SWITCH.M
    % Master control file to toggle individual degradation mechanisms.
    % This enables instantaneous isolation of specific pathways and sensitivity analysis.
    %
    % Usage:
    %   switches = mechanism_switch();
    %   % Then use switches in main simulation loop:
    %   if switches.peroxide
    %       rates = peroxide_generation(...);
    %   end
    %
    % Outputs:
    %   switches - Structure containing boolean flags for each mechanism

    %% 1. Chemical Degradation Mechanisms
    switches.peroxide  = true;   % Mechanism 1: Peroxide & radical generation
    switches.fenton    = true;   % Fenton chemistry (Fe2+/Fe3+ redox cycling)
    switches.sidechain = true;   % Mechanism 2: Side-chain attack and shedding
    switches.unzipping = true;   % Mechanism 3: Polymer backbone unzipping

    %% 2. Pt Kinetics
    switches.pt         = true;   % Pt dissolution/redeposition kinetics
    switches.pt_diss    = true;   % Pt dissolution specific toggle
    switches.pt_redep   = true;   % Pt redeposition specific toggle
    switches.pt_agglom  = true;   % Pt particle agglomeration

    %% 3. Physical Property Evolution
    switches.bidirectional = true;   % Bidirectional coupling between chemical and physical states
    switches.conductivity  = true;   % Conductivity degradation due to IEC loss
    switches.thickness     = true;   % Membrane thinning due to mass loss
    switches.IEC           = true;   % Ion Exchange Capacity degradation
    switches.water_uptake  = true;   % Water uptake changes with degradation
    switches.crossover     = true;   % Gas crossover evolution with degradation
    switches.porosity      = true;   % Porosity evolution
    switches.mechanical    = true;   % Mechanical property degradation

    %% 4. Performance Degradation
    switches.performance   = true;   % Performance loss due to degradation
    switches.ECSA_loss     = true;   % ECSA loss due to Pt dissolution
    switches.ohmic_loss    = true;   % Ohmic resistance increase
    switches.mass_transfer = true;   % Mass transfer limitations

    %% 5. Environmental Coupling
    switches.temperature   = true;   % Temperature dependence of all reactions
    switches.RH            = true;   % Relative humidity effects
    switches.pressure      = true;   % Pressure effects on crossover

    %% 6. Crossover Mechanisms
    switches.H2_crossover  = true;   % Hydrogen crossover
    switches.O2_crossover  = true;   % Oxygen crossover
    switches.H2O2_crossover= true;   % Peroxide crossover from cathode to anode

    %% 7. Membrane Degradation Pathways
    switches.radical_attack   = true;   % OH radical attack on membrane
    switches.chain_scission   = true;   % Chain scission reactions
    switches.crosslinking     = true;   % Crosslinking reactions (if applicable)
    switches.hydrolysis       = true;   % Hydrolysis reactions

    %% 8. Catalyst Degradation
    switches.catalyst_oxidation = true;   % Catalyst oxidation
    switches.catalyst_migration = true;   % Catalyst migration in membrane
    switches.catalyst_poisoning = true;   % Catalyst poisoning by impurities

    %% 9. ODE Solver Controls
    switches.adaptive_timestep = true;   % Adaptive timestep for stiffness
    switches.event_detection   = true;   % Event detection for end-of-life
    switches.sensitivity       = false;  % Sensitivity analysis mode (off by default)

    %% 10. Diagnostic Outputs
    switches.debug_output   = false;  % Enable debug printing
    switches.verbose        = true;   % Verbose console output
    switches.plot_progress  = false;  % Live plotting during simulation

    %% Quick Toggle Groups (For Sensitivity Studies)
    % Presets for common scenarios
    
    % Scenario 1: Fenton chemistry only (no Pt or membrane degradation)
    switches.scenario_fenton = false;
    if switches.scenario_fenton
        switches.peroxide  = true;
        switches.fenton    = true;
        switches.sidechain = false;
        switches.unzipping = false;
        switches.pt        = false;
        switches.bidirectional = false;
        switches.conductivity  = false;
    end
    
    % Scenario 2: Membrane degradation only (no Pt effects)
    switches.scenario_membrane = false;
    if switches.scenario_membrane
        switches.peroxide  = true;
        switches.fenton    = true;
        switches.sidechain = true;
        switches.unzipping = true;
        switches.pt        = false;
        switches.bidirectional = true;
        switches.conductivity  = true;
    end
    
    % Scenario 3: Pt degradation only (no chemical degradation)
    switches.scenario_pt_only = false;
    if switches.scenario_pt_only
        switches.peroxide  = false;
        switches.fenton    = false;
        switches.sidechain = false;
        switches.unzipping = false;
        switches.pt        = true;
        switches.bidirectional = false;
        switches.conductivity  = false;
    end
    
    % Scenario 4: Full degradation (all mechanisms active - DEFAULT)
    switches.scenario_full = true;
    if switches.scenario_full
        switches.peroxide  = true;
        switches.fenton    = true;
        switches.sidechain = true;
        switches.unzipping = true;
        switches.pt        = true;
        switches.bidirectional = true;
        switches.conductivity  = true;
        switches.thickness     = true;
        switches.IEC           = true;
        switches.water_uptake  = true;
        switches.crossover     = true;
        switches.performance   = true;
        switches.ECSA_loss     = true;
    end
    
    % Scenario 5: No degradation (baseline - all false)
    switches.scenario_baseline = false;
    if switches.scenario_baseline
        switches.peroxide  = false;
        switches.fenton    = false;
        switches.sidechain = false;
        switches.unzipping = false;
        switches.pt        = false;
        switches.bidirectional = false;
        switches.conductivity  = false;
        switches.thickness     = false;
        switches.IEC           = false;
        switches.water_uptake  = false;
        switches.crossover     = false;
        switches.performance   = false;
        switches.ECSA_loss     = false;
    end

    %% 11. Stiffness and Solver Settings
    % Controls for numerical stability
    switches.stiff_solver = true;    % Use stiff ODE solver (ode15s)
    switches.nonnegative  = true;    % Enforce non-negative species concentrations
    switches.max_timestep = 3600;    % Maximum timestep (seconds) - 1 hour
    
    %% 12. End-of-Life Detection
    switches.EOL_detection = true;   % Enable end-of-life detection
    switches.EOL_criteria = 'performance'; % 'performance', 'thickness', 'IEC', 'crossover'
    
    %% 13. Coupling Strength Controls
    % Adjust coupling strength between mechanisms
    switches.coupling_factor = 1.0;  % Multiplicative factor for coupling terms
    switches.feedback_loop   = true; % Enable feedback loops between mechanisms
    
    %% 14. Numerical Safety
    switches.clipping = true;        % Clip rates to prevent numerical instability
    switches.clipping_threshold = 1e-30; % Minimum concentration floor
    
    %% 15. Output Control
    switches.save_frequency = 3600;   % Save solution every N seconds
    switches.save_all       = true;   % Save all species or only selected
    
    %% 16. Model Validation
    switches.validation_mode = false; % Validate against experimental data
    switches.compare_exp     = false; % Compare with experimental results
    
    %% 17. User-Defined Customization
    % User can set custom toggles for specific applications
    switches.custom1 = true;
    switches.custom2 = true;
    switches.custom3 = false;
    
    %% 18. Mechanism Priority Order
    % Defines which mechanisms are computed first (affects coupling)
    switches.priority_order = {'peroxide', 'fenton', 'sidechain', 'unzipping', 'pt'};
    
    %% 19. Memory and Performance
    switches.use_parallel = false;   % Parallel computing if available
    switches.use_sparse   = true;    % Use sparse matrices for Jacobian
    switches.vectorized   = true;    % Use vectorized operations
    
    %% 20. Deprecated/Backward Compatibility
    % Maintain old field names for compatibility with existing code
    % Mechanism 1: Peroxide & radical generation (old name)
    if ~isfield(switches, 'peroxide')
        switches.peroxide = true;
    end
    
    % Mechanism 2: Side-chain attack (old name)
    if ~isfield(switches, 'attack') && isfield(switches, 'sidechain')
        switches.attack = switches.sidechain;
    end
    
    % Mechanism 3: Unzipping (old name)
    if ~isfield(switches, 'unzip') && isfield(switches, 'unzipping')
        switches.unzip = switches.unzipping;
    end
    
    %% 21. Display Summary (Optional)
    if switches.verbose
        fprintf('=== PEMFC Degradation Mechanism Switches ===\n');
        fprintf('Peroxide Generation:  %d\n', switches.peroxide);
        fprintf('Fenton Chemistry:     %d\n', switches.fenton);
        fprintf('Side-chain Attack:    %d\n', switches.sidechain);
        fprintf('Backbone Unzipping:   %d\n', switches.unzipping);
        fprintf('Pt Kinetics:          %d\n', switches.pt);
        fprintf('Conductivity Deg:     %d\n', switches.conductivity);
        fprintf('Bidirectional Coupling: %d\n', switches.bidirectional);
        fprintf('===========================================\n');
    end
end