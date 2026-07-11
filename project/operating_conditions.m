function op = operating_conditions()
    % OPERATING_CONDITIONS.M
    % Defines the environment and load profile for the PEMFC simulation.
    % Pt dissolution and degradation are highly sensitive to these inputs.
    
    %% 1. TEMPERATURE & PRESSURE
    op.Temperature = 363.15;        % Operating Temperature (K) - updated from original
    op.pressure    = 2.3;           % Operational Pressure (bar) - updated from original
    
    %% 2. LOAD CYCLE (Dynamic Profiling)
    % Define the target load/current density
    op.current_density = 1000;      % Operating Current Density (A/m²)
    op.voltage         = 1.0;       % Operating Cell Voltage (V) - updated from original
    op.load_cycle      = 'steady';  % Cycle type: 'steady', 'dynamic', 'startup-shutdown'
    
    %% 3. RELATIVE HUMIDITY (RH)
    op.RH_anode   = 0.30;           % Anode side RH (fraction) - updated from original
    op.RH_cathode = 0.30;           % Cathode side RH (fraction) - updated from original
    
    %% 4. GAS COMPOSITION
    op.H2 = 0.010;                  % Anode H2 fraction - updated from original (mol/m³)
    op.O2 = 0.0075;                 % Cathode O2 fraction - updated from original (mol/m³)
    
    %% 5. DEGRADATION-SENSITIVE KINETICS
    % Pt dissolution is driven by local potential (Voltage) and T
    % Redeposition is driven by Pt2+ concentration and local RH
    op.V_ocv = 0.95;                % Open Circuit Voltage (V)
    op.i0_ref = 1e-4;               % Reference current density (A/m²)
end