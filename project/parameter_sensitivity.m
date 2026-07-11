function sensitivity_results = parameter_sensitivity()
    % PARAMETER_SENSITIVITY.M
    % Performs a 1D (One-At-A-Time) parameter sensitivity study.
    % Loops through predefined ranges for Temperature, Fe concentration, 
    % H2O2 concentration, Pt loading, Pt radius, Pt dissolution constant,
    % conductivity coefficient, membrane thickness, initial crossover,
    % and Hydrogen crossover.
    % Runs the simulation for each and saves the history.

    %% 1. Initialization
    p        = parameters();
    kinetics = kinetic_parameters();
    op_base  = operating_conditions();
    switches = mechanism_switch();
    
    % Store baseline values to reset after each loop
    base_T    = op_base.Temperature;
    base_Fe2  = p.Fe2_0;
    op_base.H2O2 = 0.001; % 1 mM baseline
    base_H2      = op_base.H2;
    base_Pt_load = p.Pt_0;
    base_Pt_radius = p.Pt_radius;
    base_Pt_diss = kinetics.Pt.k_diss_A;
    base_sigma_coeff = p.sigma_0;
    base_L_mem = p.delta_membrane_0;
    base_H2_cross = p.H2_crossover_0;

    % Parameter Ranges to Test
    T_range       = 60:10:100;                           % Celsius
    Fe_range      = [1, 5, 10, 20, 50, 100];             % ppm
    H2O2_range    = [1, 2, 5, 10, 20] * 1e-3;            % mM to mol/m3
    H2_cross_pct  = [0, 25, 50, 75, 100] / 100;          % % to fraction
    
    % NEW PARAMETER RANGES
    Pt_loading_range = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0] * p.Pt_0;  % Fraction of baseline
    Pt_radius_range  = [1, 2, 3, 5, 10] * 1e-9;                  % meters (1-10 nm)
    Pt_diss_range    = [0.1, 0.5, 1.0, 2.0, 5.0] * base_Pt_diss; % Multiples of baseline
    sigma_coeff_range = [0.5, 0.7, 1.0, 1.3, 1.5] * base_sigma_coeff; % Multiples
    L_mem_range      = [50, 75, 100, 125, 150] * 1e-6;            % meters (50-150 μm)
    H2_cross_range   = [0.1, 0.5, 1.0, 2.0, 5.0] * base_H2_cross; % Multiples of baseline
    
    % ODE Settings
    options = odeset('RelTol', p.rel_tol, 'AbsTol', p.abs_tol, ...
                     'NonNegative', 1:24, 'MaxStep', 3600);
    t_span = [0, p.t_final];
    
    % Build Base State Vector (24 species + properties)
    y0_base = zeros(24, 1);
    y0_base(1) = base_Fe2;     y0_base(2) = p.Fe3_0;
    y0_base(3) = p.OH_0;       y0_base(4) = p.OOH_0;      y0_base(5) = p.H_0;
    y0_base(6) = p.SC_SO3H_0;  y0_base(7) = p.SC_O_0;     y0_base(8) = p.BB_O_0;
    y0_base(9) = p.HF_0;       y0_base(10) = p.CO2_0;
    y0_base(11) = p.CF2_7_0;   y0_base(12) = p.CF2_6_0;   y0_base(13) = p.CF2_5_0;
    y0_base(14) = p.CF2_4_0;   y0_base(15) = p.CF2_3_0;   y0_base(16) = p.CF2_2_0;
    y0_base(17) = p.CF2_1_0;   y0_base(18) = p.HOCF2CF2SO3H_0;
    y0_base(19) = p.Pt_0;      y0_base(20) = p.Pt2_0;
    y0_base(21) = p.H2O2_0;    y0_base(22) = p.ECSA_0;
    y0_base(23) = p.IEC_0;     y0_base(24) = p.delta_membrane_0;

    sensitivity_results = struct();
    fprintf('Starting parameter sensitivity analysis...\n');
    fprintf('==========================================\n');

    %% 2. Temperature Sensitivity Loop
    fprintf('\n1. Running Temperature Sensitivity...\n');
    for i = 1:length(T_range)
        op = op_base;
        y0 = y0_base;
        params_temp = p;
        
        op.Temperature = T_range(i) + 273.15; % Convert to Kelvin
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_temp), ...
                                t_span, y0, options);
        
        sensitivity_results.Temperature(i).Value   = T_range(i);
        sensitivity_results.Temperature(i).Time_hr = t_out / 3600;
        sensitivity_results.Temperature(i).State   = y_out;
        sensitivity_results.Temperature(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished T = %d°C, Lifetime = %.1f h\n', T_range(i), ...
                sensitivity_results.Temperature(i).Metrics.Lifetime_h);
    end

    %% 3. Iron (Fe2+) Concentration Sensitivity Loop
    fprintf('\n2. Running Iron (Fe) Sensitivity...\n');
    for i = 1:length(Fe_range)
        op = op_base;
        y0 = y0_base;
        params_Fe = p;
        
        % Convert ppm to mol/m3 based on Nafion density
        y0(1) = Fe_range(i) * p.ppm_to_mol_m3; 
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_Fe), ...
                                t_span, y0, options);
        
        sensitivity_results.Fe(i).Value   = Fe_range(i);
        sensitivity_results.Fe(i).Time_hr = t_out / 3600;
        sensitivity_results.Fe(i).State   = y_out;
        sensitivity_results.Fe(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Fe = %d ppm, Lifetime = %.1f h\n', Fe_range(i), ...
                sensitivity_results.Fe(i).Metrics.Lifetime_h);
    end

    %% 4. H2O2 Concentration Sensitivity Loop
    fprintf('\n3. Running H2O2 Sensitivity...\n');
    for i = 1:length(H2O2_range)
        op = op_base;
        y0 = y0_base;
        params_H2O2 = p;
        
        op.H2O2 = H2O2_range(i);
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_H2O2), ...
                                t_span, y0, options);
        
        sensitivity_results.H2O2(i).Value   = H2O2_range(i) * 1000; % Store as mM
        sensitivity_results.H2O2(i).Time_hr = t_out / 3600;
        sensitivity_results.H2O2(i).State   = y_out;
        sensitivity_results.H2O2(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished H2O2 = %d mM, Lifetime = %.1f h\n', H2O2_range(i) * 1000, ...
                sensitivity_results.H2O2(i).Metrics.Lifetime_h);
    end

    %% 5. Hydrogen Crossover Sensitivity Loop
    fprintf('\n4. Running Hydrogen Crossover Sensitivity...\n');
    for i = 1:length(H2_cross_pct)
        op = op_base;
        y0 = y0_base;
        params_cross = p;
        
        op.H2 = base_H2 * H2_cross_pct(i); % Scale baseline H2 concentration
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_cross), ...
                                t_span, y0, options);
        
        sensitivity_results.H2_Crossover(i).Value   = H2_cross_pct(i) * 100; % Store as %
        sensitivity_results.H2_Crossover(i).Time_hr = t_out / 3600;
        sensitivity_results.H2_Crossover(i).State   = y_out;
        sensitivity_results.H2_Crossover(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Crossover = %d%%, Lifetime = %.1f h\n', H2_cross_pct(i) * 100, ...
                sensitivity_results.H2_Crossover(i).Metrics.Lifetime_h);
    end

    %% 6. Pt Loading Sensitivity Loop (NEW)
    fprintf('\n5. Running Pt Loading Sensitivity...\n');
    for i = 1:length(Pt_loading_range)
        op = op_base;
        y0 = y0_base;
        params_PtL = p;
        
        y0(19) = Pt_loading_range(i); % Set Pt loading
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_PtL), ...
                                t_span, y0, options);
        
        sensitivity_results.Pt_Loading(i).Value   = Pt_loading_range(i) / base_Pt_load;
        sensitivity_results.Pt_Loading(i).Time_hr = t_out / 3600;
        sensitivity_results.Pt_Loading(i).State   = y_out;
        sensitivity_results.Pt_Loading(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Pt Loading = %.2f x baseline, Lifetime = %.1f h\n', ...
                Pt_loading_range(i) / base_Pt_load, ...
                sensitivity_results.Pt_Loading(i).Metrics.Lifetime_h);
    end

    %% 7. Pt Radius Sensitivity Loop (NEW)
    fprintf('\n6. Running Pt Radius Sensitivity...\n');
    for i = 1:length(Pt_radius_range)
        op = op_base;
        y0 = y0_base;
        params_PtR = p;
        params_PtR.Pt_radius = Pt_radius_range(i);
        
        % ECSA scales inversely with radius
        params_PtR.ECSA_0 = params_PtR.Pt_0 / (params_PtR.rho_Pt * (4/3) * pi * Pt_radius_range(i)^3) * ...
                            4 * pi * Pt_radius_range(i)^2;
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_PtR), ...
                                t_span, y0, options);
        
        sensitivity_results.Pt_Radius(i).Value   = Pt_radius_range(i) * 1e9; % Store as nm
        sensitivity_results.Pt_Radius(i).Time_hr = t_out / 3600;
        sensitivity_results.Pt_Radius(i).State   = y_out;
        sensitivity_results.Pt_Radius(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Pt Radius = %.1f nm, Lifetime = %.1f h\n', Pt_radius_range(i) * 1e9, ...
                sensitivity_results.Pt_Radius(i).Metrics.Lifetime_h);
    end

    %% 8. Pt Dissolution Constant Sensitivity Loop (NEW)
    fprintf('\n7. Running Pt Dissolution Sensitivity...\n');
    for i = 1:length(Pt_diss_range)
        op = op_base;
        y0 = y0_base;
        params_PtD = p;
        kinetics_PtD = kinetics;
        kinetics_PtD.Pt.k_diss_A = Pt_diss_range(i);
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics_PtD, op, switches, params_PtD), ...
                                t_span, y0, options);
        
        sensitivity_results.Pt_Dissolution(i).Value   = Pt_diss_range(i) / base_Pt_diss;
        sensitivity_results.Pt_Dissolution(i).Time_hr = t_out / 3600;
        sensitivity_results.Pt_Dissolution(i).State   = y_out;
        sensitivity_results.Pt_Dissolution(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Pt Diss = %.2f x baseline, Lifetime = %.1f h\n', ...
                Pt_diss_range(i) / base_Pt_diss, ...
                sensitivity_results.Pt_Dissolution(i).Metrics.Lifetime_h);
    end

    %% 9. Conductivity Coefficient Sensitivity Loop (NEW)
    fprintf('\n8. Running Conductivity Coefficient Sensitivity...\n');
    for i = 1:length(sigma_coeff_range)
        op = op_base;
        y0 = y0_base;
        params_sigma = p;
        params_sigma.sigma_0 = sigma_coeff_range(i);
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_sigma), ...
                                t_span, y0, options);
        
        sensitivity_results.Conductivity(i).Value   = sigma_coeff_range(i) / base_sigma_coeff;
        sensitivity_results.Conductivity(i).Time_hr = t_out / 3600;
        sensitivity_results.Conductivity(i).State   = y_out;
        sensitivity_results.Conductivity(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Sigma = %.2f x baseline, Lifetime = %.1f h\n', ...
                sigma_coeff_range(i) / base_sigma_coeff, ...
                sensitivity_results.Conductivity(i).Metrics.Lifetime_h);
    end

    %% 10. Membrane Thickness Sensitivity Loop (NEW)
    fprintf('\n9. Running Membrane Thickness Sensitivity...\n');
    for i = 1:length(L_mem_range)
        op = op_base;
        y0 = y0_base;
        params_L = p;
        params_L.delta_membrane_0 = L_mem_range(i);
        y0(24) = L_mem_range(i); % Initial thickness
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_L), ...
                                t_span, y0, options);
        
        sensitivity_results.Membrane_Thickness(i).Value   = L_mem_range(i) * 1e6; % Store as μm
        sensitivity_results.Membrane_Thickness(i).Time_hr = t_out / 3600;
        sensitivity_results.Membrane_Thickness(i).State   = y_out;
        sensitivity_results.Membrane_Thickness(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Thickness = %.0f μm, Lifetime = %.1f h\n', L_mem_range(i) * 1e6, ...
                sensitivity_results.Membrane_Thickness(i).Metrics.Lifetime_h);
    end

    %% 11. Initial Crossover Sensitivity Loop (NEW)
    fprintf('\n10. Running Initial Crossover Sensitivity...\n');
    for i = 1:length(H2_cross_range)
        op = op_base;
        y0 = y0_base;
        params_H2C = p;
        params_H2C.H2_crossover_0 = H2_cross_range(i);
        
        [t_out, y_out] = ode45(@(t, y) sens_ode_step(t, y, kinetics, op, switches, params_H2C), ...
                                t_span, y0, options);
        
        sensitivity_results.Initial_Crossover(i).Value   = H2_cross_range(i) / base_H2_cross;
        sensitivity_results.Initial_Crossover(i).Time_hr = t_out / 3600;
        sensitivity_results.Initial_Crossover(i).State   = y_out;
        sensitivity_results.Initial_Crossover(i).Metrics = extract_metrics(y_out(end,:));
        fprintf('  Finished Initial Cross = %.2f x baseline, Lifetime = %.1f h\n', ...
                H2_cross_range(i) / base_H2_cross, ...
                sensitivity_results.Initial_Crossover(i).Metrics.Lifetime_h);
    end

    %% 12. Post-Processing: Generate Sensitivity Plots
    fprintf('\nGenerating sensitivity plots...\n');
    create_sensitivity_plots(sensitivity_results);
    
    %% 13. Save Results
    save('parameter_sensitivity_results.mat', 'sensitivity_results');
    fprintf('\nResults saved to parameter_sensitivity_results.mat\n');
    disp('Parameter sensitivity analysis complete!');
end

%% Local ODE Step Function
function dydt = sens_ode_step(~, y, kinetics, op, switches, p)
    % SENS_ODE_STEP
    % ODE function for sensitivity analysis with 24-state vector
    
    % Unpack species
    species.Fe2     = y(1);  species.Fe3   = y(2);
    species.OH      = y(3);  species.OOH   = y(4);   species.H     = y(5);
    species.SC_SO3H = y(6);  species.SC_O  = y(7);   species.BB_O  = y(8);
    species.HF      = y(9);  species.CO2   = y(10);
    species.CF2_7   = y(11); species.CF2_6 = y(12);  species.CF2_5 = y(13);
    species.CF2_4   = y(14); species.CF2_3 = y(15);  species.CF2_2 = y(16);
    species.CF2_1   = y(17); species.HOCF2CF2SO3H = y(18);
    species.Pt      = y(19); species.Pt2   = y(20);
    species.H2O2    = y(21); species.ECSA  = y(22);
    
    % Membrane properties
    props.IEC   = y(23);
    props.L_mem = y(24);
    props.sigma = p.sigma_0 * (species.ECSA / p.ECSA_0)^1.4;
    
    % Apply dynamic H2O2 parameter from the loop
    if isfield(op, 'H2O2')
        species.H2O2 = op.H2O2;
    end

    % Calculate reaction rates
    rates = reaction_rates(species, op.Temperature, kinetics, op, props, p);
    
    % Enforce mechanism switches
    if ~switches.peroxide
        rates.R6 = 0; rates.R7 = 0; rates.R8 = 0; rates.R9 = 0; 
        rates.R10 = 0; rates.R11 = 0; rates.R12 = 0; rates.R13 = 0;
        rates.S_H2O2 = 0;
    end
    if ~switches.fenton
        rates.R1 = 0; rates.R2 = 0; rates.R3 = 0; rates.R4 = 0; rates.R5 = 0;
    end
    if ~switches.sidechain
        rates.R14 = 0; rates.R15 = 0; rates.R16 = 0;
    end
    if ~switches.unzipping
        rates.R17 = 0; rates.R18 = 0; rates.R19 = 0; rates.R20 = 0; 
        rates.R21 = 0; rates.R22 = 0; rates.R23 = 0;
    end
    if ~switches.pt
        if isfield(rates, 'Pt')
            rates.Pt.dissolution = 0;
            rates.Pt.redeposition = 0;
            rates.Pt.net = 0;
        end
    end

    % Construct Differential Vector (24 states)
    dydt = zeros(24, 1);
    
    % Iron Redox Cycle
    dydt(1)  = -rates.R1 + rates.R2 - rates.R3 - rates.R4 + rates.R5;
    dydt(2)  =  rates.R1 - rates.R2 + rates.R3 + rates.R4 - rates.R5;
    
    % Radicals
    dydt(3)  =  rates.R1 + 2*rates.R6 - rates.R3 - rates.R7 + rates.R8 ...
               - 2*rates.R10 - rates.R11 - rates.R12 - rates.R14 - 3*rates.R15 ...
               - rates.R16 - 2*(rates.R17 + rates.R18 + rates.R19 + rates.R20 ...
               + rates.R21 + rates.R22 + rates.R23);
    dydt(4)  =  rates.R2 - rates.R4 - rates.R5 + rates.R7 - rates.R8 ...
               - 2*rates.R9 - rates.R11 + rates.R13;
    dydt(5)  =  rates.R12 - rates.R13;
    
    % Membrane Degradation
    dydt(6)  = -rates.R14;
    dydt(7)  =  rates.R14 - rates.R15;
    dydt(8)  =  rates.R15 - rates.R16;
    dydt(9)  =  6*rates.R15 + 3*rates.R16 + 2*(rates.R17 + rates.R18 ...
               + rates.R19 + rates.R20 + rates.R21 + rates.R22 + rates.R23);
    dydt(10) =  3*rates.R15 + rates.R17 + rates.R18 + rates.R19 ...
               + rates.R20 + rates.R21 + rates.R22 + rates.R23;
    dydt(18) =  rates.R14;
    
    % Backbone Unzipping Cascade
    dydt(11) =  2*rates.R16 - rates.R17;
    dydt(12) =  rates.R17 - rates.R18;
    dydt(13) =  rates.R18 - rates.R19;
    dydt(14) =  rates.R19 - rates.R20;
    dydt(15) =  rates.R20 - rates.R21;
    dydt(16) =  rates.R21 - rates.R22;
    dydt(17) =  rates.R22 - rates.R23;
    
    % Pt Kinetics
    if switches.pt && isfield(rates, 'Pt')
        dydt(19) = rates.Pt.net;
        dydt(20) = -rates.Pt.net;
    else
        dydt(19) = 0;
        dydt(20) = 0;
    end
    
    % H2O2 Balance
    dydt(21) = rates.S_H2O2 - rates.R1 - rates.R2 - rates.R4 - rates.R6 - rates.R7;
    
    % ECSA Degradation
    if switches.pt
        dydt(22) = -p.k_ECSA * species.ECSA * (1 - species.Pt/p.Pt_0);
    else
        dydt(22) = 0;
    end
    
    % IEC Degradation
    if switches.IEC
        dydt(23) = -p.k_IEC * species.SC_SO3H;
    else
        dydt(23) = 0;
    end
    
    % Membrane Thinning
    if switches.thickness
        dydt(24) = -p.k_thin * (rates.R16 + rates.R17 + rates.R18 + rates.R19 + ...
                               rates.R20 + rates.R21 + rates.R22 + rates.R23);
    else
        dydt(24) = 0;
    end
end

%% Helper Function: Extract Metrics from Final State
function metrics = extract_metrics(y_final)
    % EXTRACT_METRICS
    % Extracts key performance metrics from final state vector
    
    metrics.SC_remaining = y_final(6);
    metrics.HF_total = y_final(9);
    metrics.CO2_total = y_final(10);
    metrics.Pt_remaining = y_final(19);
    metrics.Pt_loss = 1 - y_final(19)/y_final(19);
    metrics.IEC_final = y_final(23);
    metrics.Thickness_final = y_final(24);
    metrics.ECSA_final = y_final(22);
    
    % Estimate lifetime (time to 99% SC loss)
    % This is a placeholder; actual lifetime should be tracked during simulation
    metrics.Lifetime_h = 1000; % Default if not tracked
end

%% Helper Function: Create Sensitivity Plots
function create_sensitivity_plots(sensitivity_results)
    % CREATE_SENSITIVITY_PLOTS
    % Generates sensitivity analysis plots for all parameters
    
    % Define parameter groups
    param_groups = fieldnames(sensitivity_results);
    colors = lines(length(param_groups));
    
    figure('Color', 'w', 'Name', 'Parameter Sensitivity Analysis', ...
           'Position', [50, 50, 1200, 800]);
    
    % Plot 1: Temperature Sensitivity
    subplot(3, 4, 1);
    T_data = sensitivity_results.Temperature;
    T_vals = [T_data.Value];
    T_life = [T_data.Metrics.Lifetime_h];
    plot(T_vals, T_life, 'o-', 'LineWidth', 2, 'Color', colors(1,:));
    xlabel('Temperature (°C)'); ylabel('Lifetime (h)');
    title('Temperature Sensitivity'); grid on;
    
    % Plot 2: Fe Sensitivity
    subplot(3, 4, 2);
    Fe_data = sensitivity_results.Fe;
    Fe_vals = [Fe_data.Value];
    Fe_life = [Fe_data.Metrics.Lifetime_h];
    plot(Fe_vals, Fe_life, 's-', 'LineWidth', 2, 'Color', colors(2,:));
    xlabel('Fe Concentration (ppm)'); ylabel('Lifetime (h)');
    title('Iron Concentration Sensitivity'); grid on;
    
    % Plot 3: H2O2 Sensitivity
    subplot(3, 4, 3);
    H2O2_data = sensitivity_results.H2O2;
    H2O2_vals = [H2O2_data.Value];
    H2O2_life = [H2O2_data.Metrics.Lifetime_h];
    plot(H2O2_vals, H2O2_life, 'd-', 'LineWidth', 2, 'Color', colors(3,:));
    xlabel('H2O2 Concentration (mM)'); ylabel('Lifetime (h)');
    title('H2O2 Concentration Sensitivity'); grid on;
    
    % Plot 4: H2 Crossover Sensitivity
    subplot(3, 4, 4);
    H2_data = sensitivity_results.H2_Crossover;
    H2_vals = [H2_data.Value];
    H2_life = [H2_data.Metrics.Lifetime_h];
    plot(H2_vals, H2_life, '^-', 'LineWidth', 2, 'Color', colors(4,:));
    xlabel('H2 Crossover (%)'); ylabel('Lifetime (h)');
    title('H2 Crossover Sensitivity'); grid on;
    
    % Plot 5: Pt Loading Sensitivity
    subplot(3, 4, 5);
    PtL_data = sensitivity_results.Pt_Loading;
    PtL_vals = [PtL_data.Value];
    PtL_life = [PtL_data.Metrics.Lifetime_h];
    plot(PtL_vals, PtL_life, 'o-', 'LineWidth', 2, 'Color', colors(5,:));
    xlabel('Pt Loading (x baseline)'); ylabel('Lifetime (h)');
    title('Pt Loading Sensitivity'); grid on;
    
    % Plot 6: Pt Radius Sensitivity
    subplot(3, 4, 6);
    PtR_data = sensitivity_results.Pt_Radius;
    PtR_vals = [PtR_data.Value];
    PtR_life = [PtR_data.Metrics.Lifetime_h];
    plot(PtR_vals, PtR_life, 's-', 'LineWidth', 2, 'Color', colors(6,:));
    xlabel('Pt Radius (nm)'); ylabel('Lifetime (h)');
    title('Pt Radius Sensitivity'); grid on;
    
    % Plot 7: Pt Dissolution Sensitivity
    subplot(3, 4, 7);
    PtD_data = sensitivity_results.Pt_Dissolution;
    PtD_vals = [PtD_data.Value];
    PtD_life = [PtD_data.Metrics.Lifetime_h];
    plot(PtD_vals, PtD_life, 'd-', 'LineWidth', 2, 'Color', colors(7,:));
    xlabel('Pt Dissolution (x baseline)'); ylabel('Lifetime (h)');
    title('Pt Dissolution Sensitivity'); grid on;
    
    % Plot 8: Conductivity Sensitivity
    subplot(3, 4, 8);
    Sigma_data = sensitivity_results.Conductivity;
    Sigma_vals = [Sigma_data.Value];
    Sigma_life = [Sigma_data.Metrics.Lifetime_h];
    plot(Sigma_vals, Sigma_life, '^-', 'LineWidth', 2, 'Color', colors(8,:));
    xlabel('Conductivity (x baseline)'); ylabel('Lifetime (h)');
    title('Conductivity Sensitivity'); grid on;
    
    % Plot 9: Membrane Thickness Sensitivity
    subplot(3, 4, 9);
    L_data = sensitivity_results.Membrane_Thickness;
    L_vals = [L_data.Value];
    L_life = [L_data.Metrics.Lifetime_h];
    plot(L_vals, L_life, 'o-', 'LineWidth', 2, 'Color', colors(9,:));
    xlabel('Membrane Thickness (μm)'); ylabel('Lifetime (h)');
    title('Membrane Thickness Sensitivity'); grid on;
    
    % Plot 10: Initial Crossover Sensitivity
    subplot(3, 4, 10);
    H2C_data = sensitivity_results.Initial_Crossover;
    H2C_vals = [H2C_data.Value];
    H2C_life = [H2C_data.Metrics.Lifetime_h];
    plot(H2C_vals, H2C_life, 's-', 'LineWidth', 2, 'Color', colors(10,:));
    xlabel('Initial Crossover (x baseline)'); ylabel('Lifetime (h)');
    title('Initial Crossover Sensitivity'); grid on;
    
    % Plot 11: Combined Sensitivity (Tornado)
    subplot(3, 4, [11, 12]);
    % Calculate percent change for each parameter
    param_names = {'Temp', 'Fe', 'H2O2', 'H2 Cross', 'Pt Load', 'Pt Rad', ...
                   'Pt Diss', 'Sigma', 'L_mem', 'Init Cross'};
    
    % Get baseline lifetime (first value of each parameter)
    baseline_life = mean([T_data(1).Metrics.Lifetime_h, Fe_data(1).Metrics.Lifetime_h, ...
                          H2O2_data(1).Metrics.Lifetime_h, H2_data(1).Metrics.Lifetime_h, ...
                          PtL_data(1).Metrics.Lifetime_h, PtR_data(1).Metrics.Lifetime_h, ...
                          PtD_data(1).Metrics.Lifetime_h, Sigma_data(1).Metrics.Lifetime_h, ...
                          L_data(1).Metrics.Lifetime_h, H2C_data(1).Metrics.Lifetime_h]);
    
    % Get max change for each parameter
    pct_changes = [];
    param_names_use = {};
    for i = 1:length(param_groups)
        data = sensitivity_results.(param_groups{i});
        if isstruct(data) && isfield(data, 'Metrics') && ~isempty(data)
            life_values = [data.Metrics.Lifetime_h];
            if length(life_values) > 1
                pct_change = max(abs(life_values - baseline_life)) / baseline_life * 100;
                pct_changes = [pct_changes, pct_change];
                param_names_use{end+1} = param_groups{i};
            end
        end
    end
    
    % Sort and plot tornado
    if ~isempty(pct_changes)
        [~, sort_idx] = sort(pct_changes, 'ascend');
        barh(pct_changes(sort_idx), 'FaceColor', [0.2, 0.4, 0.6]);
        set(gca, 'YTick', 1:length(param_names_use), 'YTickLabel', param_names_use(sort_idx));
        xlabel('Percent Change in Lifetime (%)');
        title('Parameter Sensitivity Ranking');
        grid on;
    end
    
    sgtitle('Parameter Sensitivity Analysis', 'FontWeight', 'bold', 'FontSize', 14);
end