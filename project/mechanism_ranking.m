function results = mechanism_ranking()
    % MECHANISM_RANKING.M
    % Loops through 7 mechanisms, disables them one by one, runs the ODE,
    % measures key metrics (Lifetime, HF, FER, SC, IEC, Conductivity, Pt loss),
    % and outputs a Tornado chart and comprehensive results table.
    %
    % Mechanisms evaluated:
    %   1. Baseline (all active)
    %   2. Peroxide Generation (R6-R13)
    %   3. Fenton Chemistry (R1-R5)
    %   4. Pt Kinetics (dissolution/redeposition)
    %   5. Bidirectional Coupling (chemical-physical feedback)
    %   6. Conductivity Degradation
    %   7. Unzipping (R17-R23)

    %% 1. Initialization
    p        = parameters();
    kinetics = kinetic_parameters();
    op       = operating_conditions();
    
    % Define the 7 mechanisms to evaluate
    mechs = {'baseline', 'peroxide', 'fenton', 'pt', 'bidirectional', 'conductivity', 'unzipping'};
    num_runs = length(mechs);
    
    % Storage arrays for metrics
    m_Lifetime = zeros(num_runs, 1);
    m_HF       = zeros(num_runs, 1);
    m_FER      = zeros(num_runs, 1);
    m_SC       = zeros(num_runs, 1);
    m_IEC      = zeros(num_runs, 1);
    m_Conductivity = zeros(num_runs, 1);
    m_Pt_loss  = zeros(num_runs, 1);
    m_CO2      = zeros(num_runs, 1);
    m_Thickness = zeros(num_runs, 1);
    m_ECSA     = zeros(num_runs, 1);

    % ODE Solver Settings
    options = odeset('RelTol', p.rel_tol, 'AbsTol', p.abs_tol, ...
                     'NonNegative', 1:24, 'MaxStep', 3600);
    t_span = [0, p.t_final];
    
    % Initial State Vector (24 species + properties)
    y0 = zeros(24, 1);
    y0(1) = p.Fe2_0;          y0(2) = p.Fe3_0;
    y0(3) = p.OH_0;           y0(4) = p.OOH_0;        y0(5) = p.H_0;
    y0(6) = p.SC_SO3H_0;      y0(7) = p.SC_O_0;       y0(8) = p.BB_O_0;
    y0(9) = p.HF_0;           y0(10) = p.CO2_0;
    y0(11) = p.CF2_7_0;       y0(12) = p.CF2_6_0;     y0(13) = p.CF2_5_0;
    y0(14) = p.CF2_4_0;       y0(15) = p.CF2_3_0;     y0(16) = p.CF2_2_0;
    y0(17) = p.CF2_1_0;       y0(18) = p.HOCF2CF2SO3H_0;
    y0(19) = p.Pt_0;          y0(20) = p.Pt2_0;
    y0(21) = p.H2O2_0;        y0(22) = p.ECSA_0;
    y0(23) = p.IEC_0;         y0(24) = p.delta_membrane_0;

    disp('Starting mechanism ranking loop...');
    disp('=======================================');

    %% 2. Loop Through Mechanisms
    for i = 1:num_runs
        current_run = mechs{i};
        fprintf('Running scenario: %s disabled...\n', current_run);
        
        % Reset all switches to true
        switches = mechanism_switch();
        switches.peroxide  = true;
        switches.fenton    = true;
        switches.pt        = true;
        switches.bidirectional = true;
        switches.conductivity  = true;
        switches.unzipping = true;
        switches.sidechain = true; % Keep sidechain active for all runs
        switches.thickness  = true;
        switches.IEC        = true;
        switches.water_uptake = true;
        switches.crossover  = true;
        
        % Disable the targeted mechanism (skip for baseline)
        if i > 1
            switches.(current_run) = false;
            
            % Special handling for bidirectional coupling
            if strcmp(current_run, 'bidirectional')
                switches.conductivity = false;
                switches.thickness = false;
                switches.IEC = false;
                switches.water_uptake = false;
                switches.crossover = false;
            end
            
            % Special handling for conductivity
            if strcmp(current_run, 'conductivity')
                switches.IEC = false;
                switches.water_uptake = false;
            end
        end
        
        % Run ODE with current switches
        try
            [t_out, y_out] = ode45(@(t, y) ranked_ode_step(t, y, kinetics, op, p, switches), ...
                                   t_span, y0, options);
        catch ME
            warning('ODE solver failed for %s: %s', current_run, ME.message);
            t_out = t_span(2);
            y_out = y0';
            continue;
        end
        
        % Extract Metrics
        % SC: Remaining Side Chain at the end
        SC_final = y_out(end, 6);
        m_SC(i) = SC_final;
        
        % HF: Total HF accumulated at the end
        HF_final = y_out(end, 9);
        m_HF(i) = HF_final;
        
        % CO2: Total CO2 produced
        CO2_final = y_out(end, 10);
        m_CO2(i) = CO2_final;
        
        % Pt loss: Fraction of Pt lost
        Pt_final = y_out(end, 19);
        m_Pt_loss(i) = 1 - (Pt_final / p.Pt_0);
        
        % IEC: Final IEC value
        if size(y_out, 2) >= 23
            IEC_final = y_out(end, 23);
            m_IEC(i) = IEC_final;
        else
            m_IEC(i) = 0;
        end
        
        % Conductivity: Final conductivity (estimated from IEC)
        if size(y_out, 2) >= 23
            acid_fraction = y_out(end, 23) / p.IEC_0;
            m_Conductivity(i) = p.sigma_0 * (acid_fraction)^(1/1.4);
        else
            m_Conductivity(i) = 0;
        end
        
        % Membrane thickness: Final thickness
        if size(y_out, 2) >= 24
            m_Thickness(i) = y_out(end, 24);
        else
            m_Thickness(i) = p.delta_membrane_0;
        end
        
        % ECSA: Final ECSA
        if size(y_out, 2) >= 22
            m_ECSA(i) = y_out(end, 22);
        else
            m_ECSA(i) = p.ECSA_0;
        end
        
        % Lifetime: Time to 99% SC loss (or end of simulation if not reached)[cite: 5]
        SC_threshold = 0.01 * p.SC_SO3H_0;
        idx_fail = find(y_out(:, 6) <= SC_threshold, 1);
        if isempty(idx_fail)
            lifetime_s = t_out(end);
        else
            lifetime_s = t_out(idx_fail);
        end
        m_Lifetime(i) = lifetime_s / 3600; % Convert to hours
        
        % FER: Fluoride Emission Rate (HF / Lifetime)[cite: 5]
        if m_Lifetime(i) > 0
            m_FER(i) = HF_final / m_Lifetime(i);
        else
            m_FER(i) = 0;
        end
        
        fprintf('  Lifetime: %.1f h, HF: %.3f mol/m³, SC: %.3f mol/m³\n', ...
                m_Lifetime(i), HF_final, SC_final);
    end

    %% 3. Store Results in Struct
    results.Mechanisms = mechs(2:end);
    results.Baseline.Lifetime = m_Lifetime(1);
    results.Baseline.HF       = m_HF(1);
    results.Baseline.FER      = m_FER(1);
    results.Baseline.SC       = m_SC(1);
    results.Baseline.IEC      = m_IEC(1);
    results.Baseline.Conductivity = m_Conductivity(1);
    results.Baseline.Pt_loss  = m_Pt_loss(1);
    results.Baseline.CO2      = m_CO2(1);
    results.Baseline.Thickness = m_Thickness(1);
    results.Baseline.ECSA     = m_ECSA(1);
    
    % Store all results in table format
    results_table = table();
    results_table.Mechanism = mechs(2:end)';
    results_table.Lifetime_h = m_Lifetime(2:end);
    results_table.HF_mol_m3 = m_HF(2:end);
    results_table.FER_mol_m3h = m_FER(2:end);
    results_table.SC_mol_m3 = m_SC(2:end);
    results_table.IEC_meq_g = m_IEC(2:end);
    results_table.Conductivity_S_m = m_Conductivity(2:end);
    results_table.Pt_loss_frac = m_Pt_loss(2:end);
    results_table.CO2_mol_m3 = m_COI(2:end);
    results_table.Thickness_um = m_Thickness(2:end) * 1e6;
    results_table.ECSA_m2_m3 = m_ECSA(2:end);
    
    % Calculate % change from baseline for Tornado chart
    pct_change_HF = ((m_HF(2:end) - m_HF(1)) / m_HF(1)) * 100;
    pct_change_LT = ((m_Lifetime(2:end) - m_Lifetime(1)) / m_Lifetime(1)) * 100;
    pct_change_SC = ((m_SC(2:end) - m_SC(1)) / m_SC(1)) * 100;
    pct_change_IEC = ((m_IEC(2:end) - m_IEC(1)) / m_IEC(1)) * 100;
    pct_change_Pt = ((m_Pt_loss(2:end) - m_Pt_loss(1)) / max(m_Pt_loss(1), 1e-10)) * 100;
    
    % Store in results
    results.PercentChange.HF = pct_change_HF;
    results.PercentChange.Lifetime = pct_change_LT;
    results.PercentChange.SC = pct_change_SC;
    results.PercentChange.IEC = pct_change_IEC;
    results.PercentChange.Pt = pct_change_Pt;
    results.Table = results_table;

    %% 4. Output Results Table
    fprintf('\n===========================================\n');
    fprintf('MECHANISM RANKING RESULTS\n');
    fprintf('===========================================\n');
    disp(results_table);
    fprintf('===========================================\n');
    
    %% 5. Create Tornado Chart
    create_tornado_chart(mechs(2:end), pct_change_LT, pct_change_HF, ...
                        pct_change_SC, pct_change_IEC, pct_change_Pt);
    
    %% 6. Create Sensitivity Heatmap
    create_sensitivity_heatmap(mechs(2:end), results_table);
    
    %% 7. Save Results
    save('mechanism_ranking_results.mat', 'results');
    writetable(results_table, 'mechanism_ranking_results.csv');
    disp('Results saved to mechanism_ranking_results.mat and .csv');
end

%% Local ODE Function with Toggle Logic
function dydt = ranked_ode_step(~, y, kinetics, op, p, switches)
    % RANKED_ODE_STEP
    % ODE function with mechanism toggles for ranking analysis
    
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
    
    % Get base rates from reaction_rates
    rates = reaction_rates(species, op.Temperature, kinetics, op, props, p);
    
    % Apply Mechanism Switches
    if ~switches.peroxide
        % Disable peroxide/radical external formation (R6-R13)
        rates.R6 = 0; rates.R7 = 0; rates.R8 = 0; rates.R9 = 0; 
        rates.R10 = 0; rates.R11 = 0; rates.R12 = 0; rates.R13 = 0;
        rates.S_H2O2 = 0;
    end
    
    if ~switches.fenton
        % Disable Fenton catalysis (R1 to R5)
        rates.R1 = 0; rates.R2 = 0; rates.R3 = 0; rates.R4 = 0; rates.R5 = 0;
    end
    
    if ~switches.sidechain
        % Disable side-chain attack
        rates.R14 = 0; rates.R15 = 0; rates.R16 = 0;
    end
    
    if ~switches.unzipping
        % Disable backbone unzipping
        rates.R17 = 0; rates.R18 = 0; rates.R19 = 0; rates.R20 = 0; 
        rates.R21 = 0; rates.R22 = 0; rates.R23 = 0;
    end
    
    if ~switches.pt
        % Disable Pt kinetics
        if isfield(rates, 'Pt')
            rates.Pt.dissolution = 0;
            rates.Pt.redeposition = 0;
            rates.Pt.net = 0;
        end
    end
    
    if ~switches.bidirectional || ~switches.conductivity
        % Disable conductivity coupling (use initial conductivity)
        props.sigma = p.sigma_0;
    end
    
    % Construct Derivatives
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

%% Helper Function: Create Tornado Chart
function create_tornado_chart(mech_names, pct_LT, pct_HF, pct_SC, pct_IEC, pct_Pt)
    % CREATE_TORNADO_CHART
    % Generates a multi-panel tornado chart for all metrics
    
    figure('Color', 'w', 'Name', 'Mechanism Ranking: Tornado Chart', ...
           'Position', [100, 100, 1200, 600]);
    
    num_mechs = length(mech_names);
    [~, sort_idx_lt] = sort(abs(pct_LT), 'ascend');
    
    % Subplot 1: Lifetime Impact
    subplot(2, 3, 1);
    barh(pct_LT(sort_idx_lt), 'FaceColor', [0.2 0.4 0.6]);
    set(gca, 'YTick', 1:num_mechs, 'YTickLabel', mech_names(sort_idx_lt), 'FontSize', 9);
    xlabel('% Change in Lifetime');
    title('Impact on Lifetime');
    grid on;
    xlim([-max(abs(pct_LT))*1.2, max(abs(pct_LT))*1.2]);
    
    % Subplot 2: HF Impact
    subplot(2, 3, 2);
    [~, sort_idx_hf] = sort(abs(pct_HF), 'ascend');
    barh(pct_HF(sort_idx_hf), 'FaceColor', [0.8 0.2 0.2]);
    set(gca, 'YTick', 1:num_mechs, 'YTickLabel', mech_names(sort_idx_hf), 'FontSize', 9);
    xlabel('% Change in HF');
    title('Impact on HF Production');
    grid on;
    
    % Subplot 3: SC Impact
    subplot(2, 3, 3);
    [~, sort_idx_sc] = sort(abs(pct_SC), 'ascend');
    barh(pct_SC(sort_idx_sc), 'FaceColor', [0.2 0.8 0.2]);
    set(gca, 'YTick', 1:num_mechs, 'YTickLabel', mech_names(sort_idx_sc), 'FontSize', 9);
    xlabel('% Change in SC');
    title('Impact on Side Chain');
    grid on;
    
    % Subplot 4: IEC Impact
    subplot(2, 3, 4);
    [~, sort_idx_iec] = sort(abs(pct_IEC), 'ascend');
    barh(pct_IEC(sort_idx_iec), 'FaceColor', [0.8 0.6 0.2]);
    set(gca, 'YTick', 1:num_mechs, 'YTickLabel', mech_names(sort_idx_iec), 'FontSize', 9);
    xlabel('% Change in IEC');
    title('Impact on IEC');
    grid on;
    
    % Subplot 5: Pt Impact
    subplot(2, 3, 5);
    [~, sort_idx_pt] = sort(abs(pct_Pt), 'ascend');
    barh(pct_Pt(sort_idx_pt), 'FaceColor', [0.6 0.2 0.8]);
    set(gca, 'YTick', 1:num_mechs, 'YTickLabel', mech_names(sort_idx_pt), 'FontSize', 9);
    xlabel('% Change in Pt Loss');
    title('Impact on Pt Degradation');
    grid on;
    
    sgtitle('Tornado Chart: Mechanism Impact Analysis on PEMFC Degradation', ...
            'FontWeight', 'bold', 'FontSize', 14);
end

%% Helper Function: Create Sensitivity Heatmap
function create_sensitivity_heatmap(mech_names, results_table)
    % CREATE_SENSITIVITY_HEATMAP
    % Generates a heatmap of mechanism sensitivities
    
    % Extract data for heatmap
    metrics = {'Lifetime_h', 'HF_mol_m3', 'SC_mol_m3', 'IEC_meq_g', 'Pt_loss_frac'};
    metric_labels = {'Lifetime (h)', 'HF (mol/m³)', 'SC (mol/m³)', 'IEC (meq/g)', 'Pt Loss'};
    
    data = zeros(length(mech_names), length(metrics));
    for i = 1:length(mech_names)
        for j = 1:length(metrics)
            data(i, j) = results_table.(metrics{j})(i);
        end
    end
    
    % Normalize each metric to [0,1] for heatmap
    data_norm = zeros(size(data));
    for j = 1:length(metrics)
        min_val = min(data(:, j));
        max_val = max(data(:, j));
        if max_val > min_val
            data_norm(:, j) = (data(:, j) - min_val) / (max_val - min_val);
        else
            data_norm(:, j) = 0.5;
        end
    end
    
    % Create heatmap
    figure('Color', 'w', 'Name', 'Mechanism Sensitivity Heatmap', ...
           'Position', [200, 200, 800, 400]);
    
    imagesc(data_norm);
    colormap('hot');
    colorbar('Label', 'Normalized Sensitivity');
    
    set(gca, 'XTick', 1:length(metrics), 'XTickLabel', metric_labels, ...
             'YTick', 1:length(mech_names), 'YTickLabel', mech_names);
    xlabel('Performance Metrics');
    ylabel('Mechanisms');
    title('Mechanism Sensitivity Heatmap', 'FontWeight', 'bold');
    
    % Add text annotations
    for i = 1:length(mech_names)
        for j = 1:length(metrics)
            text(j, i, sprintf('%.2f', data_norm(i, j)), ...
                 'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
        end
    end
end