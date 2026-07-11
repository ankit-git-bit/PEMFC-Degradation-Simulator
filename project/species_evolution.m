function species_evolution()
    % SPECIES_EVOLUTION.M
    % Runs the chemical degradation simulation and automatically plots 
    % the time evolution of the requested chemical species[cite: 5].

    %% 1. Initialization
    p        = parameters();
    kinetics = kinetic_parameters();
    op       = operating_conditions();
    switches = mechanism_switch();
    
    % Set a fixed H2O2 pool for consistent baseline behavior
    op.H2O2 = 0.001; 

    % Initial State Vector
    y0 = zeros(18, 1);
    y0(1) = p.Fe2_0;       y0(2) = p.Fe3_0;
    y0(3) = p.OH_0;        y0(4) = p.OOH_0;      y0(5) = p.H_0;
    y0(6) = p.SC_SO3H_0;   y0(7) = p.SC_O_0;     y0(8) = p.BB_O_0;
    y0(9) = p.HF_0;        y0(10) = p.CO2_0;
    y0(11) = p.CF2_7_0;    y0(12) = p.CF2_6_0;   y0(13) = p.CF2_5_0;
    y0(14) = p.CF2_4_0;    y0(15) = p.CF2_3_0;   y0(16) = p.CF2_2_0;
    y0(17) = p.CF2_1_0;    y0(18) = p.HOCF2CF2SO3H_0;

    % Run ODE Solver
    options = odeset('RelTol', p.rel_tol, 'AbsTol', p.abs_tol);
    t_span = [0, p.t_final];
    disp('Running simulation for species evolution plot...');
    [t_out, y_out] = ode45(@(t, y) evo_ode_step(t, y, kinetics, op, switches), t_span, y0, options);

    % Convert time to hours for plotting
    t_hr = t_out / 3600;

    %% 2. Plotting Species Evolution
    % Create a 2x2 subplot grid to handle different concentration scales
    figure('Color', 'w', 'Name', 'Species Evolution Dynamics', 'Position', [100, 100, 1000, 700]);
    
    % Subplot 1: Iron Species (Linear Scale)
    subplot(2, 2, 1);
    plot(t_hr, y_out(:, 1), 'b', 'LineWidth', 2); hold on;
    plot(t_hr, y_out(:, 2), 'r', 'LineWidth', 2);
    title('Iron Redox Cycle', 'FontWeight', 'bold');
    xlabel('Time (Hours)'); ylabel('Concentration (M)');
    legend('Fe^{2+}', 'Fe^{3+}', 'Location', 'best');
    grid on;

    % Subplot 2: Radicals (Logarithmic Scale due to tiny steady-state values)[cite: 5]
    subplot(2, 2, 2);
    semilogy(t_hr, y_out(:, 3), 'k', 'LineWidth', 2); hold on;
    semilogy(t_hr, y_out(:, 4), 'm', 'LineWidth', 2);
    semilogy(t_hr, y_out(:, 5), 'g', 'LineWidth', 2);
    title('Radical Species', 'FontWeight', 'bold');
    xlabel('Time (Hours)'); ylabel('Concentration (M) [Log Scale]');
    legend('OH^{\bullet}', 'OOH^{\bullet}', 'H^{\bullet}', 'Location', 'best');
    grid on;

    % Subplot 3: Ionomer / Membrane Species (Linear Scale)
    subplot(2, 2, 3);
    plot(t_hr, y_out(:, 6), 'LineWidth', 2, 'Color', [0 0.5 0]); hold on;
    plot(t_hr, y_out(:, 7), 'LineWidth', 2, 'Color', [0.8 0.4 0]);
    plot(t_hr, y_out(:, 8), 'LineWidth', 2, 'Color', [0.4 0.2 0.6]);
    title('Ionomer Degradation Intermediates', 'FontWeight', 'bold');
    xlabel('Time (Hours)'); ylabel('Concentration (M)');
    legend('SC-SO_3H', 'SC-O^{\bullet}', 'BB-O^{\bullet}', 'Location', 'best');
    grid on;

    % Subplot 4: Accumulating Degradation Products (Linear Scale)
    subplot(2, 2, 4);
    plot(t_hr, y_out(:, 9), 'c', 'LineWidth', 2); hold on;
    plot(t_hr, y_out(:, 10), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    title('Degradation Products', 'FontWeight', 'bold');
    xlabel('Time (Hours)'); ylabel('Concentration (M)');
    legend('HF', 'CO_2', 'Location', 'best');
    grid on;

    sgtitle('PEMFC Chemical Degradation: Species Evolution', 'FontSize', 16, 'FontWeight', 'bold');
    disp('Species evolution plot generated.');
end

%% Local ODE Step Function
function dydt = evo_ode_step(~, y, kinetics, op, switches)
    
    species.Fe2     = y(1);  species.Fe3   = y(2);
    species.OH      = y(3);  species.OOH   = y(4);   species.H     = y(5);
    species.SC_SO3H = y(6);  species.SC_O  = y(7);   species.BB_O  = y(8);
    species.HF      = y(9);  species.CO2   = y(10);
    species.CF2_7   = y(11); species.CF2_6 = y(12);  species.CF2_5 = y(13);
    species.CF2_4   = y(14); species.CF2_3 = y(15);  species.CF2_2 = y(16);
    species.CF2_1   = y(17); species.HOCF2CF2SO3H = y(18);
    
    species.H2O2 = op.H2O2; 

    % Get base rates
    rates = reaction_rates(species, op.Temperature, kinetics, op);
    
    % Apply Mechanism Switches
    if ~switches.peroxide
        rates.R6 = 0; rates.R7 = 0; rates.R8 = 0; rates.R9 = 0; 
        rates.R10 = 0; rates.R11 = 0; rates.R12 = 0; rates.R13 = 0;
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

    % Construct Derivatives
    dydt = zeros(18, 1);
    
    dydt(1)  = -rates.R1 + rates.R2 - rates.R3 - rates.R4 + rates.R5;
    dydt(2)  =  rates.R1 - rates.R2 + rates.R3 + rates.R4 - rates.R5;
    dydt(3)  =  rates.R1 + 2*rates.R6 - rates.R3 - rates.R7 + rates.R8 ...
               - 2*rates.R10 - rates.R11 - rates.R12 - rates.R14 - 3*rates.R15 ...
               - rates.R16 - 2*(rates.R17 + rates.R18 + rates.R19 + rates.R20 ...
               + rates.R21 + rates.R22 + rates.R23);
    dydt(4)  =  rates.R2 - rates.R4 - rates.R5 + rates.R7 - rates.R8 ...
               - 2*rates.R9 - rates.R11 + rates.R13;
    dydt(5)  =  rates.R12 - rates.R13;
    
    dydt(6)  = -rates.R14;
    dydt(7)  =  rates.R14 - rates.R15;
    dydt(8)  =  rates.R15 - rates.R16;
    dydt(9)  =  6*rates.R15 + 3*rates.R16 + 2*(rates.R17 + rates.R18 ...
               + rates.R19 + rates.R20 + rates.R21 + rates.R22 + rates.R23);
    dydt(10) =  3*rates.R15 + rates.R17 + rates.R18 + rates.R19 ...
               + rates.R20 + rates.R21 + rates.R22 + rates.R23;
           
    dydt(11) =  2*rates.R16 - rates.R17;
    dydt(12) =  rates.R17 - rates.R18;
    dydt(13) =  rates.R18 - rates.R19;
    dydt(14) =  rates.R19 - rates.R20;
    dydt(15) =  rates.R20 - rates.R21;
    dydt(16) =  rates.R21 - rates.R22;
    dydt(17) =  rates.R22 - rates.R23;
    
    dydt(18) =  rates.R14;
end