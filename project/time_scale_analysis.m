function time_scale_analysis()
    % TIME_SCALE_ANALYSIS.M
    % Calculates the relative contribution of each reaction at every timestep
    % and plots the reaction contribution dynamics over time[cite: 5].

    %% 1. Initialization and ODE Solve
    p        = parameters();
    kinetics = kinetic_parameters();
    op       = operating_conditions();
    switches = mechanism_switch();
    
    % Use a fixed H2O2 pool to match previous operating conditions
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
    disp('Running simulation for time-scale analysis...');
    [t_out, y_out] = ode45(@(t, y) ts_ode_step(t, y, kinetics, op, switches), t_span, y0, options);

    %% 2. Calculate Reaction Contributions
    num_steps = length(t_out);
    num_reactions = 23;
    contributions = zeros(num_steps, num_reactions);

    disp('Calculating relative contributions at each timestep...');
    for i = 1:num_steps
        % Map current state to species struct
        species.Fe2     = y_out(i, 1);  species.Fe3   = y_out(i, 2);
        species.OH      = y_out(i, 3);  species.OOH   = y_out(i, 4);   
        species.H       = y_out(i, 5);  species.SC_SO3H = y_out(i, 6);  
        species.SC_O    = y_out(i, 7);  species.BB_O  = y_out(i, 8);
        species.HF      = y_out(i, 9);  species.CO2   = y_out(i, 10);
        species.CF2_7   = y_out(i, 11); species.CF2_6 = y_out(i, 12);  
        species.CF2_5   = y_out(i, 13); species.CF2_4 = y_out(i, 14); 
        species.CF2_3   = y_out(i, 15); species.CF2_2 = y_out(i, 16);
        species.CF2_1   = y_out(i, 17); species.HOCF2CF2SO3H = y_out(i, 18);
        species.H2O2    = op.H2O2;

        % Calculate raw reaction rates
        rates = reaction_rates(species, op.Temperature, kinetics, op);
        
        % Extract into array
        r_array = [rates.R1, rates.R2, rates.R3, rates.R4, rates.R5, rates.R6, ...
                   rates.R7, rates.R8, rates.R9, rates.R10, rates.R11, rates.R12, ...
                   rates.R13, rates.R14, rates.R15, rates.R16, rates.R17, rates.R18, ...
                   rates.R19, rates.R20, rates.R21, rates.R22, rates.R23];
        
        % Calculate fractional contribution
        total_rate = sum(r_array);
        if total_rate > 0
            contributions(i, :) = r_array / total_rate;
        end
    end

    %% 3. Plot Reaction Contribution Dynamics
    figure('Color', 'w', 'Name', 'Time-Scale Analysis: Reaction Contributions', 'Position', [100, 100, 1000, 500]);
    
    % Use an area chart to visualize the 100% stacked contribution over time
    area(t_out / 3600, contributions * 100, 'EdgeColor', 'none');
    
    xlabel('Time (Hours)', 'FontWeight', 'bold');
    ylabel('Relative Contribution (%)', 'FontWeight', 'bold');
    title('Dynamic Time-Scale Analysis of Chemical Degradation Mechanisms', 'FontWeight', 'bold', 'FontSize', 14);
    
    % Define generic legend labels for the 23 reactions
    labels = cell(1, num_reactions);
    for k = 1:num_reactions
        labels{k} = sprintf('R%d', k);
    end
    
    % Format Plot
    ylim([0 100]);
    colormap(parula(num_reactions));
    lgd = legend(labels, 'Location', 'eastoutside', 'NumColumns', 2);
    title(lgd, 'Reactions');
    grid on;
    
    disp('Time-scale analysis complete. Area chart generated.');
end

%% Local ODE Step Function
function dydt = ts_ode_step(~, y, kinetics, op, switches)
    
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