function history = lifetime_simulation(t_total_hours)
    % LIFETIME_SIMULATION.M
    % Master closed-loop PEMFC membrane degradation simulation. Implements
    % the full feedback path described in the project architecture:
    %
    %   Operating Conditions -> Hydrogen crossover -> Peroxide generation
    %   -> Fenton reactions -> Pt dissolution -> Membrane attack ->
    %   Polymer unzipping -> HF release -> IEC loss -> Conductivity loss
    %   -> Hydrogen crossover increase -> Peroxide increase (closes loop)
    %
    % Two-timescale integration strategy:
    %   (1) FAST chemistry (Fenton/radical kinetics, membrane attack,
    %       backbone unzipping, Pt dissolution/redeposition) is integrated
    %       with a stiff ODE solver (ode15s) over each 1-hour macro-step,
    %       via ode_system.m (21-state vector).
    %   (2) SLOW physical properties (IEC, conductivity, thickness,
    %       hydrogen crossover, ECSA) are updated algebraically ONCE per
    %       macro-step by degradation_update.m (the sole authorized
    %       location for that logic).
    %   (3) The updated crossover state is fed back into the operating
    %       conditions for the NEXT macro-step by bidirectional_coupling.m,
    %       closing the feedback loop.
    %   (4) performance_update.m converts the current physical state into
    %       cell voltage/power for that hour.
    %
    % Input (optional):
    %   t_total_hours - simulation duration in hours (default: 5000)
    %
    % Output (struct history):
    %   time_hr, species (Nx21), L_mem, sigma, IEC, crossover_H2, ECSA,
    %   voltage, power, efficiency, HF, Pt, Pt2

    if nargin < 1 || isempty(t_total_hours)
        t_total_hours = 5000;
    end

    %% 1. Initialization
    p        = parameters();
    kinetics = kinetic_parameters();
    op       = operating_conditions();
    switches = mechanism_switch();

    dt_seconds = 3600; % 1-hour macro-timestep

    % Initial 21-state chemistry vector
    y0 = zeros(21, 1);
    y0(1)  = p.Fe2_0;             y0(2)  = p.Fe3_0;
    y0(3)  = p.OH_0;              y0(4)  = p.OOH_0;      y0(5)  = p.H_0;
    y0(6)  = p.SC_SO3H_0;         y0(7)  = p.SC_O_0;     y0(8)  = p.BB_O_0;
    y0(9)  = p.HF_0;              y0(10) = p.CO2_0;
    y0(11) = p.CF2_7_0;           y0(12) = p.CF2_6_0;    y0(13) = p.CF2_5_0;
    y0(14) = p.CF2_4_0;           y0(15) = p.CF2_3_0;    y0(16) = p.CF2_2_0;
    y0(17) = p.CF2_1_0;           y0(18) = p.HOCF2CF2SO3H_0;
    y0(19) = p.Pt_0;              y0(20) = p.Pt2_0;      y0(21) = p.H2O2_0;

    % Initial physical-property state (owned by degradation_update.m)
    props.L_mem       = p.delta_membrane_0;
    props.sigma       = p.sigma_0;
    props.IEC         = p.IEC_0;
    props.ECSA        = p.ECSA_0;
    props.crossover_H2 = p.permeability_H2;

    %% 2. History Preallocation
    history.time_hr      = zeros(t_total_hours, 1);
    history.species      = zeros(t_total_hours, 21);
    history.L_mem        = zeros(t_total_hours, 1);
    history.sigma        = zeros(t_total_hours, 1);
    history.IEC          = zeros(t_total_hours, 1);
    history.crossover_H2 = zeros(t_total_hours, 1);
    history.ECSA         = zeros(t_total_hours, 1);
    history.voltage      = zeros(t_total_hours, 1);
    history.power        = zeros(t_total_hours, 1);
    history.efficiency   = zeros(t_total_hours, 1);
    history.HF           = zeros(t_total_hours, 1);
    history.Pt           = zeros(t_total_hours, 1);
    history.Pt2          = zeros(t_total_hours, 1);

    y_current = y0;
    op_current = op; % mutated each hour by bidirectional_coupling.m
    options_stiff = odeset('RelTol', p.rel_tol, 'AbsTol', p.abs_tol);
    options_nonstiff = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);

    %% 3. Lifetime Simulation Loop
    fprintf('Starting lifetime simulation for %d hours...\n', t_total_hours);

    for t_hr = 1:t_total_hours
        %% 3a. Fast chemistry: integrate 21-state ODE over this hour
        % NOTE ON SOLVER CHOICE: the QSSA reduction for OH/OOH/H (see
        % reaction_rates.m) removes the worst of the radical-chemistry
        % stiffness, but the system can still be locally stiff during
        % rapid H2O2/Fe2+/Fe3+ transients. ode15s (a proper variable-order
        % BDF method) is tried first as it is generally the fastest choice
        % for this class of problem in MATLAB. If it fails to converge
        % (observed intermittently with GNU Octave's bundled SUNDIALS/IDA
        % backend during development/testing of this project), the solver
        % automatically falls back to ode45 (explicit adaptive
        % Dormand-Prince), which was empirically verified to integrate
        % this system reliably, just at a higher step-count/runtime cost.
        t_span = [0, dt_seconds];
        ode_func = @(t, y) chemistry_step(t, y, kinetics, op_current, p, props, switches);
        try
            [~, y_out] = ode15s(ode_func, t_span, y_current, options_stiff);
        catch
            [~, y_out] = ode45(ode_func, t_span, y_current, options_nonstiff);
        end
        y_current = max(y_out(end, :)', 0);

        % Unpack for property update / performance
        current_species = unpack_species(y_current, props.ECSA);

        %% 3b. Slow physical properties: sole update via degradation_update.m
        if switches.bidirectional
            props = degradation_update(current_species, p, kinetics, op_current, props, dt_seconds);
        end

        %% 3c. Feedback: crossover state -> next hour's operating conditions
        if isfield(switches, 'H2_crossover') && switches.H2_crossover
            op_current = bidirectional_coupling(props, op_current, p);
        end

        %% 3d. Performance conversion
        current_perf = performance_update(props, current_species, p, op_current);

        %% 3e. Store history
        history.time_hr(t_hr)      = t_hr;
        history.species(t_hr, :)   = y_current';
        history.L_mem(t_hr)        = props.L_mem;
        history.sigma(t_hr)        = props.sigma;
        history.IEC(t_hr)          = props.IEC;
        history.crossover_H2(t_hr) = props.crossover_H2;
        history.ECSA(t_hr)         = props.ECSA;
        history.voltage(t_hr)      = current_perf.voltage;
        history.power(t_hr)        = current_perf.power;
        history.efficiency(t_hr)   = current_perf.efficiency;
        history.HF(t_hr)           = y_current(9);
        history.Pt(t_hr)           = y_current(19);
        history.Pt2(t_hr)          = y_current(20);

        if mod(t_hr, 50) == 0
            fprintf('Completed %d / %d hours... (V=%.4f, sigma=%.3f S/m, ECSA=%.2f m2/g)\n', ...
                t_hr, t_total_hours, current_perf.voltage, props.sigma, props.ECSA);
        end
    end
    fprintf('Simulation complete.\n');
end

%% ------------------------------------------------------------------
function species = unpack_species(y, ECSA)
    species.Fe2 = y(1);   species.Fe3 = y(2);   species.OH  = y(3);
    species.OOH = y(4);   species.H   = y(5);   species.SC_SO3H = y(6);
    species.SC_O = y(7);  species.BB_O = y(8);  species.HF  = y(9);
    species.CO2 = y(10);  species.CF2_7 = y(11); species.CF2_6 = y(12);
    species.CF2_5 = y(13); species.CF2_4 = y(14); species.CF2_3 = y(15);
    species.CF2_2 = y(16); species.CF2_1 = y(17); species.HOCF2CF2SO3H = y(18);
    species.Pt = y(19);   species.Pt2 = y(20);  species.H2O2 = y(21);
    species.ECSA = ECSA;
end

%% ------------------------------------------------------------------
function dydt = chemistry_step(t, y, kinetics, op, p, props, switches)
    % Thin wrapper around ode_system.m that applies mechanism_switch.m
    % on/off flags by zeroing the relevant rates when any mechanism is
    % disabled. Delegates all rate-constant/Arrhenius/ODE-mapping logic to
    % reaction_rates.m and ode_system.m (no duplicated physics) when the
    % default all-mechanisms-on configuration is used.

    all_on = switches.peroxide && switches.fenton && switches.sidechain && ...
             switches.unzipping && switches.pt;

    if all_on
        dydt = ode_system(t, y, kinetics, op, p, props);
        return;
    end

    % Mechanism-isolation path: recompute rates with switches applied.
    species = unpack_species(y, props.ECSA);
    rates = reaction_rates(species, op.Temperature, kinetics, op, props, p);

    if ~switches.fenton
        rates.R1 = 0; rates.R2 = 0; rates.R3 = 0; rates.R4 = 0; rates.R5 = 0;
    end
    if ~switches.peroxide
        rates.S_H2O2 = 0;
    end
    if ~switches.sidechain
        rates.R14 = 0; rates.R15 = 0; rates.R16 = 0;
    end
    if ~switches.unzipping
        rates.R17 = 0; rates.R18 = 0; rates.R19 = 0; rates.R20 = 0;
        rates.R21 = 0; rates.R22 = 0; rates.R23 = 0;
    end
    if ~switches.pt && isfield(rates, 'Pt')
        rates.Pt.Pt = 0; rates.Pt.Pt2 = 0;
    end

    dydt = zeros(21, 1);
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
    dydt(21) =  rates.S_H2O2 - rates.R1 - rates.R2 - rates.R4 - rates.R6 - rates.R7;
    if isfield(rates, 'Pt')
        dydt(19) = rates.Pt.Pt;
        dydt(20) = rates.Pt.Pt2;
    end
    dydt = max(dydt, -y/10);
end
