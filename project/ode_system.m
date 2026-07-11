function dydt = ode_system(t, y, kinetics, op, p, props)
    % ODE_SYSTEM.M
    % Defines the system of Ordinary Differential Equations (ODEs) for the
    % FAST chemical degradation state (Fenton/radical chemistry, membrane
    % attack, backbone unzipping, and Pt dissolution/redeposition).
    %
    % This file ONLY maps pre-calculated reaction rates (reaction_rates.m)
    % to state derivatives. It NEVER calculates Arrhenius rate constants
    % directly, and it NEVER updates the macroscopic physical properties
    % (IEC, membrane thickness, conductivity, hydrogen crossover, ECSA) -
    % those are SLOW variables updated once per macro-timestep by
    % degradation_update.m (see lifetime_simulation.m), which is the sole
    % authorized location for that logic per the project's module contract.
    % This two-timescale split (stiff sub-second chemistry integrated with
    % ode15s, algebraic hourly property update) is a standard strategy for
    % multi-timescale degradation models and avoids solving a physically
    % stiff, numerically awkward 24-state ODE with rate constants spanning
    % >10 orders of magnitude in timescale.
    %
    % State vector y layout (21 states):
    %  1:Fe2   2:Fe3   3:OH   4:OOH   5:H   6:SC_SO3H   7:SC_O   8:BB_O
    %  9:HF   10:CO2  11:CF2_7 12:CF2_6 13:CF2_5 14:CF2_4 15:CF2_3
    % 16:CF2_2 17:CF2_1 18:HOCF2CF2SO3H 19:Pt 20:Pt2 21:H2O2
    %
    % props (input) supplies the CURRENT (held-fixed-over-this-macrostep)
    % membrane state: props.L_mem, props.sigma, props.ECSA, props.IEC,
    % props.crossover_H2 - all owned/updated by degradation_update.m between
    % calls to this integrator.

    %% 1. Unpack State Vector (y) into Species Struct
    % Ensure no species concentration goes below a physical floor
    y = max(y, 0);

    species.Fe2          = y(1);
    species.Fe3          = y(2);
    species.OH           = y(3);
    species.OOH          = y(4);
    species.H            = y(5);
    species.SC_SO3H      = y(6);
    species.SC_O         = y(7);
    species.BB_O         = y(8);
    species.HF           = y(9);
    species.CO2          = y(10);
    species.CF2_7        = y(11);
    species.CF2_6        = y(12);
    species.CF2_5        = y(13);
    species.CF2_4        = y(14);
    species.CF2_3        = y(15);
    species.CF2_2        = y(16);
    species.CF2_1        = y(17);
    species.HOCF2CF2SO3H = y(18);
    species.Pt            = y(19);  % Platinum metal
    species.Pt2           = y(20);  % Platinum ions
    species.H2O2          = y(21);  % Hydrogen peroxide

    % ECSA is a slow state carried in props (owned by degradation_update.m),
    % but is needed here as a READ-ONLY input to reaction_rates.m (Pt
    % kinetics + peroxide generation depend on current catalyst area).
    species.ECSA = props.ECSA;

    %% 2. Fetch Reaction Rates
    % All Arrhenius and kinetic calculations are strictly isolated in
    % reaction_rates.m (full/coupled mode: crossover-driven H2O2 source +
    % Pt kinetics via pt_dissolution.m).
    rates = reaction_rates(species, op.Temperature, kinetics, op, props, p);

    %% 3. Construct ODE System (dydt) using ONLY Reaction Rates
    dydt = zeros(length(y), 1);

    % Iron Redox Cycle (R1 - R5)
    dydt(1)  = -rates.R1 + rates.R2 - rates.R3 - rates.R4 + rates.R5;              % d[Fe2+]/dt
    dydt(2)  =  rates.R1 - rates.R2 + rates.R3 + rates.R4 - rates.R5;              % d[Fe3+]/dt

    % Radicals (OH, OOH, H): NOT integrated as differential states. They
    % are solved algebraically via quasi-steady-state approximation inside
    % reaction_rates.m (see the QSSA rationale note there) because their
    % intrinsic reaction timescale (~ns-us) is many orders of magnitude
    % faster than the hour-scale macro-timestep, which would otherwise make
    % this system numerically intractable for a stiff BDF solver. Their
    % state-vector slots (y(3),y(4),y(5)) are held constant by the
    % integrator (dydt=0) and overwritten with the fresh QSSA solution by
    % lifetime_simulation.m after each macro-step, purely for history/
    % reporting purposes; reaction_rates.m already substitutes the correct
    % QSSA values into R3, R7, R10-R23 etc. regardless of what y(3:5)
    % currently holds.
    dydt(3)  = 0;
    dydt(4)  = 0;
    dydt(5)  = 0;

    % Ionomer Side-Chain Degradation (R14 - R16)
    dydt(6)  = -rates.R14;                                                         % d[SC-SO3H]/dt
    dydt(7)  =  rates.R14 - rates.R15;                                             % d[SC-O]/dt
    dydt(8)  =  rates.R15 - rates.R16;                                             % d[BB-O]/dt

    % Cumulative Degradation Products
    dydt(9)  =  6*rates.R15 + 3*rates.R16 + 2*(rates.R17 + rates.R18 ...
               + rates.R19 + rates.R20 + rates.R21 + rates.R22 + rates.R23);       % d[HF]/dt
    dydt(10) =  3*rates.R15 + rates.R17 + rates.R18 + rates.R19 ...
               + rates.R20 + rates.R21 + rates.R22 + rates.R23;                    % d[CO2]/dt
    dydt(18) =  rates.R14;                                                         % d[HOCF2CF2SO3H]/dt

    % Backbone Unzipping Cascade (R17 - R23)
    dydt(11) =  2*rates.R16 - rates.R17;                                           % d[CF2_7]/dt
    dydt(12) =  rates.R17 - rates.R18;                                             % d[CF2_6]/dt
    dydt(13) =  rates.R18 - rates.R19;                                             % d[CF2_5]/dt
    dydt(14) =  rates.R19 - rates.R20;                                             % d[CF2_4]/dt
    dydt(15) =  rates.R20 - rates.R21;                                             % d[CF2_3]/dt
    dydt(16) =  rates.R21 - rates.R22;                                             % d[CF2_2]/dt
    dydt(17) =  rates.R22 - rates.R23;                                             % d[CF2_1]/dt

    %% 4. Peroxide Balance
    % H2O2 generation from crossover minus consumption in Fenton reactions
    dydt(21) = rates.S_H2O2 - rates.R1 - rates.R2 - rates.R4 ...
               - rates.R6 - rates.R7;                                              % d[H2O2]/dt

    %% 5. Pt Kinetics (Dissolution/Redeposition), via pt_dissolution.m
    if isfield(rates, 'Pt')
        dydt(19) = rates.Pt.Pt;       % d[Pt]/dt
        dydt(20) = rates.Pt.Pt2;      % d[Pt2+]/dt
    end

    %% 6. Ensure Physical Bounds
    % Soft floor guard: prevents the stiff solver from overshooting a
    % depleting species below zero between accepted steps.
    dydt = max(dydt, -y/10);
end
