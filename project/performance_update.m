function perf = performance_update(props, species, p, op)
    % PERFORMANCE_UPDATE.M
    % Converts physical degradation state into macroscopic PEMFC
    % performance metrics (voltage, power, efficiency) via a simplified
    % lumped Butler-Volmer + Ohmic loss model.
    %
    % FIX (previous bug): the original version read R, F, alpha, ECSA_0,
    % i0_ref, V_ocv from `op`, but operating_conditions.m never defines
    % those fields (they live in parameters.m / kinetic_parameters.m).
    % This version sources them consistently from `p`.
    %
    % Inputs:
    %   props   - membrane property struct (from degradation_update.m),
    %             uses L_mem, sigma, ECSA
    %   species - current chemical-state struct (unused directly here,
    %             kept for interface symmetry / future mass-transfer terms)
    %   p       - parameters.m struct (R, F, alpha, i0_ref, V_ocv, ECSA_0)
    %   op      - operating_conditions.m struct (Temperature, current_density)
    %
    % Outputs (struct perf):
    %   voltage         - cell voltage (V)
    %   power           - power density (W/m^2)
    %   efficiency      - voltage efficiency relative to H2 HHV (1.48 V)
    %   current_density - pass-through operating current density (A/m^2)
    %   V_ohmic, V_act  - loss breakdown (V), for diagnostics

    i_op = op.current_density; % A/m^2

    %% 1. Ohmic Loss
    % V_ohmic = i * (L_mem / sigma)   [Ohm's law for ionic resistance]
    R_ohmic = props.L_mem / max(props.sigma, 1e-6); % Ohm*m^2
    V_ohmic = i_op * R_ohmic;

    %% 2. Activation (Butler-Volmer / Tafel) Loss
    % FIX (previous bug): p.i0_ref (kinetics.perf.i0_ref) is a per-REAL-
    % catalyst-surface-area exchange current density (~1e-4 A/m^2_Pt, a
    % representative literature order of magnitude, see e.g. Neyerlin, Gu,
    % Jorne & Gasteiger, J. Electrochem. Soc. 153 (2006) A1955). The
    % previous version compared it directly against the GEOMETRIC operating
    % current density (op.current_density, ~1000 A/m^2_geometric), which is
    % a roughness-factor mismatch of ~2 orders of magnitude and produced an
    % activation loss exceeding the open-circuit voltage (i.e. voltage
    % floored at 0 for all conditions). The fix converts the per-real-area
    % exchange current density to a per-geometric-area value using the
    % catalyst-layer roughness factor (real Pt area per unit geometric
    % area = ECSA [m^2/g_Pt] x Pt areal loading [g_Pt/m^2_geometric]),
    % which is standard PEMFC electrode-kinetics practice.
    roughness_factor = props.ECSA * (p.Pt_loading_mg_cm2 * 1e-3 * 1e4); % m^2_Pt / m^2_geometric
    i0_eff = p.i0_ref * roughness_factor; % A/m^2_geometric
    V_act = (p.R * op.Temperature / (p.alpha * p.F)) * log(max(i_op, 1e-8) / max(i0_eff, 1e-12));
    V_act = max(V_act, 0); % activation loss cannot be negative in this simplified model

    %% 3. Cell Voltage, Power, Efficiency
    perf.V_ohmic = V_ohmic;
    perf.V_act   = V_act;
    perf.voltage = max(p.V_ocv - V_act - V_ohmic, 0);
    perf.power   = perf.voltage * i_op; % W/m^2
    perf.efficiency = perf.voltage / 1.48; % relative to H2 HHV (1.48 V thermoneutral)
    perf.current_density = i_op;
end
