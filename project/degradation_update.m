function props = degradation_update(species, p, kinetics, op, current_props, dt)
    % DEGRADATION_UPDATE.M
    % Centralizes ALL macroscopic physical-property evolution logic. Per the
    % project's module contract, this is the ONLY function that updates:
    %   IEC, membrane thickness (L_mem), conductivity (sigma),
    %   hydrogen crossover (crossover_H2), and ECSA.
    % No other module writes to these fields; ode_system.m only reads them
    % (as fixed "slow" inputs) and pt_dissolution.m/hydrogen_crossover.m
    % only compute local rates/fluxes, not the state itself.
    %
    % Update path (per project architecture diagram):
    %   HF / side-chain & backbone mass loss -> IEC -> Conductivity
    %   -> Thickness -> Hydrogen Crossover ; Pt loss -> ECSA
    %
    % Inputs:
    %   species        - struct with current chemical-state concentrations
    %                    (SC_SO3H, CF2_7..CF2_1, HF, Pt, ...)
    %   p              - parameters.m struct
    %   kinetics       - kinetic_parameters.m struct (uses kinetics.coupling)
    %   op             - operating_conditions.m struct (unused directly here,
    %                    kept for interface consistency / future RH-coupling)
    %   current_props  - props struct from the previous macro-timestep
    %                    (must contain ECSA; other fields optional)
    %   dt             - macro-timestep length in seconds (default: 3600 s,
    %                    i.e. the 1-hour step used by lifetime_simulation.m)
    %
    % Output:
    %   props - updated struct with fields IEC, sigma, L_mem, crossover_H2,
    %           ECSA (plus pass-through diagnostic fields)

    if nargin < 6 || isempty(dt)
        dt = 3600; % seconds (default macro-timestep)
    end

    props = current_props; % carry forward any fields not touched below

    %% 1. Ion Exchange Capacity (IEC)
    % IEC is directly proportional to surviving sulfonic-acid side chains
    % (Reference: qualitative IEC-vs-sulfonate-loss relationship,
    % Kreuer, Chem. Rev. 104 (2004) 4637; exact proportionality is a
    % simplifying ASSUMPTION appropriate for a lumped 0D model).
    sc_fraction = max(species.SC_SO3H / p.SC_SO3H_0, 0);
    props.IEC = p.IEC_0 * sc_fraction;

    %% 2. Conductivity (sigma) - percolation-theory power law vs. IEC
    tau = kinetics.coupling.tau; % see documented ASSUMPTION in kinetic_parameters.m
    props.sigma = max(p.sigma_0 * (sc_fraction)^(1/tau), 1e-6); % S/m, floored

    %% 3. Membrane Thickness (L_mem)
    % FIX (previous bug): the surviving backbone mass must include BOTH the
    % dynamically tracked CF2_7..CF2_1 unzipping-cascade pools AND the
    % fixed (non-dynamically-tracked) CF2_backbone_0 pool (positions 8-13,
    % which are not consumed until the cascade reaches position 7), not a
    % nonexistent aggregate field. See parameters.m Section 4.
    sum_cf2_tracked = species.CF2_7 + species.CF2_6 + species.CF2_5 + ...
                       species.CF2_4 + species.CF2_3 + species.CF2_2 + species.CF2_1;
    structural_mass       = sum_cf2_tracked + p.CF2_backbone_0 + species.SC_SO3H;
    total_initial_mass    = p.CF2_total_0 + p.SC_SO3H_0;
    mass_retention         = max(structural_mass / total_initial_mass, 0);
    props.L_mem = p.delta_membrane_0 * mass_retention;

    %% 4. Hydrogen Crossover
    % Crossover increases as the membrane thins and HF-induced pinholes
    % accumulate (damage multiplier keyed to HF_critical threshold - see
    % parameters.m Section 14; functional form is an engineering
    % simplification, not a fitted literature correlation).
    if props.L_mem > 0
        thickness_ratio = props.L_mem / p.delta_membrane_0; % 1 = BOL, ->0 = fully thinned
        hf_damage_multiplier = 1 + (species.HF / p.HF_critical);
        props.crossover_H2 = (p.permeability_H2 / max(thickness_ratio, 1e-3)) * hf_damage_multiplier;
    else
        props.crossover_H2 = p.permeability_H2 * 10; % Catastrophic-failure floor state
    end

    %% 5. ECSA (electrochemically active surface area)
    % Explicit first-order decay driven by Pt dissolution, integrated over
    % this macro-timestep. Consistent with the rate expression used by
    % pt_dissolution.m's d_pt.ECSA term (single physical model, evaluated
    % here as the authorized state-update location).
    if isfield(species, 'Pt') && isfield(p, 'Pt_0') && p.Pt_0 > 0 && isfield(current_props, 'ECSA')
        pt_loss_fraction = max(1 - species.Pt / p.Pt_0, 0);
        ecsa_decay_rate = -kinetics.coupling.k_ECSA * current_props.ECSA * (1 + pt_loss_fraction);
        props.ECSA = max(current_props.ECSA + ecsa_decay_rate * dt, 0);
    elseif isfield(current_props, 'ECSA')
        props.ECSA = current_props.ECSA; % No Pt tracked in this call context - hold constant
    else
        props.ECSA = p.ECSA_0;
    end
end
