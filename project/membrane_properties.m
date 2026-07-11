function props = membrane_properties(species, p, op, authoritative_props)
    % MEMBRANE_PROPERTIES.M
    % Diagnostic/derived-quantity module: water uptake, porosity, and
    % mechanical/thermal property estimates as functions of the CURRENT
    % degradation state.
    %
    % IMPORTANT (module contract): per the project's architecture,
    % degradation_update.m is the SOLE authority for IEC, sigma, L_mem,
    % crossover_H2 and ECSA. This function does NOT recompute those values
    % independently (an earlier version of this file did, which produced
    % inconsistent duplicate physics against degradation_update.m). Instead
    % it takes them as read-only inputs (authoritative_props, i.e. the
    % struct already produced by degradation_update.m for the current
    % timestep) and only adds NEW derived diagnostic quantities.
    %
    % Inputs:
    %   species             - current chemical-state struct
    %   p                   - parameters.m struct
    %   op                  - operating_conditions.m struct
    %   authoritative_props - props struct from degradation_update.m
    %                         (must contain IEC, sigma, L_mem, crossover_H2)
    %
    % Outputs (props): pass-through of authoritative_props fields, plus:
    %   lambda, water_volume_fraction, water_uptake_pct, porosity,
    %   free_volume, E_modulus, yield_strength, elongation_break, Tg,
    %   k_thermal, degradation_index, remaining_life, HF_fraction

    props = authoritative_props; % pass-through: IEC, sigma, L_mem, crossover_H2, ECSA

    mass_fraction = max(min(props.L_mem / p.delta_membrane_0, 1), 0.01);
    acid_fraction = max(species.SC_SO3H / p.SC_SO3H_0, 0);

    %% 1. Water Uptake
    % Water-content correlation follows the qualitative form of Springer,
    % Zawodzinski & Gottesfeld (J. Electrochem. Soc. 138 (1991) 2334) for
    % RH- and IEC-dependence; the specific 5*RH sorption-isotherm shape and
    % Arrhenius temperature factor are simplifying ASSUMPTIONS for this
    % lumped 0D model, not a direct literature fit.
    lambda_ref = 14; % water molecules per -SO3H at full hydration (Springer et al. 1991)
    RH_factor = 1 - exp(-5 * op.RH_cathode);
    T_factor  = exp(-3000 * (1/op.Temperature - 1/353.15));
    props.lambda = lambda_ref * acid_fraction * RH_factor * T_factor;

    V_w = 18e-6;    % m^3/mol, molar volume of water
    V_mem = 0.0005; % m^3/mol, ASSUMPTION: approximate ionomer molar volume per -SO3H site
    props.water_volume_fraction = (props.lambda * V_w) / (V_mem + props.lambda * V_w);
    props.water_uptake_pct = (props.lambda * 0.018) * 100; % g water / g dry ionomer equivalent, x100

    %% 2. Porosity and Free Volume
    % TODO/ASSUMPTION: porosity growth with mass loss and water content is
    % a simplified engineering estimate, not fit to a specific porosimetry
    % dataset.
    porosity_0 = 0.1;
    degradation_porosity = (1 - mass_fraction) * 0.5;
    props.porosity = min(porosity_0 + degradation_porosity + props.water_volume_fraction, 0.8);
    props.free_volume = props.porosity * (1 - props.water_volume_fraction);

    %% 3. Mechanical Properties (order-of-magnitude, TODO/ASSUMPTION)
    % Baseline values representative of Nafion-type PFSA membranes
    % (typical literature ranges: E ~ 100-300 MPa, yield ~ 10-25 MPa,
    % elongation ~ 200-300%; see e.g. Kusoglu & Weber, Chem. Rev. 117
    % (2017) 987 for a review of Nafion mechanical properties). Degradation
    % scaling with mass_fraction/porosity below is a simplifying
    % ASSUMPTION, not fit to a specific mechanical-test dataset.
    E0 = 250e6;         % Pa
    sigma_y0 = 15e6;     % Pa
    elon_0 = 2.5;        % 250%
    props.E_modulus       = E0 * mass_fraction * (1 - 0.5 * props.porosity);
    props.yield_strength  = sigma_y0 * mass_fraction;
    props.elongation_break = elon_0 * mass_fraction;

    %% 4. Thermal Properties (order-of-magnitude, TODO/ASSUMPTION)
    Tg_0 = 383.15; % K, ~110 C, typical Nafion dry Tg
    k0 = 0.25;     % W/(m K), typical Nafion thermal conductivity
    props.Tg = Tg_0 - 50 * props.water_volume_fraction; % plasticization effect
    props.k_thermal = k0 * (1 - 2 * props.porosity / 3);

    %% 5. Diagnostic Metrics
    props.degradation_index = 1 - mass_fraction;
    props.remaining_life = mass_fraction;
    if isfield(species, 'HF')
        props.HF_fraction = min(species.HF / p.HF_critical, 1);
    end
end
