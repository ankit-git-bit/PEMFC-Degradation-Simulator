function J = hydrogen_crossover(species, props, op, p)
    % HYDROGEN_CROSSOVER.M
    % Calculates dynamic gas crossover flux (J_H2, J_O2) based on
    % membrane conductivity, thickness, and operating conditions.
    
    %% 1. Calculate Permeability (P)
    % Permeability is a function of Temperature (T) and membrane state.
    % Uses an Arrhenius relationship for temperature dependence.
    P_H2_ref = p.permeability_H2; 
    P_H2 = P_H2_ref * exp(-(p.Ea_perm / p.R) * (1/op.Temperature - 1/p.T_ref));
    
    % Adjust permeability for conductivity/hydration state
    % As sigma decreases (degraded membrane), porosity/voids may increase
    hydration_factor = (props.sigma / p.sigma_0); 
    P_H2 = P_H2 * (1 + (1 - hydration_factor)); 
    
    %% 2. Calculate Driving Force (Pressure/Concentration)
    % J = P * (Delta_P / L_mem)
    % Delta_P is the partial pressure difference between anode and cathode
    delta_P = op.pressure * 1e5; % Convert bar to Pa
    
    %% 3. Calculate Fluxes (mol/m²·s)
    % FIX (previous bug): props.L_mem is already stored in SI units (meters,
    % e.g. delta_membrane_0 = 50e-6 m); the original code multiplied by an
    % additional 1e-6 here as if L_mem were still in micrometers, which
    % inflated the crossover flux by a factor of ~1e6 and made the coupled
    % ODE system numerically unsolvable (H2O2 generation source term became
    % unphysically large, ~1e7 mol/(m^3 s), versus realistic PEMFC
    % literature crossover-flux-equivalent current densities on the order
    % of 1-10 mA/cm^2, i.e. molar flux ~1e-8 to 1e-6 mol/(m^2 s); see e.g.
    % Kreuer, Chem. Rev. 104 (2004) 4637 and Karpenko-Jereb et al.,
    % Int. J. Hydrogen Energy).
    % J_H2: Hydrogen diffusing through the membrane
    J.H2 = (P_H2 * delta_P) / props.L_mem;

    % J_O2: Oxygen diffusing through the membrane
    % Oxygen permeability is related to H2 via gas diffusion ratios
    J.O2 = J.H2 * (sqrt(2 / 32)); % Molecular weight scaling
    
    %% 4. Apply Damage Multiplier
    % Localized pinhole damage (HF-induced) accelerates flux beyond diffusion
    damage_multiplier = 1 + (species.HF / p.HF_critical);
    J.H2 = J.H2 * damage_multiplier;
    J.O2 = J.O2 * damage_multiplier;
end