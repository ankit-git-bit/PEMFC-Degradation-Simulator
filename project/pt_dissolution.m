function d_pt = pt_dissolution(species, kinetics, p, op)
    % PT_DISSOLUTION.M
    % Mechanism: Platinum catalyst degradation via potential-driven
    % dissolution and re-deposition (Butler-Volmer type kinetics), and the
    % resulting evolution of electrochemically active surface area (ECSA).
    %
    % This is the SINGLE authoritative implementation of Pt kinetics in the
    % project; reaction_rates.m calls this function rather than duplicating
    % the dissolution/redeposition formulas inline.
    %
    % Functional form after the widely used Pt-dissolution frameworks of:
    %   - Darling & Meyers, J. Electrochem. Soc. 150 (2003) A1523
    %   - Rinaldo, Lee, Stumper & Eikerling, Electrochim. Acta 57 (2011) 273
    % Rate CONSTANTS (kinetics.Pt.k_diss_A, k_redep_A, Ea) are engineering
    % placeholders, NOT taken from a single literature source for this
    % project - see the TODO/ASSUMPTION note in kinetic_parameters.m.
    %
    % Inputs:
    %   species  - struct with fields Pt (metallic, mol/m^3), Pt2 (Pt2+
    %              ions, mol/m^3), ECSA (m^2/g_Pt)
    %   kinetics - struct from kinetic_parameters.m; uses kinetics.Pt and
    %              kinetics.coupling
    %   p        - struct from parameters.m; uses p.R, p.F, p.Pt_0, p.ECSA_0
    %   op       - struct from operating_conditions.m; uses op.Temperature,
    %              op.voltage (local electrode potential, V)
    %
    % Outputs (struct d_pt):
    %   dissolution  - forward dissolution rate (mol/(m^3 s))
    %   redeposition - reverse redeposition rate (mol/(m^3 s))
    %   net          - net dissolution rate (dissolution - redeposition)
    %   Pt           - d[Pt]/dt  (mol/(m^3 s))
    %   Pt2          - d[Pt2+]/dt (mol/(m^3 s))
    %   ECSA         - d[ECSA]/dt (m^2/g_Pt per s)

    T = op.Temperature;
    E = op.voltage;

    k_diss  = kinetics.Pt.k_diss_A  * exp(-kinetics.Pt.Ea / (p.R * T));
    k_redep = kinetics.Pt.k_redep_A * exp(-kinetics.Pt.Ea / (p.R * T));

    % Butler-Volmer potential dependence (forward = dissolution favored at
    % high anodic potential; reverse = redeposition of dissolved Pt2+)
    R_diss  = k_diss  * species.Pt  * exp( kinetics.Pt.alpha * kinetics.Pt.n * p.F * E / (p.R * T));
    R_redep = k_redep * species.Pt2 * exp(-(1 - kinetics.Pt.alpha) * kinetics.Pt.n * p.F * E / (p.R * T));

    d_pt.dissolution  = R_diss;
    d_pt.redeposition = R_redep;
    d_pt.net          = R_diss - R_redep;   % net loss rate of metallic Pt

    % Mass balance
    d_pt.Pt  = -d_pt.net;
    d_pt.Pt2 =  d_pt.net;

    % ECSA evolution: active surface area declines with net Pt dissolution,
    % scaled by the fraction of Pt already lost (captures the slowing of
    % ECSA loss as particles coarsen / surface-limited sites are consumed -
    % a simplified proxy for Ostwald-ripening-driven ECSA loss; see
    % Rinaldo et al. 2011 for the full particle-size-distribution treatment,
    % which is NOT implemented here - TODO for a future revision).
    if isfield(species, 'ECSA') && isfield(p, 'Pt_0') && p.Pt_0 > 0
        pt_loss_fraction = max(1 - species.Pt / p.Pt_0, 0);
        d_pt.ECSA = -kinetics.coupling.k_ECSA * species.ECSA * (1 + pt_loss_fraction);
    else
        d_pt.ECSA = 0;
    end
end
