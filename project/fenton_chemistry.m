function d_fenton = fenton_chemistry(rates)
    % FENTON_CHEMISTRY.M
    % Implements the radical propagation and termination reactions (R8-R13).
    % Purpose: Calculate the derivative contributions for the ROS species 
    % based on the kinetic framework from Frühwirt et al. (2020).
    
    % R8: OOH + H2O2 -> OH + H2O + O2
    % R9: 2OOH -> H2O2
    % R10: 2OH -> H2O2
    % R11: OOH + OH -> H2O + O2
    % R12: OH + H2 -> H2O + H
    % R13: H + O2 -> OOH

    %% 1. Contribution to H2O2
    d_fenton.H2O2 = -rates.R8 + rates.R9 + rates.R10;

    %% 2. Contribution to OH
    d_fenton.OH = rates.R8 - 2*rates.R10 - rates.R11 - rates.R12;

    %% 3. Contribution to OOH[cite: 5]
    d_fenton.OOH = -rates.R8 - 2*rates.R9 - rates.R11 + rates.R13;

    %% 4. Contribution to H (Hydrogen Radical)[cite: 5]
    d_fenton.H = rates.R12 - rates.R13;

    %% 5. Iron Species 
    % Note: Fenton cycles (R1-R5) are handled in peroxide_generation.m.
    % This module only handles the specific propagation rates requested.
    d_fenton.Fe2 = 0; 
    d_fenton.Fe3 = 0;
end