function out = peroxide_generation(varargin)
    % PEROXIDE_GENERATION.M
    % Calculates peroxide generation from gas crossover, and (in full mode)
    % the associated radical/iron species derivatives for Mechanism 1.
    %
    % This function supports TWO calling conventions:
    %
    %   (A) FLUX MODE (4 args) - used internally by reaction_rates.m:
    %       S_H2O2 = peroxide_generation(species, props, op, p)
    %       Returns a SCALAR volumetric H2O2 generation rate [mol/(m^3 s)]
    %       computed from hydrogen_crossover() and catalyst ECSA.
    %
    %   (B) MECHANISM-AGGREGATOR MODE (5 args) - used for isolated
    %       Mechanism-1 analysis/diagnostics:
    %       d_mech1 = peroxide_generation(rates, species, props, op, p)
    %       Returns a STRUCT with derivatives for Fe2, Fe3, OH, OOH, H,
    %       H2O2, and diagnostic fields, given already-computed reaction
    %       rates R1-R13.
    %
    % Splitting the function this way avoids duplicating the crossover-flux
    % calculation logic (single source of truth: hydrogen_crossover.m) while
    % keeping the file backward compatible with both call sites.
    %
    % References: crossover functional form after Kreuer (Chem. Rev. 2004)
    % and Karpenko-Jereb et al. (Int. J. Hydrogen Energy) gas-permeation
    % framework; see hydrogen_crossover.m for details.

    if nargin == 4
        [species, props, op, p] = deal(varargin{1}, varargin{2}, varargin{3}, varargin{4});
        out = local_flux(species, props, op, p);
    elseif nargin == 5
        [rates, species, props, op, p] = deal(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
        out = local_aggregate(rates, species, props, op, p);
    else
        error('peroxide_generation:nargin', ...
            'peroxide_generation requires 4 (flux mode) or 5 (aggregator mode) arguments, got %d.', nargin);
    end
end

%% ------------------------------------------------------------------
function S_H2O2 = local_flux(species, props, op, p)
    % Volumetric H2O2 generation rate from H2/O2 crossover reacting at Pt.

    flux = hydrogen_crossover(species, props, op, p);

    % Stoichiometric efficiency factor for H2 + O2 -> H2O2 on Pt surfaces.
    % TODO/ASSUMPTION: the 2-electron (peroxide-forming) ORR pathway is a
    % minor side-channel relative to the dominant 4-electron pathway; real
    % H2O2 selectivity at Pt is typically low (see Liu & Zuckerbrod,
    % J. Electrochem. Soc. 152 (2005) A1165). The value below is an
    % ORDER-OF-MAGNITUDE calibration chosen so that the resulting
    % steady-state membrane H2O2 concentration falls in the literature-
    % reported micromolar range (see e.g. reviews of PEMFC chemical
    % degradation, Zhao, Adzic and coworkers), rather than a value fit to
    % a specific measured H2O2 dataset for this project. An earlier,
    % un-calibrated value (0.5) produced steady-state H2O2 several orders
    % of magnitude above any literature-reported concentration, which in
    % turn produced unphysically fast (second-scale) membrane side-chain
    % loss - flagged here as a corrected engineering assumption, not a
    % validated literature constant.
    k_gen = 5e-6;

    % Scale by available catalyst surface area (normalized to BOL ECSA)
    if isfield(species, 'ECSA') && isfield(p, 'ECSA_0') && p.ECSA_0 > 0
        ecsa_ratio = species.ECSA / p.ECSA_0;
    else
        ecsa_ratio = 1.0; % No ECSA state tracked in this call context
    end

    % Volumetric generation rate [mol/(m^3 s)]: convert areal flux to
    % volumetric rate by dividing by membrane thickness.
    S_H2O2 = k_gen * (flux.H2 + flux.O2) * ecsa_ratio / max(props.L_mem, eps);

    % Higher cathode RH helps wash away radicals / limits local H2O2
    % accumulation (TODO/ASSUMPTION: linear (1-RH) scaling is a simplified
    % engineering approximation, not a fitted literature correlation).
    if isfield(op, 'RH_cathode')
        S_H2O2 = S_H2O2 * (1 - op.RH_cathode);
    end
end

%% ------------------------------------------------------------------
function d_mech1 = local_aggregate(rates, species, props, op, p)
    % Full Mechanism-1 derivative struct (iron redox + radicals + H2O2),
    % for standalone/diagnostic use. Reuses local_flux() for S_H2O2 so the
    % crossover-generation term is never duplicated.

    %% 1. Iron Redox Cycle (R1 - R5)
    d_mech1.Fe2 = -rates.R1 + rates.R2 - rates.R3 - rates.R4 + rates.R5;
    d_mech1.Fe3 =  rates.R1 - rates.R2 + rates.R3 + rates.R4 - rates.R5;

    %% 2. Radical Generation and Consumption (R1 - R13)
    d_mech1.OH  =  rates.R1 + 2*rates.R6 - rates.R3 - rates.R7 + rates.R8 ...
                  - 2*rates.R10 - rates.R11 - rates.R12;

    d_mech1.OOH =  rates.R2 - rates.R4 - rates.R5 + rates.R7 - rates.R8 ...
                  - 2*rates.R9 - rates.R11 + rates.R13;

    d_mech1.H   =  rates.R12 - rates.R13;

    %% 3. Peroxide Generation from Gas Crossover
    d_mech1.S_H2O2 = local_flux(species, props, op, p);

    %% 4. Peroxide Consumption Terms (Fenton reactions R1, R2, R4, R6, R7)
    d_mech1.H2O2 = d_mech1.S_H2O2 - rates.R1 - rates.R2 - rates.R4 ...
                   - rates.R6 - rates.R7;

    %% 5. Diagnostics
    d_mech1.total_radicals = d_mech1.OH + d_mech1.OOH + d_mech1.H;
end
