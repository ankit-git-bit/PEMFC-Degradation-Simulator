function op_updated = bidirectional_coupling(props, op, p)
    % BIDIRECTIONAL_COUPLING.M
    % Closes the bidirectional degradation feedback loop described in the
    % project architecture diagram:
    %
    %   ... -> HF release -> IEC loss -> Conductivity loss -> membrane
    %   thinning -> Hydrogen crossover increase -> Peroxide increase -> ...
    %
    % degradation_update.m computes the CURRENT physical state (props,
    % including props.crossover_H2) FROM the chemistry. This function does
    % the other half of the loop: it feeds the updated crossover state back
    % INTO the operating-condition boundary values (op.H2) that drive the
    % NEXT macro-timestep's chemistry (R12 in reaction_rates.m, and the
    % crossover-driven H2O2 source term in peroxide_generation.m). This is
    % the ONLY module that writes to op.H2/op.O2 based on degradation
    % state; reaction_rates.m and peroxide_generation.m only read them.
    %
    % Inputs:
    %   props - current membrane property struct (from degradation_update.m),
    %           must contain crossover_H2
    %   op    - current operating-condition struct
    %   p     - parameters.m struct (baseline permeability reference)
    %
    % Output:
    %   op_updated - op with H2 (and O2, scaled by the same relative
    %                permeability increase) rescaled by the ratio of
    %                current-to-baseline crossover permeability.

    op_updated = op;

    if isfield(props, 'crossover_H2') && p.permeability_H2 > 0
        crossover_ratio = props.crossover_H2 / p.permeability_H2;
        crossover_ratio = max(crossover_ratio, 0); % physical floor

        op_updated.H2 = op.H2 * crossover_ratio;
        % TODO/ASSUMPTION: O2 crossover is assumed to scale with the same
        % relative degradation factor as H2 (both driven by the same
        % membrane thinning / pinhole damage mechanism). A rigorous
        % treatment would track O2 permeability independently via its own
        % Arrhenius/damage model (see hydrogen_crossover.m J.O2 for the
        % BOL molecular-weight-scaling estimate) - simplified here to avoid
        % double-counting damage multipliers.
        op_updated.O2 = op.O2 * crossover_ratio;
    end
end
