function d_mech3 = polymer_unzipping(rates, species, kinetics, op)
    % POLYMER_UNZIPPING.M
    % Mechanism 3: Backbone unzipping sequence[cite: 5].
    % Uses ONLY reaction rates R17 to R23.
    % Outputs the partial derivatives for HF, CO2, chain length, 
    % and Fluorine Evolution Rate (FER).
    %
    % Inputs:
    %   rates   - Structure containing reaction rates R17 to R23
    %   species - Structure containing species concentrations
    %   kinetics- Structure containing kinetic parameters
    %   op      - Operating conditions structure
    %
    % Outputs:
    %   d_mech3 - Structure containing derivatives for:
    %             HF, CO2, CF2 chain sequence, chain length,
    %             FER, degradation fraction, and diagnostic metrics

    %% 1. Cumulative Degradation Products (Mechanism 3 contributions only)[cite: 5]
    % Note: HF and CO2 are also produced in Mechanism 2 (R15, R16)
    d_mech3.HF  = 2 * (rates.R17 + rates.R18 + rates.R19 + rates.R20 ...
                     + rates.R21 + rates.R22 + rates.R23);
                 
    d_mech3.CO2 = rates.R17 + rates.R18 + rates.R19 + rates.R20 ...
                + rates.R21 + rates.R22 + rates.R23;

    %% 2. CF2 Chain Unzipping Cascade[cite: 5]
    % Note: CF2_7 formation from R16 is handled in Mechanism 2.
    % This section maps the consumption and propagation down the chain.
    d_mech3.CF2_7 = -rates.R17;
    d_mech3.CF2_6 =  rates.R17 - rates.R18;
    d_mech3.CF2_5 =  rates.R18 - rates.R19;
    d_mech3.CF2_4 =  rates.R19 - rates.R20;
    d_mech3.CF2_3 =  rates.R20 - rates.R21;
    d_mech3.CF2_2 =  rates.R21 - rates.R22;
    d_mech3.CF2_1 =  rates.R22 - rates.R23;
    
    %% 3. Chain Length Evolution
    % Calculate total number of CF2 units remaining
    if isfield(species, 'CF2_7') && isfield(species, 'CF2_1')
        % Current total CF2 units
        current_CF2_total = species.CF2_7 + species.CF2_6 + species.CF2_5 + ...
                            species.CF2_4 + species.CF2_3 + species.CF2_2 + species.CF2_1;
        
        % Initial total CF2 units (from parameters)
        if isfield(species, 'CF2_total_0')
            initial_CF2_total = species.CF2_total_0;
        else
            initial_CF2_total = current_CF2_total + 0.1; % Estimate if not provided
        end
        
        % Average chain length (number of CF2 units per chain)
        % Assume each chain starts with 7 CF2 units (from CF2_7)
        if isfield(species, 'n_chains')
            n_chains = species.n_chains;
        else
            % Estimate number of chains from initial CF2_7 concentration
            n_chains = species.CF2_7 / 0.01; % Rough estimate
            n_chains = max(n_chains, 1e-10);
        end
        
        % Average chain length
        d_mech3.chain_length_avg = current_CF2_total / n_chains;
        
        % Chain length fraction (remaining relative to initial)
        d_mech3.chain_length_fraction = current_CF2_total / initial_CF2_total;
        
        % Rate of chain length change
        d_mech3.chain_length_derivative = (d_mech3.CF2_7 + d_mech3.CF2_6 + ...
            d_mech3.CF2_5 + d_mech3.CF2_4 + d_mech3.CF2_3 + ...
            d_mech3.CF2_2 + d_mech3.CF2_1) / n_chains;
    else
        % Default values if species not available
        d_mech3.chain_length_avg = 7; % Maximum chain length
        d_mech3.chain_length_fraction = 1.0;
        d_mech3.chain_length_derivative = 0;
    end
    
    %% 4. Fluorine Evolution Rate (FER)
    % FER is the rate of fluoride ion release from the membrane
    % Each HF molecule released corresponds to one fluoride ion
    
    % Instantaneous FER (mol/m³·s)
    d_mech3.FER_instantaneous = d_mech3.HF;
    
    % Cumulative FER (total fluoride released over time)
    if isfield(species, 'FER_cumulative')
        d_mech3.FER_cumulative = species.FER_cumulative + d_mech3.HF;
    else
        d_mech3.FER_cumulative = d_mech3.HF;
    end
    
    % Normalized FER (per unit membrane volume)
    if isfield(species, 'L_mem') && species.L_mem > 0
        d_mech3.FER_per_volume = d_mech3.HF / species.L_mem;
    else
        d_mech3.FER_per_volume = d_mech3.HF / 1e-4; % Default thickness
    end
    
    % FER per unit area (mol/m²·s)
    if isfield(species, 'L_mem') && species.L_mem > 0
        d_mech3.FER_per_area = d_mech3.HF * species.L_mem;
    else
        d_mech3.FER_per_area = d_mech3.HF * 1e-4;
    end
    
    % FER as fraction of total fluorine content
    if isfield(species, 'fluorine_total_0') && species.fluorine_total_0 > 0
        d_mech3.FER_fraction = d_mech3.HF / species.fluorine_total_0;
    else
        % Estimate total fluorine from initial CF2 concentration
        initial_CF2_total = species.CF2_7 + species.CF2_6 + species.CF2_5 + ...
                            species.CF2_4 + species.CF2_3 + species.CF2_2 + species.CF2_1;
        total_fluorine = initial_CF2_total * 2; % Each CF2 has 2 F atoms
        d_mech3.FER_fraction = d_mech3.HF / total_fluorine;
    end

    %% 5. Backbone Degradation Metrics
    % Total backbone degradation rate (sum of all unzipping steps)
    d_mech3.total_unzipping_rate = rates.R17 + rates.R18 + rates.R19 + rates.R20 + ...
                                   rates.R21 + rates.R22 + rates.R23;
    
    % Weighted average unzipping rate (accounting for chain length)
    if d_mech3.total_unzipping_rate > 0
        weighted_sum = 7*rates.R17 + 6*rates.R18 + 5*rates.R19 + 4*rates.R20 + ...
                       3*rates.R21 + 2*rates.R22 + 1*rates.R23;
        d_mech3.avg_unzipping_step = weighted_sum / d_mech3.total_unzipping_rate;
    else
        d_mech3.avg_unzipping_step = 0;
    end
    
    % Backbone degradation fraction
    if isfield(species, 'CF2_total_0') && species.CF2_total_0 > 0
        current_CF2_total = species.CF2_7 + species.CF2_6 + species.CF2_5 + ...
                            species.CF2_4 + species.CF2_3 + species.CF2_2 + species.CF2_1;
        d_mech3.backbone_degradation_fraction = 1 - (current_CF2_total / species.CF2_total_0);
    else
        d_mech3.backbone_degradation_fraction = 0;
    end

    %% 6. Radical Attack Efficiency
    % Efficiency of OH radicals in causing unzipping
    if isfield(rates, 'OH_total') && rates.OH_total > 0
        d_mech3.unzipping_efficiency = d_mech3.total_unzipping_rate / rates.OH_total;
    else
        d_mech3.unzipping_efficiency = 0.1; % Default assumption
    end
    
    % OH consumption in unzipping reactions
    % Each unzipping step consumes 2 OH radicals
    d_mech3.OH_consumption = 2 * d_mech3.total_unzipping_rate;

    %% 7. Chain Scission Frequency
    % Number of chain scission events per unit time
    % Each unzipping step represents one scission event
    d_mech3.scission_frequency = d_mech3.total_unzipping_rate;
    
    % Scission frequency per chain
    if isfield(species, 'n_chains') && species.n_chains > 0
        d_mech3.scissions_per_chain = d_mech3.total_unzipping_rate / species.n_chains;
    else
        d_mech3.scissions_per_chain = d_mech3.total_unzipping_rate / 1e-6;
    end

    %% 8. Molecular Weight Evolution
    % Estimate molecular weight reduction due to unzipping
    if isfield(species, 'MW_initial') && species.MW_initial > 0
        % Initial molecular weight (g/mol)
        MW_0 = species.MW_initial;
        
        % Current molecular weight based on chain length fraction
        if isfield(d_mech3, 'chain_length_fraction')
            d_mech3.MW_current = MW_0 * d_mech3.chain_length_fraction;
        else
            d_mech3.MW_current = MW_0;
        end
        
        % Rate of molecular weight change
        d_mech3.MW_derivative = -MW_0 * (1 - d_mech3.chain_length_fraction);
    else
        d_mech3.MW_current = 100000; % Default MW for Nafion (g/mol)
        d_mech3.MW_derivative = 0;
    end

    %% 9. Cumulative Degradation Products
    % Total HF produced from unzipping (cumulative)
    if isfield(species, 'HF_total_unzipping')
        d_mech3.HF_cumulative = species.HF_total_unzipping + d_mech3.HF;
    else
        d_mech3.HF_cumulative = d_mech3.HF;
    end
    
    % Total CO2 produced from unzipping (cumulative)
    if isfield(species, 'CO2_total_unzipping')
        d_mech3.CO2_cumulative = species.CO2_total_unzipping + d_mech3.CO2;
    else
        d_mech3.CO2_cumulative = d_mech3.CO2;
    end

    %% 10. Diagnostic Metrics
    % Polymer degradation index (0 = intact, 1 = fully degraded)
    if isfield(d_mech3, 'chain_length_fraction')
        d_mech3.degradation_index = 1 - d_mech3.chain_length_fraction;
    else
        d_mech3.degradation_index = 0;
    end
    
    % Average chain length (number of CF2 units)
    d_mech3.avg_chain_length = d_mech3.chain_length_avg;
    
    % Chain length distribution (if available)
    if isfield(species, 'CF2_7') && isfield(species, 'CF2_1')
        d_mech3.chain_distribution = [
            species.CF2_7, species.CF2_6, species.CF2_5, ...
            species.CF2_4, species.CF2_3, species.CF2_2, species.CF2_1
        ];
    end
    
    % Total CF2 loss rate
    d_mech3.CF2_loss_rate = -(d_mech3.CF2_7 + d_mech3.CF2_6 + d_mech3.CF2_5 + ...
                               d_mech3.CF2_4 + d_mech3.CF2_3 + d_mech3.CF2_2 + d_mech3.CF2_1);
    
    % Ensure all fields are properly assigned
    if ~isfield(d_mech3, 'FER_instantaneous')
        d_mech3.FER_instantaneous = 0;
    end
    if ~isfield(d_mech3, 'FER_cumulative')
        d_mech3.FER_cumulative = 0;
    end
    if ~isfield(d_mech3, 'degradation_index')
        d_mech3.degradation_index = 0;
    end
end