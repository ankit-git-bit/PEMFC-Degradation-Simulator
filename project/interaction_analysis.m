function interaction_analysis()
    % INTERACTION_ANALYSIS.M
    % Creates an interaction matrix and directed graph (digraph) 
    % mapping the cause-and-effect chain of chemical degradation[cite: 1, 5].

    %% 1. Define Cause-and-Effect Pairs
    % Each source-target pair represents a directed mechanism link[cite: 1, 5].
    
    sources = {
        'H2_Crossover', ...
        'O2_Crossover', ...
        'H2O2', ...
        'Fe2_Ion', ...
        'OH_Radical', ...
        'OH_Radical', ...
        'SC_SO3H', ...
        'SC_O_Radical', ...
        'BB_O_Radical', ...
        'CF2_Chain', ...
        'CF2_Chain', ...
        'SC_O_Radical', ...
        'HF'
    };

    targets = {
        'H2O2', ...             % Reactant crossover forms peroxide[cite: 1, 4]
        'H2O2', ...             % Reactant crossover forms peroxide[cite: 1, 4]
        'OH_Radical', ...       % Fenton reaction splits H2O2[cite: 1, 5]
        'OH_Radical', ...       % Catalyzed by Fe2+ impurity[cite: 1, 5]
        'SC_SO3H', ...          % Hydroxyl attack initiates side-chain split[cite: 1, 5]
        'CF2_Chain', ...        % Hydroxyl attack drives backbone unzipping[cite: 1, 5]
        'SC_O_Radical', ...     % Side-chain yields SC radical[cite: 5]
        'BB_O_Radical', ...     % SC breakdown leaves backbone radical[cite: 5]
        'CF2_Chain', ...        % BB radical initiates unzipping cascade[cite: 1, 5]
        'HF', ...               % Unzipping emits Hydrogen Fluoride[cite: 1, 5]
        'CO2', ...              % Unzipping emits Carbon Dioxide[cite: 1, 5]
        'HF', ...               % Side-chain fragmenting emits HF[cite: 5]
        'FER'                   % HF accumulation dictates Fluoride Emission Rate[cite: 3, 5]
    };

    %% 2. Generate Directed Graph (digraph)
    G = digraph(sources, targets);

    %% 3. Plot Interaction Graph
    figure('Color', 'w', 'Name', 'Chemical Degradation Interaction Matrix');
    
    % Layered layout naturally visualizes the chronological cascade
    p = plot(G, 'Layout', 'layered', ...
                'NodeColor', [0.2 0.4 0.6], ...
                'EdgeColor', [0.8 0.2 0.2], ...
                'LineWidth', 1.5, ...
                'MarkerSize', 7, ...
                'ArrowSize', 12);
            
    % Formatting properties for presentation
    title('Interaction Matrix: PEMFC Chemical Degradation Mechanisms', 'FontSize', 14);
    p.NodeFontSize = 10;
    p.NodeFontWeight = 'bold';
    axis off;

    %% 4. Extract and Display Mathematical Interaction Matrix (Adjacency)
    % This binary matrix represents 1 for an active interaction, 0 for no interaction.
    interaction_matrix = full(adjacency(G));
    node_names = G.Nodes.Name;
    
    disp('======================================================');
    disp('            CHEMICAL INTERACTION MATRIX');
    disp('======================================================');
    disp('Rows = Mechanism Source | Columns = Mechanism Target');
    disp(' ');
    
    matrix_table = array2table(interaction_matrix, ...
                               'VariableNames', node_names, ...
                               'RowNames', node_names);
    disp(matrix_table);
end