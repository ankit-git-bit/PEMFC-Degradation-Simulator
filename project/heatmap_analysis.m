function heatmap_analysis(sens)
    % Dynamically extract the final HF concentration from each sensitivity run
    % to create a meaningful matrix for visualization.
    
    % Get the number of runs for each parameter
    n_T  = length(sens.Temperature);
    n_Fe = length(sens.Fe);
    
    % Initialize a grid (assuming we want to correlate T vs Fe)
    % If ranges are different lengths, we use a shared grid or zero-pad
    heatmap_data = zeros(n_T, n_Fe);
    
    for i = 1:n_T
        for j = 1:n_Fe
            % Extract final HF concentration (index 9 in the state vector)
            % This assumes the state vector index 9 is HF
            hf_at_end = sens.Fe(j).State(end, 9);
            heatmap_data(i, j) = hf_at_end;
        end
    end
    
    % Create the heatmap
    figure('Color', 'w');
    h = heatmap(heatmap_data);
    h.Title = 'Sensitivity Heatmap: HF Accumulation (T vs Fe)';
    h.XLabel = 'Fe Concentration (ppm)';
    h.YLabel = 'Temperature (K)';
    h.XDisplayLabels = string({sens.Fe.Value});
    h.YDisplayLabels = string({sens.Temperature.Value});
end