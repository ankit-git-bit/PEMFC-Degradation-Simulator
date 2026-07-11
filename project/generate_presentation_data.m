function generate_presentation_data()
    % GENERATE_PRESENTATION_DATA.M
    % Automates the export of all 10 project visual outputs.
    
    mkdir('presentation_exports'); % Create export folder
    
    %% 1. Species Evolution
    species_evolution(); 
    saveas(gcf, 'presentation_exports/1_Species_Evolution.png');
    
    %% 2. Mechanism Contribution (Time-scale Analysis)
    time_scale_analysis(); 
    saveas(gcf, 'presentation_exports/2_Mechanism_Contribution.png');
    
    %% 3. Lifetime Simulation
    history = lifetime_simulation(); % Ensure this returns result data
    save('presentation_exports/3_Lifetime_Data.mat', 'history');
    
    %% 4. Parameter Sensitivity
    sens = parameter_sensitivity();
    saveas(gcf, 'presentation_exports/4_Sensitivity.png');
    
    %% 5. Mechanism Ranking
    results = mechanism_ranking();
    saveas(gcf, 'presentation_exports/5_Mechanism_Ranking.png');
    
    %% 6. Interaction Network
    interaction_analysis();
    saveas(gcf, 'presentation_exports/6_Interaction_Network.png');
    
    %% 7. Timeline (Custom Timeline Plot)
    figure('Color', 'w');
    timeline(history); % Custom function defined below
    saveas(gcf, 'presentation_exports/7_Timeline.png');
    
    %% 8. Waterfall Chart
    figure('Color', 'w');
    waterfall_chart(results); % Custom function defined below
    saveas(gcf, 'presentation_exports/8_Waterfall.png');
    
    %% 9. Heatmap
    figure('Color', 'w');
    heatmap_analysis(sens); % Custom function defined below
    saveas(gcf, 'presentation_exports/9_Heatmap.png');
    
    %% 10. Tornado Chart
    % Generated within mechanism_ranking() but saved here for consistency
    saveas(gcf, 'presentation_exports/10_Tornado.png');
    
    close all;
    disp('All 10 visual outputs exported to /presentation_exports/');
end

function timeline(history)
    plot(history.time_hr, history.species(:, 9), 'LineWidth', 2);
    title('Degradation Timeline: HF Accumulation', 'FontWeight', 'bold');
    ylabel('HF [mol/m^3]'); xlabel('Operational Time (Hours)');
    grid on;
end

function waterfall_chart(results)
    % Cumulative sum of impacts
    data = [results.Baseline.SC, results.Baseline.HF]; 
    waterfall(data); % Uses standard plotting or 'waterfall' tool
    title('Waterfall: Cumulative Mass Loss Impacts', 'FontWeight', 'bold');
end

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
