%%% Phase Averaging Code
% Zein Sadek, 1/23

function output = phaseaverage_Polytechnique(constants, data)


    % Reference Wave
    percent_phase_tolerance = constants.percent_phase_tolerance;
    wave_amplitude  = constants.amplitude;
    wavelength      = constants.wavelength;
    phase_offset    = constants.phase_offset;    
    phase_tolerance = percent_phase_tolerance * wavelength;
    
    % Load Wave Profiles    
    wave_profiles = data.waves;
    phase_fits = data.phases;
    D = data.D;
    X = data.X;
    x = X(1,:);
    
    % Saves
    phase_average_idx = zeros(1, D);
    
    fprintf('<phaseaverage> PROGRESS: \n');
    for frame_number = 1:D
        
        % Print Progress.
        progressbarText(frame_number/D);
        
        % Fitted Phase Offest
        fitted_phase = phase_fits(frame_number);
        
        % Bin frames to phases
        for phase = 1:length(phase_offset)   
            lower_phase = phase_offset(phase) - phase_tolerance;
            upper_phase = phase_offset(phase) + phase_tolerance;

            if fitted_phase >= lower_phase && fitted_phase <= upper_phase
                phase_average_idx(frame_number) = phase;
            else
            end
        end 
    end
    
    %%% OUTPUT
    output = phase_average_idx;


    %% Plot
    clc;
    fprintf('\nTolerance Range: %4.2f mm\n', 2 * phase_tolerance)

    % Create tiled layout plot with all 4 phases
    figure('color', 'white')
    fig = tiledlayout(length(phase_offset),1, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Title include recording name and phase tolerance
    main_title = strcat(constants.recording_name, ': Tolerance +/- ', num2str(percent_phase_tolerance * 100), '% of Wavelength ');
    title(fig, main_title, 'Interpreter', 'none')
    
    % Iterate through phases
    for phase = 1:length(phase_offset)          
        ax(phase) = nexttile;
        hold on 
        
        % Search for frames with match with current phase
        idx = find(phase_average_idx == phase);
        fprintf('\nPhase %.f: %.f Images\n', phase, length(idx))
        
        % Loop through mateched frames and plot wave profiles
        for i = 1:length(idx)
            p = plot(x, wave_profiles(idx(i), :), 'Color', 'black', 'HandleVisibility','off');
            p.Color(4) = 0.25;
        end

        % Reference Wave
        reference_profile = wave_amplitude * cos(2 * pi * (x - (range(x) / 2) - phase_offset(phase)) / wavelength);
        plot(x, reference_profile, 'Color', 'red', 'LineWidth', 2);

         % Plot Max Wave Profile
        plot(x, max(wave_profiles(idx, :), [], 1), 'Color', 'green', 'LineWidth', 2)

        xline(118, 'LineStyle', '--')
        hold off  
        title(strcat('Phase', {' '}, num2str(phase),':', {' '}, num2str(length(idx)), ' Images'))
        ylabel('y [mm]')
        axis equal
    end
    
    % Match axes
    linkaxes(ax, 'xy')
    xlim([0, max(x)])
    ylim([-20, 20])

    if phase == length(phase_offset)
        xlabel('x [mm]')
    end
end


