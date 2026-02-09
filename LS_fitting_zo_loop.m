%% Calculating Z0 Looped

%% Imports
clear; close all; clc;

% Paths
wave_params = readcell('Offshore_Waves.xlsx');
means_path  = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means';
BL_path     = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer';
save_path   = '/Users/zeinsadek/Desktop/Experiments/Offshore/Hopkins/Z0';

%% Setup Cases

% Where to Fit
upper_idx = 1;
lower_idx = 6;

% Gravity
g = 9.81;

% Von Kármán constant
kappa = 0.40;

% Viscosity
kinematic_viscosity = 1.48E-5;

% Formulate the Law of the Wall Equation
equation = @(zo, z) (1 / kappa) * log(z / zo);

% Initial guess for the reference height (y_0)
initial_guess_z0 = 1.5e-4; 


% Case Conditions
WTs = {'4', '6', '8'};
WVs = {'0', 'A', 'B', 'C', 'D'};

for i = 1:length(WTs)
    % Wind Speed
    WT = WTs{i};

    for j = 1:length(WVs)
        % Wave Type
        WV = WVs{j};

        % Active Grid
        if strcmp(WV, '0') == 1
            AG = '_AGP';
        else
            AG = '_AG0';
        end
        
        % Case Data
        caze = strcat('WT', WT, '_WV', WV, AG);
        means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
        means = means.output;
        BL = load(fullfile(BL_path, strcat(caze, '_BL.mat')));
        BL = BL.output;

        % Wave Parameters
        if strcmp(WV, '0') == 0
            wave_length = wave_params{find(strcmp(wave_params, WV) == 1), 2} * 1E-3;
            wave_amplitude = wave_params{find(strcmp(wave_params, WV) == 1), 3} * 1E-3;
            k  = 1 / wave_length;
            ak = wave_amplitude / wave_length;
            output.(caze).ak = ak;
        else
            wave_length = nan;
            wave_amplitude = nan;
            k  = 1 / wave_length;
            ak = nan;
            output.(caze).ak = ak;
        end

        % Coordinates
        x  = unique(means.X) * 1E-3;
        z  = flipud(unique(means.Y)) * 1E-3;
        u  = mean(means.ensemble.u, 2, 'omitnan');
        uv = mean(means.ensemble.uv, 2, 'omitnan');
        delta = mean(BL.ensemble.thickness, 'omitnan') * 1E-3;

        % Normalize
        u_star = sqrt(max(-1 * uv));
        u_plus = u / u_star;
        z_delta = z / delta;
        output.(caze).u_star = u_star;

        % C+
        c = sqrt(g / k);
        c_plus = c / u_star;
        output.(caze).c_plus = c_plus;
        
        % Reynolds number
        Re = (u_star * delta / kinematic_viscosity);
        output.(caze).Re = Re;
        
        % Mask outside of boundary layer
        BL_mask = z_delta < 1;
        u_plus = u_plus(BL_mask);
        z_delta = z_delta(BL_mask);
        
        % Mask nans
        nan_mask = ~isnan(u_plus);
        u_plus = u_plus(nan_mask);
        z_delta = z_delta(nan_mask);

        % Save Masked Values
        output.(caze).u_plus = u_plus;
        output.(caze).z_delta = z_delta;

        % Fitting
        uplus = double(u_plus(upper_idx:lower_idx)); 
        zz = double(z_delta(upper_idx:lower_idx));
        coefficients = lsqcurvefit(equation, initial_guess_z0, zz, uplus);
        ensemble_z0 = coefficients(1);
        output.(caze).z0 = ensemble_z0;
        fprintf('\nz0 = %3.2e\n', ensemble_z0)

        % Plot original data
        ax = figure('name', 'caze');
        hold on;
        plot(zz, uplus, 'ko-', 'DisplayName', 'Fit Data Points');
        
        % Plot fitted curve
        fitted_curve = equation(ensemble_z0, zz);
        plot(zz, fitted_curve, 'r-', 'DisplayName', 'Fitted Curve');
        
        % Add labels and legend
        title(caze, 'Interpreter', 'none')
        xlabel('$y / \delta$', 'FontSize', 14, 'interpreter', 'latex');
        ylabel('U^+', 'FontSize', 14);
        legend('show');
        
        % Show grid and hold off
        set(gca, 'XScale', 'log')
        % xlim([0.9, 1])
        % ylim([13, 27])
        grid on;
        hold off;
        
        % Add z0 to plot
        dim = [.2 .5 .3 .3];
        str = sprintf('$z_0 = %3.2e$', ensemble_z0);
        annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'Interpreter', 'latex', 'fontsize', 16);
        % exportgraphics(ax, fullfile(save_path, 'figures', strcat(caze, '_z0.png')), 'Resolution', 200);
        close all

    end
end

save(fullfile(save_path, 'PSU_z0.mat'), 'output');


 
























