%% Calculating Z0

clear; close all; clc;

% Paths
wave_params = readcell('Offshore_Waves.xlsx');
means_path  = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means';
BL_path     = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer';

%%
% Case
caze  = 'WT6_WV0_AGP';
means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
means = means.output;
BL = load(fullfile(BL_path, strcat(caze, '_BL.mat')));
BL = BL.output;

% Get wave parameters
% wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
% wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
% ak = wave_amplitude / wave_length;

% Coordinates
x  = unique(means.X) * 1E-3;
z  = flipud(unique(means.Y)) * 1E-3;
u  = mean(means.ensemble.u, 2, 'omitnan');
uv = mean(means.ensemble.uv, 2, 'omitnan');
delta = mean(BL.ensemble.thickness, 'omitnan') * 1E-3;

% Normalize
kinematic_viscosity = 1.48E-5;
u_star = sqrt(max(-1 * uv));
u_plus = u / u_star;
z_delta = z / delta;

% Reynolds number
Re = (u_star * delta / kinematic_viscosity);

% Mask outside of boundary layer
BL_mask = z_delta < 1;
u_plus = u_plus(BL_mask);
z_delta = z_delta(BL_mask);

% Mask nans
nan_mask = ~isnan(u_plus);
u_plus = u_plus(nan_mask);
z_delta = z_delta(nan_mask);

% Which Points to Fit to
upper_idx = 1;
lower_idx = 6;


% Fitting
uplus = double(u_plus(upper_idx:lower_idx));  % this has to be in plus units
zz = double(z_delta(upper_idx:lower_idx));    %vertical location of the velocities between 0.03h - 0.1h. this has to be in any units

% Define the von Kármán constant
kappa = 0.40;

% Formulate the law of the wall equation
equation = @(zo, z) (1 / kappa) * log(z / zo);

% Initial guess for the reference height (y_0)
% Adjust this based on your knowledge or estimation
initial_guess_z0 = 1.5e-4; 

% Perform the least squares fit
coefficients = lsqcurvefit(equation, initial_guess_z0, zz, uplus);

% Extract the reference height (y_0) from the coefficients
ensemble_z0 = coefficients(1);
fprintf('\nz0 = %3.2e\n', ensemble_z0)

% Plot original data
ax = figure();
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
a  = gca;
set(gca, 'XScale', 'log')
grid on;
hold off;

% Add z0 to plot
dim = [.2 .5 .3 .3];
str = sprintf('$z_0 = %3.2e$', ensemble_z0);
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'Interpreter', 'latex', 'fontsize', 16);

%%




% %% Test across different x locations
% 
% % Coordinates
% x  = unique(means.X) * 1E-3;
% z  = flipud(unique(means.Y)) * 1E-3;
% u  = means.ensemble.u;
% uv = means.ensemble.uv;
% delta = BL.ensemble.thickness * 1E-3;
% kinematic_viscosity = 1.48E-5;
% 
% c = 1;
% for i = 3:length(x)-3
%     disp(i)
% 
%     % Index Profiles
%     u_profile = u(:,i);
%     uv_profile = uv(:,i);
%     u_star = sqrt(max(-1 * uv_profile));
%     u_plus = u_profile / u_star;
%     z_delta = z / delta(i);
% 
%     % Mask outside of boundary layer
%     BL_mask = z_delta < 1;
%     u_plus = u_plus(BL_mask);
%     z_delta = z_delta(BL_mask);
% 
%     % Mask nans
%     nan_mask = ~isnan(u_plus);
%     u_plus = u_plus(nan_mask);
%     z_delta = z_delta(nan_mask);
% 
%     % Which Points to Fit to
%     upper_idx = 1;
%     lower_idx = 6;
% 
%     % Fitting
%     uplus = double(u_plus(upper_idx:lower_idx));
%     zz = double(z_delta(upper_idx:lower_idx));
%     coefficients = lsqcurvefit(equation, initial_guess_z0, zz, uplus);
%     z0 = coefficients(1);
% 
%     z0_values(c) = z0;
%     c = c + 1;
% end
% 
% figure()
% plot(x(3:end-3), z0_values)
% yline(ensemble_z0)
























