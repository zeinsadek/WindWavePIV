%% Compute the integrated boundary layer equaiton for phase averages
% Zein Sadek
% 11/2023

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
% wave_parameters = readcell('Offshore_Waves.xlsx');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case
experiment_name = 'WT4_WV0_AGP';

% Experiment Specifics
tunnel_freq = experiment_name(strfind(experiment_name, 'WT') + 2);
wave_type   = experiment_name(strfind(experiment_name, 'WV') + 2);

% Save paths
results_path  = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/';
mean_file     = strcat(results_path, 'means'  , '/', experiment_name, '_MEANS.mat');
BL_file       = strcat(results_path, 'boundary_layer/', experiment_name, '_BL.mat');

% Wave parameters
clear wave_parameters

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

means = load(mean_file); 
means = means.output;

BL_data = load(BL_file);
BL_data = BL_data.output;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PICK LOCATION FOR SLICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinates
Y = means.Y;
y = Y(:,1);
X = means.X;
x = X(1,:);

% Constants
nu = 1.48E-5;
dy = ((y(3) - y(2)) * 10^(-3));
dx = ((x(3) - x(2)) * 10^(-3));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPAIR CAT-SCRATCH IN DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

components = {'u', 'v', 'uu', 'vv', 'uv'};

clc;
for c = 1:length(components)
    component = components{c};

    % Monitor
    fprintf("%s\n", component)

    % Load component
    tmp = means.ensemble.(component);

    % Shift RHS
    LHS_values = tmp(:,86);
    RHS_values = tmp(:,87);
    differences = RHS_values - LHS_values;
    differences(isnan(differences)) = 0;
    corrected = tmp;
    corrected(:,87:end) = corrected(:,87:end) - differences;

    % Interpolate across strip
    [~, left] = min(abs(x - 115));
    [~, right] = min(abs(x - 120));
    corrected(:, left:right) = nan;
    corrected = fillmissing(corrected, 'makima', 2);
    corrected(Y < means.ensemble.max_wave_profile) = nan;

    % Resave data
    means_corrected.(component) = corrected;
    clear corrected

    fprintf("\n")
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING TO CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

levels = 10;

% Contours
for c = 1:length(components)
    component = components{c};
    figure()
    title(component)
    contourf(X, Y, means_corrected.(component), levels, 'linestyle', 'none')
    axis equal
end

% Profiles
for c = 1:length(components)
    component = components{c};
    
    figure()
    title(component)
    hold on
    for i = 1:10:length(x)
        plot(x, means.ensemble.(component)(i,:) + 0.01 * i, 'color', 'black')
        plot(x, means_corrected.(component)(i,:) + 0.01 * i, 'color', 'red')
    end
    hold off
    xlim([0, 235])
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTEGRATED BOUNDARY LAYER EQUATION (PHASE AVERAGE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Get velocities
u = means_corrected.u;
v = means_corrected.v;
uu = means_corrected.uu;
uv = means_corrected.uv;
vv = means_corrected.vv;

%% Viscous Shear Stress Term
[~, dudy] = gradient(u, dy);
term_1 = nu * dudy;
output.term_1 = term_1;

%% Turbulent Shear Stress
term_2 = -1 * uv;
output.term_2 = term_2;

%% Mean Momentum Flux
% Compute integrad
[du2dx, ~] = gradient(u.^2, dx);
% Create mask to remove nans
du2dx_mask   = isnan(du2dx);
du2dx_masked = du2dx;
du2dx_masked(du2dx_mask) = 0; 
% Compute integral
int_du2dx_dy = -1 * flipud(cumtrapz(flipud(y * 1E-3), flipud(du2dx_masked), 1));
int_du2dx_dy(du2dx_mask) = 0;

% Second half of term: as shown
[dudx, ~] = gradient(u, dx);
% Create mask to remove nans
dudx_mask   = isnan(dudx);
dudx_masked = dudx;
dudx_masked(dudx_mask) = 0;
% Compute integral
int_dudx_dy = flipud(cumtrapz(flipud(y * 1E-3), flipud(dudx_masked), 1));
int_dudx_dy(dudx_mask) = 0;

% Save term
term_3 = int_du2dx_dy + (int_dudx_dy .* u);
output.term_3 = term_3;


%% Turbulent Momentum Flux
% Compute integrad
[duudx, ~] = gradient(uu, dx);
% Create mask to remove nans
duudx_mask   = isnan(duudx);
duudx_masked = duudx;
duudx_masked(duudx_mask) = 0;
% Compute integrals
int_duudx_dy = -1 * flipud(cumtrapz(flipud(y * 1E-3), flipud(duudx_masked), 1));
int_duudx_dy(duudx_mask) = nan;
% Save
term_4 = int_duudx_dy;
output.term_4 = int_duudx_dy;


%% Pressure Correction
% Compute integrad
[dvvdx, ~] = gradient(vv, dx);
% Create mask to remove nans
dvvdx_mask   = isnan(dvvdx);
dvvdx_masked = dvvdx;
dvvdx_masked(dvvdx_mask) = 0;
% Compute integrals
int_dvvdx_dy = flipud(cumtrapz(flipud(y * 1E-3), flipud(dvvdx_masked), 1));  
int_dvvdx_dy(dvvdx_mask) = nan;
% Save
term_5 = int_dvvdx_dy;
output.term_5 = term_5;

%% Summation of all terms
output.sum = sum(cat(3, term_1, term_2, term_3, term_4, term_5), 3);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING INTEGRATED BOUNDARY LAYER EQUATION TERMS (CONTOURS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

levels = 50;
terms = fields(output);
figure()
tiledlayout(1,length(terms))

clc;
for t = 1:5
    term = terms{t};
    disp(term)
    nexttile
    hold on
    contourf(X,Y,output.(term), levels, 'linestyle', 'none')
    axis equal
    colorbar()
    title(term, 'interpreter', 'none')
    xlim([0, 235])
    ylim([-15, 200])
end
clc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING INTEGRATED BOUNDARY LAYER EQUATION TERMS (PROFILES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 86;
lw = 3;

terms = fields(output);
figure()
tiledlayout(1,4)

clc; close all;
hold on
for t = 1:5
    term = terms{t};
    disp(term)
    plot(output.(term)(:, idx), y, 'linewidth', lw, 'DisplayName', term)
    xlim([-0.2, 0.2])
    ylim([-15, 200])
end
plot(output.sum(:, idx), y, 'linewidth', lw, 'color', 'black', 'DisplayName', 'Sum')
hold off
legend('location', 'northwest', 'interpreter', 'none')

clc;

%% Plot term profiles on top of the contour

scale = 50;
component = 'uv';
lw = 2;

figure()
hold on

% Contour
contourf(X, Y, means_corrected.(component), 50, 'linestyle', 'none')
axis equal

% Profiles
for i = 1:20:length(x)
    for t = 1:length(terms)
        plot(scale * output.(terms{t})(:,i) + x(i), y, 'linewidth', lw)
    end
end
hold off
ylim([-20, 60])
xlim([0, 235])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING MEAN MOMENTUM FLUX CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
%%% Get velocities
u = means_corrected.u;
v = means_corrected.v;
uu = means_corrected.uu;
uv = means_corrected.uv;
vv = means_corrected.vv;

%%% Mean Momentum Flux
% Compute integrad
[du2dx, ~] = gradient(u.^2, dx);
% Create mask to remove nans
du2dx_mask   = isnan(du2dx);
du2dx_masked = du2dx;
du2dx_masked(du2dx_mask) = 0; 
% Compute integral
int_du2dx_dy = -1 * flipud(cumtrapz(flipud(y * 1E-3), flipud(du2dx_masked), 1));
int_du2dx_dy(du2dx_mask) = 0;

% Second half of term: as shown
[dudx, ~] = gradient(u, dx);
% Create mask to remove nans
dudx_mask   = isnan(dudx);
dudx_masked = dudx;
dudx_masked(dudx_mask) = 0;
% Compute integral
int_dudx_dy = flipud(cumtrapz(flipud(y * 1E-3), flipud(dudx_masked), 1));
int_dudx_dy(dudx_mask) = 0;

% With continuity substitution
term_3_written = int_du2dx_dy + (int_dudx_dy .* u);
term_3_continuity = int_du2dx_dy - (u .* v);

% Plot to compare
figure()
tiledlayout(1,3)
nexttile
contourf(X, Y, term_3_written, 100, 'linestyle', 'none')
axis equal
colorbar()
title("As written")

nexttile
contourf(X, Y, term_3_continuity, 100, 'linestyle', 'none')
axis equal
colorbar()
title("Continuity substitution")

nexttile
contourf(X, Y, term_3_written - term_3_continuity, 100, 'linestyle', 'none')
axis equal
colorbar()
title("Difference")


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTEGRATED BOUNDARY LAYER EQUATION (CONTINUTIY SUB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Get velocities
u = means_corrected.u;
v = means_corrected.v;
uu = means_corrected.uu;
uv = means_corrected.uv;
vv = means_corrected.vv;

%% Viscous Shear Stress Term
[~, dudy] = gradient(u, dy);
term_1 = nu * dudy;
output_sub.term_1 = term_1;

%% Turbulent Shear Stress
term_2 = -1 * uv;
output_sub.term_2 = term_2;

%% Mean Momentum Flux
% Compute integrad
[du2dx, ~] = gradient(u.^2, dx);
% Create mask to remove nans
du2dx_mask   = isnan(du2dx);
du2dx_masked = du2dx;
du2dx_masked(du2dx_mask) = 0; 
% Compute integral
int_du2dx_dy = -1 * flipud(cumtrapz(flipud(y * 1E-3), flipud(du2dx_masked), 1));
int_du2dx_dy(du2dx_mask) = 0;

% Save term
term_3 = int_du2dx_dy - (u .* v);
output_sub.term_3 = term_3;


%% Turbulent Momentum Flux
% Compute integrad
[duudx, ~] = gradient(uu, dx);
% Create mask to remove nans
duudx_mask   = isnan(duudx);
duudx_masked = duudx;
duudx_masked(duudx_mask) = 0;
% Compute integrals
int_duudx_dy = -1 * flipud(cumtrapz(flipud(y * 1E-3), flipud(duudx_masked), 1));
int_duudx_dy(duudx_mask) = nan;
% Save
term_4 = int_duudx_dy;
output_sub.term_4 = int_duudx_dy;


%% Pressure Correction
% Compute integrad
[dvvdx, ~] = gradient(vv, dx);
% Create mask to remove nans
dvvdx_mask   = isnan(dvvdx);
dvvdx_masked = dvvdx;
dvvdx_masked(dvvdx_mask) = 0;
% Compute integrals
int_dvvdx_dy = flipud(cumtrapz(flipud(y * 1E-3), flipud(dvvdx_masked), 1));  
int_dvvdx_dy(dvvdx_mask) = nan;
% Save
term_5 = int_dvvdx_dy;
output_sub.term_5 = term_5;
    

% Summation of all terms
output_sub.sum = sum(cat(3, term_1, term_2, term_3, term_4, term_5), 3);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING INTEGRATED BOUNDARY LAYER EQUATION TERMS (CONTOURS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


levels = 50;

terms = fields(output_sub);
figure()
tiledlayout(1,length(terms))

clc;
for t = 1:5
    term = terms{t};
    disp(term)
    nexttile
    hold on
    contourf(X,Y,output_sub.(term), levels, 'linestyle', 'none')
    axis equal
    colorbar()
    title(term, 'interpreter', 'none')
    xlim([0, 235])
    ylim([-15, 200])
end
clc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING INTEGRATED BOUNDARY LAYER EQUATION TERMS (PROFILES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 86;
lw = 3;

terms = fields(output_sub);
figure()
tiledlayout(1,4)

clc; close all;

hold on
for t = 1:5
    term = terms{t};
    disp(term)
    plot(output_sub.(term)(:, idx), y, 'linewidth', lw, 'DisplayName', term)
    xlim([-0.25, 0.25])
    ylim([-15, 200])
end
plot(output_sub.sum(:, idx), y, 'linewidth', lw, 'color', 'black', 'DisplayName', 'Sum')
hold off
legend('location', 'northeast', 'interpreter', 'none')
clc;



