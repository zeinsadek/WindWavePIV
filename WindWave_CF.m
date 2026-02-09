%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');


%% Import Data

% Cases with waves
cases = {'0', 'A', 'B', 'C', 'D'};
waves = {'A', 'B', 'C', 'D'};
names = {'none', 'A', 'B', 'C', 'D'};

for i = 1:length(cases)
    
    c          = strcat('WT4_WV', cases{i}, '_AG0');
    LL_path    = strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law/', c,'_LL.mat');
    BL_path    = strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer/', c,'_BL.mat');
    means_path = strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/', c,'_MEANS.mat');

    % BL Stuff
    BL_temp    = load(BL_path);
    BL_temp    = BL_temp.output;
    
    BL_data.(names{i}) = BL_temp;
    
    % LL Stuff
    LL_temp    = load(LL_path);
    LL_temp    = LL_temp.out;
    
    LL_data.(names{i}) = LL_temp;
    
    % Means Stuff
    means_temp = load(means_path);
    means_temp = means_temp.output;
    means.(names{i}) = means_temp;
end

%% Brzek 2007

wave          = 'D';
wave_length   = wave_parameters{find(strcmp(wave_parameters, wave) == 1), 2};
phase_offsets = [0, wave_length/4, wave_length/2, 3*wave_length/4];

% cos_fit       = @(b, v) wave_amplitude * cos(2 * pi * (v - (max(x)/2) - b(1)) / wave_length);



for p = 1:1
    
    deriv_cos_fit = @(b, v) -1 * b(2) * sin(2 * pi * (v - (max(x)/2) - b(1)) / wave_length);
    cos_fit       = @(b, v) b(2) * cos(2 * pi * (v - (max(x)/2) - b(1)) / wave_length);
    
%     figure()
%     hold on
%     plot(x, cos_fit(phase_offsets(p), x))
%     plot(x, deriv_cos_fit(phase_offsets(p), x))
%     hold off
    
    delta      = BL_data.(wave).phase(p).thickness;
    delta_star = BL_data.(wave).phase(p).displacement;
    theta      = BL_data.(wave).phase(p).momentum;
    H          = delta_star./ theta;
    u_inf      = BL_data.(wave).u_inf;
    
    
%     figure()
%     hold on
%     plot(x, theta)
%     plot(x, cos_fit(phase_offsets(p), x))
%     hold off

    mean_U     = means.(wave).ensemble.u;
    mean_U_row = mean_U(20,:);
    mean_U_row = mean_U_row.';
    
    x  = unique(means.D.X);
    dx = x(3) - x(2);
    
    
    
    % %%% Repair middle of image
    OG        = means.(wave).phase(p).u;
    data_mask = ~isnan(OG);
    [horz_boarders, vert_boarders]  = gradient(data_mask);
    [~, vert_locs]                  = findpeaks(abs(vert_boarders(:, 87)));
    
    means.(wave).phase(p).u(:, 86 - 1:87 + 1) = nan;
    
    for i = vert_locs(1):vert_locs(2)
        means.(wave).phase(p).u(i, 3:169) = fillmissing(means.(wave).phase(p).u(i,3:169), 'linear');
    end


    %%% Von Karaman Equation
    % dtheta / dx term
    dthetadx = diff(sgolayfilt(theta(2:end-1), 3, 31),1,1)/dx;
    sgolay_dthetadx = sgolayfilt(dthetadx, 3, 31);

    % Other term
    du_infdx   = diff(mean_U_row, 1, 1) / dx;
    other_term = (2 + H(1:end-1)) .* (theta(1:end-1) ./ mean_U_row(1:end-1)) .* du_infdx;
    other_term = other_term(1:end-2);
    

    
    figure()
    hold on
    plot(x(2:end-2), sgolay_dthetadx, 'LineWidth', 2, 'DisplayName', '$d \theta / dx$')
    plot(x, -1 * deriv_cos_fit([phase_offsets(p), 0.075], x), 'LineWidth', 2, 'DisplayName', '$\sin \left(\frac{2 \pi}{\lambda}(x - \phi) \right)$')
    plot(x, cos_fit([phase_offsets(p), 0.02], x), 'LineWidth', 2, 'Color', 'red', 'DisplayName', 'Wave')
    hold off
    xlim([0, max(x) - 5])
    xlabel('x [mm]')
    ylabel('$d \theta / dx$', 'Interpreter', 'Latex')
    legend('Interpreter', 'Latex', 'FontSize', 18)
   
    
   
%     % Plot
%     figure()
%     tiledlayout(2,1);
% 
%     nexttile
%     hold on
%     plot(x(2:end-2), sgolay_dthetadx, 'DisplayName', '$d \theta / dx$')
% %     plot(x(2:end-2), other_term, 'DisplayName', '$(2 + H) \frac{\theta}{u_{\infty}} \frac{d u_{\infty}}{dx}$');
% %     plot(x(2:end-2), sgolay_dthetadx + sgolayfilt(double(other_term), 3, 21), 'DisplayName', 'Sum')
%     xline(LL_data.(wave).profile_x_value)
%     
%     plot(x, -1 * deriv_cos_fit(phase_offsets(p), x))
%     
%     
%     hold off
%     legend('Interpreter', 'latex', 'FontSize', 16)
%     xlim([0, max(x)])
% 
% 
%     % Skin friction & friction velocity
%     nexttile
    figure()
    cf               = 2 * (sgolay_dthetadx);
    vk_u_star        = mean_U_row(2:end-2) .* sqrt(0.5 * cf);
    sgolay_vk_u_star = sgolayfilt(double(vk_u_star), 3, 11);

    hold on
    plot(x(2:end-2), vk_u_star, 'DisplayName', '$u^* = u_{\infty} \sqrt{\frac{1}{2} C_f}$')
    plot(x, sqrt(-1 * deriv_cos_fit([phase_offsets(p), 0.5], x)), 'DisplayName', '$u^* = \sqrt{\sin \left(\frac{2 \pi}{\lambda}(x - \phi) \right)}$')
    plot(x, cos_fit([phase_offsets(p), 0.075], x), 'LineWidth', 2, 'Color', 'red', 'DisplayName', 'Wave')
    hold off
    
    xlim([0, max(x) - 5])
    xlabel('x [mm]')
    ylabel('$u^*$', 'Interpreter', 'Latex')
    legend('Interpreter', 'Latex', 'FontSize', 18)
    
    
%     
%     
%     xline(LL_data.(wave).profile_x_value)
%     hold off
%     xlim([0, max(x)])
% 
%     % Compare u* values
%     fprintf('IBL u*: %3.4f m/s\n', LL_data.(wave).phase(p).u_star);
%     fprintf('VK u*: %3.4f m/s\n\n', sgolay_vk_u_star(LL_data.(wave).profile_idx));

end





%% Compare IBL and VK Equations, Ensemble Average


for i = 1:length(names)
    
    name = names{i};
    
    % BL Thicknesses
    delta      = BL_data.(name).ensemble.thickness;
    delta_star = BL_data.(name).ensemble.displacement;
    theta      = BL_data.(name).ensemble.momentum;
    H          = delta_star./ theta;
    u_inf      = BL_data.(name).u_inf;

    % Freestream
    mean_U     = means.(name).ensemble.u;
    mean_U_row = mean_U(20,:);
    mean_U_row = mean_U_row.';

    x  = unique(means.D.X);
    dx = x(3) - x(2);

    % %%% Repair middle of image
    % OG        = means.none.ensemble.u;
    % data_mask = ~isnan(OG);
    % [horz_boarders, vert_boarders]  = gradient(data_mask);
    % [~, vert_locs]                  = findpeaks(abs(vert_boarders(:, 87)));
    % 
    % means.none.ensemble.u(:, 86 - 1:87 + 1) = nan;
    % 
    % for i = vert_locs(1):vert_locs(2)
    %     means.none.ensemble.u(i, 3:169) = fillmissing(means.none.ensemble.u(i,3:169), 'linear');
    % end



    %%% Von Karaman Equation
    % dtheta / dx term
    sgolay_theta    = sgolayfilt(theta(2:end-1), 3, 31);
    dthetadx        = diff(sgolay_theta, 1, 1)/dx;
    sgolay_dthetadx = sgolayfilt(dthetadx, 3, 31);

    % Other term
    du_infdx          = diff(mean_U_row, 1, 1) / dx;
    other_term        = (2 + H(1:end-1)) .* (theta(1:end-1) ./ mean_U_row(1:end-1)) .* du_infdx;
    other_term        = other_term(1:end-2);
    sgolay_other_term = sgolayfilt(double(other_term), 3, 21);


    % Plot
    figure(i)
    t = tiledlayout(3,1);
    title(t, name)
    
    nexttile
    hold on
%     plot(x, delta)
    plot(x, theta)
    plot(x(2:end-1), sgolay_theta)
    hold off
    ylabel('$\theta$', 'Interpreter', 'latex')
    xlim([0, max(x)-2])
    title('Momentum Thickness')
    

    nexttile
    hold on
    plot(x(2:end-2), sgolay_dthetadx, 'DisplayName', '$d \theta / dx$')
    plot(x(2:end-2), other_term, 'DisplayName', '$(2 + H) \frac{\theta}{u_{\infty}} \frac{d u_{\infty}}{dx}$');
%     plot(x(2:end-2), sgolay_other_term)
    plot(x(2:end-2), dthetadx + other_term, 'DisplayName', 'Sum')
    xline(LL_data.D.profile_x_value)
    hold off
    legend('Interpreter', 'latex', 'FontSize', 16)
    xlim([0, max(x)-2])
    title('Terms for Friction Factor')


    % Skin friction & friction velocity
    cf        = 2 * (sgolay_dthetadx);% + other_term);
    vk_u_star = mean_U_row(2:end-2) .* sqrt(0.5 * cf);
    sgolay_vk_u_star = sgolayfilt(double(vk_u_star), 3, 11);

    nexttile
    hold on
    plot(x(2:end-2), vk_u_star)
%     plot(x(2:end-2), sgolay_vk_u_star)
    xline(LL_data.D.profile_x_value)
    hold off
    xlim([0, max(x)-2])
    title('Friction Velocity')
    xlabel('x [mm]')
    
    
    % Print u* values
    fprintf('%s\nIBL u*: %3.4f m/s\n', name, LL_data.(name).ensemble.u_star);
    fprintf('VK u*: %3.4f m/s\n\n', sgolay_vk_u_star(LL_data.none.profile_idx));

end












%%


% figure()
% scatter(x, means.D.ensemble.u(10,:))
% xline(x(86))
% xline(x(87))

%%% Crop boarder
OG = means.D.ensemble.u;
% Find region of good data
data_mask = ~isnan(OG);
[horz_boarders, vert_boarders]  = gradient(data_mask);

%%% Determine size of image height
[pks, vert_locs] = findpeaks(abs(vert_boarders(:, 87)));
[pks, horz_locs] = findpeaks(abs(horz_boarders(50, :)));

%%% Repair middle section
test = OG;
test(:, 86 - 1:87 + 1) = nan;

% Crop below water
OG(means.D.Y < 0) = nan;
test(means.D.Y < 0) = nan;

% figure()
% contourf(means.D.X, means.D.Y, vert_boarders, 500, 'LineStyle', 'none')
% colorbar()
% 
% figure()
% plot(abs(vert_boarders(:,50)))



for i = vert_locs(1):vert_locs(2)
    test(i, 3:169) = fillmissing(test(i,3:169), 'linear');
end

% figure()
% contourf(means.D.X, means.D.Y, test, 100, 'LineStyle', 'none')
% title('Cleaned')
% 
% figure()
% contourf(means.D.X, means.D.Y, OG, 100, 'LineStyle', 'none')
% title('Original')


figure()
hold on
for i = 10:10
%     plot(x, OG(i,:), 'red')
    plot(x, test(i,:), 'black')
    
    du_infdx = (test(i,end-1) - test(i,2)) / ((x(end-1) - x(2)) * 10^-3);
    
    plot(x, test(i,2) + (x .* 10^-3 .* du_infdx))
    
end
hold off
xlim([0, max(x)])
xlabel('x [mm]')
ylabel('u [m/s]')

























