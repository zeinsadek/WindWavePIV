%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

%% Import Log Law Values

LL_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/Hopkins/Initial_Results/log_law/';
cases   = {'WT4_WVD_AG0', 'WT6_WVD_AG0', 'WT8_WVD_AG0'};

for i = 1:length(cases)
    experiment      = cases{i};
    temp            = load(strcat(LL_path, experiment, '_shear_log_law.mat'));
    data.(experiment) = temp.out;
end


%% Ensemble BL Profiles

figure('Name', 'Norm Ens BL Profiles')
sgtitle('Ensemble Average Velocity Profiles')
hold on
for i = 1:length(cases)
    experiment = cases{i};
    
    profile = data.(experiment).ensemble.profile;
    u_inf   = data.(experiment).u_inf;
    u_star  = data.(experiment).ensemble.u_star;
    y       = data.(experiment).y;
    
    subplot(1,3,1)
    hold on
    plot(profile, y, 'DisplayName', experiment, 'LineWidth', 2)
    ylabel('y [mm]', 'Interpreter', 'Latex')
    xlabel('$u$', 'Interpreter', 'latex')
    ylim([-10, 210])
    title('Non-Normalized', 'Interpreter', 'latex')
    
    subplot(1,3,2)
    hold on
    plot(profile/u_inf, y, 'DisplayName', experiment, 'LineWidth', 2)
    xlabel('$\frac{u}{u_{\infty}}$', 'Interpreter', 'latex', 'FontSize', 18)
    ylim([-10, 210])
    title('Normalized by $$u_{\infty}$$', 'Interpreter', 'latex')
    
    subplot(1,3,3)
    hold on
    plot(profile/u_star, y, 'DisplayName', experiment, 'LineWidth', 2)
    xlabel('$\frac{u}{u^*}$', 'Interpreter', 'latex', 'FontSize', 18)
    ylim([-10, 210])
    title('Normalized by $$u^*$$', 'Interpreter', 'latex')
end
hold off
hold off
hold off
legend(cases, 'Interpreter', 'none', 'Location', 'NorthWest', 'Orientation', 'Vertical');
set(gcf, 'Units', 'Inches', 'Position', [1, 5, 15, 5])

%% Phase Average Profiles: Inflow Sorting

figure('Name', 'Norm Ens BL Profiles')

hold on
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact')
sgtitle('Phase Average Velocity Profiles: Sorted by Inflow')

for i = 1:length(cases)

    experiment  = cases{i};
    ens_profile = data.(experiment).ensemble.profile;
    u_inf       = data.(experiment).u_inf;
    ens_u_star  = data.(experiment).ensemble.u_star;
    y           = data.(experiment).y;
    
    ax(i) = nexttile;
    title(experiment, 'Interpreter', 'none')
    hold on
    for j = 1:4
        profile = data.(experiment).phase(j).profile;
        plot(profile, y, 'DisplayName', strcat('Phase', num2str(j)), 'LineWidth', 2)
        xlabel('$u [m/s]$', 'Interpreter', 'latex')
        ylim([-10, 60])
        
        if j == 1
            ylabel('y [mm]', 'Interpreter', 'Latex')
        end

    end
    p = plot(ens_profile, y, 'Color', 'black', 'DisplayName', 'Ensemble', 'LineWidth', 2);
    p.Color(4) = 0.5;
    hold off
end
hold off
linkaxes(ax, 'xy')
lh = legend({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Ensemble'}, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
lh.Layout.Tile = 'North';
set(gcf, 'Units', 'Inches', 'Position', [1, 5, 15, 5])

%% Phase Average Profiles: Inflow Sorting & Subtracting Ensemble

figure('Name', 'Norm Ens BL Profiles')

hold on
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact')
sgtitle('Phase Average  - Ensemble Velocity Profiles: Sorted by Inflow')

for i = 1:length(cases)

    experiment  = cases{i};
    ens_profile = data.(experiment).ensemble.profile;
    u_inf       = data.(experiment).u_inf;
    ens_u_star  = data.(experiment).ensemble.u_star;
    y           = data.(experiment).y;
    
    ax(i) = nexttile;
    title(experiment, 'Interpreter', 'none')
    hold on
    for j = 1:4
        profile = data.(experiment).phase(j).profile;
        phase_minus_ensemble = profile - ens_profile;
 
        plot(phase_minus_ensemble, y, 'DisplayName', strcat('Phase', num2str(j)), 'LineWidth', 2)
        xlabel('$u [m/s]$', 'Interpreter', 'latex')
        ylim([-10, 210])
        
        if j == 1
            ylabel('y [mm]', 'Interpreter', 'Latex')
        end

    end
    hold off
end
hold off
linkaxes(ax, 'xy')
lh = legend({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4'}, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
lh.Layout.Tile = 'North';
set(gcf, 'Units', 'Inches', 'Position', [1, 5, 15, 5])



%% Phase Average BL Profiles: Phase Sorting

scaling_names = {'Nothing', 'Freestream', 'Ensemble u*', 'Phase u*'};

for s = 1:4
    
    figure('Name', strcat('Phase Avg Profiles:', scaling_names{s}))
    hold on
    tiledlayout(1,4, 'Padding', 'compact', 'TileSpacing', 'compact')
    sgtitle(strcat('Phase Average Profiles: Normalized by', {' '}, scaling_names{s}))
    
    for j = 1:4 
        ax(j) = nexttile;
        for i = 1:length(cases)
            
            experiment  = cases{i};
            ens_profile = data.(experiment).ensemble.profile;
            u_inf       = data.(experiment).u_inf;
            ens_u_star  = data.(experiment).ensemble.u_star;
            y           = data.(experiment).y;

            profile = data.(experiment).phase(j).profile;
            phase_u_star = data.(experiment).phase(j).u_star;

            % Normalization
            if s == 1
                norm = 1;
                xlab = '%u%';
                
            elseif s == 2
                norm = u_inf;
                xlab = '$u / u_{\infty}$';
                
            elseif s == 3
                norm = ens_u_star;
                xlab = '$u / u^*$';
                
            elseif s == 4
                norm = phase_u_star;
                xlab = '$u / {u_{\phi}}^*$';
            end
                 
            hold on
            plot(profile / norm, y, 'DisplayName', experiment, 'LineWidth', 2)
            ylim([-10, 210])
            xlabel(xlab, 'Interpreter', 'latex', 'FontSize', 14)

            if j == 1
                ylabel('y [mm]', 'Interpreter', 'Latex', 'FontSize', 14)
            end
            title(strcat('Phase', {' '},  num2str(j)))
        end
        hold off
    end

    hold off
    linkaxes(ax, 'xy')
    lh = legend(cases, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
    lh.Layout.Tile = 'North';
    set(gcf, 'Units', 'Inches', 'Position', [1, 5, 15, 5])
end


% figure('Name', 'Norm Ens BL Profiles')
% hold on
% tiledlayout(1,4, 'Padding', 'compact', 'TileSpacing', 'compact')
% sgtitle('Phase Average Velocity Profiles: Sorted by Phase')
% for j = 1:4 
%     ax(j) = nexttile;
%     for i = 1:length(cases)
%         experiment  = cases{i};
%         ens_profile = data.(experiment).ensemble.profile;
%         u_inf       = data.(experiment).u_inf;
%         ens_u_star  = data.(experiment).ensemble.u_star;
%         y           = data.(experiment).y;
%     
%         profile = data.(experiment).phase(j).profile;
%         phase_u_star = data.(experiment).phase(j).u_star;
%         
%         % normalizing factor
%         norm = u_inf;
% 
%         hold on
%         plot(profile / norm, y, 'DisplayName', experiment, 'LineWidth', 2)
%         xlabel('$\frac{u}{u_{\infty}}$', 'Interpreter', 'latex')
%         ylim([-10, 210])
%         
%         if j == 1
%             ylabel('y [mm]', 'Interpreter', 'Latex')
%         end
%         title(strcat('Phase', {' '},  num2str(j)))
%     end
%     hold off
%   
% end
% 
% hold off
% linkaxes(ax, 'xy')
% lh = legend(cases, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
% lh.Layout.Tile = 'North';
% set(gcf, 'Units', 'Inches', 'Position', [1, 5, 15, 5])

%% All Curves Plotted Togeher

colors = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250]];
scaling_names = {'Nothing', 'Freestream', 'Ensemble u*', 'Phase u*'};

for s = 1:4
    h = zeros(1,3);
    figure('Name', strcat('All Profiles:', scaling_names{s}))
    hold on
    sgtitle(strcat('All Profiles: Normalized by', {' '}, scaling_names{s}))
    for i = 1:length(cases)

        experiment  = cases{i};
        ens_profile = data.(experiment).ensemble.profile;
        u_inf       = data.(experiment).u_inf;
        ens_u_star  = data.(experiment).ensemble.u_star;
        y           = data.(experiment).y;

        for j = 1:4
            profile      = data.(experiment).phase(j).profile;
            phase_u_star = data.(experiment).phase(j).u_star;

            % Normalization
            if s == 1
                norm = 1;
                xlab = '%u%';
                
            elseif s == 2
                norm = u_inf;
                xlab = '$u / u_{\infty}$';
                
            elseif s == 3
                norm = ens_u_star;
                xlab = '$u / u^*$';
                
            elseif s == 4
                norm = phase_u_star;
                xlab = '$u / {u_{\phi}}^*$';
            end
            
            temp = plot(profile / norm, y, 'Color', colors(i,:), 'DisplayName', experiment, 'LineStyle', '--', 'LineWidth', 2);
            temp.Color(4) = 0.5;

            xlabel('$u [m/s]$', 'Interpreter', 'latex')
            ylabel('y [mm]', 'Interpreter', 'Latex')
            ylim([-10, 210])
        end
        
        if s ~= 4
            h(i) = plot(ens_profile / norm, y, 'Color', colors(i,:), 'LineStyle', '-', 'DisplayName', experiment, 'LineWidth', 2);
        else 
            h(i) = plot(ens_profile / ens_u_star, y, 'Color', colors(i,:), 'LineStyle', '-', 'DisplayName', experiment, 'LineWidth', 2);
        end
        
    end
    hold off
    legend(h, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
    set(gcf, 'Units', 'Inches', 'Position', [1, 5, 8, 7])
end



% %% Ensemble Viscous Term
% 
% figure('Name', 'Ens nududy')
% sgtitle('Ensemble Average, Fixed Location Viscous Shear Profiles')
% hold on
% for i = 1:length(cases)
%     experiment = cases{i};
%     
%     nu_dudy = data.(experiment).ensemble.nu_dudy;
%     u_inf   = data.(experiment).u_inf;
%     y       = data.(experiment).y;
%     
%     subplot(1,2,1)
%     hold on
%     plot(nu_dudy, y(1:end-1), 'DisplayName', experiment, 'LineWidth', 2)
%     ylabel('y [mm]', 'Interpreter', 'Latex')
%     xlabel('$\nu \frac{\partial u}{\partial y}$', 'Interpreter', 'latex')
%     ylim([-10, 210])
% %     yline(0, 'Color', 'red', 'LineWidth', 2)
%     
%     subplot(1,2,2)
%     hold on
%     plot(nu_dudy/(u_inf^2), y(1:end-1), 'DisplayName', experiment, 'LineWidth', 2)
%     xlabel('$\frac{\nu \frac{\partial u}{\partial y}}{{u_{\infty}}^2}$', 'Interpreter', 'latex', 'FontSize', 18)
%     ylim([-10, 210])
% %     yline(0, 'Color', 'red', 'LineWidth', 2)
% end
% hold off
% hold off
% hold off
% legend('Interpreter', 'none', 'Location', 'NorthEast')

% %% Phase Viscous Term
% 
% figure('Name', 'Phase nududy')
% sgtitle('Phase Average, Fixed Location Viscous Shear')
% hold on
% for i = 1:length(cases)
%     experiment = cases{i};
%     
%     for j = 1:4
%         nu_dudy = data.(experiment).phase(j).nu_dudy;
%         u_inf   = data.(experiment).u_inf;
%         y       = data.(experiment).y;
% 
%         subplot(1,4,j)
%         hold on
%         plot(nu_dudy, y(1:end-1), 'DisplayName', experiment, 'LineWidth', 2)
%         xlabel('$\nu \frac{\partial u}{\partial y}$', 'Interpreter', 'latex')
%         %xlabel('$\frac{\nu \frac{\partial u}{\partial y}}{{u_{\infty}}^2}$', 'Interpreter', 'latex', 'FontSize', 18)
%         ylim([-10, 210])
%         xlim([0, 15E-4])
%         title(strcat('Phase ', num2str(j)))
% 
%         if j == 1
%             ylabel('y [mm]', 'Interpreter', 'Latex')
%         end
% 
%     end
%     
%     
% end
% hold off
% hold off
% hold off
% legend('Interpreter', 'none', 'Location', 'NorthEast')
% 

%% Ensemble Reynolds Stress

figure('Name', 'Ens uv')
sgtitle('Ensemble Average, Fixed Location Reynolds Shear Stress')
hold on
for i = 1:length(cases)
    experiment = cases{i};
    
    uv      = data.(experiment).ensemble.uv;
    u_inf   = data.(experiment).u_inf;
    ens_u_star  = data.(experiment).ensemble.u_star;
    y       = data.(experiment).y;
    
    subplot(1,3,1)
    hold on
    plot(uv, y, 'DisplayName', experiment, 'LineWidth', 2)
    ylabel('y [mm]', 'Interpreter' ,'Latex')
    xlabel('$\overline{uv}$', 'Interpreter', 'latex', 'FontSize', 18)
    ylim([-10, 210])
    
    subplot(1,3,2)
    hold on
    plot(uv/(u_inf^2), y, 'DisplayName', experiment, 'LineWidth', 2)
    xlabel('$\overline{uv} / u_{\infty}^2$', 'Interpreter', 'latex', 'FontSize', 18)
    ylim([-10, 210])
    
    subplot(1,3,3)
    hold on
    plot(uv/(ens_u_star^2), y, 'DisplayName', experiment, 'LineWidth', 2)
    xlabel('$\overline{uv} / \left(u^*\right)^2$', 'Interpreter', 'latex', 'FontSize', 18)
    ylim([-10, 210])
    

end
hold off
hold off
hold off
legend('Interpreter', 'none', 'Location', 'NorthEast')
set(gcf, 'Units', 'Inches', 'Position', [1, 5, 15, 5])

%% Ensemble & Phase Reynolds Stress

figure('Name', strcat('Phase Avg Profiles:', scaling_names{s}))
hold on
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle('Ensemble Average, Fixed Location Reynolds Shear Stress')
   
for i = 1:length(cases)
    
    h = zeros(1,3);
    experiment = cases{i};
    uv         = data.(experiment).ensemble.uv;
    u_inf      = data.(experiment).u_inf;
    ens_u_star = data.(experiment).ensemble.u_star;
    y          = data.(experiment).y;
    
    ax(i) = nexttile;
    hold on
    h(1) = plot(uv, y, 'DisplayName', 'Ensemble', 'LineWidth', 2);
    
    for j = 1:4
        phase_uv = data.(experiment).phase(j).uv;
        label = strcat('Phase', {' '}, num2str(j));
        h(j + 1) = plot(phase_uv, y, 'LineWidth', 2, 'DisplayName', label{1});
    end
    
    ylabel('y [mm]', 'Interpreter' ,'Latex')
    xlabel('$\overline{uv}$', 'Interpreter', 'latex', 'FontSize', 18)
    ylim([-10, 210])
    
    hold off

end
linkaxes(ax, 'xy')
lh = legend(h, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
lh.Layout.Tile = 'North';
set(gcf, 'Units', 'Inches', 'Position', [1, 5, 15, 5])

%% Phase Reynolds Stress

figure('Name', 'Phase uv')
sgtitle('Phase Average, Fixed Location Normalized Reynolds Shear Stress')
hold on
for i = 1:length(cases)
    experiment = cases{i};
    
    for j = 1:4
        uv = data.(experiment).phase(j).uv;
        u_inf   = data.(experiment).u_inf;
        y       = data.(experiment).y;

        subplot(1,4,j)
        hold on
        plot(uv / (u_inf^2), y, 'DisplayName', experiment, 'LineWidth', 2)
        xlabel('$\frac{uv}{{u_{\infty}}^2}$', 'Interpreter', 'latex')
        ylim([-10, 210])
        title(strcat('Phase ', num2str(j)))
        xlim([0, 5E-3])
%         yline(0, 'Color', 'red', 'LineWidth', 2, 'DisplayName', 'water')

        if j == 1
            ylabel('y [mm]', 'Interpreter', 'Latex')
        end

    end
end
hold off
hold off
hold off
legend('Interpreter', 'none', 'Location', 'NorthEast')

%% Ensemble u* Profile

figure('Name', 'Ens u*')
sgtitle('Ensemble Average, Fixed Location u* Profile')
hold on
for i = 1:length(cases)
    experiment = cases{i};
    
    u_star_profile = data.(experiment).ensemble.u_star_profile;
    y       = data.(experiment).y;
    
    plot(u_star_profile, y(1:end-1), 'DisplayName', experiment, 'LineWidth', 2)
    ylabel('y [mm]', 'Interpreter', 'Latex')
    xlabel('$u^*$', 'Interpreter', 'latex')
    ylim([-10, 210])
  
end

hold off
legend('Interpreter', 'none', 'Location', 'NorthEast')

%% Phase u* Profile

figure('Name', 'Phase u*')
sgtitle('Phase Average, Fixed Location u* Profile')
hold on
for i = 1:length(cases)
    experiment = cases{i};
    
    for j = 1:4
        u_star_profile = data.(experiment).phase(j).u_star_profile;
        y       = data.(experiment).y;

        subplot(1,4,j)
        hold on
        plot(u_star_profile, y(1:end-1), 'DisplayName', experiment, 'LineWidth', 2)
        xlabel('$u^*$', 'Interpreter', 'latex')
        ylim([-10, 210])
        xlim([0, 0.4])
        
        title(strcat('Phase ', num2str(j)))
        
        if j == 1
            ylabel('y [mm]', 'Interpreter', 'Latex')
        end

    end
end
hold off
hold off
hold off
legend('Interpreter', 'none', 'Location', 'NorthEast')

%% Ensemble Log Law

figure('Name', 'Ens LL')
sgtitle('Ensemble Average, Fixed Location Log Law')
hold on
for i = 1:length(cases)
    experiment = cases{i};
    
    u_plus = data.(experiment).ensemble.u_plus;
    y_plus = data.(experiment).ensemble.y_plus;
    
    scatter(y_plus, u_plus, 'filled', 'DisplayName', experiment)
  
end
hold off
grid on
set(gca,'xscale','log')
ylabel('$u^+$', 'Interpreter', 'latex')
xlabel('$y^+$', 'Interpreter', 'latex')
legend('Interpreter', 'none', 'Location', 'NorthWest')


%% Phase Log Law: Phase Sorted

figure('Name', 'Norm Ens BL Profiles')

hold on
tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'compact')
sgtitle('Phase Average  - Ensemble Velocity Profiles: Sorted by Inflow')


for j = 1:4
    h = zeros(1,3);
    ax(j) = nexttile;
    title(strcat('Phase', {' '}, num2str(j)))
    hold on
    for i = 1:length(cases)
        
        experiment = cases{i};
        
        u_plus = data.(experiment).phase(j).u_plus;
        y_plus = data.(experiment).phase(j).y_plus;
        
        h(i) = scatter(y_plus, u_plus, 'filled', 'DisplayName', experiment);

    end
    hold off
    grid on
    set(gca,'xscale','log')
    xlabel('$y^+$', 'Interpreter', 'latex')
    
    if j == 1
        ylabel('$u^+$', 'Interpreter', 'latex')
    end
    
end
linkaxes(ax, 'xy')
set(gcf, 'Units', 'Inches', 'Position', [1, 5, 12, 6])
lh = legend(h, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
lh.Layout.Tile = 'North';

%% Phase Log Law: Inflow Sorted

figure('Name', 'Norm Ens BL Profiles')

hold on
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact')
sgtitle('Phase Average  - Ensemble Velocity Profiles: Sorted by Inflow')


for i = 1:3
    
    h = zeros(1,5);
    
    experiment  = cases{i};
    ens_u_plus = data.(experiment).ensemble.u_plus;
    ens_y_plus = data.(experiment).ensemble.y_plus;
        
    ax(i) = nexttile;
    title(strcat(experiment, ' Log Law'), 'Interpreter', 'none')
    hold on
    
    h(1) = scatter(ens_y_plus, ens_u_plus, 'filled', 'DisplayName', 'Ensemble');
    
    for j = 1:4

        u_plus = data.(experiment).phase(j).u_plus;
        y_plus = data.(experiment).phase(j).y_plus;
        
        label = strcat('Phase', {' '}, num2str(j));
        h(j + 1) = scatter(y_plus, u_plus, 'filled', 'DisplayName', label{1});

    end
    hold off
    grid on
    set(gca,'xscale','log')
    xlabel('$y^+$', 'Interpreter', 'latex')
    
    if i == 1
       ylabel('$u^+$', 'Interpreter', 'latex') 
    end
    
end

linkaxes(ax, 'xy')
set(gcf, 'Units', 'Inches', 'Position', [1, 5, 12, 6])
lh = legend(h, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
lh.Layout.Tile = 'North';



