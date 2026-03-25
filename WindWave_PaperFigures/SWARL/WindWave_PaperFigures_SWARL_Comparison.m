clear; clc;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/SWARL')
%% ==================== WAVE PARAMETERS FROM TABLE ====================
wave = struct();

wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};

wave(1).name = '$\lambda_{410}$, $ak_{0.180}$';
wave(1).lambda = 410.8e-3;   wave(1).a = 11.8e-3;
wave(1).c = 0.81;            wave(1).ak = 0.180;
% wave(1).color = [0.9 0.1 0.1];  % Red
wave(1).color = '#FE6202';

wave(2).name = '$\lambda_{313}$, $ak_{0.211}$';
wave(2).lambda = 313.3e-3;   wave(2).a = 10.5e-3;
wave(2).c = 0.71;            wave(2).ak = 0.211;
% wave(2).color = [1.0 0.75 0.0]; % Gold/Yellow
wave(2).color = '#DC2680';

wave(3).name = '$\lambda_{189}$, $ak_{0.305}$';
wave(3).lambda = 189.6e-3;   wave(3).a = 9.2e-3;
wave(3).c = 0.57;            wave(3).ak = 0.305;
% wave(3).color = [0.0 0.6 0.2]; % Green
wave(3).color = '#775EEF';

wave(4).name = '$\lambda_{124}$, $ak_{0.267}$';
wave(4).lambda = 124.3e-3;   wave(4).a = 5.3e-3;
wave(4).c = 0.49;            wave(4).ak = 0.267;
% wave(4).color = [0.0 0.8 0.9]; % Cyan
wave(4).color = '#648FFF';

% Velocity markers
u_inf = [2.42, 3.87, 5.43];  % m/s
markers = {'^', 's', 'o'};   % Triangle, Square, Circle
marker_size = 120;

%% ==================== DIGITIZED DATA ====================
raw_data = [
    1.2707774798927614, 10.638394670566253;
    1.463806970509383,  11.06409943943456;
    1.6782841823056298, 10.641481842554228;
    1.9892761394101877, 10.37111056950199;
    2.0,                9.219676659354944;
    2.6005361930294906, 6.95149890324153;
    3.083109919571046,  5.955154764806238;
    3.1260053619302943, 6.864570639369568;
    3.4798927613941015, 8.079372816638232;
    3.812332439678284,  5.233406450564626;
    5.109919571045577,  4.818994231862863;
    6.514745308310992,  2.9508489722966935;
];

assignments = [
%   case  vel
%   case_id: 1=Red(λ410), 2=Yellow(λ313), 3=Green(λ189), 4=Cyan(λ124)
%   vel_id:  1=Triangle(2.42m/s), 2=Square(3.87m/s), 3=Circle(5.43m/s)
    4,    3;    % Row 1:  Cyan circle
    3,    3;    % Row 2:  Green circle
    4,    2;    % Row 3:  Cyan square
    3,    2;    % Row 4:  Green square
    2,    3;    % Row 5:  Yellow circle
    1,    3;    % Row 6:  Red circle
    2,    2;    % Row 7:  Yellow square
    4,    1;    % Row 8:  Cyan triangle
    3,    1;    % Row 9:  Green triangle
    1,    2;    % Row 10: Red square
    2,    1;    % Row 11: Yellow triangle
    1,    1;    % Row 12: Red triangle
];

c_ustar = raw_data(:,1);
delta_U = raw_data(:,2);
case_id = assignments(:,1);
vel_id = assignments(:,2);

%% ==================== PHYSICAL CONSTANTS ====================
nu = 1.46e-5;    % kinematic viscosity [m^2/s]
kappa = 0.4;    % von Karman constant

%% ==================== CALCULATE MODEL PREDICTIONS ====================
delta_U_model = zeros(size(c_ustar));

for i = 1:length(c_ustar)
    % Get wave parameters for this case
    Hs_ref = wave(case_id(i)).a;
    lambdap = wave(case_id(i)).lambda;
    kp = 2*pi/lambdap;
    cp = wave(case_id(i)).c;
    
    % Calculate u_star from c/u_star ratio
    u_star = cp / c_ustar(i);
    
    % Model parameters
    Ret = lambdap * u_star / nu;
    delta = 3 * Hs_ref / lambdap;
    psi = 0;
    
    % Velocity-dependent zo_filt
    switch vel_id(i)
        case 1  % u_inf = 2.42 m/s (Triangle)
            zo_filt = (7.8e-5/30)/lambdap;
        case 2  % u_inf = 3.87 m/s (Square)
            zo_filt = (2.57e-4/30)/lambdap;
        case 3  % u_inf = 5.43 m/s (Circle)
            zo_filt = (7.07e-4/30)/lambdap;
        otherwise
            error('Unknown vel_id: %d', vel_id(i));
    end
    
    flag_sgs = 1;
    max_iterations = 10^6;
    tolerance = 1e-6;
    initial_guess = 30;
    
    % Call your Newton-Raphson solver
    try
        [U] = solve_U_NewtonRaphson_simp_stokes(psi, Hs_ref*kp, cp/u_star, delta, Ret, zo_filt, flag_sgs, max_iterations, tolerance, initial_guess);
        
        % Calculate delta_U from your model
        lambda_coef = 1/U^2;
        zo_wave = delta * exp(-(kappa/sqrt(lambda_coef) + psi)) * lambdap;
        delta_U_model(i) = -(1/kappa)*log(nu/(zo_wave*u_star)) + 5;
    catch ME
        warning('Solver failed for point %d: %s', i, ME.message);
        delta_U_model(i) = NaN;
    end
end

%% ==================== CREATE FIGURE WITH TWO SUBPLOTS ====================



figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test7';


tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
annotationFontSize = 8;
marker_size = 25;

markerLinewidth = 1;
trendLinewidth = 1;

clc; close all
figure('Color','w','Units','centimeters','Position',[10 10 13 5.5])
tile = tiledlayout(1,2, 'padding', 'compact', 'TileSpacing', 'loose');


ax1 = nexttile;
box off
hold on
x_fit = linspace(1, 7, 200);
y_fit = 15.5 * x_fit.^(-0.77);
plot(x_fit, y_fit, 'k-', 'LineWidth', trendLinewidth);

% Experiments
for i = 1:length(c_ustar)
    scatter(c_ustar(i), delta_U(i), marker_size, ...
        'Marker', markers{vel_id(i)}, ...
        'MarkerFaceColor', wave(case_id(i)).color, ...
        'MarkerEdgeColor', 'none', ...
        'LineWidth', markerLinewidth);
end

% SWARL
for i = 1:length(c_ustar)
    if ~isnan(delta_U_model(i))
        scatter(c_ustar(i), delta_U_model(i), marker_size, ...
            'Marker', markers{vel_id(i)}, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', wave(case_id(i)).color, ...
            'LineWidth', markerLinewidth);
    end
end
hold off
box off

axis square
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% set(ax1, 'TickDir', 'none')
xlabel('$c \, / \, \overline{(u_\xi^*)}_{_\varphi}$', 'Interpreter', 'latex', 'FontSize', labelFontSize)
ylabel('$\overline{(\Delta U_\xi^+)}_{_\varphi}$', 'Interpreter', 'latex', 'FontSize', labelFontSize)
xlim([0.5 7])
ylim([0 13])



% Model vs Experiments
ax2 = nexttile;
hold on

% Plot data points
for i = 1:length(delta_U)
    if ~isnan(delta_U_model(i))

        scatter(delta_U(i), delta_U_model(i), marker_size, ...
            'Marker', markers{vel_id(i)}, ...
            'MarkerFaceColor', wave(case_id(i)).color, ...
            'MarkerEdgeColor', 'none', ...
            'LineWidth', markerLinewidth);
    end
end

% Diagonal line (perfect agreement)
delta_U_range = [0, 12];
P = plot(delta_U_range, delta_U_range, 'k-', 'LineWidth', trendLinewidth);
uistack(P, 'bottom')


legend_marker_size = 4;
% Legend handles - Wave cases (colors)
h1 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', wave(1).color, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', legend_marker_size);
h2 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', wave(2).color, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', legend_marker_size);
h3 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', wave(3).color, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', legend_marker_size);
h4 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', wave(4).color, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', legend_marker_size);


% Add gap in legend
hgap1 = plot(nan, nan, 'color', 'white', 'displayname', ' ');

% Wind speed markers
h5 = plot(nan, nan, '^', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', legend_marker_size);
h6 = plot(nan, nan, 's', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', legend_marker_size);
h7 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', legend_marker_size);

% Add gap in legend
hgap2 = plot(nan, nan, 'color', 'white', 'displayname', ' ');

% Legend for model vs experiments
% h8 = scatter(NaN, NaN, legend_marker_size, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% h9 = scatter(NaN, NaN, legend_marker_size, 'o', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'k', 'linewidth', 1.5);

h8 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerSize', legend_marker_size);
h9 = plot(nan, nan, 'o', ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', 'k', ...
    'linewidth', markerLinewidth, ...
    'MarkerSize', legend_marker_size);

hold off





box off
axis square
set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize);
xlabel('$\overline{(\Delta U_\xi^+)}_{_\varphi}$ (Experiments)', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$\overline{(\Delta U_\xi^+)}_{_\varphi}$ (SWARL)', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', labelFontSize);
xticks(ax2, 4:4:12)
yticks(ax2, 4:4:12)
xlim([2 13]);
ylim([2 13]);







% Calculate statistics
valid_idx = ~isnan(delta_U_model);
[corr_val, e1, e2, e3] = calculate_stats(delta_U_model(valid_idx), delta_U(valid_idx));

% Display statistics on plot
text(0.25, 0.92, ['$\rho = ' num2str(corr_val,'%.3f') '$'], 'FontSize', annotationFontSize, 'FontWeight', 'bold', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Units', 'normalized');
text(0.25, 0.82, ['$e = ' num2str(e2,'%.3f') '$'], 'FontSize', annotationFontSize, 'FontWeight', 'bold', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Units', 'normalized');

leg = legend([h1 h2 h3 h4 hgap1 h5 h6 h7 hgap2 h8 h9], ...
    {wave(1).name, wave(2).name, wave(3).name, wave(4).name, ...
     ' ', ...
     '$u_\infty = 2.42$ m/s', '$u_\infty = 3.87$ m/s', '$u_\infty = 5.43$ m/s', ...
     ' ', ...
     'Experiments', 'SWARL'}, ...
    'Interpreter', 'latex', 'FontSize', legendFontSize, 'orientation', 'vertical', 'Box','off');

leg.Layout.Tile = 'east';
leg.ItemTokenSize(1) = 10;
% leg.NumColumns = 4;


addPanelLabels([ax1, ax2], {'a', 'b'}, 'Offset', [-0.22, 1.12])


% Save figure
% pause(3)
% figure_name = 'LogLaw_Curvilinear_SWARL_Comparison.pdf';
% exportgraphics(tile, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 300, 'Resolution', 600, 'ContentType', 'image');
% close all
% fprintf('Generated figure: %s\n\n', figure_name)






%% ==================== PRINT COMPARISON TABLE ====================
fprintf('\n===== MODEL vs DATA COMPARISON =====\n')
fprintf('%-4s %-8s %-10s %-12s %-10s\n', 'Row', 'c/u*', 'dU+(data)', 'dU+(model)', 'Error(%)')
fprintf('%s\n', repmat('-', 1, 55))

for i = 1:length(c_ustar)
    if ~isnan(delta_U_model(i))
        err = 100*(delta_U_model(i) - delta_U(i))/delta_U(i);
        fprintf('%-4d %-8.2f %-10.2f %-12.2f %-10.1f\n', ...
            i, c_ustar(i), delta_U(i), delta_U_model(i), err)
    else
        fprintf('%-4d %-8.2f %-10.2f %-12s %-10s\n', ...
            i, c_ustar(i), delta_U(i), 'FAILED', '-')
    end
end

fprintf('\n===== STATISTICS =====\n')
fprintf('Correlation coefficient: %.4f\n', corr_val)
fprintf('e1 (mean log error): %.4f\n', e1)
fprintf('e2 (mean relative error): %.4f\n', e2)
fprintf('e3 (mean absolute error): %.4f\n', e3)
fprintf('\nTotal points: %d data + %d model = %d\n', ...
    length(c_ustar), sum(~isnan(delta_U_model)), length(c_ustar) + sum(~isnan(delta_U_model)))


%% ==================== STATISTICS FUNCTION ====================
function [corr_val, e1, e2, e3] = calculate_stats(predicted, measured)
    % Calculate correlation coefficient from scratch
    % Step 1: Calculate the means
    mean_predicted = mean(predicted);
    mean_measured = mean(measured);
    
    % Step 2: Calculate the deviations from the means
    dev_predicted = predicted - mean_predicted;
    dev_measured = measured - mean_measured;
    
    % Step 3: Calculate the numerator (covariance)
    numerator = sum(dev_predicted .* dev_measured);
    
    % Step 4: Calculate the denominator (product of standard deviations)
    std_predicted = sqrt(sum(dev_predicted.^2));
    std_measured = sqrt(sum(dev_measured.^2));
    denominator = std_predicted * std_measured;
    
    % Step 5: Calculate the correlation coefficient
    corr_val = numerator / denominator;
    
    % Calculate error metrics
    e1 = mean(abs(log10(predicted ./ measured)));
    e2 = mean(abs(predicted - measured) ./ measured);
    e3 = mean(abs(predicted - measured));
end





function addPanelLabels(ax, labels, varargin)
% addPanelLabels(ax, labels) adds (a),(b),... just OUTSIDE top-left of each axes.
% ax     : array of axes handles (e.g., from tiledlayout / findall)
% labels : cellstr like {'a','b','c'} or string array ["a" "b" "c"]
%
% Optional name-value:
% 'Offset'   : [dx dy] in normalized axes units (default [-0.10 1.02])
% 'FontSize' : default 12
% 'FontName' : default 'Times New Roman'

p = inputParser;
addParameter(p,'Offset',[-0.10 1.1]);
addParameter(p,'FontSize', 10);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});
off = p.Results.Offset;

labels = string(labels);
for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Plain parentheses + italic letter:
    % TeX interpreter: \it turns italic ON, \rm returns to roman.
    s = sprintf('(\\ita\\rm)');              % placeholder
    s = sprintf('(\\it%s\\rm)', labels(k));  % actual label

    text(ax(k), off(1), off(2), s, ...
        'Units','normalized', ...
        'Interpreter','tex', ...           % keeps italics control simple
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');                 % critical: allow outside axes
    end
end




