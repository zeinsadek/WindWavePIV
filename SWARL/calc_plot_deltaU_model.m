clear; clc;
addpath('/Users/zeinsadek/Downloads')
%% ==================== WAVE PARAMETERS FROM TABLE ====================
wave = struct();

wave(1).name = '$\lambda_{410}$, $ak_{0.180}$';
wave(1).lambda = 410.8e-3;   wave(1).a = 11.8e-3;
wave(1).c = 0.81;            wave(1).ak = 0.180;
wave(1).color = [0.9 0.1 0.1];  % Red

wave(2).name = '$\lambda_{313}$, $ak_{0.211}$';
wave(2).lambda = 313.3e-3;   wave(2).a = 10.5e-3;
wave(2).c = 0.71;            wave(2).ak = 0.211;
wave(2).color = [1.0 0.75 0.0]; % Gold/Yellow

wave(3).name = '$\lambda_{189}$, $ak_{0.305}$';
wave(3).lambda = 189.6e-3;   wave(3).a = 9.2e-3;
wave(3).c = 0.57;            wave(3).ak = 0.305;
wave(3).color = [0.0 0.6 0.2]; % Green

wave(4).name = '$\lambda_{124}$, $ak_{0.267}$';
wave(4).lambda = 124.3e-3;   wave(4).a = 5.3e-3;
wave(4).c = 0.49;            wave(4).ak = 0.267;
wave(4).color = [0.0 0.8 0.9]; % Cyan

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
figure(1)
set(gcf,'Color','w','Units','centimeters','Position',[20 20 40 15])

% Subplot positions
xloc1 = 4;  yloc = 4;  xsize = 12; ysize = 10;
xloc2 = xloc1 + xsize + 5;  % Second plot to the right

ax1 = axes;
ax1.Units = 'Centimeters';
hold on

x_fit = linspace(1, 7, 200);
y_fit = 15.5 * x_fit.^(-0.77);
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5);

for i = 1:length(c_ustar)
    scatter(c_ustar(i), delta_U(i), marker_size, ...
        'Marker', markers{vel_id(i)}, ...
        'MarkerFaceColor', wave(case_id(i)).color, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.5);
end

for i = 1:length(c_ustar)
    if ~isnan(delta_U_model(i))
        scatter(c_ustar(i), delta_U_model(i), marker_size, ...
            'Marker', markers{vel_id(i)}, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', wave(case_id(i)).color, ...
            'LineWidth', 1.5);
    end
end



hold off

box on
set(ax1, 'LineWidth', 1.5);
ax1.Layer = 'top';
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', 20)
set(ax1, 'TickDir', 'none')
xlabel('$c \, / \, \overline{(u_\xi^*)}_{_\varphi}$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$\overline{(\Delta U_\xi^+)}_{_\varphi}$', 'Interpreter', 'latex', 'FontSize', 20)
xlim([1 7])
ylim([0 16])
ax1.Position = [xloc1 yloc xsize ysize];

ax2 = axes;
ax2.Units = 'Centimeters';
hold on

% Plot data points
for i = 1:length(delta_U)
    if ~isnan(delta_U_model(i))
        plot(delta_U(i), delta_U_model(i), markers{vel_id(i)}, ...
            'LineWidth', 1.5, ...
            'MarkerSize', 10, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', wave(case_id(i)).color);
    end
end

% Diagonal line (perfect agreement)
delta_U_range = [0, 12];
P = plot(delta_U_range, delta_U_range, 'k-', 'LineWidth', 1.5);
uistack(P, 'bottom')


% Legend handles - Wave cases (colors)
h1 = scatter(NaN, NaN, 80, 'o', 'MarkerFaceColor', wave(1).color, 'MarkerEdgeColor', wave(1).color);
h2 = scatter(NaN, NaN, 80, 'o', 'MarkerFaceColor', wave(2).color, 'MarkerEdgeColor', wave(2).color);
h3 = scatter(NaN, NaN, 80, 'o', 'MarkerFaceColor', wave(3).color, 'MarkerEdgeColor', wave(3).color);
h4 = scatter(NaN, NaN, 80, 'o', 'MarkerFaceColor', wave(4).color, 'MarkerEdgeColor', wave(4).color);
% Velocity markers (shapes)
h5 = scatter(NaN, NaN, 80, '^', 'k', 'MarkerFaceColor', 'k');
h6 = scatter(NaN, NaN, 80, 's', 'k', 'MarkerFaceColor', 'k');
h7 = scatter(NaN, NaN, 80, 'o', 'k', 'MarkerFaceColor', 'k');

hold off

box on
set(ax2, 'LineWidth', 1.5);
set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', 20);

xlabel('$\overline{(\Delta U_\xi^+)}_{_\varphi}$ (Ref)', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\overline{(\Delta U_\xi^+)}_{_\varphi}$ (SWARL)', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 20);

xlim([2 12]);
ylim([2 12]);

% Calculate statistics
valid_idx = ~isnan(delta_U_model);
[corr_val, e1, e2, e3] = calculate_stats(delta_U_model(valid_idx), delta_U(valid_idx));

% Display statistics on plot
text(0.25, 0.92, ['$\rho = ' num2str(corr_val,'%.4f') '$'], 'FontSize', 20, 'FontWeight', 'bold', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Units', 'normalized');
text(0.25, 0.82, ['$e = ' num2str(e2,'%.4f') '$'], 'FontSize', 20, 'FontWeight', 'bold', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Units', 'normalized');

legend([h1 h2 h3 h4 h5 h6 h7 ], ...
    {wave(1).name, wave(2).name, wave(3).name, wave(4).name, ...
     '$u_\infty = 2.42$ m/s', '$u_\infty = 3.87$ m/s', '$u_\infty = 5.43$ m/s'}, ...
    'Interpreter', 'latex', 'FontSize', 18, 'Location', 'northeastoutside','Box','off')

ax2.Position = [xloc2 yloc xsize ysize];

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