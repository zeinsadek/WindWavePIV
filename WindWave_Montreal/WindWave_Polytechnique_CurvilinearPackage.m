% Wrap up the curvilinear-projected instantaneous and phase averages for
% Montreal Polytechnique


clc; close all; clear

% Get case to process
wind_speed = 'WT6';
wave_type = 'C';
caze = [wind_speed, '_WV', wave_type, '_AG0'];

% Load all useful data
data   = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/data/', caze, '_DATA.mat'));
waves  = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/', caze, '_WAVE.mat'));
means  = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/', caze, '_MEANS.mat'));
phases = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/phase/', caze, '_PHASE.mat'));

% Rename structure
data   = data.output;
waves  = waves.output;
means  = means.output;
phases = phases.output;
clc; fprintf('Data loaded\n\n')

% Get wave parameters
wave_parameters = readcell("/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2}; %#ok<*FNDSB>
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
frequency       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
phase_offset    = [0, wavelength/4, wavelength/2, 3*wavelength/4];
wavenumber      = (2 * pi) / wavelength;

% Get freestream
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

clear wave_type


%% Load relevent data

% Coordinates
X = data.X;
Y = data.Y;
x = X(1,:);
y = Y(:,1);

% Differential elements
dx = mean(diff(x));
dy = mean(diff(y));

% Bin instantaneous snapshots by phase
for phase = 1:4

    % Indicies of images (out of the 6000 frames) that belong to 'phase'
    indicies = find(phases.phase_average_idx == phase);
    fprintf('Phase %1.0f: %4.0f Images\n\n', phase, length(indicies))

    % Save Cartesian, phase-binned velocities
    cartesian_phaseBinnedInstantaneous(phase).U     = data.U(:,:,indicies); %#ok<*SAGROW>
    cartesian_phaseBinnedInstantaneous(phase).V     = data.V(:,:,indicies);
    cartesian_phaseBinnedInstantaneous(phase).waves = imresize(waves.wave_profiles(indicies,:), [length(indicies), length(x)]);
end

clear phases wave_phases_indicies waves phase indicies


%% Repair cat-scratch in Cartesian MEANS

% Components to fix
components = {'u', 'v', 'uu', 'vv', 'uv'};

clc;
% Loop through phases
for phase = 1:4

    % Loop through components
    for c = 1:length(components)
        component = components{c};
    
        % Load component
        tmp = means.phase(phase).(component);
        max_wave_profile = means.phase(phase).max_wave_profile;
    
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
        corrected(Y < max_wave_profile) = nan;
        
        % Make sure original nans are in place
        corrected(isnan(tmp)) = nan;
    
        % Resave data
        cartesian_repairedMeans(phase).(component) = corrected;
    end
end

clear tmp LHS_values RHS_values differences left right 
clear max_wave_profile phase corrected c component


%% Compute curvilinear grid pher phase to be used in each step

% Centering wave profile around data
x_shift = mean(x);
top_crop = 210;

% Loop through different phases
for phase = 1:4
    
    % Specific wave shape at each phase
    max_wave_profile = means.phase(phase).max_wave_profile;
    wave_profile = amplitude * cos(2 * pi * (x - x_shift - phase_offset(phase)) / wavelength);
    wave_orthogonal = (2 * pi / wavelength) * amplitude * sin(2 * pi * (x - x_shift - phase_offset(phase)) / wavelength);

    clc;
    % Output size will match input PIV grid
    [rows, columns]  = size(X);
    horizontal_lines = nan(rows, columns);
    vertical_lines   = nan(rows, columns);
    
    % Horizontal (Zeta) Lines
    for r = 1:rows
        y0 = y(r);
        horizontal_lines(r, :) = y0 + wave_profile .* exp(-wavenumber * y0);
    end
    
    % Build Vertical (Xi) Lines 
    for c = 1:columns
        x0 = x(c);
        % Zero to wave surface
        y_zeroed = y - wave_profile(c);
        vertical_lines(:, c) = (y_zeroed .* wave_orthogonal(c) .* exp(-wavenumber * y_zeroed)) + x0;
    end
    
    % Crop below wave profile
    vertical_lines(horizontal_lines < wave_profile) = nan;
    horizontal_lines(horizontal_lines < wave_profile) = nan;

     % Compute slope of tangent lines across curvilinear space
    [dz_dx, dz_dy] = gradient(horizontal_lines, dx, dy);
    dz_dx(horizontal_lines < wave_profile) = nan;
    dz_dy(horizontal_lines < wave_profile) = nan;
    angles = atan2(dz_dx, dz_dy);
    angles(horizontal_lines < wave_profile) = nan;

    % Save
    curvilinear_mesh(phase).vertical_lines = vertical_lines;
    curvilinear_mesh(phase).horizontal_lines = horizontal_lines;
    curvilinear_mesh(phase).wave_profile = wave_profile;
    curvilinear_mesh(phase).max_wave_profile = max_wave_profile;
    curvilinear_mesh(phase).angles = angles;

end

clear r c x0 y0 y_zeroed wave_orthogonal wave_profile phase x_shift
clear angles dz_dx dz_dy vertical_lines horizontal_lines max_wave_profile
clear rows columns dx dy phase_offset


%% Plot curvilinear grids to verify they are correct

% % Plot settings
% pad = 10;
% sz = 10;
% spacing = 5;
% colors = parula(numel(vertical_lines(1:spacing:end, 1:spacing:end)));
% 
% % Figure
% figure('color', 'white')
% tile = tiledlayout(1,4);
% 
% % Loop through phases
% for phase = 1:4
% 
%     h(phase) = nexttile;
%     hold on
% 
%     % Plot curves
%     for i = 1:spacing:rows
%         plot(x, curvilinear_mesh(phase).horizontal_lines(i,:), 'color', 'black')
%     end
%     for i = 1:spacing:columns
%         plot(curvilinear_mesh(phase).vertical_lines(:,i), y + curvilinear_mesh(phase).wave_profile(i), 'color', 'black')
%     end
% 
%     % Plot scatter of points that we will interpolate at
%     c = 1;
%     for i = 1:spacing:rows
%         for j = 1:spacing:columns
%             scatter(curvilinear_mesh(phase).vertical_lines(i,j), curvilinear_mesh(phase).horizontal_lines(i,j), sz, 'filled', 'MarkerFaceColor', colors(c,:))
%             c = c + 1;
%         end
%     end
% 
%     % Plot wave surface
%     plot(x, curvilinear_mesh(phase).max_wave_profile, 'linewidth', 3, 'Color', 'blue')
%     plot(x, curvilinear_mesh(phase).wave_profile, 'linewidth', 3)
%     hold off
% 
%     title(sprintf('Phase %1.0f', phase))
%     axis equal
%     xline([min(x), max(x)])
%     yline([min(y), max(y)])
% end
% 
% linkaxes(h, 'xy')
% xlim([min(x) - pad, max(x) + pad])
% ylim([min(y) - pad, max(y) + pad])
% 
% clear c i j pad sz colors spacing 


%% Compute MEAN curvilinear velocities and stresses

% Interpolate Cartesian quantities onto curvilinear grid
for phase = 1:4

    % Load data
    u  = cartesian_repairedMeans(phase).u;
    v  = cartesian_repairedMeans(phase).v;
    uu = cartesian_repairedMeans(phase).uu;
    vv = cartesian_repairedMeans(phase).vv;
    uv = cartesian_repairedMeans(phase).uv;

    % Load curvilinear mesh
    vertical_lines = curvilinear_mesh(phase).vertical_lines;
    horizontal_lines = curvilinear_mesh(phase).horizontal_lines;
    max_wave_profile = curvilinear_mesh(phase).max_wave_profile;

    % Fill NaNs in original data just in case
    u_filled  = u;
    v_filled  = v;
    uu_filled = uu;
    vv_filled = vv;
    uv_filled = uv;

    u_filled(isnan(u_filled))   = 0;
    v_filled(isnan(v_filled))   = 0;
    uu_filled(isnan(uu_filled)) = 0;
    vv_filled(isnan(vv_filled)) = 0;
    uv_filled(isnan(uv_filled)) = 0;
    
    % Interpolate
    u_interp  = interp2(X, Y, u_filled,  vertical_lines, horizontal_lines, 'cubic', nan);
    v_interp  = interp2(X, Y, v_filled,  vertical_lines, horizontal_lines, 'cubic', nan);
    uu_interp = interp2(X, Y, uu_filled, vertical_lines, horizontal_lines, 'cubic', nan);
    vv_interp = interp2(X, Y, vv_filled, vertical_lines, horizontal_lines, 'cubic', nan);
    uv_interp = interp2(X, Y, uv_filled, vertical_lines, horizontal_lines, 'cubic', nan);

    u_interp(u_interp == 0)   = nan;
    v_interp(v_interp == 0)   = nan;
    uu_interp(uu_interp == 0) = nan;
    vv_interp(vv_interp == 0) = nan;
    uv_interp(vv_interp == 0) = nan;

    % Crop wave
    u_interp(horizontal_lines  < max_wave_profile) = nan;
    v_interp(horizontal_lines  < max_wave_profile) = nan;
    uu_interp(horizontal_lines < max_wave_profile) = nan;
    vv_interp(horizontal_lines < max_wave_profile) = nan;
    uv_interp(horizontal_lines < max_wave_profile) = nan;
    
    % Crop top
    u_interp(horizontal_lines  > top_crop) = nan;
    v_interp(horizontal_lines  > top_crop) = nan;
    uu_interp(horizontal_lines > top_crop) = nan;
    vv_interp(horizontal_lines > top_crop) = nan;
    uv_interp(horizontal_lines > top_crop) = nan;

    % Save interpolated values
    cartesian_interpolatedMeans(phase).u  = u_interp;
    cartesian_interpolatedMeans(phase).v  = v_interp;
    cartesian_interpolatedMeans(phase).uu = uu_interp;
    cartesian_interpolatedMeans(phase).vv = vv_interp;
    cartesian_interpolatedMeans(phase).uv = uv_interp;


end

clear u_filled v_filled uu_filled vv_filled uv_filled vertical_lines horizontal_lines
clear u_interp v_interp uu_interp vv_interp uv_interp u v uu vv uv phase



% Project velocities and stresses onto curvilinear grid
for phase = 1:4

    % Load interpolated velocities
    u_interp  = cartesian_interpolatedMeans(phase).u;
    v_interp  = cartesian_interpolatedMeans(phase).v;
    uu_interp = cartesian_interpolatedMeans(phase).uu;
    vv_interp = cartesian_interpolatedMeans(phase).vv;
    uv_interp = cartesian_interpolatedMeans(phase).uv;

    % Load curvilinear mesh
    horizontal_lines = curvilinear_mesh(phase).horizontal_lines;
    max_wave_profile = curvilinear_mesh(phase).max_wave_profile;
    angles = curvilinear_mesh(phase).angles;
    
    % Compute unit vectors
    xi_hat_x = cos(angles);    
    xi_hat_y = sin(angles);
    zeta_hat_x = -sin(angles);   
    zeta_hat_y =  cos(angles);

    % Project Cartesian velocities into curvilinear space
    u_tangent = u_interp .* xi_hat_x + v_interp .* xi_hat_y;
    u_normal  = (u_interp .* zeta_hat_x + v_interp .* zeta_hat_y);

    % Project Cartesian stresses into curvilinear space
    tau_xixi   = uu_interp .* xi_hat_x.^2 + 2 * uv_interp .* xi_hat_x .* xi_hat_y + vv_interp .* xi_hat_y.^2;
    tau_zetazeta = uu_interp .* zeta_hat_x.^2 + 2 * uv_interp .* zeta_hat_x .* zeta_hat_y + vv_interp .* zeta_hat_y.^2;
    tau_xizeta = uu_interp .* xi_hat_x .* zeta_hat_x + ...
                 uv_interp .* (xi_hat_x .* zeta_hat_y + xi_hat_y .* zeta_hat_x) + ...
                 vv_interp .* xi_hat_y .* zeta_hat_y;

    % Crop below wave
    u_tangent(horizontal_lines    < max_wave_profile) = nan;
    u_normal(horizontal_lines     < max_wave_profile) = nan;
    tau_xixi(horizontal_lines     < max_wave_profile) = nan;
    tau_zetazeta(horizontal_lines < max_wave_profile) = nan;
    tau_xizeta(horizontal_lines   < max_wave_profile) = nan;
    
    % Crop top
    u_tangent(horizontal_lines    > top_crop) = nan;
    u_normal(horizontal_lines     > top_crop) = nan;
    tau_xixi(horizontal_lines     > top_crop) = nan;
    tau_zetazeta(horizontal_lines > top_crop) = nan;
    tau_xizeta(horizontal_lines   > top_crop) = nan;

    % Save projected curvilinear mean values
    curvilinear_projectedMeans(phase).u  = u_tangent;
    curvilinear_projectedMeans(phase).v  = u_normal;
    curvilinear_projectedMeans(phase).uu = tau_xixi;
    curvilinear_projectedMeans(phase).vv = tau_zetazeta;
    curvilinear_projectedMeans(phase).uv = tau_xizeta;

end

clear u_interp v_interp uu_interp vv_interp uv_interp phase
clear u_tangent u_normal tau_xixi tau_zetazeta tau_xizeta
clear angles xi_hat_x xi_hat_y zeta_hat_x zeta_hat_y dz_dx dz_dy
clear wave_profile max_wave_profile horizontal_lines vertical_lines



% Plot curvilinear velocities and stresses to make sure we are okay
% Loop through components
for c = 1:length(components)
    component = components{c};

    % Figure
    figure('color', 'white')
    tiledlayout(1,4)
    sgtitle(component)

    % Loop through phases
    for phase = 1:4

        % Get coordinates
        vertical_lines = curvilinear_mesh(phase).vertical_lines;
        horizontal_lines = curvilinear_mesh(phase).horizontal_lines;
        wave_profile = curvilinear_mesh(phase).wave_profile;
        max_wave_profile = curvilinear_mesh(phase).max_wave_profile;

        % Plot
        h(phase) = nexttile;
        hold on
        contourf(vertical_lines, horizontal_lines, curvilinear_projectedMeans(phase).(component), 20, 'linestyle', 'none')
        plot(x, wave_profile, 'linewidth', 2, 'color', 'black')
        plot(x, max_wave_profile, 'linewidth', 2, 'color', 'red')
        hold off
        axis equal
        title(sprintf('Phase %1.0f', phase))
        colorbar()
    end
    linkaxes(h, 'xy')
    ylim([-20, 200])
end

clear c component phase vertical_lines horizontal_lines wave_profile
clear max_wave_profile tmp h 


%% Compute INSTANTANEOUS curvilinear velocities


% [~, left] = min(abs(x - 115));
% [~, right] = min(abs(x - 120));

clc; close all
% Loop through phases
for phase = 1:4

    % Get coordinates
    vertical_lines   = curvilinear_mesh(phase).vertical_lines;
    horizontal_lines = curvilinear_mesh(phase).horizontal_lines;
    max_wave_profile = curvilinear_mesh(phase).max_wave_profile;
    angles           = curvilinear_mesh(phase).angles;
    
    % Get instantaneous snapshots
    u_phase_snapshots = cartesian_phaseBinnedInstantaneous(phase).U;
    v_phase_snapshots = cartesian_phaseBinnedInstantaneous(phase).V;
    wave_snapshots    = cartesian_phaseBinnedInstantaneous(phase).waves;
    num_images        = size(u_phase_snapshots, 3);

    % Compute unit vectors
    xi_hat_x = cos(angles);    
    xi_hat_y = sin(angles);
    zeta_hat_x = -sin(angles);   
    zeta_hat_y =  cos(angles);

    fprintf('Phase %1.0f\n\n', phase)
    % Loop through snapshots
    for i = 1:num_images
    % for i = 1

        % Progress bar
        progressbarText(i / num_images);

        %%% Get instantanrous u, v, and wave profile
        u_inst    = u_phase_snapshots(:,:,i);
        v_inst    = v_phase_snapshots(:,:,i);
        wave_inst = wave_snapshots(i,:);

        % Turn NaNs we want to keep into 0
        u_inst(Y < wave_inst) = 0;
        v_inst(Y < wave_inst) = 0;

        % Fill in NaNs in the flow
        u_inst_NaN_filled = inpaint_nans(double(u_inst));
        v_inst_NaN_filled = inpaint_nans(double(v_inst));

        % Return below the wave to NaNs
        u_inst_NaN_filled(Y < wave_inst) = nan;
        v_inst_NaN_filled(Y < wave_inst) = nan;


        %%% Fix cat scratch in the instantaneous
        % u velocity snapshot
        % u_tmp = u_inst_NaN_filled;
        % LHS_values = u_tmp(:,86);
        % RHS_values = u_tmp(:,87);
        % differences = RHS_values - LHS_values;
        % differences(isnan(differences)) = 0;
        % u_corrected = u_tmp;
        % u_corrected(:,87:end) = u_corrected(:,87:end) - differences;
        % 
        % % v velocity snapshot
        % v_tmp = v_inst_NaN_filled;
        % LHS_values = v_tmp(:,86);
        % RHS_values = v_tmp(:,87);
        % differences = RHS_values - LHS_values;
        % differences(isnan(differences)) = 0;
        % v_corrected = v_tmp;
        % v_corrected(:,87:end) = v_corrected(:,87:end) - differences;
        % 
        % 
        % % Interpolate across strip
        % u_corrected(:, left:right) = nan;
        % u_corrected = fillmissing(u_corrected, 'makima', 2);
        % u_corrected(Y < wave_inst) = nan;
        % v_corrected(:, left:right) = nan;
        % v_corrected = fillmissing(v_corrected, 'makima', 2);
        % v_corrected(Y < wave_inst) = nan;
        % 
        % % Make sure original nans are in place
        % u_corrected(isnan(u_tmp)) = nan;
        % v_corrected(isnan(v_tmp)) = nan;
        % u_inst_NaN_filled = u_corrected;
        % v_inst_NaN_filled = v_corrected;


        %%% Interpolate onto curvilinear grid
        u_inst_filled = u_inst_NaN_filled;
        v_inst_filled = v_inst_NaN_filled;
        u_inst_filled(isnan(u_inst_filled)) = 0;
        v_inst_filled(isnan(v_inst_filled)) = 0;

        u_inst_interp = interp2(X, Y, u_inst_filled,  vertical_lines, horizontal_lines, 'cubic', nan);
        v_inst_interp = interp2(X, Y, v_inst_filled,  vertical_lines, horizontal_lines, 'cubic', nan);

        % Return NaNs
        u_inst_interp(u_inst_interp == 0) = nan;
        v_inst_interp(v_inst_interp == 0) = nan;

        % Crop to instantaneous wave profile
        u_inst_interp(horizontal_lines < wave_inst) = nan;
        v_inst_interp(horizontal_lines < wave_inst) = nan;

        % Crop to max wave profile in phase-bin
        u_inst_interp_max_crop = u_inst_interp;
        v_inst_interp_max_crop = v_inst_interp;
        u_inst_interp_max_crop(horizontal_lines < max_wave_profile) = nan;
        v_inst_interp_max_crop(horizontal_lines < max_wave_profile) = nan;

        

        %%% Project velocities
        % Instantaneous crop 
        u_tangent = u_inst_interp .* xi_hat_x + v_inst_interp .* xi_hat_y;
        u_normal  = u_inst_interp .* zeta_hat_x + v_inst_interp .* zeta_hat_y;

        % Max wave profile crop
        u_tangent_max_crop = u_inst_interp_max_crop .* xi_hat_x + v_inst_interp_max_crop .* xi_hat_y;
        u_normal_max_crop  = u_inst_interp_max_crop .* zeta_hat_x + v_inst_interp_max_crop .* zeta_hat_y;

        % Save
        % Instantaneous crop 
        curvilinear_instantaneousProjected(phase).instantaneous_crop.u(:,:,i) = u_tangent;
        curvilinear_instantaneousProjected(phase).instantaneous_crop.v(:,:,i) = u_normal;

        % Max wave profile crop
        curvilinear_instantaneousProjected(phase).maximum_crop.u(:,:,i) = u_tangent_max_crop;
        curvilinear_instantaneousProjected(phase).maximum_crop.v(:,:,i) = u_normal_max_crop; 

        % Save waves
        curvilinear_instantaneousProjected(phase).waves(i,:) = wave_inst;



        %%% Plot for checking
        % figure('color', 'white')
        % tiledlayout(1,2)
        % sgtitle(sprintf('Phase %1.0f', phase))
        % 
        % h(1) = nexttile;
        % hold on
        % contourf(X, Y, v_inst_NaN_filled)
        % plot(x, wave_profile, 'linewidth', 2, 'color', 'black')
        % plot(x, wave_inst', 'linewidth', 2, 'color', 'red')
        % hold off
        % axis equal
        % 
        % h(2) = nexttile;
        % hold on
        % contourf(vertical_lines, horizontal_lines, u_normal)
        % plot(x, wave_profile, 'linewidth', 2, 'color', 'black')
        % plot(x, wave_inst', 'linewidth', 2, 'color', 'red')
        % hold off
        % axis equal
        % linkaxes(h, 'xy')

    end

    % Save number of images
    curvilinear_instantaneousProjected(phase).num_images = num_images;

    clc;
end
clc;

clear left right phase i vertical_lines horizontal_lines wave_profile max_wave_profile tmp h 
clear angles xi_hat_x xi_hat_y zeta_hat_x zeta_hat_y u_phase_snapshots v_phase_snapshots
clear wave_snapshots num_images u_inst_filled u_inst u_inst_interp u_inst_interp_max_crop
clear u_inst_NaN_filled u_normal u_normal_max_crop u_tangent u_tangent_max_crop v_inst v_inst_filled
clear v_inst_interp v_inst_interp_max_crop v_inst_NaN_filled wave_inst wave_snapshots


%% Compute means and stresses from curvilinear instantaneous to compare

% Which crop method to compute
crop_method = 'instantaneous_crop';

% Loop through phases
for phase = 1:4

    % Compute means and fluctuations from instantaneous curvilinear fields
    mean_u = mean(curvilinear_instantaneousProjected(phase).(crop_method).u, 3);
    mean_v = mean(curvilinear_instantaneousProjected(phase).(crop_method).v, 3);
    u_flucs = curvilinear_instantaneousProjected(phase).(crop_method).u - mean_u;
    v_flucs = curvilinear_instantaneousProjected(phase).(crop_method).v - mean_v;
    
    % Compute stress
    uu = mean(u_flucs .* u_flucs, 3);
    vv = mean(v_flucs .* v_flucs, 3);
    uv = mean(u_flucs .* v_flucs, 3);

    % Save
    curvilinear_meansFromInstantaneous(phase).u = mean_u;
    curvilinear_meansFromInstantaneous(phase).v = mean_v;
    curvilinear_meansFromInstantaneous(phase).uu = uu;
    curvilinear_meansFromInstantaneous(phase).vv = vv;
    curvilinear_meansFromInstantaneous(phase).uv = uv;
end


% % Compare to means/stresses transformed from Cartesian means
% diff_uu = curvilinear_projectedMeans(phase).uu - uu;
% diff_vv = curvilinear_projectedMeans(phase).vv - vv;
% diff_uv = curvilinear_projectedMeans(phase).uv - uv;

clear phase crop_method mean_u mean_v u_flucs v_flucs uu vv uv h tile

%% Plot means and stresses from curvilinear instantaneous to compare

% Component to plot
component = 'vv';

% Figure
figure('color', 'white')
tiledlayout(1,4);

% Loop through phases
for phase = 1:4
  
    h(phase) = nexttile;
    hold on
    contourf(curvilinear_mesh(phase).vertical_lines, ...
             curvilinear_mesh(phase).horizontal_lines, ...
             curvilinear_meansFromInstantaneous(phase).(component), ...
             100, 'linestyle', 'none')
    plot(x, curvilinear_mesh(phase).wave_profile, 'linewidth', 2, 'color', 'black')
    plot(x, curvilinear_mesh(phase).max_wave_profile, 'linewidth', 2, 'color', 'red')
    hold off
    title(sprintf('Phase %1.0f', phase))
    axis equal
    colorbar()
    ylim([-20, 200])
end

linkaxes(h, 'xy')

clear phase component h


%% Collect all useful quantities to save into a matfile

% Curvilinear quantities
output.curvilinear_grid = curvilinear_mesh;
output.curvilinear_means = curvilinear_projectedMeans;
output.curvilinear_instantaneous = curvilinear_instantaneousProjected;

% Cartesian quantities
output.Cartesian_grid.X = X;
output.Cartesian_grid.Y = Y;
output.Cartesian_grid.x = x;
output.Cartesian_grid.y = y;

% Wave properties
output.constants.wavelength_mm = wavelength;
output.constants.amplitude_mm = amplitude;
output.constants.frequency_hz = frequency;
output.constants.steepness = (2 * pi * amplitude) / wavelength;
output.constants.wave_speed = frequency * wavelength * 1E-3;

% Freestream
output.constants.freestream = freestreams.(wind_speed);



%%% Save to matfile
save_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/curvilinear_montreal';
caze_name   = [caze, '_CurvilinearInstantaneous_Polytechnique.mat'];
clc; fprintf('Saving Matfile...\n')
pause(3)
save(fullfile(save_folder, caze_name), 'output', '-v7.3')
clc; fprintf('Matfile saved!\n')

