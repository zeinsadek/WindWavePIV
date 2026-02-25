%% Plotting z0 Parameters

clear; close all; clc;
data = load('/Users/zeinsadek/Desktop/Experiments/Offshore/Hopkins/Z0/PSU_z0.mat');
data = data.output;

cases = fields(data);
cases_w_waves = cases(find(contains(cases, 'AG0')));


% figure()
% hold on
for i = 1:length(cases_w_waves)
    caze = cases{i};
    Re(i) = data.(caze).Re;
    ak(i) = data.(caze).ak;
    z0(i) = data.(caze).z0;
    % scatter3(data.(caze).Re, data.(caze).ak, data.(caze).z0, 'filled')
end
% hold off
% view(40,35)
% grid on


xlin  = linspace(min(Re), max(Re), 100);
ylin  = linspace(min(ak), max(ak), 100);
[X,Y] = meshgrid(xlin, ylin);
Z     = griddata(Re,ak,z0,X,Y,'v4');
figure()
mesh(X,Y,Z)






