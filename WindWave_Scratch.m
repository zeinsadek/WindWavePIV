% %% Compare Dewraped high-res image to PIV to check vertical misalignment
% 
% clear;
% addpath('C:/Users/Zein/Documents/MATLAB/readimx-v2.1.8-win64/');
% %%
% 
% raw_path = ['N:\WT7_WVD_AGS30\WT7_WVD_AGS30\PIV_MP(2x24x24_50%ov)\PostProc\B05997.vc7'];
% % PIV_path = 'N:\WT15_WVA_AGS30\WT15_WVA_AGS30\PIV_MP(2x24x24_50%ov)\PostProc\B00001.vc7';
% 
% 
% % Define Image Dimensions
% raw       = readimx(raw_path);
% raw_image = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
% 
% % PIV       = readimx(PIV_path);
% % PIV_image = PIV.Frames{1,1}.Components{1,1}.Planes{1,1};
% 
% 
% nf_raw = size(raw_image);
% x_raw  = raw.Frames{1,1}.Scales.X.Slope.*linspace(1, nf_raw(1), nf_raw(1)).* raw.Frames{1,1}.Grids.X + raw.Frames{1,1}.Scales.X.Offset;
% y_raw  = raw.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf_raw(2), nf_raw(2)).* raw.Frames{1,1}.Grids.Y + raw.Frames{1,1}.Scales.Y.Offset;
% [X_raw, Y_raw] = meshgrid(x_raw, y_raw);
% 
% % nf_PIV = size(PIV_image);
% % x_PIV  = PIV.Frames{1,1}.Scales.X.Slope.*linspace(1, nf_PIV(1), nf_PIV(1)).* PIV.Frames{1,1}.Grids.X + PIV.Frames{1,1}.Scales.X.Offset;
% % y_PIV  = PIV.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf_PIV(2), nf_PIV(2)).* PIV.Frames{1,1}.Grids.Y + PIV.Frames{1,1}.Scales.Y.Offset;
% % [X_PIV, Y_PIV] = meshgrid(x_PIV, y_PIV);
% 
% raw_image = rot90(raw_image, -1);
% % PIV_image = rot90(PIV_image, -1);
% 
% % vertical_shift = abs(min(Y_raw, [], 'all') - min(Y_PIV, [], 'all'));
% 
% figure()
% contourf(X_raw, Y_raw, raw_image, 100, 'linestyle', 'none')
% yline(-104.1)
% axis equal
% colorbar()
% xlim([-150, 150])
% ylim([-150, 150])
% title('A')
% 
% % t =tiledlayout(1,2);
% % nexttile()
% % contourf(X_raw, Y_raw, raw_image, 100, 'linestyle', 'none')
% % yline(-104.1)
% % axis equal
% % colorbar()
% % xlim([-150, 150])
% % ylim([-150, 150])
% % title('A')
% % 
% % nexttile()
% % contourf(X_PIV, Y_PIV, PIV_image, 100, 'linestyle', 'none')
% % yline(-104.1)
% % axis equal
% % colorbar()
% % xlim([-150, 150])
% % ylim([-150, 150])
% % title('B')
% 
% 
% % figure()
% % hold on
% % contourf(X_raw, Y_raw, raw_image, 100, 'linestyle', 'none')
% % colormap('parula');
% % [c, h] = contourf(X_PIV, Y_PIV - vertical_shift, PIV_image, 10, '--w');
% % h.FaceAlpha = 0.0;
% % hold off
% % axis equal
% 
% % figure()
% % hold on
% % contourf(X_raw, Y_raw, raw_image)
% % % contourf(X_raw_resize, Y_raw_resize, raw_image_resize)
% % [c, h] = contourf(X_PIV, Y_PIV - vertical_shift, PIV_image, 10, '--w');
% % % [c,h] = contourf(X_raw_resize, Y_raw_resize, raw_image_resize, 20, '--w');
% % colormap('parula');
% % h.FaceAlpha = 0.5;
% % hold off
% % axis equal


clear; close all; clc;
data = load('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/WT4_WVA_AG0_WAVE.mat');
data = data.output;

%%

x = data.x;
colors = turbo(25);

wv1 = data.wave_profiles(1,:);
wv2 = data.wave_profiles(2,:);

figure()
hold on
plot(x, wv1, 'color', 'black')
plot(x, wv2, 'color', 'red')
hold off





























