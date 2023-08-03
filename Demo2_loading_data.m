%% ----------- Acoustic Beamforming Demo for Loading Data (Scan-frequency Version)
%
% Input:
%   - data: signal collected by mircophone array or simulated by software
%   - fs: sampling frequency
%   - z_source: distance from scanning plane to microphone array plane
%   - scan_freq: scanning frequency band
%   - scan_x: range of scanning plane (x-axis)
%   - scan_y: range of scanning plane (y-axis)
%   - scan_resolution: scan plane resolution
%
% This demo includes the following methods:
%
% -- DAS
% -- MUSIC
% -- DAMAS
% -- DAMAS2
% -- DAMAS-FISTA
% -- CLEAN-PSF
% -- CLEAN-SC
% -- FFT-NNLS 
% -- FFT-FISTA
%
% Author: Hao Liang 
% Last modified by: 23/08/03
%


clear; clc; close all;

% Parameter setting
fs = 51200;     % sampling frequency
z_source = [4.94 3.82 3.36];        % see './data_label/label.txt'
scan_freq = [1950 2050];            % the range of scanning frequency
scan_x = [-2,2]; scan_y = [-2,2];   % the range of scanning plane 
scan_resolution = 0.1;              % resolution of scanning plane

% Add data path
source_path = './data';
source_dir = dir(source_path);

% Select the type of steering vector and beamforming method
steering_type = 2;            % see'./Preprocess/steerVector.m'
beamform_type = 'Fast-DAS';   % see './Algorithm'

% Process loading data
for i = 1:length(source_dir)-2
    
    % Load data
    process_data_file = source_dir(i+2).name;
    disp(['Processing data ' process_data_file])
    data = h5read([source_path,'/', source_dir(i+2).name],'/time_data');
        
    % Acoustic beamforming procedure
    z_dist = z_source(i);
    [X, Y, B] = Acoustic_data_process(data, steering_type, beamform_type, fs, z_dist, scan_freq, scan_x, scan_y, scan_resolution);
    
    % Convert source power to sound pressure level (SPL) 
    B(B<0)=0;
    SPL = 20*log10((eps+sqrt(real(B)))/2e-5);

    % Plot beamforming map, and see './data_label/label.txt' for ground truth
    figure; BF_dr = 6; maxSPL = ceil(max(SPL(:)));
    contourf(X, Y, SPL, (maxSPL-BF_dr):1:maxSPL); colorbar; clim([maxSPL-BF_dr maxSPL])
    set(gcf,'Position',[20 100 640 500]);	 
    set(gca,'FontSize',24); set(gca,'linewidth',2); set(gcf,'Color','w');	
    xlim([-2 2]); ylim([-2 2]);

end

