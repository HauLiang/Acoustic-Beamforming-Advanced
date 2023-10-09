function [X, Y, B] = Acoustic_data_process(data, steering_type, beamform_type, fs, z_source, scan_freq, scan_x, scan_y, scan_resolution)
%
% This code implements the acoustic beamforming demo for loading data
%
%
% Inputs:
%    data: .h5 file or others, signal collected by mircophone array or simulated by software
%    steering_type:  type of steering vecotor (1, 2, 3, or 4), see './Preprocess/steerVector.m'
%    beamform_type:  type of beamforming algorithm, see './Algorithm'
%    fs:   ampling frequency
%    z_source:    distance from scanning plane to microphone array plane
%    scan_freq:   scanning frequency band
%    scan_x:   range of scanning plane (x-axis)
%    scan_y:   range of scanning plane (y-axis)
%    scan_resolution:   scan resolution
%
% Outputs:
%    X & Y:  Two-dimensional coordinates corresponding to the beamforming map
%    B:  beamforming map
%
% Author: Hao Liang 
% Last modified by: 23/08/03
%


%% Preliminary

% Add path
addpath('.\Algorithm')
addpath('.\Toolbox')
addpath('.\Preprocess')
addpath('.\MicArray')

% Set speed of sound
c = 343; 

% Load coordinates of microphone sensors
load 56_spiral_array.mat
mic_x_axis = array(:,1); mic_y_axis = array(:,2); mic_z_axis = 0;
mic_pos = [mic_x_axis mic_y_axis]; mic_pos(:,3) = mic_z_axis;
mic_centre = mean(mic_pos);  % coordinates of the center of the microphone array


%% Calculate cross-spectrum matrix (CSM)

% Determine the scanning frequency band (Hz)
search_freql = scan_freq(1); search_frequ = scan_freq(end);
[CSM, freqs] = developCSM(data.', search_freql, search_frequ, fs);


%% Calculate the steering vector

if steering_type == 1
    [g, w] = steerVector1(z_source, freqs, [scan_x scan_y], scan_resolution, mic_pos.', c, mic_centre);
elseif steering_type == 2
    [g, w] = steerVector2(z_source, freqs, [scan_x scan_y], scan_resolution, mic_pos.', c, mic_centre);
elseif steering_type == 3
    [g, w] = steerVector3(z_source, freqs, [scan_x scan_y], scan_resolution, mic_pos.', c, mic_centre);
elseif steering_type == 4
    [g, w] = steerVector4(z_source, freqs, [scan_x scan_y], scan_resolution, mic_pos.', c, mic_centre);
else
    error('There are only four types of steering vectors!')
end


%% Acoustic beamforming procedure

if (strcmp(beamform_type, 'DAS'))
    % Implementation of DAS algorithm
    [X, Y, B] = DAS(CSM, w, freqs, [scan_x scan_y], scan_resolution);
  
elseif (strcmp(beamform_type, 'Fast-DAS'))
    % Implementation of Fast-DAS algorithm
    [X, Y, B] = Fast_DAS(CSM, w, freqs, [scan_x scan_y], scan_resolution);
    
elseif (strcmp(beamform_type, 'MUSIC'))
    % Implementation of MUSIC algorithm
    nSources = 2;      % number of sources
    [X, Y, B] = MUSIC(CSM, w, freqs, [scan_x scan_y], scan_resolution, nSources);
    
elseif (strcmp(beamform_type, 'DAMAS'))
    % Implementation of DAMAS algorithm
    deps = 0.1;        % tolerance of convergence criterion
    maxIter = 1000;    % number of iterations
    [X, Y, B] = DAMAS(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, deps, maxIter);

elseif (strcmp(beamform_type, 'DAMAS2'))
    % Implementation of DAMAS2 algorithm
    tol = 1e-4;        % tolerance of convergence criterion
    maxIter = 1000;    % number of iterations
    [X, Y, B] = DAMAS2(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, tol, maxIter);

elseif (strcmp(beamform_type, 'CLEAN-PSF'))
    % Implementation of CLEAN-PSF algorithm
    loopgain = 0.99;   % damping parameter
    maxIter = 50;      % number of iterations
    [X, Y, B] = CLEAN_PSF(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, loopgain, maxIter);
    
elseif (strcmp(beamform_type, 'CLEAN-SC'))
    % Implementation of CLEAN-SC algorithm
    loopgain = 0.99;   % damping parameter
    maxIter = 50;      % number of iterations
    [X, Y, B] = CLEAN_SC(CSM, w, freqs, [scan_x scan_y], scan_resolution, loopgain, maxIter);
    
elseif (strcmp(beamform_type, 'FFT-NNLS'))
    % Implementation of FFT-NNLS algorithm
    tol = 1e-4;        % tolerance of convergence criterion
    maxIter = 1000;    % number of iterations
    [X, Y, B] = FFT_NNLS(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, tol, maxIter);
    
elseif (strcmp(beamform_type, 'FFT-DFISTA'))
    % Implementation of FFT-DFISTA algorithm
    lambda = 1e-3;    % penalty parameter
    tol = 1e-4;       % tolerance of convergence criterion
    maxIter = 1000;   % number of iterations
    [X, Y, B] = FFT_DFISTA(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, lambda, tol, maxIter);

elseif (strcmp(beamform_type, 'DAMAS-FISTA'))
    % Implementation of DAMAS-FISTA algorithm
    tol = 1e-4;        % tolerance of convergence criterion
    maxIter = 1000;    % number of iterations
    [X, Y, B] = DAMAS_FISTA(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, tol, maxIter);

else
    error('This algorithm does not exist!')
    
end

end

