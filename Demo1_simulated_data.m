%% ----------- Acoustic Beamforming Demo for Simulated Data (Scan-frequency Version)
% ------ Support for setting multiple sources with different frequencies and SPL
% ------ Support for setting scanning planes with different distances
% ------ Support for selecting different microphone arrays
% ------ Support for setting the scanning frequency band
% ------ Support for setting the grid resolution
%
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
% The preliminary version can be found at: 
% -- https://github.com/HauLiang/Acoustic-Beamforming-Methods
%
% Author: Hao Liang 
% Last modified by: 23/08/03
%

close all; clear; clc;


%% Parameter settings for signal and beamforming

% Add path
addpath('.\Algorithm')
addpath('.\Toolbox')
addpath('.\Preprocess')
addpath('.\MicArray')

% Set scanning plane
scan_x = [-3 3];  % x-limit
scan_y = [-3 3];  % y-limit

% Set distance from scanning plane to microphone array plane
z_source = 3;

% Set speed of sound
c = 343; 

% Set sampling frequency
fs = 8000;

% Load coordinates of microphone sensors
load 56_spiral_array.mat
mic_x_axis = array(:,1); mic_y_axis = array(:,2); mic_z_axis = 0;
mic_pos = [mic_x_axis mic_y_axis]; mic_pos(:,3) = mic_z_axis;
mic_centre = mean(mic_pos);  % coordinates of the center of the microphone array

% Plot the microphone array
figure;
plot(mic_x_axis, mic_y_axis,'k.', 'MarkerSize',20);
xlim([min(mic_x_axis)-0.1, max(mic_x_axis)+0.1])
ylim([min(mic_y_axis)-0.1, max(mic_y_axis)+0.1])
set(gcf,'Position',[20 100 640 500]);	 
set(gca,'FontSize',24); set(gca,'linewidth',2); set(gcf,'Color','w');	

% Set coordinates of sources
source_x = [-1,0.5]';   % x-coordinates
source_y = [0,1]';      % y-coordinates

% Set the source frequency
source1_freq = 1950;  source2_freq = 2000;
sources_freq = [source1_freq, source2_freq]';

% Set signal duration
t_start = 0;  t_end = 0.02;
source_duration = t_end*ones(length(source_x),1);

% Set the sound pressure level (SPL) of sources
source1_spl = 100; source2_spl = 100;
sources_spl = [source1_spl, source2_spl].';

% Integrate source information to the vector "source_info"
% - x-coordinate / y-coordinate / z-coordinate / source frequency / SPL / source duration
source_info = [source_x, source_y, z_source*ones(length(source_x),1), ...
    sources_freq, sources_spl, source_duration];


%% Simulate the signal collected by the microphone array

mic_signal = simulateArraydata(source_info, mic_pos, c, fs, source_duration, mic_centre);


%% Set the scanning frequency band and calculate the cross spectrum matrix (CSM)

% Set the scanning frequency band (Hz)
search_freql = 1950;   % lowest scanning frequency
search_frequ = 2050;   % upper scanning frequency

% Calculate the CSM
[CSM, freqs] = developCSM(mic_signal.', search_freql, search_frequ, fs, t_start, t_end);


%% Calculate the steering vector

scan_resolution = 0.1;  % scan resolution setting
[g, w] = steerVector2(z_source, freqs, [scan_x scan_y], scan_resolution, mic_pos.', c, mic_centre);


%% Acoustic beamforming procedure

tic;

% Choose one algorithm to perform acoustic imaging:
% - Just delete the symbol "%" in front of the corresponding algorithm

% Implementation of DAS algorithm
[X, Y, B] = DAS(CSM, w, freqs, [scan_x scan_y], scan_resolution);

% Implementation of Fast-DAS algorithm
% [X, Y, B] = Fast_DAS(CSM, w, freqs, [scan_x scan_y], scan_resolution);

% Implementation of MUSIC algorithm
% nSources = 2;     % number of sources
% [X, Y, B] = MUSIC(CSM, w, freqs, [scan_x scan_y], scan_resolution, nSources);

% Implementation of DAMAS algorithm
% deps = 0.1;       % tolerance of convergence criterion
% maxIter = 1000;   % number of iterations
% [X, Y, B] = DAMAS(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, deps, maxIter);

% Implementation of DAMAS2 algorithm
% tol = 1e-4;       % tolerance of convergence criterion
% maxIter = 1000;   % number of iterations
% [X, Y, B] = DAMAS2(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, tol, maxIter);

% Implementation of CLEAN-PSF algorithm
% loopgain = 0.99;  % damping parameter
% maxIter = 50;     % number of iterations
% [X, Y, B] = CLEAN_PSF(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, loopgain, maxIter);

% Implementation of CLEAN-SC algorithm
% loopgain = 0.99;  % damping parameter
% maxIter = 50;     % number of iterations
% [X, Y, B] = CLEAN_SC(CSM, w, freqs, [scan_x scan_y], scan_resolution, loopgain, maxIter);

% Implementation of FFT-NNLS algorithm
% tol = 1e-4;       % tolerance of convergence criterion
% maxIter = 1000;   % number of iterations
% [X, Y, B] = FFT_NNLS(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, tol, maxIter);

% Implementation of FFT-FISTA algorithm
% tol = 1e-4;       % tolerance of convergence criterion
% maxIter = 1000;   % number of iterations
% [X, Y, B] = FFT_FISTA(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, tol, maxIter);

% Implementation of DAMAS-FISTA algorithm
% tol = 1e-4;       % tolerance of convergence criterion
% maxIter = 1000;   % number of iterations
% [X, Y, B] = DAMAS_FISTA(CSM, g, w, freqs, [scan_x scan_y], scan_resolution, tol, maxIter);
% 

toc;


%% Numerical result display 

% Converting source power to sound pressure level (SPL)
B(B<0)=0;    % delete negative power
SPL = 20*log10((eps+sqrt(real(B)))/2e-5);

% Plot the beamforming map
figure;
BF_dr = 6; maxSPL = ceil(max(SPL(:)));
contourf(X, Y, SPL, (maxSPL-BF_dr):1:maxSPL,'LineStyle','none'); colorbar; clim([maxSPL-BF_dr maxSPL]);
hold on; plot(source_x(:),source_y(:),'r*');   % true location
set(gcf,'Position',[20 100 640 500]);	 
set(gca,'FontSize',24); set(gca,'linewidth',2); set(gcf,'Color','w');	
xlim([-3 3]); ylim([-3 3]);

x1 = source_x(1)-0.2; y1 = source_y(1)-0.2; 
x2 = source_x(1)+0.2; y2 = source_y(1)+0.2;
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','k','Linewidth',1);

xx1 = source_x(2)-0.2; yy1 = source_y(2)-0.2; 
xx2 = source_x(2)+0.2; yy2 = source_y(2)+0.2;
hold on; rectangle('Position',[xx1 yy1 xx2-xx1 yy2-yy1],'EdgeColor','k','Linewidth',1);

h1 = axes('position',[0.2 0.15 0.25 0.25]); 
axis(h1);
contourf(X, Y, SPL, (maxSPL-BF_dr):1:maxSPL,'LineStyle','none'); 
hold on; plot(source_x(:),source_y(:),'r*'); 
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'fontsize',12)


h2 = axes('position',[0.5 0.15 0.25 0.25]);
axis(h2);
contourf(X, Y, SPL, (maxSPL-BF_dr):1:maxSPL,'LineStyle','none'); 
hold on; plot(source_x(:),source_y(:),'r*'); 
xlim([xx1 xx2]);ylim([yy1 yy2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'fontsize',12)

