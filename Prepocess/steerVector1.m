function [g, w] = steerVector1(plane_distance, frequencies, scan_limits, scan_resolution, mic_positions, c, mic_centre)
%
% This code implements the generation of the steering vector -- Formulation I
%
% More information about the steering vector can be found in the paper:
%    Sarradj, Ennes, 
%    "Three-dimensional acoustic source mapping with different beamforming steering vector formulations", 
%    Advances in Acoustics and Vibration, 2012.
%
%
% Inputs:
%    plane_distance:    distance from scanning plane to microphone array plane
%    frequencies:   scan-frequency band
%    scan_limits:   scanning plane
%    scan_resolution:   scan resolution
%    mic_positions:     coordinates of the microphone array
%    c:    speed of sound
%    mic_centre:    coordinates of the center of the microphone array
%    
% Outputs:
%    g:    steering vector
%    w:    weighted steering vector
%
% Author: Hao Liang 
% Last modified by: 23/07/29
%


% Number of microphones and the number of scanning frequency points
N_mic = size(mic_positions, 2);
N_freqs = length(frequencies);

% Define the scanning plane
x = scan_limits(1):scan_resolution:scan_limits(2); 
y = scan_limits(3):scan_resolution:scan_limits(4); 
z = plane_distance;
N_X = length(x); N_Y = length(y);
X = repmat(x,N_X,1); Y = repmat(y.',1,N_Y);

% Initialize the steering vector and the weighted steering vector
g = zeros(N_X, N_Y, N_mic, N_freqs);
w = zeros(N_X, N_Y, N_mic, N_freqs);

% Calculate the distance from the scanning plane to the center of the microphone array
r_scan_to_mic_centre = sqrt((X-mic_centre(1)).^2 + (Y-mic_centre(2)).^2 + (z-mic_centre(3))^2);  

% Initialize the parameter
r_scan_to_mic = zeros(N_X, N_Y, N_mic);

% Calculate the steering vector and the weighted steering vector
for K = 1:N_freqs

    % Angular frequency, w
    omega = 2*pi*frequencies(K);  

    for m = 1:N_mic
        
        % Calculate the distance from the scan plane to the m-th microphone
        r_scan_to_mic(:,:,m) = sqrt((X-mic_positions(1,m)).^2+(Y-mic_positions(2,m)).^2 + z^2); 
        
        % Calculate the weighted steering vecor (formulation I)
        w(:,:,m, K) = exp(-1j*omega.*(r_scan_to_mic(:,:,m)-r_scan_to_mic_centre)./c);

        % Calculate the steering vector (formulation I)
        g(:,:,m, K) = exp(-1j*omega.*(r_scan_to_mic(:,:,m)-r_scan_to_mic_centre)./c);

    end

end

end