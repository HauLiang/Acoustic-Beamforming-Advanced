function mic_signal = simulateArraydata(source_info, mic_info, c, fs, duration, mic_centre)
%
% This code implements the generation of the simulated time-domain signal
%
%
% Inputs:
%    source_info:  the vector incorporating source information
%    mic_info:     the vector incorporating microphone array information
%    c:    speed of sound
%    fs:   sampling frequency
%    duration:     signal duration
%    mic_centre:   coordinates of the center of the microphone array
%    
% Outputs:
%    mic_signal:   simulated time-domain signal   
%
% Author: Hao Liang 
% Last modified by: 23/07/29
%


% Number of sources and microphone sensors
N_source = size(source_info, 1); 
N_mic = size(mic_info, 1); 

% Calculate the number of sample points 
t = 0:1/fs:(duration-1/fs);
N_samples = length(t);

% Initialize the simulated signal
mic_signal = zeros(N_mic, N_samples);

% Accumulate the signal for each source (here based on the assumption of incoherent sources)
for I = 1:N_source
    
    % Calculate the distance from the source to the center of the microphone array
    r_source_to_centre = norm(mic_centre-source_info(I, 1:3));
    
    % Calculate the amplitude of the signal (power = amp^2)
    % - Note that the formulation of sound pressure level (SPL) is SPL = 20*log10(amp/2e-5)
    % - Here the amplitude is calculated by the SPL of the source
    amp = 2e-5*10^(source_info(I, 5)/20);
    
    % Simulate the multi-channel signal of the I-th source
    for J = 1:N_mic

        % Calculate the distance from the J-th microphone to the simulated source 
        r_source_to_mic = sqrt(dot(mic_info(J, :) - source_info(I, 1:3), mic_info(J, :) - source_info(I, 1:3)) );
        
        % Simulate the signal collected by the J-th microphone
        delay_time = (r_source_to_mic-r_source_to_centre)/c;  % Calculate the delay time
        mic_signal(J, :) = mic_signal(J, :) + sqrt(2)*amp*cos(2*pi*source_info(I, 4)*(t-delay_time))*r_source_to_centre/r_source_to_mic;  
    
    end 
    
end    

end