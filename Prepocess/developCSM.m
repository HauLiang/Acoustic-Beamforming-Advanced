function [CSM, freqs] = developCSM(mic_signal, search_freql, search_frequ, fs, t_start, t_end)
%
% This code implements the generation of the cross-spectrum matrix (CSM)
%
%
% Inputs:
%    mic_signal:     time-domain signal collected by the microphone array
%    search_freql:   lowest scanning frequency
%    search_frequ:   upper scanning frequency
%    fs:    sampling frequency
%    t_start:   signal start time  
%    t_end:     signal termination time
%    
% Outputs:
%    CSM:     cross-spectrum matrix (CSM)
%    freqs:   scan-frequency band 
%
% Author: Hao Liang 
% Last modified by: 23/07/29
%


if nargin < 6
    t_start = 0;
    t_end = size(mic_signal, 1)/fs;
end

% Number of microphone sensors
N_mic = size(mic_signal, 2);

% Calculate the starting and end sample points
start_sample = floor(t_start*fs) + 1;
end_sample = ceil(t_end*fs);  

% Select the frequency points in the scanning frequency band
x_fr = fs / end_sample * (0:floor(end_sample/2)-1);
freq_sels = find((x_fr>=search_freql).*(x_fr<=search_frequ));

% Number of scanning frequency points
N_freqs = length(freq_sels);

% Initialize the cross-spectrum matrix (CSM)
CSM = zeros(N_mic, N_mic, N_freqs);

% Perform Fourier transform
mic_signal_fft = sqrt(2)*fft(mic_signal(start_sample:end_sample,:))/(end_sample-start_sample+1);

% Develop CSM
for K = 1:N_freqs

    % Calculate the CSM corresponding to the frequency K
    CSM(:,:,K) = mic_signal_fft(freq_sels(K),:).'*conj(mic_signal_fft(freq_sels(K),:));
    % CSM(:,:,K) = CSM(:,:,K) - diag(diag(CSM(:,:,K)));  %  Diagonal removal process

end
    
% Frequency points to be scanned
freqs = x_fr(freq_sels);

end