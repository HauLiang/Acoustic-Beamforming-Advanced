function [X, Y, MUSIC_result] = MUSIC(CSM, w, frequencies, scan_limits, scan_resolution, nSources)
%
% This code implements the MUSIC algorithm
%
% More information about MUSIC can be found in the paper:
%    Schmidt, Ralph, 
%    "Multiple emitter location and signal parameter estimation", 
%    IEEE transactions on antennas and propagation, 1986.
%
%
% Inputs:
%    CSM:  cross-spectrum matrix (CSM)
%    w:    weighted steering vector
%    frequencies:   scan-frequency band
%    scan_limits:   scanning plane
%    scan_resolution:   scan resolution
%    nSources:   number of sources
%
% Outputs:
%    X & Y:  Two-dimensional coordinates corresponding to the beamforming map
%    MUSIC_result:  beamforming map, obtained by MUSIC
%
% Author: Hao Liang 
% Last modified by: 23/07/28
%


% Scanning plane setting
X = scan_limits(1):scan_resolution:scan_limits(2);
Y = scan_limits(3):scan_resolution:scan_limits(4);
N_X = length(X); N_Y = length(Y); N_mic = size(w,3);

% Parameter initialization
MUSIC_result = zeros(N_X, N_Y);
N_freqs = length(frequencies);


% Start scan-frequency beamforming
for K = 1:N_freqs

    % Calculate the weighted steering vector and CSM corresponding to the frequency K
    wk = w(:,:,:,K);
    CSM_k = CSM(:,:,K);
   
    % Cross spectral matrix with diagonal loading
%     CSM_k = CSM_k + trace(CSM_k)/(N_mic^2)*eye(N_mic, N_mic);  
    CSM_k = CSM_k/N_mic;

    % Eigenvectors of CSM
    [Vec, Val] = eig(CSM_k);                           
    [~, Seq] = sort(max(Val));
    
    % Noise eigenvectors
    Vn = Vec(:,Seq(1:end-nSources));    

    % Spatial spectrum (pseudo spectral) imaging
    MUSIC_result_K = zeros(N_X, N_Y); 
    for i = 1:N_X 
        for j = 1:N_Y    
            wk_reshape = reshape(wk(i,j,:), N_mic, 1);
            MUSIC_result_K(i, j) = 1./(wk_reshape'*(Vn*Vn')*wk_reshape);
        end
    end
       
    % Accumulate the frequency components
    MUSIC_result = MUSIC_result + MUSIC_result_K;

end

end