function [X, Y, CLEAN_SC_result] = CLEAN_SC(CSM, w, frequencies, scan_limits, scan_resolution, loopgain, maxIter)
%
% This code implements the CLEAN-SC algorithm
%
% More information about CLEAN-SC can be found in the paper:
%    Sijtsma, Pieter, 
%    "CLEAN based on spatial source coherence", 
%    International journal of aeroacoustics, 2007.
%
%
% Inputs:
%    CSM:  cross-spectrum matrix (CSM)
%    w:    weighted steering vector
%    frequencies:   scan-frequency band
%    scan_limits:   scanning plane
%    scan_resolution:   scan resolution
%    loopgain:  loop gain
%    maxIter:   the maximum allowable iterations
%
% Outputs:
%    X & Y:  Two-dimensional coordinates corresponding to the beamforming map
%    CLEAN_SC_result:  beamforming map, obtained by CLEAN-SC
%
% Author: Hao Liang 
% Last modified by: 23/07/28
%


% Scanning plane setting
X = scan_limits(1):scan_resolution:scan_limits(2);
Y = scan_limits(3):scan_resolution:scan_limits(4);
N_X = length(X); N_Y = length(Y); N_mic = size(w,3);

% Parameter initialization       
CLEAN_SC_result = zeros(N_X, N_Y);
N_scanpoints = N_X*N_Y; N_freqs = length(frequencies);


% Start scan-frequency beamforming
for K = 1:N_freqs 
    
    % Calculate the weighted steering vector corresponding to the frequency K
    wk = w(:,:,:,K); wk_reshape = reshape(wk,[], N_mic).'./N_mic;
    
    % Start CLEAN-SC procedure
    Clean_map = zeros(1, N_scanpoints); Degraded_CSM = CSM(:,:,K); Dirty_map = sum(conj(wk_reshape).*(Degraded_CSM*wk_reshape), 1);
    Dcurr = sum(abs(Degraded_CSM(:))); count = 0; Dprev = 1e10;

    while ( Dcurr < Dprev ) && (count < maxIter)

        % Determine peak source
        [Map_max, index_max] = max(abs(Dirty_map)); 
        Map_maxCleanBeam = zeros(size(Clean_map)); Map_maxCleanBeam(index_max) = 1;
        ispositiv = real(Dirty_map(index_max)) > 0;  % Delete negative pressure maps
        hmax = Degraded_CSM*wk_reshape(:,index_max)/Map_max;

        % Update degraded CSM
        Cource = Map_max*(hmax)*hmax';
        Degraded_CSM = Degraded_CSM - loopgain*Cource;

        % Update dirty map 
        Dirty_map = sum(conj(wk_reshape).*(Degraded_CSM*wk_reshape), 1);

        % Update clean map 
        Clean_map = Clean_map + loopgain*Map_max*Map_maxCleanBeam*ispositiv;

        % Stopping criteria
        Dprev = Dcurr; Dcurr = sum(abs(Degraded_CSM(:)));
        count = count + 1;

    end
    
    if count == maxIter
        disp(['Stopped after maximum iterations (' num2str(maxIter) ')'])
    else
        disp(['Converged after ' num2str(count) ' iterations'])
    end

    % Delete negative power
    Dirty_map(real(Dirty_map)<0) = 0;

    % Final beamforming map equal to the sum of clean map and dirty map
    CLEAN_SC_map = Clean_map + Dirty_map;
    CLEAN_SC_map = reshape(CLEAN_SC_map, N_X, N_Y);
    
    % Accumulate the frequency components
    CLEAN_SC_result = CLEAN_SC_result + CLEAN_SC_map;
    
end

end
