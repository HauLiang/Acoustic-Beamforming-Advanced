function [X, Y, DAMAS_result] = DAMAS(CSM, g, w, frequencies, scan_limits, scan_resolution, deps, maxIter)
%
% This code implements the DAMAS algorithm
%
% More information about DAMAS can be found in the paper:
%    Brooks, Thomas F and Humphreys, William M, 
%    "A deconvolution approach for the mapping of acoustic sources (DAMAS) determined from phased microphone arrays", 
%    Journal of sound and vibration, 2006.
%
%
% Inputs:
%    CSM:  cross-spectrum matrix (CSM)
%    g:    steering vector
%    w:    weighted steering vector
%    frequencies:   scan-frequency band
%    scan_limits:   scanning plane
%    scan_resolution:   scan resolution
%    deps:    stopping threshold
%    maxIter: the maximum allowable iterations
%
% Outputs:
%    X & Y:  Two-dimensional coordinates corresponding to the beamforming map
%    DAMAS_result:  beamforming map, obtained by DAMAS
%
% Author: Hao Liang 
% Last modified by: 23/07/28
%


% Scanning plane setting
X = scan_limits(1):scan_resolution:scan_limits(2);
Y = scan_limits(3):scan_resolution:scan_limits(4);
N_X = length(X); N_Y = length(Y); N_mic = size(g,3);

% Parameter initialization
DAMAS_result = zeros(N_X, N_Y);
N_freqs = length(frequencies);


% Start scan-frequency beamforming
for K = 1:N_freqs  

    % Calculate the steering vector and weighted steering vector corresponding to the frequency K
    wk = w(:,:,:,K); gk = g(:,:,:,K);
    wk_reshape = reshape(wk, [], N_mic).';

    % Calculate the DAS beamforming map corresponding the K frequency
    DAS_freqK = sum(conj(wk_reshape).*(CSM(:,:,K)*wk_reshape), 1)./N_mic^2;
    DAS_freqK = reshape(DAS_freqK, N_X, N_Y);

    % Handle the dirty map by DAMAS
    temp = real(DAS_freqK);
    
    % Form the dictionary A
    w_reshape = wk_reshape.'; g_reshape = reshape(gk, [], N_mic);
    A = (abs(conj(w_reshape)*g_reshape.').^2)./N_mic^2;
    
    % Initialization
    Q = zeros(size(temp)); Q0 = temp;
    
    % Solve the inverse problem B_freqK = AQ for Q using Gauss-Seidel iteration
    % -- temp: DAS result (at frequency K)
    % -- A: reconstruction dictionary (at frequency K)
    % -- Q: DAMAS_result (at frequency K)
    for i = 1:maxIter

        for n = 1:N_X*N_Y
            Q(n) = max(0, temp(n) - A(n, 1:n-1)*Q(1:n-1)' ...
                - A(n, n+1:end)*Q0(n+1:end)');
        end
        
        dX = (Q - Q0);
        maxd = max(abs(dX(:)))/mean(Q0(:));
        
        % Stopping criteria
        if  maxd < deps
            disp(['Converged after ' num2str(i) ' iterations'])
            break;
        end
        
        Q0 = Q;

    end
    
    % Accumulate the frequency components
    DAMAS_result = DAMAS_result + Q;
    
end

end
