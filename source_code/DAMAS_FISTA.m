function [X, Y, DAMAS_FISTA_result] = DAMAS_FISTA(CSM, g, w, frequencies, scan_limits, scan_resolution, tol, maxIter)
%
% This code implements the DAMAS-FISTA algorithm
%
% More information about DAMAS-FISTA can be found in the paper:
%    Liang, Hao and Zhou, Guanxing and Tu, Xiaotong and Jakobsson, Andreas and Ding, Xinghao and Huang, Yue,
%    "Learning an Interpretable End-to-End Network for Real-Time Acoustic Beamforming", 
%    arXiv:2306.10772, 2023.
%
%
% Inputs:
%    CSM:  cross-spectrum matrix (CSM)
%    g:    steering vector
%    w:    weighted steering vector
%    frequencies:   scan-frequency band
%    scan_limits:   scanning plane
%    scan_resolution:   scan resolution
%    tol:    stopping threshold
%    maxIter: the maximum allowable iterations
%
% Outputs:
%    X & Y:  Two-dimensional coordinates corresponding to the beamforming map
%    DAMAS_FISTA_result:  beamforming map, obtained by DAMAS-FISTA
%
% Author: Hao Liang 
% Last modified by: 23/07/28
%


% Scanning plane setting
X = scan_limits(1):scan_resolution:scan_limits(2);
Y = scan_limits(3):scan_resolution:scan_limits(4);
N_X = length(X); N_Y = length(Y); N_mic = size(g,3);

% Parameter initialization
DAMAS_FISTA_result = zeros(N_X*N_Y, 1);
N_freqs = length(frequencies);


% Start scan-frequency beamforming
for K = 1:N_freqs  
    
    % Calculate the steering vector and weighted steering vector corresponding to the frequency K
    wk = w(:,:,:,K); gk = g(:,:,:,K);
    wk_reshape = reshape(wk, [], N_mic).';

    % Calculate the DAS beamforming map corresponding the K frequency
    DAS_freqK = sum(conj(wk_reshape).*(CSM(:,:,K)*wk_reshape), 1)./N_mic^2;

    % Handle the dirty map by DAMAS-FISTA
    temp = real(DAS_freqK);
    
    % Form the dictionary A
    w_reshape = wk_reshape.'; g_reshape = reshape(gk, [], N_mic);
    A = (abs(conj(w_reshape)*g_reshape.').^2)./N_mic^2;
    
    % Initialization
    x0 = zeros(size(temp.')); b = temp.';
    
    % Initialize variables
    x = x0; xold = x;
    y = x; t = 1;       

    % Precompute ATA and ATb
    ATA = A.'*A; ATb = A.'*b;
    
    % Compute Lipschitz constant
    L = max(eig(ATA));

    % Calculate gradient
    grad = ATA*y - ATb;
    
    % Start iteration
    for n = 1:maxIter   

        % Update x
        x = max(0,y - (1/L)*grad);
        
        xDif = norm(x-xold)/norm(xold);
        if xDif < tol
            disp(['Converged after ' num2str(n) ' iterations'])
            break
        end
        
        % Update step
        tnew = (1+sqrt(1+4*t*t))/2;
        
        % Update y
        y = x + ((t-1)/tnew)*(x-xold);
        
        % Re-calculate gradient
        grad = ATA*y - ATb;
        xold = x;
        t = tnew;  

    end

    % Accumulate the frequency components
    DAMAS_FISTA_result = DAMAS_FISTA_result + x;

end

% Reshape to size of N_X*N_Y
DAMAS_FISTA_result = reshape(DAMAS_FISTA_result, N_X, N_Y);

end
