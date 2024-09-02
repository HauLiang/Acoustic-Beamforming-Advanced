function [X, Y, DAMAS2_result] = DAMAS2(CSM, g, w, frequencies, scan_limits, scan_resolution, tol, maxIter)
%
% This code implements the DAMAS2 algorithm
%
% More information about DAMAS2 can be found in the paper:
%    Dougherty, Robert, 
%    "Extensions of DAMAS and benefits and limitations of deconvolution in beamforming", 
%    11th AIAA/CEAS aeroacoustics conference, 2005.
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
%    DAMAS2_result:  beamforming map, obtained by DAMAS2
%
% Author: Hao Liang 
% Last modified by: 23/07/28
%


% Scanning plane setting
X = scan_limits(1):scan_resolution:scan_limits(2);
Y = scan_limits(3):scan_resolution:scan_limits(4);
N_X = length(X); N_Y = length(Y); N_mic = size(g,3);

% Parameter initialization
DAMAS2_result = zeros(N_X, N_Y);
N_freqs = length(frequencies);


% Start scan-frequency beamforming
for K = 1:N_freqs  
    
    % Calculate the steering vector and weighted steering vector corresponding to the frequency K
    wk = w(:,:,:,K); gk = g(:,:,:,K);
    wk_reshape = reshape(wk, [], N_mic).';

    % Calculate the DAS beamforming map corresponding the K frequency
    DAS_freqK = sum(conj(wk_reshape).*(CSM(:,:,K)*wk_reshape), 1)./N_mic^2;
    DAS_freqK = reshape(DAS_freqK, N_X, N_Y);
    
    % Calculate the steering vector corresponding to the source position in the centre of the region of interest
    X_centre_index = round(length(X)/2);
    Y_centre_index = round(length(Y)/2);
    g_centre = gk(X_centre_index, Y_centre_index, :); 
    
    % Calculate the shift-invariant point spread function (PSF)
    w_reshape = wk_reshape.'; g_reshape = reshape(g_centre, [], N_mic);
    PSF = (abs(conj(w_reshape)*g_reshape.').^2)./N_mic^2;
    PSF = reshape(PSF, N_X, N_Y);
    
    % To avoid wraparound effects, zero padding is used
    DAS_freqK_zeropad = real(zeropad(DAS_freqK));
    PSF_zeropad = zeropad(PSF);
    x0 = zeros(2*N_X);
    
    % Calculate the discrete integration constant
    a = sum(sum(PSF));
    
    % Design Gaussian regularization filter function
    k = 2*pi*frequencies(K);
    k_c = 1.2*k;
    psi = exp(-k^2/2/k_c^2);  
    
    % Start DAMAS2
    x = x0;
    
    % Precompute fft of PSF
    fft_psf = fft2(PSF_zeropad);
    
    % Start iteration
    for n = 1:maxIter

        % In the first step the residual vector r_nnls is calculated, see eq.(27)
        r = fftshift(ifft2(fft2(x).*fft_psf.*psi));
        
        % Finally the new solution is calculated, see eq.(25)
        xold = x; x = max(0, xold + (DAS_freqK_zeropad - r)/a);  
        
        xDif = norm(x-xold)/norm(xold);
        if xDif < tol
            disp(['Converged after ' num2str(n) ' iterations'])
            break
        end

    end
    
    % Remove zero-padding
    x = x(int64(N_X/2)+1:int64(N_X/2 + N_X),int64(N_X/2)+1:int64(N_X/2 + N_X));
    
    % Accumulate the frequency components
    DAMAS2_result = DAMAS2_result + x;
    
end

end
