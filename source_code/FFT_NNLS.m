function [X, Y, FFT_NNLS_result] = FFT_NNLS(CSM, g, w, frequencies, scan_limits, scan_resolution, tol, maxIter)
%
% This code implements the FFT-NNLS algorithm
%
% More information about FFT-NNLS can be found in the paper:
%    Ehrenfried, Klaus and Koop, Lars,
%    "Comparison of iterative deconvolution algorithms for the mapping of acoustic sources",
%    AIAA journal, 2007.
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
%    FFT_NNLS_result:  beamforming map, obtained by FFT-NNLS
%
% Author: Hao Liang 
% Last modified by: 23/07/28
%


% Scanning plane setting
X = scan_limits(1):scan_resolution:scan_limits(2);
Y = scan_limits(3):scan_resolution:scan_limits(4);
N_X = length(X); N_Y = length(Y); N_mic = size(g,3);

% Parameter initialization
FFT_NNLS_result = zeros(N_X, N_Y);
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
    
    % Start FFT-NNLS
    x = x0;

    % Precompute fft of PSF
    fft_psf = fft2(PSF_zeropad);
    
    % Function to calculate the residual vector and gradient
    nnls_handle = @(x) nnls(DAS_freqK_zeropad,x,fft_psf); 

    % Start iteration
    for n = 1:maxIter

        % In the first step the residual vector r_nnls is calculated, see eq.(27)
        [r_nnls, grad] = nnls_handle(x);  
        
        % Then the vector w_nnls = -grad is calculated, see eq.(28)
        w_nnls = -grad;                   
        w_nnls(x == 0 & w_nnls < 0) = 0;      
        
        % The auxiliary vector is obtained, see eq.(29)
        g_nnls = fftshift(ifft2(fft2(w_nnls).*fft_psf));   
        
        % An optimal step vector lambda is estimated, see eq.(23)
        lambda = -dot(g_nnls(:),r_nnls(:))/dot(g_nnls(:),g_nnls(:));    
        
        % Finally the new solution is calculated, see eq.(25)
        xold = x; x = max(0, x + lambda*w_nnls);  
        
        xDif = norm(x-xold)/norm(xold);
        if xDif < tol
            disp(['Converged after ' num2str(n) ' iterations'])
            break
        end

    end
    
    % Remove zero-padding
    x = x(int64(N_X/2)+1:int64(N_X/2 + N_X),int64(N_X/2)+1:int64(N_X/2 + N_X));
    
    % Accumulate the frequency components
    FFT_NNLS_result = FFT_NNLS_result + x;

end

end
