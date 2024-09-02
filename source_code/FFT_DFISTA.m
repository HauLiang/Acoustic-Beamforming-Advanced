function [X, Y, FFT_DFISTA_result] = FFT_DFISTA(CSM, g, w, frequencies, scan_limits, scan_resolution, lambda, tol, maxIter)
%
% This code implements the FFT-DFISTA algorithm
%
% More information about FFT-NNLS can be found in the paper:
%    Ding, Xinghao and Liang, Hao and Jakobsson, Andreas and Tu, Xiaotong and Huang, Yue,
%    "High-Resolution Source Localization Exploiting the Sparsity of the Beamforming Map",
%    Signal Processing, 2022.
%
%
% Inputs:
%    CSM:  cross-spectrum matrix (CSM)
%    g:    steering vector
%    w:    weighted steering vector
%    frequencies:   scan-frequency band
%    scan_limits:   scanning plane
%    scan_resolution:   scan resolution
%    lambda: penalty parameter
%    tol:    stopping threshold
%    maxIter: the maximum allowable iterations
%
% Outputs:
%    X & Y:  Two-dimensional coordinates corresponding to the beamforming map
%    FFT_DFISTA_result:  beamforming map, obtained by FFT-NNLS
%
% Author: Hao Liang 
% Last modified by: 23/10/09
%


% Scanning plane setting
X = scan_limits(1):scan_resolution:scan_limits(2);
Y = scan_limits(3):scan_resolution:scan_limits(4);
N_X = length(X); N_Y = length(Y); N_mic = size(g,3);

% Parameter initialization
FFT_DFISTA_result = zeros(N_X, N_Y);
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
    
    % Construct the modified first-order difference matrix using sparse matrix
    D = zeros(2*size(PSF));
    for i = size(D,1)-1
        D(i,i) = 1; D(i,i+1) = -1;
    end
    D(size(D,1),size(D,1)) = 1;
    D = sparse(D);

    % To avoid wraparound effects, zero padding is used
    DAS_freqK_zeropad = real(zeropad(DAS_freqK));
    PSF_zeropad = zeropad(PSF);
    x0 = zeros(2*N_X);
    
    % Start FFT-DFISTA
    x = x0; xold = x;
    y = x; t = 1; 

    % Precompute fft of PSF
    fft_psf = fft2(PSF_zeropad);
    
    % Compute Lipschitz constant
    L = dfista_lipschitz(PSF_zeropad, fft_psf, D, lambda);

    % For n = 1
    gradient = dfista_gradient(fft_psf, DAS_freqK_zeropad, y, D, lambda);

    
    % Start iteration
    for n = 1:maxIter   
        % Step1. Update x
        x = max(0,y - (1/L)*gradient);
        
        % Convergence judgment
        xDif = norm(x-xold)/norm(xold);
        if xDif < tol
            disp(['Converged after ' num2str(n) ' iterations'])
            break
        end
        
        % Step2. Update t
        tnew = (1+sqrt(1+4*t*t))/2;
        
        % UStep3. Update y
        y = x + ((t-1)/tnew)*(x-xold);
        
        % Calculate gradient
        gradient = dfista_gradient(fft_psf, DAS_freqK_zeropad, y, D, lambda);
        xold = x;
        t = tnew;  
    end
    
    % Remove zero-padding
    x = x(int64(N_X/2)+1:int64(N_X/2 + N_X),int64(N_X/2)+1:int64(N_X/2 + N_X));
    
    % Accumulate the frequency components
    FFT_DFISTA_result = FFT_DFISTA_result + x;

end

end


