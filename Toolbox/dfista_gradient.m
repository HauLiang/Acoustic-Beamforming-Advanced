function dfista_gradient = dfista_gradient(fft_psf,b,x,D,lambda)
% Compute gradient for FFT-DFISTA algorithm
%
% Input:
%   fft_psf: point spread function (PSF) after Fourier transform
%   b:   beamforming map 
%   x:   beamforming map after deconvolution
%   D:   modified first-order difference matrix
%   lambda: penalty parameter
%
% Output: 
%   L:  estimated Lipschitz constant
%

% Residual vector
r = fftshift(ifft2(fft2(x).*fft_psf)) - b;


% Calculate gradient
dfista_gradient = fftshift(ifft2(fft2(r).*fft_psf))+lambda*(D.'*D)*x;

end