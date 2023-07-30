function [r, grad] = nnls(b, x, fft_psf)
% Nonnegative least squares objective function
%
% Input:
%   b : beamforming map 
%   x: 	beamforming map after deconvolution
%   fft_PSF:  point spread function (PSF) after Fourier transform
%
% Output: 
%   r:  residual  
%   g:  gradient 
%

% Calculate residual 
r = fftshift(ifft2(fft2(x).*fft_psf)) - b;

% Calculate gradient
grad = fftshift(ifft2(fft2(r).*fft_psf));


end