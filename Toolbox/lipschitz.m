function L = lipschitz(PSF, fft_psf, D, lambda)
% Calculate Lipschitz constant for FFT-DFISTA algorithm
%
% Input:
%   PSF :  point spread function
%   fft_psf: point spread function (PSF) after Fourier transform
%   D:     modified first-order difference matrix
%   lambda: penalty parameter
%
% Output: 
%   L:  estimated Lipschitz constant
%

% Using power iteration
xx = rand(size(PSF));
for k = 1:10
    xx = (fftshift(ifft2(fft2(xx).*fft_psf))+lambda*D*xx)/norm(xx,'fro');
end

L = norm(xx,'fro')^2;    

end
