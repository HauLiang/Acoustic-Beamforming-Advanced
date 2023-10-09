function L = lipschitz(PSF, fft_psf, D, lambda)

% Calculate Lipschitz constant by power iteration for FFT-DFISTA algorithm
xx = rand(size(PSF));
for k = 1:10
    xx = (fftshift(ifft2(fft2(xx).*fft_psf))+lambda*D*xx)/norm(xx,'fro');
end

L = norm(xx,'fro')^2;    

end
