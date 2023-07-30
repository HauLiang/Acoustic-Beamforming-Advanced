function L = lipschitz(PSF, fft_psf)

% Estimate Lipschitz constant by power iteration
xx = rand(size(PSF));
for k = 1:10
    xx = fftshift(ifft2(fft2(xx).*fft_psf))/norm(xx,'fro');
end

L = norm(xx,'fro')^2;    % lambda(A'A) Assuming a symmetric matric A

end
