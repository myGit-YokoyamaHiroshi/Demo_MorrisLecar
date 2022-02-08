function [f, fftp] = calc_fft_spectrum(X, fs, nfft)
    if size(X,2)>size(X,1)
        X = X.';
    end
    timesteps = size(X, 1);
    srate = fs;
    
    if nargin < 3
        NFFT = 2^nextpow2(timesteps);
    else
        NFFT = nfft;
    end
    
    f = srate/2*linspace(0, 1, NFFT/2+1);
    fft_spectrum = fft(X, NFFT).';
    fftp = abs(fft_spectrum(:,1:NFFT/2+1)).^2;
end