% Title: Digital Filter Convolution
% 
% Author: Abhiroop Paul
% Date: 27th January, 2026
%
% This script takes two vectors as input: an impulse response vector "h"
% and an audio input vector "x". It convolves x with h and returns the
% output in a vector y.
%
% It does so by using the 'fast convolution' algorithm where it multiplies
% the DFT of the signal and impulse response rather than convolving
% directly.
% ______________________________________________________________________

function out = AP_myfastconv(h,x)

    %Number of taps N
    N = length(h) - 1;

    %Length of DFT L
    L = N + length(x) - 1;

    LFFT = 2^(ceil(log2(L)));

    % Compute the DFT of the input signal and impulse response
    X = fft(x, LFFT);
    H = fft(h, LFFT);


    %Multiplication (Convolution in frequency domain)
    Y = H .* X;
    
    %IFFT
    y = ifft(Y, LFFT, "symmetric");

    % Normalisation
    y = y / max(abs(y));

    % Output truncation to match input length
    out = y(1:length(x));

end

