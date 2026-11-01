function y_fo = AP_fadeOut(x_in,tFo,Fs)
% Author: Abhiroop Paul
% Date: 4th Aug 2026
%
% Title: Fadeout Adder
% 
% Description:
% This is a function that applies a fadeout to the input audio. The
% function needs three arguments: the input audio specified as a vector
% "x_in", the fadeout time in seconds "tFo" and the sample rate "Fs".
% Please note this only works for mono audio.
% =======================================================================
N       =   length(x_in);                   % Num. samples in input audio

nFo     =   min(N,floor(Fs*tFo));           % Fadeout length in samples
fade    =   flip([0:nFo-1]')/nFo;           % Fadeout multiplier   

y_fo    =   x_in;
y_fo(N - nFo + 1:end)    =   x_in(N - nFo + 1:end) .* fade;    % Faded output

end