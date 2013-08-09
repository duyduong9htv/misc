function [out_t, out_f] = BandPassFilter(signal_t, f_i, f_f, f_samp)
% function [out_t, out_f] = BandPassFilter(signal_t, f_i, f_f, f_samp)
% bandpasses the time-domain signal signal_t using the specified initial
% and final frequencies f_i and f_f. A rectangular 
% INPUTS        signal_t: time domain input signal
%                     f_i: initial frequency
%                     f_f: final frequency
%                  f_samp: sampling frequency, typically 8000 Hz 
% OUTPUTS       out_t, out_f: time and frequency domain output signals 

                
    
if size(signal_t, 1) < size(signal_t, 2)
    signal_t = transpose(signal_t);
end


T = size(signal_t, 1)/f_samp; 
signal_f = ifft(signal_t)*T; 
df = 1/T; 
for ff = 1:size(signal_f, 1)
    if (ff < f_i/df) || (ff > f_f/df)
        signal_f(ff, :) = zeros(1, size(signal_f, 2)); 
    end
end
out_f = signal_f; 
out_t = fft(signal_f)/T; 
end 


