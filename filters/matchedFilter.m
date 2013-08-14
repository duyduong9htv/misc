function [MFt, MFf] = matchedFilter(signal_f,  Hf, fs) 
% function [MFt, MFf] = matchedFilter(signal_f,  Hf, fs) 
% calculates the (complex) matched filter output in time and frequency for
% an input signal in the frequency domain, when matched with a system
% function H. 
% INPUTS    :-signal_f, input signal in frequency domain, obtained by
%             taking the FFT (in our acoustic software suite: IFFT) of the
%             input time series. 
%            -Hf: Complex frequency spectrum (or IFFT/FFT) of the matched
%            signal. Must be of the same length as signal input 
%            -fs: sampling frequency. 
% OUTPUTS   :-MFt, MFf: time and frequency representation of the
%               matched-filter output. 
% Example: 
% fs = 8000; T = 4; t_p = 1; t = 1/fs*[1:1:fs*T]; 
% [st1, Sf1, rst1, rSf1] = getLFM(fs, 1, t_p, T, 390, 440);
% [Ht, Hf, ~, ~] = getLFM(8000, 0, t_p, T, 390, 440); 
% [MFt, MFf] = matchedFilter(Sf1, Hf, 8000); 
% figure; subplot(2, 1, 1); plot(t, real(st1)); title('original signal'); 
% subplot(2, 1, 2); plot(t, real(MFt)); title('matched-filtered signal showing peak at arrival time of original signal');
% xlabel('Time (seconds)'); 
%
% Written by DD Tran. Last updated Aug 9, 2013. 

if nargin > 2 
    T = length(Hf)/fs; %signal duration; 
else 
    T = length(Hf);  
end 

MFf = signal_f.*conj(Hf); 
MFt = fft(MFf)/T; %divided by T for proper energy normalization and satisfying the PArseval's theorem. 

end 
