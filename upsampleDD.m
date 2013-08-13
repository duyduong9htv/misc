function out = upsampleDD(s, fs, fu)
%function out = upsampleDD(s, fs, fu)
%upsamples a signal s, sampled at frequency fs to the upsampling frequency
%of fu (Hz) by interpolating between all measurement instances. 

t = [1:1:length(s)]/fs; 
tu = linspace(t(1), t(end), fu*(t(end) - t(1))); 
out = interp1(t, s, tu); 


end 