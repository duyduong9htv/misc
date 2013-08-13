% step 2: 
% zoom, get out whale call on beam

sn = linspace(-1, 1, 64); %sn = sin(theta), scanning angle. -1 for - 90 degrees and 1 for 90 degrees 

zoom = axis; 
t_start = zoom(3);
t_end = zoom(4);
f_samp = 8e3; 
rcv_t = input_series(round(t_start*f_samp):round(t_end*f_samp), :); 
rcv_f = ifft(rcv_t)*(t_end-t_start); %Fourier transform. Note the factor t_end - t_start. 
% This has to be done such that the PArseval theorem is satisfied, because
% IFFT/FFT function in Matlab is not giving the correct total signal
% energy after transformed. 

 % bandpass filter and beamform the whale call 
f_i = zoom(1); f_f = zoom(2); 
T_duration = t_end - t_start; 

for ii = 1:max(size(rcv_t))
    if (ii-f_i*T_duration)*(ii-f_f*T_duration)>0
        rcv_f(ii,:) = zeros(1, 64); 
    end        
end

final_sum = beamform2(rcv_f, 'hanning', aperture); 
time_vector = linspace(t_start, t_end, size(rcv_t, 1));
figure; 
imagesc(sn,time_vector, 10*log10(abs(final_sum').^2));


