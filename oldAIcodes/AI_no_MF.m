% this code aims to test the array invariant (or waveguide invariant for a
% horizontal line array) when no matched-filtering is employed. Results are
% from data collected during the GOM 2006 experiment. 

%% Initialization 
close all; clear; 
pingname = 'fora2006jd276t002500.dat'; 
aperture1 = 'lf'; %adjust accordingly 


%% data read 
input_series = read_input_series(pingname, aperture1);
t1 = 1; 
t2 = 15; 


%% spectrogram
P1 = 0; 
for ii = 1:1:64
    disp(ii)
    [S, F, T, P] = spectrogram(input_series(round(8000*linspace(t1, t2, round((t2-t1)*8000))), ii), 512, 256, 8192, 8000); 
    P1 = P1+P; 
end
% select whale call of interest 
figure; imagesc(F, T+t1, 10*log10(P1')); 
setFont(18); setFigureAuto; 
% caxis([40 90])
title(pingname)
pause; 

%% beamform 

sn = linspace(-1, 1, 401); %sn = sin(theta), scanning angle. -1 for - 90 degrees and 1 for 90 degrees 
zoom = axis; 
t_start = zoom(3);
t_end = zoom(4);
f_samp = 8e3; 
rcv_t = input_series(round(t_start*f_samp):round(t_end*f_samp), :); 
rcv_f = ifft(rcv_t)*(t_end-t_start); %Fourier transform. Note the factor t_end - t_start. 
% This has to be done such that the PArseval theorem is satisfied, because
% IFFT/FFT function in Matlab is not giving the correct total signal
% energy after transformed. 
% bandpass filter and beamform the signal�i
f_i = zoom(1); f_f = zoom(2); 
T_duration = t_end - t_start; 

for ii = 1:max(size(rcv_t))
    if (ii-f_i*T_duration)*(ii-f_f*T_duration)>0
        rcv_f(ii,:) = zeros(1, 64); 
    end        
end

final_sum = beamform2(rcv_f, 'rect', aperture1); 
time_vector = linspace(t_start, t_end, size(rcv_t, 1));
figure; 
imagesc(sn,time_vector, 10*log10(abs(final_sum').^2));
axis xy; 
setFont(18); setFigureAuto; 
xlabel('sin \theta'); ylabel('Time (seconds)'); 
%% 
hold on; 

for k = 1:size(final_sum, 2)
    maxval = max(abs(final_sum(:, k))); 
    ind1 = find(abs(final_sum(:, k))==maxval); 
    plot(sn(ind1), time_vector(k), 'k*'); 
end

%% starting from peak 
max0 = max(max(abs(final_sum))); 
[m, n] = find(abs(final_sum) == max0); 
S = []; 
T = []; 
th = 35; %20 dB down threshold; 

for n1 = n:(n+8e3)
    max1 = max(abs(final_sum(:, n1))); 
    if 20*log10(max1) > (20*log10(max0) - th)
        m1 = find(abs(final_sum(:, n1))==max1); 
        plot(sn(m1), time_vector(n1), 'ko'); 
        S = [S; sn(m1)]; 
        T = [T; time_vector(n1)]; 
    end 
end

%% finish up the tracking of beam intensity  migration before matched filtering - Do evaluation 



figure; plot(S, T, '*'); hold on; 
P1 = polyfit(S, T, 1); 
Tfit = polyval(P1, S); 
plot(S, Tfit, '--r'); 
