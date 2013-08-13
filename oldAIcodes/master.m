%master program. Calls out 4 different steps to 
%   -step 1     : plots out spectrogram of the received whale call
%   -step 2     : beamforms and gets the whale call with higher SNR to
%                 design the matched filter
%   -step 3     : apply matched filter on the beamformed data and plots out
%                 the migration
%   -step 4     : get range estimate using the array invariant 

%note: can use table : whaleSounds_MF_571_5.csv or any of those .csv table
%from Dave & Hari to dig out the time for t1 and t2 containing the whale
%calls.
addpath /Users/dtran/Research/whale_localization/ % ping_number = 39
clear; close all; 
ping_name = '223000'; 
jd = '276'
if ping_name(end)=='0'
aperture = 'lf'; 
else
    aperture = 'hf'; 
end
%% step 1
t1 =0.1%start time 
t2 =5%end time for spectrogram. 

% ping_name = num2str(track_whale_table(ping_number,2)); 
while length(ping_name)<6
    ping_name = ['0' ping_name]; 
end
if aperture == 'hf' 
    aperture1 = 'mf';
else 
    aperture1 = aperture; 
end
eval(['input_series = read_input_series(''fora2006jd' jd 't' ping_name '.dat'', aperture1);']);

%% spectrogram
P1 = 0; 
for ii = 1:64
    disp(ii)
    [S, F, T, P] = spectrogram(input_series(round(8000*linspace(t1, t2, round((t2-t1)*8000))), ii), 512, 384, 8192, 8000); 
    P1 = P1+P; 
end
% select whale call of interest 
figure; imagesc(T+t1, F, 10*log10(P1)); 

caxis([40 90])
title(ping_name)

%% step 2
step2; %beamform and get out whale call to design matched filter 
% zoom in to the whale beam after this step

%% step 2a
zoom =  axis; 
sn_inds = find((sn(:)>zoom(1))&(sn(:)<zoom(2))); 
sn_zoom = (sn(sn_inds))'; 
time_inds = find((time_vector(:)>zoom(3))&(time_vector(:)<zoom(4))); 
time_zoom = time_vector(time_inds); 
zoomed_btdata = 10*log10(abs(final_sum).^2);
zoomed_btdata = zoomed_btdata(sn_inds, time_inds); 
zoomed_btdata = zoomed_btdata'; 
max_value = max(max(zoomed_btdata));
[m,n] = find(zoomed_btdata == max_value);
whale_beamform = final_sum(round((sn_zoom(n)+1)/0.005), :);
step3; % apply matched filter, beamform, and plot out beam migration
%% step 4
step4; % get range estimate using the array invariant

eval(['save whale_' jd '_' ping_name '_' num2str(round(whale_time)) '_' num2str(round(abs(whale_bearing))) '.mat whale_time whale_bearing r_est r_esth']); 

%% 
% add whale heading column to it 
load whale_table_track_570_1.csv; 
ping_list = whale_table_track_570_1(:, 2)*10000 + whale_table_track_570_1(:, 3)*100 + whale_table_track_570_1(:, 4); 
load array_heading
array_heading_whale = [];
for ii = 1:length(ping_list)
    ind = find(array_heading(:, 1) == ping_list(ii)); 
    array_heading_whale = [array_heading_whale; array_heading(ind, 2)];
end
whale_table_track_570_1 = [whale_table_track_570_1 array_heading_whale];
whale_heading = whale_table_track_570_1(:, 9) + 90 - whale_table_track_570_1(:, 6); 
whale_table_track_570_1 = [whale_table_track_570_1 whale_heading]; 
whale_table_track_570_1(:, 10) = whale_table_track_570_1(:, 10) - 180; 

 %%  add another 2 columes with  
 % a. range from source
 % b. angle from source, w.r.t true North. 
 % i.e. convert all the whale location into polar coordinates 
 
 


