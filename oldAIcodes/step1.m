%step 1
% ping_name = num2str(track_whale_table(ping_number,2)); 

while length(ping_name)<6
    ping_name = ['0' ping_name]; 
end
if aperture == 'hf' 
    aperture1 = 'mf';
else 
    aperture1 = aperture; 
end
tic
eval(['input_series = read_input_series(''fora2006jd' jd 't' ping_name '.dat'', aperture1);']);
toc

%% spectrogram
input_series = input_series(1:4:end, :); 
P1 = 0; 
tic
for ii = 1:1:64
    disp(ii)
    [S, F, T, P] = spectrogram(input_series(round(2000*linspace(t1, t2, round((t2-t1)*2000))), ii), 512, 256, 1024, 2000); 
    P1 = P1+P; 
end
toc
% select whale call of interest 
figure; imagesc(F, T+t1, 10*log10(P1')); 
% caxis([40 90])
title(ping_name)