%step3


figure; h =plot_spectrogram(real(whale_beamform), 512, 384, 8192, 8000);
threshold_dB_down = 15; 
f_samp = 8e3; 
rcv_t = input_series(round(t_start*f_samp):round(t_end*f_samp), :); 
[replica_t, replica_f] =  replica_creator2(h, threshold_dB_down, f_samp);

replica_t = replica_t(1:size(rcv_t, 1));
replica_f = ifft(replica_t)*length(replica_t)/8000; 

h = plot_spectrogram(real(replica_t), 512, 384, 8192, 8000);


f = linspace(0, f_samp, length(rcv_t)); 
T_w = length(replica_t)/f_samp; 
MF_t = replica_t;
T_m = 0; 
MF_f =(T_w.*ifft(MF_t))'.*exp(+j*2*pi*f*T_m); 

mf_spectrum = zeros(size(rcv_f)); 
mf_out = zeros(size(rcv_t)); 

for ii = 1:size(rcv_f, 2)
    mf_spectrum(:, ii) = rcv_f(:, ii).*transpose(MF_f); 
    mf_out(:, ii) = fft(mf_spectrum(:, ii))./T_w;
end


figure; plot(real(mf_out(:, 1)))

% matched filter with all the beams: 
final_sum_f = ifft(transpose(final_sum)); 
MF_array = repmat(transpose(MF_f), 1, size(final_sum, 1)); 
MF_out_bf = MF_array.*final_sum_f; 
MF_out_bf_t = fft(MF_out_bf); 

figure; imagesc(sn, time_vector, 20*log10(abs(MF_out_bf_t))); 
axis xy; 
final_sum = MF_out_bf_t; 


