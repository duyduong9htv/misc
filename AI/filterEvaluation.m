%this code does a formal evaluation of the bearing measurements and SNR
%obtained by the non-linear and linear MF 



cd /Volumes/scratch/duong/whale_localization_data/Tracks_data/track570_4/DAT

% input_series = read_input_series('fora2006jd275t040730.dat', 'lf'); 


!ls fora*0.dat > fora_list_415
fid = fopen('fora_list_415'); 
while ~feof(fid)
    %read file, find peak arrival time, silly things, etc. 
    filename = fgetl(fid); 
    input_series = read_input_series(filename, 'lf'); 
    T = 20; 
    t_start = 0.1; t_end = 20;  f_samp = 8e3; 
    rcv_t1 = input_series(round(t_start*f_samp):round(t_end*f_samp), :); 
    rcv_f1 = ifft(rcv_t1)*(t_end-t_start); 
    fStart = 350; fEnd = 500;
    [rcv_t, rcv_f] = BandPassFilter(rcv_t1, fStart, fEnd, 8000);   
    ind = find(rcv_t(:, 1)==max(rcv_t(:, 1))); 
    t_peak = ind/f_samp; 
    t_start = t_peak - 1; 
    t_end = t_peak +2; 
    rcv_t1 = input_series(round(t_start*f_samp):round(t_end*f_samp), :); 
    rcv_f1 = ifft(rcv_t1)*(t_end-t_start);      
    final_sum = beamform2(rcv_f1,'hanning', 'lf'); 
    temp = sum(abs(final_sum').^2); 
    ind = find(temp(:) == max(temp)); 
    sn = linspace(-1, 1, 401); 
    whale_beamform = steer_beam(rcv_f1, 'hann', 'lf', sn(ind)); 
    whaleF = ifft(whale_beamform)*size(whale_beamform, 1)/f_samp; 

    figure; h =plot_spectrogram(real(whale_beamform), 512, 384, 8192, 8000);
    xlim([fStart fEnd])

    % ylim(timeLimSet)
    threshold_dB_down = 20; 
    f_samp = 8e3; 
    rcv_t = input_series(round(t_start*f_samp):round(t_end*f_samp), :); 
    [replica_t, replica_f, fTrace, tTrace] =  replica_creator2(h, threshold_dB_down, f_samp);

    %% plot received spectrogram with peak tracing 
    figure(54); clf;
    hold on; 
    [S, F, T, P] = spectrogram(real(whale_beamform), 512, 384, 8192, 8000); 
    s(1) = subplot(2, 3, 1); 
    P = P./max(max(abs(P))); 
    imagesc(T, F, 10*log10(P)); 
    hold on; 
    l1 = plot(tTrace, fTrace, '--', 'color', [0.7 0.7 0.7], 'linewidth', 2);
    l2 = plot(linspace(tTrace(1)-0.025, tTrace(1)-0.025+1, 100), linspace(390, 440, 100), '-', ...
              'linewidth', 2, 'color', [0.7 0.7 0.7]);
    legend([l2 l1], {'Linear kernel', 'Non-linear kernel'})
    for adjustfigure = 1    
    axis xy; 
    caxis([-40 0])
    colormap(flipud(gray))
    ylim([fStart fEnd])    
    set(gca, 'xtick', [0:0.5:5])
    xlim([tTrace(1)-0.05 tTrace(end) + 0.05]); 
    xWidth = tTrace(end) - tTrace(1); 
    setFont(22); 
    setFigureAuto; 
    xlabel('Reduced travel time (s)'); 
    ylabel('Frequency (Hz)'); 
    c = colorbar('location', 'south'); 
    caxis([-40 0])
    set(c, 'position', [0.14 0.5918 0.1 0.018])
    set(c, 'fontsize', 16, 'fontweight', 'bold'); 
    set(c, 'xcolor', 'k', 'ycolor', 'k')
    ti = get(c, 'title'); 
    set(ti, 'string', 'dB', 'color', 'k', 'fontsize', 16, 'fontweight', 'bold')
    end 

    
    %% matched-filter spectrogram 
    s(2) = subplot(2, 3, 2); 
    [S, F, T, P] = spectrogram(real(replica_t), 512, 384, 8192, 8000); 
    P = P./max(max(abs(P))); 
    imagesc(T, F, 10*log10(P)); 

    for adjustFigure = 1
    axis xy; 
    caxis auto
    colormap(flipud(gray))
    ylim([fStart fEnd])
    set(gca, 'xtick', [0.02:0.5:5])
    xtl = get(gca, 'xtick'); 
    set(gca, 'xticklabel', xtl-0.02); 
    xlim([0.02 0.02+xWidth])

    % xlim([0 2])
    setFont(22); 
    setFigureAuto; 
    xlabel('Time (s)'); 

    c1 = colorbar('location', 'south'); 
    caxis([-40 0])
    set(c1, 'position', [0.42 0.5918 0.1 0.018])
    set(c1, 'fontsize', 16, 'fontweight', 'bold');
    set(c1, 'xcolor', 'k', 'ycolor', 'k')
    ti = get(c1, 'title'); 
    set(ti, 'string', 'dB', 'color', 'k', 'fontsize', 16, 'fontweight', 'bold')
    end 
    title(filename); 
    
    
    %% comparison received spectrum and matched-filter spectrum 
    s(3) = subplot(2, 3, 3); 

    plot(linspace(0, 8000, size(whaleF, 1)),...
                    20*log10(abs(whaleF)./max(abs(whaleF))),...
                    '--', 'linewidth', 2, 'color', [0.6 0.6 0.6]); 
    hold on; 
    plot(linspace(0, 8000, size(replica_f, 1)), ...
            20*log10(abs(replica_f)./max(abs(replica_f))), ...
            '-k', 'linewidth', 2); 

    for adjust = 1 
    xlim([fStart fEnd]); 
        ylim([-45 10]);

    setFont(22); 
    xlabel('Frequency (Hz)'); 
    ylabel('dB'); 
    legend('Received signal', 'Matched-filter'); 
    setFigureAuto; 
    end 
    
    %% received signal in time
    s(4) = subplot(2, 3, 4); 
    whale_beamform = steer_beam(rcv_f1, 'hann', 'lf', -1); 
    [whale_beamform1, ~] = BandPassFilter(whale_beamform, fStart, fEnd, 8000); 
    plot(1/8000*[1:1:length(whale_beamform1)], real(whale_beamform1)/(10.^(208/20)), '-k'); 

    setFigureAuto
    setFont(22); 
    xlabel('Reduced travel time (s)'); 
    ylabel('Normalized received pressure level ( \mu Pa)'); 
    legend('Received signal in time'); 

    %% matched-filter kernel in time 
    s(5) = subplot(2, 3, 5); 
    plot([linspace(0, 1, 8000) 1+1/8000*[1:1:(length(replica_t)-8000)]],...
        [zeros(8000, 1); real(replica_t(1:end-8000))], '-k'); 
    setFont(22); 
    legend('Matched-filter kernel'); 
    xlabel('Time (s)'); 
    ylabel('$h(t)$', 'interpreter', 'latex'); 

    %% MF output 

    f2 = linspace(0, 8000, length(rcv_f)); 
    f1 = linspace(0, 8000, length(replica_f)); 
    T_w = length(replica_f)/f_samp; 
    MF_f = interp1(f1, replica_f, f2); 
    T_m = -1; 
    f = linspace(0, f_samp, length(rcv_f)); 
    MF_f = conj(MF_f.*exp(+j*2*pi*f*T_m)); 

    mf_spectrum = zeros(size(rcv_f)); 
    mf_out = zeros(size(rcv_f)); 

    for ii = 1:size(rcv_f, 2)
        mf_spectrum(:, ii) = rcv_f(:, ii).*transpose(MF_f); 
        mf_out(:, ii) = fft(mf_spectrum(:, ii))./(length(rcv_f(:, ii))/f_samp);
    end

    s(6) = subplot(2, 3, 6); 
    plot(1/8000*[1:1:size(mf_out, 1)], (abs(mf_out(:, 4)).^2)/(10^20.8), '-k', 'linewidth', 2)
    hold on; 
    setFont(22); 
    xlabel('Reduced travel time (s)'); 
    ylabel('$|{\rm \it MF}(t|ti)|^2$', 'interpreter', 'latex'); 
    legend('Matched-filter output'); 
    
    %%linear MF 

    [st1, Sf1, rst1, rSf1] = getLFM(8000, 0, 1, length(rcv_f)/f_samp, 390, 440); 
    % Sf1 = [Sf1 0]; 

    for ii = 1:size(rcv_f, 2)
        mf_spectrum_linear(:, ii) = rcv_f(:, ii).*transpose(conj(Sf1)); 
        mf_out_linear(:, ii) = fft(mf_spectrum_linear(:, ii))./(length(rcv_f(:, ii))/f_samp);
    end
    
    plot(0.5+1/8000*[1:1:size(mf_out, 1)], (abs(mf_out_linear(:, 4)).^2)/(10^20.8), ...
        '--', 'linewidth', 2, 'color', [0.5 0.5 0.5])
    xlim([5 8])
    
    xtl = get(gca, 'xtick'); 
    set(gca, 'xticklabel', xtl-5); 
    legend('Non-linear', 'Linear'); 
    
    eval(['print -dtiff MF_evaluation_' filename(1:(end-4)) '.tif']); 


end 



