% function  [replica_t, replica_f] =  replica_creator2(threshold_dB_down, f_samp);
% creates the non linear replica of a whale call based on current spectrogram figure of the whale call. 
% last modified Dec 02 '10: for flexibility of f_samp
% input : threshold_dB_down : how much below the maximum amplitude you want the peak detector to trace 
% output:  replica_t, replica_f : time and frequency domain replica signal (normalized)

function [replica_t, replica_f, fTrace, tTrace]=replica_creator2(h, threshold_dB_down, f_samp)


%current figure has to be one that is showing the spectrogram of whale call
%for peak detector to work
zoom = axis; 
% line_handles = get(gca, 'children'); 
% F = get(line_handles(1),'XData'); 
% T = get(line_handles(1),'ydata'); 
% S1 = get(line_handles(1), 'Cdata'); %get image indices, already plotted
%  in dB 
F = get(h, 'xdata');
T = get(h, 'ydata'); 
S1 = get(h, 'Cdata'); 

% S1 = get(h, 'Cdata'); 
S1 = S1'; 

%zoom to a window containing the whale call 
t_zoom_inds = find((T(:)-zoom(3)).*(T(:)-zoom(4))<0); 
f_zoom_inds = find((F(:)-zoom(1)).*(F(:)-zoom(2))<0); 

S_zoom = S1(f_zoom_inds,t_zoom_inds); 
T_zoom = T(t_zoom_inds); 
F_zoom = F(f_zoom_inds); 


call_duration = T_zoom(end)  - T_zoom(1);
S_max = max(max(S_zoom)); 
% 
% figure; imagesc(F_zoom, T_zoom, S_zoom'); 
% axis xy
% set(gca, 'fontsize', 16); 
% xlabel('Frequency (Hz)'); 
% ylabel('Time (seconds)'); 

f_zoom = F_zoom; 
time_zoom = T_zoom; 

S_log = S_zoom'; 


threshold = S_max - threshold_dB_down; %set the threshold for peak tracer to work

replica_t_f = []; 

for tt = 1:length(time_zoom)
    max_S = max(S_log(tt, :)); 
    ff = find(S_log(tt, :)==max_S); 
    if max_S > threshold
        replica_t_f = [replica_t_f; time_zoom(tt)  f_zoom(ff) sqrt(10^(max_S/10))];
    end  
end
replica_t_f_stored = replica_t_f; 

% figure(1); hold on; 
% plot(replica_t_f(:, 2), replica_t_f(:, 1), '-k', 'linewidth', 2); 
% % set(gcf, 'position', position)
% set(gcf, 'paperPositionMode', 'auto'); 
% title('Spectrogram of received whale call'); 

% 
NI = 10; %interpolation factor if the replica time and frequency samples are too sparse
a = interp(replica_t_f_stored(:,1), NI);  %time
b = interp(replica_t_f_stored(:,2), NI);  %frequency
c = interp(replica_t_f_stored(:,3), NI);  %amplitude 

%output 
fTrace = b; 
tTrace = a; 

replica_t_f = [a b c]; 

% 
% hold on; 
% plot(replica_t_f(:, 2), replica_t_f(:, 1), '-k', 'linewidth', 2); 
set(gca,'position',[0.09 0.11 0.76 0.76])
% set(gcf, 'position', position)
set(gcf, 'paperPositionMode', 'auto'); 
c = colorbar;
set(c,'linewidth',1.5,'fontsize',12,'fontweight','bold')
set(c,'position',[0.88 0.115 0.03 0.5])
ti = get(c,'title');
set(ti,'string','          dB re 1{\mu}Pa','fontsize',12,'fontweight','bold')
title('Matched filter designing '); 

%% generate the whale call replica 

t_start = min(replica_t_f(:, 1)); 
t_end = max(replica_t_f(:, 1)); 

T_w = call_duration + 1.5; %extend the whale window to 1second before and 1 after 
t = linspace(t_start-t_start, T_w, round(f_samp*T_w)); 
dt = t(2) - t(1); 
df = 1/T_w; 
ind = 0; 

temp = zeros(length(replica_t_f(:,1)), length(t)); 
A = zeros(length(t), 1); 
synthetic_signal = zeros(length(t), 1); 
f = zeros(size(A));


%interpolate for a frequency sampling of 8000 Hz
for kk = 1:length(replica_t_f)-1
    window_start = replica_t_f(kk, 1)-t_start;
    window_end = replica_t_f(kk+1, 1)-t_start; 
    f_start = replica_t_f(kk, 2); 
    f_end = replica_t_f(kk+1, 2); 
    A_start = replica_t_f(kk, 3); 
    A_end = replica_t_f(kk+1, 3); 
    t_inds = find((window_start<=t)&(t<=window_end));
    f(t_inds) = f_start + (t(t_inds)-window_start)*(f_end - f_start)/(window_end - window_start); 
    A(t_inds) = A_start + (t(t_inds)-window_start)*(A_end - A_start)/(window_end - window_start); 
end


phi = zeros(size(A)); 
temp= 0; 

%generate synthetic signal. because d\phi/dt = f, the phase \phi = \int f*dt 

for ii = 1:length(t)
    temp = temp + f(ii)*dt; 
    phi(ii) = temp; 
    synthetic_signal(ii) = A(ii)*exp(-j*2*pi*phi(ii)); 
end




replica_t = synthetic_signal ;
replica_f = ifft(replica_t)*T_w; 

%normalization for unit filter energy

signal_energy  = sum(abs(replica_f).^2*df)

replica_f = replica_f./sqrt(signal_energy); 
replica_t = replica_t./sqrt(signal_energy); 
% 
% 
% figure; plot(t, real(replica_t));xlabel('Time (s)'); ylabel('s(t)');  title('replica in time'); 
% figure; plot(linspace(0, f_samp, round(f_samp*T_w)), abs(replica_f)); xlabel('Frequency'); ylabel('|S(f)|'); title('replica in frequency'); 
% %% check : plot spectrogram of synthetic whale call 
% % 
% % [s, f, t, p] = spectrogram(real(synthetic_signal), 512, 384, 8192, 8000); 
% % figure; imagesc(f, t, 10*log10(p')); 
% axis xy
% xlim([200 1000])
% % ylim([0 1.2])

