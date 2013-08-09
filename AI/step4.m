% step 4, range estimate

zoom =  axis; 
sn_inds = find((sn(:)>zoom(1))&(sn(:)<zoom(2))); 
sn_zoom = (sn(sn_inds))'; 

time_inds = find((time_vector(:)>zoom(3))&(time_vector(:)<zoom(4))); 
time_zoom = time_vector(time_inds); 

zoomed_btdata = 10*log10(abs(final_sum').^2);

zoomed_btdata = zoomed_btdata(sn_inds, time_inds); 
zoomed_btdata = zoomed_btdata'; 


max_value = max(max(zoomed_btdata));

crit = max_value - 25; %the optimum threshold level for array invariant is about 30 dB below the peak value of the matched filter, for SNR from 5dB onwards.
%for the whale: set it lower 


[m,n] = find(zoomed_btdata == max_value);
whale_time = time_zoom(m)
whale_bearing = asind(sn_zoom(n))
t_0 = time_zoom(m);
idx_t = [];
idx_smax = [];
jj = 1;

for ii = m:length(time_zoom)
% for ii = 1:m
    temp = max(zoomed_btdata(ii,:));
    if temp>=crit
        [int_max(jj) idx_smax(jj)] = max(zoomed_btdata(ii,:));
        jj = jj + 1;
        idx_t = [idx_t;ii];
    end
end

T2 = [time_zoom(idx_t)' ones(length(idx_smax),1)];
temp2 = inv(transpose(T2)*T2)*transpose(T2)*(sn_zoom(idx_smax));
temp3 = inv(transpose(T2)*T2)*transpose(T2)*(1./(sn_zoom(idx_smax)));
chi_h = temp3(1);
d_h = temp3(2);

chi_l = temp2(1);
d_l = temp2(2);
hat_s = chi_l*time_zoom(idx_t)+d_l;
hat_sh = 1./(chi_h*time_zoom(idx_t)+d_h);
adjust = hat_s(1)-sn_zoom(n);
adjust1 = hat_sh(1)-sn_zoom(n);
% adjust = hat_s(end)-sn_zoom(n);
% adjust1 = hat_sh(end)-sn_zoom(n);
s_est = hat_s-adjust;
s_est1 = hat_sh-adjust1;
r_est = -1500/chi_l*sn_zoom(n)
r_esth = 1500/(chi_h*sn_zoom(n))
hold on; 
 plot(s_est,time_zoom(idx_t),'k-','linewidth',2)
plot([sn_zoom(n) sn_zoom(n)],[time_zoom(1) time_zoom(end)],'k-.','linewidth',1.5)
h = legend('$$\hat{s}_l(t)$$','$$sin\hat{\theta}_0$$');
set(h,'interpreter','latex','fontsize',16,'fontweight','bold','location','northeast')