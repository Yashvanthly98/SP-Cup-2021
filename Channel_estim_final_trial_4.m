% O.N.N % S.K.P
tic;
clear all;
%clc;
load('dataset2.mat');

%% OFDM - Flat Channel estimate - after nulling out insignificant taps
for l=1:50
H_estimate(:,:,l)=diag(transmitSignal)\receivedSignal(:,1:end,l);
H_irs_est_1(:,:,l) = (1/4096)*H_estimate(:,:,l)*pilotMatrix';
H_irs_est_time = ifft(H_irs_est_1(:,:,l));
H_irs_time_new=H_irs_est_time(1:20,:);
H_t=fft(H_irs_time_new,500);
h_eff_final(:,:,l)=H_t.';
end
toc;

save('h_eff_final.mat','h_eff_final');