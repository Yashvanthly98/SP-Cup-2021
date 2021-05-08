tic;
clear all;
load('dataset2.mat');

 %% OFDM - Flat Channel estimate - after nulling out insignificant taps - sorted version
for l=1:50
H_estimate(:,:,l)=diag(transmitSignal)\receivedSignal(:,1:end,l);
H_irs_est_1(:,:,l) = (1/4096)*H_estimate(:,:,l)*pilotMatrix';
H_irs_est_time = ifft(H_irs_est_1(:,:,l));
%h_sorted = sort(ifft(H_irs_est_1(:,:,l)),'descend'); % For comparison
[~,max_index] = maxk(abs(H_irs_est_time),20);  % search the largest 20 entries in each column
H_sort=zeros(500,4096);  % Nulling matrix
    for o=1:4096
         H_sort(max_index(:,o),o) = H_irs_est_time(max_index(:,o),o); % Assign the largest 20 coefficients to nulling matrix retaining the same indices as in ifft - time matrix
    end
H_t=fft(H_sort,500);
h_eff_final(:,:,l)=H_t.';
end

save('h_eff_final.mat','h_eff_final');