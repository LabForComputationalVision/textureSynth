function X=SNR(s,n);

% Compute the signal-to-noise ratio in dB
%  	X=SNR(signal,noise);
% (it does not subtract the means).

es=sum(sum(abs(s).^2));
en=sum(sum(abs(n).^2));
X=10*log10(es/en);

