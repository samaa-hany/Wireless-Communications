%% Samaa Hany Seif Elyazal
%% Wireless Communication, Intake 42
%% EGC Diversty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;
%% Diversity QPSK
SNRV = -5:3:35; %SIGNAL TO NOISE RATIO IN DB
for(K=1:length(SNRV))
SNR = SNRV(K);
snr = 10^(SNR / 10);
snrv(K) = snr;
P = 1; %RX POWER
No = P / snr;
NS = 100000; %NO OF QPSK SYBOLS
M = 4; %M ARRAY OF QSPK
NB = NS*log2(M); %NO OF BITS PER SYMOBOLS
%% THEORITICAL ERROR RATE
TH_ERROR(K) = qfunc(sqrt(snr));
%% GENTERATE SYMBOLS
ber = 0;
ber_MRC = 0;
ber_EGC = 0;
for(j = 1:1:NS)
I = randi([0 1], 1, 1); 
Q = randi([0 1], 1, 1);
S = ((2*I - 1) + 1i*(2*Q - 1))*sqrt(1/2);
%abs(sqrt(P)*S).^2 %check tx power
%% AWGN
W = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);
YAWGN = sqrt(P)*S + W;
%% DECODING AWGN
I_HAT = real(YAWGN) > 0;
Q_HAT = imag(YAWGN) > 0;
ber = ber + (((I ~= I_HAT)+(Q ~=Q_HAT)));
%% Equal Gain Combiner
W1 = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);
W2 = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);

h1 = 1/sqrt(2)*(randn(1,1) + 1i*randn(1,1));
h2 = 1/sqrt(2)*(randn(1,1) + 1i*randn(1,1));

hc1=conj(h1)/abs(h1);
hc2=conj(h2)/abs(h2);

yfading1 = h1*S + W1;
yfading2 = h2*S + W2;
%% Decoding EGC
S_HAT_MCR = yfading1*conj(h1)+yfading2*conj(h2);
S_HAT_EGC = yfading1*hc1+yfading2*hc2;
%% Decoding AWGN
I_HAT_MCR = real(S_HAT_MCR) > 0;
Q_HAT_MCR = imag(S_HAT_MCR) > 0;
I_HAT_EGC = real(S_HAT_EGC) > 0;
Q_HAT_EGC = imag(S_HAT_EGC) > 0;
ber_MRC = ber_MRC + (((I ~= I_HAT_MCR)+(Q ~=Q_HAT_MCR)));
ber_EGC = ber_EGC + (((I ~= I_HAT_EGC)+(Q ~=Q_HAT_EGC)));
end
BER(K) = ber / NB;
BER_MRC(K) = ber_MRC / NB;
BER_EGC(K) = ber_EGC / NB;

end
%% Ploting
semilogy(SNRV, BER, '--r*', SNRV, TH_ERROR, 'b-o', SNRV, BER_MRC, '--m*',SNRV, BER_EGC, '--g*', SNRV, 1./snrv.^2, 'black-s' )
legend('Monte Carlo', 'Theoritical', 'MRC Fading', 'EGC Fading','TH Scale')
title('Samaa Hany')
xlabel('SNR')
ylabel('BER')

%% Diversty Order
d_hat_MRC=(log(BER_MRC(end))-log(BER_MRC(end-5)))/(-log(snrv(end))+log(snrv(end-5)))
d_hat_EGC=(log(BER_EGC(end))-log(BER_EGC(end-5)))/(-log(snrv(end))+log(snrv(end-5)))
axis([min(SNRV), max(SNRV), 1e-4,1])