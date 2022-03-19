%% Samaa Hany Seif Elyazal
%% Wireless Communication, Intake 42
%% Selective Diversty
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
ber_f = 0;
ber_f1 = 0;
ber_f2 = 0;
for(j = 1:1:NS)
I = randi([0 1], 1, 1); 
Q = randi([0 1], 1, 1);
S = ((2*I - 1) + 1i*(2*Q - 1))*sqrt(1/2);
%abs(sqrt(P)*S).^2 %check tx power
%% AWGN
W = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);
YAWGN = sqrt(P)*S + W;
%% Decoding AWGN
I_HAT = real(YAWGN) > 0;
Q_HAT = imag(YAWGN) > 0;
ber = ber + (((I ~= I_HAT)+(Q ~=Q_HAT)));
%% Selective
W1 = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);
W2 = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);

h1 = 1/sqrt(2)*(randn(1,1) + 1i*randn(1,1));
h2 = 1/sqrt(2)*(randn(1,1) + 1i*randn(1,1));

yfading1 = h1*S + W1;
yfading2 = h2*S + W2;
%% Decoding Selective
S_HAT1 = yfading1*conj(h1);
S_HAT2 = yfading2*conj(h2);
%% Decoding AWGN
I_HAT1 = real(S_HAT1) > 0;
Q_HAT1 = imag(S_HAT1) > 0;
I_HAT2 = real(S_HAT2) > 0;
Q_HAT2 = imag(S_HAT2) > 0;
ber_f = ber_f + (((I ~= I_HAT)+(Q ~=Q_HAT)));
ber_f1 = ber_f1 + (((I ~= I_HAT1)+(Q ~=Q_HAT1)));
ber_f2 = ber_f2 + (((I ~= I_HAT2)+(Q ~=Q_HAT2)));

end
ber_SC=min(ber_f1,ber_f2);
BER(K) = ber / NB;
BER_F(K) = ber_f1 / NB;
BER_SEL1(K) = ber_f1 / NB;
BER_SEL2(K) = ber_f2 / NB;
BER_SC(K) = ber_SC / NB;

end
%% Ploting

semilogy(SNRV, BER,'--r*',SNRV, TH_ERROR,'b-o',SNRV, BER_SC,'--g*', SNRV, 1./snrv.^2,'black-s' )
legend('Monte Carlo', 'Theoritical','Selective Diversty','TH Scale')
title('Samaa Hany')
xlabel('SNR')
ylabel('BER')

%% Diversty Order
d_hat=(log(BER_F(end))-log(BER_F(end-1)))/(-log(snrv(end))+log(snrv(end-1)))
d_hat_SC=(log(BER_SC(end))-log(BER_SC(end-1)))/(-log(snrv(end))+log(snrv(end-1)))

axis([min(SNRV), max(SNRV), 1e-4,1])