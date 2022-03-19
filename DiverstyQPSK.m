%% Samaa Hany Seif Elyazal
%% Wireless Communication, Intake 42
%% Diversty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;
%diversity QPSK
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
%%THEORITICAL ERROR RATE
TH_ERROR(K) = qfunc(sqrt(snr));
%%GENTERATE SYMBOLS
ber = 0;
ber_f = 0;
for j = 1:1:NS
I = randi([0 1], 1, 1); 
Q = randi([0 1], 1, 1);
S = ((2*I - 1) + 1i*(2*Q - 1))*sqrt(1/2);
%abs(sqrt(P)*S).^2 %check tx power
%%AWGN
W = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);
YAWGN = sqrt(P)*S + W;
%%DECODING AWGN
I_HAT = real(YAWGN) > 0;
Q_HAT = imag(YAWGN) > 0;
ber = ber + (((I ~= I_HAT)+(Q ~=Q_HAT)));
%% Fading
h = 1/sqrt(2)*(randn(1,1) + 1i*randn(1,1));
yfading = h*S + W;
%% Decoding fading
S_HAT = yfading /h;
%% Decoding AWGN
I_HAT = real(S_HAT) > 0;
Q_HAT = imag(S_HAT) > 0;
ber_f = ber_f + (((I ~= I_HAT)+(Q ~=Q_HAT)));
end
BER(K) = ber / NB;
BERF(K) = ber_f / NB;

end
%%plot
semilogy(SNRV, BER, '--r*', SNRV, TH_ERROR, 'b-o', SNRV, BERF, '--m*', SNRV, 1./snrv.^3, 'black-s' )
legend('Monte Carlo', 'Theoritical', 'MC Fading', 'TH Scale')
title('Samaa Hany')
xlabel('SNR')
ylabel('BER')
dhat=(log(BERF(end))-log(BERF(end-1)))/(-log(snrv(end))+log(snrv(end-1)))
axis([min(SNRV), max(SNRV), 1e-4,1])