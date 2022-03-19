%% Samaa Hany Seif Elyazal
%% Wireless Communication, Intake 42
%% Transmit Diversty (MRT) With CSI
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
ber_MRT = 0;
for(j = 1:1:NS)
I = randi([0 1], 1, 1); 
Q = randi([0 1], 1, 1);
S = ((2*I - 1) + 1i*(2*Q - 1))*sqrt(1/2);
I1 = randi([0 1], 1, 1); 
Q1 = randi([0 1], 1, 1);
S1 = ((2*I1 - 1) + 1i*(2*Q1 - 1))*sqrt(1/2);
I2 = randi([0 1], 1, 1); 
Q2 = randi([0 1], 1, 1);
S2 = ((2*I2 - 1) + 1i*(2*Q2 - 1))*sqrt(1/2);
%abs(sqrt(P)*S).^2 %check tx power
%%AWGN
W = (randn(1, 1) + 1i*randn(1, 1))*sqrt(No/2);
YAWGN = sqrt(P)*S + W;
%%DECODING AWGN
I_HAT = real(YAWGN) > 0;
Q_HAT = imag(YAWGN) > 0;
ber = ber + (((I ~= I_HAT)+(Q ~=Q_HAT)));
%% Fading
h1 = 1/sqrt(2)*(randn(1,1) + 1i*randn(1,1));
h2 = 1/sqrt(2)*(randn(1,1) + 1i*randn(1,1));

a=conj(h1)/sqrt(abs(h1).^2+abs(h2).^2);
b=conj(h2)/sqrt(abs(h1).^2+abs(h2).^2);

%% Decoding fading
y_received = a*h1*S + b*h2*S +W;
%% Decoding AWGN
I_HAT = real(y_received) > 0;
Q_HAT = imag(y_received) > 0;
ber_MRT = ber_MRT + (((I ~= I_HAT)+(Q ~=Q_HAT)));
end
BER(K) = ber / NB;
BER_MRT(K) = ber_MRT / NB;

end
%%plot
semilogy(SNRV, BER, '--r*', SNRV, TH_ERROR, 'b-o', SNRV, BER_MRT, '--m*', SNRV, 1./snrv.^3, 'black-s' )
legend('Monte Carlo', 'Theoritical', 'MRT Signal', 'TH Scale')
title('Samaa Hany')
xlabel('SNR')
ylabel('BER')
dhat=(log(BER_MRT(end))-log(BER_MRT(end-1)))/(-log(snrv(end))+log(snrv(end-1)))
axis([min(SNRV), max(SNRV), 1e-4,1])