close all
N=1000;
SNR=-5:3:10; %Single to Noise Ratio in dB
snr=10.^(SNR/10);
P=1;    %transmiter Power
Ns=1000;   %Number of QPSK symbols
M=4;    %M-ary for QPSK M=4
Nb=Ns*log2(M); %Number of bits per symbols
%% Generate symbols
ber=0;
ber_S=0;
for x=1:1:Ns
    I=randi([0 1], 1,N);
    Q=randi([0 1], 1,N);
    S=((2*I-1)+(1i*Q-1));
%% AWGN
No=P./snr;
for x=1:1:length(SNR)
W=(sqrt(No(x)/2)).*(randn(1,Ns)+1i*randn(1,Ns));
R_awgn=(sqrt(P)*S)+W;
%% Decoding AWGN
R_RealHAT=real(R_awgn)>0;
R_ImgHAT=imag(R_awgn)>0;
ber=ber+((R_RealHAT~=I)+(R_ImgHAT~=Q));
%% Fading
h=(sqrt(1/2)).*(randn(1,N)+1i*randn(1,N));
R_Fad=h.*S+W;
%% Decoding Fading
S_HAT=R_Fad/h;
S_RealHAT=real(S_HAT)>0;
S_ImgHAT=imag(S_HAT)>0;
ber_S=ber_S+((S_RealHAT~=I)+(S_ImgHAT~=Q));
ber_th=qfunc(sqrt(2/No(x)));
end
BER(x)=ber/N;
BER_S(x)=ber_S/N;
BER_th(x)=ber_th/N;
end

semilogy(SNR,BER,'--r*',SNR,BER_S,'--g*',SNR,BER_th,'--b*')