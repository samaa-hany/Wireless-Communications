%% WIRELESS COMMUNICATION FUNDAMENTALS PROJECT
%% SAMAR ASHRAF ABDELAZIZ

%% Diversity Order BPSK
%% INTIALIZATION
clc;
clear all;
close all;
SNRv=[-5:5:50];
for n=1:length(SNRv)  
SNR=SNRv(n); %Signal to Noise Ration in dB
snr=10^(SNR/10);
snrv(n)=snr;
P=1; %Transmit Power
Ns=100000; %Number of QPSK Symbol
M=2; %M-ary for QPSK M=4;
Nb= Ns*log2(M); %Number of bits per symbol
N0=P/snr;
%% Therotical Error Rate
theortical_error(n)=qfunc(sqrt(snr));
%% Generate Symbols
ber=0;
berf=0;
berfMRC=0;
berfMRC1=0;
berfA=0;
berfA1=0;
berfA2=0;
berfA3=0;
for k=1:1:Ns
I=randi([0 1],1,1);
I_1=randi([0 1],1,1);
S=((2*I-1)).*sqrt(1/2);
S_1=((2*I_1-1)).*sqrt(1/2);
%abs(sqrt(P)*S).^2 check the transmit code
%% AWGN
W=sqrt(N0/2)*(randn(1,1)+j*randn(1,1));
W1=sqrt(N0/2)*(randn(1,1)+j*randn(1,1));
W2=sqrt(N0/2)*(randn(1,1)+j*randn(1,1));
W3=sqrt(N0/2)*(randn(1,1)+j*randn(1,1));
yawgn=(sqrt(P)*S)+W;
%% Decoding AWGN
I_hat=real(yawgn)>0;
ber=ber+((I~=I_hat));
%% Fading
h=sqrt(0.5)*(randn(1,1)+j*randn(1,1));
h1=sqrt(0.5)*(randn(1,1)+j*randn(1,1));
h2=sqrt(0.5)*(randn(1,1)+j*randn(1,1));
h3=sqrt(0.5)*(randn(1,1)+j*randn(1,1));
yfading=h*S+W;
yfading1=h1*S+W1;
yfading2=h2*S+W2;
yfading3=h3*S+W3;
%% MRC Fading 
yMRC1_2=conj(h).*yfading+conj(h1).*yfading1; %(1 TX to 1RX)
yMRC1_4=conj(h).*yfading+conj(h1).*yfading1+conj(h2).*yfading2+conj(h3).*yfading3; %(1 TX to 4 RX)
%% Alamouti Fading
yfadingA=sqrt(P/2).*S*h+sqrt(P/2).*S_1*h1+W;
yfadingA1=-sqrt(P/2).*conj(S_1)*h+sqrt(P/2).*conj(S)*h1+W1;
yfadingA2=sqrt(P/2).*S*h2+sqrt(P/2).*S_1*h3+W2;
yfadingA3=-sqrt(P/2).*conj(S_1)*h2+sqrt(P/2).*conj(S)*h3+W3;
%% Alamouti Fading (2 TX to 2 RX)
S_hatAlamouti=(conj(h).*yfadingA+h1.*conj(yfadingA1));
S_hatAlamouti1=(conj(h1).*yfadingA-h.*conj(yfadingA1));
%% Alamouti Fading (2 TX to 4 RX)
S_hatAlamouti2=(conj(h).*yfadingA+h1.*conj(yfadingA1)+conj(h2).*yfadingA2+h3.*conj(yfadingA3));
S_hatAlamouti3=(conj(h1).*yfadingA-h.*conj(yfadingA1)+conj(h3).*yfadingA2-h2.*conj(yfadingA3));
%% Decoding Fading using (Rayleigh,MRC,Alamouti)
S_hat=yfading./h;
I_hat=real(S_hat)>0;
I_hatMRC=real(yMRC1_2)>0;
I_hatMRC1=real(yMRC1_4)>0;
I_hatA=real(S_hatAlamouti)>0;
I_hatA1=real(S_hatAlamouti1)>0;
I_hatA2=real(S_hatAlamouti2)>0;
I_hatA3=real(S_hatAlamouti3)>0;
berf=berf+((I~=I_hat));
berfMRC=berfMRC+((I~=I_hatMRC));
berfMRC1=berfMRC1+((I~=I_hatMRC1));
berfA=berfA+((I~=I_hatA));
berfA1=berfA1+((I_1~=I_hatA1));
berfA2=berfA2+((I~=I_hatA2));
berfA3=berfA3+((I_1~=I_hatA3));
end
SNRv(n)
BER(n)=ber/Nb;
BERF(n)=berf/Nb;
BERFMRC(n)=berfMRC/Nb;
BERFMRC1(n)=berfMRC1/Nb;
BERFA(n)=(berfA1+berfA)/(2*Nb);
BERFA1(n)=(berfA2+berfA3)/(2*Nb);
end
%% Plotting
semilogy(SNRv,BERF,'-g',SNRv,BERFMRC,'-k',SNRv,BERFMRC1,'-p',SNRv,BERFA,'-d',SNRv,BERFA1,'-^');
legend('no diversity (1 TX,1 RX)','MRRC (1 TX,2 RX)','MRRC (1 TX,4 RX)',' MIMO(2 TX,1 RX)','MIMO(2 TX,2 RX)');
xlabel('SNR(dB)');
ylabel('BER');
axis([min(SNRv),max(SNRv),1e-4,1]);

%% Diversity Order Estimation
d_hat=(log(BERF(end))-log(BERF(end-1)))/(-log(snrv(end))+log(snrv(end-1)));