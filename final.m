%%  Samaa Hany Seif Elyazal            %%%%%%%%%%%%%%%%%
%%  Wireless Communication, Intake 42  %%%%%%%%%%%%%%%%%
%%  Diversity Order BPSK               %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear
clc;
%% Intialization
SNRV =-5:5:25; %SIGNAL TO NOISE RATIO IN DB
for(K=1:length(SNRV))
SNR = SNRV(K);
snr = 10^(SNR / 10);
snrv(K) = snr;
P=1;
N0 = P / snr;
M = 2; %bits per sympol OF QSPK
NS = 10000;
NB = NS*log2(M); %NO OF BITS PER SYMBOLS
 %NO OF QPSK SYMBOLS

%% Therotical Error Rate
theortical_error(K)=Q_function(sqrt(snr));

%% Generate Symbols
ber=0;
ber_SIC=0;
bermmse=0;
bermmse_SIC=0;
ber=0;
berf=0;
ber_MRC0=0;
ber_MRC1=0;
ber_A0=0;
ber_A1=0;
ber_A2=0;
ber_A3=0;
for k=1:1:NS
I1=randi([0 1],1,1); I2=randi([0 1],1,1);
S1=((2*I1-1))*sqrt(1/2); S2=((2*I2-1))*sqrt(1/2);
%% AWGN
W0=sqrt(N0/2)*(randn(1,1)+1i*randn(1,1));
W1=sqrt(N0/2)*(randn(1,1)+1i*randn(1,1));
W2=sqrt(N0/2)*(randn(1,1)+1i*randn(1,1));
W3=sqrt(N0/2)*(randn(1,1)+1i*randn(1,1));
y_AWGN=(sqrt(P)*S1)+W0;
%% Decoding AWGN
I_HAT=real(y_AWGN)>0;
ber=ber+((I1~=I_HAT));
%% Fading
h0=sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
h1=sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
h2=sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
h3=sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
y_fading0=h0*S1+W0;
y_fading1=h1*S1+W1;
y_fading2=h2*S1+W2;
y_fading3=h3*S1+W3;
%% MRC Fading 
y_MRC1_2=conj(h0)*y_fading0+conj(h1)*y_fading1; %(1 TX to 1RX)
y_MRC1_4=conj(h0)*y_fading0+conj(h1)*y_fading1+conj(h2)*y_fading2+conj(h3).*y_fading3; %(1 TX to 4 RX)
%% Alamouti Fading
y_fadingA0=sqrt(P/2)*S1*h0+sqrt(P/2)*S2*h1+W0;
y_fadingA1=-sqrt(P/2)*conj(S2)*h0+sqrt(P/2)*conj(S1)*h1+W1;
y_fadingA2=sqrt(P/2)*S1*h2+sqrt(P/2)*S2*h3+W2;
y_fadingA3=-sqrt(P/2)*conj(S2)*h2+sqrt(P/2)*conj(S1)*h3+W3;
%% Alamouti Fading (2 TX to 2 RX)
S_HAT_Alamouti0=(conj(h0)*y_fadingA0+h1*conj(y_fadingA1));
S_HAT_Alamouti1=(conj(h1)*y_fadingA0-h0*conj(y_fadingA1));
%% Alamouti Fading (2 TX to 4 RX)
S_HAT_Alamouti2=(conj(h0)*y_fadingA0+h1*conj(y_fadingA1)+conj(h2)*y_fadingA2+h3*conj(y_fadingA3));
S_HAT_Alamouti3=(conj(h1)*y_fadingA0-h0*conj(y_fadingA1)+conj(h3)*y_fadingA2-h2*conj(y_fadingA3));
%% Decoding Fading using (Rayleigh,MRC,Alamouti)
S_HAT=y_fading0/h0;
I_HAT=real(S_HAT)>0;
I_HAT_MRC=real(y_MRC1_2)>0;
I_HAT_MRC1=real(y_MRC1_4)>0;
I_HAT_A0=real(S_HAT_Alamouti0)>0;
I_HAT_A1=real(S_HAT_Alamouti1)>0;
I_HAT_A2=real(S_HAT_Alamouti2)>0;
I_HAT_A3=real(S_HAT_Alamouti3)>0;
berf=berf+((I1~=I_HAT));
ber_MRC0=ber_MRC0+((I1~=I_HAT_MRC));
ber_MRC1=ber_MRC1+((I1~=I_HAT_MRC1));
ber_A0=ber_A0+((I1~=I_HAT_A0));
ber_A1=ber_A1+((I2~=I_HAT_A1));
ber_A2=ber_A2+((I1~=I_HAT_A2));
ber_A3=ber_A3+((I2~=I_HAT_A3));
end
BER(K)=ber/NB;
BERF(K)=berf/NB;
BER_MRC0(K)=ber_MRC0/NB;
BER_MRC1(K)=ber_MRC1/NB;
BER_A0(K)=(ber_A1+ber_A0)/(2*NB);
BER_A1(K)=(ber_A2+ber_A3)/(2*NB);
end
%% Plotting
semilogy(SNRV,BER,'black-o',SNRV,BER_MRC0,'m-*',SNRV,BER_MRC1,'r-s',SNRV,BER_A0,'b-o',SNRV,BER_A1,'g-^','linewidth',2);
legend('no diversity (1 TX,1 RX)','MRRC (1 TX,2 RX)','MRRC (1 TX,4 RX)',' MIMO(2 TX,1 RX)','MIMO(2 TX,2 RX)');
xlabel('SNR(dB)');
ylabel('BER');
axis([min(SNRV),max(SNRV),1e-4,1]);


%% Diversity Order Estimation
d_hat=(log(BERF(end))-log(BERF(end-1)))/(-log(snrv(end))+log(snrv(end-1)));





