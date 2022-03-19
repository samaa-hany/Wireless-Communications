%% Samaa Hany Seif Elyazal
%% Wireless Communication, Intake 42
%% Transmit Diversty without CSI (Alamouti)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial Multiplexing QPSk
%% Intialization
clc;
clear ;
close all;
%% Intialization
SNRV=-5:5:25;
N=6;
M=4;
Nb=2; % no of Bits Per QPSK SymbolNs=2; % no of QPSK Symbols
Ns=4;
NI=2000;
P=1;
x=1;
for x=1:length(SNRV)
SNR=SNRV(x);%Signal To Noise Ratio in dB 
Nerrors=0;
NerrorsSIC=0;
Nerrorsmmse=0;
iteration=0;
snr=10^(SNR/10);
N0=P/snr;

while Nerrors<800 && iteration<10000
%% Generate Symbols
I=randi([0 1],M,1);
Q=randi([0 1],M,1);
i=I;
q=Q;
S=((2*I-1)+1i*(2*Q-1))./sqrt(2);

%% Channel Estimation
H=(1/sqrt(2)).*(randn(N,M)+1i.*randn(N,M));
%% Noise
n=sqrt(N0/2)*(randn(N,1)+1i*randn(N,1));
%% Received Signal
y=H*S*sqrt(P/M)+n;
%% GZF Detector
Gzf=pinv(H);
%% ZF Decoding
S_hat=Gzf*y;
I_real=real(S_hat)>0;
Q_imag=imag(S_hat)>0;
Nerrors=Nerrors+sum(I_real~=I)+sum(Q_imag~=Q);
%% ZF SIC Deccoding
while length(i)~=M/2
norm_G=abs(Gzf).^2;
D=sum(norm_G,2);
index=find(D==min(D));
y=y-H(:,index)*S_hat(index,:)*sqrt(P/M);
Gzf(index,:)=[];
R=Gzf*y;
I_real=real(R)>0;
Q_imag=imag(R)>0;
i(index)=[];
q(index)=[];
end
NerrorsSIC=NerrorsSIC+sum(I_real~=i)+sum(Q_imag~=q);
iteration=iteration+1;
end
SNRV(x)
Nerrors
NerrorsSIC
iteration
berzf(x)=Nerrors/(Nb*Ns*iteration);
berzfSIC(x)=NerrorsSIC/(Nb*Ns*iteration);
end
%% Plotting
semilogy(SNRV,berzf,'r-*',SNRV,berzfSIC,'g-^','linewidth',2);
legend('ZF','ZF SIC');
xlabel('SNR (dB)');
ylabel('BER')
%dhat=(log(BER_AL22(end))-log(BER_AL22(end-1)))/(-log(snrv(end))+log(snrv(end-1)))
axis([min(SNRV), max(SNRV), 1e-4,1])