%%  Samaa Hany Seif Elyazal            %%%%%%%%%%%%%%%%%
%%  Wireless Communication, Intake 42  %%%%%%%%%%%%%%%%%
%%  Spathial Multipexing               %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear ;
close all;
%% Intialization
SNRv=-5:5:25;
N=6;
M=4;
Nb=2; % no of Bits Per QPSK SymbolNs=2; % no of QPSK Symbols
Ns=4;
NI=2000;
P=1;
for x=1:length(SNRv)
SNR=SNRv(x);%Signal To Noise Ratio in dB 
Nerrors=0;
NerrorsSIC=0;
Nerrorsmmse=0;
NerrorsmmseSIC=0;
iteration=0;
snr=10^(SNR/10);
N0=P/snr;

while Nerrors<800 && iteration<10000
%% Generate Symbols
I=randi([0 1],M,1);
Q=randi([0 1],M,1);
S=((2.*I-1)+1i.*(2.*Q-1))./sqrt(2);

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
r=y;
h=H;
i=I;
q=Q;
R_hat=S_hat;
while length(i)~=1
norm_G=abs(Gzf).^2;
D=sum(norm_G,2);
index=find(D==min(D));
r=r-h(:,index)*R_hat(index,:)*sqrt(P/M);
Gzf(index,:)=[];
R=Gzf*r;
I_real=real(R)>0;
Q_imag=imag(R)>0;
i(index)=[];
q(index)=[];
h(:,index)=[];
R_hat(index,:)=[];
end
NerrorsSIC=NerrorsSIC+sum(I_real~=i)+sum(Q_imag~=q);
%% GMMSE Detector
Gmmse=inv(N0*eye(M)+(P/M)*(H'*H))*H';
%% MMSE Decoding
S_hat=Gmmse*y;
I_real=real(S_hat)>0;
Q_imag=imag(S_hat)>0;
Nerrorsmmse=Nerrorsmmse+sum(I_real~=I)+sum(Q_imag~=Q);
%% MMSE SIC Deccoding
r=y;
h=H;
i=I;
q=Q;
R_hat=S_hat;
while length(i)~=1  
norm_G=abs(Gmmse).^2;
D=sum(norm_G,2);
index=find(D==min(D));
r=r-h(:,index)*R_hat(index,:)*sqrt(P/M);
Gmmse(index,:)=[];
R=Gmmse*r;
I_real=real(R)>0;
Q_imag=imag(R)>0;
i(index)=[];
q(index)=[];
h(:,index)=[];
R_hat(index,:)=[];
end
NerrorsmmseSIC=NerrorsmmseSIC+sum(I_real~=i)+sum(Q_imag~=q);
iteration=iteration+1;
end
SNRv(x)
Nerrors
NerrorsSIC
%iteration
berzf(x)=Nerrors/(Nb*Ns*iteration);
berzfSIC(x)=NerrorsSIC/(Nb*Ns*iteration);
bermmse(x)=Nerrorsmmse/(Nb*Ns*iteration);
bermmseSIC(x)=NerrorsmmseSIC/(Nb*Ns*iteration);
theoryBer_nRx1(x) = 0.5.*(1-1*(1+1./snr).^(-0.5)); 
end
%% Plotting
semilogy(SNRv,berzf,'r-*',SNRv,berzfSIC,'g-o',SNRv,theoryBer_nRx1,'b-^',SNRv,bermmse,'m-*',SNRv,bermmseSIC,'y--','linewidth',2);
legend('Zero Forcing','Zero Forcing with SIC','Theortical','MMSE','MMSE with SIC');
xlabel('SNR (dB)');
ylabel('BER')
axis([min(SNRv),max(SNRv),1e-4,1]);