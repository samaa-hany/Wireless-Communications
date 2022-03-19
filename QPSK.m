%% Samaa Hany Seif Elyazal
%% Wireless Communication, Intake 42
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulatical
N=10000;
Es=1;
Real=randsrc(1,N,[-1,1]);
Imag=randsrc(1,N,[-1,1]);
c=sqrt(Es/2);
S=c*(Real+j*Imag);
SNR=[-10:2:10];
No=Es*10.^(-SNR/10);
W_real=randn(1,N);
W_imag=j*randn(1,N);
w=W_real+W_imag;
l=length(S);
for k=1:length(SNR)
    W=sqrt(No(k)/2).*w;
    Rs=W+S;
    Rs_real(find(real(Rs)>0))=c;
    Rs_imag(find(imag(Rs)>0))=c;
    Rs_real(find(real(Rs)<=0))=-c;
    Rs_imag(find(imag(Rs)<=0))=-c;
    Rf=Rs_real+j*Rs_imag;
    
    error=sum(Rf~=S);
    BER_S(k)=error/l;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theoritical
BER_th=[];
for(i=No)
    ber_th=2*Q_function(sqrt(Es/i));
    BER_th=[BER_th ber_th];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ploting
semilogy(SNR,BER_S,'r*')
title('BER vs SNR for QPSK');
hold on;
semilogy(SNR,BER_th,'bo')
legend('Experimental','Theoretical');
hold off;