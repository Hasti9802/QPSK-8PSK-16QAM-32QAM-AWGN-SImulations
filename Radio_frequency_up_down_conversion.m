clear all;
close all;

n=100;
Lpacket=100;
p=16;
k=2;
load srrcFilter_v2;

Nsim=10000;
union_upper_bound=zeros(1,16);
ber=zeros(1,16);
snr_db=0:1:15;
isRFsimulationEnabled=1;
coordiante_real=[1,1,-1,-1];
coordinate_imag=[1,-1,1,-1];

systenm_avg = (coordiante_real.*coordiante_real + coordinate_imag.*coordinate_imag )/(2^k) ;
sym_e=0;
d_temp=coordiante_real(1,1);
dis=zeros(2,2^k-1);
for i=2:1:2^k
 dis(1,i-1)=abs(d_temp-coordiante_real(1,i));
end
d_temp=coordinate_imag(1,1);
for i=2:1:2^k
 dis(2,i-1)=abs(d_temp-coordinate_imag(1,i));
end

d_sym=zeros(1,2^k-1);
for j=1:1:2^k-1
for i=1:1:2
    d_sym(1,j)= d_sym(1,j) + dis(i,j)*dis(i,j);
end
d_sym(1,j)=sqrt(d_sym(1,j));
end
for i=1:1:2^k
    sym_e=sym_e +systenm_avg(1,i);
end

max=zeros(1,n);
c=1;



Rs = 50e3; % Symbol rate in symbols/sec (arbitrarily set)
Ts = 1/Rs; % Symbol duration in seconds Lpacket = 100;
% Packet length in symbols 
P_sampPerSym = 16; 
Fs=1/Ts;
% SRRC oversample rate in samples/symbol
Lsamp = Lpacket*P_sampPerSym + 96; % Number of samples in the packet
Fsamp = Rs*P_sampPerSym; Tsamp = 1/Fsamp;
t0 = (0:Lsamp-1)*Tsamp; t0 = t0(:); % Packet duration in seconds
fc = 3*Rs; % Carrier frequency in Hz (150 kHz in this case) 
cosineSignal =cos(2*pi*fc*t0); % Electromagnetic signals 
sineSignal =sin(2*pi*fc*t0); % Electromagnetic signals 
count=zeros(17,Nsim);
Nsamp=1696;
for N=1:1:Nsim

input_bit=(randi([0,1],k,n));


c=1;
while c~=17

i=1;



snr_lin=10^(snr_db(1,c)/10);
sigma=sqrt((16*sym_e)/(2*snr_lin));
No=2*sigma*sigma;



%ylim([-2,2]);
%plot(transmitted_bit,'r.');
distance = zeros(2^k,n);
% for symbol 1  
i=1;
transmitted_bit=zeros(2,n);
transmitted_bit_with_zeros=zeros(2,n*p);
transmitted_symbol_with_zero=zeros(2,1792);
transmitted_symbol_with_zeros=zeros(2,1696);
transmitted_bit_with_zero=zeros(2,Nsamp);
transmitted_symbol=zeros(2,n);
while i~=n+1   
     if (input_bit(1,i)==0) && (input_bit(2,i)==0) 
        transmitted_bit(1,i)=coordiante_real(1,1);
      transmitted_bit(2,i)=coordinate_imag(1,1);
    end
     if (input_bit(1,i)==0) && (input_bit(2,i)==1) 
       transmitted_bit(1,i)=coordiante_real(1,2);
      transmitted_bit(2,i)=coordinate_imag(1,2);
     end
     if (input_bit(1,i)==1) && (input_bit(2,i)==0)
      transmitted_bit(1,i)=coordiante_real(1,3);
      transmitted_bit(2,i)=coordinate_imag(1,3);
     end
     if (input_bit(1,i)==1) && (input_bit(2,i)==1) 
      transmitted_bit(1,i)=coordiante_real(1,4);
      transmitted_bit(2,i)=coordinate_imag(1,4);
     end
        i=i+1;
end
   
  i=1;
   for i=1:1:n
       transmitted_bit_with_zeros(1,16*(i-1)+1)=transmitted_bit(1,i);
          
   end
    for i=1:1:n
       transmitted_bit_with_zeros(2,16*(i-1)+1)=transmitted_bit(2,i);     
   end
n_1=sigma*(randn(1,n*p));
n_2=sigma*(randn(1,n*p));
transmitted_bit_with_zero(1,:)=conv(transmitted_bit_with_zeros(1,:),(srrcImpulseResponse_alpha03_P16(1,:)));
transmitted_bit_with_zero(2,:)=conv(transmitted_bit_with_zeros(2,:),(srrcImpulseResponse_alpha03_P16(1,:)));


s_tx_inphase=(transmitted_bit_with_zero(1,:)).';
s_tx_quadphase=transmitted_bit_with_zero(2,:).';
% The inphase and the quadrature branches' output after the pulse shaping at the transmitter are denoted as s_tx_inphase and s_tx_quadphase, respectively. The following is the baseband complex-envelope of the transmitted signal
s_tx = s_tx_inphase + 1i*s_tx_quadphase;

if isRFsimulationEnabled % if this flag is set to 1, the RF part is enabled
    
    % Here we upconvert the baseband modulated pulse-shaped signal to radio frequency (RF)
    % cosSignal = cos(2*pi*fc*n/Fs) and sinSignal = sin(2*pi*fc*n/Fs) need
    % to be defined as described in the Lab Manual
    
    % Implement the quadrature notation of the transmitted signal as defined in
    % the lecture slides; note s_tx_RF is not complex-valued anymore
    % Verify mathematically that sqrt(2) multiplier is needed to preserve the symbol
    % energy Es
    
    s_tx_RF = sqrt(2)*(s_tx_inphase.*cosineSignal + s_tx_quadphase.*sineSignal);
    
    % Observe the power spectral density of the transmitted signal
    
   % [P_tx,f_tx] = pwelch(s_tx_RF,[],[],[],Fs,'twosided');
    %plot(f_tx-Fs/2,10*log10(fftshift(P_tx)),'linewidth',2);
   % xlabel('Frequency in Hertz'); ylabel('dB'); 
    %title('Power Spectral Density of the Transmitted Signal');
    
    % Add the AWGN at RF. Note that the AWGN is real-valued at the RF
    % noiseSTD has the same value as in the baseband simulator with pulse-shaping
    
    nsig_RF = sigma*randn(Nsamp,1);
    r_RF = s_tx_RF + nsig_RF;
    
    % Observe the power spectral density of the received signal in the presence of the AWGN
    
   % P_rx = pwelch(r_RF,[],[],[],Fs,'twosided');
    %figure; plot(f_tx-Fs/2,10*log10(fftshift(P_rx)),'linewidth',2);
 %   xlabel('Frequency in Hertz'); ylabel('dB'); 
 %   title('Power Spectral Density of the Received Signal in the AWGN');
       
    % Downconvert the received signal at RF to obtain the baseband
    % received signal
    
    r_inphase = r_RF.*cosineSignal;
    r_quadphase = r_RF.*sineSignal;
    
    r = r_inphase - 1i*r_quadphase;
   
else % the RF simulation is disabled, the baseband simulation with pulse-shaping is enabled
    
    nsig = sigma*(randn(Nsamp,1)+1i*randn(Nsamp,1));
    r = s_tx + nsig;
    
end

transmitted_symbol_with_zeros(1,:)=r_inphase.';
transmitted_symbol_with_zeros(2,:)=r_quadphase.';

%transmitted_symbol_with_zeros(1,:)=transmitted_bit_with_zeros(1,:) + n_1;
%transmitted_symbol_with_zeros(2,:)=transmitted_bit_with_zeros(2,:) + n_2 ;
transmitted_symbol_with_zero(1,:)=conv(transmitted_symbol_with_zeros(1,:),(srrcImpulseResponse_alpha03_P16(1,:)));
transmitted_symbol_with_zero(2,:)=conv(transmitted_symbol_with_zeros(2,:),(srrcImpulseResponse_alpha03_P16(1,:)));
  for i=1:1:n
       transmitted_symbol(1,i)=transmitted_symbol_with_zero(1,16*(i-1)+ 96+1);
          
   end
    for i=1:1:n
       transmitted_symbol(2,i)=transmitted_symbol_with_zero(2,16*(i-1)+96+1);     
   end
distance=zeros(2^k,n);
max=zeros(1,n);
i=1;
while i~=n+1   
    
    for j=1:1:2^k
      distance(j,i)=(abs(transmitted_symbol(1,i)-coordiante_real(1,j))*abs(transmitted_symbol(1,i)-coordiante_real(1,j)))  +  (abs(transmitted_symbol(2,i)-coordinate_imag(1,j))*abs(transmitted_symbol(2,i)-coordinate_imag(1,j)))  ;
    end
    max(1,i)=distance(1,i);
    min_t=1;
    for j=2:1:2^k
        if max(1,i)>distance(j,i)
            max(1,i)=distance(j,i);
            min_t=j;
        end
    end
    transmitted_symbol(1,i)=coordiante_real(1,min_t);
     transmitted_symbol(2,i)= coordinate_imag(1,min_t);
   
    i=i+1;
end
i=1;
while i~=n+1
    if transmitted_symbol(1,i)~=transmitted_bit(1,i) || (transmitted_symbol(2,i)~=transmitted_bit(2,i))   
        count(c,N)=count(c,N) + 1;
    end
    i=i+1;
end

count(c,N)=count(c,N)/n;
ber(1,c)=ber(1,c) + count(c,N);
union_upper_bound(1,c)=2*qfunc(min(d_sym)/(2*sigma/4));
c=c+1;
end

end
ber=ber./Nsim;
figure(2);



semilogy(snr_db,ber,'bo:','linewidth',2,'markerfacecolor','b'); 
hold on;
semilogy(snr_db,union_upper_bound,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error probablity for Radio Frequency Up-conversion and down conversion for QPSK(k=2)');
ylabel('Probability of Symbol Error'); 
xlabel('E_s/N_0 in dB');
grid 