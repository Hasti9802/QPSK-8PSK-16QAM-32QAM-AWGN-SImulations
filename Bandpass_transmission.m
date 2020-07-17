clearvars; 
close all;
clc;

k=2;
M=2^k;
l=100;

coordinates=[1+0j,0+1j,-1+0j,0-1j];
origin=0+0j;

Es= (sum(abs(sqrt(coordinates.*coordinates)))/M);

SNRdb = 0:15;
SNRlin=10.^(SNRdb/10);
variance=16*0.5*(Es./SNRlin);
sigma=sqrt(variance);
Nsim=10000;
% load srrcFilter;
load srrcFilter_v2;

symbol_error_rate=zeros(1,length(SNRdb));
P=16;

for i=1:16
    error=0;
    for x=1:Nsim
        my_sig= randi(2,k,l)-1;
        my_final = signal_generator(my_sig,l,coordinates);
        after_zero = z_inserter(my_final,P); % Inserts P-1 zeros after every sample to increase the sample rate for pulse shaping
        after_nois = conv(after_zero,srrcImpulseResponse_alpha03_P16,'same');
%         eyediagram(transpose(final_trans),16);
        rcvd=noise(after_nois,l,sigma(i));
        final_rcvd = conv(rcvd,srrcImpulseResponse_alpha03_P16,'same');
%         eyediagram(transpose(final_rcvd),16);
        answ = zero_rem(final_rcvd,l);% Selects every Pth data point from the given stream of 1600 data points
        decoded=demodulated(answ,coordinates,l,M);
        error=error+ser_calc(decoded,l,my_final);
    end
    symbol_error_rate(i)=error/(l*Nsim);
end




figure(1);
semilogy(SNRdb,symbol_error_rate,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
es_no_db = 0:0.1:15;
es_no_lin = 10.^(es_no_db/10);
p=2*qfunc(sqrt(es_no_lin));

semilogy(es_no_db,p,'r','linewidth',2);

legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for QPSK Modulation');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); 
set(gca,'xtick',0:1:15);
grid on;



function [euclid_dist] = distance(a,b)
    dist = sqrt((real(a)-real(b))*(real(a)-real(b))+(imag(a)-imag(b))*(imag(a)-imag(b)));
    euclid_dist = dist;
end

function [final_return] = signal_generator(input,l,cordinates)
    final_return = zeros(1,l);
    for i=1:l
        if(input(1,i)==0 && input(2,i)==0)
            final_return(i)=cordinates(1);
        elseif(input(1,i)==0 && input(2,i)==1)
            final_return(i)=cordinates(2);
        elseif(input(1,i)==1 && input(2,i)==0)
            final_return(i)=cordinates(3);
        elseif(input(1,i)==1 && input(2,i)==1)
            final_return(i)=cordinates(4);
        end
    end
end

function [returned_sig] = z_inserter(my_fin,P)
    returned_sig = zeros(1,length(my_fin)*P);
    for i=1:length(my_fin)
        returned_sig(P*(i-1) + 1) = my_fin(i); 
    end
end


function [error] = ser_calc(decoded,l,final_sig)
    error=0;
    for i=1:l
        if decoded(i) == final_sig(i)
            error = error + 0;
        else
            error=error + 1;
        end
    end
end
function [rcvd] = noise(final_sig,l,sigma)
    noise = genearte_noise(sigma,l);
    rcvd = final_sig + noise;
end

function [resultant] = zero_rem(final_rcvd,l)
    resultant = zeros(1,l);
    for i =1: length(final_rcvd)
        if mod(i,16) == 1
            resultant(ceil(i/16)) = final_rcvd(i);
        end
    end
end

function [demodulated_sig] = demodulated(resultant,points,l,M)
    demodulated_sig=zeros(1,l);
    dist = zeros(1,4);
    for i=1:l
        for k=1:M
            dist(k) = distance(resultant(i),points(k));
        end
        d = min(dist);
        if d == dist(1)
            demodulated_sig(i) = points(1);
        elseif d == dist(2)
            demodulated_sig(i) = points(2);
        elseif d == dist(3)
            demodulated_sig(i) = points(3);
        else
            demodulated_sig(i) = points(4);
        end
    end
end

function [noise] = genearte_noise(sigma,l)
    a=sigma.*randn(1,l*16);
    b=sigma.*randn(1,l*16);
    noise=complex(a,b);
end