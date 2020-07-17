%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AUTHOR:- HASTI DOSHI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
L=100;
k=3;

Nsim=10000;
union_upper_bound=zeros(1,16);
sum=zeros(1,16);
snr_db=0:1:15;

coordinate_real=[1/sqrt(2),1/2,  0,        -1/2, -1/sqrt(2),-1/2,0,         1/2];
coordinate_imag=[0,       1/2,   1/sqrt(2),-1/2,  0,        1/2,-1/sqrt(2),-1/2];

%avd symbol energy

sym_e_ar= (coordinate_real.*coordinate_real + coordinate_imag.*coordinate_imag )/8 ;
sym_e=0;
ds=coordinate_real(1,1);

%dis has symbol distance

dis=zeros(2,2^k-1);
for i=2:1:2^k
 dis(1,i-1)=abs(ds-coordinate_real(1,i));
end
ds=coordinate_imag(1,1);
for i=2:1:2^k
 dis(2,i-1)=abs(ds-coordinate_imag(1,i));
end

d_sym=zeros(1,2^k-1);%this has minimum symbol distance

for j=1:1:2^k-1
for i=1:1:2
    d_sym(1,j)= d_sym(1,j) + dis(i,j)*dis(i,j);
end
d_sym(1,j)=sqrt(d_sym(1,j));
end
for i=1:1:2^k
    sym_e=sym_e +sym_e_ar(1,i);
end

max=zeros(1,L);
c=1;

count=zeros(17,Nsim);

for N=1:1:Nsim

input_bit=(randi([0,1],k,L));


c=1;
while c~=17
i=1;
snr_lin=10^(snr_db(1,c)/10);
sigma=sqrt(sym_e/(2*snr_lin));
No=2*sigma*sigma;

%ylim([-2,2]);
%plot(transmitted_bit,'r.');

distance = zeros(2^k,L);

% for symbol 1  
i=1;
transmitted_bit=zeros(2,L);
transmitted_symbol=zeros(2,L);
while i~=L+1   
    if (input_bit(1,i)==0) && (input_bit(2,i)==0) &&  (input_bit(3,i)==0)
        transmitted_bit(1,i)=coordinate_real(1,1);
      transmitted_bit(2,i)=coordinate_imag(1,1);
    end
     if (input_bit(1,i)==0) && (input_bit(2,i)==0)  && (input_bit(3,i)==1)
       transmitted_bit(1,i)=coordinate_real(1,2);
      transmitted_bit(2,i)=coordinate_imag(1,2);
     end
     if (input_bit(1,i)==0) && (input_bit(2,i)==1) && (input_bit(3,i)==0)
      transmitted_bit(1,i)=coordinate_real(1,3);
      transmitted_bit(2,i)=coordinate_imag(1,3);
     end
     if (input_bit(1,i)==0) && (input_bit(2,i)==1) && (input_bit(3,i)==1)
  transmitted_bit(1,i)=coordinate_real(1,4);
      transmitted_bit(2,i)=coordinate_imag(1,4);
     end
        if (input_bit(1,i)==1) && (input_bit(2,i)==0) && (input_bit(3,i)==0)
       transmitted_bit(1,i)=coordinate_real(1,5);
      transmitted_bit(2,i)=coordinate_imag(1,5);
        end
         if (input_bit(1,i)==1) && (input_bit(2,i)==0) && (input_bit(3,i)==1)
     transmitted_bit(1,i)=coordinate_real(1,6);
      transmitted_bit(2,i)=coordinate_imag(1,6);
         end
         if (input_bit(1,i)==1) && (input_bit(2,i)==1) && (input_bit(3,i)==0)
         transmitted_bit(1,i)=coordinate_real(1,7);
      transmitted_bit(2,i)=coordinate_imag(1,7);
         end
         if (input_bit(1,i)==1) && (input_bit(2,i)==1) && (input_bit(3,i)==1)
        transmitted_bit(1,i)=coordinate_real(1,8);
      transmitted_bit(2,i)=coordinate_imag(1,8);
         end
        i=i+1;
end
   
i=1;

n_1=sigma*(randn(1,L));
n_2=sigma*(randn(1,L));

transmitted_symbol(1,:)=transmitted_bit(1,:) + n_1;
transmitted_symbol(2,:)=transmitted_bit(2,:) + n_2;

distance=zeros(2^k,L);
max=zeros(1,L);
while i~=L+1   
    
    for j=1:1:2^k
      distance(j,i)=(abs(transmitted_symbol(1,i)-coordinate_real(1,j))*abs(transmitted_symbol(1,i)-coordinate_real(1,j)))  +  (abs(transmitted_symbol(2,i)-coordinate_imag(1,j))*abs(transmitted_symbol(2,i)-coordinate_imag(1,j)))  ;
    end
    max(1,i)=distance(1,i);
    min_t=1;
    for j=2:1:2^k
        if max(1,i)>distance(j,i)
            max(1,i)=distance(j,i);
            min_t=j;
        end
    end
    transmitted_symbol(1,i)=coordinate_real(1,min_t);
     transmitted_symbol(2,i)= coordinate_imag(1,min_t);
   
    i=i+1;
end
i=1;
while i~=L+1
    if transmitted_symbol(1,i)~=transmitted_bit(1,i) || (transmitted_symbol(2,i)~=transmitted_bit(2,i))   
        count(c,N)=count(c,N) + 1;
    end
    i=i+1;
end

count(c,N)=count(c,N)/L;
sum(1,c)=sum(1,c) + count(c,N);
union_upper_bound(1,c)=2*qfunc(min(d_sym)/(2*sigma));
c=c+1;
end

end
sum=sum./Nsim;
figure(2);



semilogy(snr_db,sum,'bo:','linewidth',2,'markerfacecolor','b'); 
hold on;
semilogy(snr_db,union_upper_bound,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for 8-PSK Modulation');
ylabel('Probability of Symbol Error'); 
xlabel('E_s/N_0 in dB');
grid;
