function [F] = Generate_Channel_F(N,M,dis,large_fading,f_c)
% N number of receiver
% M number of transmitter
f_c=f_c*10^9;
Lambda=3*10^8/f_c;

F = zeros(N,M);
for aa=1:N
    for bb=1:M
       F(aa,bb) =( randn()+1j*randn())/sqrt(2);
    end
end
A=Lambda^2/4;
a=sqrt(A/4/pi)/dis;

% a = 10^(-3.53)/(dis^large_fading);
% a= sqrt(a);

F = a*F;
end

