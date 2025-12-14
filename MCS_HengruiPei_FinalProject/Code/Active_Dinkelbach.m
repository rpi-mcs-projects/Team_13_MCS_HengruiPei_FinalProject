function [R, P] = Active_Dinkelbach(T,C,C_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G)

dink = 10;
EE = 0.5*ones(dink, 1);

for d = 2:dink
    eta = EE(d-1);
    [R, P] = Acive_Precoding(T,C,C_mat,eta,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
    EE(d) = R/P;
    obj = R-eta*P;
    if abs(obj)<0.02
        break;
    end
end

end