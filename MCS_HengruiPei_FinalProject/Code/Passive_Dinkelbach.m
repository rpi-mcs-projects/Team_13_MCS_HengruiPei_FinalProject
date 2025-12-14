function [R, P] = Passive_Dinkelbach(xi,M,K,N,PBS,WBS,WUE,WRIS,sigma1,Theta,W,h_k,f_k,G)

dink = 10;
EE = 0.5*ones(dink, 1);

for d = 2:dink
    eta = EE(d-1);
    [R, P] = Passive_Precoding(eta,xi,M,K,N,PBS,WBS,WUE,WRIS,sigma1,Theta,W,h_k,f_k,G);
    EE(d) = R/P;
    obj = R-eta*P;
    if abs(obj)<0.02
        break;
    end
end

end