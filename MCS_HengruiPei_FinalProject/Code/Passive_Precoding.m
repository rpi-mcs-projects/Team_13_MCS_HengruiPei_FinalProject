function [R, P] = Passive_Precoding(EE,xi,M,K,N,PBS,WBS,WUE,WRIS,sigma1,Theta,W,h_k,f_k,G)
iteration = 30;
obj = zeros(2*iteration, 1);
nu = ones(K, 1);
theta = diag(Theta);
for i = 1:iteration
    w_k = Update_wk(K,M,W);
    mu = Update_mu(N,M,K,h_k,f_k,G,Theta,nu,w_k); % mu
    nu = Update_Passive_nu(K,M,N,mu,h_k,f_k,G,Theta,w_k,sigma1); % nu
    [V,A] = Update_Passive_S_u(K,M,N,mu,nu,h_k,f_k,G,Theta,EE,xi);
    W = Update_W(K,M,w_k);
    W = cvx_solve_W_passive(M,K,V,A,W,PBS,WBS,xi);
    w_k = Update_wk(K,M,W);
    [~, ~, obj(2*i-1)] = Update_Passive_SINR(K,M,N,h_k,f_k,G,Theta,w_k,sigma1,EE,WBS,WUE,WRIS,xi);
    [v,Q] = Update_Passive_Q_v(K,M,N,mu,nu,h_k,f_k,G,w_k);
    theta = cvx_solve_theta_passive(N,v,Q,theta);
    Theta = diag(theta);
    [R, P, obj(2*i)] = Update_Passive_SINR(K,M,N,h_k,f_k,G,Theta,w_k,sigma1,EE,WBS,WUE,WRIS,xi);
    
    if i>1
        if (obj(2*i)-obj(2*i-1))/obj(2*i-1)<0.01
            break;
        end
    end
end
end