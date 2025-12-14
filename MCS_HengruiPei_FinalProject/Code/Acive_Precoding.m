function [R, P] = Acive_Precoding(T,C,C_mat,EE,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G)
iteration = 30;
obj = zeros(2*iteration, 1);
nu = ones(K, 1);
for i = 1:iteration
    w_k = Update_wk(K,M,W);
    mu = Update_mu(N,M,K,h_k,f_k,G,Theta,nu,w_k); % mu
    nu = Update_Active_nu(K,M,N,mu,h_k,f_k,G,Theta,w_k,sigma1,sigma2); % nu
    [u,S,T_mat] = Update_Active_S_u(K,M,N,mu,nu,h_k,f_k,G,Theta,EE,xi,zeta);
    W = Update_W(K,M,w_k);
    W = cvx_solve_W_active(M,T,N,K,Theta,u,S,W,PBS,PRIS,sigma2,T_mat,WBS,WRIS,WRA,xi,zeta);
    w_k = Update_wk(K,M,W);
    [~, ~, obj(2*i-1)] = Update_Active_SINR(T,K,M,N,h_k,f_k,G,Theta,w_k,sigma1,sigma2,EE,WBS,WUE,WRIS,WRA,xi,zeta);
    [v,Q] = Update_Active_Q_v(K,M,N,mu,nu,h_k,f_k,G,w_k,sigma2,EE,zeta);
    theta = cvx_solve_theta_active(N,K,M,v,Q,w_k,G,PRIS,sigma2,WRIS,WRA,T,zeta);
    theta_phase = angle(theta);
    theta_amp = abs(theta); theta_amp = C_mat*theta_amp;
    theta = C*theta_amp.*exp(1j*theta_phase);
    Theta = diag(theta);
    [R, P, obj(2*i)] = Update_Active_SINR(T,K,M,N,h_k,f_k,G,Theta,w_k,sigma1,sigma2,EE,WBS,WUE,WRIS,WRA,xi,zeta);
    
    if i>1
        if (obj(2*i)-obj(2*i-1))/obj(2*i-1)<0.01
            break;
        end
    end
end
end