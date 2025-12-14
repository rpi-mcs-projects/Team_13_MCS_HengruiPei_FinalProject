function [nu,Lam] = Update_Active_Q_v(K,M,N,Rho_k,eps_k,h_k,f_k,G,w_k,sigma2,EE,zeta)

nu_k=zeros(K,N);
Lam_k=zeros(K,N,N);

for k=1:K
    temp_f_k=reshape(f_k(k,:),N,1);
    temp_h_k=reshape(h_k(k,:),M,1);
    temp_w_k=reshape(w_k(k,:),M,1);
    
    temp1=2*sqrt(1+Rho_k(k))*eps_k(k)'*diag(temp_f_k')*G*temp_w_k;
    
    temp2=zeros(M,M);
    for j=1:K
        temp3=reshape(w_k(j,:),M,1);
        temp3=temp3*temp3';
        temp2=temp2+temp3;
    end
    
    beta_k = G*temp_w_k;
    
    nu_k(k,:)=temp1-abs(eps_k(k))^2*diag(temp_f_k')*G*temp2*temp_h_k;
    
    Lam_k(k,:,:)=abs(eps_k(k))^2*diag(abs(temp_f_k).^2)*sigma2+abs(eps_k(k))^2*diag(temp_f_k')*G*temp2*G'*diag(temp_f_k)+EE*zeta*diag(abs(beta_k).^2);
end

nu=zeros(N,1);
Lam=zeros(N,N);

for k=1:K
    nu=nu+reshape(nu_k(k,:),N,1);
    Lam=Lam+reshape(Lam_k(k,:,:),N,N);
end

Lam = Lam+EE*zeta*sigma2*eye(N);

end

