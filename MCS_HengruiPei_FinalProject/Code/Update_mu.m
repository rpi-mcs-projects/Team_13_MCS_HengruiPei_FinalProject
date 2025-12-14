function [temp1, Rho_k] = Update_mu(N,M,K,h_k,f_k,G,Theta,eps_k,w_k)
    Rho_k=zeros(K,1);
    for k=1:K
        h_k_temp=reshape(h_k(k,:),M,1);
        f_k_temp=reshape(f_k(k,:),N,1);
        temp1=h_k_temp+G'*Theta*f_k_temp;
        
        temp2=reshape(w_k(k,:),M,1);
        temp=real(eps_k(k)'*temp1'*temp2);
        Rho_k(k)=(temp^2+temp*sqrt(temp^2+4))/2;
    end
    
end

