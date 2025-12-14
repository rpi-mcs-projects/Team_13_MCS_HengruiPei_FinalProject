function [V,A,T_mat]=Update_Active_S_u(K,M,N,Rho_k,eps_k,h_k,f_k,G,Theta,EE,xi,zeta)
v_k=zeros(K,M);
A_k=zeros(K,M,M);
    
for k=1:K
    h_k_temp=reshape(h_k(k,:),M,1);
    f_k_temp=reshape(f_k(k,:),N,1);
    temp1=h_k_temp+G'*Theta*f_k_temp;
    
    v_k(k,:)=2*sqrt(1+Rho_k(k))*eps_k(k)*temp1;
    
    A_k(k,:,:)=abs(eps_k(k))^2*(temp1*temp1');
end    

V=zeros(K*M,1);

for k=1:K
    V(M*(k-1)+1:1:M*k)=reshape(v_k(k,:),M,1);
end

temp=0;
for k=1:K
    temp=temp+A_k(k,:,:);
end

temp=reshape(temp,M,M);

%%%%% new %%%%%
theta = diag(Theta);
T_mat = G'*diag(abs(theta).^2)*G;
temp = temp+EE*xi*eye(M)+EE*zeta*T_mat;
%%%%% new %%%%%

A=kron(eye(K),temp);

end

