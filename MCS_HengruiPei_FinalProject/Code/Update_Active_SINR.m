function [R, P, obj] = Update_Active_SINR(T,K,M,N,h_k,f_k,G,Theta,w_k,sigma1,sigma2,EE,WBS,WUE,WRIS,WRA,xi,zeta)
L = N/T;
gamma_k=zeros(K,1);

for k=1:K
    h_k_temp=reshape(h_k(k,:),M,1);
    f_k_temp=reshape(f_k(k,:),N,1);
    temp1=h_k_temp+G'*Theta*f_k_temp;
    
    temp2=reshape(w_k(k,:),M,1);
    temp3=temp1'*temp2;
    
    temp4=0;
    for j=1:K
        w_j=reshape(w_k(j,:),M,1);
        temp4=temp4+norm(temp1'*w_j)^2;
    end
    temp4=temp4-norm(temp3)^2+norm(f_k_temp'*Theta')^2*sigma2+sigma1;    
    gamma_k(k)=norm(temp3)^2/temp4;
end
R=0;
for k=1:K
    R=R+log2(1+gamma_k(k));
end

%%%%% new %%%%%
temp_sum = 0;
for k = 1:K
    wk_temp = reshape(w_k(k,:),M,1);
    temp_sum = temp_sum+norm(Theta*G*wk_temp)^2;
end
theta = diag(Theta);
P = xi*norm(w_k, 'fro')^2+zeta*(temp_sum+norm(theta)^2*sigma2)+K*WUE+WBS+N*WRIS+L*WRA;
obj = R-EE*P;
%%%%% new %%%%%
end