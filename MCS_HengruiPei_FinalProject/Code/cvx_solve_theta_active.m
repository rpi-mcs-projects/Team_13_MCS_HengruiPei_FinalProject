function theta = cvx_solve_theta_active(N,K,M,nu,Lam,w_k,G,PRIS,sigma2,WRIS,WRA,T,zeta)
L = N/T;

U=zeros(N,N);
for k=1:K
    w_k_temp=reshape(w_k(k,:),M,1);
    beta_k = G*w_k_temp;
    U=U+diag(abs(beta_k).^2);
end
U=U+sigma2*eye(N);

Lam=0.5*(Lam+Lam');
U=0.5*(U+U');

PRIS = (PRIS-N*WRIS-L*WRA)/zeta;

lambda=0;
% cvx_begin quiet
%     variable theta(N,1) complex;
%     maximize(-(theta')*Lam*theta+2*real((theta')*nu))
%     subject to
%     theta'*U*theta<=Pr_max;
% cvx_end    
nu = nu/2;
theta=inv(Lam+lambda*U)*nu;

lambda_left=0;
lambda_right=10;

theta_left=inv(Lam+lambda_left*U)*nu;
theta_right=inv(Lam+lambda_right*U)*nu;
P_left=theta_left'*U*theta_left;
P_right=theta_right'*U*theta_right;

if P_left<PRIS
    theta=theta_left;
    return;
end

while P_left>PRIS && P_right>PRIS
	lambda_left=lambda_right;
	lambda_right=lambda_right*10;
	theta_left=inv(Lam+lambda_left*U)*nu;
	theta_right=inv(Lam+lambda_right*U)*nu;
	P_left=theta_left'*U*theta_left;
	P_right=theta_right'*U*theta_right;
end  

while abs(lambda_left-lambda_right)>0.0001

lambda_middle=(lambda_left+lambda_right)/2;
theta_middle=inv(Lam+lambda_middle*U)*nu;
P_middle=theta_middle'*U*theta_middle;

if P_middle<PRIS
    lambda_right=lambda_middle;
    theta_right=inv(Lam+lambda_right*U)*nu;
    P_right=theta_right'*U*theta_right;
else if P_middle>PRIS
    lambda_left=lambda_middle;
	theta_left=inv(Lam+lambda_left*U)*nu;
    P_left=theta_left'*U*theta_left;
    else 
        theta=theta_middle;
        return;
end
    
end

theta=theta_right;

%%%%% new %%%%%
%theta_phase = 
%%%%% new %%%%%
end
end