function W=cvx_solve_W_active(M,T,N,K,Theta,V,A,W,Ps_max,Pr_max,sigmar2,T_mat,WBS,WRIS,WRA,xi,zeta)
L = N/T;
theta=diag(Theta);

D=kron(eye(K),T_mat);
D=0.5*(D+D');

A=0.5*(A+A');

cvx_begin quiet
    variable W(M*K,1) complex;
    minimize((W')*A*W-real((V')*W))
    subject to
    W'*W<=(Ps_max-WBS)/xi;
    W'*D*W<=(Pr_max-norm(theta)^2*sigmar2-N*WRIS-L*WRA)/zeta;

cvx_end

end