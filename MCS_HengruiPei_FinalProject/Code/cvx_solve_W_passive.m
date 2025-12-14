function W=cvx_solve_W_passive(M,K,V,A,W,Ps_max,WBS,xi)


A=0.5*(A+A');

cvx_begin quiet
    variable W(M*K,1) complex;
    minimize((W')*A*W-real((V')*W))
    subject to
    W'*W<=(Ps_max-WBS)/xi;
cvx_end

end