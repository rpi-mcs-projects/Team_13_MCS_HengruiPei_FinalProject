function theta = cvx_solve_theta_passive(N,nu,Lam,theta)

Lam=0.5*(Lam+Lam');

cvx_begin quiet
    variable theta(N,1) complex;
    minimize ((theta')*Lam*theta-real((theta')*nu))
    subject to
    for n = 1:N
        abs(theta(n))<=1;
    end
cvx_end

end