function [x,t,time] = wclo_sca(A,c,Q,b0,err,x0)
%      min -t
% s.t. g_z(t,y) := (t-b0^T*y)^2-2*z^T*Q^T*y+||z||_2^2 <= 0;
%      A^T*y <= c;
%      y <= 0
% x below replace y above.

[m,n] = size(A);
MaxIter = 100;
MaxTime = 300;
tstart = tic;
%t = b0'*x0;
z = Q'*x0;
iter = 1;

while iter < MaxIter && toc(tstart) < MaxTime
    %% CVX code for test
%     cvx_begin
%     variable t
%     variable x(m)
%     cvx_quiet true
%     minimize -t
%     subject to
%     (t-b0'*x)^2-2*z'*Q'*x+z'*z <= 0;
%     A'*x <= c;
%     x <= 0;
%     cvx_end
    %%
    %% cplex code
    H = zeros(m+1);
    f = [-1,zeros(1,m)];
    Aineq = [zeros(n,1),A'];
    bineq = c;
    Aeq = [];
    beq = [];
    lb = [];
    l = [0,-2*z'*Q']';
    cQ = [1,-b0';-b0,b0*b0'];  %% Quadratic term of QC, equals to [1;-b0]*[1,-b0'];
    r = -z'*z;
    ub = [inf;zeros(m,1)];
    tx = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,cQ,r,lb,ub);
    t = tx(1);
    x = tx(2:m+1);
    %%
    zp = Q'*x;
    if norm(z-zp) < err
        time = toc(tstart);        
        return;
    else
        iter = iter+1;
        z = zp;
    end
end
