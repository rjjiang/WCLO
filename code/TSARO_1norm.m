function [v_opt,time,status] = TSARO_1norm(xi,zeta,c,eta,A,D,Q,d0,bigM)
%% Proposition 14: TSARO_1 has the same optimal value with the following LO problem
%         min       \xi'*y+\zeta'*z+w
% {y,z,w,x^1,...,x^{2r}}
%         s.t. w >= c'*x^i,             i = 1,...,2r,
%              A*x^i+D*z <= Q*e_i+d0,   i = 1,...,2r,
%               x^i >= 0,               i = 1,...,2r,
%               z <= diag(\eta)*y,
%               z >= 0,
%               y \in {0,1}^m.

tstart = tic;
MaxTime = 600;
m = length(xi);
[n,r] = size(Q);
n = n-m;
% x=[y;  z;  w; x^1,...,x^{2r}];
f = [xi;zeta;1;zeros(2*r*m*n,1)];

%z <= diag(\eta)*y,
Aineq = [-diag(eta),eye(m),zeros(m,1),zeros(m,2*r*m*n)];
bineq = [zeros(m,1)];
lb = [zeros(m,1);zeros(m,1);-inf;zeros(2*r*m*n,1)];
ub = [ones(m,1);inf*ones(m,1);inf;inf*ones(2*r*m*n,1)];
%ctype = char(m+m+1+2*r*m*n);
ctype(1:m) = 'B';           % y \in {0,1}^m.
ctype(m+1:m+m) = 'C';
ctype(m+m+1) = 'C';
ctype(m+m+2:m+m+1+2*r*m*n) = 'C';
options = cplexoptimset('MaxTime',600,'TolXInteger',1e-7,'TolRLPFun',1e-6,'TolFun',1e-6,'Display','off');
E = [eye(r),-eye(r)];
for i = 1:2*r
    % w >= c'*x^i,
    Aineq(m+i,m+m+1) = -1;
    Aineq(m+i,m+m+1+(i-1)*m*n+1:m+m+1+i*m*n) = c';
    bineq(m+i) = 0;
    
    % A*x^i+D*z <= Q*e_i+d0
    Aineq(m+2*r+1:m+2*r+m+n,1:m+m) = D;
    Aineq(m+2*r+1:m+2*r+m+n,m+m+1+(i-1)*m*n+1:m+m+1+i*m*n) = A;
    bineq(m+2*r+1:m+2*r+m+n) = Q*E(:,i)+d0;
end
[x,v_opt,exitflag,output] = cplexmilp(f,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],options);
status = output.cplexstatus;
time = toc(tstart);
    
    