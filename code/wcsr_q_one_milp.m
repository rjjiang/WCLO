function [z,y,u,v,w,fval,time,status] = wcsr_q_one_milp(A,c,Q,b0)
%% Call the cplexmiqcp solver to solve wcsr problem with p=\infty. 
%
% The wcsr problem with p=\infty:
%     max        min     c'*z
% ||u||_\infty <= 1   z\in R^m 
% s.t.  A*z <= Q*u+b0
%         z >= 0.
%
% Using KKT conditions, wcsr is equivalent to the following MILP:
%    max        c'*z
%  z,y,u,v,w
%  s.t.  A*z <= Q*u+b0,
%        A'*y <= c,
%        -y <= M*w,
%        Q*u+b0-A*z <= M*(e-w), e \in R^n is the vector of all ones,
%        z <= M*v,
%        c-A'*y <= M*(e-v), e \in R^m is the vector of all ones,
%        z >= 0, z \in R^m,
%        y <= 0, y \in R^n,
%        ||u||_\infty <= 1, u \in R^r,
%        v \in {0,1}^m,
%        w \in {0,1}^n.
tstart = tic;
[n,m] = size(A);
[~,r] = size(Q);
bigM = 100000;  % M in MILP
f = [-c;zeros(n+r+m+n,1)]; % min -c'*z
Aineq = [A,zeros(n,n),-Q,zeros(n,m),zeros(n,n)]; % A*z <= Q*u+b0
Aineq = [Aineq;zeros(m,m),A',zeros(m,r),zeros(m,m),zeros(m,n)]; % A'*y <= c
Aineq = [Aineq;zeros(n,m),-eye(n),zeros(n,r),zeros(n,m),-bigM*eye(n)]; % -y <= M*w
Aineq = [Aineq;-A,zeros(n,n),Q,zeros(n,m),bigM*eye(n)]; % Q*u+b0-A*z <= M*(e-w)
Aineq = [Aineq;eye(m),zeros(m,n),zeros(m,r),-bigM*eye(m),zeros(m,n)]; % z <= M*v
Aineq = [Aineq;zeros(m,m),-A',zeros(m,r),bigM*eye(m),zeros(m,n)]; % c-A'*y <= M*(e-v)
bineq = [b0;c;zeros(n,1);bigM*ones(n,1)-b0;zeros(m,1);bigM*ones(m,1)-c];

% z >= 0, y <= 0, ||u||_\infty <= 1
lb = [zeros(m,1);-inf*ones(n,1);-ones(r,1);-inf*ones(m,1);-inf*ones(n,1)];
ub = [ones(m,1);zeros(n,1);ones(r,1);inf*ones(m,1);inf*ones(n,1)];

ctype = char(m+n+r+m+n);
ctype(1:m+n+r) = 'C'; % z \in R^m, y \in R^n, u \in R^r,
ctype(m+n+r+1:m+n+r+m+n) = 'B'; % v \in {0,1}^m, w \in {0,1}^n,

options = cplexoptimset('MaxTime',600,'TolXInteger',1e-7,'TolRLPFun',1e-6,'TolFun',1e-6);
[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],options);
status = output.cplexstatus;
z = x(1:m);
y = x(m+1:m+n);
u = x(m+n+1:m+n+r);
v = x(m+n+r+1:m+n+r+m);
w = x(m+n+r+m+1:m+n+r+m+n);
time = toc(tstart);


