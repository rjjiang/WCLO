function [z,y,u,v,w,fval,time,status] = wcsr_q_2_miqcp(A,c,Q,b0)
%% Call the cplexmiqcp solver to solve wcsr problem with p=2. 
%
% The wcsr problem with p=2:
%     max        min     c'*z
% ||u||_2 <= 1   z\in R^m 
% s.t.  A*z <= Q*u+b0
%         z >= 0.
%
% Using KKT conditions, wcsr is equivalent to the following MIQCP:
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
%        ||u||_2 <= 1, u \in R^r,
%        v \in {0,1}^m,
%        w \in {0,1}^n.
tstart = tic;
[n,m] = size(A);
[~,r] = size(Q);
bigM = 100000; % M in MIQCP
H = zeros(m+n+r+m+n);
f = [-c;zeros(n+r+m+n,1)]; % min -c'*z
Aineq = [A,zeros(n,n),-Q,zeros(n,m),zeros(n,n)]; % A*z <= Q*u+b0
Aineq = [Aineq;zeros(m,m),A',zeros(m,r),zeros(m,m),zeros(m,n)]; % A'*y <= c
Aineq = [Aineq;zeros(n,m),-eye(n),zeros(n,r),zeros(n,m),-bigM*eye(n)]; % -y <= M*w
Aineq = [Aineq;-A,zeros(n,n),Q,zeros(n,m),bigM*eye(n)]; % Q*u+b0-A*z <= M*(e-w)
Aineq = [Aineq;eye(m),zeros(m,n),zeros(m,r),-bigM*eye(m),zeros(m,n)]; % z <= M*v
Aineq = [Aineq;zeros(m,m),-A',zeros(m,r),bigM*eye(m),zeros(m,n)]; % c-A'*y <= M*(e-v)
bineq = [b0;c;zeros(n,1);bigM*ones(n,1)-b0;zeros(m,1);bigM*ones(m,1)-c];

% z >= 0, y <= 0, -e <= u <= e according to ||u||_2<=1.
lb = [zeros(m,1);-inf*ones(n,1);-ones(r,1);-inf*ones(m,1);-inf*ones(n,1)];
ub = [ones(m,1);zeros(n,1);ones(r,1);inf*ones(m,1);inf*ones(n,1)];

ctype = char(m+n+r+m+n);
ctype(1:m+n+r) = 'C'; % z \in R^m, y \in R^n, u \in R^r,
ctype(m+n+r+1:m+n+r+m+n) = 'B'; % v \in {0,1}^m, w \in {0,1}^n,

% ||u||_2<=1,
cl = zeros(m+n+r+m+n,1);
cQ = zeros(m+n+r+m+n);
cQ(m+n+1:m+n+r,m+n+1:m+n+r) = eye(r);
cr = 1;

options = cplexoptimset('MaxTime',600,'TolXInteger',1e-7,'TolRLPFun',1e-6,'TolFun',1e-6);
[x,fval,exitflag,output] = cplexmiqcp(H,f,Aineq,bineq,[],[],cl,cQ,cr,[],[],[],lb,ub,ctype,[],options);
status = output.cplexstatus;
z = x(1:m);
y = x(m+1:m+n);
u = x(m+n+1:m+n+r);
v = x(m+n+r+1:m+n+r+m);
w = x(m+n+r+m+1:m+n+r+m+n);
time = toc(tstart);