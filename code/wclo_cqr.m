function [t,y,s,z,v] = wclo_cqr(A,c,Q,b0,zlb,zub)
%% This m-file calls cplexqcp to solve the Convex Quadratic Relaxation problem (15):
%        min -t
%   s.t. (t-b0^T*y)^2-sum{s_i} <= 0
%        Q^T*y = z;
%        A^T*y <= c;
%        y <= 0;
%        zlb <= z <= zub;
%        z_i^2 <= s_i, i = 1,...,r;
%        s_i <= (zlb_i+zub_i)*z_i-zlb_i*zub_i, i = 1,...,r;
%
% x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,x0,options) solve:
%  min  1/2*x'*H*x + f*x
%  s.t. Aineq*x <= bineq
%       Aeq*x = beq;
%       l(i,:)'*x + x'*Q{i}*x <= r(i); i = 1:m
%       lb <= x <= ub
% Note: m is the number of qudratic constraints;
%       x0 (starting point) must be n-by-1 matrix as x; 
%       options replace the default optimization options;
%       l must be m-by-n matrix;
%       Q must be cell(1,m) with cell{i} be n-by-n PSD matrix;
% See also \IBM\ILOG\CPLEX_Studio126\cplex\examples\src\matlab\qcpdual.m

[m,n] = size(A);
[~,r] = size(Q);
% cvx code for test:
% cvx_begin
% variable t;
% variable y(m);
% variable s(r);
% variable z(r);
% cvx_quiet true
% minimize -t
% subject to
% (t-b0'*y)^2-sum(s) <= 0;
% Q'*y == z;
% A'*y <= c;
% y <= 0;
% z >= zlb;
% z <= zub;
% for i = 1:r
%     z(i)^2 <= s(i);
%     s(i) <= (zlb(i)+zub(i))*z(i)-zlb(i)*zub(i);
% end
% cvx_end
% v = norm(Q'*y)+b0'*y;


%% for test
% Q'*y - z
% (t-b0'*y)^2-sum(s)
% A'*y - c
% y
% zlb-z
% z-zub
% z.^2-s
% s-(zlb+zub).*z+zlb.*zub

%% cplex code
H = zeros(1+m+r+r);
f = [-1,zeros(1,m+r+r)];
Aeq = [zeros(r,1),Q',zeros(r),-eye(r)];
beq = zeros(r,1);
Aineq = [zeros(n,1),A',zeros(n,r),zeros(n,r)];
bineq = c;
for i = 1:r
    Aineq(n+i,1+m+i) = 1;
    Aineq(n+i,1+m+r+i) = -zlb(i)-zub(i);
    bineq(n+i) = -zlb(i)*zub(i);
end
lb = [-inf*ones(1+m+r,1);zlb];
ub = [inf;zeros(m,1);inf*ones(r,1);zub];
cl = zeros(1+m+r+r,1+r);
cQ = cell(1,r+1);
cr = zeros(1,1+r);
cl(:,1) = [zeros(1+m,1);-ones(r,1);zeros(r,1)];
tmp = [1,-b0';-b0,b0*b0'];
tmp = [tmp,zeros(1+m,r+r);zeros(r+r,1+m),zeros(r+r)];
cQ{1} = tmp;
cr(1) = 0;
for i = 1:r
    cl(:,1+i) = zeros(1+m+r+r,1);
    cl(1+m+i,1+i) = -1;
    cQ{1+i} = zeros(1+m+r+r);
    cQ{1+i}(1+m+r+i,1+m+r+i) = 1;
    cr(1+i) = 0;
end
H = sparse(H);
f = sparse(f);
Aeq = sparse(Aeq);
beq = sparse(beq);
Aineq = sparse(Aineq);
bineq = sparse(bineq);
cl = sparse(cl);
cr = sparse(cr);
for i = 1:r
    cQ{1+i} = sparse(cQ{1+i});
end
[tysz,~,exitflag,output] = cplexqcp(H,f,Aineq,bineq,Aeq,beq,cl,cQ,cr,lb,ub);
% output.cplexstatus
% exitflag
if exitflag == 1 | exitflag == 5
    t = tysz(1);
    y = tysz(2:m+1);
    s = tysz(m+2:1+m+r);
    z = tysz(1+m+r+1:end);
    v = norm(Q'*y)+b0'*y;
elseif exitflag == -2
    t = [];
    y = [];
    s = [];
    z = [];
    v = [];
else
    disp('Error in solving convex ralexation problem!')
end

