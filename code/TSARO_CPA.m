function [LB,UB,time,iter] = TSARO_CPA(xi,zeta,c,eta,A,D,Q,d0,bigM,Utype)
%% The Cutting Plane Algorithm (CPA), which is Algorithm 5 in the paper, is 
%% designed to solve the following problem named TSARO:
%
%   min  {\xi'*y+\zeta'*z +    max       min       Tr(CX')}
%  {y,z}                    d\in D_p   X\in S(d,z)
%   s.t.  z <= diag(\eta)*y,
%         y \in {0,1}^m,
%         z >= 0.
% where D_p = {d=Qu+d0:||u||_p<=1,u\in R^r} with p = 2,\infty denotes the
% uncertainty set based on the l_p-norm,
% S(d,z) = {X\in R_+^m*n: X*\bar e <= z, X'*\hat e >= d}.

% bigM is similar to \bar M in equation (30) to tackle the unboundness of set C.
% Utype : two  ---- p=2
%         inf  ---- p=\infty 

tstart = tic;
MaxTime = 600;
m = length(xi);
[n,r] = size(Q);
n = n-m;
%bigM = (n+m)*1000; % According to the constraints A'*y<=c,y<=0
k = 0;
gap2 = inf;
% max xi'*y +zeta'*z
% s.t. z <= diag(eta)y,
%      y_i binary, i = 1,...,m
%      z >= 0.
% f = -[xi;zeta];
% Aineq = [-diag(eta),eye(m)];
% bineq = zeros(m,1);
% lb = [-inf*ones(m,1);zeros(m,1)];
% ub = [inf*ones(m,1);inf*ones(m,1)];
% ctype = char(m+m);
% ctype(1:m) = 'B';
% ctype(m+1:m+m) = 'C';
% options = cplexoptimset('MaxTime',600,'TolXInteger',1e-7,'TolRLPFun',1e-6,'TolFun',1e-6);
% [z,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],options);
z = [ones(m,1);eta]; % Initial feasible z^0 used in step 0 (i) in Alg. 5
if Utype == 'inf'
    [y_opt,v_opt,time,gap,iter] = wclo_q_one([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z,1e-5,1);
    %[y_opt,v_opt,time,gap,iter] = wclo_q_one([A,-ones(m+n,1)],[c;100],Q,d0-D*z,1e-5,1)

elseif Utype == 'two'
    [y_opt,v_opt,scobound,time,gap,iter,nslo,iscot] = wcsr_sca_cqr([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z,1e-5,1,'cqr');
    v_opt = -v_opt;
end
% [~,y_opt,~,~,~,v_opt,time,status] = wcsr_q_one_milp([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z)
% v_opt = -v_opt;
%[~,y_opt,~,~,~,v_opt,time,status] = wcsr_q_one_milp(A,c,Q,d0-D*z)

UB = [xi;zeta]'*z+v_opt;

f = [xi;zeta;1];
Aineq = [-diag(eta),eye(m),zeros(m,1)];
bineq = zeros(m,1);
lb = [zeros(m,1);zeros(m,1);-inf];
ub = [ones(m,1);inf*ones(m,1);inf];
ctype = char(m+m+1);
ctype(1:m) = 'B';
ctype(m+1:m+m) = 'C';
ctype(m+m+1) = 'C';
options = cplexoptimset('MaxTime',600,'TolXInteger',1e-7,'TolRLPFun',1e-6,'TolFun',1e-6,'Display','off');

while toc(tstart) < MaxTime && gap2 > 1e-5 
    k = k+1;
    Aineq = [Aineq;-y_opt'*D,-1];
    if Utype == 'inf'
        bineq = [bineq;-norm(Q'*y_opt,1)-d0'*y_opt];
    elseif Utype == 'two'
        bineq = [bineq;-norm(Q'*y_opt,2)-d0'*y_opt];
    end 
    [z_eta,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],options);
    LB = fval;
    %output.cplexstatus
    z = z_eta(1:end-1);
    if Utype == 'inf'
        [y_opt,v_opt,time,gap,iter] = wclo_q_one([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z,1e-5,1);
    elseif Utype == 'two'
        [y_opt,v_opt,scobound,time,gap,iter,nslo,iscot] = wcsr_sca_cqr([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z,1e-5,1,'cqr');
        v_opt = -v_opt;
    end
    UB = min([UB,[xi;zeta]'*z+v_opt]);
    %[~,y_opt,~,~,~,v_opt,time,status] = wcsr_q_one_milp([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z)
    %[~,y_opt,~,~,~,v_opt,time,status] = wcsr_q_one_milp(A,c,Q,d0-D*z)
    %v_opt = -v_opt
    gap2 = UB-LB;
    fprintf(1,'Iter=%d,LB=%f,UB=%f,gap=%f\n',k,LB,UB,UB-LB);

end
opt_val = fval;
iter = k;
time = toc(tstart);


