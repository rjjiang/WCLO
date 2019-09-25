function [LB,UB,time,iter] = TSARO_CCG(xi,zeta,c,eta,A,D,Q,d0,bigM,Utype,SP2solver)
%% The Column-and-constraint generation (C&CG) algorithm framework in Zeng and Zhao
%% (2013) to solve the following problem named TSARO:
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
% SP2solver: MILP  ---- 0-1 mixed integer program subproblem solver for p=\infty in C&CG
%            FBBA  ---- a finite B&B algorithm subproblem solver for p=\infty in algorithm 5
%            MIQCP ---- 0-1 mixed integer quadratic constrained program subproblem solver for p=2 in C&CG
%            SCOBB ---- the SCOBB algorithm subproblem solver for p=2 in algorithm 5
%
% Note: If SP2solver is set to MILP or MIQCP, it is the C&CG algorithm in Zeng
% and Zhao (2013), and if SP2solver is set to be FBBA or SCOBB, it is Algorithm 6 (the ICP
% Algorithm) actually. 

tstart = tic;
MaxTime = 600;
m = length(xi);
[n,r] = size(Q);
n = n-m;
%bigM = (n+m)*1000; % According to the constraints A'*y<=c,y<=0
LB = -inf;
UB = inf;
k = 0;
O = [];
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
%[z_eta,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],options)
z_eta = [ones(m,1);eta;-inf];
fval = -inf;
LB = fval;
if Utype == 'inf'
    if SP2solver == 'MILP'
        %[y,v_opt,~,~,~] = wclo_q_one([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m),1e-5,1)
        [z,y,u,v,w,v_opt,time,status] = TSARO_SP2_infnorm([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m));
        
%     u,
%     [A,-ones(m+n,1)]*z-Q*u-(d0-D*z_eta(1:2*m))
%     A'*y-c
        %[z,y,u,v,w,v_opt,time,status] = TSARO_SP2_infnorm(A,c,Q,d0-D*z_eta(1:2*m));
    elseif SP2solver == 'FBBA'
        [y,v_opt,~,~,~] = wclo_q_one([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m),1e-5,1);
        u = ones(r,1);
        for i = 1:r
            if Q(:,i)'*y < 0       u(i) = -1; end
        end
    
%     u,
%     [x_opt1,fval1,exitflag1] = cplexlp(c',A,Q*u+d0-D*z_eta(1:2*m),[],[],zeros(m*n,1),[]);
%     A*x_opt1-Q*u-(d0-D*z_eta(1:2*m))
%     A'*y-c
    end
elseif Utype == 'two'
    if SP2solver == 'MIQCP'
        [z,y,u,v,w,v_opt,time,status] = TSARO_SP2_2norm([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m));
%         v_opt
%         status
%         [y,v_opt,~,~,~,~,~,~] = wcsr_sca_cqr([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m),1e-5,1,'cqr')
%         v_opt
    elseif SP2solver == 'SCOBB'
        [y,v_opt,~,time,~,iter,~,~] = wcsr_sca_cqr([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m),1e-5,1,'cqr');
        v_opt = -v_opt;
        u = Q'*y;
        u = u/norm(u,2);
    end
end
%[y_opt,v_opt,time,gap,iter] =
%wclo_q_one([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z,1e-5,1);
UB = min([UB,[xi;zeta]'*z_eta(1:2*m)+v_opt]);
gap = UB-LB;
%fprintf(1,'Iter=%d,LB=%f,UB=%f,gap=%f\n',k,LB,UB,gap);
while toc(tstart) < MaxTime && gap > 1e-5
    k = k+1;
    if v_opt < inf
        f = [f;zeros(m*n,1)];
        [row,col] = size(Aineq);
        % add variable x^{k+1}
        Aineq = [Aineq,zeros(row,m*n)];
        % add eta >= b^T*x^{k+1}
        Aineq(row+1,m+m+1) = -1;
        Aineq(row+1,m+m+2+(k-1)*(m*n):m+m+1+k*(m*n)) = c';
        bineq = [bineq;0];
        % add Ey+Gx^{k+1}>=h-Mu_{k+1}^*
        Aineq(row+2:row+1+m+n,1:m+m) = D;
        Aineq(row+2:row+1+m+n,m+m+2+(k-1)*(m*n):m+m+1+k*(m*n)) = A;
        bineq = [bineq;Q*u+d0];
        lb = [lb;zeros(m*n,1)];
        ub = [ub;inf*ones(m*n,1)];
        ctype(col+1:col+m*n) = 'C'; 
        options = cplexoptimset('MaxTime',600,'TolXInteger',1e-7,'TolRLPFun',1e-6,'TolFun',1e-6,'Display','off');
        [z_eta,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,[],[],[],[],[],lb,ub,ctype,[],options);
        %output.cplexstatus
        
        LB = fval;
        if Utype == 'inf'
            if SP2solver == 'MILP'
                %[y,v_opt,~,~,~] = wclo_q_one([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m),1e-5,1)
                [z,y,u,v,w,v_opt,time,status] = TSARO_SP2_infnorm([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m));
%             u,
            %[z,y,u,v,w,v_opt,time,status] = TSARO_SP2_infnorm(A,c,Q,d0-D*z_eta(1:2*m));
            elseif SP2solver == 'FBBA'
                [y,v_opt,~,~,~] = wclo_q_one([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m),1e-5,1);
                u = ones(r,1);
                for i = 1:r
                    if Q(:,i)'*y < 0       u(i) = -1; end
                end
            end
        elseif Utype == 'two'
            if SP2solver == 'MIQCP'
                [z,y,u,v,w,v_opt,time,status] = TSARO_SP2_2norm([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m));
            elseif SP2solver == 'SCOBB'
                [y,v_opt,~,time,~,iter,~,~] = wcsr_sca_cqr([A,-ones(m+n,1)],[c;bigM],Q,d0-D*z_eta(1:2*m),1e-5,1,'cqr');
                v_opt = -v_opt;
                u = Q'*y;
                u = u/norm(u,2);
            end
        end
    end
    UB = min([UB,[xi;zeta]'*z_eta(1:2*m)+v_opt]);
    gap = UB-LB;
    fprintf(1,'Iter=%d,LB=%f,UB=%f,gap=%f\n',k,LB,UB,gap);

end
opt_val = LB;
iter = k;
time = toc(tstart);
