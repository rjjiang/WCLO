function [x_opt,v_opt,scobound,time,gap,iter,nslo,iscot] = wcsr_sca_cqr(A,c,Q,b0,err,kk,relax,MaxTime)
%% Find a globally optimal solution of WCLO problem:
%             max_{u} min_{x} c'*x 
%        s.t. A*x <= Q*u + b0, x >= 0, ||u||_2 <= 1
%  by Integrated SCO-B&B Algorithm.
%% It is equivalent to the dual problem:
%             min_{y} f(y):= -||Q^T*y||_2-b0'*y
%        s.t. A^T*y <= c, y <= 0;
% kk:       1 ---- restart SCA
%           0 ---- do not restart SCA
% relax:  cqr ---- Convex Quadratic Relaxation (15)
%        cqcr ---- Conditionally Quasi-Convex Relaxation (a SDP Relaxation)

if nargin < 8
    MaxTime = 600;
end
print = 0;
[nm,nn] = size(A);
[~,nr] = size(Q);
zl0 = zeros(nr,1);
zu0 = zeros(nr,1);
for i=1:nr
    [~,zl0(i),exitflag] = cplexlp(Q(:,i)',A',c,[],[],[],zeros(nm,1));
    [~,zu0(i),exitflag] = cplexlp(-Q(:,i)',A',c,[],[],[],zeros(nm,1));
end
zu0 = -zu0;
%
if nr == 2
    [x1,x2] = ndgrid([-1,1],[-1 1]);
    rho = [vec(x1),vec(x2)]';
elseif nr == 3
    [x1,x2,x3] = ndgrid([-1,1],[-1 1],[-1 1]);
    rho = [vec(x1),vec(x2),vec(x3)]';
elseif nr == 4
    [x1,x2,x3,x4] = ndgrid([-1,1],[-1 1],[-1 1],[-1 1]);
    rho = [vec(x1),vec(x2),vec(x3),vec(x4)]';
elseif nr == 5
    [x1,x2,x3,x4,x5] = ndgrid([-1,1],[-1 1],[-1 1],[-1 1],[-1 1]);
    rho = [vec(x1),vec(x2),vec(x3),vec(x4),vec(x5)]';
else
    %rho = [-ones(nr,1),ones(nr,1)];
    rho = [eye(nr),-eye(nr),-ones(nr,1),ones(nr,1)];
end

MaxIter = 100000;
epsilon = 1e-5;
if print == 1
    fprintf(1,'========= Solving by Branch and Bound =========\n');
end
tstart = tic;

nrho = size(rho,2);
x_all = zeros(nm,nrho);
v_all = zeros(nrho,1);
for ii = 1:nrho
%     f = [zeros(nm,1);rho(:,ii)];
%     Aeq = [Q',-eye(nr)];
%     beq = zeros(nr,1);
%     Aineq = [A',zeros(nn,nr)];
%     bineq = c;
%     xlb = [-inf*ones(nm,1);zl0];
%     xub = [zeros(nm,1);zu0];
%     [x0,~,exitflag] = cplexlp(f,Aineq,bineq,Aeq,beq,xlb,xub);
%     [x1,t] = wclo_sca(A,c,Q,b0,err,x0(1:nm));
    [x0,~,exitflag] = cplexlp(Q*rho(:,ii),[A';-Q';Q'],[c;-zl0;zu0],[],[],[],zeros(nm,1));
    [x1,t] = wclo_sca(A,c,Q,b0,err,x0);
    v1 = -t; % or v1 = -norm(Q'*x1)-b0'*x1; 
    v_all(ii) = v1;
    x_all(:,ii) = x1;
end
[v_opt,imax] = min(v_all); %% best upper bound 
v0 = v_opt;
scobound = v_opt;
%fprintf(1,'scobound = %f\n',scobound);
x_opt = x_all(:,imax);
iscot = toc(tstart);
nslo = 0;
clear x_all v_all;

% solve the relaxation at root nodes
if strcmp(relax,'cqr')
    [t,x,s,z,v] = wclo_cqr(A,c,Q,b0,zl0,zu0);
end
if strcmp(relax,'cqcr')
    [x,Y,v,t] = wclo_cqcr(A,c,Q,b0,zl0,zu0,epsilon);
    s = zeros(nr,1);
    for i = 1:nr
        s(i) = trace(Q(:,i)*Q(:,i)'*Y);
    end
    z = Q'*x;
end
prob.zlb = zl0;
prob.zub = zu0;
prob.dist = s-z.^2;
prob.lb = -t;
prob.s = s;
prob.z = z;
v = -v;
if kk == 1 && v < v_opt-0.001*abs(v_opt)
    nslo = nslo+1;    
    [x1,t] = wclo_sca(A,c,Q,b0,err,x);    
    v1 = -t; %or v1 = -norm(Q'*x1)-b0'*x1;
    fprintf(2,'v1 - v = %.10f\n',v1-v);
    if v1 < v_opt-1.0e-4
        fprintf(1,'Call SCA successfully!\n');
        v_opt = v1;
        x_opt = x1;
    elseif abs(v1-v_opt) <= 1.0e-4
        x_opt = [x_opt,x1];
    end
else
    if v < v_opt-1.0e-4
        v_opt = v;
        x_opt = x;
    elseif abs(v-v_opt) <= 1.0e-4
        x_opt = [x_opt,x];
    end
end

%%
LB = prob.lb;
AllNodes = [];
iter = 0;
if print == 1
    fprintf(1,'Iter = %3d,  LB = %12.8f,  UB = %12.8f, gap = %5.2f%%\n',...
        iter,LB,v_opt,100*abs(v_opt-LB)/(1.0e-10+abs(v_opt)));
end
subprob = cell(2,1);
while iter < MaxIter && toc(tstart) < MaxTime  && abs(v_opt-LB)/(1.0e-10+abs(v_opt))>5e-5 && abs(v_opt-LB)>5e-5
    iter = iter+1;
    subprob{1} = prob;
    subprob{2} = prob;    
    [~,i_max] = max(prob.dist);
%     fprintf(1,'dist = %d,  %.8f\n',i_max,max(prob.dist));
    midpoint = (prob.zub(i_max)+prob.zlb(i_max))/2;
    if prob.s(i_max) <= (prob.zlb(i_max)+midpoint)*prob.z(i_max)-prob.zlb(i_max)*midpoint...
            || prob.s(i_max) <= (prob.zub(i_max)+midpoint)*prob.z(i_max)-prob.zub(i_max)*midpoint        
         midpoint = prob.z(i_max);
    end
    subprob{1}.zub(i_max) = midpoint;
    subprob{2}.zlb(i_max) = midpoint;
    
    for ii=1:2
        if strcmp(relax,'cqr')
            [t,x,s,z,v] = wclo_cqr(A,c,Q,b0,subprob{ii}.zlb,subprob{ii}.zub);
        end
        if strcmp(relax,'cqcr')
            [x,Y,v,t,bsa_iter] = wclo_cqcr(A,c,Q,b0,subprob{ii}.zlb,subprob{ii}.zub,epsilon);
            s = zeros(nr,1);
            for i = 1:nr
                s(i) = trace(Q(:,i)*Q(:,i)'*Y);
            end
            z = Q'*x;
        end
        if ~isempty(x)    
            subprob{ii}.dist = s-z.^2;
            subprob{ii}.lb = -t;
            subprob{ii}.s = s;
            subprob{ii}.z = z;
            v = -v;           
            if kk == 1 && v < v_opt-0.001*abs(v_opt)
                nslo = nslo+1;
                [x1,t] = wclo_sca(A,c,Q,b0,err,x);
                v1 = -t;               
                fprintf(2,'v1 - v = %.10f\n',v1-v); 
                if v1 < v_opt-1.0e-4
                    fprintf(1,'Call SCA successfully!\n');
                    v_opt = v1;
                    x_opt = x1;
                elseif abs(v1-v_opt) <= 1.0e-4
                    x_opt = [x_opt,x1];
                end
            else
                if v < v_opt-1.0e-4
                    v_opt = v;
                    x_opt = x;         
                elseif abs(v-v_opt) <= 1.0e-4
                    x_opt = [x_opt,x];
                end
            end
            if sum(subprob{ii}.dist) > err
                AllNodes = [AllNodes,subprob{ii}];
            end
        end
    end  
%     fprintf(1,'lb1=%.8f, lb2=%.8f\n',subprob(1).lb,subprob(2).lb);
    % cut branch
    if ~isempty(AllNodes)
        ii = 0;
        while ii < length(AllNodes)
            ii = ii+1;
            if AllNodes(ii).lb >= v_opt-err
                AllNodes(ii) = [];
                ii = ii-1;
            end
        end
    end
    if isempty(AllNodes) 
        LB = v_opt;
        break;
    else
        LB = AllNodes(1).lb;
        index_LB = 1;
        ii = 1;
        while ii < length(AllNodes)
            ii = ii+1;
            if AllNodes(ii).lb < LB
                LB = AllNodes(ii).lb;
                index_LB = ii;
            end
        end
        prob = AllNodes(index_LB);
%         fprintf(1,'norm = %.8f, err = %.8f\n',norm(prob.tub-prob.tlb)^2/4,err);
        AllNodes(index_LB) = [];
    end
    if mod(iter,100) == 0 & print == 1
        fprintf(1,'Iter = %3d,  LB = %12.8f,  UB = %12.8f, gap = %5.2f%%\n',...
            iter,LB,v_opt,100*abs(v_opt-LB)/(1.0e-10+abs(v_opt)));
    end
end
numsol = size(x_opt,2);
x_opt = x_opt(:,1);
if print == 1
    fprintf(1,'Iter = %3d,  LB = %12.8f,  UB = %12.8f, gap = %5.2f%%\n',...
        iter,LB,v_opt,100*abs(v_opt-LB)/(1.0e-10+abs(v_opt)));
end
gap = abs(v_opt-LB);
time = toc(tstart);
end
