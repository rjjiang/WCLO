function [y_opt,v_opt,time,gap,iter] = wclo_q_one(A,c,Q,b0,err,kk)
%% Find a globally optimal solution of WCLO problem with p = \infty:
%             max_{u} min_{x} c'*x 
%        s.t. A*x <= Q*u + b0, x >= 0, ||u||_p <= 1
%  by Integrated SCO-B&B Algorithm
%% It is equivalent to the dual problem (q = 1):
%             min_{y} f(y):= -||Q^T*y||_q-b0'*y
%        s.t. A^T*y <= c, y <= 0;


tstart = tic;
[nm,nn] = size(A);
[~,nr] = size(Q);
zl0 = zeros(nr,1);
zu0 = zeros(nr,1);
for i=1:nr
    [~,zl0(i),exitflag] = cplexlp(Q(:,i)',A',c,[],[],[],zeros(nm,1));
    %linprog(-Q(:,i)',A',c,[],[],[],zeros(nm,1))
    [~,zu0(i),exitflag] = cplexlp(-Q(:,i)',A',c,[],[],[],zeros(nm,1));
end
zu0 = -zu0;

% fprintf(1,'========= Solving by Branch and Bound =========\n');
%% get initial lower bound
u = zeros(2*nr,1);
% max b0'*y + e'*z
% s.t. A'*y <= c,
%      y <= 0,
%      Q'*y = z.
yz = cplexlp(-[b0',ones(1,nr)],[A',zeros(nn,nr)],c,[Q',-eye(nr)],zeros(nr,1),-inf*ones(nm+nr,1),[zeros(nm,1);inf*ones(nr,1)]);
u(1:nr) = max(yz(nm+1:nm+nr),0);
u(nr+1:2*nr) = -min(yz(nm+1:nm+nr),0);
yu_opt = [yz(1:nm);u];
v_opt = [b0',ones(1,2*nr)]*yu_opt; % Initial lower bound
y_opt = yz(1:nm);

%% solve the relaxation at root nodes
uub1 = max(zu0,0);
uub2 = -min(zl0,0);
uub = [uub1;uub2];
% max b0'*y + e'*u
% s.t. A'*y <= c,
%      y <= 0,
%      Q'*y = (I,-I)*u,
%      0 <= u <= uub.
[yu,fval,exitflag] = cplexlp(-[b0',ones(1,2*nr)],[A',zeros(nn,2*nr)],c,[Q',-eye(nr),eye(nr)],zeros(nr,1),[-inf*ones(nm,1);zeros(2*nr,1)],[zeros(nm,1);uub]);
if yu(nm+1:nm+nr)'*yu(nm+nr+1:nm+2*nr) == 0
    if -fval > v_opt
        v_opt = -fval;
        y_opt = yu(1:nm);
    end
end

prob.ub = -fval;   % Initial upper bound obtained from the root node
prob.yu = yu;
prob.I = zeros(nr,1);
for i = 1:nr
    if yu(nm+i) == 0
        prob.I(i) = 1;
    elseif yu(nm+nr+i) == 0
        prob.I(i) = 2;
    end
end

UB = prob.ub;
AllNodes = [];
iter = 0;
% fprintf(1,'Iter = %3d,  LB = %12.8f,  UB = %12.8f, gap = %5.2f%%\n',...
%     iter,v_opt,UB,100*abs(UB-v_opt)/(1.0e-10+abs(UB)));
gap = 100*abs(UB-v_opt)/(1.0e-10+abs(UB));
subprob = cell(2,1);
status = 0;

while gap > 0
    iter = iter+1;
    subprob{1} = prob;
    subprob{2} = prob;
    %[~,i_branch] = min(prob.I);
    [~,i_branch] = max(prob.yu(nm+1:nm+nr).*prob.yu(nm+nr+1:nm+2*nr));
    subprob{1}.I(i_branch) = 1;
    subprob{2}.I(i_branch) = 2;
    
    for ii = 1:2
        uub1 = max(zu0,0);
        uub2 = -min(zl0,0);
        uub = [uub1;uub2];
        for i = 1:nr % nr can be changed to i_branch
            if subprob{ii}.I(i) == 1
                uub(i) = 0;
            elseif subprob{ii}.I(i) == 2
                uub(nr+i) = 0;
            end
        end
        
        [yu,fval,exitflag] = cplexlp(-[b0',ones(1,2*nr)],[A',zeros(nn,2*nr)],c,[Q',-eye(nr),eye(nr)],zeros(nr,1),[-inf*ones(nm,1);zeros(2*nr,1)],[zeros(nm,1);uub]);
        subprob{ii}.ub = -fval;
        subprob{ii}.yu = yu;
        AllNodes = [AllNodes,subprob{ii}];
        %yu(nm+1:nm+nr)'*yu(nm+nr+1:nm+2*nr),
        if yu(nm+1:nm+nr)'*yu(nm+nr+1:nm+2*nr) == 0
            if -fval > v_opt
                v_opt = -fval;
                y_opt = yu(1:nm);
                %disp('-----------Lower bound updated-----------')
            end
        end
    end
    % cut branch
    if ~isempty(AllNodes)
        ii = 0;
        while ii < length(AllNodes)
            ii = ii+1;
            if AllNodes(ii).ub <= v_opt+err
                AllNodes(ii) = [];
                ii = ii-1;
            end
        end
    end
    if isempty(AllNodes) 
        UB = v_opt;
        break;
    else
        UB = AllNodes(1).ub;
        index_UB = 1;
        ii = 1;
        while ii < length(AllNodes)
            ii = ii+1;
            if AllNodes(ii).ub > UB
                UB = AllNodes(ii).ub;
                index_UB = ii;
            end
        end
        prob = AllNodes(index_UB);
%         fprintf(1,'norm = %.8f, err = %.8f\n',norm(prob.tub-prob.tlb)^2/4,err);
        AllNodes(index_UB) = [];
    end
    gap = 100*abs(UB-v_opt)/(1.0e-10+abs(UB));
%     if mod(iter,1) == 0
%         fprintf(1,'Iter = %3d,  LB = %12.8f,  UB = %12.8f, gap = %5.2f%%\n',iter,v_opt,UB,100*abs(UB-v_opt)/(1.0e-10+abs(UB)));
%     end
end
% fprintf(1,'Iter = %3d,  LB = %12.8f,  UB = %12.8f, gap = %5.2f%%\n',iter,v_opt,UB,100*abs(UB-v_opt)/(1.0e-10+abs(UB)));
gap = 100*abs(UB-v_opt)/(1.0e-10+abs(UB));
time = toc(tstart);
        
    

