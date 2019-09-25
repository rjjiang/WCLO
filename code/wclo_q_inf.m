function [x_opt,v_opt,time,maxindex,sgn] = wclo_q_inf(A,c,Q,b0,err,kk)
%% Solve the global solution of WCLO problem:
%             max_{u} min_{x} c'*x 
%        s.t. A*x <= Q*u + b0, x >= 0, ||u||_p <= 1
%  by Integrated SCO-B&B Algorithm
%% that is equivalent to the dual problem:
%             min_{y} f(y):= -||Q^T*y||_q-b0'*y
%        s.t. A^T*y <= c, y<= 0;
% By solving 2*r linear programming
tstart = tic;
[nm,nn] = size(A);
[~,nr] = size(Q);
i = 1;
Yopt = zeros(nm+1,nr);
sgn = ones(nr,1);
while i <= nr
    [yu1,fval1,exitflag] = cplexlp(-[b0',1],[A',zeros(nn,1)],c,[Q(:,i)',1],0,[-inf*ones(nm,1);0],[zeros(nm,1);inf]);
    [yu2,fval2,exitflag] = cplexlp(-[b0',1],[A',zeros(nn,1)],c,[Q(:,i)',-1],0,[-inf*ones(nm,1);0],[zeros(nm,1);inf]);
    if fval1<fval2
        Yopt(1:nm,i) = yu1(1:nm);
        Yopt(nm+1,i) = -fval1;
        sgn(i) = -1;
    else
        Yopt(1:nm,i) = yu2(1:nm);
        Yopt(nm+1,i) = -fval2;
    end
    i = i+1;
end
[v_opt,maxindex] = max(Yopt(nm+1,:));
x_opt = Yopt(1:nm,maxindex);
time = toc(tstart);

    