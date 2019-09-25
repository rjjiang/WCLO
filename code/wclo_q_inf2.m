function [v_opt,time] = wclo_q_inf2(A,c,Q,b0)
%% Solve the global solution of WCLO problem:
%             max_{u} min_{x} c'*x 
%        s.t. A*x <= Q*u + b0, x >= 0, ||u||_p <= 1
%  by Integrated SCO-B&B Algorithm
%% that is equivalent to the dual problem:
%             min_{y} f(y):= -||Q^T*y||_q-b0'*y
%        s.t. A^T*y <= c, y<= 0;
% By solving an linear programming
tstart = tic;
[m,n] = size(A);
[~,r] = size(Q);
nn = 2*r*n+1;
mm = 2*r+2*r*m;
AA = zeros(mm,nn);
bb = zeros(mm,1);
EE = [eye(r),-eye(r)];
for i = 1:2*r
    AA(i,(i-1)*n+1:i*n) = c';
    AA(i,nn) = -1;
    bb(i) = 0;
    AA(2*r+(i-1)*m+1:2*r+i*m,(i-1)*n+1:i*n) = A;
    AA(2*r+(i-1)*m+1:2*r+i*m,nn) = zeros(m,1);
    bb(2*r+(i-1)*m+1:2*r+i*m) = Q*EE(:,i)+b0;
end
lb = zeros(nn,1);
lb(nn) = -inf;
f = zeros(nn,1);
f(nn) = 1;
[x,fval,exitflag] = cplexlp(f,AA,bb,[],[],lb,[]);
v_opt = fval;
time = toc(tstart);