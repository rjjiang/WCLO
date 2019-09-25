function [lbqp,ub,i,time,status,flag,nfe] = nlsdp(A,Q,b0,c)
% solving problem: 
%             max_{u} min_{x} c'*x 
%                      s.t. A*x <= Q*u + b0 
%                           x >= 0
%                           ||u||_2 <= 1
% status =1:converged; =0:infeasible; =2:not converged
% flag =1:solved no gap; =0: solved with gap
tstart = tic;
N_loops = 30;
err = 1e-5; % accuracy
MaxTime = 600;

% Step 0:
[m,n]=size(A);
[m,l]=size(Q);
%A = [A;-eye(n,n)];
%Q = [Q;zeros(n,l)];
%b0 = [b0;zeros(n,1)];

% Step 1: obtain intial bounds
%[m,n]=size(A);
%[m,l]=size(Q);
i=1;
nfe = 0;
cvx_begin
% from cvxopt import solvers
%cvx_solver sedumi
cvx_quiet true
variable Y(m,m) symmetric
variable y(m)
variable yt(m)
variable t
variable s
maximize t
subject to
trace([(Q*Q'-b0*b0'),b0;b0', -1]*[Y, yt;yt', s]) >= 0;
[1,y',t;y,Y,yt;t,yt',s] == semidefinite((m+2),(m+2));
% diag(A'*Y*A) == diag(c*c');
Y >= 0;
A'*y <= c;
y <= 0;
A'*Y >= c*y';
% additional contraints to help the solver
yt <= 0;
A'*yt <= c*t;
cvx_end
if (cvx_status(1:6) == 'Inaccu'|cvx_status(1:6) == 'Solved')
    ub = t;
    lb = b0'*y + norm(Q'*y);
    lbqp = b0'*y + norm(Q'*y);
    if (abs(ub-lb) <= err)
        status = 1;
        flag = 1;
    else
        status = 2;
        flag = 0;
    end
else
    ub = 0;
    lb = 0;
    lbqp = 0;
    status = 0;
    flag = 0;
end
delta = (lb+ub)/2;


% Step 2: loop
while status == 2 && i <= N_loops && toc(tstart) < MaxTime
    [u,lbqp,stat,flag,cut_val,eig_max] = sdpqp(A,Q,b0,c,delta);
    nfe = nfe+1;
 %   cut_val1=cut_val;
 %   cut_val1=eig_max;
    if (stat == 0)
        status = 0;
    else
        if (u > err)
            lb = max(lbqp,delta);
            if ((flag == 1)&&(abs(u) <= 0.1*delta)) 
                delta = lbqp; % fast update if get rank one solution
            else
                delta = (lb+ub)/2;
            end
        elseif (u < -err)
            ub = delta;
            lb = max(lbqp,lb);
            if ((flag == 1)&&(abs(u) <= 0.1*delta)) 
                delta = lbqp; % fast update if get rank one solution
            else
                delta = (lb+ub)/2;
            end
        else
            status = 1;
            ub = delta;
      %      feas1_val=cut_val;
      %      eig_val_max=eig_max;
        end
    end
    i = i+1;
end
time = toc(tstart);
% while ((status == 2)&&(i <= 100))
%     [u,lbqp,stat,flag] = sdpqp(A,Q,b0,c,delta);
%     if (stat == 0)
%         status = 0;
%     else %status = 1:successfully solved
%         if (u > 0.00001)
%             lb = max(lbqp,delta);
%             delta = (lb+ub)/2;
%         elseif (u < -0.00001)
%             ub = delta;
%             lb = max(lbqp,lb);
%             delta = (lb+ub)/2;
%         else
%             status = 1;
%         end
%     end
%     i = i+1;
% end
end


function [u,lbqp,stat,flg,cut_val,eig_max] = sdpqp(A,Q,b0,c,delta)
% stat = 1:solved; =0:failed
% flg =1:rank1; =0:non-rank1
[m,n]=size(A);
[m,l]=size(Q);
cvx_begin
%cvx_solver sedumi
%cvx_solver_settings( 'maxiters', '10' );  % wrong commands
cvx_quiet true
variable Y(m,m) symmetric
variable y(m)
variable c1
maximize c1
subject to
c1 == trace((Q*Q'-b0*b0')*Y)+2*b0'*y*delta - delta^2;
[1,y';y,Y;] == semidefinite((m+1),(m+1));
%diag(A'*Y*A) == diag(c*c');
Y >= 0;
A'*y <= c;
y <= 0;
% additional contraints to help the solver
A'*Y >= c*y';
cvx_end
eig_max=max(eig(Y-y*y'));
cut_val=b0'*(Y-y*y')*b0;

if (cvx_status(1:6) == 'Inaccu'|cvx_status(1:6) == 'Solved')
    u = c1;
    lbqp = b0'*y + norm(Q'*y);
    stat = 1;
    if (norm(Y-y*y') <= 0.00001)
        flg = 1;
    else
        flg = 0;
    end
else
    u = c1;
    lbqp = 0;
    stat = 0;
    flg = 0;
end
end