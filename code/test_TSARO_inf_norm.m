function test_TSARO_inf_norm
%% Test the performance of Alg. 5, Alg. 6 and C&CG for TSRLTP with p=\infty
NN = 5; % total number of tested problem
R = [];
for i = 1:13
    for j = 1:NN
        dataname = strcat('.\data_TSARO\tsaro',num2str((i-1)*NN+j));
        load(dataname,'xi','zeta','c','eta','A','D','Q','d0','bigM','scale');
        m = length(xi);
        [n,r] = size(Q);
        n = n-m;
        fprintf(1,'n=%d,m=%d,r=%d\n',n,m,r);
        % Alg. 5
        [LB1,UB1,time1,iter1] = TSARO_CPA(xi,zeta,c,eta,A,D,Q,d0,bigM,'inf');
        
        % Alg. 6
        [LB2,UB2,time2,iter2] = TSARO_CCG(xi,zeta,c,eta,A,D,Q,d0,bigM,'inf','FBBA');
        
        % C&CG
        [LB3,UB3,time3,iter3] = TSARO_CCG(xi,zeta,c,eta,A,D,Q,d0,bigM,'inf','MILP');
        
        R = [R;(i-1)*NN+j,m,n,r,LB1,UB1,time1,iter1,LB2,UB2,time2,iter2,LB3,UB3,time3,iter3];
        save tsaro_results_inf_norm.mat R;
    end
end