function test_TSARO_1_norm
%% Test the performance of Alg. C for TSRLTP with p=1
NN = 5; % total number of tested problem
R = [];
for i = 6:13
    for j = 1:NN
        dataname = strcat('.\data_TSARO\tsaro',num2str((i-1)*NN+j));
        load(dataname,'xi','zeta','c','eta','A','D','Q','d0','bigM','scale');
        m = length(xi);
        [n,r] = size(Q);
        n = n-m;
        fprintf(1,'n=%d,m=%d,r=%d\n',n,m,r);
        [v_opt,time,status] = TSARO_1norm(xi,zeta,c,eta,A,D,Q,d0,bigM);
        status = double(status);
        R = [R;(i-1)*NN+j,m,n,r,v_opt,time,status];
        save tsaro_results_1_norm.mat R;
    end
end