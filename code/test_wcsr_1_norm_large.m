function test_wcsr_1_norm_large
%% Test the performance of Alg. A and B for WCSR with p=1 and n>100
%% to get table 7 in the paper

NN = 5; %  number of tested problems for each [n,m,l]
R = [];
 for i = 1:8
    for j = 1:NN
        dataname = strcat('.\data_wcsr1028_allnorm_40ex\wcsr',num2str((i-1)*NN+j));
        load(dataname, 'b0', 'c','Q','A');
        [m,n] = size(A);
        [~,l] = size(Q);
        fprintf(1,'n=%d,m=%d,r=%d\n',n,m,l);
        err = 1e-5;
        [y_opt,v_opt,time,maxindex,sgn] = wclo_q_inf(A,c,Q,b0,err,1);
        [x_opt3,fval3,exitflag3] = cplexlp(c',A,sgn(maxindex)*Q(:,maxindex)+b0,[],[],zeros(n,1),[]);
        [v_opt2,time2] = wclo_q_inf2(A,c,Q,b0);
        x_opt = 1-x_opt3;
        R = [R;(i-1)*NN+j,n,m,l,v_opt,time,mean(x_opt),v_opt2,time2];
        save wcsr_1_norm_large_results.mat R;
    end
end