function test_wcsr_inf_norm_large
%% Test the performance of Algorithm 4 for WCSR with p=\infty and n>=200
NN = 5; %  number of tested problems for each [n,m,l]
R = [];
 for i = 1:8
    for j = 1:NN
        dataname = strcat('.\data_wcsr1028_allnorm_40ex\wcsr',num2str((i-1)*NN+j));
        load(dataname, 'b0', 'c','Q','A');
        [m,n] = size(A);
        [~,l] = size(Q);
        fprintf(1,'n=%d,m=%d,r=%d\n',n,m,l);
        bigM = 10;
        err = 1e-5;
        [y_opt1,v_opt1,time1,gap1,iter1] = wclo_q_one([A,-ones(m,1)],[c;bigM],Q,b0,err,1);
        u1 = ones(l,1);
        for k = 1:l
            if Q(:,k)'*y_opt1 < 0       u1(k) = -1; end
        end
        [z_opt1,fval1,exitflag1] = cplexlp(c',A,Q*u1+b0,[],[],zeros(n,1),[]);
        x_opt1 = 1-z_opt1;
        R = [R;(i-1)*NN+j,n,m,l,v_opt1,time1,gap1,iter1,mean(x_opt1)];
        save wcsr_inf_norm_large_results.mat R;
    end
end