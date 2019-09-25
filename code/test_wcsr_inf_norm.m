function test_wcsr_inf_norm
%% Test the performance of Algorithm 4 and MIPR for WCSR with p=\infty
%% and n<=150
NN = 5; %  number of tested problems for each [n,m,l]
R = [];
 for i = 1:10
    for j = 1:NN
        dataname = strcat('.\data_wcsr1015_2norm_60ex\wcsr',num2str((i-1)*NN+j));
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
        [z_opt1,fval1,exitflag1] = cplexlp(c',A,Q*u1+b0,[],[],zeros(n,1),[])
        x_opt1 = 1-z_opt1;
        [z,y,u,v,w,fval2,time2,status] = wcsr_q_one_milp([A,-ones(m,1)],[c;bigM],Q,b0)
        status = double(status);
        fval2 = -fval2;
        R = [R;(i-1)*NN+j,n,m,l,v_opt1,time1,gap1,iter1,mean(x_opt1),fval2,time2,status];
        save wcsr_inf_norm_results.mat R;
    end
end