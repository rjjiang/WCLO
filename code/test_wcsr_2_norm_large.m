function test_wcsr_2_norm_large
%% Test the performance of SCOBB for WCSR with p=2 and n>100 to get table 4 in paper

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
        % SCOBB
        [y_opt1,v_opt1,scobound1,time1,gap1,iter1,nslo1,iscot1] = wcsr_sca_cqr([A,-ones(m,1)],[c;bigM],Q,b0,err,1,'cqr');
        v_opt1 = -v_opt1;
        scobound1 = -scobound1;
        
        u1 = Q'*y_opt1;
        u1 = u1/norm(u1,2);
        [x_opt1,fval2,exitflag2] = cplexlp(c',A,Q*u1+b0,[],[],zeros(n,1),[]);
        x = 1-x_opt1;
        
        
        R = [R;(i-1)*NN+j,n,m,l,v_opt1,scobound1,time1,gap1,iter1,nslo1,iscot1,mean(x)];
        save wcsr_2_norm_large_results.mat R;
    end
 end