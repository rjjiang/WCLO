function test_wcsr_2_norm
%% Test the performance of SCOBB, BARON, MIPR and NLSDP for WCSR with p=2 and
%% n<=100 to get table 2 in paper

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
        % SCOBB
        [x_opt1,v_opt1,scobound1,time1,gap1,iter1,nslo1,iscot1] = wcsr_sca_cqr([A,-ones(m,1)],[c;bigM],Q,b0,err,1,'cqr');
        v_opt1 = -v_opt1;
        scobound1 = -scobound1;
        
        % BARON
        fun2 = @(y) (-norm(Q'*y,2)-b0'*y);
        % For win64 system: (The path below should be changed according to your computer system)
        opts = baronset('CplexLibName','C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab\x64_win64\cplex1260.dll','MaxTime',600);
        % For win32 system: (The path below should be changed according to your computer system)
        % opts = baronset('CplexLibName','C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\bin\x86_win32\cplex1260.dll','MaxTime',600);
        tstart = tic;
        [x2,fval2,exitflag2,info] = baron(fun2,[A';-ones(1,m)],-inf*ones(n+1,1),[c;bigM],-inf*ones(m,1),zeros(m,1),[],[],[],[],[],opts);
        fval2 = -fval2;
        time2 = toc(tstart);
        iter2 = info.BaR_Iterations;
        gap2 = info.Upper_Bound-info.Lower_Bound;
        
        % MIPR
        [z3,y3,u3,v3,w3,fval3,time3,status3] = wcsr_q_2_miqcp([A,-ones(m,1)],[c;bigM],Q,b0);
        x_opt3 = 1 - z3;
        status3 = double(status3);
        
        % NLSDP
        [lb4,ub4,bisiter4,time4,status4,flag4] = nlsdp([A,-ones(m,1)],Q,b0,[c;bigM]);
        
        R = [R;(i-1)*NN+j,n,m,l,fval2,time2,gap2,iter2,v_opt1,scobound1,time1,gap1,iter1,nslo1,iscot1,...
            lb4,ub4,bisiter4,time4,status4,flag4,-fval3,time3,mean(x_opt3),status3];
        save wcsr_2_norm_results.mat R;
    end
 end

 
 
 
 
 
 
 