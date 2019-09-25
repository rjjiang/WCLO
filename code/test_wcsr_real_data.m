function test_wcsr_real_data(period)
%% Test the performance of SCOBB, Alg. 4 and Alg. A for WCSR with partial
%% real data.
% period: 1   ----   Dec. 2009
%         2   ----   Jun. 2010

err = 1e-5;
R = [];
for i = 1:5
    if period == 1
        dataname = strcat('.\real_data\wcsr_Dec_2009_',num2str(i));
    elseif period == 2
        dataname = strcat('.\real_data\wcsr_Jun_2010_',num2str(i));
    end
    load(dataname, 'b0', 'c','Q','A','bigM');

    [~,l] = size(Q);
    [m,n] = size(A);

% 2-norm solution
    [y_opt2,v_opt2,scobound2,time2,gap2,iter2,nslo2,iscot2] = wcsr_sca_cqr([A,-ones(m,1)],[c;bigM],Q,b0,err,1,'cqr');
    v_opt2 = -v_opt2;
    u2 = Q'*y_opt2;
    u2 = u2/norm(u2,2);
    [z_opt2,fval2,exitflag2] = cplexlp(c',A,Q*u2+b0,[],[],zeros(n,1),[]);
    x_opt2 = 1-z_opt2;

% inf-norm solution
    [y_opt3,v_opt3,time3,gap3,iter3] = wclo_q_one([A,-ones(m,1)],[c;bigM],Q,b0,err,1);
    u3 = ones(l,1);
    for k = 1:l
        if Q(:,k)'*y_opt3 < 0       u3(k) = -1; end
    end
    [z_opt3,fval3,exitflag3,output3] = cplexlp(c',A,Q*u3+b0,[],[],zeros(n,1),[]);
    x_opt3 = 1-z_opt3;

% 1-norm solution
    [y_opt1,v_opt1,time1,maxindex,sgn] = wclo_q_inf([A,-ones(m,1)],[c;bigM],Q,b0,err,1);
    [z_opt1,fval1,exitflag1] = cplexlp(c',A,sgn(maxindex)*Q(:,maxindex)+b0,[],[],zeros(n,1),[]);
    u1 = zeros(l,1);
    u1(maxindex) = sgn(maxindex);
    x_opt1 = 1-z_opt1;

% fprintf(1,'        u1           u2          uinf\r\n')
% fprintf(1,'EXP: %f,   %f,   %f\r\n',mean(x_opt1),mean(x_opt2),mean(x_opt3))
% fprintf(1,'FVL: %f,   %f,   %f\r\n',v_opt1,v_opt2,v_opt3)
% fprintf(1,'TIM: %f,   %f,   %f\r\n',time1,time2,time3)
    R = [R;i,v_opt1,time1,mean(x_opt1),v_opt2,time2,mean(x_opt2),v_opt3,time3,mean(x_opt3)];
    if period == 1
        save wcsr_real_data_Dec_2009_results.mat R;
    elseif period == 2
        save wcsr_real_data_Jun_2010_results.mat R;
    end
end
