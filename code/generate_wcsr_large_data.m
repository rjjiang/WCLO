function generate_wcsr_large_data
%% Generate 40 test problems for WCSR with n = 200,300,400,500 and r = 5,10 
mkdir('.\data_wcsr1028_allnorm_40ex')
NN = 5; % total number of tested problem
T = [200 5; 300 5; 400 5; 500 5; 200 10; 300 10; 400 10; 500 10];
for i = 1:8
    for j = 1:NN
        N = T(i,1);
        l = T(i,2);
        % generate the liability matrix L (N by N)
        Cmin = 0.5; % initial capital lower bound
        Cmax = 5; % initial capital upper bound
        muL = 0.03; % mean of liabilities of all banks
        sigL = 1; % sigma of liabilities of all banks

        % Initial capital: sample from uniform in [Cmin,Cmax]
        C = Cmin+(Cmax-Cmin).*rand(N,1);
        % Liability matrix L: sample from LOGN(muL,sigL)
        L = lognrnd(muL,sigL,[N,N]);
        for k = 1:N
            L(k,k)=0;
        end
        % construct data
        b0 = C - L*ones(N,1) + L'*ones(N,1);
        c = ones(N,1);
        Q = [6*(rand(N,l)-0.5)];
        A = L' - diag(L*ones(N,1));

        % save data
        dataname = strcat('.\data_wcsr1028_allnorm_40ex\wcsr',num2str((i-1)*NN+j));
        save(dataname,'A','c','Q','b0');
    end
end