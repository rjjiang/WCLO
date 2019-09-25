function generate_wcsr_data
%% Generate 60 test problems for WCSR with n = 20,30,50,80,100,150 and r = 5,10 
mkdir('.\data_wcsr1015_2norm_60ex')
NN = 5; % total number of tested problem
T = [20 5; 30 5; 50 5; 80 5; 100 5; 20 10; 30 10; 50 10; 80 10;100 10;150 5;150 10];
for i = 1:12
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
        dataname = strcat('.\data_wcsr1015_2norm_60ex\wcsr',num2str((i-1)*NN+j));
        save(dataname,'A','c','Q','b0');
    end
end