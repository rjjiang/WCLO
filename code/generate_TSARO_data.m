function generate_TSARO_data
%% Generate 65 test problems for TSRLTP.
mkdir('.\data_TSARO')
T = [3 2; 5 3; 7 4; 10 5; 15 5; 20 5; 30 5; 50 5; 100 5; 20 10; 30 10; 50 10; 100 10];
NN = 5;

for i = 1:13
    for j = 1:NN
        % generate data (A,c,Q,d0)
        m = T(i,1);
        n = T(i,1);
        r = T(i,2);
        xi = 100+900*rand(m,1);
        zeta = 10+90*rand(m,1);
        C = 1+999*rand(m,n);
        c = reshape(C,m*n,1);
        eta = 200+500*rand(m,1);
        Q = -100+200*rand(n,r);    % n*r in this problem
        d0 = 10+490*rand(n,1); % n-dim in this problem
        while sum(eta)-sum(d0)<0 % Make sure that the transportation problem is feasible.
            d0 = 10+490*rand(n,1);
        end

        A = zeros(m+n,m*n);
        for k = 1:m
            A(k,n*(k-1)+1:n*k) = ones(1,n);
            A(m+1:m+n,n*(k-1)+1:n*k) = -eye(n);
        end
        D = [zeros(m),-eye(m);zeros(n,m),zeros(n,m)];
        Q = [zeros(m,r);-Q];
        d0 = [zeros(m,1);-d0];
        bigM = (n+m)*1000; % According to the constraints A'*y<=c,y<=0
        
        % Do scaling
        scale = 1e5;
        xi = xi/scale;
        zeta = zeta/scale;
        c = c/scale;
        bigM = bigM/scale;
        dataname = strcat('.\data_TSARO\tsaro',num2str((i-1)*NN+j));
        save(dataname,'xi','zeta','c','eta','A','D','Q','d0','bigM','scale');
    end
end