function generate_wcsr_real_data
%% Generate 10 test problems for WCSR with partial real data.
mkdir('.\real_data')
NN = 5;
N = 8;
l = 5;
%% Dec. 2009
L = [0.00 500.62 341.62 409.36 189.95 231.97 36.22 10.43;
     172.97 0.00 292.94 51.02 176.58 36.35 20.52 4.62;
     239.17 195.64 0.00 50.42 92.73 20.60 32.57 8.08;
     114.14 237.98 219.64 0.00 119.73 30.23 26.56 28.08;
     96.69 155.65 150.57 22.82 0.00 15.47 28.11 11.39;
     187.51 183.76 60.33 15.66 30.82 0.00 64.50 21.52;
     30.72 40.68 301.37 9.42 131.55 6.11 0.00 1.17;
     24.26 47.38 44.74 86.08 12.41 5.43 3.14 0.00];

for i = 1:NN
    % Initial capital: sample from uniform in [Cmin,Cmax]
    Cmin = 50;
    Cmax = 200;
    C = Cmin+(Cmax-Cmin).*rand(N,1);
    % construct data
    b0 = C - L*ones(N,1) + L'*ones(N,1);
    c = ones(N,1);
    Q = [30*(rand(N,l)-0.5)];
    A = L' - diag(L*ones(N,1));
    bigM = 10;
    dataname = strcat('.\real_data\wcsr_Dec_2009_',num2str(i));
    save(dataname,'A','c','Q','b0','L','C','bigM');
end

%% Jun. 2010
L = [0.00 462.07 327.72 386.37 135.37 208.97 43.14 7.72;
     172.18 0.00 255.00 39.08 149.82 32.11 20.93 3.93;
     257.11 196.84 0.00 26.26 80.84 18.11 29.70 8.21;
     110.85 181.65 162.44 0.00 72.67 25.34 18.75 23.09;
     141.39 148.62 126.38 20.66 0.00 12.45 23.14 11.11;
     148.51 138.57 50.08 13.98 21.20 0.00 53.99 19.38;
     29.15 35.14 253.13 5.67 108.68 5.32 0.00 0.39;
     22.39 37.24 41.90 78.29 5.13 5.15 2.57 0.00];

for i = 1:NN
    % Initial capital: sample from uniform in [Cmin,Cmax]
    Cmin = 50;
    Cmax = 200;
    C = Cmin+(Cmax-Cmin).*rand(N,1);
    % construct data
    b0 = C - L*ones(N,1) + L'*ones(N,1);
    c = ones(N,1);
    Q = [30*(rand(N,l)-0.5)];
    A = L' - diag(L*ones(N,1));
    bigM = 10;
    dataname = strcat('.\real_data\wcsr_Jun_2010_',num2str(i));
    save(dataname,'A','c','Q','b0','L','C','bigM');
end
 
 