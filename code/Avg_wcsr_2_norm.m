function Avg_wcsr_2_norm
load('wcsr_2_norm_results.mat')
[row,column] = size(R);
R = [R,zeros(row,2)];
T = R;
for i = 1:row       % BARON succeeded
    if T(i,7) < 1e-2 && T(i,6) < 600
        T(i,column+1) = 1;
    else            % BARON failed
        T(i,5) = 0;
        T(i,6) = 0;
        T(i,7) = 1;
    end
    
    if T(i,25) < 103  % MIPR succeeded
        T(i,column+2) = 1;
    else             % MIPR failed
        T(i,22) = 0;
        T(i,23) = 0;
        T(i,24) = 0;
    end
end

dir = pwd;
txtname = strcat(dir,'/table2.txt');
fid = fopen(txtname,'a');
fprintf(fid,'-------------- Comparison of the average performance of SCOBB, BARON, MIPR and NLSDP for WCSR with p=2 --------------\r\n\r\n');
fprintf(fid,'\\hline                      SCOBB                      &      BARON          &             MIPR   &              NLSDP            \\\\\r\n');
fprintf(fid,'\\hline   n  &  r & Time  & Opt.val &  Iter   & Val_SCO &   Time &  Opt.val   &  Time &  Opt.val   &  Time &    LB   &    UB   & T \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline %3d & %3d & %5.1f & %7.4f & %6.1f  & %7.4f & %6.1f & %7.4f(%d) & %5.1f & %7.4f(%d) & %5.1f & %7.4f & %7.4f & %d \\\\\r\n',T(k,2),T(k,4),...
        mean(T(k:k+4,11)),mean(T(k:k+4,9)),mean(T(k:k+4,13)),mean(T(k:k+4,10)),...
        sum(T(k:k+4,6))/sum(T(k:k+4,column+1)),sum(T(k:k+4,5))/sum(T(k:k+4,column+1)),sum(T(k:k+4,column+1)),...
        sum(T(k:k+4,23))/sum(T(k:k+4,column+2)),sum(T(k:k+4,22))/sum(T(k:k+4,column+2)),sum(T(k:k+4,column+2)),...
        mean(T(k:k+4,19)),mean(T(k:k+4,16)),mean(T(k:k+4,17)),sum(T(k:k+4,21)));
end
fclose(fid);