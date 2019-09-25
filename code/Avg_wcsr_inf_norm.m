function Avg_wcsr_inf_norm
load('wcsr_inf_norm_results.mat')
[row,column] = size(R);
R = [R,zeros(row,1)];
T = R;
for i = 1:row
    if T(i,12) < 103
        T(i,column+1) = 1;
        R(i,column+1) = 1;
    else
        T(i,10) = 0;
        T(i,11) = 0;
    end
end

dir = pwd;
txtname = strcat(dir,'/table5.txt');
fid = fopen(txtname,'a');
fprintf(fid,'-------------- Table 5 Average numerical results of Algorithm 4 and MIPR for WCSR with p=infty    ---------------\r\n\r\n');
fprintf(fid,'\\hline      Size   &            Alg4                      &             MILP          \\\\\r\n');
fprintf(fid,'\\hline   n   &   r &  time  &     vopt   & Iter &   avg.x &   time     &      vopt    \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline%4d  & %3d & %7.4f & %9.6f  &%5.1f & %7.4f & %7.4f(%d) & %9.6f(%d) \\\\\r\n',T(k,2),T(k,4),...
        mean(T(k:k+4,6)),mean(T(k:k+4,5)),mean(T(k:k+4,8)),mean(T(k:k+4,9)),sum(T(k:k+4,11))/sum(T(k:k+4,column+1)),sum(T(k:k+4,column+1)),...
        sum(T(k:k+4,10))/sum(T(k:k+4,column+1)),sum(T(k:k+4,column+1)));
end
fclose(fid);