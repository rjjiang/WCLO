function Avg_wcsr_inf_norm_large
load('wcsr_inf_norm_large_results.mat')
[row,column] = size(R);
dir = pwd;
txtname = strcat(dir,'/table6.txt');
fid = fopen(txtname,'a');
fprintf(fid,'-------------- Average numerical results of Algorithm 4 for WCSR with p=infty    ---------------\r\n\r\n');
fprintf(fid,'\\hline     Size   &                  Alg4                    \\\\\r\n');
fprintf(fid,'\\hline   n  &   r &    time  &     vopt  & Iter &   avg.x    \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline%4d  & %3d & %7.4f & %9.6f  &%5.1f & %7.4f  \\\\\r\n',R(k,2),R(k,4),...
        mean(R(k:k+4,6)),mean(R(k:k+4,5)),mean(R(k:k+4,8)),mean(R(k:k+4,9)));
end
fclose(fid);