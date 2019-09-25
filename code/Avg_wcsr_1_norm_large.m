function Avg_wcsr_1_norm_large
load('wcsr_1_norm_large_results.mat')
[row,column] = size(R);
dir = pwd;
txtname = strcat(dir,'/table7.txt');
fid = fopen(txtname,'a');
fprintf(fid,'-------------- Average numerical results of Algorithm A and B for WCSR with p=1    ---------------\r\n\r\n');
fprintf(fid,'\\hline     Size   &                 Alg. A       &            Alg. B              \\\\\r\n');
fprintf(fid,'\\hline   n  &   r &   Time  &   Opt.val  & avg.x &  Time   &   Opt.val  & avg.x   \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline%4d  & %3d & %7.4f & %9.6f  &%6.4f & %7.4f & %9.6f  &%6.4f   \\\\\r\n',R(k,2),R(k,4),...
        mean(R(k:k+4,6)),mean(R(k:k+4,5)),mean(R(k:k+4,7)),mean(R(k:k+4,9)),mean(R(k:k+4,8)),mean(R(k:k+4,7)));
end
fclose(fid);