function Avg_wcsr_2_norm_large
load('wcsr_2_norm_large_results.mat')
[row,column] = size(R);
dir = pwd;
txtname = strcat(dir,'/table4.txt');
fid = fopen(txtname,'a');
fprintf(fid,'-------------- Average numerical results of SCOBB for WCSR with p=2 --------------\r\n\r\n');
fprintf(fid,'\\hline      Size               &    SCOBB                          \\\\\r\n');
fprintf(fid,'\\hline   n  &  r & Time  &  Opt.val  &  Iter   &  Val_SCO  &   avg.x  \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline %3d & %3d & %5.1f & %9.6f & %6.1f  & %9.6f & %6.4f\\\\\r\n',...
    R(k,2),R(k,4),mean(R(k:k+4,7)),mean(R(k:k+4,5)),mean(R(k:k+4,9)),mean(R(k:k+4,6)),mean(R(k:k+4,12)));
end
fclose(fid);