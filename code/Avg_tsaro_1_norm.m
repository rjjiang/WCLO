function Avg_tsaro_1_norm
load('tsaro_results_1_norm.mat')
[row,col] = size(R);
dir = pwd;
txtname = strcat(dir,'/table12.txt');
fid = fopen(txtname,'a');
fprintf(fid,'--------- Average numerical results of Alg. C for TSRLTP with p=1 -----------\r\n\r\n');

fprintf(fid,'\\hline   m  &   n  &   r  &  Opt.val  &  Time    \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline %4d & %4d & %4d & %9.6f & %9.6f\\\\\r\n',R(k,2),R(k,3),R(k,4),...
        mean(R(k:k+4,5)),mean(R(k:k+4,6)));
end
fclose(fid);