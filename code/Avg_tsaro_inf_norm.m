function Avg_tsaro_inf_norm
load('tsaro_results_inf_norm.mat')
[row,col] = size(R);
R = [R,zeros(row,3)];
T = R;
for i = 1:row
    if T(i,6)-T(i,5) < 1e-5      T(i,col+1) = 1; end
    if T(i,10)-T(i,9)< 1e-5      T(i,col+2) = 1; end
    if T(i,14)-T(i,13)<1e-5      T(i,col+3) = 1; end
end
dir = pwd;
txtname = strcat(dir,'/table11.txt');
fid = fopen(txtname,'a');
fprintf(fid,'--------------------------------   compared with Alg.5, Alg.6 and CCG    ---------------------------------\r\n\r\n');
fprintf(fid,'\\hline                       Alg.5                   &         Alg.6              &           CCG              \\\\\r\n');
fprintf(fid,'\\hline  m  &  n  &   r  &   ral. gap     &    Time   &  iter  &   ral. gap   &    Time   &  iter &   ral. gap   &    Time   &  iter \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline %4d & %4d & %4d & %3.1e(%d) &  %7.3f  & %5.1f& %3.1e(%d) & %7.3f  & %5.1f & %3.1e(%d) & %7.3f  & %5.1f\\\\\r\n',T(k,2),T(k,3),T(k,4),...
        mean((T(k:k+4,6)-T(k:k+4,5))./max(1,T(k:k+4,5))),sum(T(k:k+4,col+1)),mean(T(k:k+4,7)),mean(T(k:k+4,8)),...
        mean((T(k:k+4,10)-T(k:k+4,9))./max(1,T(k:k+4,9))),sum(T(k:k+4,col+2)),mean(T(k:k+4,11)),mean(T(k:k+4,12)),...
        mean((T(k:k+4,14)-T(k:k+4,13))./max(1,T(k:k+4,13))),sum(T(k:k+4,col+3)),mean(T(k:k+4,15)),mean(T(k:k+4,16)));
end
fclose(fid);