function Avg_tsaro_2_norm
load('tsaro_results_2_norm.mat')
[row,col] = size(R);
R = [R,zeros(row,3)];
T = R;
for i = 1:row
    if T(i,6)-T(i,5) < 1e-5      T(i,col+1) = 1; end
    if T(i,10)-T(i,9)< 1e-5      T(i,col+2) = 1; end
    if T(i,14)-T(i,13)<1e-5      T(i,col+3) = 1; end
end
G = [(T(:,6)-T(:,5))./max(1,T(:,5)),(T(:,10)-T(:,9))./max(1,T(:,9)),(T(:,14)-T(:,13))./max(1,T(:,13))]; % relative gap
dir = pwd;
txtname = strcat(dir,'/table10.txt');
fid = fopen(txtname,'a');
fprintf(fid,'------------ Average numerical results of Alg.5, Alg.6 and CCG for TSRLTP with p=2 -------------\r\n\r\n');
fprintf(fid,'\\hline       Size      &      Alg.5              &         Alg.6          &           CCG          \\\\\r\n');
fprintf(fid,'\\hline  m  &  n  &   r &  gap  &  Time  &  iter  &  gap  &  Time  &  iter &  gap  &  Time  &  iter \\\\\r\n');
for k = 1:5:row
    fprintf(fid,'\\hline %3d & %3d & %3d & %3.1e(%d) & %7.3f  & %5.1f& %3.1e(%d) & %7.3f  & %5.1f & %3.1e(%d) & %7.3f  & %5.1f\\\\\r\n',T(k,2),T(k,3),T(k,4),...
        mean(G(k:k+4,1)),sum(T(k:k+4,col+1)),mean(T(k:k+4,7)),mean(T(k:k+4,8)),...
        mean(G(k:k+4,2)),sum(T(k:k+4,col+2)),mean(T(k:k+4,11)),mean(T(k:k+4,12)),...
        mean(G(k:k+4,3)),sum(T(k:k+4,col+3)),mean(T(k:k+4,15)),mean(T(k:k+4,16)));
end
fclose(fid);