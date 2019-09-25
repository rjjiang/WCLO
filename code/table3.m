function table3
load('wcsr_2_norm_results.mat')
dir = pwd;
txtname = strcat(dir,'/table3.txt');
fid = fopen(txtname,'a');
fprintf(fid,'-------------- Numerical results of SCOBB and NLSDP for three specifical instances of WCSR with p =2 --------------\r\n\r\n');
fprintf(fid,'\\hline      Instance   &         SCOBB            &              NLSDP            \\\\\r\n');
fprintf(fid,'\\hline  ID &  n  &  r  & Time  & Opt.val   & Iter &  Time &    LB      &    UB    \\\\\r\n');
for i = [11 26 45]
    fprintf(fid,'\\hline %3d & %3d & %3d & %5.1f & %9.6f & %4d & %5.1f & %9.6f  & %9.6f\\\\\r\n',...
        R(i,1),R(i,2),R(i,4),R(i,11),R(i,9),R(i,13),R(i,19),R(i,16),R(i,17));
end
fclose(fid)