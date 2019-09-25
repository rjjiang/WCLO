function Avg_wcsr_real_data
% To get table 9
load('wcsr_real_data_Dec_2009_results.mat')
dir = pwd;
txtname = strcat(dir,'/table9.txt');
fid = fopen(txtname,'a');
fprintf(fid,'-------------- Average performance of SCOBB, Alg. 4 and Alg. A for WCSR with partial real data    ---------------\r\n\r\n');
fprintf(fid,'\\hline                            SCOBB       &              Alg. 4        &            Alg. A           \\\\\r\n');
fprintf(fid,'\\hline              Opt.val  &  Time  & avg.x &   Opt.val &  Time  & avg.x &   Opt.val &  Time  & avg.x  \\\\\r\n');
fprintf(fid,'\\hline Dec. 2009 & %9.6f & %5.3f  & %5.3f & %9.6f & %5.3f  & %5.3f & %9.6f & %5.3f  & %5.3f  \\\\\r\n',...
        mean(R(:,5)),mean(R(:,6)),mean(R(:,7)),mean(R(:,8)),mean(R(:,9)),mean(R(:,10)),mean(R(:,2)),mean(R(:,3)),mean(R(:,4)));
    
load('wcsr_real_data_Jun_2010_results.mat')
fprintf(fid,'\\hline Jun. 2010 & %9.6f & %5.3f  & %5.3f & %9.6f & %5.3f  & %5.3f & %9.6f & %5.3f  & %5.3f  \\\\\r\n',...
        mean(R(:,5)),mean(R(:,6)),mean(R(:,7)),mean(R(:,8)),mean(R(:,9)),mean(R(:,10)),mean(R(:,2)),mean(R(:,3)),mean(R(:,4)));

fclose(fid);