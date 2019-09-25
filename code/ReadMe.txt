The codes are used to test the algorithms proposed in the paper entitled "Complexity Results and Effective Algorithms for Worst-case Linear Optimization under Uncertainties" (JOC-2018-08-OA-147.R2), which has been submitted to INFORMS Journal on Computing.

To run the codes successfully, MATLAB, BARON, CVX and CPLEX should be installed in advance.

Beyond that, the path of the interface of BARON, CVX and Cplex for Matlab should be added to the search path of Matlab. Additionally, "CplexLibName" for BARON should be set according to your computer system, which has been marked in "test_wcsr_2_norm.m".


The codes are organized as follows.

1.Solver:

"wcsr_sca_cqr.m" is used to find a globally optimal solution for WCSR with p=2 by SCOBB algorithm, which is Algorithm 2 in the paper.

"wclo_cqr.m" is used to solve the convex quadratic relaxation problem to get the upper bound for SCOBB algorithm.

"wclo_sca.m" is used to get the lower bound for SCOBB algorithm. 

"nlsdp.m" is used to get lower bound and upper bound for WCSR with p=2 by solving an nonlinear semi-definite programming.

"wcsr_q_2_miqcp.m" is used to find a globally optimal solution for WCSR with p=2 by MIPR, which solves a 0-1 mixed integer quadratic constrained programming.

"wclo_q_one.m" is used to find a globally optimal solution for WCSR with p=infty by Algorithm 4 in the paper.

"wcsr_q_one_milp.m" is used to find a globally optimal solution for WCSR with p=infty by MIPR, which solve a 0-1 mixed integer linear programming.

"wclo_q_inf.m" is used to find a globally optimal solution for WCSR with p=1 by Algorithm A in the paper.

"wclo_q_inf2.m" is used to find a globally optimal solution for WCSR with p=1 by Algorithm B in the paper.

"TSARO_CPA.m" is used to find a globally optimal solution for TSARO problems with p=2 or p=infty by a cutting plane algorithm, which is Algorithm 5 in the paper.

"TSARO_CCG.m" is used to find a globally optimal solution for TSARO problems with p=2 or p=infty by the C&CG algorithm framework in Zeng and Zhao (2013). If the parameter "SP2solver" is set to MILP or MIQCP, it is indeed the C&CG algorithm in Zeng and Zhao(2013). If the parameter "SP2solver" is set to FBBA or SCOBB, it is indeed the Algorithm 6 in the paper. See "TSARO_CCG.m" for more details.

"TSARO_SP2_2norm.m" is used to solve the SP2 subproblem with p=2 in C&CG.

"TSARO_SP2_infnorm.m" is used to solve the SP2 subproblem with p=infty in C&CG.

"TSARO_1norm.m" is used to find a globally optimal solution for TSARO problems with p=1 by solving a linear programming, which is Algorithm C in the paper.


2. Data generation:

"generate_wcsr_data.m" is used to generate 60 test problems for WCSR with n = 20,30,50,80,100,150 and r = 5,10.

"generate_wcsr_large_data.m" is used to generate 40 test problems for WCSR with n = 200,300,400,500 and r = 5,10.

"generate_wcsr_real_data.m" is used to generate 10 test problems for WCSR with partial real data.

"generate_TSARO_data.m" is used to generate 65 test problems for TSRLTP.

In order to get better results, some scaling techniques to the orginal data should be done. Indeed, the parameters "xi", "zeta" and "c" has been devided by 10000 for TSARO data. See "generate_TSARO_data.m" for more details.

If you want test these codes, you can run the above 4 m-files to generate the data, or you can download the data we generated in the paper. Then put the data and the codes in the same directory and do as follows.  

3. Test and get the average results:

"test_wcsr_2_norm.m" is used to test the performance of SCOBB, BARON, MIPR and NLSDP for WCSR with p=2 and n<=100. After it is run, you can run "Avg_wcsr_2_norm.m" directly to get table 2 in the paper. The average results will be displayed in LaTex format.
 
"test_wcsr_2_norm_large.m" is used to test the performance of SCOBB for WCSR with p=2 and n>100. After it is run, you can run "Avg_wcsr_2_norm_large.m" directly to get table 4 in the paper. The average results will be displayed in LaTex format.

"test_wcsr_inf_norm.m" is used to test the performance of Algorithm 4 and MIPR for WCSR with p=infty and n<=150. After it is run, you can run "Avg_wcsr_inf_norm.m" directly to get table 5 in the paper. The average results will be displayed in LaTex format.

"test_wcsr_inf_norm_large" is used to test the performance of Algorithm 4 for WCSR with p=infty and n>=200. After it is run, you can run "Avg_wcsr_inf_norm_large.m" directly to get table 6 in the paper. The average results will be displayed in LaTex format.

"test_wcsr_1_norm_large.m" is used to test the performance of Algorithm A and B for WCSR with p=1 and n>100. After it is run, you can run "Avg_wcsr_1_norm_large.m" directly to get table 7 in the paper. The average results will be displayed in LaTex format.

"test_wcsr_real_data.m" is used to test the performance of SCOBB, Algorithm 4 and Algorithm A for WCSR with partial real data. Run "test_wcsr_real_data(1)" to get the results for December 2009 and "test_wcsr_real_data(2)" to get the results for June 2010. Then you can run "Avg_wcsr_real_data.m" directly to get table 9 in the paper. The average results will be displayed in LaTex format.

"test_TSARO_2_norm.m" is used to test the performance of Algorithm 5, Algorithm 6 and C&CG for TSRLTP with p=2. After it is run, you can run "Avg_tsaro_2_norm.m" directly to get table 10 in the paper. The average results will be displayed in LaTex format.

"test_TSARO_inf_norm.m" is used to test the performance of Algorithm 5, Algorithm 6 and C&CG for TSRLTP with p=infinity. After it is run, you can run "Avg_tsaro_inf_norm.m" directly to get table 11 in the paper. The average results will be displayed in LaTex format.

"test_TSARO_1_norm.m" is used to test the performance of Algorithm C for TSRLTP with p=1. After it is run, you can run "Avg_tsaro_1_norm.m" directly to get table 12 in the paper. The average results will be displayed in LaTex format.

"table3.m" is used to list the three specifical instances of WCSR with p=2 to get table 3.

Most of the m-files have necessary annotation, see them for more details.



