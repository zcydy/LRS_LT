The project is for the "paper Lunar Radar Sounder-Based Detection of Lunar Subsurface Anomalies from Time-Frequency and Robust Principal Component Analysis Approaches"
Currently, the project is at first version, the readability and efficiency requires further modification. 


There are 3 .txt files contains gprMax .in file, all  simulation parameters are included

1.rough_3d_lt is a 3d simulation of a lava tube underneath multiple craters(sim-rough_lt.m).
2.rough_lt is a 2d simulation of a lava tube underneath multiple craters.
3.rough_flat is a 2d simulation of a flat surface with craters and slop(sim-rough.m). 

There is a main.m file, 9 .mat  files , and a RPCA pack.


main.m file is all in one file that contains plot generation and lava tube detection. It can be directly run. 
1.echo_ln is the linear scaled LRS data read from LRS_SAR10KM_20080824184044
2.echo_db is the dB scaled LRS data read from LRS_SAR10KM_20080824184044
3.sim_rough.mat was introduced above 
4.sim_rough_lt.mat was introduced above 
5.sparsity.mat is the results of sparsity test(step1)
6.ranky.mat  is the results of rank test(step2)
7.ano.mat is the results of RPCA test(step3)
8.lon. mat is longitude information
9.lat. mat is latitude information (for marking the lat and lon in figure3 (d) )


RPCA pack contains ( code below is Copyright 2009 Perception and Decision Lab, University of Illinois at Urbana-Champaign, and Microsoft Research Asia, Beijing. ):
1.Propack： prerequisite functions 
2.DS. store file 
3.Choosesvd: effecient svd algorithms
4.inexact_alm_rpaca: RPCA used in the main.m
(special thanks to Z. Lin, M. Chen, L. Wu, and Y. Ma to share their code)

