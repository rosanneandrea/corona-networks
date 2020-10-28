

Using BFGS method to optimize GTR rate parameters, to disable this specify "--no-bfgs" 



This is RAxML version 8.2.12 released by Alexandros Stamatakis on May 2018.

With greatly appreciated code contributions by:
Andre Aberer      (HITS)
Simon Berger      (HITS)
Alexey Kozlov     (HITS)
Kassian Kobert    (HITS)
David Dao         (KIT and HITS)
Sarah Lutteropp   (KIT and HITS)
Nick Pattengale   (Sandia)
Wayne Pfeiffer    (SDSC)
Akifumi S. Tanabe (NRIFS)
Charlie Taylor    (UF)


Alignment has 20 distinct alignment patterns

Proportion of gaps and completely undetermined characters in this alignment: 0.87%

RAxML rapid bootstrapping and subsequent ML search

Using 1 distinct models/data partitions with joint branch length optimization



Executing 100 rapid bootstrap inferences and thereafter a thorough ML search 

All free model parameters will be estimated by RAxML
ML estimate of 25 per site rate categories

Likelihood of final tree will be evaluated and optimized under GAMMA

GAMMA Model parameters will be estimated up to an accuracy of 0.1000000000 Log Likelihood units

Partition: 0
Alignment Patterns: 20
Name: No Name Provided
DataType: DNA
Substitution Matrix: GTR




RAxML was called as follows:

./raxmlHPC-PTHREADS-AVX -f a -m GTRCAT -s ../../data/Corona_edited_seqsel10.ph.genes_seqsel10.E.txt -n genes_seqsel10.E -p 123 -T 8 -x 123 -# 100 -w /home/rosanne/Documents/minor_research_project/results/raxml 



Time for BS model parameter optimization 0.011511
Bootstrap[0]: Time 0.020460 seconds, bootstrap likelihood -388.696275, best rearrangement setting 5
Bootstrap[1]: Time 0.023492 seconds, bootstrap likelihood -358.085064, best rearrangement setting 6
Bootstrap[2]: Time 0.008827 seconds, bootstrap likelihood -363.171290, best rearrangement setting 7
Bootstrap[3]: Time 0.017976 seconds, bootstrap likelihood -365.154273, best rearrangement setting 11
Bootstrap[4]: Time 0.015847 seconds, bootstrap likelihood -362.811665, best rearrangement setting 15
Bootstrap[5]: Time 0.013347 seconds, bootstrap likelihood -352.963124, best rearrangement setting 15
Bootstrap[6]: Time 0.011772 seconds, bootstrap likelihood -359.632130, best rearrangement setting 5
Bootstrap[7]: Time 0.008779 seconds, bootstrap likelihood -334.566191, best rearrangement setting 12
Bootstrap[8]: Time 0.008477 seconds, bootstrap likelihood -346.404529, best rearrangement setting 11
Bootstrap[9]: Time 0.014749 seconds, bootstrap likelihood -386.516043, best rearrangement setting 5
Bootstrap[10]: Time 0.016990 seconds, bootstrap likelihood -364.192273, best rearrangement setting 13
Bootstrap[11]: Time 0.006759 seconds, bootstrap likelihood -348.423080, best rearrangement setting 7
Bootstrap[12]: Time 0.006772 seconds, bootstrap likelihood -335.116798, best rearrangement setting 6
Bootstrap[13]: Time 0.005366 seconds, bootstrap likelihood -340.264012, best rearrangement setting 5
Bootstrap[14]: Time 0.008559 seconds, bootstrap likelihood -359.265271, best rearrangement setting 5
Bootstrap[15]: Time 0.014154 seconds, bootstrap likelihood -337.744494, best rearrangement setting 8
Bootstrap[16]: Time 0.005121 seconds, bootstrap likelihood -351.147565, best rearrangement setting 14
Bootstrap[17]: Time 0.006668 seconds, bootstrap likelihood -350.907096, best rearrangement setting 14
Bootstrap[18]: Time 0.008535 seconds, bootstrap likelihood -354.933975, best rearrangement setting 12
Bootstrap[19]: Time 0.018968 seconds, bootstrap likelihood -361.745237, best rearrangement setting 5
Bootstrap[20]: Time 0.012795 seconds, bootstrap likelihood -367.615816, best rearrangement setting 5
Bootstrap[21]: Time 0.005418 seconds, bootstrap likelihood -355.149580, best rearrangement setting 7
Bootstrap[22]: Time 0.007116 seconds, bootstrap likelihood -357.038726, best rearrangement setting 9
Bootstrap[23]: Time 0.021247 seconds, bootstrap likelihood -356.059264, best rearrangement setting 14
Bootstrap[24]: Time 0.007181 seconds, bootstrap likelihood -346.991703, best rearrangement setting 9
Bootstrap[25]: Time 0.015349 seconds, bootstrap likelihood -372.035637, best rearrangement setting 8
Bootstrap[26]: Time 0.006203 seconds, bootstrap likelihood -352.116633, best rearrangement setting 12
Bootstrap[27]: Time 0.007274 seconds, bootstrap likelihood -344.033591, best rearrangement setting 8
Bootstrap[28]: Time 0.006741 seconds, bootstrap likelihood -362.338801, best rearrangement setting 11
Bootstrap[29]: Time 0.005936 seconds, bootstrap likelihood -355.428773, best rearrangement setting 15
Bootstrap[30]: Time 0.028824 seconds, bootstrap likelihood -352.788557, best rearrangement setting 8
Bootstrap[31]: Time 0.006240 seconds, bootstrap likelihood -360.002929, best rearrangement setting 5
Bootstrap[32]: Time 0.008647 seconds, bootstrap likelihood -354.925564, best rearrangement setting 10
Bootstrap[33]: Time 0.011019 seconds, bootstrap likelihood -343.665475, best rearrangement setting 11
Bootstrap[34]: Time 0.011129 seconds, bootstrap likelihood -377.118831, best rearrangement setting 14
Bootstrap[35]: Time 0.014848 seconds, bootstrap likelihood -339.569928, best rearrangement setting 10
Bootstrap[36]: Time 0.006889 seconds, bootstrap likelihood -352.558625, best rearrangement setting 8
Bootstrap[37]: Time 0.011559 seconds, bootstrap likelihood -366.913734, best rearrangement setting 6
Bootstrap[38]: Time 0.008732 seconds, bootstrap likelihood -363.929208, best rearrangement setting 10
Bootstrap[39]: Time 0.008248 seconds, bootstrap likelihood -367.188477, best rearrangement setting 5
Bootstrap[40]: Time 0.020821 seconds, bootstrap likelihood -346.688065, best rearrangement setting 6
Bootstrap[41]: Time 0.014829 seconds, bootstrap likelihood -338.149413, best rearrangement setting 9
Bootstrap[42]: Time 0.009961 seconds, bootstrap likelihood -362.866073, best rearrangement setting 14
Bootstrap[43]: Time 0.009518 seconds, bootstrap likelihood -372.704765, best rearrangement setting 10
Bootstrap[44]: Time 0.017005 seconds, bootstrap likelihood -353.188867, best rearrangement setting 8
Bootstrap[45]: Time 0.009891 seconds, bootstrap likelihood -332.469285, best rearrangement setting 14
Bootstrap[46]: Time 0.008552 seconds, bootstrap likelihood -351.032203, best rearrangement setting 14
Bootstrap[47]: Time 0.009953 seconds, bootstrap likelihood -349.432970, best rearrangement setting 7
Bootstrap[48]: Time 0.021034 seconds, bootstrap likelihood -356.498960, best rearrangement setting 12
Bootstrap[49]: Time 0.034061 seconds, bootstrap likelihood -372.934884, best rearrangement setting 11
Bootstrap[50]: Time 0.004381 seconds, bootstrap likelihood -342.762990, best rearrangement setting 12
Bootstrap[51]: Time 0.006269 seconds, bootstrap likelihood -375.291318, best rearrangement setting 13
Bootstrap[52]: Time 0.014078 seconds, bootstrap likelihood -369.194861, best rearrangement setting 7
Bootstrap[53]: Time 0.008305 seconds, bootstrap likelihood -347.149828, best rearrangement setting 6
Bootstrap[54]: Time 0.005475 seconds, bootstrap likelihood -362.081988, best rearrangement setting 10
Bootstrap[55]: Time 0.008469 seconds, bootstrap likelihood -338.336477, best rearrangement setting 11
Bootstrap[56]: Time 0.004978 seconds, bootstrap likelihood -369.629769, best rearrangement setting 12
Bootstrap[57]: Time 0.008981 seconds, bootstrap likelihood -339.919945, best rearrangement setting 9
Bootstrap[58]: Time 0.006748 seconds, bootstrap likelihood -360.771083, best rearrangement setting 11
Bootstrap[59]: Time 0.004204 seconds, bootstrap likelihood -338.209062, best rearrangement setting 8
Bootstrap[60]: Time 0.007735 seconds, bootstrap likelihood -361.635833, best rearrangement setting 11
Bootstrap[61]: Time 0.015959 seconds, bootstrap likelihood -371.807528, best rearrangement setting 12
Bootstrap[62]: Time 0.014673 seconds, bootstrap likelihood -342.794520, best rearrangement setting 15
Bootstrap[63]: Time 0.011056 seconds, bootstrap likelihood -347.783654, best rearrangement setting 6
Bootstrap[64]: Time 0.006755 seconds, bootstrap likelihood -344.447933, best rearrangement setting 15
Bootstrap[65]: Time 0.007458 seconds, bootstrap likelihood -335.112232, best rearrangement setting 15
Bootstrap[66]: Time 0.005891 seconds, bootstrap likelihood -359.848577, best rearrangement setting 9
Bootstrap[67]: Time 0.005666 seconds, bootstrap likelihood -380.566181, best rearrangement setting 7
Bootstrap[68]: Time 0.017271 seconds, bootstrap likelihood -379.138237, best rearrangement setting 6
Bootstrap[69]: Time 0.016999 seconds, bootstrap likelihood -366.065406, best rearrangement setting 10
Bootstrap[70]: Time 0.009171 seconds, bootstrap likelihood -338.471921, best rearrangement setting 15
Bootstrap[71]: Time 0.009642 seconds, bootstrap likelihood -365.435565, best rearrangement setting 8
Bootstrap[72]: Time 0.005105 seconds, bootstrap likelihood -357.837812, best rearrangement setting 8
Bootstrap[73]: Time 0.008547 seconds, bootstrap likelihood -370.483363, best rearrangement setting 14
Bootstrap[74]: Time 0.007647 seconds, bootstrap likelihood -354.099971, best rearrangement setting 11
Bootstrap[75]: Time 0.007419 seconds, bootstrap likelihood -353.659538, best rearrangement setting 11
Bootstrap[76]: Time 0.009283 seconds, bootstrap likelihood -343.606196, best rearrangement setting 7
Bootstrap[77]: Time 0.004517 seconds, bootstrap likelihood -347.867376, best rearrangement setting 7
Bootstrap[78]: Time 0.006028 seconds, bootstrap likelihood -337.360160, best rearrangement setting 15
Bootstrap[79]: Time 0.006988 seconds, bootstrap likelihood -363.032507, best rearrangement setting 14
Bootstrap[80]: Time 0.004135 seconds, bootstrap likelihood -342.418560, best rearrangement setting 13
Bootstrap[81]: Time 0.023846 seconds, bootstrap likelihood -365.028072, best rearrangement setting 12
Bootstrap[82]: Time 0.004218 seconds, bootstrap likelihood -377.807135, best rearrangement setting 10
Bootstrap[83]: Time 0.004981 seconds, bootstrap likelihood -361.127793, best rearrangement setting 8
Bootstrap[84]: Time 0.013919 seconds, bootstrap likelihood -360.181903, best rearrangement setting 15
Bootstrap[85]: Time 0.005366 seconds, bootstrap likelihood -347.621346, best rearrangement setting 13
Bootstrap[86]: Time 0.006525 seconds, bootstrap likelihood -356.515424, best rearrangement setting 10
Bootstrap[87]: Time 0.004968 seconds, bootstrap likelihood -332.830584, best rearrangement setting 11
Bootstrap[88]: Time 0.005901 seconds, bootstrap likelihood -366.587755, best rearrangement setting 11
Bootstrap[89]: Time 0.007345 seconds, bootstrap likelihood -354.932007, best rearrangement setting 6
Bootstrap[90]: Time 0.017396 seconds, bootstrap likelihood -328.951564, best rearrangement setting 10
Bootstrap[91]: Time 0.005328 seconds, bootstrap likelihood -345.245663, best rearrangement setting 9
Bootstrap[92]: Time 0.005328 seconds, bootstrap likelihood -373.858548, best rearrangement setting 15
Bootstrap[93]: Time 0.005722 seconds, bootstrap likelihood -337.003076, best rearrangement setting 14
Bootstrap[94]: Time 0.004602 seconds, bootstrap likelihood -352.865433, best rearrangement setting 6
Bootstrap[95]: Time 0.004689 seconds, bootstrap likelihood -361.662227, best rearrangement setting 9
Bootstrap[96]: Time 0.005273 seconds, bootstrap likelihood -346.065442, best rearrangement setting 15
Bootstrap[97]: Time 0.014496 seconds, bootstrap likelihood -349.937725, best rearrangement setting 6
Bootstrap[98]: Time 0.006561 seconds, bootstrap likelihood -363.412422, best rearrangement setting 6
Bootstrap[99]: Time 0.007091 seconds, bootstrap likelihood -372.038713, best rearrangement setting 12


Overall Time for 100 Rapid Bootstraps 1.031365 seconds
Average Time per Rapid Bootstrap 0.010314 seconds

Starting ML Search ...

Fast ML optimization finished

Fast ML search Time: 0.222527 seconds

Slow ML Search 0 Likelihood: -414.067696
Slow ML Search 1 Likelihood: -414.067696
Slow ML Search 2 Likelihood: -414.067696
Slow ML Search 3 Likelihood: -414.067696
Slow ML Search 4 Likelihood: -414.067696
Slow ML Search 5 Likelihood: -414.067696
Slow ML Search 6 Likelihood: -414.067696
Slow ML Search 7 Likelihood: -414.067696
Slow ML Search 8 Likelihood: -414.067696
Slow ML Search 9 Likelihood: -414.067696
Slow ML optimization finished

Slow ML search Time: 0.343235 seconds
Thorough ML search Time: 0.010203 seconds

Final ML Optimization Likelihood: -414.067696

Model Information:

Model Parameters of Partition 0, Name: No Name Provided, Type of Data: DNA
alpha: 0.020000
Tree-Length: 0.096679
rate A <-> C: 0.000100
rate A <-> G: 7.981953
rate A <-> T: 3.639653
rate C <-> G: 0.000100
rate C <-> T: 4.089802
rate G <-> T: 1.000000

freq pi(A): 0.221252
freq pi(C): 0.201601
freq pi(G): 0.189956
freq pi(T): 0.387191


ML search took 0.577062 secs or 0.000160 hours

Combined Bootstrap and ML search took 1.608490 secs or 0.000447 hours

Drawing Bootstrap Support Values on best-scoring ML tree ...



Found 1 tree in File /home/rosanne/Documents/minor_research_project/results/raxml/RAxML_bestTree.genes_seqsel10.E



Found 1 tree in File /home/rosanne/Documents/minor_research_project/results/raxml/RAxML_bestTree.genes_seqsel10.E

Program execution info written to /home/rosanne/Documents/minor_research_project/results/raxml/RAxML_info.genes_seqsel10.E
All 100 bootstrapped trees written to: /home/rosanne/Documents/minor_research_project/results/raxml/RAxML_bootstrap.genes_seqsel10.E

Best-scoring ML tree written to: /home/rosanne/Documents/minor_research_project/results/raxml/RAxML_bestTree.genes_seqsel10.E

Best-scoring ML tree with support values written to: /home/rosanne/Documents/minor_research_project/results/raxml/RAxML_bipartitions.genes_seqsel10.E

Best-scoring ML tree with support values as branch labels written to: /home/rosanne/Documents/minor_research_project/results/raxml/RAxML_bipartitionsBranchLabels.genes_seqsel10.E

Overall execution time for full ML analysis: 1.611384 secs or 0.000448 hours or 0.000019 days

