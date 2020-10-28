The 'blocks' and 'genes' folders both contain the following folders:

  -- raxml contains the separate output trees from RAxML before any preprocessing and folders with additional output files from raxml that were not used in this project
  -- original-rooted contains the original tree sets (all binary) after rooting, including branch lengths and bootstrap support values
  -- nonbinary-original contains the not necessarily binary tree sets after edge contraction
  -- nonbinary-nexus-mpl contains the same tree sets as nonbinary-original, but then in NEXUS format, which were used as input for Maximum Pseudo-Likelihood algorithm
  -- nonbinary-dated-icass contains the same tree sets as nonbinary-original, but then with collection dates and without bootstrap support values, as suitable for ICass
  -- binary-treechild-semitemp contains binary tree sets, which were used as input for the Tree-Child and Semi-Temporal algorithms. Files with 'nolabels' in their name are the same tree sets as original-rooted, but without branch lengths and bootstrap support values. Files with 'binary' in their name resulted from resolving multifurcations using the Data Transformer for Real-World Data with the nonbinary-original tree sets as input
  -- additional contains binary tree sets resulting from using original-rooted tree sets as input for the Data Transformer for Real-World Data, resulting in a set of only unique trees (without branch lengths and bootstrap support values). These were not used as input for any of the algorithms, but could have been used instead of the files with 'nolabels' in their name in the binary-treechild-semitemp folder 
  -- data-transformer-logs contains logs of the output from the Data Transformer for Real-World Data, contanining for example the number of input and output trees and the number of resolved multifurcations
  
Files are named as follows:
  -- blocks/genes indicate the breakpoints
  -- A, B, C, Aminus, Bminus, Cminus indicate taxon selections
  -- len indicates branch length threshold for edge contraction (if applied)
  -- supp indicates bootstrap support value threshold for edge contraction (if applied)
