Folders contain output files (or files created from output) from the different phylogenetic network algorithms

Files are named as follows:
  -- blocks/genes indicate the breakpoints
  -- A, B, C, Aminus, Bminus, Cminus indicate taxon selections
  -- len indicates branch length threshold for edge contraction (if applied)
  -- supp indicates bootstrap support value threshold for edge contraction (if applied)
  
Different file types:
  -- .log is a log of the output from the algorithm
  -- .net or .eNewick is the constructed network in eNewick format
  -- .png is the visualisation by Dendroscope of the constructed network (note that not all networks could be visualised correctly with Dendroscope!)
  -- .dot is an additional output file from TriLoNet which contains the constructed network to visualise in GraphViz
  -- .tnet are the trinets extracted from the TriLoNet output, which were used as input for TriL2Net
