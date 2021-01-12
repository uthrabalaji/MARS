========================================================================
Demo3: Normalization & Differential expression analysis - RNaseR dataset
========================================================================

Step1 
-----

Read in demo data for ribominus/RNaseR- BSJ ; ribominus/RNaseR- gene counts ; ribinominus/RNaseR+ BSJ::   
	
	setwd("/test_repo/demo_data/demo3_data")

	rminus.bsj <- read.csv("demo3_rminus_BSJ.csv",row.names = 1)

	rplus.bsj  <- read.csv("demo3_rplus_BSJ.csv",row.names = 1)

	rminus.gc  <- read.csv("demo3_rminus_GC.csv",row.names = 1)
	
Read in design file::

	designFile <- read.csv("demo3_design.csv")

**Note:** The design file must contain 2 columns called **'Sample'** (matched with column names in data files) and **'Group'** (conditions/grouping for all samples in the data)


Step2
-----

Load MARS::

	setwd("/test_repo")

	source("runMARS.R")

Run MARS method::

	mars.out <- runMARS(Rminus_BSJ = rminus.bsj, Rplus_BSJ = rplus.bsj, Rminus_GC = rminus.gc, design_file = designFile)


Step3
-----

Write MARS output::

	setwd("_path_to_output")

	write.csv(mars.out,file="mars_results.csv",row.names=F)


MARS output format
------------------

The output file of MARS contains the follwing columns
  	
	+--------+----------------------------------+------------------------------+
	| column |        name        		    |     description              |
	+========+==================================+==============================+
      	|   1    |       circRNA      		    |     circular RNA             |
	+--------+----------------------------------+------------------------------+
 	|   2 	 |Ribominus/RNaseR- baseMean        | base mean in totalRNA  	   |
        +--------+----------------------------------+------------------------------+
	|   3    |Ribominus/RNaseR- log2FoldChange  | log2 fold change in total RNA|
        +--------+----------------------------------+------------------------------+
 	|   4    |Ribominus/RNaseR- pvalue          | p-value in total RNA         |
        +--------+----------------------------------+------------------------------+
	|   5    |Ribominus/RNaseR+ baseMean        | base mean in RNaseR          |
        +--------+----------------------------------+------------------------------+
        |   6    |Ribominus/RNaseR+ log2FoldChange  | log2 fold change in RNaseR   |
        +--------+----------------------------------+------------------------------+
        |   7    |Ribominus/RNaseR+ pvalue          | p-value in RNaseR            |
        +--------+----------------------------------+------------------------------+
	|   8    |meta.pvalue		            | meta analysis p-value        |
        +--------+----------------------------------+------------------------------+
