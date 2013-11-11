\name{mig_rsq}
\alias{mig_rsq}
\title{Memory-efficient Implementation of Gabriel et al. 2002 haplotype block definition using the r^2 coefficient}
\description{
	Function for the efficient whole-genome haplotype block partitioning.
	Haplotype blocks are defined based on r^2 coefficient of linkage disequilibrium.
}
\usage{
	mig_rsq(phase_file, output_file, phase_file_format = "VCF", 
	map_file = NULL, region = NULL, maf = 0.0, 
	weak_rsq = 0.5, strong_rsq = 0.8, fraction = 0.95, pruning_method = "MIG++", window = NULL)
}
\arguments{
	\item{phase_file}{
		Name of the input file with phased genotypes in VCF or HAPMAP2 format.
	}
	\item{output_file}{
		Name of the output file where to store the haplotype blocks.
	}
	\item{phase_file_format}{
		Format of the phase_file: VCF (default) or HAPMAP2.
		If VCF, then only SNPs with "PASS" or "." in the FILTER field are considered.
	}
	\item{map_file}{
		Name of the map file with base-pair positions of each SNP.
		Mandatory when file_format = HAPMAP2.
	}
	\item{region}{
		Numeric vector with start and end positions (in base-pairs) of the chromosomal region to be partitioned.
		If NULL (default), then the whole chromosome is processed.
	}
	\item{maf}{
		Minor Allele Frequency (MAF) threshold: SNPs with MAF <= maf will not be considered.
		The threshold may vary from 0 (default) to 0.5.
	}
	\item{weak_rsq}{
		r^2 threshold value for the SNP pairs in weak LD.
		A SNP pair is in weak LD if r^2 < weak_rsq.
		By default, weak_rsq = 0.5.
	}
	\item{strong_rsq}{
		r^2 threshold value for the SNP pairs in strong LD.
		A SNP pair is in strong LD if r^2 >= strong_rsq.
		By default, strong_rsq = 0.8.
	}
	\item{fraction}{
		Fraction of strong LD SNP pairs over all informative SNP pairs that is needed to classify a sequence of SNP as a haplotype block.
		Informative SNP pairs are those that are in strong LD or in weak LD.
		Default is 0.95.
	}
	\item{pruning_method}{
		Name of a search space pruning method.
		Supported  methods are MIG, MIG+ and MIG++ (default).
	}
	\item{window}{
		Number of SNPs within the window in MIG++ search space pruning method.
		If NULL (default), it is calculated on the fly based on the region length and ld_fraction.
	}
}
\section{Haplotype Blocks}{
	The haplotype blocks are defined based on r^2 coefficient of linkage disequilibrium (LD) between a pair of SNPs following the logic suggested by Gabriel et al., 2002.
	
	\subsection{SNP Pairs Classification}{
		According to  the r^2 statistic, LD between any SNP pair is classified as follows:  \cr 
		(1) strong LD if r^2 >= 0.8;  \cr
		(2) weak LD if r^2 < 0.5;
		
		All SNP pairs satisfying conditions (1) and (2) are called informative.
	}

	\subsection{Haplotype Block Definition}{
		A chromosomal region is a haplotype block if:  \cr
		(a) the two outer-most SNPs are in strong LD;  \cr
		(b) at least 95\% of informative SNP pairs within the region are in strong LD.  
	}
}
\section{Search Space Pruning}{
	\subsection{MIG}{
		The MIG is a memory efficient implementation of Gabriel et al. 2002 haplotype block definition using the r^2 coefficient.  
		It has extremely low memory requirements (linear memory complexity) and is able to represent LD structure of the entire chromosome (millions of SNPs) in the main memory.  
		The MIG doesn't prune the search space i.e. it computes LD for all possible SNP pairs.     
	}
	
	\subsection{MIG+}{
		The MIG+ identifies and omits computations in regions that can never become haplotype blocks.  
		The chromosome is processed from left to right.  
		Assuming all unprocessed SNP pairs to the right being in strong LD, an early termination point is determined using already known LD structure to the left.  
		It is theoretically guaranteed, that haplotype blocks produced by MIG and MIG+ are identical.  
		The MIG+ has low memory requirements and improved runtime compared to MIG.  
	}
	
	\subsection{MIG++}{
		The MIG++ is the improved version of MIG+ search space pruning method.  
		The chromosome is processed in several iterations.  
		At every iteration it is scanned from the very start to the very end and the computations are restricted by a user-defined window.  
		The window is expressed as a maximal number of SNPs that will be considered when computing LD with the currently considered SNP.  
		In every new iteration, the partially calculated LD structure across the whole chromosome allows to predict termination points better.  
		It is theoretically guaranteed, that haplotype blocks produced by MIG, MIG+ and MIG++ are identical.  
		The MIG++ has low memory requirements and significantly better runtime compared to MIG and MIG+.  
	}
}
\section{Output File}{
	The output file consists of the following columns:
	\tabular{ll}{
		BLOCK_NAME \tab Generated unique block name\cr
		FIRST_SNP \tab Name of the first SNP in block\cr
		LAST_SNP \tab Name of the last SNP in block\cr
		FIRST_SNP_ID \tab Index of the first SNP in block with respect to the filtered SNPs\cr
		LAST_SNP_ID \tab Index of the last SNP in block with respect to the filtered SNPs\cr
		START_BP \tab The base-pair position of the first SNP in block\cr
		END_BP \tab The base-pair position of the last SNP in block\cr
		N_SNPS \tab Number of SNPs in block\cr
		N_HAPS \tab Number of haplotypes in block\cr
		N_UNIQUE_HAPS \tab Number of unique haplotypes in block\cr
		N_COMMON_HAPS \tab Number of common (which appear more than once) haplotypes in block\cr
		N_HAPS_DIVERSITY \tab The haplotype diversity in block (Patil et al., 2001). 1 - low diversity, 0 - high diversity. 
	}
}
\references{
	Patil, N. et al. (2001) Blocks of Limited Haplotype Diversity Revealed by High-Resolution Scanning of Human Chromosome 21. \emph{Science}, \bold{294}(5547), 1719--1723.
	
	Gabriel, S. B. et al. (2002) The Structure of Haplotype Blocks in the Human Genome. \emph{Science}, \bold{296}(5576), 2225--2229.
}
\author{Daniel Taliun, Johann Gamper, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{package}
\examples{
\dontshow{
    # change the workspace
    currentWd <- getwd()
    newWd <- paste(system.file(package="LDExplorer"), "doc", sep="/")
    setwd(newWd)
}
	
    # load LDExplorer library
    library(LDExplorer)
	
    # run mig_rsq() function on HapMap phase II CEU data with default arguments.
    mig_rsq(
     phase_file = "HapMapII_r22_CEU_chr20_31767872_33700401.phase.gz", 
     output_file = "HapMapII_r22_CEU_chr20_31767872_33700401.rsq_blocks.txt",
     phase_file_format = "HAPMAP2",
     map_file = "HapMapII_r22_CEU_chr20_31767872_33700401.legend.txt.gz"
    )
	
    # show contents of the output file
    file.show(
     "HapMapII_r22_CEU_chr20_31767872_33700401.rsq_blocks.txt",
     title="HapMapII_r22_CEU_chr20_31767872_33700401.rsq_blocks.txt"
    )
	
    # run mig_rsq() function on 1000 Genomes Project CEU data with default arguments.
    mig_rsq(
     phase_file = "1000G_phase1_v3_20101123_CEU_chr2_89153688_89307566.vcf.gz", 
     output_file = "1000G_phase1_v3_20101123_CEU_chr2_89153688_89307566.rsq_blocks.txt"
    )
	
    # show contents of the output file
    file.show(
     "1000G_phase1_v3_20101123_CEU_chr2_89153688_89307566.rsq_blocks.txt",
     title="1000G_phase1_v3_20101123_CEU_chr2_89153688_89307566.rsq_blocks.txt"
    )
	
\dontshow{
    # restore previous workspace
    setwd(currentWd)
    
    # all input and output files are located in the subdirectory "doc" of the installed LDExplorer package
    message <- c("\n", rep("#", 40), "\n")
    message <- c(message, "\nAll input and output files of this example are located in directory:\n", newWd, "\n")
    message <- c(message, "\n", rep("#", 40),"\n")
    cat(message, sep="")
}
}