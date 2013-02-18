\name{mig_multi_regions}
\alias{mig_multi_regions}
\title{Memory-efficient Implementation of Gabriel et al. 2002 haplotype block definition}
\description{
	Function for the efficient whole-genome haplotype block partitioning.
	It is analogous to \code{\link{mig}} and allows to specify several chromosomal regions at once.
	Haplotype blocks are defined based on D' coefficient of linkage disequilibrium (Gabriel et al., 2002).
}
\usage{
	mig_multi_regions(phase_file, output_files, regions_start, regions_end, processes = 1, 
	phase_file_format = "VCF", map_file = NULL, maf = 0.0, ci_method = "WP",
	l_density = 100, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, 
	ld_fraction = 0.95, pruning_method = "MIG++", windows = NULL)
}
\arguments{
	\item{phase_file}{
		Name of the input file with phased genotypes in VCF or HAPMAP2 format.
	}
	\item{output_files}{
		The list of names of the output files where to store the haplotype blocks.
		One output file for every region.
	}
	\item{regions_start}{
		Numeric vector with start positions (in base-pairs) of the chromosomal regions to be partitioned.
	}
	\item{regions_end}{
		Numeric vector with end positions (in base-pairs) of the chromosomal regions to be partitioned.
	}
	\item{processes}{
		An integer >= 1, which indicates the number of parallel processes. 
		All parallel processes are created on the same machine and shares the main memory.
	}
	\item{phase_file_format}{
		Format of the phase_file: VCF (default) or HAPMAP2.
	}
	\item{map_file}{
		Name of the map file with base-pair positions of each SNP.
		Mandatory when file_format = HAPMAP2.
	}
	\item{maf}{
		Minor Allele Frequency (MAF) threshold: SNPs with MAF <= maf will not be considered.
		The threshold may vary from 0 (default) to 0.5.
	}
	\item{ci_method}{
		Confidence interval (CI) estimation method.
		Supported methods are WP (default) = Wall and Pritchard (2003) method; AV = approximate variance estimator by Zapata et al. (1997).
	}
	\item{l_density}{
		Number of points at which to evaluate the likelihood (applies only to the WP method). 
		Default is 100. 
		The higher the number the longer the runtime. 
		The lower the number the lower the precision.
	}
	\item{ld_ci}{
		Numeric vector with 2 values: thresholds for the lower bound (CL) and upper bound (CU) of the 90\% CI of D'.
		Following Gabriel et al. (2002), default is c(0.7, 0.98).
	}
	\item{ehr_ci}{
		Threshold value for the evidence of historical recombination. 
		Following Gabriel et al. (2002), default is 0.9.
	}
	\item{ld_fraction}{
		Fraction of strong LD SNP pairs over all informative pairs that is needed to classify a sequence of SNP as a haplotype block.
		Following Gabriel et al. (2002), default is 0.95.
	}
	\item{pruning_method}{
		Name of a search space pruning method.
		Supported  methods are MIG, MIG+ and MIG++ (default).
	}
	\item{windows}{
		Numeric vector where every value corresponds to the according region and specifies the number of SNPs within the window in MIG++ search space pruning method.
		If NULL (default), all values are calculated on the fly based on the corresponding region lengths and ld_fraction.
	}
}
\author{Daniel Taliun, Johann Gamper, Cristian Pattaro}
\seealso{
	See \link{mig} for the haplotype block definition, description of the D' distribution modeling and pruning methods.
}
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
	
    # run mig_multi_regions() function on HapMap phase II CEU data with default arguments.
    mig_multi_regions(
     phase_file = "HapMapII_r22_CEU_chr20_31767872_33700401.phase.gz", 
     output_file = c("HapMapII_r22_CEU_chr20_31767872_32767871.blocks.txt", "HapMapII_r22_CEU_chr20_32767872_33700401.blocks.txt"),
     regions_start = c(31767872, 32767872),
     regions_end = c(32767871, 33700401),
     processes = 2,
     phase_file_format = "HAPMAP2",
     map_file = "HapMapII_r22_CEU_chr20_31767872_33700401.legend.txt.gz"
    )
    
    # show contents of the output files
    file.show(
     "HapMapII_r22_CEU_chr20_31767872_32767871.blocks.txt",
     title="HapMapII_r22_CEU_chr20_31767872_32767871.blocks.txt"
    )
    
    file.show(
     "HapMapII_r22_CEU_chr20_32767872_33700401.blocks.txt",
     title="HapMapII_r22_CEU_chr20_32767872_33700401.blocks.txt"
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