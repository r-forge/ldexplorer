\name{mig_multi_chr}
\alias{mig_multi_chr}
\title{Memory-efficient Implementation of Gabriel et al. 2002 haplotype block definition}
\description{
	Function for the efficient whole-genome haplotype block partitioning.
	It is analogous to \code{\link{mig}} and allows to specify several input files (chromosomes) at once and process them in parallel.
	Haplotype blocks are defined based on D' coefficient of linkage disequilibrium (Gabriel et al., 2002).
}
\usage{
	mig_multi_chr(phase_files, output_files, processes = 1, phase_file_format = "VCF", 
	map_files = NULL, maf = 0.0, ci_method = "WP",
	l_density = 100, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, 
	ld_fraction = 0.95, pruning_method = "MIG++", windows = NULL)
}
\arguments{
	\item{phase_files}{
		The list of names of the input files with phased genotypes.
		All input files must be in the same format: VCF or HAPMAP2.
	}
	\item{output_files}{
		The list of names of the output files where to store the haplotype blocks.
		One output file for every input file.
	}
	\item{processes}{
		An integer >= 1, which indicates the number of parallel processes. 
		All processes are created on the \emph{localhost} and communicate through sockets.
	}
	\item{phase_file_format}{
		Format of the phase_files: VCF (default) or HAPMAP2.
	}
	\item{map_files}{
		The list of names of the map files with base-pair positions of each SNP.
		One map file for every input file.
		The order of map files must correspond to the order of input files.
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
		Numeric vector where every value corresponds to the according input file (chromosome) and specifies the number of SNPs within the window in MIG++ search space pruning method.
		If NULL (default), all values are calculated on the fly based on the corresponding chromosome lengths and ld_fraction.
	}
}
\references{
	Zapata, C., Alvarez, G., Carollo, C. (1997) Approximate variance of the standardized measure of gametic disequilibrium D'. \emph{American Journal of Human Genetics}, \bold{61}(3), 771--774.

	Gabriel, S. B. et al. (2002) The Structure of Haplotype Blocks in the Human Genome. \emph{Science}, \bold{296}(5576), 2225--2229.
	
	Wall, J. D. and Pritchard, J. K. (2003) Assessing the performance of the haplotype block model of linkage disequilibrium. \emph{American Journal of Human Genetics}, \bold{73}(3), 502--515.
	
	Barrett, J. C. et al. (2005) Haploview: analysis and visualization of LD and haplotype maps. \emph{Bioinformatics}, \bold{21}(2), 263--265.
}
\author{Daniel Taliun, Johann Gamper, Cristian Pattaro}
\seealso{
	See \link{mig} for the haplotype block definition, description of the D' distribution modeling and pruning methods.
}
\keyword{misc}
\keyword{utilities}
\keyword{package}