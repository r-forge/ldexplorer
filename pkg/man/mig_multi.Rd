\name{mig_multi}
\alias{mig_multi}
\title{Memory-efficient Implementation of Gabriel et al. 2002 haplotype block definition}
\description{
	Function for the efficient whole-genome haplotype block partitioning.
	Haplotype blocks are defined based on D' coefficient of linkage disequilibrium (Gabriel et al., 2002).
}
\usage{mig_multi(phase_files, output_files, phase_file_format = "VCF", map_files = NULL, maf = 0.0, ci_method = "WP", l_density = 100, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, ld_fraction = 0.95, pruning_method = "MIG++", windows = NULL)}
\arguments{
	\item{phase_files}{
		Name of the input file with phased genotypes in VCF or HAPMAP2 format.
	}
	\item{output_files}{
		Name of the output file where to store the haplotype blocks.
	}
	\item{phase_file_format}{
		Format of the phase_file: VCF (default) or HAPMAP2.
	}
	\item{map_files}{
		Name of the map file with base-pair position for each SNP.
		Mandatory when file_formatq = HAPMAP2.
	}
	\item{maf}{
		Minor Allele Frequency (MAF) threshold: SNPs with MAF <= maf will not be considered.
		The threshold may vary from 0 (default) to 0.5.
	}
	\item{ci_method}{
		Confidence interval (CI) estimation method.
		Supported methods are WP (default) = Wall and Pitchard (2003) method; AV = approximate variance estimator by Zapata et al (1997).
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
		Size of the window in MIG++ search space pruning method.
		During an iterative processing of a chromosome, the number of considered SNP pairs at every iteration is restricted by the window.
		For every SNP, only its correlation with SNPs, which are within the window, are computed.
		By default, the window is calculated on the fly based on the chromosome length and ld_fraction threshold.
	}
}
\references{
	Gabriel, Stacey B. et al. (2002) \emph{The Structure of Haplotype Blocks in the Human Genome}. Science, 296(5576), 2225-2229.
	
	J.C. Barrett et al. (2005) \emph{Haploview: analysis and visualization of LD and haplotype maps}. Bioinformatics, 21(2), 263-265.
}
\author{Daniel Taliun, Johann Gamper, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{package}
\examples{}