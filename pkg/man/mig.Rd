\name{mig}
\alias{mig}
\title{Memory-efficient Implementation of Gabriel et al. 2002 haplotype block definition}
\description{
	Function for the efficient whole-genome haplotype block partitioning.
	Haplotype blocks are defined based on D' coefficient of linkage disequilibrium (Gabriel et al., 2002).
}
\usage{mig(phase_file, output_file, file_format = "VCF", legend_file = NULL, region = NULL, maf = 0.0, ci_method = "AV", ci_precision = 1000, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, ld_fraction = 0.95, pruning_method = "MIG++", window = NULL)}
\arguments{
	\item{phase_file}{
		Name of an input file with phased genotypes in VCF or HAPMAP2 format.
	}
	\item{output_file}{
		Name of an output file with haplotype blocks.
	}
	\item{file_format}{
		Format of input file(s) with phased genotypes. 
		Supported formats are VCF (default) and HAPMAP2.
	}
	\item{legend_file}{
		Name of an input file with base pair positions of markers.
		Mandatory only when HAPMAP2 format is specified.
	}
	\item{region}{
		A chromosomal region to be partitioned.
		Is a vector of two values that correspond to the start and end position in base pairs.
		If NULL (default), then the whole chromosome is partitioned.
	}
	\item{maf}{
		A threshold for the minor allele frequency of SNPs.
		It varies from 0 (default) to 0.5.
		SNPs with lower (<=) minor allele frequency than the specified threshold are not considered.
	}
	\item{ci_method}{
		Name of a method to compute the confidence interval (CI) of D' between a SNP pair.
		Supported methods are WP and AV (default).
		The WP method computes CI of |D'| and is based on the work of Wall and Pitchard (2003) and the Haploview (Barrett et al., 2005) implementation.
		The AV method computes CI of D' using the approximate variance estimator of D' (Zapata, 1997).
	}
	\item{ci_precision}{
		Number of likelihood estimation points in WP method to compute confidence interval (CI) of |D'|.
		The default is 1000, the minimum recommended is 100.
		Higher number of likelihood estimation points improves the precision of CI values, but degrades the runtime.
	}
	\item{ld_ci}{
		Thresholds for the lower bound (CL) and upper bound (CU) of the 90\% confidence interval (CI) of D'.
		Must be specified as a vector of two values, where the first value corresponds to the CL and the second value corresponds to the CU.
		If WP method to compute CI of |D'| was specified, then a SNP pair is in strong LD if CL >= ld_lower_ci and CU >= ld_upper_ci.
		If AV method to compute CI of D' was specified, then a SNP pair is in strong LD if CL >= ld_lower_ci and CU >= ld_upper_ci or CL <= -ld_upper_ci and CU <= ld_lower_ci.
		By default, the threshold for the CL is 0.7 and for the CU is 0.98 (Gabriel et al., 2002).
	}
	\item{ehr_ci}{
		A threshold for the upper bound (CU) of the 90\% confidence interval (CI) of D'.
		If WP method to compute CI of |D'| was specified, then a SNP pair shows strong evidence of historical recombination (EHR) if CU <= ehr_upper_ci. 
		If AV method to compute CI of D' was specified, then a SNP pair shows strong EHR if CL >= ehr_upper_ci and CU <= ehr_upper_ci.
		By default, the threshold for the CU is 0.9 (Gabriel et al., 2002).
	}
	\item{ld_fraction}{
		A threshold for the fraction of strong LD SNP pairs within a haplotype block.
		If the ratio of strong LD SNP pairs to all informative (strong LD and strong EHR) SNP pairs is greater (>=) than ld_fraction, then a chromosomal region is a haplotype block.
		By default, the threshold is 0.95 (Gabriel et al., 2002).
	}
	\item{pruning_method}{
		Name of a search space pruning method.
		Supported methods are MIG, MIG+ and MIG++ (default).
		The MIG doesn't prune the search space and calculates correlations between all possible SNP pairs.
		The MIG+ avoids the calculation of correlations between SNP pairs that don't influence the final haplotype block partitioning.
		The MIG++ is an improved version of MIG+ that performs more precise pruning and significantly reduces the runtime.
		All pruning methods produce identical haplotype block partitioning.
	}
	\item{window}{
		Size of the window in MIG++ search space pruning method.
		During an iterative processing of a chromosome, the number of considered SNP pairs at every iteration is restricted by the window.
		For every SNP, only its correlation with SNPs, which are within the window, are computed.
		By default, the window is calculated on the fly based on the chromosome length and ld_fraction threshold.
	}
}
\author{Daniel Taliun, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{package}
\examples{}