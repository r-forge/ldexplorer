\name{mig}
\alias{mig}
\title{Memory-efficient Implementation of Gabriel et al. 2002 haplotype block definition}
\description{
	Function for the efficient whole-genome haplotype block partitioning.
	Haplotype blocks are defined based on D' coefficient of linkage disequilibrium (Gabriel et al., 2002).
}
\usage{mig(file_format, phase_file, legend_file, output_file, region, maf, ld_lower_ci, ld_upper_ci, ehr_upper_ci, ld_fraction, prune, window)}
\arguments{
	\item{file_format}{
		Format of input file(s) with phased genotypes. 
		Supported formats are VCF and HAPMAP2.
	}
	\item{phase_file}{
		Name of an input file with phased genotypes in VCF or HAPMAP2 format.
	}
	\item{legend_file}{
		Name of an input file with base pair positions of markers.
		Mandatory only when HAPMAP2 format is specified.
	}
	\item{output_file}{
		Name of an output file with haplotype blocks.
	}
	\item{region}{
		A chromosomal region for haplotype block partitioning.
		Is a vector of two values that correspond to the start and end position in base pairs.
		If NA (default), then the whole chromosome is partitioned.
	}
	\item{maf}{
		A threshold for the minor allele frequency of SNPs.
		It varies from 0 (default) to 0.5.
		SNPs with lower (<=) minor allele frequency are not considered.
	}
	\item{ld_lower_ci}{
		A threshold for the lower bound of D' confidence interval between a SNP pair that is in strong LD.
	}
	\item{ld_upper_ci}{
		A threshold for the lower upper of D' confidence interval between a SNP pair that is in strong LD.
	}
	\item{ehr_upper_ci}{
		A threshold for the lower upper of D' confidence interval between a SNP pair that shows strong evidence of historical recombination.
	}
	\item{ld_fraction}{
		A threshold for the fraction of SNP pairs in string LD within a haplotype block.
	}
	\item{prune}{
		Name of a pruning method.
		Supported methods are MIG, MIG+ and MIG++ (default).
		MIG doesn't prune the search space.
		MIG++ performs the greatest pruning and significantly improves the runtime.
	}
	\item{window}{
		Size (number of SNPs) of the window when MIG++ pruning method is specified.
	}
}
\author{Daniel Taliun}
\keyword{misc}
\keyword{utilities}
\keyword{package}
\examples{}