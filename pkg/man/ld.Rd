\name{ld}
\alias{ld}
\title{Linkage disequilibrium (LD) computation from phased data}
\description{
	Function to compute linkage disequilibrium (LD) from phased data.
	The LD is computed between a SNP of interest and all other SNPs within a specified window around it.
}
\usage{
	ld(phase_file, snps_file, output_file, window = 500000, coefficient = "dprime", maf = 0.0, gzip = TRUE)
}
\arguments{
	\item{phase_file}{
		The name of the input file with phased genotypes in the VCF format.
	}
	\item{snps_file}{
		The name of the input file with the list of SNP identifiers (one per line).
		The SNP identifier can be an rsID, chromosomal position in base-pairs or any other identifier format that is present in phase_file.
		For every SNP in the list, the LD between it and all other SNPs within a specified window around it is computed.
	}
	\item{output_file}{
		The name of the output file.
		Output file includes five tab-separated columns.
		The first four columns are: FIRST_MARKER, FIRST_BP, SECOND_MARKER, SECOND_BP.
		The last column depends on the specified coefficient and can be D, DPRIME or R2, accordingly.
	}
	\item{window}{
		The number of base pairs to consider around every SNP of interest.
	}
	\item{coefficient}{ 
		The LD coefficient to be calculated between a pair of SNPs.
		The supported LD coefficients are "d" (D), "dprime" (D') or "r2" (r^2). 
	}
	\item{maf}{
		Minor Allele Frequency (MAF) threshold: SNPs with MAF <= maf will not be considered.
		The threshold may vary from 0 (default) to 0.5.
	}
	\item{gzip}{
		TRUE if output file is in gzip format.
	}
}
\author{Daniel Taliun, Johann Gamper, Cristian Pattaro}
\keyword{misc}
\keyword{utilities}
\keyword{package}