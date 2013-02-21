#
# Copyright © 2013 Daniel Taliun, Johann Gamper and Cristian Pattaro. All rights reserved.
#
# This file is part of LDExplorer.
#
# LDExplorer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LDExplorer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LDExplorer.  If not, see <http://www.gnu.org/licenses/>.
#

mig_multi_chr <- function(phase_files, output_files, processes = 1, phase_file_format = "VCF", map_files = NULL, maf = 0.0, ci_method = "WP", l_density = 100, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, ld_fraction = 0.95, pruning_method = "MIG++", windows = NULL) {
	if (missing(phase_files)) {
		stop("The 'phase_files' argument is missing.");
	}
	
	if (is.character(phase_files)) {
		if (length(phase_files) <= 0) {
			stop("The 'phase_files' argument must be a non-empty character vector.")
		} 
	} else {
		stop("The 'phase_files' argument must be a character vector.")
	}
	
	if (missing(output_files)) {
		stop("The 'output_files' argument is missing.");
	}
	
	if (is.character(output_files)) {
		if (length(phase_files) != length(output_files)) {
			stop("The 'phase_files' and 'output_files' arguments must be character vectors of the same length.")
		}
	} else {
		stop("The 'output_files' argument must be a character vector.")
	}
	
	if ((toupper(phase_file_format) == "HAPMAP2")) {
		if (!is.null(map_files)) {
			if (is.character(map_files)) {
				if (length(phase_files) != length(map_files)) {
					stop("The 'phase_files' and 'map_files' arguments must be character vectors of the same length.")
				}
			} else {
				stop("The 'map_files' argument must be a character vector.")
			}
		} else {
			stop("The 'map_files' argument is not specified.");
		}
	}
	
	if (toupper(pruning_method) == "MIG++") {
		if (!is.null(windows)) {
			if (is.numeric(windows)) {
				if (length(phase_files) != length(windows)) {
					stop("The 'phase_files' and 'windows' arguments must be vectors of the same length.")
				}
			} else {
				stop("The 'windows' argument must be a numeric vector.")
			}
		}
	}
	
	if (is.numeric(processes)) {
		if (abs(processes - round(processes)) < .Machine$double.eps ^ 0.5) {
			processes <- as.integer(processes)
			if (is.na(processes)) {
				stop("Error while casting number of processes argument to integer.")
			} else if (processes < 1) {
				stop("The number of processes must be greater than 0.")	
			}
		} else {
			stop("The number of processes argument must be integer.")
		}
	} else {
		stop("The number of processes argument must be numeric.")
	}
	
	args_per_cluster <- list()
	
	for (i in 1:length(phase_files)) {
		phase_file <- phase_files[i]
		output_file <- output_files[i]
		
		if (!is.null(map_files)) {
			map_file <- map_files[i]
		} else {
			map_file <- NULL
		}
		
		if (!is.null(windows)) {
			window <- windows[i]
		} else {
			window <- NULL
		}
		
		args <- list(phase_file = phase_file, output_file = output_file, phase_file_format = phase_file_format, map_file = map_file, region = NULL, maf = maf, ci_method = ci_method, l_density = l_density, ld_ci = ld_ci, ehr_ci = ehr_ci, ld_fraction = ld_fraction, pruning_method = pruning_method, window = window)
		
		args_per_cluster <- c(args_per_cluster, list(args))
	}
	
	if (processes == 1) {
		for (i in 1:length(args_per_cluster)) {
			x <- args_per_cluster[[i]]
			mig(x$phase_file, x$output_file, x$phase_file_format, x$map_file, x$region, x$maf, x$ci_method, x$l_density, x$ld_ci, x$ehr_ci, x$ld_fraction, x$pruning_method, x$window)
		}
	} else {
		cat("Initializing cluster... ")
		start_time <- proc.time()
		
		require(parallel, quietly = TRUE)
	
		clusters <- makeCluster(rep("localhost", processes), type="PSOCK")
		ldexplorer_package_path <- dirname(.path.package("LDExplorer"))
		clusterExport(clusters, "ldexplorer_package_path", envir=environment())
		clusterEvalQ(clusters, .libPaths(union(ldexplorer_package_path, .libPaths())))
		clusterEvalQ(clusters, library(LDExplorer))
		
		elapsed_time <- proc.time() - start_time
		cat("Done (", elapsed_time[3], " sec).\n", sep="")
		
		cat("Processing... ")
		start_time <- proc.time()
		
		cluster_result <- clusterApply(clusters, args_per_cluster, function(x) {
					sink(paste(x$output_file, ".log", sep=""), split=F)
					mig(x$phase_file, x$output_file, x$phase_file_format, x$map_file, x$region, x$maf, x$ci_method, x$l_density, x$ld_ci, x$ehr_ci, x$ld_fraction, x$pruning_method, x$window)
				})
		
		elapsed_time <- proc.time() - start_time
		cat("Done (", elapsed_time[3], " sec).\n", sep="")
		
		cat("Stopping cluster... ")
		start_time <- proc.time()
		
		stopCluster(clusters)
		
		elapsed_time <- proc.time() - start_time
		cat("Done (", elapsed_time[3], " sec).\n", sep="")
	}
}