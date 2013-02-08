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

mig_multi <- function(phase_files, output_files, phase_file_format = "VCF", map_files = NULL, maf = 0.0, ci_method = "WP", l_density = 100, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, ld_fraction = 0.95, pruning_method = "MIG++", windows = NULL) {
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
			
		mig(phase_file, output_file, phase_file_format, map_file, NULL, maf, ci_method, l_density, ld_ci, ehr_ci, ld_fraction, pruning_method, window)
	}
}