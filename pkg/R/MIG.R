#
# Copyright © 2013 Daniel Taliun and Cristian Pattaro. All rights reserved.
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

mig <- function(phase_file, output_file, file_format = "VCF", legend_file = NULL, region = NULL, maf = 0.0, ci_method = "AV", ci_precision = 1000, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, ld_fraction = 0.95, pruning_method = "MIG++", window = NULL) {
	if (missing(phase_file)) {
		stop("The 'phase_file' argument is missing.");
	}
	
	if (missing(output_file)) {
		stop("The 'output_file' argument is missing.");
	}
	
	result <- .Call("mig", phase_file, output_file, file_format, legend_file, region, maf, ci_method, ci_precision, ld_ci, ehr_ci, ld_fraction, pruning_method, window)
}