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

mig_multi_regions <- function(phase_file, output_files, regions_start, regions_end, processes = 1, phase_file_format = "VCF", map_file = NULL, maf = 0.0, ci_method = "WP", l_density = 100, ld_ci = c(0.7, 0.98), ehr_ci = 0.9, ld_fraction = 0.95, pruning_method = "MIG++", windows = NULL) {
	if (missing(phase_file)) {
		stop("The 'phase_file' argument is missing.");
	}
	
	if (missing(output_files)) {
		stop("The 'output_files' argument is missing.");
	}
	
	if (missing(regions_start)) {
		stop("The 'regions_start' argument is missing.");
	}
	
	if (missing(regions_end)) {
		stop("The 'regions_end' argument is missing.");
	}
	
	result <- .Call("mig_multi_regions", phase_file, output_files, regions_start, regions_end, processes, phase_file_format, map_file, maf, ci_method, l_density, ld_ci, ehr_ci, ld_fraction, pruning_method, windows)
}