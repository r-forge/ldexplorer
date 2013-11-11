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

mig_rsq <- function(phase_file, output_file, phase_file_format = "VCF", map_file = NULL, region = NULL, maf = 0.0, weak_rsq = 0.5, strong_rsq = 0.8, fraction = 0.95, pruning_method = "MIG++", window = NULL) {
	if (missing(phase_file)) {
		stop("The 'phase_file' argument is missing.");
	}
	
	if (missing(output_file)) {
		stop("The 'output_file' argument is missing.");
	}
	
	result <- .Call("mig_rsq", phase_file, output_file, phase_file_format, map_file, region, maf, weak_rsq, strong_rsq, fraction, pruning_method, window)
}