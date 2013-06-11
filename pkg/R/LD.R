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

ld <- function(phase_file, snps_file, output_file, window = 500000, coefficient = "dprime", maf = 0.0, gzip = TRUE) {
	if (missing(phase_file)) {
		stop("The 'phase_file' argument is missing.");
	}
	
	if (missing(snps_file)) {
		stop("The 'snps_file' argument is missing.");
	}
	
	if (missing(output_file)) {
		stop("The 'output_file' argument is missing.");
	}
	
	result <- .Call("ld", phase_file, snps_file, output_file, window, coefficient, maf, gzip);
}