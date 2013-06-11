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

create_browser_track <- function(input_file, output_file, chromosome, strand, track_name = "", track_desc = "") {
	if (missing(input_file)) {
		stop("The 'phase_file' argument is missing.");
	}
	
	if (missing(output_file)) {
		stop("The 'output_file' argument is missing.");
	}
	
	if (missing(chromosome)) {
		stop("The 'chromosome' argument is missing.");
	}
	
	if (missing(strand)) {
		stop("The 'strand' argument is missing.");
	}
	
	if (is.character(input_file)) {
		if (length(input_file) <= 0) {
			stop("The input file name is empty.")
		} else if (length(input_file) > 1) {
			stop("The input file name has multiple values.")
		}
		input_file <- gsub("^\\s+|\\s+$", "", input_file)
		if (nchar(input_file) <= 0) {
			stop("The input file name must be a non-blank character string.")
		}
	} else {
		stop("The input file name must be a character string.")
	}
	
	if (is.character(output_file)) {
		if (length(output_file) <= 0) {
			stop("The output file name is empty.")
		} else if (length(output_file) > 1) {
			stop("The output file name has multiple values.")
		}
		output_file <- gsub("^\\s+|\\s+$", "", output_file)
		if (nchar(output_file) <= 0) {
			stop("The output file name must be a non-blank character string.")
		}
	} else {
		stop("The output file name must be a character string.")
	}
	
	if (is.character(chromosome)) {
		if (length(chromosome) <= 0) {
			stop("The chromosome name is empty.")
		} else if (length(chromosome) > 1) {
			stop("The chromosome name has multiple values.")
		}
		chromosome <- gsub("^\\s+|\\s+$", "", chromosome)
		if (nchar(chromosome) <= 0) {
			stop("The chromosome name must be a non-blank character string.")
		}
		if (!grepl("^chr", chromosome)[1]) {
			stop("The chromosome name must start with \"chr\" prefix.")
		}
	} else {
		stop("The chromosome name must be a character string.")
	}
	
	if (is.character(strand)) {
		if (length(strand) <= 0) {
			stop("The strand is empty.")
		} else if (length(strand) > 1) {
			stop("The strand has multiple values.")
		}
		strand <- gsub("^\\s+|\\s+$", "", strand)
		if (nchar(strand) <= 0) {
			stop("The strand a non-blank character.")
		}
		if ((strand != "+") && (strand != "-")) {
			stop("The strand must be equal to \"+\" or \"-\".")
		}
	} else {
		stop("The strand must be a character.")
	}
	
	t <- read.table(input_file, header=T, stringsAsFactors=F, sep="\t")
	
	definition <- "track name="
	definition <- paste(definition, "\"", track_name,"\" description=", sep="")
	definition <- paste(definition, "\"", track_desc,"\" itemRgb=\"on\"", sep="")

	colors_function <- colorRampPalette(c("grey80", "red"), space = "rgb")
	colors <- col2rgb(colors_function(10))
	colors_indices <- ceiling(t$HAPS_DIVERSITY * 10)
	
	t$chrom <- chromosome
	t$chromStart <- t$START_BP
	t$chromEnd <- t$END_BP
	t$name <- t$BLOCK_NAME
	t$score <- 0
	t$strand <- strand
	t$thickStart <- 0
	t$thickEnd <- 0
	t$itemRgb <- paste(colors[1, colors_indices], ",", colors[2, colors_indices], ",", colors[3, colors_indices], sep="")
	
	t <- t[c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")]
	
	connection = file(output_file, open="wt")
	writeLines(definition, connection)
	close(connection)
	
	connection = file(output_file, open="at")
	write.table(t, connection, col.names=F, row.names=F, quote=F, sep="\t")
	close(connection)
	
}