/*
 * Copyright © 2013 Daniel Taliun, Johann Gamper and Cristian Pattaro. All rights reserved.
 *
 * This file is part of LDExplorer.
 *
 * LDExplorer is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LDExplorer is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LDExplorer.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "include/DbView.h"

DbView::DbView(double maf_threshold, unsigned long int start_position, unsigned long int end_position) :
	maf_threshold(maf_threshold), start_position(start_position), end_position(end_position),
	n_unfiltered_markers(0u), n_haplotypes(0u), n_markers(0u), markers(NULL), positions(NULL),
	major_alleles(NULL), minor_alleles(NULL), major_allele_freqs(NULL), haplotypes(NULL) {

}

DbView::~DbView() {
	if (markers != NULL) {
		free(markers);
		markers = NULL;
	}

	if (positions != NULL) {
		free(positions);
		positions = NULL;
	}

	if (major_alleles != NULL) {
		free(major_alleles);
		major_alleles = NULL;
	}

	if (minor_alleles != NULL) {
		free(minor_alleles);
		minor_alleles = NULL;
	}

	if (major_allele_freqs != NULL) {
		free(major_allele_freqs);
		major_allele_freqs = NULL;
	}

	if (haplotypes != NULL) {
		free(haplotypes);
		haplotypes = NULL;
	}
}

double DbView::get_memory_usage() {
	double memory_usage = 0.0;

	if (markers != NULL) {
		memory_usage += (n_markers * sizeof(char*)) / 1048576.0;
	}

	if (positions != NULL) {
		memory_usage += (n_markers * sizeof(unsigned long int)) / 1048576.0;
	}

	if (major_alleles != NULL) {
		memory_usage += (n_markers * sizeof(char)) / 1048576.0;
	}

	if (minor_alleles != NULL) {
		memory_usage += (n_markers * sizeof(char)) / 1048576.0;
	}

	if (major_allele_freqs != NULL) {
		memory_usage += (n_markers * sizeof(double)) / 1048576.0;
	}

	if (haplotypes != NULL) {
		memory_usage += (n_markers * sizeof(char*)) / 1048576.0;
	}

	return memory_usage;
}
