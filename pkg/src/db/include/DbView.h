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

#ifndef DBVIEW_H_
#define DBVIEW_H_

#include <stdlib.h>

using namespace std;

class DbView {
private:
	DbView(double maf_threshold, unsigned long int start_position, unsigned long int end_position);

public:
	const char* hap_file_name;
	const char* map_file_name;

	double maf_threshold;
	unsigned long int start_position;
	unsigned long int end_position;

	unsigned int n_unfiltered_markers;

	unsigned int n_haplotypes;
	unsigned int n_markers;
	char** markers;
	unsigned long int* positions;
	char* major_alleles;
	char* minor_alleles;
	double* major_allele_freqs;
	char** haplotypes;

	virtual ~DbView();

	double get_memory_usage();

	friend class Db;
};

#endif
