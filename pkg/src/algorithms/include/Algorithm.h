/*
 * Copyright © 2013 Daniel Taliun and Cristian Pattaro. All rights reserved.
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

#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include <map>
#include <math.h>

#include "../../auxiliary/include/auxiliary.h"
#include "../../db/include/Db.h"
#include "../../writer/include/WriterFactory.h"
#include "CIFactory.h"

using namespace std;

class Algorithm {
protected:
	static const unsigned int STRONG_PAIRS_SIZE_INIT;
	static const unsigned int STRONG_PAIRS_SIZE_INCREMENT;
	static const unsigned int BLOCKS_SIZE_INIT;
	static const unsigned int BLOCKS_SIZE_INCREMENT;

	struct pair {
		unsigned int first;
		unsigned int last;
		unsigned long int distance;
	};

	Db* db;
	pair* strong_pairs;
	unsigned int n_strong_pairs;
	unsigned int strong_pairs_size;
	unsigned int* blocks;
	unsigned int n_blocks;
	unsigned int blocks_size;

	bool is_compatible_haplotype(const char* first, const char* second);
	void get_block_diversity(unsigned int block_id,
			unsigned int* n_haps, unsigned int* n_unique_haps, unsigned int* n_common_haps, double* haps_diversity) throw (Exception);

	static int paircmp(const void* first, const void* second);

public:
	static const double EPSILON;

	Algorithm(Db& db) throw (Exception);
	virtual ~Algorithm();

	virtual void compute_preliminary_blocks() throw (Exception) = 0;

	void sort_preliminary_blocks();
	void select_final_blocks() throw (Exception);

	unsigned int get_n_strong_pairs();
	unsigned int get_n_blocks();

	void write_blocks(const char* output_file_name, const char* input_phase_file_name, const char* input_map_file_name,
			double maf_threshold, bool region, unsigned long int start, unsigned long int end, const char* ci_method) throw (Exception);

	double get_max_memory_usage();
};

#endif
