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

#ifndef PARTITION_H_
#define PARTITION_H_

#include <map>
#include <vector>

#include "../../LDExplorer.h"
#include "../../writer/include/WriterFactory.h"
#include "../../db/include/DbView.h"

using namespace std;

class Partition {
private:
	static const unsigned int BLOCKS_SIZE_INIT;
	static const unsigned int BLOCKS_SIZE_INCREMENT;

	struct block {
		unsigned int start;
		unsigned int end;
	};

	const DbView* db;

	block* blocks;
	unsigned int n_blocks;
	unsigned int blocks_size;

	bool is_compatible_haplotype(const char* first, const char* second);
	void get_block_diversity(unsigned int block_id,
			unsigned int* n_haps, unsigned int* n_unique_haps, unsigned int* n_common_haps, double* haps_diversity) throw (Exception);

public:
	bool rsq_blocks;
	const char* ci_method;
	unsigned int likelihood_density;
	double strong_pair_cl;
	double strong_pair_cu;
	double recomb_pair_cu;
	double strong_pair_rsq;
	double strong_pairs_fraction;
	const char* pruning_method;
	unsigned int window;

	Partition(const DbView* db) throw (Exception);
	virtual ~Partition();

	void add_block(unsigned int start, unsigned int end) throw (Exception);
	unsigned int get_n_blocks();

	void write(const char* output_file_name) throw (Exception);

	double get_memory_usage();
};

#endif
