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

#include "include/Algorithm.h"

const char* Algorithm::ALGORITHM_MIG = "MIG";
const char* Algorithm::ALGORITHM_MIGP = "MIG+";
const char* Algorithm::ALGORITHM_MIGPP = "MIG++";

const unsigned int Algorithm::PRELIMINARY_BLOCKS_SIZE_INIT = 100000;
const unsigned int Algorithm::PRELIMINARY_BLOCKS_SIZE_INCREMENT = 10000;

const double Algorithm::EPSILON = 0.000000001;

Algorithm::Algorithm() throw (Exception) :
		db(NULL), ci_method(NULL),
		pos_strong_pair_cl(0.7), neg_strong_pair_cl(-0.7),
		pos_strong_pair_cu(0.98), neg_strong_pair_cu(-0.98),
		pos_recomb_pair_cu(0.9), neg_recomb_pair_cu(0.9),
		strong_pairs_fraction(0.95), strong_pair_weight(0.05), recomb_pair_weight(0.95),
		preliminary_blocks(NULL), n_preliminary_blocks(0u), preliminary_blocks_size(PRELIMINARY_BLOCKS_SIZE_INIT) {

	preliminary_blocks = (preliminary_block*)malloc(preliminary_blocks_size * sizeof(preliminary_block));
	if (preliminary_blocks == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}
}

Algorithm::~Algorithm() {
	db = NULL;

	free(preliminary_blocks);
	preliminary_blocks = NULL;
}

void Algorithm::set_dbview(const DbView* db) {
	this->db = db;
}

void Algorithm::set_ci_method(const char* ci_method) {
	this->ci_method = ci_method;
}

void Algorithm::set_likelihood_density(unsigned int likelihood_density) {
	this->likelihood_density = likelihood_density;
}

void Algorithm::set_strong_pair_cl(double ci_lower_bound) {
	pos_strong_pair_cl = ci_lower_bound;
	neg_strong_pair_cl = -pos_strong_pair_cl;
}

void Algorithm::set_strong_pair_cu(double ci_upper_bound) {
	pos_strong_pair_cu = ci_upper_bound;
	neg_strong_pair_cu = -pos_strong_pair_cu;
}

void Algorithm::set_recomb_pair_cu(double ci_upper_bound) {
	pos_recomb_pair_cu = ci_upper_bound;
	neg_recomb_pair_cu = -pos_recomb_pair_cu;
}

void Algorithm::set_strong_pairs_fraction(double fraction) {
	strong_pairs_fraction = fraction;
	strong_pair_weight = 1.0 - strong_pairs_fraction;
	recomb_pair_weight = strong_pairs_fraction;
}

void Algorithm::sort_preliminary_blocks() {
	qsort(preliminary_blocks, n_preliminary_blocks, sizeof(preliminary_block), preliminary_blocks_cmp);
}

Partition* Algorithm::get_block_partition() throw (Exception) {
	Partition* partition = NULL;

	bool* used_markers = NULL;

	unsigned int start = 0u;
	unsigned int end = 0u;

	partition = new Partition(db);

	partition->ci_method = ci_method;
	partition->likelihood_density = likelihood_density;
	partition->strong_pair_cl = pos_strong_pair_cl;
	partition->strong_pair_cu = pos_strong_pair_cu;
	partition->recomb_pair_cu = pos_recomb_pair_cu;
	partition->strong_pairs_fraction = strong_pairs_fraction;

	used_markers = (bool*)malloc(db->n_markers * sizeof(bool));
	if (used_markers == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < db->n_markers; ++i) {
		used_markers[i] = false;
	}

	for (unsigned int b = 0u; b < n_preliminary_blocks; ++b) {
		start = preliminary_blocks[b].start;
		end = preliminary_blocks[b].end;

		if (used_markers[start] || used_markers[end]) {
			continue;
		}

		partition->add_block(start, end);

		for (unsigned int i = start; i <= end; ++i) {
			used_markers[i] = true;
		}
	}

	free(used_markers);
	used_markers = NULL;

	return partition;
}

unsigned int Algorithm::get_n_preliminary_blocks() {
	return n_preliminary_blocks;
}

double Algorithm::get_memory_usage_preliminary_blocks() {
	return ((preliminary_blocks_size * sizeof(preliminary_block)) / 1048576.0);
}

double Algorithm::get_memory_usage() {
	return 0.0;
}

int Algorithm::preliminary_blocks_cmp(const void* first, const void* second) {
	preliminary_block* first_block = (preliminary_block*)first;
	preliminary_block* second_block = (preliminary_block*)second;

	if (first_block->length_bp > second_block->length_bp) {
		return -1;
	} else if (first_block->length_bp < second_block->length_bp) {
		return 1;
	} else {
		return first_block->start - second_block->start;
	}
}
