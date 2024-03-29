/*
 * Copyright � 2013 Daniel Taliun, Johann Gamper and Cristian Pattaro. All rights reserved.
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

#include "include/AlgorithmMIG.h"

AlgorithmMIG::AlgorithmMIG() : Algorithm() {

}

AlgorithmMIG::~AlgorithmMIG() {

}

void AlgorithmMIG::compute_preliminary_blocks() throw (Exception) {
	CI* ci = NULL;

	long double* w_values = NULL;
	long double w_values_sum = 0.0;

	double lower_ci = 0.0;
	double upper_ci = 0.0;

	preliminary_block* new_strong_pairs = NULL;

	n_preliminary_blocks = 0u;
	rsq_preliminary_blocks = false;

	ci = CIFactory::create(ci_method, likelihood_density);
	ci->set_dbview(db);

	w_values = (long double*)malloc(db->n_markers * sizeof(long double));
	if (w_values == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < db->n_markers; ++i) {
		w_values[i] = 0.0;
	}

	for (unsigned int i = 1u; i < db->n_markers; ++i) {
		w_values_sum = 0.0;
		for (long int j = i - 1u; j >= 0; --j) {
			ci->get_CI(i, j, &lower_ci, &upper_ci);
			if (!isnan(lower_ci) && !isnan(upper_ci)) {
				if (((auxiliary::fcmp(lower_ci, pos_strong_pair_cl, EPSILON) >= 0) && (auxiliary::fcmp(upper_ci, pos_strong_pair_cu, EPSILON) >= 0)) ||
						((auxiliary::fcmp(lower_ci, neg_strong_pair_cu, EPSILON) <= 0) && (auxiliary::fcmp(upper_ci, neg_strong_pair_cl, EPSILON) <= 0))) {
					w_values_sum += strong_pair_weight;
					w_values[j] += w_values_sum;
					if (auxiliary::fcmp(w_values[j], 0.0, EPSILON) >= 0) {
						if (n_preliminary_blocks >= preliminary_blocks_size) {
							preliminary_blocks_size += PRELIMINARY_BLOCKS_SIZE_INCREMENT;
							new_strong_pairs = (preliminary_block*)realloc(preliminary_blocks, preliminary_blocks_size * sizeof(preliminary_block));
							if (new_strong_pairs == NULL) {
								delete ci;
								ci = NULL;

								free(w_values);
								w_values = NULL;

								throw Exception(__FILE__, __LINE__, "Error in memory reallocation.");
							}
							preliminary_blocks = new_strong_pairs;
							new_strong_pairs = NULL;
						}

						preliminary_blocks[n_preliminary_blocks].start = j;
						preliminary_blocks[n_preliminary_blocks].end = i;
						preliminary_blocks[n_preliminary_blocks].length_bp = db->positions[i] - db->positions[j];

						++n_preliminary_blocks;
					}
				} else if ((auxiliary::fcmp(lower_ci, neg_recomb_pair_cu, EPSILON) >= 0) && (auxiliary::fcmp(upper_ci, pos_recomb_pair_cu, EPSILON) <= 0)) {
					w_values_sum -= recomb_pair_weight;
					w_values[j] += w_values_sum;
				} else {
					w_values[j] += w_values_sum;
				}
			} else {
				w_values[j] += w_values_sum;
			}
		}
	}

	delete ci;
	ci = NULL;

	free(w_values);
	w_values = NULL;
}

void AlgorithmMIG::compute_preliminary_blocks_rsq() throw (Exception) {
	CI* ci = NULL;

	long double* w_values = NULL;
	long double w_values_sum = 0.0;

	double rsq = 0.0;

	preliminary_block* new_strong_pairs = NULL;

	n_preliminary_blocks = 0u;
	rsq_preliminary_blocks = true;

	ci = CIFactory::create(CI::NONE, likelihood_density);
	ci->set_dbview(db);

	w_values = (long double*)malloc(db->n_markers * sizeof(long double));
	if (w_values == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < db->n_markers; ++i) {
		w_values[i] = 0.0;
	}

	for (unsigned int i = 1u; i < db->n_markers; ++i) {
		w_values_sum = 0.0;
		for (long int j = i - 1u; j >= 0; --j) {
			rsq = ci->get_rsq(i, j);
			if (!isnan(rsq)) {
				if (auxiliary::fcmp(rsq, strong_pair_rsq, EPSILON) >= 0) {
					w_values_sum += strong_pair_weight;
					w_values[j] += w_values_sum;
					if (auxiliary::fcmp(w_values[j], 0.0, EPSILON) >= 0) {
						if (n_preliminary_blocks >= preliminary_blocks_size) {
							preliminary_blocks_size += PRELIMINARY_BLOCKS_SIZE_INCREMENT;
							new_strong_pairs = (preliminary_block*)realloc(preliminary_blocks, preliminary_blocks_size * sizeof(preliminary_block));
							if (new_strong_pairs == NULL) {
								delete ci;
								ci = NULL;

								free(w_values);
								w_values = NULL;

								throw Exception(__FILE__, __LINE__, "Error in memory reallocation.");
							}
							preliminary_blocks = new_strong_pairs;
							new_strong_pairs = NULL;
						}

						preliminary_blocks[n_preliminary_blocks].start = j;
						preliminary_blocks[n_preliminary_blocks].end = i;
						preliminary_blocks[n_preliminary_blocks].length_bp = db->positions[i] - db->positions[j];

						++n_preliminary_blocks;
					}
				} else if (auxiliary::fcmp(rsq, weak_pair_rsq, EPSILON) < 0) {
					w_values_sum -= recomb_pair_weight;
					w_values[j] += w_values_sum;
				} else {
					w_values[j] += w_values_sum;
				}
			} else {
				w_values[j] += w_values_sum;
			}
		}
	}

	delete ci;
	ci = NULL;

	free(w_values);
	w_values = NULL;
}

Partition* AlgorithmMIG::get_block_partition() throw (Exception) {
	Partition* partition = Algorithm::get_block_partition();

	partition->pruning_method = Algorithm::ALGORITHM_MIG;
	partition->window = 0u;

	return partition;
}


double AlgorithmMIG::get_memory_usage() {
	return ((db->n_markers * sizeof(long double)) / 1048576.0);
}
