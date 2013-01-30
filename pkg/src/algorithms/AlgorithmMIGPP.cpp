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

#include "include/AlgorithmMIGPP.h"

AlgorithmMIGPP::AlgorithmMIGPP(Db& db) : Algorithm(db) {

}

AlgorithmMIGPP::~AlgorithmMIGPP() {

}

void AlgorithmMIGPP::compute_preliminary_blocks(const char* ci_method, unsigned int ci_precision, unsigned int window) throw (Exception) {
	CI* ci = NULL;

	long double* w_values = NULL;
	long double* w_values_sums = NULL;

	long double w_values_sum_left = 0.0;
	long double* w_values_sums_left = NULL;

	long double w_value_max = 0.0;
	long double* w_values_max = NULL;

	long int* breakpoints = NULL;
	long int* terminations = NULL;

	double lower_ci = 0.0;
	double upper_ci = 0.0;

	unsigned int current_window = 0u;

	long int breakpoint = 0;
	long int updated_breakpoint = 0;

	unsigned long int calculations = numeric_limits<unsigned long int>::max();;

	pair* new_strong_pairs = NULL;

	ci = CIFactory::create(*db, ci_method, ci_precision);

	w_values = (long double*)malloc(db->n_markers * sizeof(long double));
	if (w_values == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	w_values_sums = (long double*)malloc(db->n_markers * sizeof(long double));
	if (w_values_sums == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	w_values_sums_left = (long double*)malloc(db->n_markers * sizeof(long double));
	if (w_values_sums_left == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	w_values_max = (long double*)malloc(db->n_markers * sizeof(long double));
	if (w_values_max == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	terminations = (long int*)malloc(db->n_markers * sizeof(long int));
	if (terminations == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	breakpoints = (long int*)malloc(db->n_markers * sizeof(long int));
	if (breakpoints == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < db->n_markers; ++i) {
		w_values[i] = 0.0;
		w_values_sums[i] = 0.0;
		terminations[i] = i;
		breakpoints[i] = i;
	}

	for (unsigned int i = 0u; i < db->n_markers; ++i) {
		w_values_sum_left += strong_pair_weight * terminations[i];
		w_values_sums_left[i] = w_values_sum_left;
	}

	w_values_max[db->n_markers - 1u] = w_values_sums_left[db->n_markers - 1u];
	for (long int i = db->n_markers - 1u; i >= 2u; --i) {
		w_values_max[i - 1u] = w_values_max[i] > w_values_sums_left[i] ? w_values_max[i] : w_values_sums_left[i];
	}

	while (calculations > 0u) {
		current_window += window;

		calculations = 0u;

		breakpoint = 0;
		updated_breakpoint = 0;

		w_values_sum_left = 0.0;

		for (long int i = 1; i < db->n_markers; ++i) {
			if (updated_breakpoint == breakpoints[i]) {
				breakpoints[i] = breakpoint;
				breakpoint = updated_breakpoint = terminations[i];

				w_values_sum_left += strong_pair_weight * terminations[i] + w_values_sums[i];
				w_values_sums_left[i] = w_values_sum_left;

				continue;
			}

			if ((i - updated_breakpoint) > current_window) {
				breakpoints[i] = breakpoint = i - current_window;
			} else {
				breakpoints[i] = breakpoint;
				breakpoint = updated_breakpoint;
			}

			updated_breakpoint = terminations[i];

			for (long int j = terminations[i] - 1u; j >= breakpoint; --j) {
				++calculations;

				ci->get_CI(i, j, &lower_ci, &upper_ci);
				if (!isnan(lower_ci) && !isnan(upper_ci)) {
					if (((auxiliary::fcmp(lower_ci, pos_strong_pair_cl, EPSILON) >= 0) && (auxiliary::fcmp(upper_ci, pos_strong_pair_cu, EPSILON) >= 0)) ||
							((auxiliary::fcmp(lower_ci, neg_strong_pair_cu, EPSILON) <= 0) && (auxiliary::fcmp(upper_ci, neg_strong_pair_cl, EPSILON) <= 0))) {
						w_values_sums[i] += strong_pair_weight;
						w_values[j] += w_values_sums[i];
						if (auxiliary::fcmp(w_values[j], 0.0, EPSILON) >= 0) {
							if (n_strong_pairs >= strong_pairs_size) {
								strong_pairs_size += STRONG_PAIRS_SIZE_INCREMENT;
								new_strong_pairs = (pair*)realloc(strong_pairs, strong_pairs_size * sizeof(pair));
								if (new_strong_pairs == NULL) {
									delete ci;
									ci = NULL;

									free(w_values);
									w_values = NULL;

									free(w_values_sums);
									w_values_sums = NULL;

									free(w_values_sums_left);
									w_values_sums_left = NULL;

									free(w_values_max);
									w_values_max = NULL;

									free(terminations);
									terminations = NULL;

									free(breakpoints);
									breakpoints = NULL;

									throw Exception(__FILE__, __LINE__, "Error in memory reallocation.");
								}
								strong_pairs = new_strong_pairs;
								new_strong_pairs = NULL;
							}

							strong_pairs[n_strong_pairs].first = j;
							strong_pairs[n_strong_pairs].last = i;
							strong_pairs[n_strong_pairs].distance = db->positions[i] - db->positions[j];

							++n_strong_pairs;
						}
					} else if ((auxiliary::fcmp(lower_ci, neg_recomb_pair_cu, EPSILON) >= 0) && (auxiliary::fcmp(upper_ci, pos_recomb_pair_cu, EPSILON) <= 0)) {
						w_values_sums[i] -= recomb_pair_weight;
						w_values[j] += w_values_sums[i];
					} else {
						w_values[j] += w_values_sums[i];
					}
				} else {
					w_values[j] += w_values_sums[i];
				}

				/* With prior using pre-calculated sums, medium conservative. */
				w_value_max = w_values[j] + w_values_max[i] - w_values_sums_left[i];
				if (auxiliary::fcmp(w_value_max, 0.0, EPSILON) >= 0) {
					updated_breakpoint = j;
				}
			}

			terminations[i] = breakpoint;

			w_values_sum_left += strong_pair_weight * terminations[i] + w_values_sums[i];
			w_values_sums_left[i] = w_values_sum_left;
		}

		w_values_max[db->n_markers - 1u] = w_values_sums_left[db->n_markers - 1u];
		for (long int k = db->n_markers - 1u; k >= 2u; --k) {
			w_values_max[k - 1u] = w_values_max[k] > w_values_sums_left[k] ? w_values_max[k] : w_values_sums_left[k];
		}
	}

	delete ci;
	ci = NULL;

	free(w_values);
	w_values = NULL;

	free(w_values_sums);
	w_values_sums = NULL;

	free(w_values_sums_left);
	w_values_sums_left = NULL;

	free(w_values_max);
	w_values_max = NULL;

	free(terminations);
	terminations = NULL;

	free(breakpoints);
	breakpoints = NULL;
}

double AlgorithmMIGPP::get_memory_usage() {
	double memory_usage = 0.0;

	memory_usage += (4u * db->n_markers * sizeof(long double)) / 1048576.0;
	memory_usage += (2u * db->n_markers * sizeof(long int)) / 1048576.0;

	return memory_usage;
}
