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

#include "include/AlgorithmMIGP.h"

AlgorithmMIGP::AlgorithmMIGP(Db& db) : Algorithm(db) {

}

AlgorithmMIGP::~AlgorithmMIGP() {

}

void AlgorithmMIGP::compute_preliminary_blocks(const char* ci_method, unsigned int ci_precision, unsigned int window) throw (Exception) {
	CI* ci = NULL;

	long double* w_values = NULL;
	long double w_values_sum = 0.0;

	long double w_value_max = 0.0;

	double lower_ci = 0.0;
	double upper_ci = 0.0;

	long int breakpoint = 0;
	long int updated_breakpoint = 0;

	pair* new_strong_pairs = NULL;

	ci = CIFactory::create(*db, ci_method, ci_precision);

	w_values = (long double*)malloc(db->n_markers * sizeof(long double));
	if (w_values == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < db->n_markers; ++i) {
		w_values[i] = 0.0;
	}

	for (unsigned int i = 1u; i < db->n_markers; ++i) {
		w_values_sum = 0.0;
		breakpoint = updated_breakpoint;
		updated_breakpoint = i;
		for (long int j = i - 1u; j >= breakpoint; --j) {
			ci->get_CI(i, j, &lower_ci, &upper_ci);
			if (!isnan(lower_ci) && !isnan(upper_ci)) {
				if (((auxiliary::fcmp(lower_ci, 0.7, EPSILON) >= 0) && (auxiliary::fcmp(upper_ci, 0.98, EPSILON) >= 0)) ||
						((auxiliary::fcmp(lower_ci, -0.98, EPSILON) <= 0) && (auxiliary::fcmp(upper_ci, -0.7, EPSILON) <= 0))) {
					w_values_sum += 0.05;
					w_values[j] += w_values_sum;
					if (auxiliary::fcmp(w_values[j], 0.0, EPSILON) >= 0) {
						if (n_strong_pairs >= strong_pairs_size) {
							strong_pairs_size += STRONG_PAIRS_SIZE_INCREMENT;
							new_strong_pairs = (pair*)realloc(strong_pairs, strong_pairs_size * sizeof(pair));
							if (new_strong_pairs == NULL) {
								delete ci;
								ci = NULL;

								free(w_values);
								w_values = NULL;

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
				} else if ((auxiliary::fcmp(lower_ci, -0.9, EPSILON) >= 0) && (auxiliary::fcmp(upper_ci, 0.9, EPSILON) <= 0)) {
					w_values_sum -= 0.95;
					w_values[j] += w_values_sum;
				} else {
					w_values[j] += w_values_sum;
				}
			} else {
				w_values[j] += w_values_sum;
			}

			w_value_max = w_values[j] + 0.025 * ((db->n_markers - i - 1u) * (db->n_markers + i - j - j));
			if (auxiliary::fcmp(w_value_max, 0.0, EPSILON) >= 0) {
				updated_breakpoint = j;
			}
		}
	}

	delete ci;
	ci = NULL;

	free(w_values);
	w_values = NULL;
}
