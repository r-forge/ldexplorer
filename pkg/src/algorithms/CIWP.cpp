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

#include "include/CIWP.h"

CIWP::CIWP(Db& db, unsigned int precision) throw (Exception) : CI(db),
		generation_number(precision),
		generated_dprime(NULL),
		generated_freq_haplotype_ref_a_ref_b(0.0), generated_freq_haplotype_ref_a_alt_b(0.0), generated_freq_haplotype_alt_a_ref_b(0.0), generated_freq_haplotype_alt_a_alt_b(0.0),
		log_likelihood(NULL), max_log_likelihood(0.0),
		posterior_dist(NULL), total_posterior_dist_area(0.0), tail_posterior_dist_area(0.0), covered_posterior_dist_area(0.0),
		tmp_n_observed_haplotype_ref_a_ref_b(0u), tmp_n_observed_haplotype_alt_a_alt_b(0u),
		dmax(0.0) {

	generated_dprime = (double*)malloc((generation_number + 1) * sizeof(double));
	if (generated_dprime == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation");
	}

	for (unsigned int i = 0u; i <= generation_number; i++) {
		generated_dprime[i] = i / (double)generation_number;
	}

	posterior_dist = (double*)malloc((generation_number + 1) * sizeof(double));
	if (posterior_dist == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	log_likelihood = (double*)malloc((generation_number + 1) * sizeof(double));
	if (log_likelihood == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}
}

CIWP::~CIWP() {
	free(generated_dprime);
	generated_dprime = NULL;

	free(posterior_dist);
	posterior_dist = NULL;

	free(log_likelihood);
	log_likelihood = NULL;
}

void CIWP::get_CI(unsigned int marker_a, unsigned int marker_b, double* dprime_lower_ci, double* dprime_upper_ci) {
	observed_ref_allele_a = db->major_alleles[marker_a];
	observed_alt_allele_a = db->minor_alleles[marker_a];

	observed_ref_allele_b = db->major_alleles[marker_b];
	observed_alt_allele_b = db->minor_alleles[marker_b];

	observed_major_af_a = db->major_allele_freqs[marker_a];
	observed_major_af_b = db->major_allele_freqs[marker_b];

	observed_haplotype_a = db->haplotypes[marker_a];
	observed_haplotype_b = db->haplotypes[marker_b];

	n_observed_haplotype_ref_a_ref_b = n_observed_haplotype_ref_a_alt_b = n_observed_haplotype_alt_a_ref_b = n_observed_haplotype_alt_a_alt_b = 0u;

	for (unsigned int i = 0u; i < db->n_haplotypes; ++i) {
		observed_allele_a = observed_haplotype_a[i];
		observed_allele_b = observed_haplotype_b[i];
		if (observed_allele_a == observed_ref_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				++n_observed_haplotype_ref_a_ref_b;
			} else if (observed_allele_b == observed_alt_allele_b) {
				++n_observed_haplotype_ref_a_alt_b;
			}
		} else if (observed_allele_a == observed_alt_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				++n_observed_haplotype_alt_a_ref_b;
			} else if (observed_allele_b == observed_alt_allele_b) {
				++n_observed_haplotype_alt_a_alt_b;
			}
		}
	}

	observed_d = (n_observed_haplotype_ref_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b)) - (observed_major_af_a * observed_major_af_b);

	switch (auxiliary::fcmp(observed_d, 0.0, EPSILON)) {
		case -1:
			tmp_n_observed_haplotype_ref_a_ref_b = n_observed_haplotype_ref_a_ref_b;
			tmp_n_observed_haplotype_alt_a_alt_b = n_observed_haplotype_alt_a_alt_b;

			n_observed_haplotype_ref_a_ref_b = n_observed_haplotype_ref_a_alt_b;
			n_observed_haplotype_alt_a_alt_b = n_observed_haplotype_alt_a_ref_b;

			n_observed_haplotype_ref_a_alt_b = tmp_n_observed_haplotype_ref_a_ref_b;
			n_observed_haplotype_alt_a_ref_b = tmp_n_observed_haplotype_alt_a_alt_b;

			observed_major_af_b = 1.0 - observed_major_af_b;

			observed_d = (n_observed_haplotype_ref_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b)) - (observed_major_af_a * observed_major_af_b);
		case 1:
			dmax = min(observed_major_af_a * (1 - observed_major_af_b), (1 - observed_major_af_a) * observed_major_af_b);
			break;
		default:
			*dprime_lower_ci = *dprime_upper_ci = numeric_limits<double>::quiet_NaN();
			return;
	}

	max_log_likelihood = -numeric_limits<double>::infinity();
	total_posterior_dist_area = 0.0;
	for (unsigned int i = 0u; i <= generation_number; i++) {
		generated_freq_haplotype_ref_a_ref_b = generated_dprime[i] * dmax + observed_major_af_a * observed_major_af_b;
		generated_freq_haplotype_ref_a_alt_b = observed_major_af_a - generated_freq_haplotype_ref_a_ref_b;
		generated_freq_haplotype_alt_a_ref_b = observed_major_af_b - generated_freq_haplotype_ref_a_ref_b;
		generated_freq_haplotype_alt_a_alt_b = (1 - observed_major_af_a) - generated_freq_haplotype_alt_a_ref_b;

		log_likelihood[i] =
				n_observed_haplotype_ref_a_ref_b * log10(auxiliary::fcmp(generated_freq_haplotype_ref_a_ref_b, 0.0, EPSILON) <= 0 ? 0.0000000001 : generated_freq_haplotype_ref_a_ref_b) +
				n_observed_haplotype_ref_a_alt_b * log10(auxiliary::fcmp(generated_freq_haplotype_ref_a_alt_b, 0.0, EPSILON) <= 0 ? 0.0000000001 : generated_freq_haplotype_ref_a_alt_b) +
				n_observed_haplotype_alt_a_ref_b * log10(auxiliary::fcmp(generated_freq_haplotype_alt_a_ref_b, 0.0, EPSILON) <= 0 ? 0.0000000001 : generated_freq_haplotype_alt_a_ref_b) +
				n_observed_haplotype_alt_a_alt_b * log10(auxiliary::fcmp(generated_freq_haplotype_alt_a_alt_b, 0.0, EPSILON) <= 0 ? 0.0000000001 : generated_freq_haplotype_alt_a_alt_b);

		if (log_likelihood[i] > max_log_likelihood) {
			max_log_likelihood = log_likelihood[i];
		}
	}

	for (unsigned int i = 0u; i <= generation_number; i++) {
		log_likelihood[i] -= max_log_likelihood;
		posterior_dist[i] = pow(10.0, log_likelihood[i]);
		total_posterior_dist_area += posterior_dist[i];
	}

	tail_posterior_dist_area = 0.05 * total_posterior_dist_area;

	covered_posterior_dist_area = 0.0;
	for (unsigned int i = 0u; i <= generation_number; i++) {
		covered_posterior_dist_area += posterior_dist[i];
		if (covered_posterior_dist_area > tail_posterior_dist_area) {
			if (i != 0u) {
				*dprime_lower_ci = generated_dprime[i - 1u];
			} else {
				*dprime_lower_ci = generated_dprime[0u];
			}
			break;
		}
	}

	covered_posterior_dist_area = 0.0;
	for (unsigned int i = generation_number; i >= 0u; i--) {
		covered_posterior_dist_area += posterior_dist[i];
		if (covered_posterior_dist_area > tail_posterior_dist_area) {
			if (i != generation_number) {
				*dprime_upper_ci =  generated_dprime[i + 1u];
			} else {
				*dprime_upper_ci =  generated_dprime[generation_number];
			}
			break;
		}
	}

//	cout << marker_a << "\t" << marker_b << "\t" << observed_d << "\t" << (observed_d / dmax) << "\t" << "[" << *dprime_lower_ci << ", " << *dprime_upper_ci << "]" << endl;
}
