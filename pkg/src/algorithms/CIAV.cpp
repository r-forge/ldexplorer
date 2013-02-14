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

#include "include/CIAV.h"

CIAV::CIAV() :
	var_d(0.0),
	dmax_first(0.0), dmax_second(0.0), dmax(0.0),
	f(0.0),
	dprime(0.0), abs_dprime(0.0), var_dprime(0.0) {
}

CIAV::~CIAV() {

}

void CIAV::get_CI(unsigned int marker_a, unsigned int marker_b, double* dprime_lower_ci, double* dprime_upper_ci) {
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
	var_d = (observed_major_af_a * (1.0 - observed_major_af_a) * observed_major_af_b * (1.0 - observed_major_af_b) + observed_d * ((1.0 - observed_major_af_a) - observed_major_af_a) * ((1.0 - observed_major_af_b) - observed_major_af_b) - observed_d * observed_d) / db->n_haplotypes;

	switch (auxiliary::fcmp(observed_d, 0.0, EPSILON)) {
		case 1:
			dmax_first = observed_major_af_a * (1.0 - observed_major_af_b);
			dmax_second = (1.0 - observed_major_af_a) * observed_major_af_b;

			if (auxiliary::fcmp(dmax_first, dmax_second, EPSILON) <= 0) {
				f = n_observed_haplotype_ref_a_alt_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b);
				dmax = dmax_first;
			} else {
				f = n_observed_haplotype_alt_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b);
				dmax = dmax_second;
			}

			dprime = observed_d / dmax;
			abs_dprime = fabs(dprime);

			var_dprime = (1.0 / (db->n_haplotypes * dmax * dmax)) *
					((1.0 - abs_dprime) * (db->n_haplotypes * var_d - abs_dprime * dmax * (observed_major_af_a * observed_major_af_b + (1.0 - observed_major_af_a) * (1.0 - observed_major_af_b) - 2.0 * fabs(observed_d))) +
							abs_dprime * f * (1.0 - f));

			break;
		case -1:
			dmax_first = observed_major_af_a * observed_major_af_b;
			dmax_second = (1.0 - observed_major_af_a) * (1.0 - observed_major_af_b);

			if (auxiliary::fcmp(dmax_first, dmax_second, EPSILON) <= 0) {
				f = n_observed_haplotype_ref_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b);
				dmax = dmax_first;
			} else {
				f = n_observed_haplotype_alt_a_alt_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b);
				dmax = dmax_second;
			}

			dprime = observed_d / dmax;
			abs_dprime = fabs(dprime);

			var_dprime = (1.0 / (db->n_haplotypes * dmax * dmax)) *
					((1.0 - abs_dprime) * (db->n_haplotypes * var_d - abs_dprime * dmax * (observed_major_af_a * (1.0 - observed_major_af_b) + (1.0 - observed_major_af_a) * observed_major_af_b - 2.0 * fabs(observed_d))) +
							abs_dprime * f * (1.0 - f));

			break;
		default:
			*dprime_lower_ci = *dprime_upper_ci = numeric_limits<double>::quiet_NaN();
			return;
	}

	if (auxiliary::fcmp(var_dprime, 0.0, EPSILON) <= 0) {
		var_dprime = 0.0;
	}

	*dprime_lower_ci = dprime - 1.644854 * sqrt(var_dprime);
	*dprime_upper_ci = dprime + 1.644854 * sqrt(var_dprime);

	if (auxiliary::fcmp(*dprime_lower_ci, -1.0, EPSILON) < 0) {
		*dprime_lower_ci = -1.0;
	}

	if (auxiliary::fcmp(*dprime_upper_ci, 1.0, EPSILON) > 0) {
		*dprime_upper_ci = 1.0;
	}
}
