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

#include "include/CI.h"

const char* CI::NONE = "NONE";
const char* CI::CI_WP = "WP";
const char* CI::CI_AV = "AV";

const double CI::EPSILON = 0.000000001;

CI::CI() :
		db(NULL),
		observed_haplotype_a(NULL), observed_haplotype_b(NULL),
		n_observed_haplotype_ref_a_ref_b(0u), n_observed_haplotype_ref_a_alt_b(0u), n_observed_haplotype_alt_a_ref_b(0u), n_observed_haplotype_alt_a_alt_b(0u),
		observed_allele_a('\0'), observed_allele_b('\0'),
		observed_ref_allele_a('\0'), observed_alt_allele_a('\0'), observed_ref_allele_b('\0'), observed_alt_allele_b('\0'),
		observed_major_af_a(0.0), observed_major_af_b(0.0),
		observed_d(0.0) {

}

CI::~CI() {
	db = NULL;
}

void CI::set_dbview(const DbView* db) {
	this->db = db;
}

double CI::get_D(unsigned int marker_a, unsigned int marker_b) {
	observed_haplotype_a = db->haplotypes[marker_a];
	observed_haplotype_b = db->haplotypes[marker_b];

	n_observed_haplotype_ref_a_ref_b = n_observed_haplotype_ref_a_alt_b = n_observed_haplotype_alt_a_ref_b = n_observed_haplotype_alt_a_alt_b = 0u;

	observed_ref_allele_a = db->major_alleles[marker_a];
	observed_alt_allele_a = db->minor_alleles[marker_a];

	observed_ref_allele_b = db->major_alleles[marker_b];
	observed_alt_allele_b = db->minor_alleles[marker_b];

	observed_major_af_a = db->major_allele_freqs[marker_a];
	observed_major_af_b = db->major_allele_freqs[marker_b];

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

	return (n_observed_haplotype_ref_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b)) - (observed_major_af_a * observed_major_af_b);
}

double CI::get_Dprime(unsigned int marker_a, unsigned int marker_b) {
	observed_haplotype_a = db->haplotypes[marker_a];
	observed_haplotype_b = db->haplotypes[marker_b];

	n_observed_haplotype_ref_a_ref_b = n_observed_haplotype_ref_a_alt_b = n_observed_haplotype_alt_a_ref_b = n_observed_haplotype_alt_a_alt_b = 0u;

	observed_ref_allele_a = db->major_alleles[marker_a];
	observed_alt_allele_a = db->minor_alleles[marker_a];

	observed_ref_allele_b = db->major_alleles[marker_b];
	observed_alt_allele_b = db->minor_alleles[marker_b];

	observed_major_af_a = db->major_allele_freqs[marker_a];
	observed_major_af_b = db->major_allele_freqs[marker_b];

	for (unsigned int i = 0u; i < db->n_haplotypes; i++) {
		observed_allele_a = observed_haplotype_a[i];
		observed_allele_b = observed_haplotype_b[i];
		if (observed_allele_a == observed_ref_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				n_observed_haplotype_ref_a_ref_b++;
			} else if (observed_allele_b == observed_alt_allele_b) {
				n_observed_haplotype_ref_a_alt_b++;
			}
		} else if (observed_allele_a == observed_alt_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				n_observed_haplotype_alt_a_ref_b++;
			} else if (observed_allele_b == observed_alt_allele_b) {
				n_observed_haplotype_alt_a_alt_b++;
			}
		}
	}

	observed_d = (n_observed_haplotype_ref_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b)) - (observed_major_af_a * observed_major_af_b);

	if (observed_d < 0.0) {
		return observed_d / min(observed_major_af_a * observed_major_af_b, (1 - observed_major_af_a) * (1 - observed_major_af_b));
	} else if (observed_d > 0.0) {
		return observed_d / min(observed_major_af_a * (1 - observed_major_af_b), (1 - observed_major_af_a) * observed_major_af_b);
	} else {
		return numeric_limits<double>::quiet_NaN();
	}
}

double CI::get_r(unsigned int marker_a, unsigned int marker_b) {
	observed_haplotype_a = db->haplotypes[marker_a];
	observed_haplotype_b = db->haplotypes[marker_b];

	n_observed_haplotype_ref_a_ref_b = n_observed_haplotype_ref_a_alt_b = n_observed_haplotype_alt_a_ref_b = n_observed_haplotype_alt_a_alt_b = 0u;

	observed_ref_allele_a = db->major_alleles[marker_a];
	observed_alt_allele_a = db->minor_alleles[marker_a];

	observed_ref_allele_b = db->major_alleles[marker_b];
	observed_alt_allele_b = db->minor_alleles[marker_b];

	observed_major_af_a = db->major_allele_freqs[marker_a];
	observed_major_af_b = db->major_allele_freqs[marker_b];

	for (unsigned int i = 0u; i < db->n_haplotypes; i++) {
		observed_allele_a = observed_haplotype_a[i];
		observed_allele_b = observed_haplotype_b[i];
		if (observed_allele_a == observed_ref_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				n_observed_haplotype_ref_a_ref_b++;
			} else if (observed_allele_b == observed_alt_allele_b) {
				n_observed_haplotype_ref_a_alt_b++;
			}
		} else if (observed_allele_a == observed_alt_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				n_observed_haplotype_alt_a_ref_b++;
			} else if (observed_allele_b == observed_alt_allele_b) {
				n_observed_haplotype_alt_a_alt_b++;
			}
		}
	}

	observed_d = (n_observed_haplotype_ref_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b)) - (observed_major_af_a * observed_major_af_b);

	return observed_d / sqrt(observed_major_af_a * (1.0 - observed_major_af_a) * observed_major_af_b * (1.0 - observed_major_af_b));
}

double CI::get_rsq(unsigned int marker_a, unsigned int marker_b) {
	observed_haplotype_a = db->haplotypes[marker_a];
	observed_haplotype_b = db->haplotypes[marker_b];

	n_observed_haplotype_ref_a_ref_b = n_observed_haplotype_ref_a_alt_b = n_observed_haplotype_alt_a_ref_b = n_observed_haplotype_alt_a_alt_b = 0u;

	observed_ref_allele_a = db->major_alleles[marker_a];
	observed_alt_allele_a = db->minor_alleles[marker_a];

	observed_ref_allele_b = db->major_alleles[marker_b];
	observed_alt_allele_b = db->minor_alleles[marker_b];

	observed_major_af_a = db->major_allele_freqs[marker_a];
	observed_major_af_b = db->major_allele_freqs[marker_b];

	for (unsigned int i = 0u; i < db->n_haplotypes; i++) {
		observed_allele_a = observed_haplotype_a[i];
		observed_allele_b = observed_haplotype_b[i];
		if (observed_allele_a == observed_ref_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				n_observed_haplotype_ref_a_ref_b++;
			} else if (observed_allele_b == observed_alt_allele_b) {
				n_observed_haplotype_ref_a_alt_b++;
			}
		} else if (observed_allele_a == observed_alt_allele_a) {
			if (observed_allele_b == observed_ref_allele_b) {
				n_observed_haplotype_alt_a_ref_b++;
			} else if (observed_allele_b == observed_alt_allele_b) {
				n_observed_haplotype_alt_a_alt_b++;
			}
		}
	}

	observed_d = (n_observed_haplotype_ref_a_ref_b / (double)(n_observed_haplotype_ref_a_ref_b + n_observed_haplotype_ref_a_alt_b + n_observed_haplotype_alt_a_ref_b + n_observed_haplotype_alt_a_alt_b)) - (observed_major_af_a * observed_major_af_b);

	return (observed_d * observed_d) / (observed_major_af_a * (1.0 - observed_major_af_a) * observed_major_af_b * (1.0 - observed_major_af_b));
}

void CI::get_CI(unsigned int marker_a, unsigned int marker_b, double* dprime_lower_ci, double* dprime_upper_ci) {

}
