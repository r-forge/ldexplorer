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

#ifndef CIWP_H_
#define CIWP_H_

#include "CI.h"
#include "../../auxiliary/include/auxiliary.h"

class CIWP: public CI {
private:
	unsigned int generation_number;
	double* generated_dprime;
	double generated_freq_haplotype_ref_a_ref_b;
	double generated_freq_haplotype_ref_a_alt_b;
	double generated_freq_haplotype_alt_a_ref_b;
	double generated_freq_haplotype_alt_a_alt_b;

	double* log_likelihood;
	double max_log_likelihood;

	double* posterior_dist;
	double total_posterior_dist_area;
	double tail_posterior_dist_area;
	double covered_posterior_dist_area;

	unsigned int tmp_n_observed_haplotype_ref_a_ref_b;
	unsigned int tmp_n_observed_haplotype_alt_a_alt_b;

	double dmax;

public:
	CIWP(const DbView* db, unsigned int precision) throw (Exception);
	virtual ~CIWP();

	void get_CI(unsigned int marker_a, unsigned int marker_b, double* dprime_lower_ci, double* dprime_upper_ci);
};

#endif
