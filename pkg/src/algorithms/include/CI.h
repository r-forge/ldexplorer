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

#ifndef ALGORITHMCI_H_
#define ALGORITHMCI_H_

#include "../../db/include/DbView.h"
#include "../../writer/include/WriterFactory.h"

using namespace std;

class CI {
protected:
	const DbView* db;

	char* observed_haplotype_a;
	char* observed_haplotype_b;

	unsigned int n_observed_haplotype_ref_a_ref_b;
	unsigned int n_observed_haplotype_ref_a_alt_b;
	unsigned int n_observed_haplotype_alt_a_ref_b;
	unsigned int n_observed_haplotype_alt_a_alt_b;

	char observed_allele_a;
	char observed_allele_b;

	char observed_ref_allele_a;
	char observed_alt_allele_a;

	char observed_ref_allele_b;
	char observed_alt_allele_b;

	double observed_major_af_a;
	double observed_major_af_b;

	double observed_d;

public:
	static const char* NONE;
	static const char* CI_WP;
	static const char* CI_AV;

	static const double EPSILON;

	CI();
	virtual ~CI();

	void set_dbview(const DbView* db);

	double get_D(unsigned int marker_a, unsigned int marker_b);
	double get_Dprime(unsigned int marker_a, unsigned int marker_b);
	double get_r(unsigned int marker_a, unsigned int marker_b);
	double get_rsq(unsigned int marker_a, unsigned int marker_b);

	virtual void get_CI(unsigned int marker_a, unsigned int marker_b, double* dprime_lower_ci, double* dprime_upper_ci);

};

#endif
