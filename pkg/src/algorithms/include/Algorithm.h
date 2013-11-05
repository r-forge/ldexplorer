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

#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include <math.h>

#include "../../auxiliary/include/auxiliary.h"
#include "../../db/include/DbView.h"
#include "../../writer/include/WriterFactory.h"
#include "CIFactory.h"
#include "Partition.h"

using namespace std;

class Algorithm {
protected:
	static const unsigned int PRELIMINARY_BLOCKS_SIZE_INIT;
	static const unsigned int PRELIMINARY_BLOCKS_SIZE_INCREMENT;

	struct preliminary_block {
		unsigned int start;
		unsigned int end;
		unsigned long int length_bp;
	};

	const DbView* db;

	const char* ci_method;
	unsigned int likelihood_density;

	double pos_strong_pair_cl;
	double neg_strong_pair_cl;

	double pos_strong_pair_cu;
	double neg_strong_pair_cu;

	double pos_recomb_pair_cu;
	double neg_recomb_pair_cu;

	double strong_pair_rsq;

	double strong_pairs_fraction;

	double strong_pair_weight;
	double recomb_pair_weight;

	preliminary_block* preliminary_blocks;
	unsigned int n_preliminary_blocks;
	unsigned int preliminary_blocks_size;
	bool rsq_preliminary_blocks;

	static int preliminary_blocks_cmp(const void* first, const void* second);

public:
	static const char* ALGORITHM_MIG;
	static const char* ALGORITHM_MIGP;
	static const char* ALGORITHM_MIGPP;

	static const double EPSILON;

	Algorithm() throw (Exception);
	virtual ~Algorithm();

	void set_dbview(const DbView* db);

	void set_ci_method(const char* ci_method);
	void set_likelihood_density(unsigned int likelihood_density);
	void set_strong_pair_cl(double ci_lower_bound);
	void set_strong_pair_cu(double ci_upper_bound);
	void set_recomb_pair_cu(double ci_upper_bound);
	void set_strong_pair_rsq(double rsq_lower_bound);
	void set_strong_pairs_fraction(double fraction);

	virtual void compute_preliminary_blocks() throw (Exception) = 0;
	virtual void compute_preliminary_blocks_rsq() throw (Exception) = 0;
	unsigned int get_n_preliminary_blocks();
	void sort_preliminary_blocks();

	virtual Partition* get_block_partition() throw (Exception);

	double get_memory_usage_preliminary_blocks();
	virtual double get_memory_usage();
};

#endif
