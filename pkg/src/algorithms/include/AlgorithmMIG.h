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

#ifndef ALGORITHMMIG_H_
#define ALGORITHMMIG_H_

#include "Algorithm.h"

using namespace std;

class AlgorithmMIG: public Algorithm {

public:
	AlgorithmMIG();
	virtual ~AlgorithmMIG();

	void compute_preliminary_blocks() throw (Exception);
	void compute_preliminary_blocks_rsq() throw (Exception);

	Partition* get_block_partition() throw (Exception);

	double get_memory_usage();
};

#endif
