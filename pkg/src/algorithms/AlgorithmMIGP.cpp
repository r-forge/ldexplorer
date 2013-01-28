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

#include "include/AlgorithmMIGP.h"

AlgorithmMIGP::AlgorithmMIGP(Db& db) : Algorithm(db) {
	cout << "MIG+ algorithm." << endl;
}

AlgorithmMIGP::~AlgorithmMIGP() {

}

void AlgorithmMIGP::compute_preliminary_blocks(const char* ci_method, unsigned int ci_precision, unsigned int window) throw (Exception) {
	CI* ci = NULL;

	ci = CIFactory::create(*db, ci_method, ci_precision);

	cout << "MIG+ computations." << endl;
}
