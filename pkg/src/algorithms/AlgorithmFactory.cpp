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

#include "include/AlgorithmFactory.h"

const char* AlgorithmFactory::ALGORITHM_MIG = "MIG";
const char* AlgorithmFactory::ALGORITHM_MIGP = "MIG+";
const char* AlgorithmFactory::ALGORITHM_MIGPP = "MIG++";

AlgorithmFactory::AlgorithmFactory() {

}

AlgorithmFactory::~AlgorithmFactory() {

}

Algorithm* AlgorithmFactory::create(Db& db, const char* name) throw (Exception) {
	if (auxiliary::strcmp_ignore_case(name, ALGORITHM_MIG) == 0) {
		return new AlgorithmMIG(db);
	} else if (auxiliary::strcmp_ignore_case(name, ALGORITHM_MIGP) == 0) {
		return new AlgorithmMIGP(db);
	} else if (auxiliary::strcmp_ignore_case(name, ALGORITHM_MIGPP) == 0) {
		return new AlgorithmMIGPP(db);
	} else {
		throw Exception(__FILE__, __LINE__, "Unknown algorithm '%s' was specified.", name);
	}
}
