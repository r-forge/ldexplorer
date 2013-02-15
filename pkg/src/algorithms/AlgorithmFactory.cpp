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

#include "include/AlgorithmFactory.h"

AlgorithmFactory::AlgorithmFactory() {

}

AlgorithmFactory::~AlgorithmFactory() {

}

Algorithm* AlgorithmFactory::create(const char* name, unsigned int window) throw (Exception) {
	if (auxiliary::strcmp_ignore_case(name, Algorithm::ALGORITHM_MIG) == 0) {
		return new AlgorithmMIG();
	} else if (auxiliary::strcmp_ignore_case(name, Algorithm::ALGORITHM_MIGP) == 0) {
		return new AlgorithmMIGP();
	} else if (auxiliary::strcmp_ignore_case(name, Algorithm::ALGORITHM_MIGPP) == 0) {
		return new AlgorithmMIGPP(window);
	} else {
		throw Exception(__FILE__, __LINE__, "Unknown algorithm '%s' was specified.", name);
	}
}
