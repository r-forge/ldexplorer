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

#include "include/CIFactory.h"

CIFactory::CIFactory() {

}

CIFactory::~CIFactory() {

}

CI* CIFactory::create(const char* method, unsigned int likelihood_density) throw (Exception) {
	if (auxiliary::strcmp_ignore_case(method, CI::NONE) == 0) {
		return new CI();
	} else if (auxiliary::strcmp_ignore_case(method, CI::CI_WP) == 0) {
		return new CIWP(likelihood_density);
	} else if (auxiliary::strcmp_ignore_case(method, CI::CI_AV) == 0) {
		return new CIAV();
	} else {
		throw Exception(__FILE__, __LINE__, "Unknown D' CI computation method '%s' was specified.", method);
	}
}
