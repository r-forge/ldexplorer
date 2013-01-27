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

#include "include/CIFactory.h"

const char* CIFactory::NONE = "NONE";
const char* CIFactory::WP_CI = "WP";
const char* CIFactory::AV_CI = "AV";

CIFactory::CIFactory() {

}

CIFactory::~CIFactory() {

}

CI* CIFactory::create(Db& db, const char* method, unsigned int precision) throw (Exception) {
	if (auxiliary::strcmp_ignore_case(method, NONE) == 0) {
		return new CI(db);
	} else if (auxiliary::strcmp_ignore_case(method, WP_CI) == 0) {
		return new CIWP(db, precision);
	} else if (auxiliary::strcmp_ignore_case(method, AV_CI) == 0) {
		return new CIAV(db);
	} else {
		throw Exception(__FILE__, __LINE__, "Unknown D' CI computation method '%s' was specified.", method);
	}
}
