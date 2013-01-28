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

#include "include/WriterFactory.h"

const char* WriterFactory::TEXT = "TEXT";
const char* WriterFactory::GZIP = "GZIP";

WriterFactory::WriterFactory() {

}

WriterFactory::~WriterFactory() {

}

Writer* WriterFactory::create(const char* type) throw (Exception) {
	if (auxiliary::strcmp_ignore_case(type, TEXT) == 0) {
		return new TextWriter();
	} else if (auxiliary::strcmp_ignore_case(type, GZIP) == 0) {
		return new GzipWriter();
	} else {
		throw Exception(__FILE__, __LINE__, "Unrecognized file type '%s'.", type);
	}
}
