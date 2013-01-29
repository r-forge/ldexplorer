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

#include "include/Unique.h"

const unsigned int Unique::DEFAULT_BUFFER_SIZE = 8192;

Unique::Unique(unsigned int buffer_size) throw (Exception) : value(NULL), buffer_size(buffer_size), check(NULL) {
	value = (char*)malloc((buffer_size + 1) * sizeof(char));
	if (value == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	value[0] = '\0';

	check = &Unique::set_value;
}

Unique::~Unique() {
	free(value);
	value = NULL;
}

bool Unique::set_value(const char* value) {
	strncpy(this->value, value, buffer_size);
	this->value[buffer_size] = '\0';
	check = &Unique::check_value;
	return true;
}

bool Unique::check_value(const char* value) {
	return (auxiliary::strcmp_ignore_case(value, this->value) == 0);
}

const char* Unique::get_value() {
	return value;
}

void Unique::reset() {
	check = &Unique::set_value;
}
