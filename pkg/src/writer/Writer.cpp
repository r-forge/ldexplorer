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

#include "include/Writer.h"

Writer::Writer() : file_name(NULL) {

}

Writer::~Writer() {
	free(file_name);
	file_name = NULL;
}

void Writer::set_file_name(const char* file_name) throw (Exception) {
	free(this->file_name);
	this->file_name = NULL;
	this->file_name = (char*)malloc((strlen(file_name) + 1) * sizeof(char));
	if (this->file_name == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}
	strcpy(this->file_name, file_name);
}

const char* Writer::get_file_name() {
	return file_name;
}
