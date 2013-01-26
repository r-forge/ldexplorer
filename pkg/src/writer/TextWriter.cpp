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

#include "include/TextWriter.h"

const unsigned int TextWriter::DEFAULT_BUFFER_SIZE = 16777216;

TextWriter::TextWriter(unsigned int buffer_size) throw (Exception) : buffer(NULL) {
	buffer = (char*)malloc((buffer_size + 1u) * sizeof(char));
	if (buffer == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	buffer[0] = '\0';
}

TextWriter::~TextWriter() {
	free(buffer);
	buffer = NULL;
}

void TextWriter::open() throw (Exception) {
	if (ofile_stream.is_open()) {
		close();
	}

	ofile_stream.clear();
	ofile_stream.open(file_name, ios::binary);

	if (ofile_stream.fail()) {
		throw Exception(__FILE__, __LINE__, "Error while opening '%s' file.", file_name);
	}
}

void TextWriter::close() throw (Exception) {
	if (ofile_stream.is_open()) {
		ofile_stream.clear();
		ofile_stream.close();

		if (ofile_stream.fail()) {
			throw Exception(__FILE__, __LINE__, "Error while closing '%s' file.", file_name);
		}
	}
}

void TextWriter::write(const char* format, ...) throw (Exception) {
	va_list arguments;

	va_start(arguments, format);
	if (vsprintf(buffer, format, arguments) < 0) {
		throw Exception(__FILE__, __LINE__, "Error while writing '%s' file.", file_name);
	}
	va_end(arguments);

	ofile_stream << buffer;

	if (ofile_stream.fail()) {
		throw Exception(__FILE__, __LINE__, "Error while writing '%s' file.", file_name);
	}
}
