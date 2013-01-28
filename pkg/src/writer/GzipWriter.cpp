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

#include "include/GzipWriter.h"

const unsigned int GzipWriter::DEFAULT_BUFFER_SIZE = 16777216;

GzipWriter::GzipWriter(unsigned int buffer_size) throw (Exception) : buffer(NULL) {
	buffer = (char*)malloc((buffer_size + 1u) * sizeof(char));
	if (buffer == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	buffer[0] = '\0';
}

GzipWriter::~GzipWriter() {
	free(buffer);
	buffer = NULL;
}

void GzipWriter::open() throw (Exception) {
	outfile = gzopen(file_name, "wb");
	if (outfile == NULL) {
		throw Exception(__FILE__, __LINE__, "Error while opening '%s' file.", file_name);
	}
}

void GzipWriter::close() throw (Exception) {
	int gzerrno = 0;

	gzerrno = gzclose(outfile);
	if (gzerrno != Z_OK) {
		throw Exception(__FILE__, __LINE__, "Error while closing '%s' file.", file_name);
	}
}

void GzipWriter::write(const char* format, ...) throw (Exception) {
	va_list arguments;

	va_start(arguments, format);
	if (vsprintf(buffer, format, arguments) < 0) {
		throw Exception(__FILE__, __LINE__, "Error while writing '%s' file.", file_name);
	}
	va_end(arguments);

	if (gzputs(outfile, buffer) < 0) {
		throw Exception(__FILE__, __LINE__, "Error while writing '%s' file.", file_name);
	}
}
