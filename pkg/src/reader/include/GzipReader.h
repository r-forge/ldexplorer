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

#ifndef GZIPREADER_H_
#define GZIPREADER_H_

#include "Reader.h"
#include "../../zlib/zlib.h"

class GzipReader: public Reader {
private:
	gzFile infile;

	int buffer_size;
	char* buffer;

public:
	static const unsigned int DEFAULT_BUFFER_SIZE;

	GzipReader(unsigned int buffer_size = DEFAULT_BUFFER_SIZE) throw (Exception);
	virtual ~GzipReader();

	void open() throw (Exception);
	void close() throw (Exception);
	int read_line() throw (Exception);
	void reset() throw (Exception);
	bool eof();
	bool is_compressed();
};

#endif
