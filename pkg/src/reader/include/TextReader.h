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

#ifndef TEXTREADER_H_
#define TEXTREADER_H_

#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include <limits>

#include "Reader.h"

class TextReader : public Reader {
private:
	ifstream ifile_stream;

	int buffer_size;
	char* buffer;

public:
	static const unsigned int DEFAULT_BUFFER_SIZE;

	TextReader(unsigned int buffer_size = DEFAULT_BUFFER_SIZE) throw (Exception);
	virtual ~TextReader();

	void open() throw (Exception);
	void close() throw (Exception);
	int read_line() throw (Exception);
	void reset() throw (Exception);
	bool eof();
	bool is_compressed();
};

#endif
