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

#ifndef READER_H_
#define READER_H_

#include <cstdlib>
#include <cstring>
#include "../../exception/include/Exception.h"

using namespace std;

class Reader {
protected:
	char* file_name;
	bool compressed;

	Reader(char** buffer);

public:
	char* const* line;

	virtual ~Reader();

	void set_file_name(const char* file_name) throw (Exception);
	const char* get_file_name();

	virtual void open() throw (Exception) = 0;
	virtual void close() throw (Exception) = 0;
	virtual int read_line() throw (Exception) = 0;
	virtual void reset() throw (Exception) = 0;
	virtual bool eof() = 0;
	virtual bool is_compressed() = 0;
};

#endif
