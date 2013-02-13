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

#ifndef WRITER_H_
#define WRITER_H_

#include <cstdarg>
#include <cstring>
#include "../../exception/include/Exception.h"

using namespace std;

class Writer {
protected:
	char* file_name;

	Writer();

public:
	virtual ~Writer();

	void set_file_name(const char* file_name) throw (Exception);
	const char* get_file_name();

	virtual void open(bool append) throw (Exception) = 0;
	virtual void close() throw (Exception) = 0;
	virtual void write(const char* format, ...) throw (Exception) = 0;
};

#endif
