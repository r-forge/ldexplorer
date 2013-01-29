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

#ifndef UNIQUE_H_
#define UNIQUE_H_

#include "../../exception/include/Exception.h"
#include "../../auxiliary/include/auxiliary.h"

class Unique {
private:
	char* value;
	unsigned int buffer_size;

	bool set_value(const char* value);
	bool check_value(const char* value);

public:
	static const unsigned int DEFAULT_BUFFER_SIZE;

	typedef bool (Unique::*check_function)(const char* value);

	Unique(unsigned int buffer_size = DEFAULT_BUFFER_SIZE) throw (Exception);
	virtual ~Unique();

	check_function check;

	const char* get_value();
	void reset();
};

#endif
