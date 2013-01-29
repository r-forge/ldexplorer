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

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <exception>
#include <limits>
#include <list>
#include <cstdarg>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cstdio>

using namespace std;

class Exception: public exception {
private:
	static const int MESSAGE_BUFFER_SIZE;

	struct message {
		char* text;
		const char* source;
		int source_line;
		message(): text(NULL), source(NULL), source_line(numeric_limits<int>::min()) {};
	};

	list<message*> trace;

	void format_message_text(char** text, const char* text_template, va_list arguments);
	void add_message(const char* source, int source_line, const char* text_message, va_list arguments);

public:
	Exception(const char* source, int source_line, const char* text_message, ... );
	virtual ~Exception() throw();
	void add_message(const char* source, int source_line, const char* text_message, ... );
	virtual const char* what() const throw();
};

#endif
