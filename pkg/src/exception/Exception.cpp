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

#include "include/Exception.h"

const int Exception::MESSAGE_BUFFER_SIZE = 8000;

Exception::Exception(const char* source, int source_line, const char* text_message, ... ): exception() {
	va_list arguments;

	va_start(arguments, text_message);
	add_message(source, source_line, text_message, arguments);
	va_end(arguments);
}

Exception::~Exception() throw() {
	list<message*>::iterator trace_it;

	trace_it = trace.begin();
	while (trace_it != trace.end()) {
		if ((*trace_it)->text != NULL) {
			free((*trace_it)->text);
			(*trace_it)->text = NULL;
			(*trace_it)->source = NULL;
		}
		delete *trace_it;

		trace_it++;
	}

	trace.clear();
}

void Exception::format_message_text(char** text, const char* text_template, va_list arguments) {
	if (*text != NULL) {
		free(*text);
		*text = NULL;
	}

	*text = (char*)malloc(MESSAGE_BUFFER_SIZE * sizeof(char));
	if (*text != NULL) {
		if (vsprintf(*text, text_template, arguments) < 0) {
			free(*text);
			*text = NULL;
		}
	}
}

void Exception::add_message(const char* source, int source_line, const char* text_message, va_list arguments) {
	message* msg = new message();

	format_message_text(&(msg->text), text_message, arguments);
	msg->source = source;
	msg->source_line = source_line;

	trace.push_front(msg);
}

void Exception::add_message(const char* source, int source_line, const char* text_message, ... ) {
	va_list arguments;

	va_start(arguments, text_message);
	add_message(source, source_line, text_message, arguments);
	va_end(arguments);
}

const char* Exception::what() const throw() {
	stringstream string_stream;
	list<message*>::const_iterator trace_it;

	string_stream << setfill(' ');
	trace_it = trace.begin();
	while (trace_it != trace.end()) {
		if ((*trace_it)->source != NULL) {
			string_stream << (*trace_it)->source << " ";
		}

		string_stream << "(" << (*trace_it)->source_line << ") : ";

		if ((*trace_it)->text != NULL) {
			string_stream << (*trace_it)->text;
		} else {
			string_stream << "Exception.";
		}

		string_stream << endl;
		trace_it++;
	}

	return string_stream.str().c_str();
}
