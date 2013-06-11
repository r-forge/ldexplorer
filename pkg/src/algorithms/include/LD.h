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

#ifndef LD_H_
#define LD_H_

#include <iostream>
#include <vector>
#include <math.h>
#include "../../db/include/DbView.h"
#include "../../reader/include/ReaderFactory.h"
#include "../../writer/include/WriterFactory.h"
#include "CIFactory.h"
#include "../../auxiliary/include/auxiliary.h"

using namespace std;

class LD {
private:
	const DbView* db;

	struct variant {
		char* chromosome;
		unsigned long int position;
		char* name;

		variant() : chromosome(NULL), position(0ul), name(NULL) {}

		void free_variant() {
			if (chromosome != NULL) {
				free(chromosome);
				chromosome = NULL;
			}

			if (name != NULL) {
				free(name);
				name = NULL;
			}
		}
	};

	struct marker_index_entry {
		const char* marker;
		unsigned int location;

		marker_index_entry() : marker(NULL), location(0u) {}

		~marker_index_entry() {
			marker = NULL;
		}
	};

	vector<variant> variants;
	vector<variant>::iterator variants_it;

	marker_index_entry* marker_index;
	unsigned int marker_index_size;

	bool lookup_db_by_marker(marker_index_entry* query_marker, unsigned int* start_location, unsigned int* end_location);
	bool lookup_db_by_position(unsigned long int position, unsigned int* location);

	static int marker_index_entry_cmp(const void* first, const void* second);
	static int ulong_int_cmp(const void* first, const void* second);

public:
	static const char* D;
	static const char* DPRIME;
	static const char* R2;

	LD();
	virtual ~LD();

	void set_dbview(const DbView* db);

	void load_markers(const char* file_name) throw (Exception);
	void index_db_markers() throw (Exception);
	void compute_ld(const char* output_file_name, const char* coefficient, unsigned int window, bool gzip) throw (Exception);

	unsigned int get_n_snps();

	double get_used_memory();
};

#endif
