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

#include "include/LD.h"

const char* LD::D = "D";
const char* LD::DPRIME = "DPRIME";
const char* LD::R2 = "R2";

LD::LD() : db(NULL), marker_index(NULL), marker_index_size(0u) {

}

LD::~LD() {
	db = NULL;

	variants.clear();
}

void LD::set_dbview(const DbView* db) {
	this->db = db;

	for (variants_it = variants.begin(); variants_it != variants.end(); variants_it++) {
		variants_it->free_variant();
	}
	variants.clear();

	if (marker_index != NULL) {
		free(marker_index);
		marker_index = NULL;
	}
}

void LD::load_markers(const char* file_name) throw (Exception) {
	Reader* reader = NULL;

	char* line = NULL;
	int line_length = 0;

	unsigned int tokens_number = 0u;
	char* token = 0;
	char* tokens[2] = {NULL, NULL};
	char* end_ptr = NULL;

	try {
		reader = ReaderFactory::create(file_name);
		reader->open();

		variants.clear();

		variant temp_variant;

		while ((line_length = reader->read_line()) >= 0) {
			if (line_length == 0) {
				continue;
			}

			line = *(reader->line);

			auxiliary::trim(line);

			if (auxiliary::strcmp_ignore_case(line, "rs", 2) == 0) {
				temp_variant.chromosome = NULL;
				temp_variant.position = 0ul;
				temp_variant.name = (char*)malloc((strlen(line) + 1u) * sizeof(char));
				if (temp_variant.name == NULL) {
					throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
				}
				strcpy(temp_variant.name, line);
			} else {
				tokens_number = 0u;
				tokens[0u] = tokens[1u] = NULL;
				while ((tokens_number < 2u) && ((token = auxiliary::strtok(&line, ':')) != NULL)) {
					tokens[tokens_number++] = token;
				}

				if (tokens[1u] != NULL) {
					temp_variant.chromosome = (char*)malloc((strlen(tokens[0u]) + 1u) * sizeof(char));
					if (temp_variant.chromosome == NULL) {
						throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
					}
					strcpy(temp_variant.chromosome, tokens[0u]);
					temp_variant.position = strtoul(tokens[1u], &end_ptr, 10);
					if (*end_ptr != '\0') {
						throw Exception(__FILE__, __LINE__, "The SNP position '%s' could not be parsed to unsigned integer.", tokens[1u]);
					}
					temp_variant.name = NULL;
				} else {
					temp_variant.chromosome = NULL;
					temp_variant.position = strtoul(tokens[0u], &end_ptr, 10);
					if (*end_ptr != '\0') {
						temp_variant.position = 0ul;
						temp_variant.name = (char*)malloc((strlen(tokens[0u]) + 1u) * sizeof(char));
						if (temp_variant.name == NULL) {
							throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
						}
						strcpy(temp_variant.name, tokens[0u]);
					} else {
						temp_variant.name = NULL;
					}
				}
			}

			variants.push_back(temp_variant);
		}

		temp_variant.chromosome = NULL;
		temp_variant.position = 0u;
		temp_variant.name = NULL;

		reader->close();
		delete reader;
	} catch (Exception &e) {
		e.add_message(__FILE__, __LINE__, "Error while loading '%s' file.", file_name);
		throw;
	}
}

void LD::index_db_markers() throw (Exception) {
	marker_index = (marker_index_entry*)malloc(db->n_markers * sizeof(marker_index_entry));
	if (marker_index == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	marker_index_size = db->n_markers;

	for (unsigned int i = 0u; i < marker_index_size; ++i) {
		marker_index[i].marker = db->markers[i];
		marker_index[i].location = i;
	}

	qsort(marker_index, marker_index_size, sizeof(marker_index_entry), marker_index_entry_cmp);
}

void LD::compute_ld(const char* output_file_name, const char* coefficient, unsigned int window, bool gzip) throw (Exception) {
	Writer* writer = NULL;
	CI* ci = NULL;

	marker_index_entry query_marker;
	unsigned long int window_start = 0ul;
	unsigned long int window_end = 0ul;
	unsigned int start_location = 0;
	unsigned int end_location = 0;
	unsigned int location = 0u;
	long int i = 0;

	try {
		writer = WriterFactory::create(gzip == true ? Writer::GZIP : Writer::TEXT);
		writer->set_file_name(output_file_name);
		writer->open(false);

		ci = CIFactory::create(CI::NONE);
		ci->set_dbview(db);

		writer->write("FIRST_MARKER\tFIRST_BP\tSECOND_MARKER\tSECOND_BP");

		if (auxiliary::strcmp_ignore_case(coefficient, D) == 0) {
			writer->write("\t%s\n", D);

			for (variants_it = variants.begin(); variants_it != variants.end(); ++variants_it) {
				if (variants_it->position > 0ul) {
					if (!lookup_db_by_position(variants_it->position, &location)) {
						continue;
					}

					window_start = db->positions[location] > window ? db->positions[location] - window : 0ul;
					window_end = db->positions[location] + window;

					i = location;
					while ((--i >= 0) && (db->positions[i] >= window_start));

					while ((++i < db->n_markers) && (db->positions[i] <= window_end)) {
						writer->write("%s\t%lu\t%s\t%lu\t%.5f\n", db->markers[location], db->positions[location], db->markers[i], db->positions[i], ci->get_D(location, i));
					}
				} else if (variants_it->name != NULL) {
					query_marker.marker = variants_it->name;
					if (!lookup_db_by_marker(&query_marker, &start_location, &end_location)) {
						continue;
					}

					while (start_location <= end_location) {
						window_start = db->positions[start_location] > window ? db->positions[start_location] - window : 0ul;
						window_end = db->positions[start_location] + window;

						i = start_location;
						while ((--i >= 0) && (db->positions[i] >= window_start));

						while ((++i < db->n_markers) && (db->positions[i] <= window_end)) {
							writer->write("%s\t%lu\t%s\t%lu\t%.5f\n", db->markers[start_location], db->positions[start_location], db->markers[i], db->positions[i], ci->get_D(start_location, i));
						}

						++start_location;
					}
				}
			}
		} else if (auxiliary::strcmp_ignore_case(coefficient, DPRIME) == 0) {
			writer->write("\t%s\n", DPRIME);

			for (variants_it = variants.begin(); variants_it != variants.end(); ++variants_it) {
				if (variants_it->position > 0ul) {
					if (!lookup_db_by_position(variants_it->position, &location)) {
						continue;
					}

					window_start = db->positions[location] > window ? db->positions[location] - window : 0ul;
					window_end = db->positions[location] + window;

					i = location;
					while ((--i >= 0) && (db->positions[i] >= window_start));

					while ((++i < db->n_markers) && (db->positions[i] <= window_end)) {
						writer->write("%s\t%lu\t%s\t%lu\t%.5f\n", db->markers[location], db->positions[location], db->markers[i], db->positions[i], ci->get_Dprime(location, i));
					}
				} else if (variants_it->name != NULL) {
					query_marker.marker = variants_it->name;
					if (!lookup_db_by_marker(&query_marker, &start_location, &end_location)) {
						continue;
					}

					while (start_location <= end_location) {
						window_start = db->positions[start_location] > window ? db->positions[start_location] - window : 0ul;
						window_end = db->positions[start_location] + window;

						i = start_location;
						while ((--i >= 0) && (db->positions[i] >= window_start));

						while ((++i < db->n_markers) && (db->positions[i] <= window_end)) {
							writer->write("%s\t%lu\t%s\t%lu\t%.5f\n", db->markers[start_location], db->positions[start_location], db->markers[i], db->positions[i], ci->get_Dprime(start_location, i));
						}

						++start_location;
					}
				}
			}
		} else if (auxiliary::strcmp_ignore_case(coefficient, R2) == 0) {
			writer->write("\t%s\n", R2);

			for (variants_it = variants.begin(); variants_it != variants.end(); ++variants_it) {
				if (variants_it->position > 0ul) {
					if (!lookup_db_by_position(variants_it->position, &location)) {
						continue;
					}

					window_start = db->positions[location] > window ? db->positions[location] - window : 0ul;
					window_end = db->positions[location] + window;

					i = location;
					while ((--i >= 0) && (db->positions[i] >= window_start));

					while ((++i < db->n_markers) && (db->positions[i] <= window_end)) {
						writer->write("%s\t%lu\t%s\t%lu\t%.5f\n", db->markers[location], db->positions[location], db->markers[i], db->positions[i], pow(ci->get_r(location, i), 2.0));
					}
				} else if (variants_it->name != NULL) {
					query_marker.marker = variants_it->name;
					if (!lookup_db_by_marker(&query_marker, &start_location, &end_location)) {
						continue;
					}

					while (start_location <= end_location) {
						window_start = db->positions[start_location] > window ? db->positions[start_location] - window : 0ul;
						window_end = db->positions[start_location] + window;

						i = start_location;
						while ((--i >= 0) && (db->positions[i] >= window_start));

						while ((++i < db->n_markers) && (db->positions[i] <= window_end)) {
							writer->write("%s\t%lu\t%s\t%lu\t%.5f\n", db->markers[start_location], db->positions[start_location], db->markers[i], db->positions[i], pow(ci->get_r(start_location, i), 2.0));
						}

						++start_location;
					}
				}
			}
		}

		delete ci;
		ci = NULL;

		writer->close();
		delete writer;
	} catch (Exception &e) {
		if (writer != NULL) {
			delete writer;
		}
		if (ci != NULL) {
			delete ci;
		}
		throw;
	}
}

bool LD::lookup_db_by_marker(marker_index_entry* query_marker, unsigned int* start_location, unsigned int* end_location) {
	marker_index_entry* found_marker = NULL;

	long int start = 0;
	long int end = 0;

	found_marker = (marker_index_entry*)bsearch(query_marker, marker_index, marker_index_size, sizeof(marker_index_entry), marker_index_entry_cmp);
	if (found_marker == NULL) {
		return false;
	}

	start = found_marker - marker_index;
	end = start;


	while ((start >= 0) && (marker_index_entry_cmp(&(marker_index[start]), query_marker) == 0)) {
		--start;
	}
	++start;

	while ((++end < marker_index_size) && (marker_index_entry_cmp(&(marker_index[end]), query_marker) == 0)) {}
	--end;

	*start_location = marker_index[start].location;
	*end_location = marker_index[end].location;

	return true;
}

bool LD::lookup_db_by_position(unsigned long int position, unsigned int* location) {
	unsigned long int* found_position = NULL;

	found_position = (unsigned long int*)bsearch(&position, db->positions, db->n_markers, sizeof(unsigned long int), ulong_int_cmp);
	if (found_position == NULL) {
		return false;
	}

	*location = found_position - db->positions;

	return true;
}

unsigned int LD::get_n_snps() {
	return variants.size();
}

double LD::get_used_memory() {
	return (variants.size() * sizeof(variant) + marker_index_size * sizeof(marker_index_entry)) / 1048576.0;
}

int LD::marker_index_entry_cmp(const void* first, const void* second) {
	marker_index_entry* first_entry = (marker_index_entry*)first;
	marker_index_entry* second_entry = (marker_index_entry*)second;

	return auxiliary::strcmp_ignore_case(first_entry->marker, second_entry->marker);
}

int LD::ulong_int_cmp(const void* first, const void* second) {
	unsigned long int first_ulong_int = *(unsigned long int*)first;
	unsigned long int second_ulong_int = *(unsigned long int*)second;

	if (first_ulong_int > second_ulong_int) {
		return 1;
	} else if (first_ulong_int < second_ulong_int) {
		return -1;
	}

	return 0;
}
