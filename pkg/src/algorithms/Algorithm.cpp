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

#include "include/Algorithm.h"

const unsigned int Algorithm::STRONG_PAIRS_SIZE_INIT = 100000;
const unsigned int Algorithm::STRONG_PAIRS_SIZE_INCREMENT = 1000;
const unsigned int Algorithm::BLOCKS_SIZE_INIT = 10000;
const unsigned int Algorithm::BLOCKS_SIZE_INCREMENT = 1000;

const double Algorithm::EPSILON = 0.000000001;

Algorithm::Algorithm() throw (Exception) :
		db(NULL),
		pos_strong_pair_cl(0.7), neg_strong_pair_cl(-0.7),
		pos_strong_pair_cu(0.98), neg_strong_pair_cu(-0.98),
		pos_recomb_pair_cu(0.9), neg_recomb_pair_cu(0.9),
		strong_pairs_fraction(0.95), strong_pair_weight(0.05), recomb_pair_weight(0.95),
		strong_pairs(NULL), n_strong_pairs(0u), strong_pairs_size(STRONG_PAIRS_SIZE_INIT),
		blocks(NULL), n_blocks(0u), blocks_size(BLOCKS_SIZE_INIT) {

	strong_pairs = (pair*)malloc(strong_pairs_size * sizeof(pair));
	if (strong_pairs == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	blocks = (unsigned int*)malloc(blocks_size * sizeof(unsigned int));
	if (blocks == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}
}

Algorithm::~Algorithm() {
	db = NULL;

	free(strong_pairs);
	strong_pairs = NULL;

	free(blocks);
	blocks = NULL;
}

void Algorithm::set_dbview(const DbView* db) {
	this->db = db;
}

void Algorithm::set_strong_pair_cl(double ci_lower_bound) {
	pos_strong_pair_cl = ci_lower_bound;
	neg_strong_pair_cl = -pos_strong_pair_cl;
}

void Algorithm::set_strong_pair_cu(double ci_upper_bound) {
	pos_strong_pair_cu = ci_upper_bound;
	neg_strong_pair_cu = -pos_strong_pair_cu;
}

void Algorithm::set_recomb_pair_cu(double ci_upper_bound) {
	pos_recomb_pair_cu = ci_upper_bound;
	neg_recomb_pair_cu = -pos_recomb_pair_cu;
}

void Algorithm::set_strong_pairs_fraction(double fraction) {
	strong_pairs_fraction = fraction;
	strong_pair_weight = 1.0 - strong_pairs_fraction;
	recomb_pair_weight = strong_pairs_fraction;
}

void Algorithm::sort_preliminary_blocks() {
	qsort(strong_pairs, n_strong_pairs, sizeof(pair), paircmp);
}

void Algorithm::select_final_blocks() throw (Exception) {
	unsigned int* new_blocks = NULL;
	bool* used_markers = NULL;

	unsigned int first = 0u;
	unsigned int last = 0u;

	n_blocks = 0u;

	used_markers = (bool*)malloc(db->n_markers * sizeof(bool));
	if (used_markers == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < db->n_markers; ++i) {
		used_markers[i] = false;
	}

	for (unsigned int b = 0u; b < n_strong_pairs; ++b) {
		first = strong_pairs[b].first;
		last = strong_pairs[b].last;

		if (used_markers[first] || used_markers[last]) {
			continue;
		}

		if (n_blocks >= blocks_size) {
			blocks_size += BLOCKS_SIZE_INCREMENT;
			new_blocks = (unsigned int*)realloc(blocks, blocks_size * sizeof(unsigned int));
			if (new_blocks == NULL) {
				throw Exception(__FILE__, __LINE__, "Error in memory reallocation.");
			}
			blocks = new_blocks;
			new_blocks = NULL;
		}

		blocks[n_blocks] = b;
		++n_blocks;

		for (unsigned int i = first; i <= last; ++i) {
			used_markers[i] = true;
		}
	}

	free(used_markers);
	used_markers = NULL;
}

unsigned int Algorithm::get_n_strong_pairs() {
	return n_strong_pairs;
}

unsigned int Algorithm::get_n_blocks() {
	return n_blocks;
}

//void Algorithm::write_blocks(const char* output_file_name,
//		const char* input_phase_file_name, const char* input_map_file_name,
//		double maf_threshold, bool region, unsigned long int start, unsigned long int end,
//		const char* ci_method) throw (Exception) {

void Algorithm::write_blocks(const char* output_file_name) throw (Exception) {
	Writer* writer = NULL;

	pair strong_pair;

	const char* first_marker = NULL;
	const char* last_marker = NULL;
	unsigned int start_bp = 0u;
	unsigned int end_bp = 0u;
	unsigned int n_markers = 0u;
	unsigned int n_haps = 0u;
	unsigned int n_unique_haps = 0u;
	unsigned int n_common_haps = 0u;
	double haps_diversity = 0.0;

	try {
		writer = WriterFactory::create(WriterFactory::TEXT);
		writer->set_file_name(output_file_name);
		writer->open(true);

//		writer->write("#PHASE FILE: %s\n", input_phase_file_name);
//		if (input_map_file_name != NULL) {
//			writer->write("#MAP FILE: %s\n", input_map_file_name);
//		}
//		if (region) {
//			writer->write("#REGION: [%u, %u]\n", start, end);
//		}
//		writer->write("#HAPLOTYPES: %u\n", db->n_haplotypes);
//		writer->write("#MAF > %g\n", maf_threshold);
//		writer->write("#ALL SNPs: %u\n", db->all_n_markers);
//		writer->write("#FILTERED SNPs: %u\n", db->n_markers);
//		writer->write("#D' CI COMPUTATION METHOD: %s\n", ci_method);

		writer->write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				"BLOCK_NAME", "FIRST_SNP", "LAST_SNP", "FIRST_SNP_ID", "LAST_SNP_ID", "START_BP", "END_BP", "N_SNPS", "N_HAPS", "N_UNIQUE_HAPS", "N_COMMON_HAPS", "HAPS_DIVERSITY");

		for (unsigned int b = 0u; b < n_blocks; ++b) {
			strong_pair = strong_pairs[blocks[b]];

			first_marker = db->markers[strong_pair.first];
			last_marker = db->markers[strong_pair.last];
			start_bp = db->positions[strong_pair.first];
			end_bp = db->positions[strong_pair.last];
			n_markers = strong_pair.last - strong_pair.first + 1u;

			get_block_diversity(b, &n_haps, &n_unique_haps, &n_common_haps, &haps_diversity);

			writer->write("BLOCK_%07u\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%g\n",
					b + 1u, first_marker, last_marker, strong_pair.first, strong_pair.last, start_bp, end_bp, n_markers, n_haps, n_unique_haps, n_common_haps, haps_diversity);
		}

		writer->close();
		delete writer;
	} catch (Exception &e) {
		if (writer != NULL) {
			delete writer;
		}
		throw;
	}
}

double Algorithm::get_memory_usage_preliminary_blocks() {
	return ((strong_pairs_size * sizeof(pair)) / 1048576.0);
}

double Algorithm::get_memory_usage_final_blocks() {
	return ((blocks_size * sizeof(unsigned int)) / 1048576.0);
}

double Algorithm::get_memory_usage() {
	return 0.0;
}

double Algorithm::get_total_memory_usage() {
	double memory_usage = 0.0;

	memory_usage += get_memory_usage_preliminary_blocks();
	memory_usage += get_memory_usage_final_blocks();
	memory_usage += get_memory_usage();

	return memory_usage;
}

int Algorithm::paircmp(const void* first, const void* second) {
	pair* first_pair = (pair*)first;
	pair* second_pair = (pair*)second;

	if (first_pair->distance > second_pair->distance) {
		return -1;
	} else if (first_pair->distance < second_pair->distance) {
		return 1;
	} else {
		return first_pair->first - second_pair->first;
	}
}

bool Algorithm::is_compatible_haplotype(const char* first, const char* second) {
	unsigned int i = 0u;
	unsigned int length = max(strlen(first), strlen(second));
	char first_char = '\0';
	char second_char = '\0';

	while (i < length) {
		first_char = tolower(first[i]);
		if ((first_char != 'a') && (first_char != 'c') && (first_char != 'g') && (first_char != 't')) {
			++i;
			continue;
		}

		second_char = tolower(second[i]);
		if ((second_char != 'a') && (second_char != 'c') && (second_char != 'g') && (second_char != 't')) {
			++i;
			continue;
		}

		if (first_char != second_char) {
			return false;
		}

		++i;
	}

	if ((first[i] == '\0') && (second[i] == '\0')) {
		return true;
	}

	return false;
}

void Algorithm::get_block_diversity(unsigned int block_id, unsigned int* n_haps, unsigned int* n_unique_haps, unsigned int* n_common_haps, double* haps_diversity) throw (Exception) {
	pair strong_pair;

	unsigned int n_markers = 0u;
	unsigned int n_all_common_haps = 0u;

	char* hap = NULL;
	char* hap_new = NULL;

	map<char*, unsigned int, bool(*)(const char*, const char*)> haps(auxiliary::bool_strcmp_ignore_case);
	map<char*, unsigned int, bool(*)(const char*, const char*)>::iterator haps_it;

	vector<const char*> compatible_haps;

	map<const char*, unsigned int, bool(*)(const char*, const char*)> unambiguous_haps(auxiliary::bool_strcmp_ignore_case);
	map<const char*, unsigned int, bool(*)(const char*, const char*)>::iterator unambiguous_haps_it;

	vector< std::pair<const char*, unsigned int> > unique_haps;
	vector< std::pair<const char*, unsigned int> >::iterator unique_haps_it;

	try {
		*n_haps = 0u;
		*n_unique_haps = 0u;
		*n_common_haps = 0u;
		*haps_diversity = 0.0;

		strong_pair = strong_pairs[blocks[block_id]];

		n_markers = strong_pair.last - strong_pair.first + 1u;

		hap = (char*)malloc((n_markers + 1u) * sizeof(char));
		if (hap == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}
		hap[n_markers] = '\0';

		/* Read all (ambiguous & unambiguous) haplotypes */
		for (unsigned int j = 0u; j < db->n_haplotypes; ++j) {
			for (unsigned int i = strong_pair.first, k = 0u; i <= strong_pair.last; ++i, ++k) {
				hap[k] = db->haplotypes[i][j];
			}

			haps_it = haps.find(hap);
			if (haps_it == haps.end()) {
				hap_new = (char*)malloc((n_markers + 1u) * sizeof(char));
				if (hap_new == NULL) {
					throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
				}
				strcpy(hap_new, hap);

				haps.insert(std::pair<char*, unsigned int>(hap_new, 1u));
			} else {
				haps_it->second++;
			}
		}

		/* Select unambiguous haplotypes */
		haps_it = haps.begin();
		while (haps_it != haps.end()) {
			unambiguous_haps.insert(std::pair<const char*, unsigned int>(haps_it->first, haps_it->second));
			++haps_it;
		}

		haps_it = haps.begin();
		while (haps_it != haps.end()) {
			compatible_haps.clear();

			unambiguous_haps_it = unambiguous_haps.begin();
			while (unambiguous_haps_it != unambiguous_haps.end()) {
				if (is_compatible_haplotype(haps_it->first, unambiguous_haps_it->first)) {
					compatible_haps.push_back(unambiguous_haps_it->first);
				}
				++unambiguous_haps_it;
			}

			for (unsigned int j = 1u; j < compatible_haps.size(); ++j) {
				for (unsigned int i = 0u; i < j; ++i) {
					if (!is_compatible_haplotype(compatible_haps.at(i), compatible_haps.at(j))) {
						unambiguous_haps.erase(haps_it->first);
						break;
					}
				}
			}

			++haps_it;
		}

		/* Group unambiguous haplotypes by compatibility */
		unambiguous_haps_it = unambiguous_haps.begin();
		while (unambiguous_haps_it != unambiguous_haps.end()) {

			unique_haps_it = unique_haps.begin();
			while (unique_haps_it != unique_haps.end()) {
				if (is_compatible_haplotype(unambiguous_haps_it->first, unique_haps_it->first)) {
					unique_haps_it->second += unambiguous_haps_it->second;
					break;
				}
				++unique_haps_it;
			}

			if (unique_haps_it == unique_haps.end()) {
				unique_haps.push_back(std::pair<const char*, unsigned int>(unambiguous_haps_it->first, unambiguous_haps_it->second));
			}

			++unambiguous_haps_it;
		}

		unique_haps_it = unique_haps.begin();
		while (unique_haps_it != unique_haps.end()) {
			*n_haps += unique_haps_it->second;
			if (unique_haps_it->second > 1u) {
				*n_common_haps += 1u;
				n_all_common_haps += unique_haps_it->second;
			}
			++unique_haps_it;
		}

		*n_unique_haps = unique_haps.size();
		*haps_diversity = ((double)n_all_common_haps) / ((double)*n_haps);

		compatible_haps.clear();
		unambiguous_haps.clear();
		unique_haps.clear();

		for (haps_it = haps.begin(); haps_it != haps.end(); ++haps_it) {
			free(haps_it->first);
		}
		haps.clear();

		free(hap);
	} catch (Exception &e) {
		compatible_haps.clear();
		unambiguous_haps.clear();
		unique_haps.clear();

		for (haps_it = haps.begin(); haps_it != haps.end(); ++haps_it) {
			free(haps_it->first);
		}
		haps.clear();

		if (hap != NULL) {
			free(hap);
		}

		throw;
	}
}
