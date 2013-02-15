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

#include "include/Partition.h"

const unsigned int Partition::BLOCKS_SIZE_INIT = 10000;
const unsigned int Partition::BLOCKS_SIZE_INCREMENT = 1000;

Partition::Partition(const DbView* db) throw (Exception) : db(db),
		blocks(NULL), n_blocks(0u), blocks_size(BLOCKS_SIZE_INIT),
		ci_method(NULL), likelihood_density(0u),
		strong_pair_cl(numeric_limits<double>::quiet_NaN()), strong_pair_cu(numeric_limits<double>::quiet_NaN()),
		recomb_pair_cu(numeric_limits<double>::quiet_NaN()), strong_pairs_fraction(numeric_limits<double>::quiet_NaN()),
		pruning_method(NULL), window(0u) {

	blocks = (block*)malloc(blocks_size * sizeof(block));
	if (blocks == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}
}

Partition::~Partition() {
	db = NULL;

	free(blocks);
	blocks = NULL;
}

void Partition::add_block(unsigned int start, unsigned int end) throw (Exception) {

	if (n_blocks >= blocks_size) {
		blocks_size += BLOCKS_SIZE_INCREMENT;
		block* new_blocks = (block*)realloc(blocks, blocks_size * sizeof(block));
		if (new_blocks == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory reallocation.");
		}
		blocks = new_blocks;
	}

	blocks[n_blocks].start = start;
	blocks[n_blocks].end = end;

	++n_blocks;
}

unsigned int Partition::get_n_blocks() {
	return n_blocks;
}

void Partition::write(const char* output_file_name) throw (Exception) {
	Writer* writer = NULL;

	block block;

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
		writer = WriterFactory::create(Writer::TEXT);
		writer->set_file_name(output_file_name);
		writer->open(true);

		writer->write("# PHASE FILE: %s\n", db->hap_file_name);
		writer->write("# MAP FILE: %s\n", db->map_file_name == NULL ? "NA" : db->map_file_name);
		if ((db->start_position > 0u) || (db->end_position != numeric_limits<unsigned long int>::max())) {
			writer->write("# REGION: [%u, %u]\n", db->start_position, db->end_position);
		} else {
			writer->write("# REGION: NA\n");
		}
		writer->write("# MAF FILTER: > %g\n", db->maf_threshold);
		writer->write("# ALL SNPs: %u\n", db->n_unfiltered_markers);
		writer->write("# FILTERED SNPs: %u\n",db->n_markers);
		writer->write("# HAPLOTYPES: %u\n", db->n_haplotypes);
		writer->write("# D' CI COMPUTATION METHOD: %s\n", ci_method);
		if (likelihood_density > 0u) {
			writer->write("# D' LIKELIHOOD DENSITY: %u\n", likelihood_density);
		} else {
			writer->write("# D' LIKELIHOOD DENSITY: NA\n");
		}
		writer->write("# D' CI LOWER BOUND FOR STRONG LD: >= %g\n", strong_pair_cl);
		writer->write("# D' CI UPPER BOUND FOR STRONG LD: >= %g\n", strong_pair_cu);
		writer->write("# D' CI UPPER BOUND FOR RECOMBINATION: <= %g\n", recomb_pair_cu);
		writer->write("# FRACTION OF STRONG LD SNP PAIRS: >= %g\n", strong_pairs_fraction);
		writer->write("# PRUNING METHOD: %s\n", pruning_method);
		if (window > 0u) {
			writer->write("# WINDOW: %ld\n", window);
		} else {
			writer->write("# WINDOW: NA\n", window);
		}

		writer->write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				"BLOCK_NAME", "FIRST_SNP", "LAST_SNP", "FIRST_SNP_ID", "LAST_SNP_ID", "START_BP", "END_BP", "N_SNPS", "N_HAPS", "N_UNIQUE_HAPS", "N_COMMON_HAPS", "HAPS_DIVERSITY");

		for (unsigned int b = 0u; b < n_blocks; ++b) {
			block = blocks[b];

			first_marker = db->markers[block.start];
			last_marker = db->markers[block.end];
			start_bp = db->positions[block.start];
			end_bp = db->positions[block.end];
			n_markers = block.end - block.start + 1u;

			get_block_diversity(b, &n_haps, &n_unique_haps, &n_common_haps, &haps_diversity);

			writer->write("BLOCK_%07u\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%g\n",
					b + 1u, first_marker, last_marker, block.start, block.end, start_bp, end_bp, n_markers, n_haps, n_unique_haps, n_common_haps, haps_diversity);
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

double Partition::get_memory_usage() {
	return ((blocks_size * sizeof(block)) / 1048576.0);
}

bool Partition::is_compatible_haplotype(const char* first, const char* second) {
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

void Partition::get_block_diversity(unsigned int block_id, unsigned int* n_haps, unsigned int* n_unique_haps, unsigned int* n_common_haps, double* haps_diversity) throw (Exception) {
	block block;

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

		block = blocks[block_id];

		n_markers = block.end - block.start + 1u;

		hap = (char*)malloc((n_markers + 1u) * sizeof(char));
		if (hap == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}
		hap[n_markers] = '\0';

		/* Read all (ambiguous & unambiguous) haplotypes */
		for (unsigned int j = 0u; j < db->n_haplotypes; ++j) {
			for (unsigned int i = block.start, k = 0u; i <= block.end; ++i, ++k) {
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
		for (haps_it = haps.begin(); haps_it != haps.end(); ++haps_it) {
			unambiguous_haps.insert(std::pair<const char*, unsigned int>(haps_it->first, haps_it->second));
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
