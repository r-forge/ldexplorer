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

#include "include/Db.h"

const char* Db::VCF = "VCF";
const char* Db::HAPMAP2 = "HAPMAP2";

const char Db::VCF_FIELD_SEPARATOR = '\t';
const char* Db::VCF_FILE_FORMAT = "##fileformat";
const char* Db::VCF_CHROM = "#CHROM";
const char* Db::VCF_POS = "POS";
const char* Db::VCF_ID = "ID";
const char* Db::VCF_REF = "REF";
const char* Db::VCF_ALT = "ALT";
const char* Db::VCF_QUAL = "QUAL";
const char* Db::VCF_FILTER = "FILTER";
const char* Db::VCF_INFO = "INFO";
const char* Db::VCF_FORMAT = "FORMAT";
const char Db::VCF_INFO_FIELD_SEPARATOR = ';';
const char* Db::VCF_VARIANT_TYPE = "VT";
const char* Db::VCF_SNP_TYPE = "SNP";
const unsigned int Db::VCF_MANDATORY_COLUMNS_SIZE = 9u;
const char* Db::vcf_mandatory_columns[VCF_MANDATORY_COLUMNS_SIZE] = {
		VCF_CHROM, VCF_POS, VCF_ID, VCF_REF, VCF_ALT, VCF_QUAL, VCF_FILTER, VCF_INFO, VCF_FORMAT
};

const char Db::HAPMAP2_MAP_FIELD_SEPARATOR = '\t';
const char Db::HAPMAP2_HAP_FIELD_SEPARATOR = ' ';
const char* Db::HAPMAP2_MAP_RS = "rs";
const char* Db::HAPMAP2_MAP_POSITION = "position";
const char* Db::HAPMAP2_MAP_0 = "0";
const char* Db::HAPMAP2_MAP_1 = "1";
const unsigned int Db::HAPMAP2_MAP_MANDATORY_COLUMNS_SIZE = 4u;
const char* Db::hapmap2_map_mandatory_columns[HAPMAP2_MAP_MANDATORY_COLUMNS_SIZE] = {
		HAPMAP2_MAP_RS, HAPMAP2_MAP_POSITION, HAPMAP2_MAP_0, HAPMAP2_MAP_1
};

const unsigned int Db::HEAP_SIZE = 2000000;
const unsigned int Db::HEAP_INCREMENT = 100000;

const double Db::EPSILON = 0.000000001;

Db::Db() throw (Exception): n_haplotypes(0u), all_n_markers(0u), all_markers(NULL), all_positions(NULL),
		all_major_alleles(NULL), all_minor_alleles(), all_major_allele_freqs(NULL), all_haplotypes(NULL),
		current_heap_size(HEAP_SIZE) {

	all_markers = (char**)malloc(current_heap_size * sizeof(char*));
	if (all_markers == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < current_heap_size; ++i) {
		all_markers[i] = NULL;
	}

	all_positions = (unsigned long int*)malloc(current_heap_size * sizeof(unsigned long int));
	if (all_positions == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	all_major_alleles = (char*)malloc(current_heap_size * sizeof(char));
	if (all_major_alleles == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	all_minor_alleles = (char*)malloc(current_heap_size * sizeof(char));
	if (all_minor_alleles == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < current_heap_size; ++i) {
		all_major_alleles[i] = '\0';
		all_minor_alleles[i] = '\0';
	}

	all_major_allele_freqs = (double*)malloc(current_heap_size * sizeof(double));
	if (all_major_allele_freqs == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	all_haplotypes = (char**)malloc(current_heap_size * sizeof(char*));
	if (all_haplotypes == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	for (unsigned int i = 0u; i < current_heap_size; ++i) {
		all_haplotypes[i] = NULL;
	}
}

Db::~Db() {
	vector<DbView*>::iterator views_it;

	for (views_it = views.begin(); views_it != views.end(); ++views_it) {
		delete *views_it;
	}
	views.clear();

	free_markers(current_heap_size);
	free_positions(current_heap_size);
	free_alleles(current_heap_size);
	free_major_allele_freqs(current_heap_size);
	free_haplotypes(current_heap_size);
}

void Db::free_markers(unsigned int heap_size) {
	if (all_markers != NULL) {
		for (unsigned int i = 0u; i < heap_size; ++i) {
			if (all_markers[i] != NULL) {
				free(all_markers[i]);
				all_markers[i] = NULL;
			}
		}

		free(all_markers);
		all_markers = NULL;
	}
}

void Db::free_positions(unsigned int heap_size) {
	if (all_positions != NULL) {
		free(all_positions);
		all_positions = NULL;
	}
}

void Db::free_alleles(unsigned int heap_size) {
	if (all_major_alleles != NULL) {
		free(all_major_alleles);
		all_major_alleles = NULL;
	}

	if (all_minor_alleles != NULL) {
		free(all_minor_alleles);
		all_minor_alleles = NULL;
	}
}

void Db::free_major_allele_freqs(unsigned int heap_size) {
	if (all_major_allele_freqs != NULL) {
		free(all_major_allele_freqs);
		all_major_allele_freqs = NULL;
	}
}

void Db::free_haplotypes(unsigned int heap_size) {
	if (all_haplotypes != NULL) {
		for (unsigned int i = 0u; i < heap_size; ++i) {
			if (all_haplotypes[i] != NULL) {
				free(all_haplotypes[i]);
				all_haplotypes[i] = NULL;
			}
		}

		free(all_haplotypes);
		all_haplotypes = NULL;
	}
}

void Db::reallocate() throw (Exception) {
	char** new_all_markers = NULL;
	unsigned long int* new_all_positions = NULL;
	char* new_all_major_alleles = NULL;
	char* new_all_minor_alleles = NULL;
	double* new_all_major_allele_freqs = NULL;
	char** new_all_haplotypes = NULL;
	unsigned int new_heap_size = current_heap_size + HEAP_INCREMENT;

	new_all_markers = (char**)realloc(all_markers, new_heap_size * sizeof(char*));
	if (new_all_markers == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
	}
	all_markers = new_all_markers;
	new_all_markers = NULL;

	for (unsigned int i = current_heap_size; i < new_heap_size; ++i) {
		all_markers[i] = NULL;
	}

	new_all_positions = (unsigned long int*)realloc(all_positions, new_heap_size * sizeof(unsigned long int));
	if (new_all_positions == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
	}
	all_positions = new_all_positions;
	new_all_positions = NULL;

	new_all_major_alleles = (char*)realloc(all_major_alleles, new_heap_size * sizeof(char));
	if (new_all_major_alleles == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
	}
	all_major_alleles = new_all_major_alleles;
	new_all_major_alleles = NULL;

	new_all_minor_alleles = (char*)realloc(all_minor_alleles, new_heap_size * sizeof(char));
	if (new_all_minor_alleles == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
	}
	all_minor_alleles = new_all_minor_alleles;
	new_all_minor_alleles = NULL;

	for (unsigned int i = current_heap_size; i < new_heap_size; ++i) {
		all_major_alleles[i] = '\0';
		all_minor_alleles[i] = '\0';
	}

	new_all_major_allele_freqs = (double*)realloc(all_major_allele_freqs, new_heap_size * sizeof(double));
	if (new_all_major_allele_freqs == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
	}
	all_major_allele_freqs = new_all_major_allele_freqs;
	new_all_major_allele_freqs = NULL;

	new_all_haplotypes = (char**)realloc(all_haplotypes, new_heap_size * sizeof(char*));
	if (new_all_haplotypes == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
	}
	all_haplotypes = new_all_haplotypes;
	new_all_haplotypes = NULL;

	for (unsigned int i = current_heap_size; i < new_heap_size; ++i) {
		all_haplotypes[i] = NULL;
	}

	current_heap_size = new_heap_size;
}

void Db::load(const char* hap_file_name, const char* map_file_name, unsigned long int start_position, unsigned long int end_position, const char* type) throw (Exception) {
	if (auxiliary::strcmp_ignore_case(type, VCF) == 0) {
		if ((start_position != 0u) || (end_position != numeric_limits<unsigned long int>::max())) {
			load_vcf(hap_file_name, start_position, end_position);
		} else {
			load_vcf(hap_file_name);
		}
	} else if (auxiliary::strcmp_ignore_case(type, HAPMAP2) == 0) {
		if ((start_position != 0u) || (end_position != numeric_limits<unsigned long int>::max())) {
			load_hapmap2(map_file_name, hap_file_name, start_position, end_position);
		} else {
			load_hapmap2(map_file_name, hap_file_name);
		}
	} else {
		throw Exception(__FILE__, __LINE__, "Unknown file type '%s' was specified.", type);
	}
}

void Db::load_vcf(const char* file_name, unsigned long int start_position, unsigned long int end_position) throw (Exception) {
	Reader* reader = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 0u;

	char* token = NULL;
	char** tokens = NULL;

	unsigned int total_column_number = 0u;
	unsigned int column_number = 0u;
	int sample_number = 0;

	Unique chromosome;
	unsigned int length = 0u;

	unsigned int first_allele_index = 0u;
	unsigned int second_allele_index = 0u;

	unsigned int n_ref_allele = 0u;
	unsigned int n_alt_allele = 0u;
	char swap_allele = '\0';

	int vt_length = strlen(VCF_VARIANT_TYPE);

	try {
		reader = ReaderFactory::create(file_name);
		reader->open();

		/* Read the first required line with file format description. */
		if ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			if ((token = auxiliary::strtok(&line, '=')) != NULL) {
				auxiliary::trim_end(token);
				if (auxiliary::strcmp_ignore_case(token, VCF_FILE_FORMAT) != 0) {
					throw Exception(__FILE__, __LINE__, "The mandatory VCF file format information line in '%s' file is incorrect.", file_name);
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The mandatory VCF file format information line in '%s' file is incorrect.", file_name);
			}
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The mandatory VCF file format information line in '%s' file is empty.", file_name);
		}

		/* Read the mandatory header. Meta-info lines are optional. */
		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);

			if (line_length > 1) {
				if (line[0u] != '#') {
					throw Exception(__FILE__, __LINE__, "The mandatory VCF header line in '%s' file was not found.", file_name);
				}

				if (line[1u] != '#') {
					while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
						if (total_column_number < VCF_MANDATORY_COLUMNS_SIZE) {
							if (auxiliary::strcmp_ignore_case(token, vcf_mandatory_columns[total_column_number]) != 0) {
								throw Exception(__FILE__, __LINE__, "Column '%s' is missing on position %d in '%s' file.", vcf_mandatory_columns[total_column_number], total_column_number + 1u, file_name);
							}
						} else {
							/* sample columns */
						}
						++total_column_number;
					}
					break;
				} else {
					/* process meta-info line if necessary */
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The header/meta-information on line %d in '%s' file is incorrect.", line_number, file_name);
			}
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, file_name);
		}

		if ((sample_number  = total_column_number - VCF_MANDATORY_COLUMNS_SIZE) <= 0) {
			throw Exception(__FILE__, __LINE__, "No sample columns were found in '%s' file.", file_name);
		}

		n_haplotypes = 2u * ((unsigned int)sample_number);

		/* Read data. */
		tokens = (char**)malloc(total_column_number * sizeof(char*));
		if (tokens == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			column_number = 0u;

			while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
				if (column_number < total_column_number) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number != total_column_number) {
				throw Exception(__FILE__, __LINE__, "The number of columns (%d) on line %d in '%s' file is not equal to the expected (%d).", column_number, line_number, file_name, total_column_number);
			}

			/* tokens[0] -- chromosome. check if unique accross all VCF file. must be one file per chromosome. */
			if (!(chromosome.*chromosome.check)(tokens[0u])) {
				throw Exception(__FILE__, __LINE__, "Unexpected chromosome '%s' (expected chromosome is '%s') was found on line %d in '%s' file.", tokens[0u], chromosome.get_value(), line_number, file_name);
			}

			if (all_n_markers >= current_heap_size) {
				reallocate();
			}

			/* tokens[1] -- position. parse to unsigned long integer */
			if (!auxiliary::to_ulong_int(tokens[1u], &(all_positions[all_n_markers]))) {
				throw Exception(__FILE__, __LINE__, "The chromosomal position '%s' on line %d in '%s' file could not be parsed to unsigned integer.", tokens[1u], line_number, file_name);
			}

			/* tokens[7] -- check info field. if variant type is specified, then it must be SNP. */
			while ((token = auxiliary::strtok(&(tokens[7u]), VCF_INFO_FIELD_SEPARATOR)) != NULL) {
				auxiliary::trim_start(token);
				if (auxiliary::strcmp_ignore_case(token, VCF_VARIANT_TYPE, vt_length) == 0) {
					token = strchr(token, '=');
					if (token != NULL) {
						++token;

						auxiliary::trim_start(token);
						auxiliary::trim_end(token);

						if (auxiliary::strcmp_ignore_case(token, VCF_SNP_TYPE) != 0) {
							continue;
						}
					}
					break;
				}
			}

			if ((all_positions[all_n_markers] >= start_position) && (all_positions[all_n_markers] <= end_position)) {
				/* tokens[4] -- alt allele, if more than one value, then skip. if one value and more than one letter -- skip. Otherwise check if in {A,C,G,T} */
				if ((length = strlen(tokens[4u])) < 1u) {
					throw Exception(__FILE__, __LINE__, "The alternate allele value on line %d in '%s' file is empty.", line_number, file_name);
				} else if (length == 1u) {
					all_minor_alleles[all_n_markers] = toupper(tokens[4u][0u]);
					if (all_minor_alleles[all_n_markers] == '.') {
						/* monomorphic SNP, not interesting for haplotypes */
						continue;
					} else if ((all_minor_alleles[all_n_markers] != 'A') && (all_minor_alleles[all_n_markers] != 'C') &&
							(all_minor_alleles[all_n_markers] != 'G') && (all_minor_alleles[all_n_markers] != 'T')) {
						throw Exception(__FILE__, __LINE__, "The alternate allele value '%s' on line %d in '%s' file is incorrect.", tokens[4u], line_number, file_name);
					}
				} else {
					/* multi-allelic SNP or indel, deletion and etc. */
					continue;
				}

				/* tokens[3] -- ref allele, if more than one letter, then indel -- skip. Otherwise check if in {A,C,G,T}. */
				if ((length = strlen(tokens[3u])) < 1u) {
					throw Exception(__FILE__, __LINE__, "The reference allele value on line %d in '%s' file is empty.", line_number, file_name);
				} else if (length == 1u) {
					all_major_alleles[all_n_markers] = toupper(tokens[3u][0u]);
					if ((all_major_alleles[all_n_markers] != 'A') && (all_major_alleles[all_n_markers] != 'C') &&
							(all_major_alleles[all_n_markers] != 'G') && (all_major_alleles[all_n_markers] != 'T')) {
						throw Exception(__FILE__, __LINE__, "The reference allele value '%s' on line %d in '%s' file is incorrect.", tokens[3u], line_number, file_name);
					}
				} else {
					/* indel, deletion and etc. */
					continue;
				}

				/* tokens[2] -- rsId. Copy it. */
				all_markers[all_n_markers] = (char*)malloc((strlen(tokens[2u]) + 1u) * sizeof(char));
				if (all_markers[all_n_markers] == NULL) {
					throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
				}
				strcpy(all_markers[all_n_markers], tokens[2u]);

				/* tokens[i] with i = 9,...,N are samples */
				all_haplotypes[all_n_markers] = (char*)malloc(n_haplotypes * sizeof(char));
				if (all_haplotypes[all_n_markers] == NULL) {
					throw Exception(__FILE__, __LINE__, "Error in memory allocation");
				}

				n_ref_allele = 0u;
				n_alt_allele = 0u;

				for (unsigned int i = 9u; i < total_column_number; ++i) {
					if ((token = auxiliary::strtok(&(tokens[i]), ':')) != NULL) {
						if (strlen(token) != 3u) {
							throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has an incorrect genotype value '%s'.", i - 9u, line_number, file_name, token);
						}

						if (token[1u] != '|') {
							throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has UNPHASED genotype '%s'.", i - 9u, line_number, file_name, token);
						}

						first_allele_index = 2u * (i - 9u);
						second_allele_index = 2u * (i - 9u) + 1u;

						if ((token[0u] == '.') && (token[2u] == '.')) {
							all_haplotypes[all_n_markers][first_allele_index] = '.';
							all_haplotypes[all_n_markers][second_allele_index] = '.';
						}

						if (token[0u] == '0') {
							all_haplotypes[all_n_markers][first_allele_index] = all_major_alleles[all_n_markers];
							++n_ref_allele;
						} else if (token[0u] == '1') {
							all_haplotypes[all_n_markers][first_allele_index] = all_minor_alleles[all_n_markers];
							++n_alt_allele;
						} else {
							throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has unexpected first allele '%c'.", i - 9, line_number, file_name, token[0u]);
						}

						if (token[2u] == '0') {
							all_haplotypes[all_n_markers][second_allele_index] = all_major_alleles[all_n_markers];
							++n_ref_allele;
						} else if (token[2u] == '1') {
							all_haplotypes[all_n_markers][second_allele_index] = all_minor_alleles[all_n_markers];
							++n_alt_allele;
						} else {
							throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has unexpected second allele '%c'.", i - 9u, line_number, file_name, token[2u]);
						}
					} else {
						throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has an incorrect value.", i - 9u, line_number, file_name);
					}
				}

				/* Save major allele frequency. Swap. First allele in alleles array must be major. */
				if (n_ref_allele < n_alt_allele) {
					swap_allele = all_major_alleles[all_n_markers];
					all_major_alleles[all_n_markers] = all_minor_alleles[all_n_markers];
					all_minor_alleles[all_n_markers] = swap_allele;
					all_major_allele_freqs[all_n_markers] = ((double)n_alt_allele) / ((double)(n_ref_allele + n_alt_allele));
				} else {
					all_major_allele_freqs[all_n_markers] = ((double)n_ref_allele) / ((double)(n_ref_allele + n_alt_allele));
				}

				++all_n_markers;
			}
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, file_name);
		}

		reader->close();
		delete reader;

		free(tokens);
		tokens = NULL;

//		mask(numeric_limits<double>::quiet_NaN());
	} catch (Exception &e) {
		e.add_message(__FILE__, __LINE__, "Error while loading '%s' file.", file_name);
		throw;
	}

}

void Db::load_vcf(const char* file_name) throw (Exception) {
	Reader* reader = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 0u;
	char* token = NULL;
	char** tokens = NULL;

	unsigned int total_column_number = 0u;
	unsigned int column_number = 0u;
	int sample_number = 0;

	Unique chromosome;
	unsigned int length = 0u;

	unsigned int first_allele_index = 0u;
	unsigned int second_allele_index = 0u;

	unsigned int n_ref_allele = 0u;
	unsigned int n_alt_allele = 0u;
	char swap_allele = '\0';

	int vt_length = strlen(VCF_VARIANT_TYPE);

	try {
		reader = ReaderFactory::create(file_name);
		reader->open();

		/* Read the first required line with file format description. */
		if ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			if ((token = auxiliary::strtok(&line, '=')) != NULL) {
				auxiliary::trim_end(token);
				if (auxiliary::strcmp_ignore_case(token, VCF_FILE_FORMAT) != 0) {
					throw Exception(__FILE__, __LINE__, "The mandatory VCF file format information line in '%s' file is incorrect.", file_name);
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The mandatory VCF file format information line in '%s' file is incorrect.", file_name);
			}
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The mandatory VCF file format information line in '%s' file is empty.", file_name);
		}

		/* Read the mandatory header. Meta-info lines are optional. */
		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);

			if (line_length > 1) {
				if (line[0u] != '#') {
					throw Exception(__FILE__, __LINE__, "The mandatory VCF header line in '%s' file was not found.", file_name);
				}

				if (line[1u] != '#') {
					while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
						if (total_column_number < VCF_MANDATORY_COLUMNS_SIZE) {
							if (auxiliary::strcmp_ignore_case(token, vcf_mandatory_columns[total_column_number]) != 0) {
								throw Exception(__FILE__, __LINE__, "Column '%s' is missing on position %d in '%s' file.", vcf_mandatory_columns[total_column_number], total_column_number + 1u, file_name);
							}
						} else {
							/* sample columns */
						}
						++total_column_number;
					}
					break;
				} else {
					/* process meta-info line if necessary */
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The header/meta-information on line %d in '%s' file is incorrect.", line_number, file_name);
			}
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, file_name);
		}

		if ((sample_number  = total_column_number - VCF_MANDATORY_COLUMNS_SIZE) <= 0) {
			throw Exception(__FILE__, __LINE__, "No sample columns were found in '%s' file.", file_name);
		}

		n_haplotypes = 2u * ((unsigned int)sample_number);

		/* Read data. */
		tokens = (char**)malloc(total_column_number * sizeof(char*));
		if (tokens == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			column_number = 0u;

			while ((token = auxiliary::strtok(&line, VCF_FIELD_SEPARATOR)) != NULL) {
				if (column_number < total_column_number) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number != total_column_number) {
				throw Exception(__FILE__, __LINE__, "The number of columns (%d) on line %d in '%s' file is not equal to the expected (%d).", column_number, line_number, file_name, total_column_number);
			}

			/* tokens[0] -- chromosome. check if unique accross all VCF file. must be one file per chromosome. */
			if (!(chromosome.*chromosome.check)(tokens[0u])) {
				throw Exception(__FILE__, __LINE__, "Unexpected chromosome '%s' (expected chromosome is '%s') was found on line %d in '%s' file.", tokens[0u], chromosome.get_value(), line_number, file_name);
			}

			if (all_n_markers >= current_heap_size) {
				reallocate();
			}

			/* tokens[1] -- position. parse to unsigned long integer */
			if (!auxiliary::to_ulong_int(tokens[1u], &(all_positions[all_n_markers]))) {
				throw Exception(__FILE__, __LINE__, "The chromosomal position '%s' on line %d in '%s' file could not be parsed to unsigned integer.", tokens[1u], line_number, file_name);
			}

			/* tokens[7] -- check info field. if variant type is specified, then it must be SNP. */
			while ((token = auxiliary::strtok(&(tokens[7u]), VCF_INFO_FIELD_SEPARATOR)) != NULL) {
				auxiliary::trim_start(token);
				if (auxiliary::strcmp_ignore_case(token, VCF_VARIANT_TYPE, vt_length) == 0) {
					token = strchr(token, '=');
					if (token != NULL) {
						++token;

						auxiliary::trim_start(token);
						auxiliary::trim_end(token);

						if (auxiliary::strcmp_ignore_case(token, VCF_SNP_TYPE) != 0) {
							continue;
						}
					}
					break;
				}
			}

			/* tokens[4] -- alt allele, if more than one value, then skip. if one value and more than one letter -- skip. Otherwise check if in {A,C,G,T} */
			if ((length = strlen(tokens[4u])) < 1u) {
				throw Exception(__FILE__, __LINE__, "The alternate allele value on line %d in '%s' file is empty.", line_number, file_name);
			} else if (length == 1u) {
				all_minor_alleles[all_n_markers] = toupper(tokens[4u][0u]);
				if (all_minor_alleles[all_n_markers] == '.') {
					/* monomorphic SNP, not interesting for haplotypes */
					continue;
				} else if ((all_minor_alleles[all_n_markers] != 'A') && (all_minor_alleles[all_n_markers] != 'C') &&
						(all_minor_alleles[all_n_markers] != 'G') && (all_minor_alleles[all_n_markers] != 'T')) {
					throw Exception(__FILE__, __LINE__, "The alternate allele value '%s' on line %d in '%s' file is incorrect.", tokens[4u], line_number, file_name);
				}
			} else {
				/* multi-allelic SNP or indel, deletion and etc. */
				continue;
			}

			/* tokens[3] -- ref allele, if more than one letter, then indel -- skip. Otherwise check if in {A,C,G,T}. */
			if ((length = strlen(tokens[3u])) < 1u) {
				throw Exception(__FILE__, __LINE__, "The reference allele value on line %d in '%s' file is empty.", line_number, file_name);
			} else if (length == 1u) {
				all_major_alleles[all_n_markers] = toupper(tokens[3u][0u]);
				if ((all_major_alleles[all_n_markers] != 'A') && (all_major_alleles[all_n_markers] != 'C') &&
						(all_major_alleles[all_n_markers] != 'G') && (all_major_alleles[all_n_markers] != 'T')) {
					throw Exception(__FILE__, __LINE__, "The reference allele value '%s' on line %d in '%s' file is incorrect.", tokens[3u], line_number, file_name);
				}
			} else {
				/* indel, deletion and etc. */
				continue;
			}

			/* tokens[2] -- rsId. Copy it. */
			all_markers[all_n_markers] = (char*)malloc((strlen(tokens[2u]) + 1u) * sizeof(char));
			if (all_markers[all_n_markers] == NULL) {
				throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
			}
			strcpy(all_markers[all_n_markers], tokens[2u]);

			/* tokens[i] with i = 9,...,N are samples */
			all_haplotypes[all_n_markers] = (char*)malloc(n_haplotypes * sizeof(char));
			if (all_haplotypes[all_n_markers] == NULL) {
				throw Exception(__FILE__, __LINE__, "Error in memory allocation");
			}

			n_ref_allele = 0u;
			n_alt_allele = 0u;

			for (unsigned int i = 9u; i < total_column_number; ++i) {
				if ((token = auxiliary::strtok(&(tokens[i]), ':')) != NULL) {
					if (strlen(token) != 3u) {
						throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has an incorrect genotype value '%s'.", i - 9u, line_number, file_name, token);
					}

					if (token[1u] != '|') {
						throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has UNPHASED genotype '%s'.", i - 9u, line_number, file_name, token);
					}

					first_allele_index = 2u * (i - 9u);
					second_allele_index = 2u * (i - 9u) + 1u;

					if ((token[0u] == '.') && (token[2u] == '.')) {
						all_haplotypes[all_n_markers][first_allele_index] = '.';
						all_haplotypes[all_n_markers][second_allele_index] = '.';
					}

					if (token[0u] == '0') {
						all_haplotypes[all_n_markers][first_allele_index] = all_major_alleles[all_n_markers];
						++n_ref_allele;
					} else if (token[0u] == '1') {
						all_haplotypes[all_n_markers][first_allele_index] = all_minor_alleles[all_n_markers];
						++n_alt_allele;
					} else {
						throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has unexpected first allele '%c'.", i - 9u, line_number, file_name, token[0u]);
					}

					if (token[2u] == '0') {
						all_haplotypes[all_n_markers][second_allele_index] = all_major_alleles[all_n_markers];
						++n_ref_allele;
					} else if (token[2u] == '1') {
						all_haplotypes[all_n_markers][second_allele_index] = all_minor_alleles[all_n_markers];
						++n_alt_allele;
					} else {
						throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has unexpected second allele '%c'.", i - 9u, line_number, file_name, token[2u]);
					}
				} else {
					throw Exception(__FILE__, __LINE__, "Sample %d on line %d in '%s' file has an incorrect value.", i - 9u, line_number, file_name);
				}
			}

			/* Save major allele frequency. Swap. First allele in alleles array must be major. */
			if (n_ref_allele < n_alt_allele) {
				swap_allele = all_major_alleles[all_n_markers];
				all_major_alleles[all_n_markers] = all_minor_alleles[all_n_markers];
				all_minor_alleles[all_n_markers] = swap_allele;
				all_major_allele_freqs[all_n_markers] = ((double)n_alt_allele) / ((double)(n_ref_allele + n_alt_allele));
			} else {
				all_major_allele_freqs[all_n_markers] = ((double)n_ref_allele) / ((double)(n_ref_allele + n_alt_allele));
			}

			++all_n_markers;
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, file_name);
		}

		reader->close();
		delete reader;

		free(tokens);
		tokens = NULL;
	} catch (Exception &e) {
		e.add_message(__FILE__, __LINE__, "Error while loading '%s' file.", file_name);
		throw;
	}
}

void Db::load_hapmap2(const char* map_file_name, const char* hap_file_name, unsigned long int start_position, unsigned long int end_position) throw (Exception) {
	Reader* reader = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 0u;

	char* token = NULL;
	char** tokens = NULL;

	unsigned int total_column_number = 0u;
	unsigned int column_number = 0u;

	unsigned char* filter = NULL;
	unsigned char* new_filter = NULL;
	unsigned int current_filter_size = HEAP_SIZE;
	unsigned int filter_n_markers = 0u;

	unsigned long int position = 0u;
	char first_allele = '\0';
	char second_allele = '\0';

	char swap_allele = '\0';
	unsigned int* n_first_alleles = NULL;
	unsigned int* n_second_alleles = NULL;
	char* extended_haplotype = NULL;

	filter = (unsigned char*)malloc(current_filter_size * sizeof(unsigned char));
	if (filter == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	/* Read map file. */
	try {
		reader = ReaderFactory::create(map_file_name);
		reader->open();

		/* Read map header. */
		if ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);

			while ((token = auxiliary::strtok(&line, HAPMAP2_MAP_FIELD_SEPARATOR)) != NULL) {
				if (total_column_number < HAPMAP2_MAP_MANDATORY_COLUMNS_SIZE) {
					if (auxiliary::strcmp_ignore_case(token, hapmap2_map_mandatory_columns[total_column_number]) != 0) {
						throw Exception(__FILE__, __LINE__, "Column '%s' is missing on position %d in '%s' file.", hapmap2_map_mandatory_columns[total_column_number], total_column_number + 1u, map_file_name);
					}
				} else {
					/* other columns */
				}
				++total_column_number;
			}
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The mandatory header line in '%s' HAPMAP2 legend file is empty.", map_file_name);
		}

		/* Read map data. */
		tokens = (char**)malloc(total_column_number * sizeof(char*));
		if (tokens == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			column_number = 0u;

			while ((token = auxiliary::strtok(&line, HAPMAP2_MAP_FIELD_SEPARATOR)) != NULL) {
				if (column_number < total_column_number) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number != total_column_number) {
				throw Exception(__FILE__, __LINE__, "The number of columns (%d) on line %d in '%s' file is not equal to the expected (%d).", column_number, line_number, map_file_name, total_column_number);
			}

			if (!auxiliary::to_ulong_int(tokens[1u], &position)) {
				throw Exception(__FILE__, __LINE__, "The chromosomal position '%s' on line %d in '%s' file could not be parsed to unsigned integer.", tokens[1u], line_number, map_file_name);
			}

			if (strlen(tokens[2u]) == 1u) {
				first_allele = toupper(tokens[2u][0u]);
				if ((first_allele != 'A') && (first_allele != 'C') && (first_allele != 'G') && (first_allele != 'T')) {
					throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d  in '%s' file is incorrect.", tokens[2u], line_number, map_file_name);
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d in '%s' file is incorrect.", tokens[2u], line_number, map_file_name);
			}

			if (strlen(tokens[3u]) == 1u) {
				second_allele = toupper(tokens[3u][0u]);
				if ((second_allele != 'A') && (second_allele != 'C') && (second_allele != 'G') && (second_allele != 'T')) {
					throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d in '%s' file is incorrect.", tokens[3u], line_number, map_file_name);
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d in '%s' file is incorrect.", tokens[3u], line_number, map_file_name);
			}

			if (filter_n_markers >= current_filter_size) {
				current_filter_size += HEAP_INCREMENT;
				new_filter = (unsigned char*)realloc(filter, current_filter_size * sizeof(unsigned char));
				if (new_filter == NULL) {
					free(filter);
					filter = NULL;
					throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
				}
				filter = new_filter;
				new_filter = NULL;
			}

			if ((position >= start_position) && (position <= end_position)) {
				filter[filter_n_markers] = 0u;

				if (all_n_markers >= current_heap_size) {
					reallocate();
				}

				all_markers[all_n_markers] = (char*)malloc((strlen(tokens[0u]) + 1u) * sizeof(char));
				if (all_markers[all_n_markers] == NULL) {
					throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
				}
				strcpy(all_markers[all_n_markers], tokens[0u]);

				all_positions[all_n_markers] = position;
				all_major_alleles[all_n_markers] = first_allele;
				all_minor_alleles[all_n_markers] = second_allele;

				++all_n_markers;
			} else {
				filter[filter_n_markers] = 1u;
			}

			++filter_n_markers;
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, map_file_name);
		}

		reader->close();
		delete reader;

		free(tokens);
	} catch (Exception &e) {
		e.add_message(__FILE__, __LINE__, "Error while loading '%s' file.", map_file_name);
		throw;
	}

	line_number = 0u;

	/* Read haplotype file. */
	try {
		reader = ReaderFactory::create(hap_file_name);
		reader->open();

		tokens = (char**)malloc(filter_n_markers * sizeof(char*));
		if (tokens == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		n_first_alleles = (unsigned int*)malloc(all_n_markers * sizeof(unsigned int));
		if (n_first_alleles == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		n_second_alleles = (unsigned int*)malloc(all_n_markers * sizeof(unsigned int));
		if (n_second_alleles == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		for (unsigned int i = 0u; i < all_n_markers; ++i) {
			n_first_alleles[i] = 0u;
			n_second_alleles[i] = 0u;
		}

		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			column_number = 0u;

			/* eliminate extra field separator at the end of line */
			if (line[line_length - 1] == HAPMAP2_HAP_FIELD_SEPARATOR) {
				line[line_length - 1] = '\0';
			}

			while ((token = auxiliary::strtok(&line, HAPMAP2_HAP_FIELD_SEPARATOR)) != NULL) {
				if (column_number < filter_n_markers) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number != filter_n_markers) {
				throw Exception(__FILE__, __LINE__, "The number of columns (%d) on line %d in '%s' file is not equal to the expected (%d).", column_number, line_number, hap_file_name, filter_n_markers);
			}

			for (unsigned int i = 0u, j = 0u; i < filter_n_markers; ++i) {
				if (strlen(tokens[i]) != 1u) {
					throw Exception(__FILE__, __LINE__, "Sample on line %d in '%s' file has incorrect allele value '%s' for marker on position %d.", line_number, hap_file_name, tokens[i], i + 1u);
				} else if ((tokens[i][0u] != '0') && (tokens[i][0u] != '1')) {
					throw Exception(__FILE__, __LINE__, "Sample on line %d in '%s' file has unexpected allele value '%c' for marker on position %d.", line_number, hap_file_name, tokens[i][0u], i + 1u);
				}

				if (filter[i] == 0u) {
					extended_haplotype = (char*)realloc(all_haplotypes[j], (n_haplotypes + 1u) * sizeof(char));
					if (extended_haplotype == NULL) {
						throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
					}
					all_haplotypes[j] = extended_haplotype;
					extended_haplotype = NULL;

					if (tokens[i][0u] == '0') {
						++n_first_alleles[j];
						all_haplotypes[j][n_haplotypes] = all_major_alleles[j];
					} else if (tokens[i][0u] == '1') {
						++n_second_alleles[j];
						all_haplotypes[j][n_haplotypes] = all_minor_alleles[j];
					}

					++j;
				}
			}

			++n_haplotypes;
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, hap_file_name);
		}

		reader->close();
		delete reader;

		/* Calculate major allele frequencies. */
		for (unsigned int i = 0u; i < all_n_markers; ++i) {
			if (n_first_alleles[i] < n_second_alleles[i]) {
				swap_allele = all_major_alleles[i];
				all_major_alleles[i] = all_minor_alleles[i];
				all_minor_alleles[i] = swap_allele;
				all_major_allele_freqs[i] = ((double)n_second_alleles[i]) / ((double)(n_first_alleles[i] + n_second_alleles[i]));
			} else {
				all_major_allele_freqs[i] = ((double)n_first_alleles[i]) / ((double)(n_first_alleles[i] + n_second_alleles[i]));
			}
		}

		free(tokens);
		tokens = NULL;

		free(n_first_alleles);
		n_first_alleles = NULL;

		free(n_second_alleles);
		n_second_alleles = NULL;
	} catch (Exception &e) {
		e.add_message(__FILE__, __LINE__, "Error while loading '%s' file.", hap_file_name);
		throw;
	}

	free(filter);
	filter = NULL;
}

void Db::load_hapmap2(const char* map_file_name, const char* hap_file_name) throw (Exception) {
	Reader* reader = NULL;

	char* line = NULL;
	int line_length = 0;
	unsigned int line_number = 0u;

	char* token = NULL;
	char** tokens = NULL;

	unsigned int total_column_number = 0u;
	unsigned int column_number = 0u;

	char swap_allele = '\0';
	unsigned int* n_first_alleles = NULL;
	unsigned int* n_second_alleles = NULL;
	char* extended_haplotype = NULL;

	/* Read map file. */
	try {
		reader = ReaderFactory::create(map_file_name);
		reader->open();

		/* Read header. */
		if ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);

			while ((token = auxiliary::strtok(&line, HAPMAP2_MAP_FIELD_SEPARATOR)) != NULL) {
				if (total_column_number < HAPMAP2_MAP_MANDATORY_COLUMNS_SIZE) {
					if (auxiliary::strcmp_ignore_case(token, hapmap2_map_mandatory_columns[total_column_number]) != 0) {
						throw Exception(__FILE__, __LINE__, "Column '%s' is missing on position %d in '%s' file.", hapmap2_map_mandatory_columns[total_column_number], total_column_number + 1u, map_file_name);
					}
				} else {
					/* other columns */
				}
				++total_column_number;
			}
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The mandatory header line in '%s' HAPMAP2 legend file is empty.", map_file_name);
		}

		/* Read data. */
		tokens = (char**)malloc(total_column_number * sizeof(char*));
		if (tokens == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			column_number = 0u;

			while ((token = auxiliary::strtok(&line, HAPMAP2_MAP_FIELD_SEPARATOR)) != NULL) {
				if (column_number < total_column_number) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number != total_column_number) {
				throw Exception(__FILE__, __LINE__, "The number of columns (%d) on line %d in '%s' file is not equal to the expected (%d).", column_number, line_number, map_file_name, total_column_number);
			}

			if (all_n_markers >= current_heap_size) {
				reallocate();
			}

			all_markers[all_n_markers] = (char*)malloc((strlen(tokens[0u]) + 1u) * sizeof(char));
			if (all_markers[all_n_markers] == NULL) {
				throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
			}
			strcpy(all_markers[all_n_markers], tokens[0u]);

			if (!auxiliary::to_ulong_int(tokens[1u], &(all_positions[all_n_markers]))) {
				throw Exception(__FILE__, __LINE__, "The chromosomal position '%s' on line %d in '%s' file could not be parsed to unsigned integer.", tokens[1u], line_number, map_file_name);
			}

			if (strlen(tokens[2u]) == 1u) {
				all_major_alleles[all_n_markers] = toupper(tokens[2u][0u]);
				if ((all_major_alleles[all_n_markers] != 'A') && (all_major_alleles[all_n_markers] != 'C') &&
						(all_major_alleles[all_n_markers] != 'G') && (all_major_alleles[all_n_markers] != 'T')) {
					throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d in '%s' file is incorrect.", tokens[2u], line_number, map_file_name);
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d in '%s' file is incorrect.", tokens[2u], line_number, map_file_name);
			}

			if (strlen(tokens[3u]) == 1u) {
				all_minor_alleles[all_n_markers] = toupper(tokens[3u][0u]);
				if ((all_minor_alleles[all_n_markers] != 'A') && (all_minor_alleles[all_n_markers] != 'C') &&
						(all_minor_alleles[all_n_markers] != 'G') && (all_minor_alleles[all_n_markers] != 'T')) {
					throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d in '%s' file is incorrect.", tokens[3u], line_number, map_file_name);
				}
			} else {
				throw Exception(__FILE__, __LINE__, "The allele value '%s' on line %d in '%s' file is incorrect.", tokens[3u], line_number, map_file_name);
			}

			++all_n_markers;
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, map_file_name);
		}

		reader->close();
		delete reader;

		free(tokens);
	} catch (Exception &e) {
		e.add_message(__FILE__, __LINE__, "Error while loading '%s' file.", map_file_name);
		throw;
	}

	line_number = 0u;

	/* Read haplotype file. */
	try {
		reader = ReaderFactory::create(hap_file_name);
		reader->open();

		tokens = (char**)malloc(all_n_markers * sizeof(char*));
		if (tokens == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		n_first_alleles = (unsigned int*)malloc(all_n_markers * sizeof(unsigned int));
		if (n_first_alleles == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		n_second_alleles = (unsigned int*)malloc(all_n_markers * sizeof(unsigned int));
		if (n_second_alleles == NULL) {
			throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
		}

		for (unsigned int i = 0u; i < all_n_markers; ++i) {
			n_first_alleles[i] = 0u;
			n_second_alleles[i] = 0u;
		}

		while ((line_length = reader->read_line()) > 0) {
			++line_number;
			line = *(reader->line);
			column_number = 0u;

			/* eliminate extra field separator at the end of line */
			if (line[line_length - 1] == HAPMAP2_HAP_FIELD_SEPARATOR) {
				line[line_length - 1] = '\0';
			}

			while ((token = auxiliary::strtok(&line, HAPMAP2_HAP_FIELD_SEPARATOR)) != NULL) {
				if (column_number < all_n_markers) {
					tokens[column_number] = token;
				}
				++column_number;
			}

			if (column_number != all_n_markers) {
				throw Exception(__FILE__, __LINE__, "The number of columns (%d) on line %d in '%s' file is not equal to the expected (%d).", column_number, line_number, hap_file_name, all_n_markers);
			}

			for (unsigned int i = 0u; i < all_n_markers; ++i) {
				extended_haplotype = (char*)realloc(all_haplotypes[i], (n_haplotypes + 1u) * sizeof(char));
				if (extended_haplotype == NULL) {
					throw Exception(__FILE__, __LINE__, "Error in memory re-allocation.");
				}
				all_haplotypes[i] = extended_haplotype;
				extended_haplotype = NULL;

				if (strlen(tokens[i]) != 1u) {
					throw Exception(__FILE__, __LINE__, "Sample on line %d in '%s' file has incorrect allele value '%s' for marker on position %d.", line_number, hap_file_name, tokens[i], i + 1u);
				} else if (tokens[i][0u] == '0') {
					++n_first_alleles[i];
					all_haplotypes[i][n_haplotypes] = all_major_alleles[i];
				} else if (tokens[i][0u] == '1') {
					++n_second_alleles[i];
					all_haplotypes[i][n_haplotypes] = all_minor_alleles[i];
				} else {
					throw Exception(__FILE__, __LINE__, "Sample on line %d in '%s' file has unexpected allele value '%c' for marker on position %d.", line_number, hap_file_name, tokens[i][0u], i + 1u);
				}
			}

			++n_haplotypes;
		}

		if (line_length == 0) {
			throw Exception(__FILE__, __LINE__, "The line %d in '%s' file is empty.", line_number + 1u, hap_file_name);
		}

		reader->close();
		delete reader;

		/* Calculate major allele frequencies. */
		for (unsigned int i = 0u; i < all_n_markers; ++i) {
			if (n_first_alleles[i] < n_second_alleles[i]) {
				swap_allele = all_major_alleles[i];
				all_major_alleles[i] = all_minor_alleles[i];
				all_minor_alleles[i] = swap_allele;
				all_major_allele_freqs[i] = ((double)n_second_alleles[i]) / ((double)(n_first_alleles[i] + n_second_alleles[i]));
			} else {
				all_major_allele_freqs[i] = ((double)n_first_alleles[i]) / ((double)(n_first_alleles[i] + n_second_alleles[i]));
			}
		}

		free(tokens);
		tokens = NULL;

		free(n_first_alleles);
		n_first_alleles = NULL;

		free(n_second_alleles);
		n_second_alleles = NULL;
	} catch (Exception &e) {
		e.add_message(__FILE__, __LINE__, "Error while loading '%s' file.", hap_file_name);
		throw;
	}
}

const DbView* Db::create_view(double maf_threshold, unsigned long int start_position, unsigned long int end_position) throw (Exception) {
	DbView* view = NULL;

	unsigned int start_index = 0u;
	unsigned int end_index = 0u;
	unsigned int n_markers = 0u;

	if ((start_position > 0u) || (end_position != numeric_limits<unsigned long int>::max())) {
		for (unsigned int i = 0u; i < all_n_markers; ++i) {
			if (all_positions[i] >= start_position) {
				start_index = i;
				end_index = i;
				break;
			}
		}

		for (unsigned int i = all_n_markers - 1u; i >= start_index; --i) {
			if (all_positions[i] <= end_position) {
				end_index = i;
				break;
			}
		}
	} else {
		start_index = 0u;
		end_index = all_n_markers - 1u;
	}

	if ((end_index - start_index) == 0u) {
		return NULL;
	}

	if (!isnan(maf_threshold)) {
		for (unsigned int i = start_index; i <= end_index; ++i) {
			if (auxiliary::fcmp((1.0 - all_major_allele_freqs[i]), maf_threshold, EPSILON) > 0) {
				++n_markers;
			}
		}
	} else {
		n_markers = end_index - start_index + 1u;
	}

	if (n_markers == 0u) {
		return NULL;
	}

	view = new DbView();

	views.push_back(view);

	view->n_haplotypes = n_haplotypes;

	view->n_unfiltered_markers = end_index - start_index + 1u;
	view->n_markers = n_markers;

	view->markers = (char**)malloc(view->n_markers * sizeof(char*));
	if (view->markers == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	view->positions = (unsigned long int*)malloc(view->n_markers * sizeof(unsigned long int));
	if (view->positions == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	view->major_alleles = (char*)malloc(view->n_markers * sizeof(char));
	if (view->major_alleles == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	view->minor_alleles = (char*)malloc(view->n_markers * sizeof(char));
	if (view->minor_alleles == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	view->major_allele_freqs = (double*)malloc(view->n_markers * sizeof(double));
	if (view->major_allele_freqs == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	view->haplotypes = (char**)malloc(view->n_markers * sizeof(char*));
	if (view->haplotypes == NULL) {
		throw Exception(__FILE__, __LINE__, "Error in memory allocation.");
	}

	if (!isnan(maf_threshold)) {
		for (unsigned int i = start_index, j = 0u; i <= end_index; ++i) {
			if (auxiliary::fcmp((1.0 - all_major_allele_freqs[i]), maf_threshold, EPSILON) > 0) {
				if (j >= view->n_markers) {
					throw Exception(__FILE__, __LINE__, "Array index is out of range.");
				}

				view->markers[j] = all_markers[i];
				view->positions[j] = all_positions[i];
				view->major_alleles[j] = all_major_alleles[i];
				view->minor_alleles[j] = all_minor_alleles[i];
				view->major_allele_freqs[j] = all_major_allele_freqs[i];
				view->haplotypes[j] = all_haplotypes[i];

				++j;
			}
		}
	} else {
		for (unsigned int i = start_index, j = 0u; i <= end_index; ++i) {
			if (j >= view->n_markers) {
				throw Exception(__FILE__, __LINE__, "Array index is out of range.");
			}

			view->markers[j] = all_markers[i];
			view->positions[j] = all_positions[i];
			view->major_alleles[j] = all_major_alleles[i];
			view->minor_alleles[j] = all_minor_alleles[i];
			view->major_allele_freqs[j] = all_major_allele_freqs[i];
			view->haplotypes[j] = all_haplotypes[i];

			++j;
		}
	}

	return view;
}

unsigned int Db::get_n_haplotypes() {
	return n_haplotypes;
}

unsigned int Db::get_all_n_markers() {
	return all_n_markers;
}

double Db::get_memory_usage() {
	double memory_usage = 0.0;
	vector<DbView*>::iterator views_it;

	for (views_it = views.begin(); views_it != views.end(); ++views_it) {
		memory_usage += (*views_it)->get_memory_usage();
	}

	memory_usage += (current_heap_size * sizeof(char*)) / 1048576.0;
	for (unsigned int i = 0u; i < current_heap_size; ++i) {
		if (all_markers[i] != NULL) {
			memory_usage += ((strlen(all_markers[i]) + 1u) * sizeof(char)) / 1048576.0;
		}
	}

	memory_usage += (current_heap_size * sizeof(unsigned long int)) / 1048576.0;
	memory_usage += (2u * current_heap_size * sizeof(char)) / 1048576.0;
	memory_usage += (current_heap_size * sizeof(double)) / 1048576.0;

	memory_usage += (current_heap_size * sizeof(char*)) / 1048576.0;
	for (unsigned int i = 0u; i < current_heap_size; ++i) {
		if (all_haplotypes[i] != NULL) {
			memory_usage += (n_haplotypes * sizeof(char)) / 1048576.0;
		}
	}

	return memory_usage;
}
