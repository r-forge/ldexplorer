/*
 * Copyright © 2013 Daniel Taliun and Cristian Pattaro. All rights reserved.
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

#ifndef DB_H_
#define DB_H_

#include <iostream>
#include <vector>

#include "../../auxiliary/include/auxiliary.h"
#include "../../reader/include/ReaderFactory.h"
#include "Unique.h"

using namespace std;

class Db {
private:
	static const char VCF_FIELD_SEPARATOR;
	static const char* VCF_FILE_FORMAT;
	static const char* VCF_CHROM;
	static const char* VCF_POS;
	static const char* VCF_ID;
	static const char* VCF_REF;
	static const char* VCF_ALT;
	static const char* VCF_QUAL;
	static const char* VCF_FILTER;
	static const char* VCF_INFO;
	static const char* VCF_FORMAT;
	static const char VCF_INFO_FIELD_SEPARATOR;
	static const char* VCF_VARIANT_TYPE;
	static const char* VCF_SNP_TYPE;
	static const unsigned int VCF_MANDATORY_COLUMNS_SIZE;
	static const char* vcf_mandatory_columns[];

	static const char HAPMAP2_MAP_FIELD_SEPARATOR;
	static const char HAPMAP2_HAP_FIELD_SEPARATOR;
	static const char* HAPMAP2_MAP_RS;
	static const char* HAPMAP2_MAP_POSITION;
	static const char* HAPMAP2_MAP_0;
	static const char* HAPMAP2_MAP_1;
	static const unsigned int HAPMAP2_MAP_MANDATORY_COLUMNS_SIZE;
	static const char* hapmap2_map_mandatory_columns[];

	unsigned int n_haplotypes;

	unsigned int all_n_markers;
	char** all_markers;
	unsigned long int* all_positions;
	char* all_major_alleles;
	char* all_minor_alleles;
	double* all_major_allele_freqs;
	char** all_haplotypes;

	unsigned int n_markers;
	char** markers;
	unsigned long int* positions;
	char* major_alleles;
	char* minor_alleles;
	double* major_allele_freqs;
	char** haplotypes;

	unsigned int current_heap_size;

	void free_markers(unsigned int heap_size);
	void free_positions(unsigned int heap_size);
	void free_alleles(unsigned int heap_size);
	void free_major_allele_freqs(unsigned int heap_size);
	void free_haplotypes(unsigned int heap_size);

	void reallocate() throw (Exception);

public:
	static const unsigned int HEAP_SIZE;
	static const unsigned int HEAP_INCREMENT;

	static const double EPSILON;

	static const char* VCF;
	static const char* HAPMAP2;

	Db() throw (Exception);
	virtual ~Db();

	void load_vcf(const char* file_name, unsigned long int start_position, unsigned long int end_position) throw (Exception);
	void load_vcf(const char* file_name) throw (Exception);

	void load_hapmap2(const char* map_file_name, const char* hap_file_name, unsigned long int start_position, unsigned long int end_position) throw (Exception);
	void load_hapmap2(const char* map_file_name, const char* hap_file_name) throw (Exception);

	void mask(double maf_threshold) throw (Exception);

	unsigned int get_n_haplotypes();
	unsigned int get_all_n_markers();
	unsigned int get_n_markers();
	const char* get_marker(unsigned int index);
	unsigned int get_position(unsigned int index);
	const char* get_haplotype(unsigned int index);

	double get_memory_usage();

	friend class CI;
	friend class CIWP;
	friend class CIAV;
	friend class Algorithm;
	friend class AlgorithmMIG;
	friend class AlgorithmMIGP;
	friend class AlgorithmMIGPP;
};

#endif
