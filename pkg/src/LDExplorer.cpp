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

#include <iostream>
#include <limits>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "algorithms/include/CIFactory.h"
#include "algorithms/include/AlgorithmFactory.h"
#include "algorithms/include/LD.h"
#include "db/include/Db.h"

#include <R.h>
#include <Rinternals.h>

using namespace std;

extern "C" {

	int validateBoolean(SEXP value, const char* name) {
		if (!isLogical(value)) {
			error("'%s' argument is not a boolean.", name);
		}

		if (length(value) <= 0) {
			error("'%s' argument contains no values.", name);
		}

		if (length(value) > 1) {
			error("'%s' argument contains multiple values.", name);
		}

		return LOGICAL(value)[0];
	}

	const char* validateString(SEXP value, const char* name) {
		if (!isString(value)) {
			error("'%s' argument is not a string.", name);
		}

		if (length(value) <= 0) {
			error("'%s' argument contains no values.", name);
		}

		if (length(value) > 1) {
			error("'%s' argument contains multiple values.", name);
		}

		if (STRING_ELT(value, 0) == R_NaString) {
			error("'%s' argument is NA.", name);
		}

		if (STRING_ELT(value, 0) == R_BlankString) {
			error("'%s' argument is blank.", name);
		}

		return CHAR(STRING_ELT(value, 0));
	}

	void validateStringsLengthFree(SEXP value, const char* name, vector<const char*>& c_value) {
		long int size = 0;

		if (!isString(value)) {
			error("'%s' argument is not a string.", name);
		}

		size = length(value);

		for (long int i = 0; i < size; ++i) {
			if (STRING_ELT(value, i) == R_NaString) {
				error("'%s' argument contains NA values.", name);
			}

			if (STRING_ELT(value, i) == R_BlankString) {
				error("'%s' argument contains blank values.", name);
			}

			c_value.push_back(CHAR(STRING_ELT(value, i)));
		}
	}

	double validateDouble(SEXP value, const char* name) {
		double c_value = numeric_limits<double>::quiet_NaN();

		if (!isNumeric(value)) {
			error("'%s' argument is not numeric.", name);
		}

		if (isLogical(value)) {
			error("'%s' argument is logical.", name);
		}

		if (length(value) <= 0) {
			error("'%s' argument contains no values.", name);
		}

		if (length(value) > 1) {
			error("'%s' argument contains multiple values.", name);
		}

		if (isInteger(value)) {
			c_value = (double)INTEGER(value)[0];
		} else {
			c_value = REAL(value)[0];
			if (isnan(c_value)) {
				error("'%s' argument is NA/NaN.", name);
			}
		}

		return c_value;
	}

	void validateDoubles(SEXP value, const char* name, double* c_value, unsigned int length) {
		if (!isNumeric(value)) {
			error("'%s' argument is not numeric.", name);
		}

		if (isLogical(value)) {
			error("'%s' argument is logical.", name);
		}

		if (length(value) < (long int)length) {
			error("'%s' argument contains less than %u value(s).", name, length);
		}

		if (length(value) > (long int)length) {
			error("'%s' argument contains more than %u value(s).", name, length);
		}

		if (isInteger(value)) {
			for (unsigned int i = 0; i < length; ++i) {
				c_value[i] = (double)INTEGER(value)[i];
			}
		} else {
			for (unsigned int i = 0; i < length; ++i) {
				c_value[i] = REAL(value)[i];
				if (isnan(c_value[i])) {
					error("'%s' argument contains NA/NaN value(s).", name);
				}
			}
		}
	}

	long int validateInteger(SEXP value, const char* name) {
		double c_value_double = numeric_limits<double>::quiet_NaN();
		long int c_value_int = numeric_limits<long int>::min();

		if (!isNumeric(value)) {
			error("'%s' argument is not numeric.", name);
		}

		if (isLogical(value)) {
			error("'%s' argument is logical.", name);
		}

		if (length(value) <= 0) {
			error("'%s' argument contains no values.", name);
		}

		if (length(value) > 1) {
			error("'%s' argument contains multiple values.", name);
		}

		if (isInteger(value)) {
			c_value_int = INTEGER(value)[0];
		} else {
			c_value_double = REAL(value)[0];
			if (isnan(c_value_double)) {
				error("'%s' argument is NA/NaN.", name);
			}
			c_value_int = (long int)c_value_double;
		}

		return c_value_int;
	}

	void validateIntegers(SEXP value, const char* name, long int* c_value, unsigned int length) {
		double c_value_double = numeric_limits<double>::quiet_NaN();

		if (!isNumeric(value)) {
			error("'%s' argument is not numeric.", name);
		}

		if (isLogical(value)) {
			error("'%s' argument is logical.", name);
		}

		if (length(value) < (long int)length) {
			error("'%s' argument contains less than %u value(s).", name, length);
		}

		if (length(value) > (long int)length) {
			error("'%s' argument contains more than %u value(s).", name, length);
		}

		if (isInteger(value)) {
			for (unsigned int i = 0; i < length; ++i) {
				c_value[i] = INTEGER(value)[i];
			}
		} else {
			for (unsigned int i = 0; i < length; ++i) {
				c_value_double = REAL(value)[i];
				if (isnan(c_value_double)) {
					error("'%s' argument contains NA/NaN value(s).", name);
				}
				c_value[i] = (long int)c_value_double;
			}
		}
	}

	void validateIntegersLengthFree(SEXP value, const char* name, vector<long int>& c_value) {
		long int size = 0;
		double c_value_double = numeric_limits<double>::quiet_NaN();

		if (!isNumeric(value)) {
			error("'%s' argument is not numeric.", name);
		}

		if (isLogical(value)) {
			error("'%s' argument is logical.", name);
		}

		size = length(value);

		if (isInteger(value)) {
			for (long int i = 0; i < size; ++i) {
				c_value.push_back(INTEGER(value)[i]);
			}
		} else {
			for (long int i = 0; i < size; ++i) {
				c_value_double = REAL(value)[i];
				if (isnan(c_value_double)) {
					error("'%s' argument contains NA/NaN value(s).", name);
				}
				c_value.push_back((long int)c_value_double);
			}
		}
	}

	SEXP mig(SEXP phase_file, SEXP output_file, SEXP phase_file_format, SEXP map_file,
			SEXP region, SEXP maf, SEXP ci_method, SEXP l_density, SEXP ld_ci, SEXP ehr_ci, SEXP ld_fraction,
			SEXP pruning_method, SEXP window) {

		const char* c_phase_file = NULL;
		const char* c_output_file = NULL;
		const char* c_phase_file_format = NULL;
		const char* c_map_file = NULL;
		long int c_region[2] = {numeric_limits<long int>::min(), numeric_limits<long int>::min()};
		double c_maf = numeric_limits<double>::quiet_NaN();
		const char* c_ci_method = NULL;
		long int c_l_density = numeric_limits<long int>::min();
		double c_ld_ci[2] = {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
		double c_ehr_ci = numeric_limits<double>::quiet_NaN();
		double c_ld_fraction = numeric_limits<double>::quiet_NaN();
		const char* c_pruning_method = NULL;
		long int c_window = numeric_limits<long int>::min();

//		Validate phase_file argument.
		if (!isNull(phase_file)) {
			c_phase_file = validateString(phase_file, "phase_file");
		} else {
			error("'%s' argument is NULL.", "phase_file");
		}

//		Validate output_file argument.
		if (!isNull(output_file)) {
			c_output_file = validateString(output_file, "output_file");
		} else {
			error("'%s' argument is NULL.", "output_file");
		}

//		Validate file_format argument.
		if (!isNull(phase_file_format)) {
			c_phase_file_format = validateString(phase_file_format, "file_format");
			if ((auxiliary::strcmp_ignore_case(c_phase_file_format, Db::VCF) != 0) &&
					(auxiliary::strcmp_ignore_case(c_phase_file_format, Db::HAPMAP2) != 0)) {
				error("The file format, specified in '%s' argument, must be '%s' or '%s'.", "phase_file_format", Db::VCF, Db::HAPMAP2);
			}
		} else {
			error("'%s' argument is NULL.", "phase_file_format");
		}

//		Validate legend_file argument.
		if (auxiliary::strcmp_ignore_case(c_phase_file_format, Db::HAPMAP2) == 0) {
			if (!isNull(map_file)) {
				c_map_file = validateString(map_file, "map_file");
			} else {
				error("'%s' argument is NULL.", "map_file");
			}
		}

//		Validate region argument.
		if (!isNull(region)) {
			validateIntegers(region, "region", c_region, 2u);
			if (c_region[0] < 0) {
				error("The region start position, specified in '%s' argument, must be positive.", "region");
			}
			if (c_region[1] < 0) {
				error("The region end position, specified in '%s' argument, must be positive.", "region");
			}
			if (c_region[0] >= c_region[1]) {
				error("The region end position, specified in '%s' argument, must be strictly greater than the region start position.", "region");
			}
		}

//		Validate maf argument.
		if (!isNull(maf)) {
			c_maf = validateDouble(maf, "maf");
			if ((c_maf < 0.0) || (c_maf > 0.5)) {
				error("The minor allele frequency, specified in '%s' argument, must be in [0, 0.5] interval.", "maf");
			}
		} else {
			error("'%s' argument is NULL.", "maf");
		}

//		Validate ci_method argument.
		if (!isNull(ci_method)) {
			c_ci_method = validateString(ci_method, "ci_method");
			if ((auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_WP) != 0) &&
					(auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_AV) != 0)) {
				error("The method to compute the confidence interval (CI) of D', specified in '%s' argument, must be '%s' or '%s'.", "ci_method", CI::CI_WP, CI::CI_AV);
			}
		} else {
			error("'%s' argument is NULL.", "ci_method");
		}

//		Validate likelihood density argument if WP method to compute D' CI was specified.
		if (auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_WP) == 0) {
			if (!isNull(l_density)) {
				c_l_density = validateInteger(l_density, "l_density");
				if (c_l_density <= 0) {
					error("The number of likelihood estimation points to compute confidence interval, specified in '%s' argument, must be strictly greater then 0.", "l_density");
				}
			} else {
				error("'%s' argument is NULL.", "l_density");
			}
		}

//		Validate ld_ci argument.
		if (!isNull(ld_ci)) {
			validateDoubles(ld_ci, "ld_ci", c_ld_ci, 2u);
			if ((c_ld_ci[0] < 0.0) || (c_ld_ci[0] > 1.0)) {
				error("The lower bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ld_ci");
			}
			if ((c_ld_ci[1] < 0.0) || (c_ld_ci[1] > 1.0)) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ld_ci");
			}
			if (c_ld_ci[0] >= c_ld_ci[1]) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be greater than the lower bound.", "ld_ci");
			}
		} else {
			error("'%s' argument is NULL.", "ld_ci");
		}

//		Validate ehr_ci argument.
		if (!isNull(ehr_ci)) {
			c_ehr_ci = validateDouble(ehr_ci, "ehr_ci");
			if ((c_ehr_ci < 0.0) || (c_ehr_ci > 1.0)) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ehr_ci");
			}
		} else {
			error("'%s' argument is NULL.", "ehr_ci");
		}

//		Validate ld_fraction argument.
		if (!isNull(ld_fraction)) {
			c_ld_fraction = validateDouble(ld_fraction, "ld_fraction");
			if ((c_ld_fraction <= 0.0) || (c_ld_fraction > 1.0)) {
				error("The fraction of strong LD SNP pairs within a haplotype block, specified in '%s' argument, must be in (0.0, 1.0] interval.", "ld_fraction");
			}
		} else {
			error("'%s' argument is NULL.", "ld_fraction");
		}

//		Validate pruning_method argument.
		if (!isNull(pruning_method)) {
			c_pruning_method = validateString(pruning_method, "pruning_method");
			if ((auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIG) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGP) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) != 0)) {
				error("The search space pruning method, specified in '%s' argument, must be '%s', '%s' or '%s'.",
						"file_format", Algorithm::ALGORITHM_MIG, Algorithm::ALGORITHM_MIGP, Algorithm::ALGORITHM_MIGPP);
			}
		} else {
			error("'%s' argument is NULL.", "pruning_method");
		}

//		Validate window argument if MIG++ search space pruning method was specified.
		if (auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) == 0) {
			if (!isNull(window)) {
				c_window = validateInteger(window, "window");
				if (c_window <= 0) {
					error("The window size, specified in '%s' argument, must be strictly greater than 0.", "window");
				}
			}
		}

		Algorithm* algorithm = NULL;
		Partition* partition = NULL;

		try {
			clock_t start_time = 0;
			double execution_time = 0.0;

			Db db;
			const DbView* dbview = NULL;

			Rprintf("Loading data...\n");

			start_time = clock();
			db.set_hap_file(c_phase_file);
			db.set_map_file(c_map_file);
			db.load(c_region[0] == numeric_limits<long int>::min() ? 0u : (unsigned long int)c_region[0], c_region[1] == numeric_limits<long int>::min() ? numeric_limits<unsigned long int>::max() : (unsigned long int)c_region[1], c_phase_file_format);
			dbview = db.create_view(c_maf, c_region[0] == numeric_limits<long int>::min() ? 0u : (unsigned long int)c_region[0], c_region[1] == numeric_limits<long int>::min() ? numeric_limits<unsigned long int>::max() : (unsigned long int)c_region[1]);
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			if (dbview == NULL) {
				Rprintf("\tNot enough SNPs (<= 1) in the specified region.\n");
				Rprintf("Done (%.3f sec)\n", execution_time);
				return R_NilValue;
			}

			Rprintf("\tPhase file: %s\n", c_phase_file);
			Rprintf("\tMap file: %s\n", c_map_file == NULL ? "NA" : c_map_file);
			if ((c_region[0] != numeric_limits<long int>::min()) && (c_region[1] != numeric_limits<long int>::min())) {
				Rprintf("\tRegion: [%u, %u]\n", c_region[0], c_region[1]);
			} else {
				Rprintf("\tRegion: NA\n");
			}
			Rprintf("\tMAF filter: > %g\n", dbview->maf_threshold);
			Rprintf("\tAll SNPs: %u\n", dbview->n_unfiltered_markers);
			Rprintf("\tFiltered SNPs: %u\n", dbview->n_markers);
			Rprintf("\tHaplotypes: %u\n", dbview->n_haplotypes);
			Rprintf("\tUsed memory (Mb): %.3f\n", db.get_memory_usage());
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Initializing algorithm...\n");
			start_time = clock();

			Rprintf("\tD' CI computation method: %s\n", c_ci_method);
			Rprintf("\tD' likelihood density: ");
			if (auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_WP) == 0) {
				Rprintf("%u\n", c_l_density);
			} else {
				Rprintf("NA\n");
			}
			Rprintf("\tD' CI lower bound for strong LD: >= %g\n", c_ld_ci[0]);
			Rprintf("\tD' CI upper bound for strong LD: >= %g\n", c_ld_ci[1]);
			Rprintf("\tD' CI upper bound for recombination: <= %g\n", c_ehr_ci);
			Rprintf("\tFraction of strong LD SNP pairs: >= %g\n", c_ld_fraction);
			Rprintf("\tPruning method: %s\n", c_pruning_method);
			Rprintf("\tWindow: ");
			if (auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) == 0) {
				if (c_window == numeric_limits<long int>::min()) {
					c_window = (long int)(((double)dbview->n_markers * (1.0 - c_ld_fraction)) / 2.0);
					if (c_window <= 0) {
						c_window = 1;
					}
				}
				Rprintf("%ld\n", c_window);
			} else {
				Rprintf("NA\n", c_window);
			}

			algorithm = AlgorithmFactory::create(c_pruning_method, c_window);

			algorithm->set_dbview(dbview);
			algorithm->set_ci_method(c_ci_method);
			algorithm->set_likelihood_density(c_l_density);
			algorithm->set_strong_pair_cl(c_ld_ci[0]);
			algorithm->set_strong_pair_cu(c_ld_ci[1]);
			algorithm->set_recomb_pair_cu(c_ehr_ci);
			algorithm->set_strong_pairs_fraction(c_ld_fraction);

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Processing data...\n");
			start_time = clock();

			algorithm->compute_preliminary_blocks();
			Rprintf("\tPreliminary haplotype blocks: %u\n", algorithm->get_n_preliminary_blocks());

			algorithm->sort_preliminary_blocks();
			partition = algorithm->get_block_partition();

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("\tFinal haplotype blocks: %u\n", partition->get_n_blocks());

			Rprintf("\tMemory used for preliminary haplotype blocks (Mb): %.3g\n", algorithm->get_memory_usage_preliminary_blocks());
			Rprintf("\tMemory used for final haplotype blocks (Mb): %.3g\n", partition->get_memory_usage());
			Rprintf("\tMemory used by algorithm (Mb): %.3g\n", algorithm->get_memory_usage());
			Rprintf("\tTotal used memory (Mb): %.3g\n", algorithm->get_memory_usage_preliminary_blocks() + partition->get_memory_usage() + algorithm->get_memory_usage());
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Writing results...\n");
			Rprintf("\tOutput file: %s\n", c_output_file);

			start_time = clock();
			partition->write(c_output_file);
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("Done (%.3f sec)\n", execution_time);

			delete partition;
			partition = NULL;

			delete algorithm;
			algorithm = NULL;

		} catch (Exception &e) {
			delete partition;
			partition = NULL;

			delete algorithm;
			algorithm = NULL;

			error("%s", e.what());
		}

		return R_NilValue;
	}

	SEXP mig_multi_regions(SEXP phase_file, SEXP output_files, SEXP regions_start, SEXP regions_end, SEXP processes,
			SEXP phase_file_format, SEXP map_file,
			SEXP maf, SEXP ci_method, SEXP l_density, SEXP ld_ci, SEXP ehr_ci, SEXP ld_fraction,
			SEXP pruning_method, SEXP windows) {

		const char* c_phase_file = NULL;
		vector<const char*> c_output_files;
		const char* c_phase_file_format = NULL;
		const char* c_map_file = NULL;
		vector<long int> c_regions_start;
		vector<long int> c_regions_end;
		long int c_processes = numeric_limits<long int>::min();
		double c_maf = numeric_limits<double>::quiet_NaN();
		const char* c_ci_method = NULL;
		long int c_l_density = numeric_limits<long int>::min();
		double c_ld_ci[2] = {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
		double c_ehr_ci = numeric_limits<double>::quiet_NaN();
		double c_ld_fraction = numeric_limits<double>::quiet_NaN();
		const char* c_pruning_method = NULL;
		vector<long int> c_windows;
		long int c_window = numeric_limits<long int>::min();

//		Validate phase_file argument.
		if (!isNull(phase_file)) {
			c_phase_file = validateString(phase_file, "phase_file");
		} else {
			error("'%s' argument is NULL.", "phase_file");
		}

//		Validate output_file argument.
		if (!isNull(output_files)) {
			validateStringsLengthFree(output_files, "output_files", c_output_files);
		} else {
			error("'%s' argument is NULL.", "output_files");
		}

//		Validate file_format argument.
		if (!isNull(phase_file_format)) {
			c_phase_file_format = validateString(phase_file_format, "file_format");
			if ((auxiliary::strcmp_ignore_case(c_phase_file_format, Db::VCF) != 0) &&
					(auxiliary::strcmp_ignore_case(c_phase_file_format, Db::HAPMAP2) != 0)) {
				error("The file format, specified in '%s' argument, must be '%s' or '%s'.", "phase_file_format", Db::VCF, Db::HAPMAP2);
			}
		} else {
			error("'%s' argument is NULL.", "phase_file_format");
		}

//		Validate legend_file argument.
		if (auxiliary::strcmp_ignore_case(c_phase_file_format, Db::HAPMAP2) == 0) {
			if (!isNull(map_file)) {
				c_map_file = validateString(map_file, "map_file");
			} else {
				error("'%s' argument is NULL.", "map_file");
			}
		}

//		Validate regions_start and regions_end arguments
		if (!isNull(regions_start)) {
			validateIntegersLengthFree(regions_start, "regions_start", c_regions_start);
		} else {
			error("'%s' argument is NULL.", "regions_start");
		}

		if (!isNull(regions_end)) {
			validateIntegersLengthFree(regions_end, "regions_end", c_regions_end);
		} else {
			error("'%s' argument is NULL.", "regions_end");
		}

		for (unsigned int i = 0u; i < c_regions_start.size(); ++i) {
			if (c_regions_start.at(0) < 0) {
				error("The region start positions, specified in '%s' argument, must be positive.", "regions_start");
			}
		}

		for (unsigned int i = 0u; i < c_regions_end.size(); ++i) {
			if (c_regions_end.at(0) < 0) {
				error("The region end positions, specified in '%s' argument, must be positive.", "regions_end");
			}
		}

		if (c_regions_start.size() != c_regions_end.size()) {
			error("The number of region start and end positions, specified in '%s' and '%s' arguments, must be identical.", "regions_start", "regions_end");
		}

		for (unsigned int i = 0u; i < c_regions_start.size(); ++i) {
			if (c_regions_start.at(i) >= c_regions_end.at(i)) {
				error("The region end positions, specified in '%s' argument, must be strictly greater than the region start positions, specified in '%s' argument.", "regions_end", "regions_start");
			}
		}

		if (c_output_files.size() != c_regions_start.size()) {
			error("The number of the specified output files must correspond to the number of the specified regions.");
		}

//		Validate processes argument.
		if (!isNull(processes)) {
			c_processes = validateInteger(processes, "processes");
			if (c_processes < 1) {
				error("The number of processes, specified in '%s' argument, must be greater than 0.", "processes");
			}
		} else {
			error("'%s' argument is NULL.", "processes");
		}

//		Validate maf argument.
		if (!isNull(maf)) {
			c_maf = validateDouble(maf, "maf");
			if ((c_maf < 0.0) || (c_maf > 0.5)) {
				error("The minor allele frequency, specified in '%s' argument, must be in [0, 0.5] interval.", "maf");
			}
		} else {
			error("'%s' argument is NULL.", "maf");
		}

//		Validate ci_method argument.
		if (!isNull(ci_method)) {
			c_ci_method = validateString(ci_method, "ci_method");
			if ((auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_WP) != 0) &&
					(auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_AV) != 0)) {
				error("The method to compute the confidence interval (CI) of D', specified in '%s' argument, must be '%s' or '%s'.", "ci_method", CI::CI_WP, CI::CI_AV);
			}
		} else {
			error("'%s' argument is NULL.", "ci_method");
		}

//		Validate likelihood density argument if WP method to compute D' CI was specified.
		if (auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_WP) == 0) {
			if (!isNull(l_density)) {
				c_l_density = validateInteger(l_density, "l_density");
				if (c_l_density <= 0) {
					error("The number of likelihood estimation points to compute confidence interval, specified in '%s' argument, must be strictly greater then 0.", "l_density");
				}
			} else {
				error("'%s' argument is NULL.", "l_density");
			}
		}

//		Validate ld_ci argument.
		if (!isNull(ld_ci)) {
			validateDoubles(ld_ci, "ld_ci", c_ld_ci, 2u);
			if ((c_ld_ci[0] < 0.0) || (c_ld_ci[0] > 1.0)) {
				error("The lower bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ld_ci");
			}
			if ((c_ld_ci[1] < 0.0) || (c_ld_ci[1] > 1.0)) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ld_ci");
			}
			if (c_ld_ci[0] >= c_ld_ci[1]) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be greater than the lower bound.", "ld_ci");
			}
		} else {
			error("'%s' argument is NULL.", "ld_ci");
		}

//		Validate ehr_ci argument.
		if (!isNull(ehr_ci)) {
			c_ehr_ci = validateDouble(ehr_ci, "ehr_ci");
			if ((c_ehr_ci < 0.0) || (c_ehr_ci > 1.0)) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ehr_ci");
			}
		} else {
			error("'%s' argument is NULL.", "ehr_ci");
		}

//		Validate ld_fraction argument.
		if (!isNull(ld_fraction)) {
			c_ld_fraction = validateDouble(ld_fraction, "ld_fraction");
			if ((c_ld_fraction <= 0.0) || (c_ld_fraction > 1.0)) {
				error("The fraction of strong LD SNP pairs within a haplotype block, specified in '%s' argument, must be in (0.0, 1.0] interval.", "ld_fraction");
			}
		} else {
			error("'%s' argument is NULL.", "ld_fraction");
		}

//		Validate pruning_method argument.
		if (!isNull(pruning_method)) {
			c_pruning_method = validateString(pruning_method, "pruning_method");
			if ((auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIG) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGP) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) != 0)) {
				error("The search space pruning method, specified in '%s' argument, must be '%s', '%s' or '%s'.",
						"file_format", Algorithm::ALGORITHM_MIG, Algorithm::ALGORITHM_MIGP, Algorithm::ALGORITHM_MIGPP);
			}
		} else {
			error("'%s' argument is NULL.", "pruning_method");
		}

//		Validate window argument if MIG++ search space pruning method was specified.
		if (auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) == 0) {
			if (!isNull(windows)) {
				validateIntegersLengthFree(windows, "windows", c_windows);

				for (unsigned int i = 0u; i < c_windows.size(); ++i) {
					if (c_windows.at(i) <= 0) {
						error("The window sizes, specified in '%s' argument, must be strictly greater than 0.", "windows");
					}
				}

				if (c_output_files.size() != c_windows.size()) {
					error("The number of the specified output files must correspond to the number of the specified windows.");
				}
			}
		}

		Algorithm* algorithm = NULL;
		vector<Algorithm*> algorithms;

		Partition* partition = NULL;
		vector<Partition*> partitions;

		try {
			clock_t start_time = 0;
			double start_time_omp = 0.0;
			double execution_time = 0.0;

			Db db;
			const DbView* dbview = NULL;
			vector<const DbView*> dbviews;
			bool all_empty = true;
			int omp_i = 0;

			Rprintf("Loading data...\n");

			start_time = clock();
			db.set_hap_file(c_phase_file);
			db.set_map_file(c_map_file);
			db.load(0u, numeric_limits<unsigned long int>::max(), c_phase_file_format);
			for (unsigned int i = 0u; i < c_output_files.size(); ++i) {
				dbview = db.create_view(c_maf, c_regions_start.at(i), c_regions_end.at(i));
				dbviews.push_back(dbview);
				if ((all_empty == true) && (dbview != NULL)) {
					all_empty = false;
				}
			}
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			if (all_empty == true) {
				Rprintf("\tNot enough SNPs (<= 1) in any of the specified regions.\n");
				Rprintf("Done (%.3f sec)\n", execution_time);
				return R_NilValue;
			}

			Rprintf("\tPhase file: %s\n", c_phase_file);
			Rprintf("\tMap file: %s\n", c_map_file == NULL ? "NA" : c_map_file);
			Rprintf("\tRegions:\n");
			for (unsigned int i = 0u; i < dbviews.size(); ++i) {
				dbview = dbviews.at(i);
				if (dbview != NULL) {
					Rprintf("\t- Region: [%u, %u]\n", c_regions_start.at(i), c_regions_end.at(i));
					Rprintf("\t--  MAF filter: > %g\n", dbview->maf_threshold);
					Rprintf("\t--  All SNPs: %u\n", dbview->n_unfiltered_markers);
					Rprintf("\t--  Filtered SNPs: %u\n", dbview->n_markers);
					Rprintf("\t--  Haplotypes: %u\n", dbview->n_haplotypes);
				}
			}
			Rprintf("\tUsed memory (Mb): %.3f\n", db.get_memory_usage());
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Initializing algorithm...\n");

			start_time = clock();
			if (c_windows.size() == 0) {
				for (unsigned int i = 0u; i < dbviews.size(); ++i) {
					c_window = numeric_limits<long int>::min();
					dbview = dbviews.at(i);
					if (dbview != NULL) {
						c_window = (long int)(((double)dbview->n_markers * (1.0 - c_ld_fraction)) / 2.0);
						if (c_window <= 0) {
							c_window = 1;
						}
					}
					c_windows.push_back(c_window);
				}
			}

			for (unsigned int i = 0u; i < dbviews.size(); ++i) {
				algorithm = NULL;
				dbview = dbviews.at(i);
				if (dbview != NULL) {
					algorithm = AlgorithmFactory::create(c_pruning_method, c_windows.at(i));
					algorithm->set_dbview(dbview);
					algorithm->set_ci_method(c_ci_method);
					algorithm->set_likelihood_density(c_l_density);
					algorithm->set_strong_pair_cl(c_ld_ci[0]);
					algorithm->set_strong_pair_cu(c_ld_ci[1]);
					algorithm->set_recomb_pair_cu(c_ehr_ci);
					algorithm->set_strong_pairs_fraction(c_ld_fraction);
				}
				algorithms.push_back(algorithm);
			}
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("\tD' CI computation method: %s\n", c_ci_method);
			Rprintf("\tD' likelihood density: ");
			if (auxiliary::strcmp_ignore_case(c_ci_method, CI::CI_WP) == 0) {
				Rprintf("%u\n", c_l_density);
			} else {
				Rprintf("NA\n");
			}
			Rprintf("\tD' CI lower bound for strong LD: >= %g\n", c_ld_ci[0]);
			Rprintf("\tD' CI upper bound for strong LD: >= %g\n", c_ld_ci[1]);
			Rprintf("\tD' CI upper bound for recombination: <= %g\n", c_ehr_ci);
			Rprintf("\tFraction of strong LD SNP pairs: >= %g\n", c_ld_fraction);
			Rprintf("\tPruning method: %s\n", c_pruning_method);


			Rprintf("\tWindows: ");
			if (auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) == 0) {
				Rprintf("\n");
				for (unsigned int i = 0u; i < dbviews.size(); ++i) {
					dbview = dbviews.at(i);
					if (dbview != NULL) {
						Rprintf("\t- Region: [%u, %u]\n", c_regions_start.at(i), c_regions_end.at(i));
						Rprintf("\t-- Window: %ld\n", c_windows.at(i));
					}
				}
			} else {
				Rprintf("NA\n", c_window);
			}


			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Processing data (%d processes)...\n", c_processes);

#ifdef	_OPENMP
			start_time_omp = omp_get_wtime();
#else
			start_time = clock();
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(c_processes) private(omp_i, algorithm) schedule(dynamic, 1)
#endif
			for (omp_i = 0; omp_i < (int)algorithms.size(); ++omp_i) {
				algorithm = algorithms.at(omp_i);
				if (algorithm != NULL) {
					algorithm->compute_preliminary_blocks();
					algorithm->sort_preliminary_blocks();
				}
			}

			for (unsigned int i = 0; i < algorithms.size(); ++i) {
				partition = NULL;
				algorithm = algorithms.at(i);
				if (algorithm != NULL) {
					partition = algorithm->get_block_partition();
				}
				partitions.push_back(partition);
			}

#ifdef	_OPENMP
			execution_time = omp_get_wtime() - start_time_omp;
#else
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
#endif

			for (unsigned int i = 0u; i < dbviews.size(); ++i) {
				dbview = dbviews.at(i);
				if (dbview != NULL) {
					algorithm = algorithms.at(i);
					partition = partitions.at(i);
					Rprintf("\t- Region: [%u, %u]\n", c_regions_start.at(i), c_regions_end.at(i));
					Rprintf("\t-- Preliminary haplotype blocks: %u\n", algorithm->get_n_preliminary_blocks());
					Rprintf("\t-- Final haplotype blocks: %u\n", partition->get_n_blocks());
					Rprintf("\t-- Memory used for preliminary haplotype blocks (Mb): %.3g\n", algorithm->get_memory_usage_preliminary_blocks());
					Rprintf("\t-- Memory used for final haplotype blocks (Mb): %.3g\n", partition->get_memory_usage());
					Rprintf("\t-- Memory used by algorithm (Mb): %.3g\n", algorithm->get_memory_usage());
					Rprintf("\t-- Total used memory (Mb): %.3g\n", algorithm->get_memory_usage_preliminary_blocks() + partition->get_memory_usage() + algorithm->get_memory_usage());
				}
			}
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Writing results...\n");

			start_time = clock();
			for (unsigned int i = 0; i < partitions.size(); ++i) {
				Rprintf("\tOutput file: %s\n", c_output_files.at(i));

				partition = partitions.at(i);
				if (partition != NULL) {
					partition->write(c_output_files.at(i));
				}
			}
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("Done (%.3f sec)\n", execution_time);

			for (unsigned int i = 0u; i < partitions.size(); ++i) {
				partition = partitions.at(i);
				if (partition != NULL) {
					delete partition;
				}
			}
			partitions.clear();

			for (unsigned int i = 0u; i < algorithms.size(); ++i) {
				algorithm = algorithms.at(i);
				if (algorithm != NULL) {
					delete algorithm;
				}
			}
			algorithms.clear();

		} catch (Exception &e) {
			for (unsigned int i = 0u; i < partitions.size(); ++i) {
				partition = partitions.at(i);
				if (partition != NULL) {
					delete partition;
				}
			}
			partitions.clear();

			for (unsigned int i = 0u; i < algorithms.size(); ++i) {
				algorithm = algorithms.at(i);
				if (algorithm != NULL) {
					delete algorithm;
				}
			}
			algorithms.clear();

			error("%s", e.what());
		}

		return R_NilValue;
	}

	SEXP mig_rsq(SEXP phase_file, SEXP output_file, SEXP phase_file_format, SEXP map_file,
			SEXP region, SEXP maf, SEXP weak_rsq, SEXP strong_rsq, SEXP fraction,
			SEXP pruning_method, SEXP window) {

		const char* c_phase_file = NULL;
		const char* c_output_file = NULL;
		const char* c_phase_file_format = NULL;
		const char* c_map_file = NULL;
		long int c_region[2] = {numeric_limits<long int>::min(), numeric_limits<long int>::min()};
		double c_maf = numeric_limits<double>::quiet_NaN();
		double c_weak_rsq = numeric_limits<double>::quiet_NaN();
		double c_strong_rsq = numeric_limits<double>::quiet_NaN();
		double c_fraction = numeric_limits<double>::quiet_NaN();
		const char* c_pruning_method = NULL;
		long int c_window = numeric_limits<long int>::min();

//		Validate phase_file argument.
		if (!isNull(phase_file)) {
			c_phase_file = validateString(phase_file, "phase_file");
		} else {
			error("'%s' argument is NULL.", "phase_file");
		}

//		Validate output_file argument.
		if (!isNull(output_file)) {
			c_output_file = validateString(output_file, "output_file");
		} else {
			error("'%s' argument is NULL.", "output_file");
		}

//		Validate file_format argument.
		if (!isNull(phase_file_format)) {
			c_phase_file_format = validateString(phase_file_format, "file_format");
			if ((auxiliary::strcmp_ignore_case(c_phase_file_format, Db::VCF) != 0) &&
					(auxiliary::strcmp_ignore_case(c_phase_file_format, Db::HAPMAP2) != 0)) {
				error("The file format, specified in '%s' argument, must be '%s' or '%s'.", "phase_file_format", Db::VCF, Db::HAPMAP2);
			}
		} else {
			error("'%s' argument is NULL.", "phase_file_format");
		}

//		Validate legend_file argument.
		if (auxiliary::strcmp_ignore_case(c_phase_file_format, Db::HAPMAP2) == 0) {
			if (!isNull(map_file)) {
				c_map_file = validateString(map_file, "map_file");
			} else {
				error("'%s' argument is NULL.", "map_file");
			}
		}

//		Validate region argument.
		if (!isNull(region)) {
			validateIntegers(region, "region", c_region, 2u);
			if (c_region[0] < 0) {
				error("The region start position, specified in '%s' argument, must be positive.", "region");
			}
			if (c_region[1] < 0) {
				error("The region end position, specified in '%s' argument, must be positive.", "region");
			}
			if (c_region[0] >= c_region[1]) {
				error("The region end position, specified in '%s' argument, must be strictly greater than the region start position.", "region");
			}
		}

//		Validate maf argument.
		if (!isNull(maf)) {
			c_maf = validateDouble(maf, "maf");
			if ((c_maf < 0.0) || (c_maf > 0.5)) {
				error("The minor allele frequency, specified in '%s' argument, must be in [0, 0.5] interval.", "maf");
			}
		} else {
			error("'%s' argument is NULL.", "maf");
		}

//		Validate weak_rsq argument.
		if (!isNull(weak_rsq)) {
			c_weak_rsq = validateDouble(weak_rsq, "weak_rsq");
			if ((c_weak_rsq <= 0.0) || (c_weak_rsq > 1.0)) {
				error("The upper bound of the r^2 for the weak LD SNP pairs, specified in '%s' argument, must be in (0, 1] interval.", "weak_rsq");
			}
		} else {
			error("'%s' argument is NULL.", "weak_rsq");
		}

//		Validate strong_rsq argument.
		if (!isNull(strong_rsq)) {
			c_strong_rsq = validateDouble(strong_rsq, "strong_rsq");
			if ((c_strong_rsq <= 0.0) || (c_strong_rsq > 1.0)) {
				error("The lower bound of the r^2 for the strong LD SNP pairs, specified in '%s' argument, must be in (0, 1] interval.", "strong_rsq");
			}
		} else {
			error("'%s' argument is NULL.", "strong_rsq");
		}

//		Validate ld_fraction argument.
		if (!isNull(fraction)) {
			c_fraction = validateDouble(fraction, "fraction");
			if ((c_fraction <= 0.0) || (c_fraction > 1.0)) {
				error("The fraction of the strong LD SNP pairs within a haplotype block, specified in '%s' argument, must be in (0.0, 1.0] interval.", "fraction");
			}
		} else {
			error("'%s' argument is NULL.", "fraction");
		}

//		Validate pruning_method argument.
		if (!isNull(pruning_method)) {
			c_pruning_method = validateString(pruning_method, "pruning_method");
			if ((auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIG) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGP) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) != 0)) {
				error("The search space pruning method, specified in '%s' argument, must be '%s', '%s' or '%s'.",
						"file_format", Algorithm::ALGORITHM_MIG, Algorithm::ALGORITHM_MIGP, Algorithm::ALGORITHM_MIGPP);
			}
		} else {
			error("'%s' argument is NULL.", "pruning_method");
		}

//		Validate window argument if MIG++ search space pruning method was specified.
		if (auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) == 0) {
			if (!isNull(window)) {
				c_window = validateInteger(window, "window");
				if (c_window <= 0) {
					error("The window size, specified in '%s' argument, must be strictly greater than 0.", "window");
				}
			}
		}

		Algorithm* algorithm = NULL;
		Partition* partition = NULL;

		try {
			clock_t start_time = 0;
			double execution_time = 0.0;

			Db db;
			const DbView* dbview = NULL;

			Rprintf("Loading data...\n");

			start_time = clock();
			db.set_hap_file(c_phase_file);
			db.set_map_file(c_map_file);
			db.load(c_region[0] == numeric_limits<long int>::min() ? 0u : (unsigned long int)c_region[0], c_region[1] == numeric_limits<long int>::min() ? numeric_limits<unsigned long int>::max() : (unsigned long int)c_region[1], c_phase_file_format);
			dbview = db.create_view(c_maf, c_region[0] == numeric_limits<long int>::min() ? 0u : (unsigned long int)c_region[0], c_region[1] == numeric_limits<long int>::min() ? numeric_limits<unsigned long int>::max() : (unsigned long int)c_region[1]);
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			if (dbview == NULL) {
				Rprintf("\tNot enough SNPs (<= 1) in the specified region.\n");
				Rprintf("Done (%.3f sec)\n", execution_time);
				return R_NilValue;
			}

			Rprintf("\tPhase file: %s\n", c_phase_file);
			Rprintf("\tMap file: %s\n", c_map_file == NULL ? "NA" : c_map_file);
			if ((c_region[0] != numeric_limits<long int>::min()) && (c_region[1] != numeric_limits<long int>::min())) {
				Rprintf("\tRegion: [%u, %u]\n", c_region[0], c_region[1]);
			} else {
				Rprintf("\tRegion: NA\n");
			}
			Rprintf("\tMAF filter: > %g\n", dbview->maf_threshold);
			Rprintf("\tAll SNPs: %u\n", dbview->n_unfiltered_markers);
			Rprintf("\tFiltered SNPs: %u\n", dbview->n_markers);
			Rprintf("\tHaplotypes: %u\n", dbview->n_haplotypes);
			Rprintf("\tUsed memory (Mb): %.3f\n", db.get_memory_usage());
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Initializing algorithm...\n");
			start_time = clock();

			Rprintf("\tWeak LD SNP pairs r^2: < %g\n", c_weak_rsq);
			Rprintf("\tStrong LD SNP pairs r^2: >= %g\n", c_strong_rsq);
			Rprintf("\tFraction of strong LD SNP pairs: >= %g\n", c_fraction);
			Rprintf("\tPruning method: %s\n", c_pruning_method);
			Rprintf("\tWindow: ");
			if (auxiliary::strcmp_ignore_case(c_pruning_method, Algorithm::ALGORITHM_MIGPP) == 0) {
				if (c_window == numeric_limits<long int>::min()) {
					c_window = (long int)(((double)dbview->n_markers * (1.0 - c_fraction)) / 2.0);
					if (c_window <= 0) {
						c_window = 1;
					}
				}
				Rprintf("%ld\n", c_window);
			} else {
				Rprintf("NA\n", c_window);
			}

			algorithm = AlgorithmFactory::create(c_pruning_method, c_window);

			algorithm->set_dbview(dbview);
			algorithm->set_weak_pair_rsq(c_weak_rsq);
			algorithm->set_strong_pair_rsq(c_strong_rsq);
			algorithm->set_strong_pairs_fraction(c_fraction);

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Processing data...\n");
			start_time = clock();

			algorithm->compute_preliminary_blocks_rsq();
			Rprintf("\tPreliminary haplotype blocks: %u\n", algorithm->get_n_preliminary_blocks());

			algorithm->sort_preliminary_blocks();
			partition = algorithm->get_block_partition();

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("\tFinal haplotype blocks: %u\n", partition->get_n_blocks());

			Rprintf("\tMemory used for preliminary haplotype blocks (Mb): %.3g\n", algorithm->get_memory_usage_preliminary_blocks());
			Rprintf("\tMemory used for final haplotype blocks (Mb): %.3g\n", partition->get_memory_usage());
			Rprintf("\tMemory used by algorithm (Mb): %.3g\n", algorithm->get_memory_usage());
			Rprintf("\tTotal used memory (Mb): %.3g\n", algorithm->get_memory_usage_preliminary_blocks() + partition->get_memory_usage() + algorithm->get_memory_usage());
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Writing results...\n");
			Rprintf("\tOutput file: %s\n", c_output_file);

			start_time = clock();
			partition->write(c_output_file);
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("Done (%.3f sec)\n", execution_time);

			delete partition;
			partition = NULL;

			delete algorithm;
			algorithm = NULL;

		} catch (Exception &e) {
			delete partition;
			partition = NULL;

			delete algorithm;
			algorithm = NULL;

			error("%s", e.what());
		}

		return R_NilValue;
	}

	SEXP ld(SEXP phase_file, SEXP snps_file, SEXP output_file, SEXP window, SEXP coefficient, SEXP maf, SEXP gzip) {
		const char* c_phase_file = NULL;
		const char* c_snps_file = NULL;
		const char* c_output_file = NULL;
		long int c_window = numeric_limits<long int>::min();
		const char* c_coefficient = NULL;
		double c_maf = numeric_limits<double>::quiet_NaN();
		int c_gzip = 0;

//		Validate phase_file argument.
		if (!isNull(phase_file)) {
			c_phase_file = validateString(phase_file, "phase_file");
		} else {
			error("'%s' argument is NULL.", "phase_file");
		}

//		Validate phase_file argument.
		if (!isNull(snps_file)) {
			c_snps_file = validateString(snps_file, "snps_file");
		} else {
			error("'%s' argument is NULL.", "snps_file");
		}

//		Validate output_file argument.
		if (!isNull(output_file)) {
			c_output_file = validateString(output_file, "output_file");
		} else {
			error("'%s' argument is NULL.", "output_file");
		}

//		Validate window argument.
		if (!isNull(window)) {
			c_window = validateInteger(window, "window");
			if (c_window <= 0) {
				error("The window size, specified in '%s' argument, must be strictly greater than 0 base pairs.", "window");
			}
		} else {
			error("'%s' argument is NULL.", "window");
		}

//		Validate coefficient argument.
		if (!isNull(coefficient)) {
			c_coefficient = validateString(coefficient, "coefficient");
			if ((auxiliary::strcmp_ignore_case(c_coefficient, LD::D) != 0) &&	(auxiliary::strcmp_ignore_case(c_coefficient, LD::DPRIME) != 0) && (auxiliary::strcmp_ignore_case(c_coefficient, LD::R2) != 0)) {
				error("The LD coefficient, specified in '%s' argument, must be '%s', '%s' or '%s'.", "coefficient", LD::D, LD::DPRIME, LD::R2);
			}
		} else {
			error("'%s' argument is NULL.", "coefficient");
		}

//		Validate maf argument.
		if (!isNull(maf)) {
			c_maf = validateDouble(maf, "maf");
			if ((c_maf < 0.0) || (c_maf > 0.5)) {
				error("The minor allele frequency, specified in '%s' argument, must be in [0, 0.5] interval.", "maf");
			}
		} else {
			error("'%s' argument is NULL.", "maf");
		}

//		Validate gzip argument.
		if (!isNull(gzip)) {
			c_gzip = validateBoolean(gzip, "gzip");
		} else {
			error("'%s' argument is NULL.", "gzip");
		}

		try {
			clock_t start_time = 0;
			double execution_time = 0.0;

			Db db;
			const DbView* dbview = NULL;
			LD ld;

			Rprintf("Loading data...\n");

			start_time = clock();
			db.set_hap_file(c_phase_file);
			db.load(0u, numeric_limits<unsigned long int>::max(), Db::VCF);
			dbview = db.create_view(c_maf, 0u, numeric_limits<unsigned long int>::max());
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			if (dbview == NULL) {
				Rprintf("\tNot enough SNPs (<= 1) in the specified region.\n");
				Rprintf("Done (%.3f sec)\n", execution_time);
				return R_NilValue;
			}

			Rprintf("\tPhase file: %s\n", c_phase_file);
			Rprintf("\tMAF filter: > %g\n", dbview->maf_threshold);
			Rprintf("\tAll SNPs: %u\n", dbview->n_unfiltered_markers);
			Rprintf("\tFiltered SNPs: %u\n", dbview->n_markers);
			Rprintf("\tHaplotypes: %u\n", dbview->n_haplotypes);
			Rprintf("\tUsed memory (Mb): %.3f\n", db.get_memory_usage());
			Rprintf("Done (%.3f sec)\n", execution_time);

			ld.set_dbview(dbview);

			Rprintf("Calculating LD coefficients...\n");

			Rprintf("\tLD coefficient: %s\n", c_coefficient);
			Rprintf("\tWindow: +/-%d bp\n", c_window);

			start_time = clock();
			ld.load_markers(c_snps_file);
			ld.index_db_markers();
			ld.compute_ld(c_output_file, c_coefficient, c_window, c_gzip);
			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("\tInput SNPs: %u\n", ld.get_n_snps());
			Rprintf("\tOuptu file: %s\n", c_output_file);
			Rprintf("\tUsed memory (Mb): %.3f\n", ld.get_used_memory());
			Rprintf("Done (%.3f sec)\n", execution_time);

		} catch (Exception &e) {
			error("%s", e.what());
		}

		return R_NilValue;
	}
}

int main(int args, char** argv) {
	return 0;
}
