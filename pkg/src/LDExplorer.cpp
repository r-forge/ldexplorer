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

#include "algorithms/include/CIFactory.h"
#include "algorithms/include/AlgorithmFactory.h"
#include "db/include/Db.h"

#include <R.h>
#include <Rinternals.h>

using namespace std;

#define VERSION "LDExplorer 1.0.0"

extern "C" {

	struct config {
		const char* phase_file;
		const char* output_file;
		const char* phase_file_format;
		const char* map_file;
		bool is_region;
		long int region[2];
		double maf;
		const char* ci_method;
		long int l_density;
		double ld_ci[2];
		double ehr_ci;
		double ld_fraction;
		const char* pruning_method;
		long int window;

		unsigned int haplotypes;
		unsigned int all_snps;
		unsigned int filtered_snps;

		config() : phase_file(NULL), output_file(NULL), phase_file_format(NULL), map_file(NULL),
				is_region(false), maf(0.0),	ci_method(NULL), l_density(0), ehr_ci(0.0), ld_fraction(0.0),
				pruning_method(NULL), window(numeric_limits<long int>::max()),
				all_snps(0u), filtered_snps(0u) {

			region[0] = 0;
			region[1] = numeric_limits<long int>::max();

			ld_ci[0] = ld_ci[1] = 0.0;
		};
	};

	void write_info(const char* output_file, config& c) throw (Exception) {
		Writer* writer = NULL;

		try{
			writer = WriterFactory::create(WriterFactory::TEXT);
			writer->set_file_name(output_file);
			writer->open(true);

			writer->write("# VERSION: %s\n", VERSION);
			writer->write("# PHASE FILE: %s\n", c.phase_file);
			writer->write("# MAP FILE: %s\n", c.map_file == NULL ? "NA" : c.map_file);
			if ((c.region[0] > 0) || (c.region[1] != numeric_limits<long int>::max())) {
				writer->write("# REGION: [%u, %u]\n", c.region[0], c.region[1]);
			} else {
				writer->write("# REGION: NA\n");
			}
			writer->write("# MAF FILTER: > %g\n", c.maf);
			writer->write("# ALL SNPs: %u\n", c.all_snps);
			writer->write("# FILTERED SNPs: %u\n", c.filtered_snps);
			writer->write("# HAPLOTYPES: %u\n", c.haplotypes);
			writer->write("# D' CI COMPUTATION METHOD: %s\n", c.ci_method);
			if (auxiliary::strcmp_ignore_case(c.ci_method, CIFactory::CI_WP) == 0) {
				writer->write("# D' LIKELIHOOD DENSITY: %u\n", c.l_density);
			} else {
				writer->write("# D' LIKELIHOOD DENSITY: NA");
			}
			writer->write("# D' CI LOWER BOUND FOR STRONG LD: >= %g\n", c.ld_ci[0]);
			writer->write("# D' CI UPPER BOUND FOR STRONG LD: >= %g\n", c.ld_ci[1]);
			writer->write("# D' CI UPPER BOUND FOR RECOMBINATION: <= %g\n", c.ehr_ci);
			writer->write("# FRACTION OF STRONG LD SNP PAIRS: >= %g\n", c.ld_fraction);
			writer->write("# PRUNING METHOD: %s\n", c.pruning_method);
			if (auxiliary::strcmp_ignore_case(c.pruning_method, AlgorithmFactory::ALGORITHM_MIGPP) == 0) {
				writer->write("# WINDOW: %ld\n", c.window);
			} else {
				writer->write("# WINDOW: NA\n", c.window);
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
				c_value[i] = (double)INTEGER(value)[i];
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

	SEXP mig(SEXP phase_file, SEXP output_file, SEXP phase_file_format, SEXP map_file,
			SEXP region, SEXP maf, SEXP ci_method, SEXP l_density, SEXP ld_ci, SEXP ehr_ci, SEXP ld_fraction,
			SEXP pruning_method, SEXP window) {

		config c;

		bool is_default_window = false;

//		Validate phase_file argument.
		if (!isNull(phase_file)) {
			c.phase_file = validateString(phase_file, "phase_file");
		} else {
			error("'%s' argument is NULL.", "phase_file");
		}

//		Validate output_file argument.
		if (!isNull(output_file)) {
			c.output_file = validateString(output_file, "output_file");
		} else {
			error("'%s' argument is NULL.", "output_file");
		}

//		Validate file_format argument.
		if (!isNull(phase_file_format)) {
			c.phase_file_format = validateString(phase_file_format, "file_format");
			if ((auxiliary::strcmp_ignore_case(c.phase_file_format, Db::VCF) != 0) &&
					(auxiliary::strcmp_ignore_case(c.phase_file_format, Db::HAPMAP2) != 0)) {
				error("The file format, specified in '%s' argument, must be '%s' or '%s'.", "phase_file_format", Db::VCF, Db::HAPMAP2);
			}
		} else {
			error("'%s' argument is NULL.", "phase_file_format");
		}

//		Validate legend_file argument.
		if (auxiliary::strcmp_ignore_case(c.phase_file_format, Db::HAPMAP2) == 0) {
			if (!isNull(map_file)) {
				c.map_file = validateString(map_file, "map_file");
			} else {
				error("'%s' argument is NULL.", "map_file");
			}
		}

//		Validate region argument.
		if (!isNull(region)) {
			validateIntegers(region, "region", c.region, 2u);
			if (c.region[0] < 0) {
				error("The region start position, specified in '%s' argument, must be positive.", "region");
			}
			if (c.region[1] < 0) {
				error("The region end position, specified in '%s' argument, must be positive.", "region");
			}
			if (c.region[0] >= c.region[1]) {
				error("The region end position, specified in '%s' argument, must be strictly greater than the region start position.", "region");
			}
			c.is_region = true;
		}

//		Validate maf argument.
		if (!isNull(maf)) {
			c.maf = validateDouble(maf, "maf");
			if ((c.maf < 0.0) || (c.maf > 0.5)) {
				error("The minor allele frequency, specified in '%s' argument, must be in [0, 0.5] interval.", "maf");
			}
		} else {
			error("'%s' argument is NULL.", "maf");
		}

//		Validate ci_method argument.
		if (!isNull(ci_method)) {
			c.ci_method = validateString(ci_method, "ci_method");
			if ((auxiliary::strcmp_ignore_case(c.ci_method, CIFactory::CI_WP) != 0) &&
					(auxiliary::strcmp_ignore_case(c.ci_method, CIFactory::CI_AV) != 0)) {
				error("The method to compute the confidence interval (CI) of D', specified in '%s' argument, must be '%s' or '%s'.", "file_format", CIFactory::CI_WP, CIFactory::CI_AV);
			}
		} else {
			error("'%s' argument is NULL.", "ci_method");
		}

//		Validate ci_precision argument if WP method to compute D' CI was specified.
		if (auxiliary::strcmp_ignore_case(c.ci_method, CIFactory::CI_WP) == 0) {
			if (!isNull(l_density)) {
				c.l_density = validateInteger(l_density, "ci_precision");
				if (c.l_density <= 0) {
					error("The number of likelihood estimation points to compute confidence interval, specified in '%s' argument, must be strictly greater then 0.", "l_density");
				}
			} else {
				error("'%s' argument is NULL.", "l_density");
			}
		}

//		Validate ld_ci argument.
		if (!isNull(ld_ci)) {
			validateDoubles(ld_ci, "ld_ci", c.ld_ci, 2u);
			if ((c.ld_ci[0] < 0.0) || (c.ld_ci[0] > 1.0)) {
				error("The lower bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ld_ci");
			}
			if ((c.ld_ci[1] < 0.0) || (c.ld_ci[1] > 1.0)) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ld_ci");
			}
			if (c.ld_ci[0] >= c.ld_ci[1]) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be greater than the lower bound.", "ld_ci");
			}
		} else {
			error("'%s' argument is NULL.", "ld_ci");
		}

//		Validate ehr_ci argument.
		if (!isNull(ehr_ci)) {
			c.ehr_ci = validateDouble(ehr_ci, "ehr_ci");
			if ((c.ehr_ci < 0.0) || (c.ehr_ci > 1.0)) {
				error("The upper bound of confidence interval, specified in '%s' argument, must be in [0, 1] interval.", "ehr_ci");
			}
		} else {
			error("'%s' argument is NULL.", "ehr_ci");
		}

//		Validate ld_fraction argument.
		if (!isNull(ld_fraction)) {
			c.ld_fraction = validateDouble(ld_fraction, "ld_fraction");
			if ((c.ld_fraction <= 0.0) || (c.ld_fraction > 1.0)) {
				error("The fraction of strong LD SNP pairs within a haplotype block, specified in '%s' argument, must be in (0.0, 1.0] interval.", "ld_fraction");
			}
		} else {
			error("'%s' argument is NULL.", "ld_fraction");
		}

//		Validate pruning_method argument.
		if (!isNull(pruning_method)) {
			c.pruning_method = validateString(pruning_method, "pruning_method");
			if ((auxiliary::strcmp_ignore_case(c.pruning_method, AlgorithmFactory::ALGORITHM_MIG) != 0) &&
					(auxiliary::strcmp_ignore_case(c.pruning_method, AlgorithmFactory::ALGORITHM_MIGP) != 0) &&
					(auxiliary::strcmp_ignore_case(c.pruning_method, AlgorithmFactory::ALGORITHM_MIGPP) != 0)) {
				error("The search space pruning method, specified in '%s' argument, must be '%s', '%s' or '%s'.",
						"file_format", AlgorithmFactory::ALGORITHM_MIG, AlgorithmFactory::ALGORITHM_MIGP, AlgorithmFactory::ALGORITHM_MIGPP);
			}
		} else {
			error("'%s' argument is NULL.", "pruning_method");
		}

//		Validate window argument if MIG++ search space pruning method was specified.
		if (auxiliary::strcmp_ignore_case(c.pruning_method, AlgorithmFactory::ALGORITHM_MIGPP) == 0) {
			if (!isNull(window)) {
				c.window = validateInteger(window, "window");
				if (c.window <= 0) {
					error("The window size, specified in '%s' argument, must be strictly greater than 0.", "window");
				}
			} else {
				is_default_window = true;
			}
		}

		Algorithm* algorithm = NULL;

		try {
			clock_t start_time = 0;
			double execution_time = 0.0;

			Db db;
			const DbView* dbview = NULL;

			Rprintf("Loading data...\n");
			start_time = clock();

			Rprintf("\tPhase file: %s\n", c.phase_file);
			Rprintf("\tMap file: %s\n", c.map_file == NULL ? "NA" : c.map_file);
			if ((c.region[0] > 0) || (c.region[1] != numeric_limits<long int>::max())) {
				Rprintf("\tRegion: [%u, %u]\n", c.region[0], c.region[1]);
			} else {
				Rprintf("\tRegion: NA\n");
			}
			Rprintf("\tMAF filter: > %g\n", c.maf);

			db.load(c.phase_file, c.map_file, c.region[0], c.region[1], c.phase_file_format);
			dbview = db.create_view(c.maf, 0u, numeric_limits<long int>::max());

			c.all_snps = dbview->n_unfiltered_markers;
			c.filtered_snps = dbview->n_markers;
			c.haplotypes = dbview->n_haplotypes;


			Rprintf("\tAll SNPs: %u\n", c.all_snps);
			Rprintf("\tFiltered SNPs: %u\n", c.filtered_snps);
			Rprintf("\tHaplotypes: %u\n", c.haplotypes);
			Rprintf("\tUsed memory (Mb): %.3f\n", db.get_memory_usage());

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Initializing algorithm...\n");
			start_time = clock();

			Rprintf("\tD' CI computation method: %s\n", c.ci_method);
			Rprintf("\tD' likelihood density: ");
			if (auxiliary::strcmp_ignore_case(c.ci_method, CIFactory::CI_WP) == 0) {
				Rprintf("%u\n", c.l_density);
			} else {
				Rprintf("NA\n");
			}
			Rprintf("\tD' CI lower bound for strong LD: >= %g\n", c.ld_ci[0]);
			Rprintf("\tD' CI upper bound for strong LD: >= %g\n", c.ld_ci[1]);
			Rprintf("\tD' CI upper bound for recombination: <= %g\n", c.ehr_ci);
			Rprintf("\tFraction of strong LD SNP pairs: >= %g\n", c.ld_fraction);
			Rprintf("\tPruning method: %s\n", c.pruning_method);
			Rprintf("\tWindow: ");
			if (auxiliary::strcmp_ignore_case(c.pruning_method, AlgorithmFactory::ALGORITHM_MIGPP) == 0) {
				if (is_default_window) {
					c.window = (long int)(((double)dbview->n_markers * (1.0 - c.ld_fraction)) / 2.0);
					if (c.window <= 0) {
						c.window = 1;
					}
				}
				Rprintf("%ld\n", c.window);
			} else {
				Rprintf("NA\n", c.window);
			}

			algorithm = AlgorithmFactory::create(dbview, c.pruning_method);

			algorithm->set_strong_pair_cl(c.ld_ci[0]);
			algorithm->set_strong_pair_cu(c.ld_ci[1]);
			algorithm->set_recomb_pair_cu(c.ehr_ci);
			algorithm->set_strong_pairs_fraction(c.ld_fraction);

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Processing data...\n");
			start_time = clock();

			algorithm->compute_preliminary_blocks(c.ci_method, c.l_density, c.window);
			Rprintf("\tPreliminary haplotype blocks: %u\n", algorithm->get_n_strong_pairs());

			algorithm->sort_preliminary_blocks();
			algorithm->select_final_blocks();

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;

			Rprintf("\tFinal haplotype blocks: %u\n", algorithm->get_n_blocks());

			Rprintf("\tMemory used for preliminary haplotype blocks (Mb): %.3g\n", algorithm->get_memory_usage_preliminary_blocks());
			Rprintf("\tMemory used for final haplotype blocks (Mb): %.3g\n", algorithm->get_memory_usage_final_blocks());
			Rprintf("\tMemory used by algorithm (Mb): %.3g\n", algorithm->get_memory_usage());
			Rprintf("\tTotal used memory (Mb): %.3g\n", algorithm->get_total_memory_usage());
			Rprintf("Done (%.3f sec)\n", execution_time);

			Rprintf("Writing results...\n");
			start_time = clock();

			Rprintf("\tOutput file: %s\n", c.output_file);

			write_info(c.output_file, c);

			algorithm->write_blocks(c.output_file);

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
			Rprintf("Done (%.3f sec)\n", execution_time);

			delete algorithm;
			algorithm = NULL;

		} catch (Exception &e) {
			delete algorithm;
			algorithm = NULL;

			error("%s", e.what());
		}

		return R_NilValue;
	}
}

int main(int args, char** argv) {
	return 0;
}
