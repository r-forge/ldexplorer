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

#include <iostream>
#include <math.h>

#include "algorithms/include/CIFactory.h"
#include "algorithms/include/AlgorithmFactory.h"
#include "db/include/Db.h"

#include <R.h>
#include <Rinternals.h>

using namespace std;

extern "C" {

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

	SEXP mig(SEXP phase_file, SEXP output_file, SEXP file_format, SEXP legend_file,
			SEXP region, SEXP maf, SEXP ci_method, SEXP ci_precision, SEXP ld_ci, SEXP ehr_ci, SEXP ld_fraction,
			SEXP pruning_method, SEXP window) {

		const char* c_phase_file = NULL;
		const char* c_output_file = NULL;
		const char* c_file_format = NULL;
		const char* c_legend_file = NULL;
		long int c_region[2] = {0, numeric_limits<long int>::max()};
		bool is_region = false;
		double c_maf = 0.0;
		const char* c_ci_method = NULL;
		long int c_ci_precision = 0;
		double c_ld_ci[] = {0.0, 0.0};
		double c_ehr_ci = 0.0;
		double c_ld_fraction = 0.0;
		const char* c_pruning_method = NULL;
		long int c_window = numeric_limits<long int>::max();
		bool is_default_window = false;

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
		if (!isNull(file_format)) {
			c_file_format = validateString(file_format, "file_format");
			if ((auxiliary::strcmp_ignore_case(c_file_format, Db::VCF) != 0) &&
					(auxiliary::strcmp_ignore_case(c_file_format, Db::HAPMAP2) != 0)) {
				error("The file format, specified in '%s' argument, must be '%s' or '%s'.", "file_format", Db::VCF, Db::HAPMAP2);
			}
		} else {
			error("'%s' argument is NULL.", "file_format");
		}

//		Validate legend_file argument.
		if (!isNull(legend_file)) {
			c_legend_file = validateString(legend_file, "legend_file");
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
			is_region = true;
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
			if ((auxiliary::strcmp_ignore_case(c_ci_method, CIFactory::CI_WP) != 0) &&
					(auxiliary::strcmp_ignore_case(c_ci_method, CIFactory::CI_AV) != 0)) {
				error("The method to compute the confidence interval (CI) of D', specified in '%s' argument, must be '%s' or '%s'.", "file_format", CIFactory::CI_WP, CIFactory::CI_AV);
			}
		} else {
			error("'%s' argument is NULL.", "ci_method");
		}

//		Validate ci_precision argument if WP method to compute D' CI was specified.
		if (auxiliary::strcmp_ignore_case(c_ci_method, CIFactory::CI_WP) == 0) {
			if (!isNull(ci_precision)) {
				c_ci_precision = validateInteger(ci_precision, "ci_precision");
				if (c_ci_precision <= 0) {
					error("The number of likelihood estimation points to compute confidence interval, specified in '%s' argument, must be strictly greater then 0.", "ci_precision");
				}
			} else {
				error("'%s' argument is NULL.", "ci_precision");
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
			if ((auxiliary::strcmp_ignore_case(c_pruning_method, AlgorithmFactory::ALGORITHM_MIG) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, AlgorithmFactory::ALGORITHM_MIGP) != 0) &&
					(auxiliary::strcmp_ignore_case(c_pruning_method, AlgorithmFactory::ALGORITHM_MIGPP) != 0)) {
				error("The search space pruning method, specified in '%s' argument, must be '%s', '%s' or '%s'.",
						"file_format", AlgorithmFactory::ALGORITHM_MIG, AlgorithmFactory::ALGORITHM_MIGP, AlgorithmFactory::ALGORITHM_MIGPP);
			}
		} else {
			error("'%s' argument is NULL.", "pruning_method");
		}

//		Validate window argument if MIG++ search space pruning method was specified.
		if (auxiliary::strcmp_ignore_case(c_pruning_method, AlgorithmFactory::ALGORITHM_MIGPP) == 0) {
			if (!isNull(window)) {
				c_window = validateInteger(window, "window");
				if (c_window <= 0) {
					error("The window size, specified in '%s' argument, must be strictly greater than 0.", "window");
				}
			} else {
				is_default_window = true;
			}
		}

		try {
			clock_t start_time = 0;
			double execution_time = 0.0;

			Db db;
			Algorithm* algorithm = NULL;

			Rprintf("Loading data...\n");
			start_time = clock();

			Rprintf("\tPhase file: %s\n", c_phase_file);
			Rprintf("\tLegend file: %s\n", c_legend_file == NULL ? "NA" : c_legend_file);
			Rprintf("\tRegion: ");
			if (is_region) {
				Rprintf("[%u, %u]\n", c_region[0], c_region[1]);
			} else {
				Rprintf("NA\n");
			}
			Rprintf("\tMAF filter: > %g\n", c_maf);

			if (is_region) {
				if (auxiliary::strcmp_ignore_case(c_file_format, Db::VCF) == 0) {
					db.load_vcf(c_phase_file, c_region[0], c_region[1]);
				} else if (auxiliary::strcmp_ignore_case(c_file_format, Db::HAPMAP2) == 0) {
					db.load_hapmap2(c_legend_file, c_phase_file, c_region[0], c_region[1]);
				}
			} else {
				if (auxiliary::strcmp_ignore_case(c_file_format, Db::VCF) == 0) {
					db.load_vcf(c_phase_file);
				} else if (auxiliary::strcmp_ignore_case(c_file_format, Db::HAPMAP2) == 0) {
					db.load_hapmap2(c_legend_file, c_phase_file);
				}
			}

			db.mask(c_maf);

			Rprintf("\tAll SNPs: %u\n", db.get_all_n_markers());
			Rprintf("\tFiltered SNPs: %u\n", db.get_n_markers());
			Rprintf("\tHaplotypes: %u\n", db.get_n_haplotypes());
			Rprintf("\tUsed memory (Mb): %.3g\n", db.get_memory_usage());

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
			Rprintf("Done (%.3g sec)\n", execution_time);

			Rprintf("Processing data...\n");
			start_time = clock();

			Rprintf("\tD' CI computation method: %s\n", c_ci_method);
			Rprintf("\tD' CI precision: ");
			if (auxiliary::strcmp_ignore_case(c_ci_method, CIFactory::CI_WP) == 0) {
				Rprintf("%u\n", c_ci_precision);
			} else {
				Rprintf("NA\n");
			}
			Rprintf("\tD' CI lower bound for strong LD: >= %g\n", c_ld_ci[0]);
			Rprintf("\tD' CI upper bound for strong LD: >= %g\n", c_ld_ci[1]);
			Rprintf("\tD' CI upper bound for recombination: <= %g\n", c_ehr_ci);
			Rprintf("\tFraction of strong LD SNP pairs: >= %g\n", c_ld_fraction);
			Rprintf("\tPruning method: %s\n", c_pruning_method);
			Rprintf("\tWindow: ");
			if (auxiliary::strcmp_ignore_case(c_pruning_method, AlgorithmFactory::ALGORITHM_MIGPP) == 0) {
				if (is_default_window) {
					c_window = (long int)(((double)db.get_n_markers() * (1.0 - c_ld_fraction)) / 2.0);
				}
				Rprintf("%ld\n", c_window);
			} else {
				Rprintf("NA\n", c_window);
			}

			algorithm = AlgorithmFactory::create(db, c_pruning_method);
			algorithm->compute_preliminary_blocks(c_ci_method, c_ci_precision, c_window);

			execution_time = (clock() - start_time)/(double)CLOCKS_PER_SEC;
			Rprintf("Done (%.3g sec)\n", execution_time);

		} catch (Exception &e) {
			error("%s", e.what());
		}

		return R_NilValue;
	}
}

int main(int args, char** argv) {
	return 0;
}
