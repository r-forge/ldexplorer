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

#include <Rinternals.h>

using namespace std;

extern "C" {

	const char* validateString(SEXP value, const char* name, bool optional) {
		if (isNull(value)) {
			if (!optional) {
				error("%s is NULL.", name);
			} else {
				return NULL;
			}
		}

		if (!isString(value)) {
			error("%s is not a string.", name);
		}

		if (length(value) <= 0) {
			error("%s contains no values.", name);
		}

		if (length(value) > 1) {
			error("%s contains multiple values.", name);
		}

		if (STRING_ELT(value, 0) == R_NaString) {
			error("%s is NA.", name);
		}

		if (STRING_ELT(value, 0) == R_BlankString) {
			error("%s is blank.", name);
		}

		return CHAR(STRING_ELT(value, 0));
	}

	double validateDouble(SEXP value, const char* name, bool optional) {
		double c_value = numeric_limits<double>::quiet_NaN();

		if (isNull(value)) {
			if (!optional) {
				error("%s is NULL.", name);
			} else {
				return numeric_limits<double>::quiet_NaN();
			}
		}

		if (!isNumeric(value)) {
			error("%s is not numeric.", name);
		}

		if (isLogical(value)) {
			error("%s is logical.", name);
		}

		if (length(value) <= 0) {
			error("%s contains no values.", name);
		}

		if (length(value) > 1) {
			error("%s contains multiple values.", name);
		}

		if (isInteger(value)) {
			c_value = (double)INTEGER(value)[0];
		} else {
			c_value = REAL(value)[0];
			if (isnan(c_value)) {
				error("%s is NA/NaN.", name);
			}
		}

		return c_value;
	}

	void validateDoubles(SEXP value, const char* name, double* c_value, unsigned int length, bool optional) {
		if (isNull(value)) {
			if (!optional) {
				error("%s is NULL.", name);
			} else {
				return;
			}
		}

		if (!isNumeric(value)) {
			error("%s is not numeric.", name);
		}

		if (isLogical(value)) {
			error("%s is logical.", name);
		}

		if (length(value) < (long int)length) {
			error("%s contains less than %u value(s).", name, length);
		}

		if (length(value) > (long int)length) {
			error("%s contains more than %u value(s).", name, length);
		}

		if (isInteger(value)) {
			for (unsigned int i = 0; i < length; ++i) {
				c_value[i] = (double)INTEGER(value)[i];
			}
		} else {
			for (unsigned int i = 0; i < length; ++i) {
				c_value[i] = REAL(value)[i];
				if (isnan(c_value[i])) {
					error("%s contains NA/NaN value(s).", name);
				}
			}
		}
	}

	long int validateInteger(SEXP value, const char* name, bool optional) {
		double c_value_double = numeric_limits<double>::quiet_NaN();
		long int c_value_int = numeric_limits<long int>::min();

		if (isNull(value)) {
			if (!optional) {
				error("%s is NULL.", name);
			} else {
				return numeric_limits<long int>::min();
			}
		}

		if (!isNumeric(value)) {
			error("%s is not numeric.", name);
		}

		if (isLogical(value)) {
			error("%s is logical.", name);
		}

		if (length(value) <= 0) {
			error("%s contains no values.", name);
		}

		if (length(value) > 1) {
			error("%s contains multiple values.", name);
		}

		if (isInteger(value)) {
			c_value_int = INTEGER(value)[0];
		} else {
			c_value_double = REAL(value)[0];
			if (isnan(c_value_double)) {
				error("%s is NA/NaN.", name);
			}
			c_value_int = (long int)c_value_double;
		}

		return c_value_int;
	}

	void validateIntegers(SEXP value, const char* name, long int* c_value, unsigned int length, bool optional) {
		double c_value_double = numeric_limits<double>::quiet_NaN();

		if (isNull(value)) {
			if (!optional) {
				error("%s is NULL.", name);
			} else {
				return;
			}
		}

		if (!isNumeric(value)) {
			error("%s is not numeric.", name);
		}

		if (isLogical(value)) {
			error("%s is logical.", name);
		}

		if (length(value) < (long int)length) {
			error("%s contains less than %u value(s).", name, length);
		}

		if (length(value) > (long int)length) {
			error("%s contains more than %u value(s).", name, length);
		}

		if (isInteger(value)) {
			for (unsigned int i = 0; i < length; ++i) {
				c_value[i] = (double)INTEGER(value)[i];
			}
		} else {
			for (unsigned int i = 0; i < length; ++i) {
				c_value_double = REAL(value)[i];
				if (isnan(c_value_double)) {
					error("%s contains NA/NaN value(s).", name);
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
		double c_maf = 0.0;
		const char* c_ci_method = NULL;
		long int c_ci_precision = 0;
		double c_ld_ci[] = {0.0, 0.0};
		double c_ehr_ci = 0.0;
		double c_ld_fraction = 0.0;
		const char* c_pruning_method = NULL;
		long int c_window = 0;

//		Validate phase_file argument
		c_phase_file = validateString(phase_file, "phase_file", false);

//		Validate output_file argument
		c_output_file = validateString(output_file, "output_file", false);

//		Validate file_format argument
		c_file_format = validateString(file_format, "file_format", false);

//		Validate legend_file argument
		c_legend_file = validateString(legend_file, "legend_file", true);

//		Validate region argument
		validateIntegers(region, "region", c_region, 2u, true);
		if (c_region[0] < 0) {
			cout << "start must be positive" << endl;
		}
		if (c_region[1] < 0) {
			cout << "end must be positive" << endl;
		}
		if (c_region[0] >= c_region[1]) {
			cout << "end must be strictly greater than end" << endl;
		}

//		Validate maf argument
		c_maf = validateDouble(maf, "maf", false);
		if ((c_maf < 0.0) || (c_maf > 0.5)) {
			cout << "c_maf must be in [0, 0.5]" << endl;
		}

//		Validate ci_method argument
		c_ci_method = validateString(ci_method, "ci_method", false);

//		Validate ci_precision argument
		c_ci_precision = validateInteger(ci_precision, "ci_precision", false);
		if (c_ci_precision <= 0) {
			cout << "ci_precision must be greater than 0" << endl;
		}

//		Validate ld_ci argument
		validateDoubles(ld_ci, "ld_ci", c_ld_ci, 2u, false);
		if ((c_ld_ci[0] < 0.0) || (c_ld_ci[0] > 1.0)) {
			cout << "CL must be in [0.0, 1.0]" << endl;
		}

		if ((c_ld_ci[1] < 0.0) || (c_ld_ci[1] > 1.0)) {
			cout << "CU must be in [0.0, 1.0]" << endl;
		}

		if (c_ld_ci[0] >= c_ld_ci[1]) {
			cout << "CL must be greater than CU" << endl;
		}

//		Validate ehr_ci argument
		c_ehr_ci = validateDouble(ehr_ci, "ehr_ci", false);
		if ((c_ehr_ci < 0.0) || (c_ehr_ci > 1.0)) {
			cout << "c_ehr_ci must be in [0.0, 1.0]" << endl;
		}

//		Validate ld_fraction argument
		c_ld_fraction = validateDouble(ld_fraction, "ld_fraction", false);
		if ((c_ld_fraction <= 0.0) || (c_ld_fraction > 1.0)) {
			cout << "c_ld_fraction must be in (0.0, 1.0]" << endl;
		}

//		Validate pruning_method argument
		c_pruning_method = validateString(pruning_method, c_pruning_method, false);

//		Validate window argument
		c_window = validateInteger(window, "window", true);
		if (c_window <= 0) {
			cout << "window must be greater than 0" << endl;
		}

		return R_NilValue;
	}
}

int main(int args, char** argv) {
	return 0;
}
