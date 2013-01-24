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

	SEXP mig(SEXP phase_file, SEXP output_file, SEXP file_format, SEXP legend_file,
			SEXP region, SEXP maf, SEXP ci_method, SEXP ci_precision, SEXP ld_ci, SEXP ehr_ci, SEXP ld_fraction,
			SEXP pruning_method, SEXP window) {

		const char* c_phase_file = NULL;
		const char* c_output_file = NULL;
		const char* c_file_format = NULL;
		const char* c_legend_file = NULL;
		double* c_region_double = NULL;
		long int c_region_int[2] = {0, numeric_limits<long int>::max()};
		double c_maf = 0.0;
		const char* c_ci_method = NULL;
		double c_ci_precision_double = 0.0;
		long int c_ci_precision_int = 0;
		double c_ld_ci[] = {0.0, 0.0};
		double c_ehr_ci = 0.0;
		double c_ld_fraction = 0.0;
		const char* c_pruning_method = NULL;
		double c_window_double = 0.0;
		long int c_window_int = 0;

//		Validate phase_file argument
		if (isNull(phase_file)) {
			cout << "phase_file is null" << endl;
		}

		if (!isString(phase_file)) {
			cout << "phase_file is not string" << endl;
		}

		if (length(phase_file) <= 0) {
			cout << "phase_file is empty" << endl;
		}

		if (length(phase_file) > 1) {
			cout << "phase_file has multiple values." << endl;
		}

		if (STRING_ELT(phase_file, 0) == R_NaString) {
			cout << "phase_file is na" << endl;
		}

		if (STRING_ELT(phase_file, 0) == R_BlankString) {
			cout << "phase_file is blank" << endl;
		}

		c_phase_file = CHAR(STRING_ELT(phase_file, 0));

//		Validate output_file argument
		if (isNull(output_file)) {
			cout << "output_file is null" << endl;
		}

		if (!isString(output_file)) {
			cout << "output_file is not string" << endl;
		}

		if (length(output_file) <= 0) {
			cout << "output_file is empty" << endl;
		}

		if (length(output_file) > 1) {
			cout << "output_file has multiple values." << endl;
		}

		if (STRING_ELT(output_file, 0) == R_NaString) {
			cout << "output_file is na" << endl;
		}

		if (STRING_ELT(output_file, 0) == R_BlankString) {
			cout << "output_file is blank" << endl;
		}

		c_output_file = CHAR(STRING_ELT(output_file, 0));

//		Validate file_format argument
		if (isNull(file_format)) {
			cout << "file_format is null" << endl;
		}

		if (!isString(file_format)) {
			cout << "file_format is not string" << endl;
		}

		if (length(file_format) <= 0) {
			cout << "file_format is empty" << endl;
		}

		if (length(file_format) > 1) {
			cout << "file_format has multiple values." << endl;
		}

		if (STRING_ELT(file_format, 0) == R_NaString) {
			cout << "file_format is na" << endl;
		}

		if (STRING_ELT(file_format, 0) == R_BlankString) {
			cout << "file_format is blank" << endl;
		}

		c_file_format = CHAR(STRING_ELT(file_format, 0));

//		Validate legend_file argument
		if (!isNull(legend_file)) {
			if (!isString(legend_file)) {
				cout << "legend_file is not string" << endl;
			}

			if (length(legend_file) <= 0) {
				cout << "legend_file is empty" << endl;
			}

			if (length(legend_file) > 1) {
				cout << "legend_file has multiple values." << endl;
			}

			if (STRING_ELT(legend_file, 0) == R_NaString) {
				cout << "legend_file is na" << endl;
			}

			if (STRING_ELT(legend_file, 0) == R_BlankString) {
				cout << "legend_file is blank" << endl;
			}

			c_legend_file = CHAR(STRING_ELT(legend_file, 0));
		}

//		Validate region argument
		if (!isNull(region)) {
			if(!isNumeric(region)) {
				cout << "region is not numeric" << endl;
			}

			if (isLogical(region)) {
				cout << "region is logical" << endl;
			}

			if (length(region) != 2) {
				cout << "region must contain 2 values" << endl;
			}

			if (isInteger(region)) {
				c_region_int[0] = (long int)INTEGER(region)[0];
				c_region_int[1] = (long int)INTEGER(region)[1];
			} else {
				c_region_double = REAL(region);
				if (isnan(c_region_double[0])) {
					cout << "region start is nan" << endl;
				}
				if (isnan(c_region_double[1])) {
					cout << "region end is nan" << endl;
				}
				c_region_int[0] = (long int)c_region_double[0];
				c_region_int[1] = (long int)c_region_double[1];
			}

			if (c_region_int[0] < 0) {
				cout << "start must be positive" << endl;
			}

			if (c_region_int[1] < 0) {
				cout << "end must be positive" << endl;
			}

			if (c_region_int[0] >= c_region_int[1]) {
				cout << "end must be strictly greater than end" << endl;
			}

			cout << c_region_int[0] << " " << c_region_int[1] << endl;
		}

//		Validate maf argument
		if (isNull(maf)) {
			cout << "maf is null" << endl;
		}

		if (!isNumeric(maf)) {
			cout << "maf is not numeric" << endl;
		}

		if (isLogical(maf)) {
			cout << "maf is logical" << endl;
		}

		if (length(maf) <= 0) {
			cout << "maf is empty" << endl;
		}

		if (length(maf) > 1) {
			cout << "maf has multiple values" << endl;
		}

		if (isInteger(maf)) {
			c_maf = (double)INTEGER(maf)[0];
		} else {
			c_maf = REAL(maf)[0];
			if (isnan(c_maf)) {
				cout << "c_maf is not a number" << endl;
			}
		}

		if ((c_maf < 0.0) || (c_maf > 0.5)) {
			cout << "c_maf must be in [0, 0.5]" << endl;
		}

//		Validate ci_method argument
		if (isNull(ci_method)) {
			cout << "ci_method is null" << endl;
		}

		if (!isString(ci_method)) {
			cout << "ci_method is not string" << endl;
		}

		if (length(ci_method) <= 0) {
			cout << "ci_method is empty" << endl;
		}

		if (length(ci_method) > 1) {
			cout << "ci_method has multiple values." << endl;
		}

		if (STRING_ELT(ci_method, 0) == R_NaString) {
			cout << "ci_method is na" << endl;
		}

		if (STRING_ELT(ci_method, 0) == R_BlankString) {
			cout << "ci_method is blank" << endl;
		}

		c_ci_method = CHAR(STRING_ELT(ci_method, 0));

//		Validate ci_precision argument
		if (isNull(ci_precision)) {
			cout << "ci_precision is null" << endl;
		}

		if (!isNumeric(ci_precision)) {
			cout << "ci_precision is not numeric" << endl;
		}

		if (isLogical(ci_precision)) {
			cout << "ci_precision is logical" << endl;
		}

		if (length(ci_precision) <= 0) {
			cout << "ci_precision is empty" << endl;
		}

		if (length(ci_precision) > 1) {
			cout << "ci_precision has multiple values." << endl;
		}

		if (isInteger(ci_precision)) {
			c_ci_precision_int = INTEGER(ci_precision)[0];
		} else {
			c_ci_precision_double = REAL(ci_precision)[0];
			if (isnan(c_ci_precision_double)) {
				cout << "ci_precision is nan" << endl;
			}
			c_ci_precision_int = (long int)c_ci_precision_double;
		}

		if (c_ci_precision_int <= 0) {
			cout << "ci_precision must be greater than 0" << endl;
		}

//		Validate ld_ci argument
		if (isNull(ld_ci)) {
			cout << "ld_ci is null" << endl;
		}

		if (!isNumeric(ld_ci)) {
			cout << "ld_ci is not numeric" << endl;
		}

		if (isLogical(ld_ci)) {
			cout << "ld_ci is logical" << endl;
		}

		if (length(ld_ci) != 2) {
			cout << "ld_ci must contain two values" << endl;
		}

		if (isInteger(ld_ci)) {
			c_ld_ci[0] = (double)INTEGER(ld_ci)[0];
			c_ld_ci[1] = (double)INTEGER(ld_ci)[1];
		} else {
			c_ld_ci[0] = REAL(ld_ci)[0];
			if (isnan(c_ld_ci[0])) {
				cout << "CL is nan" << endl;
			}
			c_ld_ci[1] = REAL(ld_ci)[1];
			if (isnan(c_ld_ci[1])) {
				cout << "CU is nan" << endl;
			}
		}

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
		if (isNull(ehr_ci)) {
			cout << "ehr_ci is null" << endl;
		}

		if (!isNumeric(ehr_ci)) {
			cout << "ehr_ci is not numeric" << endl;
		}

		if (isLogical(ehr_ci)) {
			cout << "ehr_ci is logical" << endl;
		}

		if (length(ehr_ci) <= 0) {
			cout << "ehr_ci is empty" << endl;
		}

		if (length(ehr_ci) > 1) {
			cout << "ehr_ci has multiple values" << endl;
		}

		if (isInteger(ehr_ci)) {
			c_ehr_ci = (double)INTEGER(ehr_ci)[0];
		} else {
			c_ehr_ci = REAL(ehr_ci)[0];
			if (isnan(c_ehr_ci)) {
				cout << "c_ehr_ci is not a number" << endl;
			}
		}

		if ((c_ehr_ci < 0.0) || (c_ehr_ci > 1.0)) {
			cout << "c_ehr_ci must be in [0.0, 1.0]" << endl;
		}

//		Validate ld_fraction argument
		if (isNull(ld_fraction)) {
			cout << "ld_fraction is null" << endl;
		}

		if (!isNumeric(ld_fraction)) {
			cout << "ld_fraction is not numeric" << endl;
		}

		if (isLogical(ld_fraction)) {
			cout << "ld_fraction is logical" << endl;
		}

		if (length(ld_fraction) <= 0) {
			cout << "ld_fraction is empty" << endl;
		}

		if (length(ld_fraction) > 1) {
			cout << "ld_fraction has multiple values" << endl;
		}

		if (isInteger(ld_fraction)) {
			c_ld_fraction = (double)INTEGER(ld_fraction)[0];
		} else {
			c_ld_fraction = REAL(ld_fraction)[0];
			if (isnan(c_ld_fraction)) {
				cout << "c_ld_fraction is not a number" << endl;
			}
		}

		if ((c_ld_fraction <= 0.0) || (c_ld_fraction > 1.0)) {
			cout << "c_ld_fraction must be in (0.0, 1.0]" << endl;
		}

//		Validate pruning_method argument
		if (isNull(pruning_method)) {
			cout << "pruning_method is null" << endl;
		}

		if (!isString(pruning_method)) {
			cout << "pruning_method is not string" << endl;
		}

		if (length(pruning_method) <= 0) {
			cout << "pruning_method is empty" << endl;
		}

		if (length(pruning_method) > 1) {
			cout << "pruning_method has multiple values." << endl;
		}

		if (STRING_ELT(pruning_method, 0) == R_NaString) {
			cout << "pruning_method is na" << endl;
		}

		if (STRING_ELT(pruning_method, 0) == R_BlankString) {
			cout << "pruning_method is blank" << endl;
		}

		c_pruning_method = CHAR(STRING_ELT(pruning_method, 0));

//		Validate window argument
		if (isNull(window)) {
			cout << "window is null" << endl;
		}

		if (!isNumeric(window)) {
			cout << "window is not numeric" << endl;
		}

		if (isLogical(window)) {
			cout << "window is logical" << endl;
		}

		if (length(window) <= 0) {
			cout << "window is empty" << endl;
		}

		if (length(window) > 1) {
			cout << "window has multiple values." << endl;
		}

		if (isInteger(window)) {
			c_window_int = INTEGER(window)[0];
		} else {
			c_window_double = REAL(window)[0];
			if (isnan(c_window_double)) {
				cout << "window is nan" << endl;
			}
			c_window_int = (long int)c_window_double;
		}

		if (c_window_int <= 0) {
			cout << "window must be greater than 0" << endl;
		}


		return R_NilValue;
	}
}

int main(int args, char** argv) {
	return 0;
}
