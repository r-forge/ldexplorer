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

#include "include/auxiliary.h"

char* auxiliary::strtok(char** start, char separator) {
	if ((*start != NULL) && (**start != '\0')) {
		char* token = *start;

		do {
			if (**start == separator) {
				**start = '\0';
				(*start)++;
				return token;
			}
			(*start)++;
		} while (**start != '\0');

		return token;
	}

	return NULL;
}

int auxiliary::strcmp_ignore_case(const char* first, const char* second) {
	int i = 0;

	while ((first[i] != '\0') && (second[i] != '\0')) {
		if (tolower(first[i]) < tolower(second[i])) {
			return -1;
		}
		else if (tolower(first[i]) > tolower(second[i])) {
			return 1;
		}
		i++;
	}

	if (first[i] == '\0') {
		if (second[i] == '\0') {
			return 0;
		}
		return -1;
	}
	return 1;
}

int auxiliary::strcmp_ignore_case(const char* first, const char* second, int n) {
	int i = 0;

	while ((first[i] != '\0') && (second[i] != '\0') && (i < n)) {
		if (tolower(first[i]) < tolower(second[i])) {
			return -1;
		}
		else if (tolower(first[i]) > tolower(second[i])) {
			return 1;
		}
		i++;
	}

	if (i >= n) {
		return 0;
	}

	if (first[i] == '\0') {
		if (second[i] == '\0') {
			return 0;
		}
		return -1;
	}

	return 1;
}

double auxiliary::stats_quantile_from_sorted_data(double* data, unsigned int size, double fraction) {
	unsigned int i = (unsigned int)floor((size - 1) * fraction);
	double delta = (size - 1) * fraction - i;

	return (1 - delta) * data[i] + delta * data[i + 1];
}
