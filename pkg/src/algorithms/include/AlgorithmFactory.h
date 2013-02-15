/*
 * Copyright � 2013 Daniel Taliun, Johann Gamper and Cristian Pattaro. All rights reserved.
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

#ifndef ALGORITHMFACTORY_H_
#define ALGORITHMFACTORY_H_

#include "../../exception/include/Exception.h"
#include "../../db/include/DbView.h"
#include "Algorithm.h"
#include "AlgorithmMIG.h"
#include "AlgorithmMIGP.h"
#include "AlgorithmMIGPP.h"

class AlgorithmFactory {
public:
	AlgorithmFactory();
	virtual ~AlgorithmFactory();

	static Algorithm* create(const char* name, unsigned int window = 0u) throw (Exception);
};

#endif
