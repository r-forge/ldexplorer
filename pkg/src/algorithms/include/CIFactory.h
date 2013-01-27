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

#ifndef CIFACTORY_H_
#define CIFACTORY_H_

#include "../../exception/include/Exception.h"
#include "../../db/include/Db.h"
#include "CI.h"
#include "CIWP.h"
#include "CIAV.h"

class CIFactory {
public:
	static const char* NONE;
	static const char* WP_CI;
	static const char* AV_CI;

	CIFactory();
	virtual ~CIFactory();

	static CI* create(Db& db, const char* method, unsigned int precision = 0u) throw (Exception);
};

#endif
