/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#ifndef CGOGN_CORE_TYPES_CMAP_CELL_H_
#define CGOGN_CORE_TYPES_CMAP_CELL_H_

#include <cgogn/core/types/cmap/dart.h>

#include <cgogn/core/utils/assert.h>
#include <cgogn/core/utils/numerics.h>

#include <cgogn/geometry/types/vector_traits.h>

/**
 * \file cgogn/core/types/cmap/cell.h
 * \brief Orbit and cell definitions used in cgogn.
 */

namespace cgogn
{

enum Orbit : uint32
{
	DART = 0,
	PHI1,
	PHI2,
	PHI21,
	PHI1_PHI2,
	PHI1_PHI3,
	PHI2_PHI3,
	PHI21_PHI31,
	PHI1_PHI2_PHI3
};

static const std::size_t NB_ORBITS = Orbit::PHI1_PHI2_PHI3 + 1;

inline std::string orbit_name(Orbit orbit)
{
	switch (orbit)
	{
	case Orbit::DART:
		return "cgogn::Orbit::DART";
	case Orbit::PHI1:
		return "cgogn::Orbit::PHI1";
	case Orbit::PHI2:
		return "cgogn::Orbit::PHI2";
	case Orbit::PHI21:
		return "cgogn::Orbit::PHI21";
	case Orbit::PHI1_PHI2:
		return "cgogn::Orbit::PHI1_PHI2";
	case Orbit::PHI1_PHI3:
		return "cgogn::Orbit::PHI1_PHI3";
	case Orbit::PHI2_PHI3:
		return "cgogn::Orbit::PHI2_PHI3";
	case Orbit::PHI21_PHI31:
		return "cgogn::Orbit::PHI21_PHI31";
	case Orbit::PHI1_PHI2_PHI3:
		return "cgogn::Orbit::PHI1_PHI2_PHI3";
		//		default: cgogn_assert_not_reached("This orbit does not exist"); return "UNKNOWN";
	}
	cgogn_assert_not_reached("This orbit does not exist");
#ifdef NDEBUG
	return "UNKNOWN"; // little trick to avoid warning on VS
#endif
}

/**
 * \brief Cellular typing
 * \tparam ORBIT The type of the orbit used to create the Cell
 */
template <Orbit ORBIT_>
struct Cell
{
	static const Orbit ORBIT = ORBIT_;
	using Self = Cell<ORBIT>;

	/**
	 * \brief the dart representing this cell
	 */
	Dart dart;

	/**
	 * \brief Creates a new empty Cell as a nil dart.
	 */
	inline Cell() : dart()
	{
	}

	/**
	 * \brief Creates a new Cell with a dart.
	 * \param[in] d dart to convert to a cell of a given orbit
	 */
	inline explicit Cell(Dart d) : dart(d)
	{
	}

	/**
	 * \brief Copy constructor.
	 * Creates a new Cell from an another one.
	 * \param[in] c a cell
	 */
	inline Cell(const Self& c) : dart(c.dart)
	{
	}

	/**
	 * \brief Tests the validity of the cell.
	 * \retval true if the cell is valid
	 * \retval false otherwise
	 */
	inline bool is_valid() const
	{
		return !dart.is_nil();
	}

	/**
	 * \brief Assigns to the left hand side cell the value
	 * of the right hand side cell.
	 * \param[in] rhs the cell to assign
	 * \return The cell with the assigned value
	 */
	inline Self& operator=(Self rhs)
	{
		dart = rhs.dart;
		return *this;
	}

	/**
	 * \brief Tests equality between two cells
	 * \param[out] boolean value of comparison
	 * \param[in] rhs the cell to compare with
	 */
	inline bool operator==(const Self& rhs) const
	{
		return (this->dart.index == rhs.dart.index);
	}

	/**
	 * \brief Prints a cell to a stream.
	 * \param[out] out the stream to print on
	 * \param[in] rhs the cell to print
	 */
	inline friend std::ostream& operator<<(std::ostream& out, const Self& rhs)
	{
		return out << rhs.dart;
	}

	/**
	 * \brief Reads a cell from a stream.
	 * \param[in] in the stream to read from
	 * \param[out] rhs the cell read
	 */
	inline friend std::istream& operator>>(std::istream& in, Self& rhs)
	{
		in >> rhs.dart;
		return in;
	}
};

class Vec3Ext : public Eigen::Vector3d
{
public:
	Vec3Ext(const Eigen::Vector3d& vec) : Eigen::Vector3d(vec)
	{
	}

	inline bool operator==(const Eigen::Vector3d& vec) const
	{
		return (this->x() == vec.x() && this->y() == vec.y() && this->z() == vec.z());
	}
};

class Vec3ExtHashFunction
{
public:
	size_t operator()(const Vec3Ext& vec) const
	{
		std::hash<int> hasher;
		size_t seed = 0;
		for (int i = 0; i < vec.size(); i++)
		{
			seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CELL_H_
