﻿/*******************************************************************************
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

#ifndef CGOGN_RENDERING_SHADERS_FLAT_COLOR_PER_VERTEX_H_
#define CGOGN_RENDERING_SHADERS_FLAT_COLOR_PER_VERTEX_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(FlatColorPerVertex, CGOGN_STR(FlatColorPerVertex))

class CGOGN_RENDERING_EXPORT ShaderParamFlatColorPerVertex : public ShaderParam
{
    inline void set_uniforms() override
    {
        shader_->set_uniforms_values(ambiant_color_, light_position_, double_side_);
    }

public:
	GLColor ambiant_color_;
	GLVec3 light_position_;
	bool double_side_;


	template<typename ...Args>
	void fill(Args&&... args)
	{
        auto a = std::forward_as_tuple(args...);
        ambiant_color_ = std::get<0>(a);
        light_position_ = std::get<1>(a);
        double_side_ = std::get<2>(a);
    }

    using LocalShader = ShaderFlatColorPerVertex;

    ShaderParamFlatColorPerVertex(LocalShader *sh)
        : ShaderParam(sh)
        , ambiant_color_(0.05f, 0.05f, 0.05f, 1)
        , light_position_(10, 100, 1000)
        , double_side_(true)
    {}

    inline ~ShaderParamFlatColorPerVertex() override {}
};

} // namespace rendering

} // namespace cgogn

#endif
