/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_SURFACE_TOOLS_H_
#define CGOGN_MODULE_SURFACE_TOOLS_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/selection.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>
#include <unordered_set>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceTools : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceTools can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	enum SelectingCell
	{
		VertexSelect = 0,
		EdgeSelect,
		FaceSelect
	};

	enum SelectionMethod
	{
		SingleCell = 0,
		WithinSphere
	};

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  selected_vertices_set_(nullptr), selected_edges_set_(nullptr), selected_faces_set_(nullptr),
			  selecting_cell_(VertexSelect), selection_method_(SingleCell)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			param_edge_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_edge_->width_ = 10.0f;
			param_edge_->set_vbos({&selected_edges_vbo_});

			param_flat_ = rendering::ShaderFlat::generate_param();
			param_flat_->front_color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_flat_->back_color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_flat_->double_side_ = true;
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_flat_->set_vbos({&selected_faces_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

	public:
		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> temp_position;
				temp_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell(
					[&](Vertex v) { temp_position.push_back(value<Vec3>(*mesh_, vertex_position_, v)); });
				rendering::update_vbo(temp_position, &selected_vertices_vbo_);
				selected_vertices_position_.swap(temp_position);
			}
		}

		void update_selected_edges_vbo()
		{
			if (selected_edges_set_)
			{
				std::vector<Vec3> temp_position;
				temp_position.reserve(selected_edges_set_->size() * 2);
				selected_edges_set_->foreach_cell([&](Edge e) {
					std::vector<Vertex> vertices = incident_vertices(*mesh_, e);
					temp_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[0]));
					temp_position.push_back(value<Vec3>(*mesh_, vertex_position_, vertices[1]));
				});
				rendering::update_vbo(temp_position, &selected_edges_vbo_);
				selected_edges_position_.swap(temp_position);
			}
		}

		void update_selected_faces_vbo()
		{
			if (selected_faces_set_)
			{
				std::vector<Vec3> temp_position;
				temp_position.reserve(selected_faces_set_->size() * 3); // TODO: manage polygonal faces
				selected_faces_set_->foreach_cell([&](Face f) {
					foreach_incident_vertex(*mesh_, f, [&](Vertex v) -> bool {
						temp_position.push_back(value<Vec3>(*mesh_, vertex_position_, v));
						return true;
					});
				});
				rendering::update_vbo(temp_position, &selected_faces_vbo_);
				selected_faces_position_.swap(temp_position);
			}
		}

		MESH* mesh_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO selected_vertices_vbo_;
		rendering::VBO selected_edges_vbo_;
		rendering::VBO selected_faces_vbo_;

		CellsSet<MESH, Vertex>* selected_vertices_set_;
		CellsSet<MESH, Edge>* selected_edges_set_;
		CellsSet<MESH, Face>* selected_faces_set_;

		std::vector<Vec3> selected_vertices_position_;
		std::vector<Vec3> selected_edges_position_;
		std::vector<Vec3> selected_faces_position_;

		SelectingCell selecting_cell_;
		SelectionMethod selection_method_;
	};

public:
	SurfaceTools(const App& app)
		: ViewModule(app, "SurfaceTools (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr)
	{
	}

	~SurfaceTools()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](Attribute<Vec3>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6);
						p.update_selected_vertices_vbo();
						p.update_selected_edges_vbo();
						p.update_selected_faces_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Vertex>>(
				m, [this, m](CellsSet<MESH, Vertex>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Edge>>(
				m, [this, m](CellsSet<MESH, Edge>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_edges_set_ == set && p.vertex_position_)
					{
						p.update_selected_edges_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Face>>(
				m, [this, m](CellsSet<MESH, Face>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_faces_set_ == set && p.vertex_position_)
					{
						p.update_selected_faces_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
	}

	void set_position(std::vector<Vec3> vec, float array[3])
	{
		array[0] = 0.0f;
		array[1] = 0.0f;
		array[2] = 0.0f;

		for (Vec3 value : vec)
		{
			array[0] += (float)value.x();
			array[1] += (float)value.y();
			array[2] += (float)value.z();
		}

		int size = vec.size();
		array[0] /= size;
		array[1] /= size;
		array[2] /= size;
	}

	void translate(const MESH& m, float old_position[3], float position[3])
	{
		if (old_position[0] == position[0] && old_position[1] == position[1] && old_position[2] == position[2])
			return;

		Parameters& p = parameters_[&m];

		Vec3 translation(position[0] - old_position[0], position[1] - old_position[1], position[2] - old_position[2]);

		if (p.selecting_cell_ == VertexSelect)
		{
			p.selected_vertices_set_->foreach_cell(
				[&](Vertex v) { value<Vec3>(*p.mesh_, p.vertex_position_, v) += translation; });
		}
		else if (p.selecting_cell_ == EdgeSelect)
		{
			std::unordered_set<Vec3Ext, Vec3ExtHashFunction> to_translate_set;
			std::vector<Vertex> to_translate_vector;
			p.selected_edges_set_->foreach_cell([&](Edge e) {
				std::vector<Vertex> vertices = incident_vertices(*p.mesh_, e);
				for (Vertex v : vertices)
				{
					Vec3Ext tempVector = (Vec3Ext)value<Vec3>(*p.mesh_, p.vertex_position_, v);
					if (to_translate_set.insert(tempVector).second)
					{
						to_translate_vector.push_back(v);
					}
				}
			});

			for (Vertex v : to_translate_vector)
			{
				value<Vec3>(*p.mesh_, p.vertex_position_, v) += translation;
			}
		}
		else if (p.selecting_cell_ == FaceSelect)
		{
			std::unordered_set<Vec3Ext, Vec3ExtHashFunction> to_translate_set;
			std::vector<Vertex> to_translate_vector;
			p.selected_faces_set_->foreach_cell([&](Face f) {
				foreach_incident_vertex(*p.mesh_, f, [&](Vertex v) -> bool {
					Vec3Ext tempVector = (Vec3Ext)value<Vec3>(*p.mesh_, p.vertex_position_, v);
					if (to_translate_set.insert(tempVector).second)
					{
						to_translate_vector.push_back(v);
					}
					return true;
				});
			});

			for (Vertex v : to_translate_vector)
			{
				value<Vec3>(*p.mesh_, p.vertex_position_, v) += translation;
			}
		}

		mesh_provider_->emit_attribute_changed(selected_mesh_, p.vertex_position_.get());
	}

	void extrude(MESH& m)
	{
		Parameters& p = parameters_[&m];

		p.selected_faces_set_->foreach_cell([&](Face f) {
			Dart base[3];
			uint8 it = 0u;

			Dart d = f.dart;
			do
			{
				base[it++] = phi2(m, d);
				phi2_unsew(m, d);
				d = phi1(m, d);
			} while (d != f.dart);

			Face extruded_face = add_face(m, 3, true);
			p.selected_faces_set_->select(extruded_face);
			foreach_incident_vertex(*p.mesh_, extruded_face, [&](Vertex v) -> bool {
				phi2_unsew(m, v.dart);
				return true;
			});

			Vec3 abc[3];
			it = 0u;
			foreach_incident_vertex(*p.mesh_, f, [&](Vertex v) -> bool {
				abc[it++] = value<Vec3>(*p.mesh_, p.vertex_position_, v);
				return true;
			});

			Vec3 ab = abc[1] - abc[0];
			Vec3 ac = abc[2] - abc[0];
			Vec3 normal = ab.cross(ac);
			float height = normal.norm() / 20.0f;
			normal.normalize();

			it = 0u;
			foreach_incident_vertex(*p.mesh_, extruded_face, [&](Vertex v) -> bool {
				value<Vec3>(*p.mesh_, p.vertex_position_, v) = abc[it++] + normal * height;
				return true;
			});

			remove_face(static_cast<CMap1&>(m), CMap1::Face(f), false);

			Face faces_ring[6];
			for (int i = 0; i < 6; i++)
			{
				faces_ring[i] = add_face(m, 3, true);
				Dart tempd = faces_ring[i].dart;
				do
				{
					phi2_unsew(m, tempd);
					tempd = phi1(m, tempd);
				} while (tempd != faces_ring[i].dart);
			}

			set_index<Vertex>(m, faces_ring[0].dart, index_of(m, Vertex(base[0])));
			set_index<Vertex>(m, phi1(m, faces_ring[0].dart), index_of(m, Vertex(base[1])));
			set_index<Vertex>(m, phi1(m, phi1(m, faces_ring[0].dart)),
							  index_of(m, Vertex(phi1(m, phi1(m, extruded_face.dart)))));

			set_index<Vertex>(m, faces_ring[1].dart, index_of(m, Vertex(phi1(m, phi1(m, extruded_face.dart)))));
			set_index<Vertex>(m, phi1(m, faces_ring[1].dart), index_of(m, Vertex(phi1(m, extruded_face.dart))));
			set_index<Vertex>(m, phi1(m, phi1(m, faces_ring[1].dart)), index_of(m, Vertex(base[0])));

			set_index<Vertex>(m, faces_ring[2].dart, index_of(m, Vertex(base[1])));
			set_index<Vertex>(m, phi1(m, faces_ring[2].dart), index_of(m, Vertex(base[2])));
			set_index<Vertex>(m, phi1(m, phi1(m, faces_ring[2].dart)), index_of(m, Vertex(extruded_face.dart)));

			set_index<Vertex>(m, faces_ring[3].dart, index_of(m, Vertex(extruded_face.dart)));
			set_index<Vertex>(m, phi1(m, faces_ring[3].dart),
							  index_of(m, Vertex(phi1(m, phi1(m, extruded_face.dart)))));
			set_index<Vertex>(m, phi1(m, phi1(m, faces_ring[3].dart)), index_of(m, Vertex(base[1])));

			set_index<Vertex>(m, faces_ring[4].dart, index_of(m, Vertex(base[2])));
			set_index<Vertex>(m, phi1(m, faces_ring[4].dart), index_of(m, Vertex(base[0])));
			set_index<Vertex>(m, phi1(m, phi1(m, faces_ring[4].dart)),
							  index_of(m, Vertex(phi1(m, extruded_face.dart))));

			set_index<Vertex>(m, faces_ring[5].dart, index_of(m, Vertex(phi1(m, extruded_face.dart))));
			set_index<Vertex>(m, phi1(m, faces_ring[5].dart), index_of(m, Vertex(extruded_face.dart)));
			set_index<Vertex>(m, phi1(m, phi1(m, faces_ring[5].dart)), index_of(m, Vertex(base[2])));

			phi2_sew(m, phi1(m, faces_ring[0].dart), phi1(m, faces_ring[3].dart));
			phi2_sew(m, phi1(m, faces_ring[2].dart), phi1(m, faces_ring[5].dart));
			phi2_sew(m, phi1(m, faces_ring[4].dart), phi1(m, faces_ring[1].dart));

			phi2_sew(m, phi1(m, phi1(m, faces_ring[0].dart)), phi1(m, phi1(m, faces_ring[1].dart)));
			phi2_sew(m, phi1(m, phi1(m, faces_ring[2].dart)), phi1(m, phi1(m, faces_ring[3].dart)));
			phi2_sew(m, phi1(m, phi1(m, faces_ring[4].dart)), phi1(m, phi1(m, faces_ring[5].dart)));
			
			phi2_sew(m, base[1], faces_ring[0].dart);
			phi2_sew(m, base[2], faces_ring[2].dart);
			phi2_sew(m, base[0], faces_ring[4].dart);

			phi2_sew(m, phi1(m, extruded_face.dart), faces_ring[1].dart);
			phi2_sew(m, phi1(m, phi1(m, extruded_face.dart)), faces_ring[3].dart);
			phi2_sew(m, extruded_face.dart, faces_ring[5].dart);
		});

		mesh_provider_->emit_attribute_changed(selected_mesh_, p.vertex_position_.get());
		mesh_provider_->emit_connectivity_changed(selected_mesh_);
		mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_faces_set_);
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, p.vertex_position_.get()) / 6); // 6 ???
			p.update_selected_vertices_vbo();
			p.update_selected_edges_vbo();
			p.update_selected_faces_vbo();
		}

		for (View* v : linked_views_)
			v->request_update();
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SurfaceTools<MESH>::init_mesh));
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (selected_mesh_ && view->shift_pressed())
		{
			// MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);
			Parameters& p = parameters_[selected_mesh_];

			if (p.vertex_position_)
			{
				rendering::GLVec3d near = view->unproject(x, y, 0.0);
				rendering::GLVec3d far = view->unproject(x, y, 1.0);
				Vec3 A{near.x(), near.y(), near.z()};
				Vec3 B{far.x(), far.y(), far.z()};

				switch (p.selection_method_)
				{
				case SingleCell: {
					switch (p.selecting_cell_)
					{
					case VertexSelect:
						if (p.selected_vertices_set_)
						{
							std::vector<Vertex> picked;
							cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
							if (!picked.empty())
							{
								switch (button)
								{
								case 0:
									p.selected_vertices_set_->select(picked[0]);
									break;
								case 1:
									p.selected_vertices_set_->unselect(picked[0]);
									break;
								}
								mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_vertices_set_);
							}
						}
						break;
					case EdgeSelect:
						if (p.selected_edges_set_)
						{
							std::vector<Edge> picked;
							cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
							if (!picked.empty())
							{
								switch (button)
								{
								case 0:
									p.selected_edges_set_->select(picked[0]);
									break;
								case 1:
									p.selected_edges_set_->unselect(picked[0]);
									break;
								}
								mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_edges_set_);
							}
						}
						break;
					case FaceSelect:
						if (p.selected_faces_set_)
						{
							std::vector<Face> picked;
							cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
							if (!picked.empty())
							{
								switch (button)
								{
								case 0:
									p.selected_faces_set_->select(picked[0]);
									break;
								case 1:
									p.selected_faces_set_->unselect(picked[0]);
									break;
								}
								mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_faces_set_);
							}
						}
						break;
					}
					break;
				}
				case WithinSphere: {
					std::vector<Vertex> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
					if (!picked.empty())
					{
						CellCache<MESH> cache = geometry::within_sphere(*selected_mesh_, picked[0],
																		p.vertex_base_size_ * p.sphere_scale_factor_,
																		p.vertex_position_.get());
						switch (p.selecting_cell_)
						{
						case VertexSelect:
							if (p.selected_vertices_set_)
							{
								switch (button)
								{
								case 0:
									foreach_cell(cache, [&p](Vertex v) -> bool {
										p.selected_vertices_set_->select(v);
										return true;
									});
									break;
								case 1:
									foreach_cell(cache, [&p](Vertex v) -> bool {
										p.selected_vertices_set_->unselect(v);
										return true;
									});
									break;
								}
								mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_vertices_set_);
							}
							break;
						case EdgeSelect:
							if (p.selected_edges_set_)
							{
								switch (button)
								{
								case 0:
									foreach_cell(cache, [&p](Edge e) -> bool {
										p.selected_edges_set_->select(e);
										return true;
									});
									break;
								case 1:
									foreach_cell(cache, [&p](Edge e) -> bool {
										p.selected_edges_set_->unselect(e);
										return true;
									});
									break;
								}
								mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_edges_set_);
							}
							break;
						case FaceSelect:
							if (p.selected_faces_set_)
							{
								switch (button)
								{
								case 0:
									foreach_cell(cache, [&p](Face f) -> bool {
										p.selected_faces_set_->select(f);
										return true;
									});
									break;
								case 1:
									foreach_cell(cache, [&p](Face f) -> bool {
										p.selected_faces_set_->unselect(f);
										return true;
									});
									break;
								}
								mesh_provider_->emit_cells_set_changed(selected_mesh_, p.selected_faces_set_);
							}
							break;
						}
					}
				}
				}
			}
		}
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.selecting_cell_ == VertexSelect && p.selected_vertices_set_ && p.selected_vertices_set_->size() > 0 &&
				p.param_point_sprite_->attributes_initialized())
			{
				p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_POINTS, 0, p.selected_vertices_set_->size());
				p.param_point_sprite_->release();
			}
			else if (p.selecting_cell_ == EdgeSelect && p.selected_edges_set_ && p.selected_edges_set_->size() > 0 &&
					 p.param_edge_->attributes_initialized())
			{
				p.param_edge_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_LINES, 0, p.selected_edges_set_->size() * 2);
				p.param_edge_->release();
			}
			else if (p.selecting_cell_ == FaceSelect && p.selected_faces_set_ && p.selected_faces_set_->size() > 0 &&
					 p.param_flat_->attributes_initialized())
			{
				p.param_flat_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_TRIANGLES, 0, p.selected_faces_set_->size() * 3); // TODO: manage polygonal faces
				p.param_flat_->release();
			}
		}
	}

	void interface() override
	{
		bool need_update = false;

		imgui_mesh_selector(mesh_provider_, selected_mesh_, [&](MESH* m) {
			selected_mesh_ = m;
			mesh_provider_->mesh_data(selected_mesh_)->outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;
			Parameters& p = parameters_[selected_mesh_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_mesh_, attribute);
												});

			if (p.vertex_position_)
			{
				ImGui::Separator();
				int* ptr_sel_cell = reinterpret_cast<int*>(&p.selecting_cell_);
				need_update |= ImGui::RadioButton("Vertex", ptr_sel_cell, VertexSelect);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Edge", ptr_sel_cell, EdgeSelect);
				ImGui::SameLine();
				need_update |= ImGui::RadioButton("Face", ptr_sel_cell, FaceSelect);

				ImGui::RadioButton("Single", reinterpret_cast<int*>(&p.selection_method_), SingleCell);
				ImGui::SameLine();
				ImGui::RadioButton("Sphere", reinterpret_cast<int*>(&p.selection_method_), WithinSphere);

				if (p.selection_method_ == WithinSphere)
					ImGui::SliderFloat("Sphere radius", &(p.sphere_scale_factor_), 10.0f, 100.0f);

				MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);

				if (p.selecting_cell_ == VertexSelect)
				{
					if (ImGui::BeginCombo("Sets", p.selected_vertices_set_ ? p.selected_vertices_set_->name().c_str()
																		   : "-- select --"))
					{
						md->template foreach_cells_set<Vertex>([&](CellsSet<MESH, Vertex>& cs) {
							bool is_selected = &cs == p.selected_vertices_set_;
							if (ImGui::Selectable(cs.name().c_str(), is_selected))
							{
								p.selected_vertices_set_ = &cs;
								p.update_selected_vertices_vbo();
								for (View* v : linked_views_)
									v->request_update();
							}
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						});
						ImGui::EndCombo();
					}
					if (p.selected_vertices_set_)
					{
						ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
						if (ImGui::Button("X##selected_vertices_set"))
							p.selected_vertices_set_ = nullptr;
					}
					if (ImGui::Button("Create##vertices_set"))
						md->template add_cells_set<Vertex>();
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);

					if (p.selected_vertices_set_ && !p.selected_vertices_set_->empty())
					{
						ImGui::Separator();
						ImGui::TextUnformatted("Vertex Tools");

						float position[3];
						set_position(p.selected_vertices_position_, position);

						float old_position[3];
						std::copy(std::begin(position), std::end(position), std::begin(old_position));

						bool update_flag = false;
						update_flag |= ImGui::InputFloat3("Translation", &position[0], 2);
						need_update |= update_flag;

						if (update_flag)
						{
							translate(*p.mesh_, old_position, position);
						}
					}
				}
				else if (p.selecting_cell_ == EdgeSelect)
				{
					if (ImGui::BeginCombo("Sets", p.selected_edges_set_ ? p.selected_edges_set_->name().c_str()
																		: "-- select --"))
					{
						md->template foreach_cells_set<Edge>([&](CellsSet<MESH, Edge>& cs) {
							bool is_selected = &cs == p.selected_edges_set_;
							if (ImGui::Selectable(cs.name().c_str(), is_selected))
							{
								p.selected_edges_set_ = &cs;
								p.update_selected_edges_vbo();
								for (View* v : linked_views_)
									v->request_update();
							}
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						});
						ImGui::EndCombo();
					}
					if (p.selected_edges_set_)
					{
						ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
						if (ImGui::Button("X##selected_edges_set"))
							p.selected_edges_set_ = nullptr;
					}
					if (ImGui::Button("Create##edges_set"))
						md->template add_cells_set<Edge>();
					ImGui::TextUnformatted("Drawing parameters");
					need_update |=
						ImGui::ColorEdit3("color##edges", p.param_edge_->color_.data(), ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("width##edges", &(p.param_edge_->width_), 1.0f, 10.0f);

					if (p.selected_edges_set_ && !p.selected_edges_set_->empty())
					{
						ImGui::Separator();
						ImGui::TextUnformatted("Edge Tools");

						float position[3];
						set_position(p.selected_edges_position_, position);

						float old_position[3];
						std::copy(std::begin(position), std::end(position), std::begin(old_position));

						bool update_flag = false;
						update_flag |= ImGui::InputFloat3("Translation", &position[0], 2);
						need_update |= update_flag;

						if (update_flag)
						{
							translate(*p.mesh_, old_position, position);
						}
					}
				}
				else if (p.selecting_cell_ == FaceSelect)
				{
					if (ImGui::BeginCombo("Sets", p.selected_faces_set_ ? p.selected_faces_set_->name().c_str()
																		: "-- select --"))
					{
						md->template foreach_cells_set<Face>([&](CellsSet<MESH, Face>& cs) {
							bool is_selected = &cs == p.selected_faces_set_;
							if (ImGui::Selectable(cs.name().c_str(), is_selected))
							{
								p.selected_faces_set_ = &cs;
								p.update_selected_faces_vbo();
								for (View* v : linked_views_)
									v->request_update();
							}
							if (is_selected)
								ImGui::SetItemDefaultFocus();
						});
						ImGui::EndCombo();
					}
					if (p.selected_faces_set_)
					{
						ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
						if (ImGui::Button("X##selected_faces_set"))
							p.selected_faces_set_ = nullptr;
					}
					if (ImGui::Button("Create##faces_set"))
						md->template add_cells_set<Face>();
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("front color##flat", p.param_flat_->front_color_.data(),
													 ImGuiColorEditFlags_NoInputs);

					if (p.selected_faces_set_ && !p.selected_faces_set_->empty())
					{
						ImGui::Separator();
						ImGui::TextUnformatted("Faces Tools");

						float position[3];
						set_position(p.selected_faces_position_, position);

						float old_position[3];
						std::copy(std::begin(position), std::end(position), std::begin(old_position));

						bool update_flag = false;
						update_flag |= ImGui::InputFloat3("Translation", &position[0], 2);
						need_update |= update_flag;

						if (update_flag)
							translate(*p.mesh_, old_position, position);

						if (ImGui::Button("Extrude"))
							extrude(*p.mesh_);
					}
				}
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	const MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_TOOLS_H_

// Vertex v = ...
// value<Vec3>(*m,position,v) = Vec3(0,0,0)
// emit_attribute_changed(&m, vertex_attribute)

// cut edge
// collapse edge
// add face
// phi2_sew
