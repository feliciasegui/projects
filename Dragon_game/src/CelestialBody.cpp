#include "CelestialBody.hpp"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/trigonometric.hpp>

#include "core/helpers.hpp"
#include "core/Log.h"

CelestialBody::CelestialBody(bonobo::mesh_data const& shape,
                             GLuint const* program,
                             GLuint diffuse_texture_id)
{
	_body.node.set_geometry(shape);
	_body.node.add_texture("diffuse_texture", diffuse_texture_id, GL_TEXTURE_2D);
	_body.node.set_program(program);
}

glm::mat4 CelestialBody::render(std::chrono::microseconds elapsed_time,
                                glm::mat4 const& view_projection,
                                glm::mat4 const& parent_transform,
                                bool show_basis)
{
	// Convert the duration from microseconds to seconds.
	auto const elapsed_time_s = std::chrono::duration<float>(elapsed_time).count();
	// If a different ratio was needed, for example a duration in
	// milliseconds, the following would have been used:
	// auto const elapsed_time_ms = std::chrono::duration<float, std::milli>(elapsed_time).count();

	//_body.spin.rotation_angle = -glm::half_pi<float>() / 2.0f;

	//glm::mat4 world = parent_transform;

	// Overwriting world matrix task 3.1.1, 3.1.2 and task 3.2.2
	//world = glm::scale(glm::mat4(1.0f), _body.scale);

	// Overwriting world matrix with first rotation matrix (around y axis) task 3.2.1 and 3.2.2
	_body.spin.rotation_angle = _body.spin.rotation_angle  + _body.spin.speed * elapsed_time_s;
	_body.orbit.rotation_angle = _body.orbit.rotation_angle + _body.orbit.speed * elapsed_time_s;

	// Rotation and translation matrices
	glm::mat4 scale = glm::scale(glm::mat4(1.0f), _body.scale);
	glm::mat4 r1 = glm::rotate(glm::mat4(1.0f),_body.spin.rotation_angle, glm::vec3(0, 1, 0)); // Task 3.2.1
	glm::mat4 r2 = glm::rotate(glm::mat4(1.0f),-_body.spin.axial_tilt, glm::vec3(0, 0, 1)); // Task 3.2.3
	glm::mat4 t0 = glm::translate(glm::mat4(1.0f), _body.orbit.radius*glm::vec3(1, 0, 0)); // Task 3.3.1
	glm::mat4 r10 = glm::rotate(glm::mat4(1.0f),_body.orbit.rotation_angle, glm::vec3(0, 1, 0));
	glm::mat4 r20 = glm::rotate(glm::mat4(1.0f),-_body.orbit.inclination, glm::vec3(0, 0, 1));

	// Task 3.2.3 and 3.2.6
	glm::mat4 r21 = r2*r1*scale;

	// Task 3.3.5
	glm::mat4 r210 = r20*r10;

	// Combining all matrices 
	glm::mat4 world = r210*t0*r21; 

	// Include parent transform task 4.1.1;
	world = parent_transform * world;

	glm::mat4 child_transform = parent_transform*r210*t0;

	//
	// Rings
	//
	//glm::mat4 scale_r = glm::scale(glm::mat4(1.0f), _ring.scale); This does not work, why?
	glm::mat4 r_rings = glm::rotate(glm::mat4(1.0f), glm::half_pi<float>() , glm::vec3(1, 0, 0));
	glm::mat4 world_r = child_transform*r_rings;//*scale_r;

	if (show_basis)
	{
		bonobo::renderBasis(1.0f, 2.0f, view_projection, world);
	}

	// Note: The second argument of `node::render()` is supposed to be the
	// parent transform of the node, not the whole world matrix, as the
	// node internally manages its local transforms. However in our case we
	// manage all the local transforms ourselves, so the internal transform
	// of the node is just the identity matrix and we can forward the whole
	// world matrix.
	_body.node.render(view_projection, world);
	_ring.node.render(view_projection, world_r);




	return child_transform;
}

void CelestialBody::add_child(CelestialBody* child)
{
	_children.push_back(child);
}

std::vector<CelestialBody*> const& CelestialBody::get_children() const
{
	return _children;
}

void CelestialBody::set_orbit(OrbitConfiguration const& configuration)
{
	_body.orbit.radius = configuration.radius;
	_body.orbit.inclination = configuration.inclination;
	_body.orbit.speed = configuration.speed;
	_body.orbit.rotation_angle = 0.0f;
}

void CelestialBody::set_scale(glm::vec3 const& scale)
{
	_body.scale = scale;
}

void CelestialBody::set_spin(SpinConfiguration const& configuration)
{
	_body.spin.axial_tilt = configuration.axial_tilt;
	_body.spin.speed = configuration.speed;
	_body.spin.rotation_angle = 0.0f;
}

void CelestialBody::set_ring(bonobo::mesh_data const& shape,
                             GLuint const* program,
                             GLuint diffuse_texture_id,
                             glm::vec2 const& scale)
{
	_ring.node.set_geometry(shape);
	_ring.node.add_texture("diffuse_texture", diffuse_texture_id, GL_TEXTURE_2D);
	_ring.node.set_program(program);

	_ring.scale = scale;
}
