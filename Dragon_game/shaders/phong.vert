#version 410

layout (location = 0) in vec3 vertex;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec3 texcoord;
layout (location = 3) in vec3 tangent;
layout (location = 4) in vec3 binormal;
uniform vec3 light_position;
uniform vec3 camera_position;
uniform mat4 vertex_model_to_world;
uniform mat4 vertex_world_to_clip;


out VS_OUT {

	vec3 fN;
	vec3 fL;
	vec3 fV;
	vec2 texcoord;
	vec3 ftangent;
	vec3 fbinormal;

} vs_out;

void main()
{
	vs_out.ftangent = tangent;
	vs_out.fbinormal = binormal;

	// Vertices defined in model space, transformed to world space with matrix
	vec3 worldPos = (vertex_model_to_world*vec4(vertex,1)).xyz;  

	// ??
	mat4 WorldIT = inverse(transpose(vertex_model_to_world));
	vs_out.fN = (WorldIT*vec4(normal,0)).xyz;  

	// Camera position is already in world space, calculating the vector from the object to the camera
	vs_out.fV = camera_position  - worldPos; 

	// Light position already in world space, calculating the vector from the object to the light
	vs_out.fL = light_position - worldPos;     

	// Creating the model to clip matrix, which is needed for the rasterizer
	mat4 model_to_clip = vertex_world_to_clip * vertex_model_to_world;
	gl_Position =  model_to_clip * vec4(vertex, 1.0);

	// Texture coordinates
	vs_out.texcoord = vec2(texcoord.x, texcoord.y);

}
