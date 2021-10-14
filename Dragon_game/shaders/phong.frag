#version 410

uniform float shininess;
uniform vec3 ambient;
uniform sampler2D diffuse_text;
uniform sampler2D specular_map;
uniform sampler2D normal_map;
uniform int has_diffuse_texture;
uniform mat4 vertex_model_to_world;
uniform int use_normal_mapping;

in VS_OUT {
	vec3 fN;
	vec3 fL;
	vec3 fV;
	vec2 texcoord;
	vec3 ftangent;
	vec3 fbinormal;
} fs_in;

out vec4 frag_color;

void main()
{

	// Creating TBN Matrix
	mat3 TBN = mat3(fs_in.ftangent, fs_in.fbinormal, fs_in.fN);

	// NOrmalizing normal, light vector  and view vector
	vec3 N = normalize(fs_in.fN);
	vec3 L = normalize(fs_in.fL);
	vec3 V = normalize(fs_in.fV);
	vec3 R = reflect( -L, N);

	// Map N to [-1, 1]
	vec3 normal_map = texture(normal_map, fs_in.texcoord).xyz;
	normal_map = normal_map * 2 - 1;


	if (use_normal_mapping != 0){
		mat4 WorldIT = inverse(transpose(vertex_model_to_world));
		mat4 TBN_extended = mat4(TBN);
		N = (WorldIT * TBN_extended * vec4(normal_map, 0.0f)).xyz;
	}
	

	// Calculating terms (changed from just constant diffuse color to texture)
	vec3 diffuse = texture(diffuse_text, fs_in.texcoord).xyz * max(dot(N,L),0.0);
	vec3 specular = texture(specular_map, fs_in.texcoord).xyz * pow(max(dot(R,V),0.0 ), shininess);

	// Summing color terms
	frag_color.xyz = ambient + diffuse + specular ;
	frag_color.w = 1.0f;

}