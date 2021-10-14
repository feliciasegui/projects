#version 410

uniform float illum; // shininess?
uniform vec3 Ka; // ambient
uniform vec3 Kd; // Diffuse
uniform vec3 Ks; // Specular


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

	// NOrmalizing normal, light vector  and view vector
	vec3 N = normalize(fs_in.fN);
	vec3 L = normalize(fs_in.fL);
	vec3 V = normalize(fs_in.fV);
	vec3 R = reflect( -L, N);

	// Calculating terms (changed from just constant diffuse color to texture)
	vec3 diffuse = Kd * max(dot(N,L),0.0);
	vec3 specular = Ks * pow(max(dot(R,V),0.0 ), illum);

	// Summing color terms
	frag_color.xyz = Ka  + diffuse + specular;
	frag_color.w = 1.0f;

}