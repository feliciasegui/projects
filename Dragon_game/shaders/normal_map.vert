#version 410

layout (location = 1) in vec3 normal;
layout (location = 3) in vec3 tangent;
layout (location = 4) in vec3 binormal;

out VS_OUT {

    vec3 ftangent;
    vec3 fbinormal;
    vec3 fnormal

} vs_out;

void main()
{
    // Tangent space -> Model space 
	

}
