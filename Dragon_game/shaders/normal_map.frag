#version 410


in VS_OUT {

    vec3 normal;
    vec3 ftangent;
    vec3 fbinormal;
    vec3 fbinormal;
    vec3 fnormal

} vs_out;

out vec4 frag_color;

void main()
{
    mat3 TBN = (vs_out.ftangent, vs_out.fbinormal, vs_out.fnormal);


}