#version 410

uniform samplerCube skybox_cube_map;
uniform vec3 light_position;
uniform float ellapsed_time_s;
uniform sampler2D normal_map;
uniform mat4 vertex_model_to_world;

in VS_OUT {
	vec3 vertex;
	vec3 normal;
    vec3 fV;
    vec3 normalCoord0;
    vec3 normalCoord1;
    vec3 normalCoord2;
    vec3 tangent;
    vec3 binormal;
    vec3 texcoord;
} fs_in;

out vec4 frag_color;

void main()
{
    vec3 N = normalize(fs_in.normal);
    vec3 V = normalize(fs_in.fV);

    // Map N to [-1, 1]
	vec3 normal_map0 = texture(normal_map, fs_in.normalCoord0.xy).xyz * 2 - 1;
    vec3 normal_map1 = texture(normal_map, fs_in.normalCoord1.xy).xyz * 2 - 1;
    vec3 normal_map2 = texture(normal_map, fs_in.normalCoord2.xy).xyz * 2 - 1;

    vec3 nbump = normalize(normal_map0 + normal_map1 + normal_map2);

    mat3 TBN = mat3(fs_in.tangent, fs_in.binormal, fs_in.normal); // obs normal är inte exakt som i beskrivningen här 


    mat4 WorldIT = inverse(transpose(vertex_model_to_world));
	mat4 TBN_extended = mat4(TBN);
	N = (WorldIT * TBN_extended * vec4(nbump, 0.0f)).xyz;

    // Water Color
    float facing = 1.0 - max(dot(V, N), 0);
    vec4 color_deep = vec4(0.545f, 0.0f, 0.0f, 1.0f);
    vec4 color_shallow = vec4(0.8f, 0.0f, 0.0f, 1.0f);

    // Fresnel
    float R0 = 0.02037;
    float fastFresnel = R0 + (1 - R0)* pow((1 - dot(V,N)), 5);

    frag_color = color_shallow;

    frag_color += (mix( color_deep, color_shallow, facing) - color_shallow );

    // Refraction
    vec4 refraction = vec4(refract(-V, N, 1.33), 0.0f);
    refraction = texture(skybox_cube_map, refraction.xyz);
    frag_color += refraction* (1 - fastFresnel);

    // Reflection
    vec3 R = reflect( -V, N);
    vec4 reflection = texture(skybox_cube_map, R);
    frag_color +=  reflection*fastFresnel;

}
