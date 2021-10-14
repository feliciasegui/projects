#version 410


layout (location = 0) in vec3 vertex;
layout (location = 2) in vec3 texcoord;

uniform mat4 vertex_model_to_world;
uniform mat4 normal_model_to_world;
uniform mat4 vertex_world_to_clip;
uniform float ellapsed_time_s;
uniform vec3 camera_position;

out VS_OUT {
	vec3 vertex;
	vec3 normal;
    vec3 fV;
    vec3 normalCoord0;
    vec3 normalCoord1;
    vec3 normalCoord2;
    vec3 tangent;
    vec3 binormal;
    vec3 texcoord;
} vs_out;

float h_function(in vec2 position, in vec2 A, in mat2 D, in vec2 f, in vec2 p, in vec2 k, in float t, out float dHx, out float dHz){
    float alpha1 = sin((D[0][0] * position.x + D[0][1] * position.y) * f[0] + t * p[0]) * 0.5 + 0.5;
    float alpha2 = sin((D[1][0] * position.x + D[1][1] * position.y) * f[1] + t * p[1]) * 0.5 + 0.5;
    
    float dg1dx = 0.5 * k[0] * f[0] * A[0] * pow(alpha1, k[0]-1) * cos((D[0][0] * position.x + D[0][1] * position.y) * f[0] + t * p[0]) * D[0][0];
    float dg1dz = 0.5 * k[0] * f[0] * A[0] * pow(alpha1, k[0]-1) * cos((D[0][0] * position.x + D[0][1] * position.y) * f[0] + t * p[0]) * D[0][1];
    float dg2dx = 0.5 * k[1] * f[1] * A[1] * pow(alpha2, k[1]-1) * cos((D[1][0] * position.x + D[1][1] * position.y) * f[1] + t * p[1]) * D[1][0];
    float dg2dz = 0.5 * k[1] * f[1] * A[1] * pow(alpha2, k[1]-1) * cos((D[1][0] * position.x + D[1][1] * position.y) * f[1] + t * p[1]) * D[1][1];
    
    dHx = dg1dx + dg2dx;
    dHz = dg1dz + dg2dz;

    return A[0] * pow(alpha1, k[0]) + A[1] * pow(alpha2, k[1]);
}



void main()
{
    vec2 A = vec2(1.0f, 0.5f);
    mat2 D =  mat2(vec2( -1.0f, 0.0f), vec2( -0.7f, 0.7f ));
    vec2 f = vec2(0.2f, 0.4f);
    vec2 p = vec2(0.5f, 1.3f);
    vec2 k = vec2(2.0f, 2.0f);
    float dHx;
    float dHz;


    // Displacing vertices to create waves
    vec3 displaced_vertex = vertex; 
    displaced_vertex.y += h_function(vertex.xz, A, D, f, p, k, ellapsed_time_s, dHx, dHz);

    vec3 n = vec3(-dHx, 1.0f, -dHz);

    // Calculate V
	vec3 worldPos = (vertex_model_to_world*vec4(displaced_vertex,1)).xyz;  
	vs_out.fV = camera_position  - worldPos; 

    // From diffuse
	vs_out.vertex = vec3(vertex_model_to_world * vec4(displaced_vertex, 1.0));
	vs_out.normal = vec3(normal_model_to_world * vec4(n, 0.0));

    // Animated Normal Mapping
    vec2 texScale  = vec2(8, 4);
    float normalTime = mod(ellapsed_time_s, 100.0);
    vec2 normalSpeed = vec2(-0.05, 0);
    
    vs_out.normalCoord0.xy = texcoord.xz*texScale   + normalTime*normalSpeed;
    vs_out.normalCoord1.xy = texcoord.xz*texScale*2 + normalTime*normalSpeed*4;
    vs_out.normalCoord2.xy = texcoord.xz*texScale*4 + normalTime*normalSpeed*8;

    vs_out.tangent = vec3(1.0f, dHx, 0.0f);
    vs_out.binormal = vec3(0.0f, dHz, 1.0f);
	gl_Position = vertex_world_to_clip * vertex_model_to_world * vec4(displaced_vertex, 1.0);
}



