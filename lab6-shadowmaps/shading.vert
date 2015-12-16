#version 130
//in vec3		vertex;
//in vec3		color;
out vec4	outColor;
//in	vec2	texCoordIn;	
out	vec2	texCoord;	
//in  vec3	normalIn;		
out vec3	normal;			
out vec3	modelViewPosition;		
uniform mat4 normalMatrix;	 
uniform mat4 modelViewMatrix;
uniform mat4 modelViewProjectionMatrix; 
uniform mat4 lightMatrix; 
out vec4 shadowMapCoord;


void main() 
{
	vec3 vertex = gl_Vertex.xyz;
	vec3 color = gl_Color.xyz;
	gl_Position = modelViewProjectionMatrix * vec4(vertex,1);
	outColor = vec4(color,1); 
	//texCoord = texCoordIn; 
	texCoord = gl_MultiTexCoord0.xy; 
	//normal = (normalMatrix * vec4(normalIn,0.0)).xyz;
	normal = (mat3(normalMatrix) *gl_Normal).xyz;
	modelViewPosition = (modelViewMatrix * vec4(vertex, 1)).xyz;
	shadowMapCoord = lightMatrix * vec4(modelViewPosition, 1.0);
	shadowMapCoord.xyz *= vec3(0.5f, 0.5f, 0.5f);	
	shadowMapCoord.xyz += shadowMapCoord.w * vec3(0.5f, 0.5f, 0.5f);

}