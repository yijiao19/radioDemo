#version 130

precision highp float;
uniform sampler2D tex0;
uniform sampler2DShadow tex1;
uniform vec3 lightPosition;

in vec4 outColor;
in vec2 texCoord;
in vec3 normal; 
in vec3 eyeVector; 
in vec3 modelViewPosition; 
in vec4 shadowMapCoord; 

uniform float fu;


out vec4 fragmentColor;


void main() 
{
	//vec3 posToLight = normalize(lightPosition - modelViewPosition);

	vec3 posToLight= vec3(1.0,2.0,1.0);
	normalize(posToLight);

	float intensity = max(0, dot(posToLight, normalize(normal)));
	
    intensity *= textureProj(tex1, shadowMapCoord); 

	intensity+= 0.01;

	intensity = intensity * fu;

	//fragmentColor.rgb = normalize(normal);
	fragmentColor = texture2D(tex0, vec2(min(0.99,intensity),0.5));
   // fragmentColor = texture2D(tex0, vec2(intensity,0.5));

   //fragmentColor.rgb = vec3(intensity,intensity,intensity);
	fragmentColor.a = 1.0;
}