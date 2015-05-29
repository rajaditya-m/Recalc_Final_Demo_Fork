varying vec3 positionEye;
varying vec3 normal;
varying vec3 LightDir;
varying vec3 ViewDir;

//varying vec2 texCoord ;
//uniform float blend_coefficient;
//uniform samplerCube cubeMap;

uniform float test;
uniform sampler2D texture;
uniform sampler2D bumpMap;

void main() {
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //  vec4 view_vector = normalize(vec4(0, 0, 0, 1) - positionEye);
  //  vec4 light_vector = normalize(vec4(0, 0, 0, 1) - positionEye);
  //  vec4 halfv = normalize(light_vector + view_vector);
  //  //  vec3 diffuse = max(0.0, dot(light_vector, normalize(normalEye))) * (texture2D(texture, gl_TexCoord[0].xy) * gl_LightSource[0].diffuse).xyz;
  //  vec3 diffuse = max(0.0, dot(light_vector, normalize(normalEye))) * (texture2D(bumpMap, gl_TexCoord[0].xy) * gl_LightSource[0].diffuse).xyz;
  //  vec3 ambient = (gl_FrontMaterial.ambient * gl_LightSource[0].ambient).xyz;
  //  vec3 specular = vec3(0, 0, 0);
  //  vec4 reflection = vec4(0, 0, 0, 0);
  //  if (dot(normalize(normalEye), light_vector) > 0.0) {
  //    reflection = 2.0 * dot(normalize(normalEye), light_vector) * normalize(normalEye) - light_vector ;
  //    float d = dot(reflection, view_vector) ;
  //    if (d > 0.0) {
  //      specular = pow(d, 512.0) * (gl_FrontMaterial.specular * gl_LightSource[0].specular).xyz;
  //    }
  //  }
  //  gl_FragColor = vec4(test, test, test, 1);
  //  gl_FragColor = vec4(ambient, 1.0) + vec4(diffuse, 1.0) + 0.3 * vec4(specular, 1.0);
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // * blend_coefficient + textureColor * (1 - blend_coefficient) + vec4(specular,1.0);
  //  reflection = reflect(view_vector, normalEye) ;
  //  reflection = normalize(reflection);
  //  vec4 textureColor = textureCube(cubeMap, reflection.xyz);
  //  gl_FragColor = vec4(ambient,1.0) + vec4(diffuse,1.0) * blend_coefficient + textureColor * (1 - blend_coefficient) + vec4(specular,1.0);
  //  gl_FragColor = vec4(ambient,1.0) + vec4(diffuse,1.0) * blend_coefficient + textureColor * (1 - blend_coefficient);

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  vec3 normalEye = texture2D(bumpMap, gl_TexCoord[0].xy).xyz;
  normalEye = 2.0 * normalEye - 1.0;
  normalEye = normalize(normalEye);
  vec3 lightv = LightDir ;// modelViewMatrix*light.position
  vec3 viewv = ViewDir;
  vec3 halfv = normalize(lightv + viewv) ;
  vec4 textureColor = texture2D(texture, gl_TexCoord[0].xy);
//  vec4 textureColor = texture2D(texture, vec2(0.5, 0.5));
  //  float a0 = dot(normalize((-positionEye).xyz), normalize(normal.xyz));
  //  float a1 = dot(lightv, normalEye);
  //  vec4 diffuse = max(0.0, dot(normalize(-positionEye), normalize(normal))) * textureColor * gl_LightSource[0].diffuse; //
  //  vec4 diffuse = max(0.0, dot(lightv, normalEye)) * textureColor * gl_LightSource[0].diffuse; //
  vec4 diffuse = max(0.0, dot(lightv, normalEye)) * textureColor * gl_LightSource[0].diffuse; //
  vec4 ambient = gl_FrontMaterial.ambient * gl_LightSource[0].ambient;
  vec4 specular = vec4(0.0, 0.0, 0.0, 0.0);
  vec3 reflection = vec3(0.0, 0.0, 0.0); //= 2*dot(normal,lightv);//(2*dotProduct)*normal-L
  float cosine = dot(normalEye, lightv);
  if (cosine > 0.0) {
    //    float a = dot(normalEye, lightv);
    reflection = 2.0 * dot(normalEye, lightv) * normalEye - lightv ;
    //            reflection = reflect(lightv, normalEye);
    //            reflection.w = 0.0;
    //            reflection = normalize(reflection);
    float d = dot(reflection, viewv);
    //    float d = 2.0 * cosine * cosine - 1.0;
    //    float d = dot(reflection, normalEye);
    //    gl_FragColor = vec4(a, a, a, 1);
    //    reflection = 0.5 * reflection + vec4(0.5, 0.5, 0.5, 0);
    if (d > 0.0) {
      specular = pow(d, 80.0) * gl_FrontMaterial.specular * gl_LightSource[0].specular;
      //      gl_FragColor = specular;
      //        gl_FragColor = vec4(d, d, d, 1);
    } else {
      //      gl_FragColor = vec4(1, 0, 0, 1);
    }
  }
  vec4 color = ambient + diffuse + specular * 0.6;
  gl_FragColor = color;
//  gl_FragColor = textureColor;
//  gl_FragColor = vec4(0, 0, 0, 1);
//  gl_FragColor = vec4(0, 0, 0, color.x);

  //  gl_FragColor = diffuse;
  //  gl_FragColor = vec4(a0, a0, a0, 1);
  //  gl_FragColor = vec4(a1, a1, a1, 1);
  //  gl_FragColor = reflection;
  //      gl_FragColor = specular;
  //  lightv = normalize(lightv);
  //  lightv = 0.5 * lightv + vec4(.5, .5, .5, 0);
  //  gl_FragColor = vec4(h, h, h, 0);
  //  gl_FragColor = color * 1.5;
  //  gl_FragColor = normalEye;
  //  gl_FragColor = textureColor;
  //  reflection = reflect(viewv, normalEye) ;
  //  reflection = normalize(reflection) ;
  //  gl_FragColor = 0.4 * textureColor + 0.6 * color;
  //  gl_FragColor = texture2D(bumpMap, gl_TexCoord[0].xy);
}
