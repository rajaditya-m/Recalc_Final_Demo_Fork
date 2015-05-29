varying vec3 positionEye;
varying vec3 normalEye;
varying vec2 texCoord;
uniform samplerCube cubeMap;
uniform int has_texture;

void main() {
  vec3 normal = normalize(normalEye);
  vec3 view_vector = normalize(vec3(0, 0, 0) - positionEye);
  vec3 light_vector = normalize(vec3(0, 0, 0) - positionEye);
  vec3 halfv = normalize(light_vector + view_vector);
  vec3 diffuse = max(0.0, dot(light_vector, normal)) * gl_FrontMaterial.diffuse.xyz;
  //  vec3 diffuse = max(0, dot(light_vector, normal)) * (gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse).xyz;
  //  vec3 diffuse = max(0, dot(light_vector, normal)) * vec3(224.0 / 255.0, 223.0 / 255.0, 219.0 / 255.0);
  //  vec3 diffuse = max(0.0, dot(light_vector, normal)) * vec3(0.50754, 0.50754,	0.50754);
  //  vec3 diffuse = max(0.0, dot(light_vector, normal)) * vec3(0.75164, 0.60648, 0.22648);
  vec3 ambient = gl_FrontMaterial.ambient.xyz;
  //  vec3 ambient = (gl_FrontMaterial.ambient * gl_LightSource[0].ambient).xyz;
  //  vec3 ambient = vec3(0.19225, 0.19225, 0.19225);//(gl_FrontMaterial.ambient * gl_LightSource[0].ambient).xyz;
  //  vec3 ambient = vec3(0.24725, 0.1995, 0.0745);//(gl_FrontMaterial.ambient * gl_LightSource[0].ambient).xyz;
  vec3 specular = vec3(0, 0, 0);
  vec3 reflection = vec3(0, 0, 0);
  if (dot(normal, light_vector) > 0.0) {
    reflection = 2.0 * dot(normal, light_vector) * normal - light_vector ;
    float d = dot(reflection, view_vector);
    if (d > 0.0) {
      specular = pow(d, gl_FrontMaterial.shininess) * gl_FrontMaterial.specular.xyz;
      //      float shininess = 128.0 * 0.4;
      //      specular = pow(d, 512) * (gl_FrontMaterial.specular * gl_LightSource[0].specular).xyz;
      //      specular = pow(d, shininess) * vec3(0.508273, 0.508273, 0.508273);
      //      specular = pow(d, shininess) * vec3(0.628281, 0.555802, 0.366065);
      //      specular = pow(d, 64.0) * vec3(1, 1, 1);;
    }
  }
  reflection = reflect(view_vector, normalEye) ;
  reflection = normalize(reflection);
  vec3 textureColor = textureCube(cubeMap, reflection.xyz).xyz;
  float blend_coefficient = 1.0;
  gl_FragColor = vec4(ambient + diffuse * blend_coefficient + textureColor * (1.0 - blend_coefficient) + specular, 1.0);
  //  gl_FragColor = vec4(textureColor, 1.0);
}
