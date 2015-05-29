varying vec3 normal_eye;
varying vec3 position_eye;
varying vec3 color_out;

void main() {
  vec3 normal = normalize(normal_eye);
  vec3 view_vector = normalize(vec3(0, 0, 0) - position_eye);
  vec3 light_vector = normalize(vec3(0, 0, 0) - position_eye);
  vec3 halfv = normalize(light_vector + view_vector);
//  vec3 diffuse = max(0.0, dot(light_vector, normal)) * gl_LightSource[0].diffuse.xyz;
  vec3 diffuse = max(0.0, dot(light_vector, normal)) * color_out;// vec3(0.4, 0.5, 1.0);
//  vec3 ambient = (gl_FrontMaterial.ambient * gl_LightSource[0].ambient).xyz;
  vec3 ambient = vec3(0.0, 0.0, 0.0);
  vec3 specular = vec3(0, 0, 0);
  vec3 reflection = vec3(0, 0, 0);
  if (dot(normal, light_vector) > 0.0) {
    reflection = 2.0 * dot(normal, light_vector) * normal - light_vector;
    float d = dot(reflection, view_vector);
    if (d > 0.0) {
//      specular = pow(d, 64.0) * (gl_FrontMaterial.specular * gl_LightSource[0].specular).xyz;
      specular = pow(d, 128.0) * vec3(1, 1, 1);
    }
  }
  gl_FragColor = vec4(ambient + diffuse + 1.0 * specular, 1.0);
}
