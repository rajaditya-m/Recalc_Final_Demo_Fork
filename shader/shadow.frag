varying vec4 fragment_position;

uniform sampler2D shadow_texture;
uniform sampler2D texture;
uniform mat4 		biased_MVP;
uniform	vec3		light_position;
uniform float ap1;
uniform float ap2;

float rand1(vec2 co) {
  return fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5453);
}


float rand2(vec2 co) {
  return fract(sin(dot(co.xy, vec2(12.9898, 68.233))) * 33058.5453);
}

void main() {
  //    vec3 V = vec3(fragment_position);                                  //in the world space
  //    V = normalize(vec3(gl_ModelViewMatrixInverse * vec4(V, 0.0)));      //from eye to world
  //   vec4 texture_color=Environment(V);
  vec3 L = normalize(light_position - fragment_position.xyz);
  vec3 N = vec3(0, 1, 0);

  vec3 shadow_coord = (biased_MVP * fragment_position).xyz;

  float visibility = 0.0;

  for (int i = -3; i <=3; i++)
    for (int j = -3; j <= 3; j++) {
      vec2 pos = shadow_coord.xy + (vec2(float(i), j) + vec2(rand1(shadow_coord.xy), rand2(shadow_coord.xy)) * ap2) * ap1;

      if (texture2D(shadow_texture, pos ).z  <  shadow_coord.z)
        visibility = visibility + 0.8;
      else visibility = visibility + 1.0;

    }
    visibility=visibility/49.0;// 25.0;
  //visibility = visibility / 9.0;
  float color = dot(L, N) * pow(visibility, 2.0) * 0.8 + 0.1;
  //  float color=dot(L, N)*(visibility) *0.8 + 0.1;
  //  vec4 textureColor = texture2D(texture, vec2(0.5,0.5));
  //  vec4 textureColor = color * texture2D(texture, gl_TexCoord[0].xy);
  vec4 textureColor = color * texture2D(texture, gl_TexCoord[0].xy);
  //  vec4 textureColor = texture2D(texture, gl_TexCoord[0].xy);
  gl_FragColor = textureColor;
  //  gl_FragColor = vec4(color, color, color*0.93, 1);
  //   	gl_FragColor = vec4(texture2D(shadow_texture, vec2(0.5, 0.5) ).z, 0, 0, 1);
  //    	gl_FragColor = vec4(texture2D(env_texture, shadow_coord).xyz, 1);
  //  gl_FragColor = vec4(shadow_coord, 1);
  //  gl_FragColor = vec4(0, 0, 0, 0);

}
