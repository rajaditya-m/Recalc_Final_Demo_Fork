attribute vec3 position;
attribute vec3 normal;
attribute vec3 color;

varying vec3 normal_eye;
varying vec3 position_eye;
varying vec3 color_out;
void main()
{ 
  normal_eye = gl_NormalMatrix * normal;
  position_eye = (gl_ModelViewMatrix * vec4(position, 1.0)).xyz;
  color_out = color;
  gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 1.0);

}
