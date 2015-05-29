attribute vec3 position;
attribute vec3 normal;
attribute vec2 tex;

varying vec3 positionEye;
varying vec3 normalEye;
varying vec2 texCoord;

void main() {
  texCoord = tex;
  positionEye = (gl_ModelViewMatrix * vec4(position, 1.0)).xyz;
  normalEye = +1.0 * gl_NormalMatrix * normal;
  gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 1.0);
}
