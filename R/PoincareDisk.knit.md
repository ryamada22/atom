---
title: "ポアンカレ　ディスク"
output: html_document
---




## 双曲面模型

この節は、[Hyperboloid model(Wikipedia記事)](https://en.wikipedia.org/wiki/Hyperboloid_model) 
に従って書いてある。

双曲面模型とは、n+1次元空間にある、ある(超)双曲面上の点を、n次元空間中の単位球体(円板)にあるルールで射影するもの。

n+1次元空間は、ミンコフスキー空間になっており、内積の定義に負の固有値が入っているため、直線・平行線・距離の概念がユークリッド空間のそれとは異なっている。

ミンコフスキー空間のm+1次元(超)平面は、n次元平面に射影されてm次元(超)平面になるが、それは、n+1次元空間の超平面が、双曲面とよぎった点の写像としてあらわされる。

n=2の場合は、ミンコフスキー空間は３次元、射影は２次元円板、ミンコフスキー空間の平面が、ミンコフスキー空間の双曲面とよぎって、３次元空間上の曲線となり、その写像が、２次元円板上の曲線となる。


$$
x^2 + y^2 - z^2 \le 0; z \ge 0
$$
の領域を考える。

この空間に
$$
x^2 + y^2 - z^2 = -1
$$
という双曲面がある。


### 3次元ミンコフスキー空間座標の単位円板への射影

その領域の点(x,y,z)を、$z=0$平面上の半径１の円板$x^2 + y^2 \le 1$上に射影することとする。

点(0,0,-1)と(x,y,z)とを通る直線が、$z=0$をよぎる点の座標(x',y',0)は
以下の式で与えられる。

$$
x' = \frac{x}{z+1}\\
y' = \frac{y}{z+1}
$$

逆写像は、次の式で与えられる。

$$
x = \frac{2x'}{1-x'^2-y'^2}\\
y = \frac{2y'}{1-x'^2-y'^2}\\
z = \frac{1+x'^2+y'^2}{1-x'^2-y'^2}
$$

```r
my.pdisk.coords <- function(x){
  x. <- x[1]/(x[3]+1)
  y. <- x[2]/(x[3]+1)
  return(c(x.,y.))
}
# 円板上の点の座標をx^2 + y^2 - z^2 = -1 という双曲面上に逆写像する
my.pdisk.coords.inv <- function(x){
  z <- (1+sum(x^2))/(1-sum(x^2))
  tmp <- 2*x/(1-sum(x^2))
  return(c(tmp,z))
}
#my.pdisk.coords.inv_bk <- function(x){
#  tmp <- x/sqrt(1-sum(x^2))
#  z <- sqrt(sum(tmp^2)+1)
#  return(c(tmp,z))
#}
# x^2+y^2-z^2=0という円錐に逆写像
#my.pdisk.coords.inv.cone <- function(x){
#  k <- (-1+sqrt(sum(x^2)))/(sum(x^2)-1)
#  
#  return(k * c(x,1))
#}
```

### 双曲面上の点の射影


```r
n <- 10000
x <- runif(n,min=-10,max=10)
y <- runif(n,min=-10,max=10)
z <- sqrt(x^2 + y^2 + 1)

xyz <- cbind(x,y,z)
```


```r
plot3d(xyz)
```

<script>/*
* Copyright (C) 2009 Apple Inc. All Rights Reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY APPLE INC. ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
* PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL APPLE INC. OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
* OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* Copyright (2016) Duncan Murdoch - fixed CanvasMatrix4.ortho,
* cleaned up.
*/
/*
CanvasMatrix4 class
This class implements a 4x4 matrix. It has functions which
duplicate the functionality of the OpenGL matrix stack and
glut functions.
IDL:
[
Constructor(in CanvasMatrix4 matrix),           // copy passed matrix into new CanvasMatrix4
Constructor(in sequence<float> array)           // create new CanvasMatrix4 with 16 floats (row major)
Constructor()                                   // create new CanvasMatrix4 with identity matrix
]
interface CanvasMatrix4 {
attribute float m11;
attribute float m12;
attribute float m13;
attribute float m14;
attribute float m21;
attribute float m22;
attribute float m23;
attribute float m24;
attribute float m31;
attribute float m32;
attribute float m33;
attribute float m34;
attribute float m41;
attribute float m42;
attribute float m43;
attribute float m44;
void load(in CanvasMatrix4 matrix);                 // copy the values from the passed matrix
void load(in sequence<float> array);                // copy 16 floats into the matrix
sequence<float> getAsArray();                       // return the matrix as an array of 16 floats
WebGLFloatArray getAsCanvasFloatArray();           // return the matrix as a WebGLFloatArray with 16 values
void makeIdentity();                                // replace the matrix with identity
void transpose();                                   // replace the matrix with its transpose
void invert();                                      // replace the matrix with its inverse
void translate(in float x, in float y, in float z); // multiply the matrix by passed translation values on the right
void scale(in float x, in float y, in float z);     // multiply the matrix by passed scale values on the right
void rotate(in float angle,                         // multiply the matrix by passed rotation values on the right
in float x, in float y, in float z);    // (angle is in degrees)
void multRight(in CanvasMatrix matrix);             // multiply the matrix by the passed matrix on the right
void multLeft(in CanvasMatrix matrix);              // multiply the matrix by the passed matrix on the left
void ortho(in float left, in float right,           // multiply the matrix by the passed ortho values on the right
in float bottom, in float top,
in float near, in float far);
void frustum(in float left, in float right,         // multiply the matrix by the passed frustum values on the right
in float bottom, in float top,
in float near, in float far);
void perspective(in float fovy, in float aspect,    // multiply the matrix by the passed perspective values on the right
in float zNear, in float zFar);
void lookat(in float eyex, in float eyey, in float eyez,    // multiply the matrix by the passed lookat
in float ctrx, in float ctry, in float ctrz,    // values on the right
in float upx, in float upy, in float upz);
}
*/
CanvasMatrix4 = function(m)
{
if (typeof m == 'object') {
if ("length" in m && m.length >= 16) {
this.load(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
return;
}
else if (m instanceof CanvasMatrix4) {
this.load(m);
return;
}
}
this.makeIdentity();
};
CanvasMatrix4.prototype.load = function()
{
if (arguments.length == 1 && typeof arguments[0] == 'object') {
var matrix = arguments[0];
if ("length" in matrix && matrix.length == 16) {
this.m11 = matrix[0];
this.m12 = matrix[1];
this.m13 = matrix[2];
this.m14 = matrix[3];
this.m21 = matrix[4];
this.m22 = matrix[5];
this.m23 = matrix[6];
this.m24 = matrix[7];
this.m31 = matrix[8];
this.m32 = matrix[9];
this.m33 = matrix[10];
this.m34 = matrix[11];
this.m41 = matrix[12];
this.m42 = matrix[13];
this.m43 = matrix[14];
this.m44 = matrix[15];
return;
}
if (arguments[0] instanceof CanvasMatrix4) {
this.m11 = matrix.m11;
this.m12 = matrix.m12;
this.m13 = matrix.m13;
this.m14 = matrix.m14;
this.m21 = matrix.m21;
this.m22 = matrix.m22;
this.m23 = matrix.m23;
this.m24 = matrix.m24;
this.m31 = matrix.m31;
this.m32 = matrix.m32;
this.m33 = matrix.m33;
this.m34 = matrix.m34;
this.m41 = matrix.m41;
this.m42 = matrix.m42;
this.m43 = matrix.m43;
this.m44 = matrix.m44;
return;
}
}
this.makeIdentity();
};
CanvasMatrix4.prototype.getAsArray = function()
{
return [
this.m11, this.m12, this.m13, this.m14,
this.m21, this.m22, this.m23, this.m24,
this.m31, this.m32, this.m33, this.m34,
this.m41, this.m42, this.m43, this.m44
];
};
CanvasMatrix4.prototype.getAsWebGLFloatArray = function()
{
return new WebGLFloatArray(this.getAsArray());
};
CanvasMatrix4.prototype.makeIdentity = function()
{
this.m11 = 1;
this.m12 = 0;
this.m13 = 0;
this.m14 = 0;
this.m21 = 0;
this.m22 = 1;
this.m23 = 0;
this.m24 = 0;
this.m31 = 0;
this.m32 = 0;
this.m33 = 1;
this.m34 = 0;
this.m41 = 0;
this.m42 = 0;
this.m43 = 0;
this.m44 = 1;
};
CanvasMatrix4.prototype.transpose = function()
{
var tmp = this.m12;
this.m12 = this.m21;
this.m21 = tmp;
tmp = this.m13;
this.m13 = this.m31;
this.m31 = tmp;
tmp = this.m14;
this.m14 = this.m41;
this.m41 = tmp;
tmp = this.m23;
this.m23 = this.m32;
this.m32 = tmp;
tmp = this.m24;
this.m24 = this.m42;
this.m42 = tmp;
tmp = this.m34;
this.m34 = this.m43;
this.m43 = tmp;
};
CanvasMatrix4.prototype.invert = function()
{
// Calculate the 4x4 determinant
// If the determinant is zero,
// then the inverse matrix is not unique.
var det = this._determinant4x4();
if (Math.abs(det) < 1e-8)
return null;
this._makeAdjoint();
// Scale the adjoint matrix to get the inverse
this.m11 /= det;
this.m12 /= det;
this.m13 /= det;
this.m14 /= det;
this.m21 /= det;
this.m22 /= det;
this.m23 /= det;
this.m24 /= det;
this.m31 /= det;
this.m32 /= det;
this.m33 /= det;
this.m34 /= det;
this.m41 /= det;
this.m42 /= det;
this.m43 /= det;
this.m44 /= det;
};
CanvasMatrix4.prototype.translate = function(x,y,z)
{
if (x === undefined)
x = 0;
if (y === undefined)
y = 0;
if (z === undefined)
z = 0;
var matrix = new CanvasMatrix4();
matrix.m41 = x;
matrix.m42 = y;
matrix.m43 = z;
this.multRight(matrix);
};
CanvasMatrix4.prototype.scale = function(x,y,z)
{
if (x === undefined)
x = 1;
if (z === undefined) {
if (y === undefined) {
y = x;
z = x;
}
else
z = 1;
}
else if (y === undefined)
y = x;
var matrix = new CanvasMatrix4();
matrix.m11 = x;
matrix.m22 = y;
matrix.m33 = z;
this.multRight(matrix);
};
CanvasMatrix4.prototype.rotate = function(angle,x,y,z)
{
// angles are in degrees. Switch to radians
angle = angle / 180 * Math.PI;
angle /= 2;
var sinA = Math.sin(angle);
var cosA = Math.cos(angle);
var sinA2 = sinA * sinA;
// normalize
var length = Math.sqrt(x * x + y * y + z * z);
if (length === 0) {
// bad vector, just use something reasonable
x = 0;
y = 0;
z = 1;
} else if (length != 1) {
x /= length;
y /= length;
z /= length;
}
var mat = new CanvasMatrix4();
// optimize case where axis is along major axis
if (x == 1 && y === 0 && z === 0) {
mat.m11 = 1;
mat.m12 = 0;
mat.m13 = 0;
mat.m21 = 0;
mat.m22 = 1 - 2 * sinA2;
mat.m23 = 2 * sinA * cosA;
mat.m31 = 0;
mat.m32 = -2 * sinA * cosA;
mat.m33 = 1 - 2 * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else if (x === 0 && y == 1 && z === 0) {
mat.m11 = 1 - 2 * sinA2;
mat.m12 = 0;
mat.m13 = -2 * sinA * cosA;
mat.m21 = 0;
mat.m22 = 1;
mat.m23 = 0;
mat.m31 = 2 * sinA * cosA;
mat.m32 = 0;
mat.m33 = 1 - 2 * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else if (x === 0 && y === 0 && z == 1) {
mat.m11 = 1 - 2 * sinA2;
mat.m12 = 2 * sinA * cosA;
mat.m13 = 0;
mat.m21 = -2 * sinA * cosA;
mat.m22 = 1 - 2 * sinA2;
mat.m23 = 0;
mat.m31 = 0;
mat.m32 = 0;
mat.m33 = 1;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else {
var x2 = x*x;
var y2 = y*y;
var z2 = z*z;
mat.m11 = 1 - 2 * (y2 + z2) * sinA2;
mat.m12 = 2 * (x * y * sinA2 + z * sinA * cosA);
mat.m13 = 2 * (x * z * sinA2 - y * sinA * cosA);
mat.m21 = 2 * (y * x * sinA2 - z * sinA * cosA);
mat.m22 = 1 - 2 * (z2 + x2) * sinA2;
mat.m23 = 2 * (y * z * sinA2 + x * sinA * cosA);
mat.m31 = 2 * (z * x * sinA2 + y * sinA * cosA);
mat.m32 = 2 * (z * y * sinA2 - x * sinA * cosA);
mat.m33 = 1 - 2 * (x2 + y2) * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
}
this.multRight(mat);
};
CanvasMatrix4.prototype.multRight = function(mat)
{
var m11 = (this.m11 * mat.m11 + this.m12 * mat.m21 +
this.m13 * mat.m31 + this.m14 * mat.m41);
var m12 = (this.m11 * mat.m12 + this.m12 * mat.m22 +
this.m13 * mat.m32 + this.m14 * mat.m42);
var m13 = (this.m11 * mat.m13 + this.m12 * mat.m23 +
this.m13 * mat.m33 + this.m14 * mat.m43);
var m14 = (this.m11 * mat.m14 + this.m12 * mat.m24 +
this.m13 * mat.m34 + this.m14 * mat.m44);
var m21 = (this.m21 * mat.m11 + this.m22 * mat.m21 +
this.m23 * mat.m31 + this.m24 * mat.m41);
var m22 = (this.m21 * mat.m12 + this.m22 * mat.m22 +
this.m23 * mat.m32 + this.m24 * mat.m42);
var m23 = (this.m21 * mat.m13 + this.m22 * mat.m23 +
this.m23 * mat.m33 + this.m24 * mat.m43);
var m24 = (this.m21 * mat.m14 + this.m22 * mat.m24 +
this.m23 * mat.m34 + this.m24 * mat.m44);
var m31 = (this.m31 * mat.m11 + this.m32 * mat.m21 +
this.m33 * mat.m31 + this.m34 * mat.m41);
var m32 = (this.m31 * mat.m12 + this.m32 * mat.m22 +
this.m33 * mat.m32 + this.m34 * mat.m42);
var m33 = (this.m31 * mat.m13 + this.m32 * mat.m23 +
this.m33 * mat.m33 + this.m34 * mat.m43);
var m34 = (this.m31 * mat.m14 + this.m32 * mat.m24 +
this.m33 * mat.m34 + this.m34 * mat.m44);
var m41 = (this.m41 * mat.m11 + this.m42 * mat.m21 +
this.m43 * mat.m31 + this.m44 * mat.m41);
var m42 = (this.m41 * mat.m12 + this.m42 * mat.m22 +
this.m43 * mat.m32 + this.m44 * mat.m42);
var m43 = (this.m41 * mat.m13 + this.m42 * mat.m23 +
this.m43 * mat.m33 + this.m44 * mat.m43);
var m44 = (this.m41 * mat.m14 + this.m42 * mat.m24 +
this.m43 * mat.m34 + this.m44 * mat.m44);
this.m11 = m11;
this.m12 = m12;
this.m13 = m13;
this.m14 = m14;
this.m21 = m21;
this.m22 = m22;
this.m23 = m23;
this.m24 = m24;
this.m31 = m31;
this.m32 = m32;
this.m33 = m33;
this.m34 = m34;
this.m41 = m41;
this.m42 = m42;
this.m43 = m43;
this.m44 = m44;
};
CanvasMatrix4.prototype.multLeft = function(mat)
{
var m11 = (mat.m11 * this.m11 + mat.m12 * this.m21 +
mat.m13 * this.m31 + mat.m14 * this.m41);
var m12 = (mat.m11 * this.m12 + mat.m12 * this.m22 +
mat.m13 * this.m32 + mat.m14 * this.m42);
var m13 = (mat.m11 * this.m13 + mat.m12 * this.m23 +
mat.m13 * this.m33 + mat.m14 * this.m43);
var m14 = (mat.m11 * this.m14 + mat.m12 * this.m24 +
mat.m13 * this.m34 + mat.m14 * this.m44);
var m21 = (mat.m21 * this.m11 + mat.m22 * this.m21 +
mat.m23 * this.m31 + mat.m24 * this.m41);
var m22 = (mat.m21 * this.m12 + mat.m22 * this.m22 +
mat.m23 * this.m32 + mat.m24 * this.m42);
var m23 = (mat.m21 * this.m13 + mat.m22 * this.m23 +
mat.m23 * this.m33 + mat.m24 * this.m43);
var m24 = (mat.m21 * this.m14 + mat.m22 * this.m24 +
mat.m23 * this.m34 + mat.m24 * this.m44);
var m31 = (mat.m31 * this.m11 + mat.m32 * this.m21 +
mat.m33 * this.m31 + mat.m34 * this.m41);
var m32 = (mat.m31 * this.m12 + mat.m32 * this.m22 +
mat.m33 * this.m32 + mat.m34 * this.m42);
var m33 = (mat.m31 * this.m13 + mat.m32 * this.m23 +
mat.m33 * this.m33 + mat.m34 * this.m43);
var m34 = (mat.m31 * this.m14 + mat.m32 * this.m24 +
mat.m33 * this.m34 + mat.m34 * this.m44);
var m41 = (mat.m41 * this.m11 + mat.m42 * this.m21 +
mat.m43 * this.m31 + mat.m44 * this.m41);
var m42 = (mat.m41 * this.m12 + mat.m42 * this.m22 +
mat.m43 * this.m32 + mat.m44 * this.m42);
var m43 = (mat.m41 * this.m13 + mat.m42 * this.m23 +
mat.m43 * this.m33 + mat.m44 * this.m43);
var m44 = (mat.m41 * this.m14 + mat.m42 * this.m24 +
mat.m43 * this.m34 + mat.m44 * this.m44);
this.m11 = m11;
this.m12 = m12;
this.m13 = m13;
this.m14 = m14;
this.m21 = m21;
this.m22 = m22;
this.m23 = m23;
this.m24 = m24;
this.m31 = m31;
this.m32 = m32;
this.m33 = m33;
this.m34 = m34;
this.m41 = m41;
this.m42 = m42;
this.m43 = m43;
this.m44 = m44;
};
CanvasMatrix4.prototype.ortho = function(left, right, bottom, top, near, far)
{
var tx = (left + right) / (left - right);
var ty = (top + bottom) / (bottom - top);
var tz = (far + near) / (near - far);
var matrix = new CanvasMatrix4();
matrix.m11 = 2 / (right - left);
matrix.m12 = 0;
matrix.m13 = 0;
matrix.m14 = 0;
matrix.m21 = 0;
matrix.m22 = 2 / (top - bottom);
matrix.m23 = 0;
matrix.m24 = 0;
matrix.m31 = 0;
matrix.m32 = 0;
matrix.m33 = -2 / (far - near);
matrix.m34 = 0;
matrix.m41 = tx;
matrix.m42 = ty;
matrix.m43 = tz;
matrix.m44 = 1;
this.multRight(matrix);
};
CanvasMatrix4.prototype.frustum = function(left, right, bottom, top, near, far)
{
var matrix = new CanvasMatrix4();
var A = (right + left) / (right - left);
var B = (top + bottom) / (top - bottom);
var C = -(far + near) / (far - near);
var D = -(2 * far * near) / (far - near);
matrix.m11 = (2 * near) / (right - left);
matrix.m12 = 0;
matrix.m13 = 0;
matrix.m14 = 0;
matrix.m21 = 0;
matrix.m22 = 2 * near / (top - bottom);
matrix.m23 = 0;
matrix.m24 = 0;
matrix.m31 = A;
matrix.m32 = B;
matrix.m33 = C;
matrix.m34 = -1;
matrix.m41 = 0;
matrix.m42 = 0;
matrix.m43 = D;
matrix.m44 = 0;
this.multRight(matrix);
};
CanvasMatrix4.prototype.perspective = function(fovy, aspect, zNear, zFar)
{
var top = Math.tan(fovy * Math.PI / 360) * zNear;
var bottom = -top;
var left = aspect * bottom;
var right = aspect * top;
this.frustum(left, right, bottom, top, zNear, zFar);
};
CanvasMatrix4.prototype.lookat = function(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz)
{
var matrix = new CanvasMatrix4();
// Make rotation matrix
// Z vector
var zx = eyex - centerx;
var zy = eyey - centery;
var zz = eyez - centerz;
var mag = Math.sqrt(zx * zx + zy * zy + zz * zz);
if (mag) {
zx /= mag;
zy /= mag;
zz /= mag;
}
// Y vector
var yx = upx;
var yy = upy;
var yz = upz;
// X vector = Y cross Z
xx =  yy * zz - yz * zy;
xy = -yx * zz + yz * zx;
xz =  yx * zy - yy * zx;
// Recompute Y = Z cross X
yx = zy * xz - zz * xy;
yy = -zx * xz + zz * xx;
yx = zx * xy - zy * xx;
// cross product gives area of parallelogram, which is < 1.0 for
// non-perpendicular unit-length vectors; so normalize x, y here
mag = Math.sqrt(xx * xx + xy * xy + xz * xz);
if (mag) {
xx /= mag;
xy /= mag;
xz /= mag;
}
mag = Math.sqrt(yx * yx + yy * yy + yz * yz);
if (mag) {
yx /= mag;
yy /= mag;
yz /= mag;
}
matrix.m11 = xx;
matrix.m12 = xy;
matrix.m13 = xz;
matrix.m14 = 0;
matrix.m21 = yx;
matrix.m22 = yy;
matrix.m23 = yz;
matrix.m24 = 0;
matrix.m31 = zx;
matrix.m32 = zy;
matrix.m33 = zz;
matrix.m34 = 0;
matrix.m41 = 0;
matrix.m42 = 0;
matrix.m43 = 0;
matrix.m44 = 1;
matrix.translate(-eyex, -eyey, -eyez);
this.multRight(matrix);
};
// Support functions
CanvasMatrix4.prototype._determinant2x2 = function(a, b, c, d)
{
return a * d - b * c;
};
CanvasMatrix4.prototype._determinant3x3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3)
{
return a1 * this._determinant2x2(b2, b3, c2, c3) -
b1 * this._determinant2x2(a2, a3, c2, c3) +
c1 * this._determinant2x2(a2, a3, b2, b3);
};
CanvasMatrix4.prototype._determinant4x4 = function()
{
var a1 = this.m11;
var b1 = this.m12;
var c1 = this.m13;
var d1 = this.m14;
var a2 = this.m21;
var b2 = this.m22;
var c2 = this.m23;
var d2 = this.m24;
var a3 = this.m31;
var b3 = this.m32;
var c3 = this.m33;
var d3 = this.m34;
var a4 = this.m41;
var b4 = this.m42;
var c4 = this.m43;
var d4 = this.m44;
return a1 * this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4) -
b1 * this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4) +
c1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4) -
d1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
};
CanvasMatrix4.prototype._makeAdjoint = function()
{
var a1 = this.m11;
var b1 = this.m12;
var c1 = this.m13;
var d1 = this.m14;
var a2 = this.m21;
var b2 = this.m22;
var c2 = this.m23;
var d2 = this.m24;
var a3 = this.m31;
var b3 = this.m32;
var c3 = this.m33;
var d3 = this.m34;
var a4 = this.m41;
var b4 = this.m42;
var c4 = this.m43;
var d4 = this.m44;
// Row column labeling reversed since we transpose rows & columns
this.m11  =   this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
this.m21  = - this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
this.m31  =   this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
this.m41  = - this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
this.m12  = - this._determinant3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
this.m22  =   this._determinant3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
this.m32  = - this._determinant3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
this.m42  =   this._determinant3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);
this.m13  =   this._determinant3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
this.m23  = - this._determinant3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
this.m33  =   this._determinant3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
this.m43  = - this._determinant3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);
this.m14  = - this._determinant3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
this.m24  =   this._determinant3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
this.m34  = - this._determinant3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
this.m44  =   this._determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
};</script>
<script>// To generate the help pages for this library, use
// jsdoc --destination ../../../doc/rglwidgetClass --template ~/node_modules/jsdoc-baseline rglClass.src.js
// To validate, use
// setwd(".../inst/htmlwidgets/lib/rglClass")
// hints <- js::jshint(readLines("rglClass.src.js"))
// hints[, c("line", "reason")]
/**
* The class of an rgl widget
* @class
*/
rglwidgetClass = function() {
this.canvas = null;
this.userMatrix = new CanvasMatrix4();
this.types = [];
this.prMatrix = new CanvasMatrix4();
this.mvMatrix = new CanvasMatrix4();
this.vp = null;
this.prmvMatrix = null;
this.origs = null;
this.gl = null;
this.scene = null;
this.select = {state: "inactive", subscene: null, region: {p1: {x:0, y:0}, p2: {x:0, y:0}}};
this.drawing = false;
};
/**
* Multiply matrix by vector
* @returns {number[]}
* @param M {number[][]} Left operand
* @param v {number[]} Right operand
*/
rglwidgetClass.prototype.multMV = function(M, v) {
return [ M.m11 * v[0] + M.m12 * v[1] + M.m13 * v[2] + M.m14 * v[3],
M.m21 * v[0] + M.m22 * v[1] + M.m23 * v[2] + M.m24 * v[3],
M.m31 * v[0] + M.m32 * v[1] + M.m33 * v[2] + M.m34 * v[3],
M.m41 * v[0] + M.m42 * v[1] + M.m43 * v[2] + M.m44 * v[3]
];
};
/**
* Multiply row vector by Matrix
* @returns {number[]}
* @param v {number[]} left operand
* @param M {number[][]} right operand
*/
rglwidgetClass.prototype.multVM = function(v, M) {
return [ M.m11 * v[0] + M.m21 * v[1] + M.m31 * v[2] + M.m41 * v[3],
M.m12 * v[0] + M.m22 * v[1] + M.m32 * v[2] + M.m42 * v[3],
M.m13 * v[0] + M.m23 * v[1] + M.m33 * v[2] + M.m43 * v[3],
M.m14 * v[0] + M.m24 * v[1] + M.m34 * v[2] + M.m44 * v[3]
];
};
/**
* Euclidean length of a vector
* @returns {number}
* @param v {number[]}
*/
rglwidgetClass.prototype.vlen = function(v) {
return Math.sqrt(this.dotprod(v, v));
};
/**
* Dot product of two vectors
* @instance rglwidgetClass
* @returns {number}
* @param a {number[]}
* @param b {number[]}
*/
rglwidgetClass.prototype.dotprod = function(a, b) {
return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
};
/**
* Cross product of two vectors
* @returns {number[]}
* @param a {number[]}
* @param b {number[]}
*/
rglwidgetClass.prototype.xprod = function(a, b) {
return [a[1]*b[2] - a[2]*b[1],
a[2]*b[0] - a[0]*b[2],
a[0]*b[1] - a[1]*b[0]];
};
/**
* Bind vectors or matrices by columns
* @returns {number[][]}
* @param a {number[]|number[][]}
* @param b {number[]|number[][]}
*/
rglwidgetClass.prototype.cbind = function(a, b) {
if (b.length < a.length)
b = this.repeatToLen(b, a.length);
else if (a.length < b.length)
a = this.repeatToLen(a, b.length);
return a.map(function(currentValue, index, array) {
return currentValue.concat(b[index]);
});
};
/**
* Swap elements
* @returns {any[]}
* @param a {any[]}
* @param i {number} Element to swap
* @param j {number} Other element to swap
*/
rglwidgetClass.prototype.swap = function(a, i, j) {
var temp = a[i];
a[i] = a[j];
a[j] = temp;
};
/**
* Flatten a matrix into a vector
* @returns {any[]}
* @param a {any[][]}
*/
rglwidgetClass.prototype.flatten = function(arr, result) {
var value;
if (typeof result === "undefined") result = [];
for (var i = 0, length = arr.length; i < length; i++) {
value = arr[i];
if (Array.isArray(value)) {
this.flatten(value, result);
} else {
result.push(value);
}
}
return result;
};
/**
* set element of 1d or 2d array as if it was flattened.
* Column major, zero based!
* @returns {any[]|any[][]}
* @param {any[]|any[][]} a - array
* @param {number} i - element
* @param {any} value
*/
rglwidgetClass.prototype.setElement = function(a, i, value) {
if (Array.isArray(a[0])) {
var dim = a.length,
col = Math.floor(i/dim),
row = i % dim;
a[row][col] = value;
} else {
a[i] = value;
}
};
/**
* Transpose an array
* @returns {any[][]}
* @param {any[][]} a
*/
rglwidgetClass.prototype.transpose = function(a) {
var newArray = [],
n = a.length,
m = a[0].length,
i;
for(i = 0; i < m; i++){
newArray.push([]);
}
for(i = 0; i < n; i++){
for(var j = 0; j < m; j++){
newArray[j].push(a[i][j]);
}
}
return newArray;
};
/**
* Calculate sum of squares of a numeric vector
* @returns {number}
* @param {number[]} x
*/
rglwidgetClass.prototype.sumsq = function(x) {
var result = 0, i;
for (i=0; i < x.length; i++)
result += x[i]*x[i];
return result;
};
/**
* Convert a matrix to a CanvasMatrix4
* @returns {CanvasMatrix4}
* @param {number[][]|number[]} mat
*/
rglwidgetClass.prototype.toCanvasMatrix4 = function(mat) {
if (mat instanceof CanvasMatrix4)
return mat;
var result = new CanvasMatrix4();
mat = this.flatten(this.transpose(mat));
result.load(mat);
return result;
};
/**
* Convert an R-style numeric colour string to an rgb vector
* @returns {number[]}
* @param {string} s
*/
rglwidgetClass.prototype.stringToRgb = function(s) {
s = s.replace("#", "");
var bigint = parseInt(s, 16);
return [((bigint >> 16) & 255)/255,
((bigint >> 8) & 255)/255,
(bigint & 255)/255];
};
/**
* Take a component-by-component product of two 3 vectors
* @returns {number[]}
* @param {number[]} x
* @param {number[]} y
*/
rglwidgetClass.prototype.componentProduct = function(x, y) {
if (typeof y === "undefined") {
this.alertOnce("Bad arg to componentProduct");
}
var result = new Float32Array(3), i;
for (i = 0; i<3; i++)
result[i] = x[i]*y[i];
return result;
};
/**
* Get next higher power of two
* @returns { number }
* @param { number } value - input value
*/
rglwidgetClass.prototype.getPowerOfTwo = function(value) {
var pow = 1;
while(pow<value) {
pow *= 2;
}
return pow;
};
/**
* Unique entries
* @returns { any[] }
* @param { any[] } arr - An array
*/
rglwidgetClass.prototype.unique = function(arr) {
arr = [].concat(arr);
return arr.filter(function(value, index, self) {
return self.indexOf(value) === index;
});
};
/**
* Shallow compare of arrays
* @returns { boolean }
* @param { any[] } a - An array
* @param { any[] } b - Another array
*/
rglwidgetClass.prototype.equalArrays = function(a, b) {
return a === b || (a && b &&
a.length === b.length &&
a.every(function(v, i) {return v === b[i];}));
};
/**
* Repeat an array to a desired length
* @returns {any[]}
* @param {any | any[]} arr The input array
* @param {number} len The desired output length
*/
rglwidgetClass.prototype.repeatToLen = function(arr, len) {
arr = [].concat(arr);
while (arr.length < len/2)
arr = arr.concat(arr);
return arr.concat(arr.slice(0, len - arr.length));
};
/**
* Give a single alert message, not to be repeated.
* @param {string} msg  The message to give.
*/
rglwidgetClass.prototype.alertOnce = function(msg) {
if (typeof this.alerted !== "undefined")
return;
this.alerted = true;
alert(msg);
};
rglwidgetClass.prototype.f_is_lit = 1;
rglwidgetClass.prototype.f_is_smooth = 2;
rglwidgetClass.prototype.f_has_texture = 4;
rglwidgetClass.prototype.f_depth_sort = 8;
rglwidgetClass.prototype.f_fixed_quads = 16;
rglwidgetClass.prototype.f_is_transparent = 32;
rglwidgetClass.prototype.f_is_lines = 64;
rglwidgetClass.prototype.f_sprites_3d = 128;
rglwidgetClass.prototype.f_sprite_3d = 256;
rglwidgetClass.prototype.f_is_subscene = 512;
rglwidgetClass.prototype.f_is_clipplanes = 1024;
rglwidgetClass.prototype.f_fixed_size = 2048;
rglwidgetClass.prototype.f_is_points = 4096;
rglwidgetClass.prototype.f_is_twosided = 8192;
rglwidgetClass.prototype.f_fat_lines = 16384;
rglwidgetClass.prototype.f_is_brush = 32768;
/**
* Which list does a particular id come from?
* @returns { string }
* @param {number} id The id to look up.
*/
rglwidgetClass.prototype.whichList = function(id) {
var obj = this.getObj(id),
flags = obj.flags;
if (obj.type === "light")
return "lights";
if (flags & this.f_is_subscene)
return "subscenes";
if (flags & this.f_is_clipplanes)
return "clipplanes";
if (flags & this.f_is_transparent)
return "transparent";
return "opaque";
};
/**
* Get an object by id number.
* @returns { Object }
* @param {number} id
*/
rglwidgetClass.prototype.getObj = function(id) {
if (typeof id !== "number") {
this.alertOnce("getObj id is "+typeof id);
}
return this.scene.objects[id];
};
/**
* Get ids of a particular type from a subscene or the whole scene
* @returns { number[] }
* @param {string} type What type of object?
* @param {number} subscene  Which subscene?  If not given, find in the whole scene
*/
rglwidgetClass.prototype.getIdsByType = function(type, subscene) {
var
result = [], i, self = this;
if (typeof subscene === "undefined") {
Object.keys(this.scene.objects).forEach(
function(key) {
key = parseInt(key, 10);
if (self.getObj(key).type === type)
result.push(key);
});
} else {
ids = this.getObj(subscene).objects;
for (i=0; i < ids.length; i++) {
if (this.getObj(ids[i]).type === type) {
result.push(ids[i]);
}
}
}
return result;
};
/**
* Get a particular material property for an id
* @returns { any }
* @param {number} id  Which object?
* @param {string} property Which material property?
*/
rglwidgetClass.prototype.getMaterial = function(id, property) {
var obj = this.getObj(id),
mat = obj.material[property];
if (typeof mat === "undefined")
mat = this.scene.material[property];
return mat;
};
/**
* Is a particular id in a subscene?
* @returns { boolean }
* @param {number} id Which id?
* @param {number} subscene Which subscene id?
*/
rglwidgetClass.prototype.inSubscene = function(id, subscene) {
return this.getObj(subscene).objects.indexOf(id) > -1;
};
/**
* Add an id to a subscene.
* @param {number} id Which id?
* @param {number} subscene Which subscene id?
*/
rglwidgetClass.prototype.addToSubscene = function(id, subscene) {
var thelist,
thesub = this.getObj(subscene),
ids = [id],
obj = this.getObj(id), i;
if (typeof obj != "undefined" && typeof (obj.newIds) !== "undefined") {
ids = ids.concat(obj.newIds);
}
thesub.objects = [].concat(thesub.objects);
for (i = 0; i < ids.length; i++) {
id = ids[i];
if (thesub.objects.indexOf(id) == -1) {
thelist = this.whichList(id);
thesub.objects.push(id);
thesub[thelist].push(id);
}
}
};
/**
* Delete an id from a subscene
* @param { number } id - the id to add
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.delFromSubscene = function(id, subscene) {
var thelist,
thesub = this.getObj(subscene),
obj = this.getObj(id),
ids = [id], i;
if (typeof obj !== "undefined" && typeof (obj.newIds) !== "undefined")
ids = ids.concat(obj.newIds);
thesub.objects = [].concat(thesub.objects); // It might be a scalar
for (j=0; j<ids.length;j++) {
id = ids[j];
i = thesub.objects.indexOf(id);
if (i > -1) {
thesub.objects.splice(i, 1);
thelist = this.whichList(id);
i = thesub[thelist].indexOf(id);
thesub[thelist].splice(i, 1);
}
}
};
/**
* Set the ids in a subscene
* @param { number[] } ids - the ids to set
* @param { number } subsceneid - the id of the subscene
*/
rglwidgetClass.prototype.setSubsceneEntries = function(ids, subsceneid) {
var sub = this.getObj(subsceneid);
sub.objects = ids;
this.initSubscene(subsceneid);
};
/**
* Get the ids in a subscene
* @returns {number[]}
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.getSubsceneEntries = function(subscene) {
return this.getObj(subscene).objects;
};
/**
* Get the ids of the subscenes within a subscene
* @returns { number[] }
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.getChildSubscenes = function(subscene) {
return this.getObj(subscene).subscenes;
};
/**
* Start drawing
* @returns { boolean } Previous state
*/
rglwidgetClass.prototype.startDrawing = function() {
var value = this.drawing;
this.drawing = true;
return value;
};
/**
* Stop drawing and check for context loss
* @param { boolean } saved - Previous state
*/
rglwidgetClass.prototype.stopDrawing = function(saved) {
this.drawing = saved;
if (!saved && this.gl && this.gl.isContextLost())
this.restartCanvas();
};
/**
* Generate the vertex shader for an object
* @returns {string}
* @param { number } id - Id of object
*/
rglwidgetClass.prototype.getVertexShader = function(id) {
var obj = this.getObj(id),
userShader = obj.userVertexShader,
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
nclipplanes = this.countClipplanes(),
fixed_size = flags & this.f_fixed_size,
is_points = flags & this.f_is_points,
is_twosided = flags & this.f_is_twosided,
fat_lines = flags & this.f_fat_lines,
is_brush = flags & this.f_is_brush,
result;
if (type === "clipplanes" || sprites_3d) return;
if (typeof userShader !== "undefined") return userShader;
result = "  /* ****** "+type+" object "+id+" vertex shader ****** */\n"+
"  attribute vec3 aPos;\n"+
"  attribute vec4 aCol;\n"+
" uniform mat4 mvMatrix;\n"+
" uniform mat4 prMatrix;\n"+
" varying vec4 vCol;\n"+
" varying vec4 vPosition;\n";
if ((is_lit && !fixed_quads && !is_brush) || sprite_3d)
result = result + "  attribute vec3 aNorm;\n"+
" uniform mat4 normMatrix;\n"+
" varying vec3 vNormal;\n";
if (has_texture || type === "text")
result = result + " attribute vec2 aTexcoord;\n"+
" varying vec2 vTexcoord;\n";
if (fixed_size)
result = result + "  uniform vec2 textScale;\n";
if (fixed_quads)
result = result + "  attribute vec2 aOfs;\n";
else if (sprite_3d)
result = result + "  uniform vec3 uOrig;\n"+
"  uniform float uSize;\n"+
"  uniform mat4 usermat;\n";
if (is_twosided)
result = result + "  attribute vec3 aPos1;\n"+
"  attribute vec3 aPos2;\n"+
"  varying float normz;\n";
if (fat_lines) {
result = result +   "  attribute vec3 aNext;\n"+
"  attribute vec2 aPoint;\n"+
"  varying vec2 vPoint;\n"+
"  varying float vLength;\n"+
"  uniform float uAspect;\n"+
"  uniform float uLwd;\n";
}
result = result + "  void main(void) {\n";
if ((nclipplanes || (!fixed_quads && !sprite_3d)) && !is_brush)
result = result + "    vPosition = mvMatrix * vec4(aPos, 1.);\n";
if (!fixed_quads && !sprite_3d && !is_brush)
result = result + "    gl_Position = prMatrix * vPosition;\n";
if (is_points) {
var size = this.getMaterial(id, "size");
result = result + "    gl_PointSize = "+size.toFixed(1)+";\n";
}
result = result + "    vCol = aCol;\n";
if (is_lit && !fixed_quads && !sprite_3d && !is_brush)
result = result + "    vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n";
if (has_texture || type == "text")
result = result + "    vTexcoord = aTexcoord;\n";
if (fixed_size)
result = result + "    vec4 pos = prMatrix * mvMatrix * vec4(aPos, 1.);\n"+
"   pos = pos/pos.w;\n"+
"   gl_Position = pos + vec4(aOfs*textScale, 0.,0.);\n";
if (type == "sprites" && !fixed_size)
result = result + "    vec4 pos = mvMatrix * vec4(aPos, 1.);\n"+
"   pos = pos/pos.w + vec4(aOfs, 0., 0.);\n"+
"   gl_Position = prMatrix*pos;\n";
if (sprite_3d)
result = result + "   vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n"+
"   vec4 pos = mvMatrix * vec4(uOrig, 1.);\n"+
"   vPosition = pos/pos.w + vec4(uSize*(vec4(aPos, 1.)*usermat).xyz,0.);\n"+
"   gl_Position = prMatrix * vPosition;\n";
if (is_twosided)
result = result + "   vec4 pos1 = prMatrix*(mvMatrix*vec4(aPos1, 1.));\n"+
"   pos1 = pos1/pos1.w - gl_Position/gl_Position.w;\n"+
"   vec4 pos2 = prMatrix*(mvMatrix*vec4(aPos2, 1.));\n"+
"   pos2 = pos2/pos2.w - gl_Position/gl_Position.w;\n"+
"   normz = pos1.x*pos2.y - pos1.y*pos2.x;\n";
if (fat_lines) 
/* This code was inspired by Matt Deslauriers' code in https://mattdesl.svbtle.com/drawing-lines-is-hard */
result = result + "   vec2 aspectVec = vec2(uAspect, 1.0);\n"+
"   mat4 projViewModel = prMatrix * mvMatrix;\n"+
"   vec4 currentProjected = projViewModel * vec4(aPos, 1.0);\n"+
"   currentProjected = currentProjected/currentProjected.w;\n"+
"   vec4 nextProjected = projViewModel * vec4(aNext, 1.0);\n"+
"   vec2 currentScreen = currentProjected.xy * aspectVec;\n"+
"   vec2 nextScreen = (nextProjected.xy / nextProjected.w) * aspectVec;\n"+
"   float len = uLwd;\n"+
"   vec2 dir = vec2(1.0, 0.0);\n"+
"   vPoint = aPoint;\n"+
"   vLength = length(nextScreen - currentScreen)/2.0;\n"+
"   vLength = vLength/(vLength + len);\n"+
"   if (vLength > 0.0) {\n"+
"     dir = normalize(nextScreen - currentScreen);\n"+
"   }\n"+
"   vec2 normal = vec2(-dir.y, dir.x);\n"+
"   dir.x /= uAspect;\n"+
"   normal.x /= uAspect;\n"+
"   vec4 offset = vec4(len*(normal*aPoint.x*aPoint.y - dir), 0.0, 0.0);\n"+
"   gl_Position = currentProjected + offset;\n";
if (is_brush)
result = result + "   gl_Position = vec4(aPos, 1.);\n";
result = result + "  }\n";
// console.log(result);
return result;
};
/**
* Generate the fragment shader for an object
* @returns {string}
* @param { number } id - Id of object
*/
rglwidgetClass.prototype.getFragmentShader = function(id) {
var obj = this.getObj(id),
userShader = obj.userFragmentShader,
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
sprites_3d = flags & this.f_sprites_3d,
is_twosided = (flags & this.f_is_twosided) > 0,
fat_lines = flags & this.f_fat_lines,
is_transparent = flags & this.f_is_transparent,
nclipplanes = this.countClipplanes(), i,
texture_format, nlights,
result;
if (type === "clipplanes" || sprites_3d) return;
if (typeof userShader !== "undefined") return userShader;
if (has_texture)
texture_format = this.getMaterial(id, "textype");
result = "/* ****** "+type+" object "+id+" fragment shader ****** */\n"+
"#ifdef GL_ES\n"+
"  precision highp float;\n"+
"#endif\n"+
"  varying vec4 vCol; // carries alpha\n"+
"  varying vec4 vPosition;\n";
if (has_texture || type === "text")
result = result + "  varying vec2 vTexcoord;\n"+
" uniform sampler2D uSampler;\n";
if (is_lit && !fixed_quads)
result = result + "  varying vec3 vNormal;\n";
for (i = 0; i < nclipplanes; i++)
result = result + "  uniform vec4 vClipplane"+i+";\n";
if (is_lit) {
nlights = this.countLights();
if (nlights)
result = result + "  uniform mat4 mvMatrix;\n";
else
is_lit = false;
}
if (is_lit) {
result = result + "   uniform vec3 emission;\n"+
"   uniform float shininess;\n";
for (i=0; i < nlights; i++) {
result = result + "   uniform vec3 ambient" + i + ";\n"+
"   uniform vec3 specular" + i +"; // light*material\n"+
"   uniform vec3 diffuse" + i + ";\n"+
"   uniform vec3 lightDir" + i + ";\n"+
"   uniform bool viewpoint" + i + ";\n"+
"   uniform bool finite" + i + ";\n";
}
}
if (is_twosided)
result = result + "   uniform bool front;\n"+
"   varying float normz;\n";
if (fat_lines)
result = result + "   varying vec2 vPoint;\n"+
"   varying float vLength;\n";
result = result + "  void main(void) {\n";
if (fat_lines) {
result = result + "    vec2 point = vPoint;\n"+
"    bool neg = point.y < 0.0;\n"+
"    point.y = neg ? "+
"      (point.y + vLength)/(1.0 - vLength) :\n"+
"     -(point.y - vLength)/(1.0 - vLength);\n";
if (is_transparent && type == "linestrip")
result = result+"    if (neg && length(point) <= 1.0) discard;\n";
result = result + "    point.y = min(point.y, 0.0);\n"+
"    if (length(point) > 1.0) discard;\n";
}
for (i=0; i < nclipplanes;i++)
result = result + "    if (dot(vPosition, vClipplane"+i+") < 0.0) discard;\n";
if (fixed_quads) {
result = result +   "    vec3 n = vec3(0., 0., 1.);\n";
} else if (is_lit) {
result = result +   "    vec3 n = normalize(vNormal);\n";
}
if (is_twosided) {
result = result +   "    if ((normz <= 0.) != front) discard;\n";
}
if (is_lit) {
result = result + "    vec3 eye = normalize(-vPosition.xyz);\n"+
"   vec3 lightdir;\n"+
"   vec4 colDiff;\n"+
"   vec3 halfVec;\n"+
"   vec4 lighteffect = vec4(emission, 0.);\n"+
"   vec3 col;\n"+
"   float nDotL;\n";
if (!fixed_quads) {
result = result +   "   n = -faceforward(n, n, eye);\n";
}
for (i=0; i < nlights; i++) {
result = result + "   colDiff = vec4(vCol.rgb * diffuse" + i + ", vCol.a);\n"+
"   lightdir = lightDir" + i + ";\n"+
"   if (!viewpoint" + i +")\n"+
"     lightdir = (mvMatrix * vec4(lightdir, 1.)).xyz;\n"+
"   if (!finite" + i + ") {\n"+
"     halfVec = normalize(lightdir + eye);\n"+
"   } else {\n"+
"     lightdir = normalize(lightdir - vPosition.xyz);\n"+
"     halfVec = normalize(lightdir + eye);\n"+
"   }\n"+
"    col = ambient" + i + ";\n"+
"   nDotL = dot(n, lightdir);\n"+
"   col = col + max(nDotL, 0.) * colDiff.rgb;\n"+
"   col = col + pow(max(dot(halfVec, n), 0.), shininess) * specular" + i + ";\n"+
"   lighteffect = lighteffect + vec4(col, colDiff.a);\n";
}
} else {
result = result +   "   vec4 colDiff = vCol;\n"+
"    vec4 lighteffect = colDiff;\n";
}
if (type === "text")
result = result +   "    vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n";
if (has_texture) {
result = result + {
rgb:            "   vec4 textureColor = lighteffect*vec4(texture2D(uSampler, vTexcoord).rgb, 1.);\n",
rgba:           "   vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n",
alpha:          "   vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
"   float luminance = dot(vec3(1.,1.,1.), textureColor.rgb)/3.;\n"+
"   textureColor =  vec4(lighteffect.rgb, lighteffect.a*luminance);\n",
luminance:      "   vec4 textureColor = vec4(lighteffect.rgb*dot(texture2D(uSampler, vTexcoord).rgb, vec3(1.,1.,1.))/3., lighteffect.a);\n",
"luminance.alpha":"    vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
"   float luminance = dot(vec3(1.,1.,1.),textureColor.rgb)/3.;\n"+
"   textureColor = vec4(lighteffect.rgb*luminance, lighteffect.a*textureColor.a);\n"
}[texture_format]+
"   gl_FragColor = textureColor;\n";
} else if (type === "text") {
result = result +   "    if (textureColor.a < 0.1)\n"+
"     discard;\n"+
"   else\n"+
"     gl_FragColor = textureColor;\n";
} else
result = result +   "   gl_FragColor = lighteffect;\n";
//if (fat_lines)
//  result = result +   "   gl_FragColor = vec4(0.0, abs(point.x), abs(point.y), 1.0);"
result = result + "  }\n";
// console.log(result);
return result;
};
/**
* Call gl functions to create and compile shader
* @returns {Object}
* @param { number } shaderType - gl code for shader type
* @param { string } code - code for the shader
*/
rglwidgetClass.prototype.getShader = function(shaderType, code) {
var gl = this.gl, shader;
shader = gl.createShader(shaderType);
gl.shaderSource(shader, code);
gl.compileShader(shader);
if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS) && !gl.isContextLost())
alert(gl.getShaderInfoLog(shader));
return shader;
};
/**
* Handle a texture after its image has been loaded
* @param { Object } texture - the gl texture object
* @param { Object } textureCanvas - the canvas holding the image
*/
rglwidgetClass.prototype.handleLoadedTexture = function(texture, textureCanvas) {
var gl = this.gl || this.initGL();
gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
gl.bindTexture(gl.TEXTURE_2D, texture);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, textureCanvas);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
gl.generateMipmap(gl.TEXTURE_2D);
gl.bindTexture(gl.TEXTURE_2D, null);
};
/**
* Get maximum dimension of texture in current browser.
* @returns {number}
*/
rglwidgetClass.prototype.getMaxTexSize = function() {
var gl = this.gl || this.initGL();	
return Math.min(4096, gl.getParameter(gl.MAX_TEXTURE_SIZE));
};
/**
* Load an image to a texture
* @param { string } uri - The image location
* @param { Object } texture - the gl texture object
*/
rglwidgetClass.prototype.loadImageToTexture = function(uri, texture) {
var canvas = this.textureCanvas,
ctx = canvas.getContext("2d"),
image = new Image(),
self = this;
image.onload = function() {
var w = image.width,
h = image.height,
canvasX = self.getPowerOfTwo(w),
canvasY = self.getPowerOfTwo(h),
gl = self.gl || self.initGL(),
maxTexSize = self.getMaxTexSize();
while (canvasX > 1 && canvasY > 1 && (canvasX > maxTexSize || canvasY > maxTexSize)) {
canvasX /= 2;
canvasY /= 2;
}
canvas.width = canvasX;
canvas.height = canvasY;
ctx.imageSmoothingEnabled = true;
ctx.drawImage(image, 0, 0, canvasX, canvasY);
self.handleLoadedTexture(texture, canvas);
self.drawScene();
};
image.src = uri;
};
/**
* Draw text to the texture canvas
* @returns { Object } object with text measurements
* @param { string } text - the text
* @param { number } cex - expansion
* @param { string } family - font family
* @param { number } font - font number
*/
rglwidgetClass.prototype.drawTextToCanvas = function(text, cex, family, font) {
var canvasX, canvasY,
textY,
scaling = 20,
textColour = "white",
backgroundColour = "rgba(0,0,0,0)",
canvas = this.textureCanvas,
ctx = canvas.getContext("2d"),
i, textHeight = 0, textHeights = [], width, widths = [], 
offsetx, offsety = 0, line, lines = [], offsetsx = [],
offsetsy = [], lineoffsetsy = [], fontStrings = [],
maxTexSize = this.getMaxTexSize(),
getFontString = function(i) {
textHeights[i] = scaling*cex[i];
var fontString = textHeights[i] + "px",
family0 = family[i],
font0 = font[i];
if (family0 === "sans")
family0 = "sans-serif";
else if (family0 === "mono")
family0 = "monospace";
fontString = fontString + " " + family0;
if (font0 === 2 || font0 === 4)
fontString = "bold " + fontString;
if (font0 === 3 || font0 === 4)
fontString = "italic " + fontString;
return fontString;
};
cex = this.repeatToLen(cex, text.length);
family = this.repeatToLen(family, text.length);
font = this.repeatToLen(font, text.length);
canvasX = 1;
line = -1;
offsetx = maxTexSize;
for (i = 0; i < text.length; i++)  {
ctx.font = fontStrings[i] = getFontString(i);
width = widths[i] = ctx.measureText(text[i]).width;
if (offsetx + width > maxTexSize) {
line += 1;
offsety = lineoffsetsy[line] = offsety + 2*textHeight;
if (offsety > maxTexSize)
console.error("Too many strings for texture.");
textHeight = 0;
offsetx = 0;
}
textHeight = Math.max(textHeight, textHeights[i]);
offsetsx[i] = offsetx;
offsetx += width;
canvasX = Math.max(canvasX, offsetx);
lines[i] = line;
}
offsety = lineoffsetsy[line] = offsety + 2*textHeight;
for (i = 0; i < text.length; i++) {
offsetsy[i] = lineoffsetsy[lines[i]];
}
canvasX = this.getPowerOfTwo(canvasX);
canvasY = this.getPowerOfTwo(offsety);
canvas.width = canvasX;
canvas.height = canvasY;
ctx.fillStyle = backgroundColour;
ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);
ctx.textBaseline = "alphabetic";
for(i = 0; i < text.length; i++) {
ctx.font = fontStrings[i];
ctx.fillStyle = textColour;
ctx.textAlign = "left";
ctx.fillText(text[i], offsetsx[i],  offsetsy[i]);
}
return {canvasX:canvasX, canvasY:canvasY,
widths:widths, textHeights:textHeights,
offsetsx:offsetsx, offsetsy:offsetsy};
};
/**
* Set the gl viewport and scissor test
* @param { number } id - id of subscene
*/
rglwidgetClass.prototype.setViewport = function(id) {
var gl = this.gl || this.initGL(),
vp = this.getObj(id).par3d.viewport,
x = vp.x*this.canvas.width,
y = vp.y*this.canvas.height,
width = vp.width*this.canvas.width,
height = vp.height*this.canvas.height;
this.vp = {x:x, y:y, width:width, height:height};
gl.viewport(x, y, width, height);
gl.scissor(x, y, width, height);
gl.enable(gl.SCISSOR_TEST);
};
/**
* Set the projection matrix for a subscene
* @param { number } id - id of subscene
*/
rglwidgetClass.prototype.setprMatrix = function(id) {
var subscene = this.getObj(id),
embedding = subscene.embeddings.projection;
if (embedding === "replace")
this.prMatrix.makeIdentity();
else
this.setprMatrix(subscene.parent);
if (embedding === "inherit")
return;
// This is based on the Frustum::enclose code from geom.cpp
var bbox = subscene.par3d.bbox,
scale = subscene.par3d.scale,
ranges = [(bbox[1]-bbox[0])*scale[0]/2,
(bbox[3]-bbox[2])*scale[1]/2,
(bbox[5]-bbox[4])*scale[2]/2],
radius = Math.sqrt(this.sumsq(ranges))*1.1; // A bit bigger to handle labels
if (radius <= 0) radius = 1;
var observer = subscene.par3d.observer,
distance = observer[2],
FOV = subscene.par3d.FOV, ortho = FOV === 0,
t = ortho ? 1 : Math.tan(FOV*Math.PI/360),
near = distance - radius,
far = distance + radius,
hlen,
aspect = this.vp.width/this.vp.height,
z = subscene.par3d.zoom;
if (far < 0.0)
far = 1.0;
if (near < far/100.0)
near = far/100.0;
hlen = t*near;
if (ortho) {
if (aspect > 1)
this.prMatrix.ortho(-hlen*aspect*z, hlen*aspect*z,
-hlen*z, hlen*z, near, far);
else
this.prMatrix.ortho(-hlen*z, hlen*z,
-hlen*z/aspect, hlen*z/aspect,
near, far);
} else {
if (aspect > 1)
this.prMatrix.frustum(-hlen*aspect*z, hlen*aspect*z,
-hlen*z, hlen*z, near, far);
else
this.prMatrix.frustum(-hlen*z, hlen*z,
-hlen*z/aspect, hlen*z/aspect,
near, far);
}
};
/**
* Set the model-view matrix for a subscene
* @param { number } id - id of the subscene
*/
rglwidgetClass.prototype.setmvMatrix = function(id) {
var observer = this.getObj(id).par3d.observer;
this.mvMatrix.makeIdentity();
this.setmodelMatrix(id);
this.mvMatrix.translate(-observer[0], -observer[1], -observer[2]);
};
/**
* Set the model matrix for a subscene
* @param { number } id - id of the subscene
*/
rglwidgetClass.prototype.setmodelMatrix = function(id) {
var subscene = this.getObj(id),
embedding = subscene.embeddings.model;
if (embedding !== "inherit") {
var scale = subscene.par3d.scale,
bbox = subscene.par3d.bbox,
center = [(bbox[0]+bbox[1])/2,
(bbox[2]+bbox[3])/2,
(bbox[4]+bbox[5])/2];
this.mvMatrix.translate(-center[0], -center[1], -center[2]);
this.mvMatrix.scale(scale[0], scale[1], scale[2]);
this.mvMatrix.multRight( subscene.par3d.userMatrix );
}
if (embedding !== "replace")
this.setmodelMatrix(subscene.parent);
};
/**
* Set the normals matrix for a subscene
* @param { number } subsceneid - id of the subscene
*/
rglwidgetClass.prototype.setnormMatrix = function(subsceneid) {
var self = this,
recurse = function(id) {
var sub = self.getObj(id),
embedding = sub.embeddings.model;
if (embedding !== "inherit") {
var scale = sub.par3d.scale;
self.normMatrix.scale(1/scale[0], 1/scale[1], 1/scale[2]);
self.normMatrix.multRight(sub.par3d.userMatrix);
}
if (embedding !== "replace")
recurse(sub.parent);
};
self.normMatrix.makeIdentity();
recurse(subsceneid);
};
/**
* Set the combined projection-model-view matrix
*/
rglwidgetClass.prototype.setprmvMatrix = function() {
this.prmvMatrix = new CanvasMatrix4( this.mvMatrix );
this.prmvMatrix.multRight( this.prMatrix );
};
/**
* Count clipping planes in a scene
* @returns {number}
*/
rglwidgetClass.prototype.countClipplanes = function() {
return this.countObjs("clipplanes");
};
/**
* Count lights in a scene
* @returns { number }
*/
rglwidgetClass.prototype.countLights = function() {
return this.countObjs("light");
};
/**
* Count objects of specific type in a scene
* @returns { number }
* @param { string } type - Type of object to count
*/
rglwidgetClass.prototype.countObjs = function(type) {
var self = this,
bound = 0;
Object.keys(this.scene.objects).forEach(
function(key) {
if (self.getObj(parseInt(key, 10)).type === type)
bound = bound + 1;
});
return bound;
};
/**
* Initialize a subscene
* @param { number } id - id of subscene.
*/
rglwidgetClass.prototype.initSubscene = function(id) {
var sub = this.getObj(id),
i, obj;
if (sub.type !== "subscene")
return;
sub.par3d.userMatrix = this.toCanvasMatrix4(sub.par3d.userMatrix);
sub.par3d.listeners = [].concat(sub.par3d.listeners);
sub.backgroundId = undefined;
sub.subscenes = [];
sub.clipplanes = [];
sub.transparent = [];
sub.opaque = [];
sub.lights = [];
for (i=0; i < sub.objects.length; i++) {
obj = this.getObj(sub.objects[i]);
if (typeof obj === "undefined") {
sub.objects.splice(i, 1);
i--;
} else if (obj.type === "background")
sub.backgroundId = obj.id;
else
sub[this.whichList(obj.id)].push(obj.id);
}
};
/**
* Copy object
* @param { number } id - id of object to copy
* @param { string } reuse - Document id of scene to reuse
*/
rglwidgetClass.prototype.copyObj = function(id, reuse) {
var obj = this.getObj(id),
prev = document.getElementById(reuse);
if (prev !== null) {
prev = prev.rglinstance;
var
prevobj = prev.getObj(id),
fields = ["flags", "type",
"colors", "vertices", "centers",
"normals", "offsets",
"texts", "cex", "family", "font", "adj",
"material",
"radii",
"texcoords",
"userMatrix", "ids",
"dim",
"par3d", "userMatrix",
"viewpoint", "finite"],
i;
for (i = 0; i < fields.length; i++) {
if (typeof prevobj[fields[i]] !== "undefined")
obj[fields[i]] = prevobj[fields[i]];
}
} else
console.warn("copyObj failed");
};
/**
* Update the triangles used to display a plane
* @param { number } id - id of the plane
* @param { Object } bbox - bounding box in which to display the plane
*/
rglwidgetClass.prototype.planeUpdateTriangles = function(id, bbox) {
var perms = [[0,0,1], [1,2,2], [2,1,0]],
x, xrow, elem, A, d, nhits, i, j, k, u, v, w, intersect, which, v0, v2, vx, reverse,
face1 = [], face2 = [], normals = [],
obj = this.getObj(id),
nPlanes = obj.normals.length;
obj.bbox = bbox;
obj.vertices = [];
obj.initialized = false;
for (elem = 0; elem < nPlanes; elem++) {
//    Vertex Av = normal.getRecycled(elem);
x = [];
A = obj.normals[elem];
d = obj.offsets[elem][0];
nhits = 0;
for (i=0; i<3; i++)
for (j=0; j<2; j++)
for (k=0; k<2; k++) {
u = perms[0][i];
v = perms[1][i];
w = perms[2][i];
if (A[w] !== 0.0) {
intersect = -(d + A[u]*bbox[j+2*u] + A[v]*bbox[k+2*v])/A[w];
if (bbox[2*w] < intersect && intersect < bbox[1+2*w]) {
xrow = [];
xrow[u] = bbox[j+2*u];
xrow[v] = bbox[k+2*v];
xrow[w] = intersect;
x.push(xrow);
face1[nhits] = j + 2*u;
face2[nhits] = k + 2*v;
nhits++;
}
}
}
if (nhits > 3) {
/* Re-order the intersections so the triangles work */
for (i=0; i<nhits-2; i++) {
which = 0; /* initialize to suppress warning */
for (j=i+1; j<nhits; j++) {
if (face1[i] == face1[j] || face1[i] == face2[j] ||
face2[i] == face1[j] || face2[i] == face2[j] ) {
which = j;
break;
}
}
if (which > i+1) {
this.swap(x, i+1, which);
this.swap(face1, i+1, which);
this.swap(face2, i+1, which);
}
}
}
if (nhits >= 3) {
/* Put in order so that the normal points out the FRONT of the faces */
v0 = [x[0][0] - x[1][0] , x[0][1] - x[1][1], x[0][2] - x[1][2]];
v2 = [x[2][0] - x[1][0] , x[2][1] - x[1][1], x[2][2] - x[1][2]];
/* cross-product */
vx = this.xprod(v0, v2);
reverse = this.dotprod(vx, A) > 0;
for (i=0; i<nhits-2; i++) {
obj.vertices.push(x[0]);
normals.push(A);
for (j=1; j<3; j++) {
obj.vertices.push(x[i + (reverse ? 3-j : j)]);
normals.push(A);
}
}
}
}
obj.pnormals = normals;
};
/**
* Initialize object for display
* @param { number } id - id of object to initialize
*/
rglwidgetClass.prototype.initObj = function(id) {
var obj = this.getObj(id),
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
is_lines = flags & this.f_is_lines,
fat_lines = flags & this.f_fat_lines,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
is_transparent = obj.is_transparent,
depth_sort = flags & this.f_depth_sort,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
fixed_size = flags & this.f_fixed_size,
is_twosided = (flags & this.f_is_twosided) > 0,
is_brush = flags & this.f_is_brush,
gl = this.gl || this.initGL(),
texinfo, drawtype, nclipplanes, f, nrows, oldrows,
i,j,v,v1,v2, mat, uri, matobj, pass, passes, pmode,
dim, nx, nz, attr;
if (typeof id !== "number") {
this.alertOnce("initObj id is "+typeof id);
}
obj.initialized = true;
if (type === "bboxdeco" || type === "subscene")
return;
if (type === "light") {
obj.ambient = new Float32Array(obj.colors[0].slice(0,3));
obj.diffuse = new Float32Array(obj.colors[1].slice(0,3));
obj.specular = new Float32Array(obj.colors[2].slice(0,3));
obj.lightDir = new Float32Array(obj.vertices[0]);
return;
}
if (type === "clipplanes") {
obj.vClipplane = this.flatten(this.cbind(obj.normals, obj.offsets));
return;
}
if (type === "background" && typeof obj.ids !== "undefined") {
obj.quad = this.flatten([].concat(obj.ids));
return;
}
if (is_transparent) {
depth_sort = ["triangles", "quads", "surface",
"spheres", "sprites", "text"].indexOf(type) >= 0;
}
if (is_brush)
this.initSelection(id);
if (typeof obj.vertices === "undefined")
obj.vertices = [];
v = obj.vertices;
obj.vertexCount = v.length;
if (!obj.vertexCount) return;
if (is_twosided) {
if (typeof obj.userAttributes === "undefined")
obj.userAttributes = {};
v1 = Array(v.length);
v2 = Array(v.length);
if (obj.type == "triangles" || obj.type == "quads") {
if (obj.type == "triangles")
nrow = 3;
else
nrow = 4;
for (i=0; i<Math.floor(v.length/nrow); i++)
for (j=0; j<nrow; j++) {
v1[nrow*i + j] = v[nrow*i + ((j+1) % nrow)];
v2[nrow*i + j] = v[nrow*i + ((j+2) % nrow)];
}
} else if (obj.type == "surface") {
dim = obj.dim[0];
nx = dim[0];
nz = dim[1];
for (j=0; j<nx; j++) {
for (i=0; i<nz; i++) {
if (i+1 < nz && j+1 < nx) {
v2[j + nx*i] = v[j + nx*(i+1)];
v1[j + nx*i] = v[j+1 + nx*(i+1)];
} else if (i+1 < nz) {
v2[j + nx*i] = v[j-1 + nx*i];
v1[j + nx*i] = v[j + nx*(i+1)];
} else {
v2[j + nx*i] = v[j + nx*(i-1)];
v1[j + nx*i] = v[j-1 + nx*(i-1)];
}
}
}
}
obj.userAttributes.aPos1 = v1;
obj.userAttributes.aPos2 = v2;
}
if (!sprites_3d) {
if (gl.isContextLost()) return;
obj.prog = gl.createProgram();
gl.attachShader(obj.prog, this.getShader( gl.VERTEX_SHADER,
this.getVertexShader(id) ));
gl.attachShader(obj.prog, this.getShader( gl.FRAGMENT_SHADER,
this.getFragmentShader(id) ));
//  Force aPos to location 0, aCol to location 1
gl.bindAttribLocation(obj.prog, 0, "aPos");
gl.bindAttribLocation(obj.prog, 1, "aCol");
gl.linkProgram(obj.prog);
var linked = gl.getProgramParameter(obj.prog, gl.LINK_STATUS);
if (!linked) {
// An error occurred while linking
var lastError = gl.getProgramInfoLog(obj.prog);
console.warn("Error in program linking:" + lastError);
gl.deleteProgram(obj.prog);
return;
}
}
if (type === "text") {
texinfo = this.drawTextToCanvas(obj.texts,
this.flatten(obj.cex),
this.flatten(obj.family),
this.flatten(obj.family));
}
if (fixed_quads && !sprites_3d) {
obj.ofsLoc = gl.getAttribLocation(obj.prog, "aOfs");
}
if (sprite_3d) {
obj.origLoc = gl.getUniformLocation(obj.prog, "uOrig");
obj.sizeLoc = gl.getUniformLocation(obj.prog, "uSize");
obj.usermatLoc = gl.getUniformLocation(obj.prog, "usermat");
}
if (has_texture || type == "text") {
if (!obj.texture)
obj.texture = gl.createTexture();
obj.texLoc = gl.getAttribLocation(obj.prog, "aTexcoord");
obj.sampler = gl.getUniformLocation(obj.prog, "uSampler");
}
if (has_texture) {
mat = obj.material;
if (typeof mat.uri !== "undefined")
uri = mat.uri;
else if (typeof mat.uriElementId === "undefined") {
matobj = this.getObj(mat.uriId);
if (typeof matobj !== "undefined") {
uri = matobj.material.uri;
} else {
uri = "";
}
} else
uri = document.getElementById(mat.uriElementId).rglinstance.getObj(mat.uriId).material.uri;
this.loadImageToTexture(uri, obj.texture);
}
if (type === "text") {
this.handleLoadedTexture(obj.texture, this.textureCanvas);
}
var stride = 3, nc, cofs, nofs, radofs, oofs, tofs, vnew, fnew,
nextofs = -1, pointofs = -1, alias, colors, key, selection, filter;
obj.alias = undefined;
colors = obj.colors;
j = this.scene.crosstalk.id.indexOf(id);
if (j >= 0) {
key = this.scene.crosstalk.key[j];
options = this.scene.crosstalk.options[j];
colors = colors.slice(0); 
for (i = 0; i < v.length; i++)
colors[i] = obj.colors[i % obj.colors.length].slice(0);
if ( (selection = this.scene.crosstalk.selection) &&
(selection.length || !options.selectedIgnoreNone) )
for (i = 0; i < v.length; i++) {
if (!selection.includes(key[i])) {
if (options.deselectedColor)
colors[i] = options.deselectedColor.slice(0);
colors[i][3] = colors[i][3]*options.deselectedFade;   /* default: mostly transparent if not selected */
} else if (options.selectedColor)
colors[i] = options.selectedColor.slice(0);
}
if ( (filter = this.scene.crosstalk.filter) )
for (i = 0; i < v.length; i++) 
if (!filter.includes(key[i])) {
if (options.filteredColor)
colors[i] = options.filteredColor.slice(0);
colors[i][3] = colors[i][3]*options.filteredFade;   /* default: completely hidden if filtered */
}
}  
nc = obj.colorCount = colors.length;
if (nc > 1) {
cofs = stride;
stride = stride + 4;
v = this.cbind(v, colors);
} else {
cofs = -1;
obj.onecolor = this.flatten(colors);
}
if (typeof obj.normals !== "undefined") {
nofs = stride;
stride = stride + 3;
v = this.cbind(v, typeof obj.pnormals !== "undefined" ? obj.pnormals : obj.normals);
} else
nofs = -1;
if (typeof obj.radii !== "undefined") {
radofs = stride;
stride = stride + 1;
// FIXME:  always concat the radii?
if (obj.radii.length === v.length) {
v = this.cbind(v, obj.radii);
} else if (obj.radii.length === 1) {
v = v.map(function(row, i, arr) { return row.concat(obj.radii[0]);});
}
} else
radofs = -1;
// Add default indices
f = Array(v.length);
for (i = 0; i < v.length; i++)
f[i] = i;
obj.f = [f,f];
if (type == "sprites" && !sprites_3d) {
tofs = stride;
stride += 2;
oofs = stride;
stride += 2;
vnew = new Array(4*v.length);
fnew = new Array(4*v.length);
alias = new Array(v.length);
var rescale = fixed_size ? 72 : 1,
size = obj.radii, s = rescale*size[0]/2;
last = v.length;
f = obj.f[0];
for (i=0; i < v.length; i++) {
if (size.length > 1)
s = rescale*size[i]/2;
vnew[i]  = v[i].concat([0,0,-s,-s]);
fnew[4*i] = f[i];
vnew[last]= v[i].concat([1,0, s,-s]);
fnew[4*i+1] = last++;
vnew[last]= v[i].concat([1,1, s, s]);
fnew[4*i+2] = last++;
vnew[last]= v[i].concat([0,1,-s, s]);
fnew[4*i+3] = last++;
alias[i] = [last-3, last-2, last-1];
}
v = vnew;
obj.vertexCount = v.length;
obj.f = [fnew, fnew];
} else if (type === "text") {
tofs = stride;
stride += 2;
oofs = stride;
stride += 2;
vnew = new Array(4*v.length);
f = obj.f[0];
fnew = new Array(4*f.length);
alias = new Array(v.length);
last = v.length;
for (i=0; i < v.length; i++) {
vnew[i]  = v[i].concat([0,-0.5]).concat(obj.adj[0]);
fnew[4*i] = f[i];
vnew[last] = v[i].concat([1,-0.5]).concat(obj.adj[0]);
fnew[4*i+1] = last++;
vnew[last] = v[i].concat([1, 1.5]).concat(obj.adj[0]);
fnew[4*i+2] = last++;
vnew[last] = v[i].concat([0, 1.5]).concat(obj.adj[0]);
fnew[4*i+3] = last++;
alias[i] = [last-3, last-2, last-1];
for (j=0; j < 4; j++) {
v1 = vnew[fnew[4*i+j]];
v1[tofs+2] = 2*(v1[tofs]-v1[tofs+2])*texinfo.widths[i];
v1[tofs+3] = 2*(v1[tofs+1]-v1[tofs+3])*texinfo.textHeights[i];
v1[tofs] = (texinfo.offsetsx[i] + v1[tofs]*texinfo.widths[i])/texinfo.canvasX;
v1[tofs+1] = 1.0-(texinfo.offsetsy[i] -
v1[tofs+1]*texinfo.textHeights[i])/texinfo.canvasY;
vnew[fnew[4*i+j]] = v1;
}
}
v = vnew;
obj.vertexCount = v.length;
obj.f = [fnew, fnew];
} else if (typeof obj.texcoords !== "undefined") {
tofs = stride;
stride += 2;
oofs = -1;
v = this.cbind(v, obj.texcoords);
} else {
tofs = -1;
oofs = -1;
}
obj.alias = alias;
if (typeof obj.userAttributes !== "undefined") {
obj.userAttribOffsets = {};
obj.userAttribLocations = {};
obj.userAttribSizes = {};
for (attr in obj.userAttributes) {
obj.userAttribLocations[attr] = gl.getAttribLocation(obj.prog, attr);
if (obj.userAttribLocations[attr] >= 0) { // Attribute may not have been used
obj.userAttribOffsets[attr] = stride;
v = this.cbind(v, obj.userAttributes[attr]);
stride = v[0].length;
obj.userAttribSizes[attr] = stride - obj.userAttribOffsets[attr];
}
}
}
if (typeof obj.userUniforms !== "undefined") {
obj.userUniformLocations = {};
for (attr in obj.userUniforms)
obj.userUniformLocations[attr] = gl.getUniformLocation(obj.prog, attr);
}
if (sprites_3d) {
obj.userMatrix = new CanvasMatrix4(obj.userMatrix);
obj.objects = this.flatten([].concat(obj.ids));
is_lit = false;
}
if (is_lit && !fixed_quads) {
obj.normLoc = gl.getAttribLocation(obj.prog, "aNorm");
}
nclipplanes = this.countClipplanes();
if (nclipplanes && !sprites_3d) {
obj.clipLoc = [];
for (i=0; i < nclipplanes; i++)
obj.clipLoc[i] = gl.getUniformLocation(obj.prog,"vClipplane" + i);
}
if (is_lit) {
obj.emissionLoc = gl.getUniformLocation(obj.prog, "emission");
obj.emission = new Float32Array(this.stringToRgb(this.getMaterial(id, "emission")));
obj.shininessLoc = gl.getUniformLocation(obj.prog, "shininess");
obj.shininess = this.getMaterial(id, "shininess");
obj.nlights = this.countLights();
obj.ambientLoc = [];
obj.ambient = new Float32Array(this.stringToRgb(this.getMaterial(id, "ambient")));
obj.specularLoc = [];
obj.specular = new Float32Array(this.stringToRgb(this.getMaterial(id, "specular")));
obj.diffuseLoc = [];
obj.lightDirLoc = [];
obj.viewpointLoc = [];
obj.finiteLoc = [];
for (i=0; i < obj.nlights; i++) {
obj.ambientLoc[i] = gl.getUniformLocation(obj.prog, "ambient" + i);
obj.specularLoc[i] = gl.getUniformLocation(obj.prog, "specular" + i);
obj.diffuseLoc[i] = gl.getUniformLocation(obj.prog, "diffuse" + i);
obj.lightDirLoc[i] = gl.getUniformLocation(obj.prog, "lightDir" + i);
obj.viewpointLoc[i] = gl.getUniformLocation(obj.prog, "viewpoint" + i);
obj.finiteLoc[i] = gl.getUniformLocation(obj.prog, "finite" + i);
}
}
obj.passes = is_twosided + 1;
obj.pmode = new Array(obj.passes);
for (pass = 0; pass < obj.passes; pass++) {
if (type === "triangles" || type === "quads" || type === "surface")
pmode = this.getMaterial(id, (pass === 0) ? "front" : "back");
else pmode = "filled";
obj.pmode[pass] = pmode;
}
obj.f.length = obj.passes;
for (pass = 0; pass < obj.passes; pass++) {
f = fnew = obj.f[pass];
pmode = obj.pmode[pass];
if (pmode === "culled")
f = [];
else if (pmode === "points") {
// stay with default
} else if ((type === "quads" || type === "text" ||
type === "sprites") && !sprites_3d) {
nrows = Math.floor(obj.vertexCount/4);
if (pmode === "filled") {
fnew = Array(6*nrows);
for (i=0; i < nrows; i++) {
fnew[6*i] = f[4*i];
fnew[6*i+1] = f[4*i + 1];
fnew[6*i+2] = f[4*i + 2];
fnew[6*i+3] = f[4*i];
fnew[6*i+4] = f[4*i + 2];
fnew[6*i+5] = f[4*i + 3];
}
} else {
fnew = Array(8*nrows);
for (i=0; i < nrows; i++) {
fnew[8*i] = f[4*i];
fnew[8*i+1] = f[4*i + 1];
fnew[8*i+2] = f[4*i + 1];
fnew[8*i+3] = f[4*i + 2];
fnew[8*i+4] = f[4*i + 2];
fnew[8*i+5] = f[4*i + 3];
fnew[8*i+6] = f[4*i + 3];
fnew[8*i+7] = f[4*i];
}
}
} else if (type === "triangles") {
nrows = Math.floor(obj.vertexCount/3);
if (pmode === "filled") {
fnew = Array(3*nrows);
for (i=0; i < fnew.length; i++) {
fnew[i] = f[i];
}
} else if (pmode === "lines") {
fnew = Array(6*nrows);
for (i=0; i < nrows; i++) {
fnew[6*i] = f[3*i];
fnew[6*i + 1] = f[3*i + 1];
fnew[6*i + 2] = f[3*i + 1];
fnew[6*i + 3] = f[3*i + 2];
fnew[6*i + 4] = f[3*i + 2];
fnew[6*i + 5] = f[3*i];
}
}
} else if (type === "spheres") {
// default
} else if (type === "surface") {
dim = obj.dim[0];
nx = dim[0];
nz = dim[1];
if (pmode === "filled") {
fnew = [];
for (j=0; j<nx-1; j++) {
for (i=0; i<nz-1; i++) {
fnew.push(f[j + nx*i],
f[j + nx*(i+1)],
f[j + 1 + nx*(i+1)],
f[j + nx*i],
f[j + 1 + nx*(i+1)],
f[j + 1 + nx*i]);
}
}
} else if (pmode === "lines") {
fnew = [];
for (j=0; j<nx; j++) {
for (i=0; i<nz; i++) {
if (i+1 < nz)
fnew.push(f[j + nx*i],
f[j + nx*(i+1)]);
if (j+1 < nx)
fnew.push(f[j + nx*i],
f[j+1 + nx*i]);
}
}
}
}
obj.f[pass] = fnew;
if (depth_sort) {
drawtype = "DYNAMIC_DRAW";
} else {
drawtype = "STATIC_DRAW";
}
}
if (fat_lines) {
alias = undefined;
obj.nextLoc = gl.getAttribLocation(obj.prog, "aNext");
obj.pointLoc = gl.getAttribLocation(obj.prog, "aPoint");
obj.aspectLoc = gl.getUniformLocation(obj.prog, "uAspect");
obj.lwdLoc = gl.getUniformLocation(obj.prog, "uLwd");
// Expand vertices to turn each segment into a pair of triangles
for (pass = 0; pass < obj.passes; pass++) {
f = obj.f[pass];	
oldrows = f.length;
if (obj.pmode[pass] === "lines") 
break;
}
if (type === "linestrip") 
nrows = 4*(oldrows - 1); 
else
nrows = 2*oldrows;
vnew = new Array(nrows);
fnew = new Array(1.5*nrows);
var fnext = new Array(nrows),
fpt = new Array(nrows), 
pt, start, gap = type === "linestrip" ? 3 : 1;
// We're going to turn each pair of vertices into 4 new ones, with the "next" and "pt" attributes
// added.
// We do this by copying the originals in the first pass, adding the new attributes, then in a 
// second pass add new vertices at the end.
for (i = 0; i < v.length; i++) {
vnew[i] = v[i].concat([0,0,0,0,0]); 
}
nextofs = stride;
pointofs = stride + 3;
stride = stride + 5;
// Now add the extras
last = v.length - 1;
ind = 0;
alias = new Array(f.length);
for (i = 0; i < f.length; i++)
alias[i] = [];
for (i = 0; i < f.length - 1; i++) {
if (type !== "linestrip" && i % 2 == 1)
continue;
k = ++last;
vnew[k] = vnew[f[i]].slice();
for (j=0; j<3; j++)
vnew[k][nextofs + j] = vnew[f[i+1]][j];
vnew[k][pointofs] = -1;
vnew[k][pointofs+1] = -1;
fnew[ind] = k;
last++;
vnew[last] = vnew[k].slice();
vnew[last][pointofs] = 1;
fnew[ind+1] = last;
alias[f[i]].push(last-1, last);
last++;
k = last;
vnew[k] = vnew[f[i+1]].slice();
for (j=0; j<3; j++)
vnew[k][nextofs + j] = vnew[f[i]][j];
vnew[k][pointofs] = -1;
vnew[k][pointofs+1] = 1;
fnew[ind+2] = k;
fnew[ind+3] = fnew[ind+1];
last++;
vnew[last] = vnew[k].slice();
vnew[last][pointofs] = 1;
fnew[ind+4] = last;
fnew[ind+5] = fnew[ind+2];
ind += 6;
alias[f[i+1]].push(last-1, last);
}
vnew.length = last+1;
v = vnew;
obj.vertexCount = v.length;
if (typeof alias !== "undefined" && typeof obj.alias !== "undefined") {  // Already have aliases from previous section?
var oldalias = obj.alias, newalias = Array(obj.alias.length);
for (i = 0; i < newalias.length; i++) {
newalias[i] = oldalias[i].slice();
for (j = 0; j < oldalias[i].length; j++)
Array.prototype.push.apply(newalias[i], alias[oldalias[j]]); // pushes each element 
}
obj.alias = newalias;
} else
obj.alias = alias;
for (pass = 0; pass < obj.passes; pass++)
if (type === "lines" || type === "linestrip" || obj.pmode[pass] == "lines") {
obj.f[pass] = fnew;
}
if (depth_sort) 
drawtype = "DYNAMIC_DRAW";
else
drawtype = "STATIC_DRAW";
}
for (pass = 0; pass < obj.passes; pass++) {
if (obj.vertexCount > 65535) {
if (this.index_uint) {
obj.f[pass] = new Uint32Array(obj.f[pass]);
obj.index_uint = true;
} else
this.alertOnce("Object has "+obj.vertexCount+" vertices, not supported in this browser.");
} else {
obj.f[pass] = new Uint16Array(obj.f[pass]);
obj.index_uint = false;
}
}
if (stride !== v[0].length) {
this.alertOnce("problem in stride calculation");
}
obj.vOffsets = {vofs:0, cofs:cofs, nofs:nofs, radofs:radofs, oofs:oofs, tofs:tofs,
nextofs:nextofs, pointofs:pointofs, stride:stride};
obj.values = new Float32Array(this.flatten(v));
if (type !== "spheres" && !sprites_3d) {
obj.buf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW); //
obj.ibuf = Array(obj.passes);
obj.ibuf[0] = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[0]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[0], gl[drawtype]);
if (is_twosided) {
obj.ibuf[1] = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[1]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[1], gl[drawtype]);
}
}
if (!sprites_3d) {
obj.mvMatLoc = gl.getUniformLocation(obj.prog, "mvMatrix");
obj.prMatLoc = gl.getUniformLocation(obj.prog, "prMatrix");
}
if (fixed_size) {
obj.textScaleLoc = gl.getUniformLocation(obj.prog, "textScale");
}
if (is_lit && !sprites_3d) {
obj.normMatLoc = gl.getUniformLocation(obj.prog, "normMatrix");
}
if (is_twosided) {
obj.frontLoc = gl.getUniformLocation(obj.prog, "front");
}
};
/**
* Set gl depth test based on object's material
* @param { number } id - object to use
*/
rglwidgetClass.prototype.setDepthTest = function(id) {
var gl = this.gl || this.initGL(),
tests = {never: gl.NEVER,
less:  gl.LESS,
equal: gl.EQUAL,
lequal:gl.LEQUAL,
greater: gl.GREATER,
notequal: gl.NOTEQUAL,
gequal: gl.GEQUAL,
always: gl.ALWAYS},
test = tests[this.getMaterial(id, "depth_test")];
gl.depthFunc(test);
};
rglwidgetClass.prototype.mode4type = {points : "POINTS",
linestrip : "LINE_STRIP",
abclines : "LINES",
lines : "LINES",
sprites : "TRIANGLES",
planes : "TRIANGLES",
text : "TRIANGLES",
quads : "TRIANGLES",
surface : "TRIANGLES",
triangles : "TRIANGLES"};
/**
* Sort objects from back to front
* @returns { number[] }
* @param { Object } obj - object to sort
*/
rglwidgetClass.prototype.depthSort = function(obj) {
var n = obj.centers.length,
depths = new Float32Array(n),
result = new Array(n),
compare = function(i,j) { return depths[j] - depths[i]; },
z, w;
for(i=0; i<n; i++) {
z = this.prmvMatrix.m13*obj.centers[i][0] +
this.prmvMatrix.m23*obj.centers[i][1] +
this.prmvMatrix.m33*obj.centers[i][2] +
this.prmvMatrix.m43;
w = this.prmvMatrix.m14*obj.centers[i][0] +
this.prmvMatrix.m24*obj.centers[i][1] +
this.prmvMatrix.m34*obj.centers[i][2] +
this.prmvMatrix.m44;
depths[i] = z/w;
result[i] = i;
}
result.sort(compare);
return result;
};
/**
* Draw an object in a subscene
* @param { number } id - object to draw
* @param { number } subsceneid - id of subscene
*/
rglwidgetClass.prototype.drawObj = function(id, subsceneid) {
var obj = this.getObj(id),
subscene = this.getObj(subsceneid),
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
is_transparent = flags & this.f_is_transparent,
depth_sort = flags & this.f_depth_sort,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
is_lines = flags & this.f_is_lines,
fat_lines = flags & this.f_fat_lines,
is_points = flags & this.f_is_points,
fixed_size = flags & this.f_fixed_size,
is_twosided = (flags & this.f_is_twosided) > 0,
gl = this.gl || this.initGL(),
mat,
sphereMV, baseofs, ofs, sscale, i, count, light,
pass, mode, pmode, attr;
if (typeof id !== "number") {
this.alertOnce("drawObj id is "+typeof id);
}
if (type === "planes") {
if (obj.bbox !== subscene.par3d.bbox || !obj.initialized) {
this.planeUpdateTriangles(id, subscene.par3d.bbox);
}
}
if (!obj.initialized)
this.initObj(id);
if (type === "clipplanes") {
count = obj.offsets.length;
var IMVClip = [];
for (i=0; i < count; i++) {
IMVClip[i] = this.multMV(this.invMatrix, obj.vClipplane.slice(4*i, 4*(i+1)));
}
obj.IMVClip = IMVClip;
return;
}
if (type === "light" || type === "bboxdeco" || !obj.vertexCount)
return;
if (!is_transparent &&
obj.someHidden) {
is_transparent = true;
depth_sort = ["triangles", "quads", "surface",
"spheres", "sprites", "text"].indexOf(type) >= 0;
}        
this.setDepthTest(id);
if (sprites_3d) {
var norigs = obj.vertices.length,
savenorm = new CanvasMatrix4(this.normMatrix);
this.origs = obj.vertices;
this.usermat = new Float32Array(obj.userMatrix.getAsArray());
this.radii = obj.radii;
this.normMatrix = subscene.spriteNormmat;
for (this.iOrig=0; this.iOrig < norigs; this.iOrig++) {
for (i=0; i < obj.objects.length; i++) {
this.drawObj(obj.objects[i], subsceneid);
}
}
this.normMatrix = savenorm;
return;
} else {
gl.useProgram(obj.prog);
}
if (sprite_3d) {
gl.uniform3fv(obj.origLoc, new Float32Array(this.origs[this.iOrig]));
if (this.radii.length > 1) {
gl.uniform1f(obj.sizeLoc, this.radii[this.iOrig][0]);
} else {
gl.uniform1f(obj.sizeLoc, this.radii[0][0]);
}
gl.uniformMatrix4fv(obj.usermatLoc, false, this.usermat);
}
if (type === "spheres") {
gl.bindBuffer(gl.ARRAY_BUFFER, this.sphere.buf);
} else {
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
}
gl.uniformMatrix4fv( obj.prMatLoc, false, new Float32Array(this.prMatrix.getAsArray()) );
gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(this.mvMatrix.getAsArray()) );
var clipcheck = 0,
clipplaneids = subscene.clipplanes,
clip, j;
for (i=0; i < clipplaneids.length; i++) {
clip = this.getObj(clipplaneids[i]);
for (j=0; j < clip.offsets.length; j++) {
gl.uniform4fv(obj.clipLoc[clipcheck + j], clip.IMVClip[j]);
}
clipcheck += clip.offsets.length;
}
if (typeof obj.clipLoc !== "undefined")
for (i=clipcheck; i < obj.clipLoc.length; i++)
gl.uniform4f(obj.clipLoc[i], 0,0,0,0);
if (is_lit) {
gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(this.normMatrix.getAsArray()) );
gl.uniform3fv( obj.emissionLoc, obj.emission);
gl.uniform1f( obj.shininessLoc, obj.shininess);
for (i=0; i < subscene.lights.length; i++) {
light = this.getObj(subscene.lights[i]);
gl.uniform3fv( obj.ambientLoc[i], this.componentProduct(light.ambient, obj.ambient));
gl.uniform3fv( obj.specularLoc[i], this.componentProduct(light.specular, obj.specular));
gl.uniform3fv( obj.diffuseLoc[i], light.diffuse);
gl.uniform3fv( obj.lightDirLoc[i], light.lightDir);
gl.uniform1i( obj.viewpointLoc[i], light.viewpoint);
gl.uniform1i( obj.finiteLoc[i], light.finite);
}
for (i=subscene.lights.length; i < obj.nlights; i++) {
gl.uniform3f( obj.ambientLoc[i], 0,0,0);
gl.uniform3f( obj.specularLoc[i], 0,0,0);
gl.uniform3f( obj.diffuseLoc[i], 0,0,0);
}
}
if (fixed_size) {
gl.uniform2f( obj.textScaleLoc, 0.75/this.vp.width, 0.75/this.vp.height);
}
gl.enableVertexAttribArray( this.posLoc );
var nc = obj.colorCount;
count = obj.vertexCount;
if (type === "spheres") {
subscene = this.getObj(subsceneid);
var scale = subscene.par3d.scale,
scount = count, indices;
gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
gl.enableVertexAttribArray(obj.normLoc );
gl.vertexAttribPointer(obj.normLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
gl.disableVertexAttribArray( this.colLoc );
var sphereNorm = new CanvasMatrix4();
sphereNorm.scale(scale[0], scale[1], scale[2]);
sphereNorm.multRight(this.normMatrix);
gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(sphereNorm.getAsArray()) );
if (nc == 1) {
gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
}
if (has_texture) {
gl.enableVertexAttribArray( obj.texLoc );
gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,
4*this.sphere.vOffsets.tofs);
gl.activeTexture(gl.TEXTURE0);
gl.bindTexture(gl.TEXTURE_2D, obj.texture);
gl.uniform1i( obj.sampler, 0);
}
if (depth_sort)
indices = this.depthSort(obj);
for (i = 0; i < scount; i++) {
sphereMV = new CanvasMatrix4();
if (depth_sort) {
baseofs = indices[i]*obj.vOffsets.stride;
} else {
baseofs = i*obj.vOffsets.stride;
}
ofs = baseofs + obj.vOffsets.radofs;
sscale = obj.values[ofs];
sphereMV.scale(sscale/scale[0], sscale/scale[1], sscale/scale[2]);
sphereMV.translate(obj.values[baseofs],
obj.values[baseofs+1],
obj.values[baseofs+2]);
sphereMV.multRight(this.mvMatrix);
gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(sphereMV.getAsArray()) );
if (nc > 1) {
ofs = baseofs + obj.vOffsets.cofs;
gl.vertexAttrib4f( this.colLoc, obj.values[ofs],
obj.values[ofs+1],
obj.values[ofs+2],
obj.values[ofs+3] );
}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.sphere.ibuf);
gl.drawElements(gl.TRIANGLES, this.sphere.sphereCount, gl.UNSIGNED_SHORT, 0);
}
return;
} else {
if (obj.colorCount === 1) {
gl.disableVertexAttribArray( this.colLoc );
gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
} else {
gl.enableVertexAttribArray( this.colLoc );
gl.vertexAttribPointer(this.colLoc, 4, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.cofs);
}
}
if (is_lit && obj.vOffsets.nofs > 0) {
gl.enableVertexAttribArray( obj.normLoc );
gl.vertexAttribPointer(obj.normLoc, 3, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.nofs);
}
if (has_texture || type === "text") {
gl.enableVertexAttribArray( obj.texLoc );
gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.tofs);
gl.activeTexture(gl.TEXTURE0);
gl.bindTexture(gl.TEXTURE_2D, obj.texture);
gl.uniform1i( obj.sampler, 0);
}
if (fixed_quads) {
gl.enableVertexAttribArray( obj.ofsLoc );
gl.vertexAttribPointer(obj.ofsLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.oofs);
}
if (typeof obj.userAttributes !== "undefined") {
for (attr in obj.userAttribSizes) {  // Not all attributes may have been used
gl.enableVertexAttribArray( obj.userAttribLocations[attr] );
gl.vertexAttribPointer( obj.userAttribLocations[attr], obj.userAttribSizes[attr],
gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.userAttribOffsets[attr]);
}
}
if (typeof obj.userUniforms !== "undefined") {
for (attr in obj.userUniformLocations) {
var loc = obj.userUniformLocations[attr];
if (loc !== null) {
var uniform = obj.userUniforms[attr];
if (typeof uniform.length === "undefined")
gl.uniform1f(loc, uniform);
else if (typeof uniform[0].length === "undefined") {
uniform = new Float32Array(uniform);
switch(uniform.length) {
case 2: gl.uniform2fv(loc, uniform); break;
case 3: gl.uniform3fv(loc, uniform); break;
case 4: gl.uniform4fv(loc, uniform); break;
default: console.warn("bad uniform length");
}
} else if (uniform.length == 4 && uniform[0].length == 4)
gl.uniformMatrix4fv(loc, false, new Float32Array(uniform.getAsArray()));
else
console.warn("unsupported uniform matrix");
}
}
}
for (pass = 0; pass < obj.passes; pass++) {
pmode = obj.pmode[pass];
if (pmode === "culled")
continue;
mode = fat_lines && (is_lines || pmode == "lines") ? "TRIANGLES" : this.mode4type[type];
if (depth_sort && pmode == "filled") {// Don't try depthsorting on wireframe or points
var faces = this.depthSort(obj),
nfaces = faces.length,
frowsize = Math.floor(obj.f[pass].length/nfaces);
if (type !== "spheres") {
var f = obj.index_uint ? new Uint32Array(obj.f[pass].length) : new Uint16Array(obj.f[pass].length);
for (i=0; i<nfaces; i++) {
for (j=0; j<frowsize; j++) {
f[frowsize*i + j] = obj.f[pass][frowsize*faces[i] + j];
}
}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, f, gl.DYNAMIC_DRAW);
}
}
if (is_twosided)
gl.uniform1i(obj.frontLoc, pass !== 0);
if (type !== "spheres") 
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
if (type === "sprites" || type === "text" || type === "quads") {
count = count * 6/4;
} else if (type === "surface") {
count = obj.f[pass].length;
}
count = obj.f[pass].length;
if (!is_lines && pmode === "lines" && !fat_lines) {
mode = "LINES";
} else if (pmode === "points") {
mode = "POINTS";
}
if ((is_lines || pmode === "lines") && fat_lines) {
gl.enableVertexAttribArray(obj.pointLoc);
gl.vertexAttribPointer(obj.pointLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.pointofs);
gl.enableVertexAttribArray(obj.nextLoc );
gl.vertexAttribPointer(obj.nextLoc, 3, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.nextofs);
gl.uniform1f(obj.aspectLoc, this.vp.width/this.vp.height);
gl.uniform1f(obj.lwdLoc, this.getMaterial(id, "lwd")/this.vp.height);
}
gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*obj.vOffsets.stride,  4*obj.vOffsets.vofs);
gl.drawElements(gl[mode], count, obj.index_uint ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT, 0);
}
};
/**
* Draw the background for a subscene
* @param { number } id - id of background object
* @param { number } subsceneid - id of subscene
*/
rglwidgetClass.prototype.drawBackground = function(id, subsceneid) {
var gl = this.gl || this.initGL(),
obj = this.getObj(id),
bg, i;
if (!obj.initialized)
this.initObj(id);
if (obj.colors.length) {
bg = obj.colors[0];
gl.clearColor(bg[0], bg[1], bg[2], bg[3]);
gl.depthMask(true);
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
}
if (typeof obj.quad !== "undefined") {
this.prMatrix.makeIdentity();
this.mvMatrix.makeIdentity();
gl.disable(gl.BLEND);
gl.disable(gl.DEPTH_TEST);
gl.depthMask(false);
for (i=0; i < obj.quad.length; i++)
this.drawObj(obj.quad[i], subsceneid);
}
};
/**
* Draw a subscene
* @param { number } subsceneid - id of subscene
* @param { boolean } opaquePass - is this the opaque drawing pass?
*/
rglwidgetClass.prototype.drawSubscene = function(subsceneid, opaquePass) {
var gl = this.gl || this.initGL(),
sub = this.getObj(subsceneid),
objects = this.scene.objects,
subids = sub.objects,
subscene_has_faces = false,
subscene_needs_sorting = false,
flags, i, obj;
if (sub.par3d.skipRedraw)
return;
for (i=0; i < subids.length; i++) {
obj = objects[subids[i]];
flags = obj.flags;
if (typeof flags !== "undefined") {
subscene_has_faces |= (flags & this.f_is_lit)
& !(flags & this.f_fixed_quads);
obj.is_transparent = (flags & this.f_is_transparent) || obj.someHidden;
subscene_needs_sorting |= (flags & this.f_depth_sort) || obj.is_transparent;
}
}
this.setViewport(subsceneid);
if (typeof sub.backgroundId !== "undefined" && opaquePass)
this.drawBackground(sub.backgroundId, subsceneid);
if (subids.length) {
this.setprMatrix(subsceneid);
this.setmvMatrix(subsceneid);
if (subscene_has_faces) {
this.setnormMatrix(subsceneid);
if ((sub.flags & this.f_sprites_3d) &&
typeof sub.spriteNormmat === "undefined") {
sub.spriteNormmat = new CanvasMatrix4(this.normMatrix);
}
}
if (subscene_needs_sorting)
this.setprmvMatrix();
var clipids = sub.clipplanes;
if (typeof clipids === "undefined") {
console.warn("bad clipids");
}
if (clipids.length > 0) {
this.invMatrix = new CanvasMatrix4(this.mvMatrix);
this.invMatrix.invert();
for (i = 0; i < clipids.length; i++)
this.drawObj(clipids[i], subsceneid);
}
subids = sub.opaque.concat(sub.transparent);
if (opaquePass) {
gl.enable(gl.DEPTH_TEST);
gl.depthMask(true);
gl.disable(gl.BLEND);
for (i = 0; i < subids.length; i++) {
if (!this.getObj(subids[i]).is_transparent)	
this.drawObj(subids[i], subsceneid);
}
} else {
gl.depthMask(false);
gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA,
gl.ONE, gl.ONE);
gl.enable(gl.BLEND);
for (i = 0; i < subids.length; i++) {
if (this.getObj(subids[i]).is_transparent)
this.drawObj(subids[i], subsceneid);
}
}
subids = sub.subscenes;
for (i = 0; i < subids.length; i++) {
this.drawSubscene(subids[i], opaquePass);
}
}
};
/**
* Respond to brush change
*/
rglwidgetClass.prototype.selectionChanged = function() {
var i, j, k, id, subid = this.select.subscene, subscene,
objids, obj,
p1 = this.select.region.p1, p2 = this.select.region.p2,
filter, selection = [], handle, keys, xmin, x, xmax, ymin, y, ymax, z, v,
someHidden;
if (!subid)
return;
subscene = this.getObj(subid);
objids = subscene.objects;
filter = this.scene.crosstalk.filter;
this.setmvMatrix(subid);
this.setprMatrix(subid);
this.setprmvMatrix();
xmin = Math.min(p1.x, p2.x);
xmax = Math.max(p1.x, p2.x);
ymin = Math.min(p1.y, p2.y);
ymax = Math.max(p1.y, p2.y);
for (i = 0; i < objids.length; i++) {
id = objids[i];
j = this.scene.crosstalk.id.indexOf(id);
if (j >= 0) {
keys = this.scene.crosstalk.key[j];
obj = this.getObj(id);
someHidden = false;
for (k = 0; k < keys.length; k++) {
if (filter && filter.indexOf(keys[k]) < 0) {
someHidden = true;
continue;
}
v = [].concat(obj.vertices[k]).concat(1.0);
v = this.multVM(v, this.prmvMatrix);
x = v[0]/v[3];
y = v[1]/v[3];
z = v[2]/v[3];
if (xmin <= x && x <= xmax && ymin <= y && y <= ymax && -1.0 <= z && z <= 1.0) {
selection.push(keys[k]);
} else
someHidden = true;
}
obj.someHidden = someHidden && (filter || selection.length);
obj.initialized = false;
/* Who should we notify?  Only shared data in the current subscene, or everyone? */
if (!this.equalArrays(selection, this.scene.crosstalk.selection)) {
handle = this.scene.crosstalk.sel_handle[j];
handle.set(selection, {rglSubsceneId: this.select.subscene});
}
}
}
};
/**
* Respond to selection or filter change from crosstalk
* @param { Object } event - crosstalk event
* @param { boolean } filter - filter or selection?
*/
rglwidgetClass.prototype.selection = function(event, filter) {
var i, j, ids, obj, keys, crosstalk = this.scene.crosstalk,
selection, someHidden;
// Record the message and find out if this event makes some objects have mixed values:
crosstalk = this.scene.crosstalk;
if (filter) {
filter = crosstalk.filter = event.value;
selection = crosstalk.selection;
} else {  
selection = crosstalk.selection = event.value;
filter = crosstalk.filter;
}
ids = crosstalk.id;
for (i = 0; i < ids.length ; i++) {
obj = this.getObj(ids[i]);
obj.initialized = false;
keys = crosstalk.key[i];
someHidden = false;
for (j = 0; j < keys.length && !someHidden; j++) {
if ((filter && filter.indexOf(keys[j]) < 0) ||
(selection.length && selection.indexOf(keys[j]) < 0))
someHidden = true;
}
obj.someHidden = someHidden;
}
this.drawScene();
};
/**
* Clear the selection brush
* @param { number } except - Subscene that should ignore this request
*/
rglwidgetClass.prototype.clearBrush = function(except) {
if (this.select.subscene != except) {
this.select.state = "inactive";
this.delFromSubscene(this.scene.brushId, this.select.subscene);
}
this.drawScene();
};
/**
* Compute mouse coordinates relative to current canvas
* @returns { Object }
* @param { Object } event - event object from mouse click
*/
rglwidgetClass.prototype.relMouseCoords = function(event) {
var totalOffsetX = 0,
totalOffsetY = 0,
currentElement = this.canvas;
do {
totalOffsetX += currentElement.offsetLeft;
totalOffsetY += currentElement.offsetTop;
currentElement = currentElement.offsetParent;
}
while(currentElement);
var canvasX = event.pageX - totalOffsetX,
canvasY = event.pageY - totalOffsetY;
return {x:canvasX, y:canvasY};
};
/**
* Set mouse handlers for the scene
*/
rglwidgetClass.prototype.setMouseHandlers = function() {
var self = this, activeSubscene, handler,
handlers = {}, drag = 0;
handlers.rotBase = 0;
this.screenToVector = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height,
radius = Math.max(width, height)/2.0,
cx = width/2.0,
cy = height/2.0,
px = (x-cx)/radius,
py = (y-cy)/radius,
plen = Math.sqrt(px*px+py*py);
if (plen > 1.e-6) {
px = px/plen;
py = py/plen;
}
var angle = (Math.SQRT2 - plen)/Math.SQRT2*Math.PI/2,
z = Math.sin(angle),
zlen = Math.sqrt(1.0 - z*z);
px = px * zlen;
py = py * zlen;
return [px, py, z];
};
handlers.trackballdown = function(x,y) {
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
handlers.rotBase = this.screenToVector(x, y);
this.saveMat = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
}
};
handlers.trackballmove = function(x,y) {
var rotCurrent = this.screenToVector(x,y),
rotBase = handlers.rotBase,
dot = rotBase[0]*rotCurrent[0] +
rotBase[1]*rotCurrent[1] +
rotBase[2]*rotCurrent[2],
angle = Math.acos( dot/this.vlen(rotBase)/this.vlen(rotCurrent) )*180.0/Math.PI,
axis = this.xprod(rotBase, rotCurrent),
objects = this.scene.objects,
activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
l = activeModel.par3d.listeners,
i;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.userMatrix.load(objects[l[i]].saveMat);
activeSub.par3d.userMatrix.rotate(angle, axis[0], axis[1], axis[2]);
}
this.drawScene();
};
handlers.trackballend = 0;
this.clamp = function(x, lo, hi) {
return Math.max(lo, Math.min(x, hi));
};
this.screenToPolar = function(x,y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height,
r = Math.min(width, height)/2,
dx = this.clamp(x - width/2, -r, r),
dy = this.clamp(y - height/2, -r, r);
return [Math.asin(dx/r), Math.asin(-dy/r)];
};
handlers.polardown = function(x,y) {
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
handlers.dragBase = this.screenToPolar(x, y);
this.saveMat = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
activeSub.camBase = [-Math.atan2(activeSub.saveMat.m13, activeSub.saveMat.m11),
Math.atan2(activeSub.saveMat.m32, activeSub.saveMat.m22)];
}
};
handlers.polarmove = function(x,y) {
var dragCurrent = this.screenToPolar(x,y),
activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
objects = this.scene.objects,
l = activeModel.par3d.listeners,
i, changepos = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
for (j=0; j<2; j++)
changepos[j] = -(dragCurrent[j] - handlers.dragBase[j]);
activeSub.par3d.userMatrix.makeIdentity();
activeSub.par3d.userMatrix.rotate(changepos[0]*180/Math.PI, 0,-1,0);
activeSub.par3d.userMatrix.multRight(objects[l[i]].saveMat);
activeSub.par3d.userMatrix.rotate(changepos[1]*180/Math.PI, -1,0,0);
}
this.drawScene();
};
handlers.polarend = 0;
handlers.axisdown = function(x,y) {
handlers.rotBase = this.screenToVector(x, this.canvas.height/2);
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
}
};
handlers.axismove = function(x,y) {
var rotCurrent = this.screenToVector(x, this.canvas.height/2),
rotBase = handlers.rotBase,
angle = (rotCurrent[0] - rotBase[0])*180/Math.PI,
rotMat = new CanvasMatrix4();
rotMat.rotate(angle, handlers.axis[0], handlers.axis[1], handlers.axis[2]);
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.userMatrix.load(activeSub.saveMat);
activeSub.par3d.userMatrix.multLeft(rotMat);
}
this.drawScene();
};
handlers.axisend = 0;
handlers.y0zoom = 0;
handlers.zoom0 = 0;
handlers.zoomdown = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
handlers.y0zoom = y;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.zoom0 = Math.log(activeSub.par3d.zoom);
}
};
handlers.zoommove = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.zoom = Math.exp(activeSub.zoom0 + (y-handlers.y0zoom)/this.canvas.height);
}
this.drawScene();
};
handlers.zoomend = 0;
handlers.y0fov = 0;
handlers.fovdown = function(x, y) {
handlers.y0fov = y;
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.fov0 = activeSub.par3d.FOV;
}
};
handlers.fovmove = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.FOV = Math.max(1, Math.min(179, activeSub.fov0 +
180*(y-handlers.y0fov)/this.canvas.height));
}
this.drawScene();
};
handlers.fovend = 0;
handlers.selectingdown = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height, 
p = {x: 2.0*x/width - 1.0, y: 2.0*y/height - 1.0};
this.select.region = {p1: p, p2: p};
if (this.select.subscene && this.select.subscene != activeSubscene)
this.delFromSubscene(this.scene.brushId, this.select.subscene);
this.select.subscene = activeSubscene;
this.addToSubscene(this.scene.brushId, activeSubscene);
this.select.state = "changing";
if (typeof this.scene.brushId !== "undefined")
this.getObj(this.scene.brushId).initialized = false;
this.drawScene();
};
handlers.selectingmove = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height;
if (this.select.state === "inactive") 
return;
this.select.region.p2 = {x: 2.0*x/width - 1.0, y: 2.0*y/height - 1.0};
if (typeof this.scene.brushId !== "undefined")
this.getObj(this.scene.brushId).initialized = false;
this.drawScene();
};
handlers.selectingend = 0;
this.canvas.onmousedown = function ( ev ){
if (!ev.which) // Use w3c defns in preference to MS
switch (ev.button) {
case 0: ev.which = 1; break;
case 1:
case 4: ev.which = 2; break;
case 2: ev.which = 3;
}
drag = ["left", "middle", "right"][ev.which-1];
var coords = self.relMouseCoords(ev);
coords.y = self.canvas.height-coords.y;
activeSubscene = self.whichSubscene(coords);
var sub = self.getObj(activeSubscene), f;
handler = sub.par3d.mouseMode[drag];
switch (handler) {
case "xAxis":
handler = "axis";
handlers.axis = [1.0, 0.0, 0.0];
break;
case "yAxis":
handler = "axis";
handlers.axis = [0.0, 1.0, 0.0];
break;
case "zAxis":
handler = "axis";
handlers.axis = [0.0, 0.0, 1.0];
break;
}
f = handlers[handler + "down"];
if (f) {
coords = self.translateCoords(activeSubscene, coords);
f.call(self, coords.x, coords.y);
ev.preventDefault();
} else
console.warn("Mouse handler '" + handler + "' is not implemented.");
};
this.canvas.onmouseup = function ( ev ){
if ( drag === 0 ) return;
var f = handlers[handler + "end"];
if (f) {
f.call(self);
ev.preventDefault();
}
drag = 0;
};
this.canvas.onmouseout = this.canvas.onmouseup;
this.canvas.onmousemove = function ( ev ) {
if ( drag === 0 ) return;
var f = handlers[handler + "move"];
if (f) {
var coords = self.relMouseCoords(ev);
coords.y = self.canvas.height - coords.y;
coords = self.translateCoords(activeSubscene, coords);
f.call(self, coords.x, coords.y);
}
};
handlers.wheelHandler = function(ev) {
var del = 1.02, i;
if (ev.shiftKey) del = 1.002;
var ds = ((ev.detail || ev.wheelDelta) > 0) ? del : (1 / del);
if (typeof activeSubscene === "undefined")
activeSubscene = self.scene.rootSubscene;
var activeSub = self.getObj(activeSubscene),
activeProjection = self.getObj(self.useid(activeSub.id, "projection")),
l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = self.getObj(l[i]);
activeSub.par3d.zoom *= ds;
}
self.drawScene();
ev.preventDefault();
};
this.canvas.addEventListener("DOMMouseScroll", handlers.wheelHandler, false);
this.canvas.addEventListener("mousewheel", handlers.wheelHandler, false);
};
/**
* Find a particular subscene by inheritance
* @returns { number } id of subscene to use
* @param { number } subsceneid - child subscene
* @param { string } type - type of inheritance:  "projection" or "model"
*/
rglwidgetClass.prototype.useid = function(subsceneid, type) {
var sub = this.getObj(subsceneid);
if (sub.embeddings[type] === "inherit")
return(this.useid(sub.parent, type));
else
return subsceneid;
};
/**
* Check whether point is in viewport of subscene
* @returns {boolean}
* @param { Object } coords - screen coordinates of point
* @param { number } subsceneid - subscene to check
*/
rglwidgetClass.prototype.inViewport = function(coords, subsceneid) {
var viewport = this.getObj(subsceneid).par3d.viewport,
x0 = coords.x - viewport.x*this.canvas.width,
y0 = coords.y - viewport.y*this.canvas.height;
return 0 <= x0 && x0 <= viewport.width*this.canvas.width &&
0 <= y0 && y0 <= viewport.height*this.canvas.height;
};
/**
* Find which subscene contains a point
* @returns { number } subscene id
* @param { Object } coords - coordinates of point
*/
rglwidgetClass.prototype.whichSubscene = function(coords) {
var self = this,
recurse = function(subsceneid) {
var subscenes = self.getChildSubscenes(subsceneid), i, id;
for (i=0; i < subscenes.length; i++) {
id = recurse(subscenes[i]);
if (typeof(id) !== "undefined")
return(id);
}
if (self.inViewport(coords, subsceneid))
return(subsceneid);
else
return undefined;
},
rootid = this.scene.rootSubscene,
result = recurse(rootid);
if (typeof(result) === "undefined")
result = rootid;
return result;
};
/**
* Translate from window coordinates to viewport coordinates
* @returns { Object } translated coordinates
* @param { number } subsceneid - which subscene to use?
* @param { Object } coords - point to translate
*/
rglwidgetClass.prototype.translateCoords = function(subsceneid, coords) {
var viewport = this.getObj(subsceneid).par3d.viewport;
return {x: coords.x - viewport.x*this.canvas.width,
y: coords.y - viewport.y*this.canvas.height};
};
/**
* Initialize the sphere object
*/
rglwidgetClass.prototype.initSphere = function() {
var verts = this.scene.sphereVerts,
reuse = verts.reuse, result;
if (typeof reuse !== "undefined") {
var prev = document.getElementById(reuse).rglinstance.sphere;
result = {values: prev.values, vOffsets: prev.vOffsets, it: prev.it};
} else
result = {values: new Float32Array(this.flatten(this.cbind(this.transpose(verts.vb),
this.transpose(verts.texcoords)))),
it: new Uint16Array(this.flatten(this.transpose(verts.it))),
vOffsets: {vofs:0, cofs:-1, nofs:-1, radofs:-1, oofs:-1,
tofs:3, nextofs:-1, pointofs:-1, stride:5}};
result.sphereCount = result.it.length;
this.sphere = result;
};
/**
* Set the vertices in the selection box object
*/
rglwidgetClass.prototype.initSelection = function(id) {
if (typeof this.select.region === "undefined")
return;
var obj = this.getObj(id),
width = this.canvas.width,
height = this.canvas.height, 
p1 = this.select.region.p1,
p2 = this.select.region.p2;
obj.vertices = [[p1.x, p1.y, 0.0],
[p2.x, p1.y, 0.0],
[p2.x, p2.y, 0.0],
[p1.x, p2.y, 0.0],
[p1.x, p1.y, 0.0]];
};
/**
* Do the gl part of initializing the sphere
*/
rglwidgetClass.prototype.initSphereGL = function() {
var gl = this.gl || this.initGL(), sphere = this.sphere;
if (gl.isContextLost()) return;
sphere.buf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, sphere.buf);
gl.bufferData(gl.ARRAY_BUFFER, sphere.values, gl.STATIC_DRAW);
sphere.ibuf = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sphere.ibuf);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, sphere.it, gl.STATIC_DRAW);
return;
};
/**
* Initialize the DOM object
* @param { Object } el - the DOM object
* @param { Object } x - the scene data sent by JSON from R
*/
rglwidgetClass.prototype.initialize = function(el, x) {
this.textureCanvas = document.createElement("canvas");
this.textureCanvas.style.display = "block";
this.scene = x;
this.normMatrix = new CanvasMatrix4();
this.saveMat = {};
this.distance = null;
this.posLoc = 0;
this.colLoc = 1;
if (el) {
el.rglinstance = this;
this.el = el;
this.webGLoptions = el.rglinstance.scene.webGLoptions;
this.initCanvas();
}
};
/**
* Restart the WebGL canvas
*/
rglwidgetClass.prototype.restartCanvas = function() {
var newcanvas = document.createElement("canvas"),
self = this;
newcanvas.width = this.el.width;
newcanvas.height = this.el.height;
newcanvas.addEventListener("webglcontextrestored",
this.onContextRestored, false);
newcanvas.addEventListener("webglcontextlost",
this.onContextLost, false);
while (this.el.firstChild) {
this.el.removeChild(this.el.firstChild);
}
this.el.appendChild(newcanvas);
this.canvas = newcanvas;
this.setMouseHandlers();
if (this.gl) 
Object.keys(this.scene.objects).forEach(function(key){
self.getObj(parseInt(key, 10)).texture = undefined; 
});
this.gl = null;
};
/**
* Initialize the WebGL canvas
*/
rglwidgetClass.prototype.initCanvas = function() {
this.restartCanvas();
var objs = this.scene.objects,
self = this;
Object.keys(objs).forEach(function(key){
var id = parseInt(key, 10),
obj = self.getObj(id);
if (typeof obj.reuse !== "undefined")
self.copyObj(id, obj.reuse);
});
Object.keys(objs).forEach(function(key){
self.initSubscene(parseInt(key, 10));
});
this.setMouseHandlers();
this.initSphere();
this.onContextRestored = function(event) {
self.initGL();
self.drawScene();
};
this.onContextLost = function(event) {
if (!self.drawing)
this.gl = null;
event.preventDefault();
};
this.initGL0();
this.lazyLoadScene = function() {
if (typeof self.slide === "undefined")
self.slide = self.getSlide();
if (self.isInBrowserViewport()) {
if (!self.gl || self.gl.isContextLost())
self.initGL();
self.drawScene();
}
};
window.addEventListener("DOMContentLoaded", this.lazyLoadScene, false);
window.addEventListener("load", this.lazyLoadScene, false);
window.addEventListener("resize", this.lazyLoadScene, false);
window.addEventListener("scroll", this.lazyLoadScene, false);
this.slide = this.getSlide();
if (this.slide) {
if (typeof this.slide.rgl === "undefined")
this.slide.rgl = [this];
else
this.slide.rgl.push(this);
if (this.scene.context.rmarkdown) 
if (this.scene.context.rmarkdown === "ioslides_presentation") {
this.slide.setAttribute("slideenter", "this.rgl.forEach(function(scene) { scene.lazyLoadScene.call(window);})");
} else if (this.scene.context.rmarkdown === "slidy_presentation") {
// This method would also work in ioslides, but it gets triggered
// something like 5 times per slide for every slide change, so
// you'd need a quicker function than lazyLoadScene.
var MutationObserver = window.MutationObserver || window.WebKitMutationObserver || window.MozMutationObserver,
observer = new MutationObserver(function(mutations) {
mutations.forEach(function(mutation) {
self.slide.rgl.forEach(function(scene) { scene.lazyLoadScene.call(window); });});});
observer.observe(this.slide, { attributes: true, attributeFilter:["class"] });
}
}
};
/**
* Start the writeWebGL scene. This is only used by writeWebGL; rglwidget has
no debug element and does the drawing in rglwidget.js.
*/
rglwidgetClass.prototype.start = function() {
if (typeof this.prefix !== "undefined") {
this.debugelement = document.getElementById(this.prefix + "debug");
this.debug("");
}
this.drag = 0;
this.drawScene();
};
/**
* Display a debug message
* @param { string } msg - The message to display
* @param { Object } [img] - Image to insert before message
*/
rglwidgetClass.prototype.debug = function(msg, img) {
if (typeof this.debugelement !== "undefined" && this.debugelement !== null) {
this.debugelement.innerHTML = msg;
if (typeof img !== "undefined") {
this.debugelement.insertBefore(img, this.debugelement.firstChild);
}
} else if (msg !== "")
alert(msg);
};
/**
* Get the snapshot image of this scene
* @returns { Object } The img DOM element
*/
rglwidgetClass.prototype.getSnapshot = function() {
var img;
if (typeof this.scene.snapshot !== "undefined") {
img = document.createElement("img");
img.src = this.scene.snapshot;
img.alt = "Snapshot";
}
return img;
};
/**
* Initial test for WebGL
*/
rglwidgetClass.prototype.initGL0 = function() {
if (!window.WebGLRenderingContext){
alert("Your browser does not support WebGL. See http://get.webgl.org");
return;
}
};
/**
* If we are in an ioslides or slidy presentation, get the
* DOM element of the current slide
* @returns { Object }
*/
rglwidgetClass.prototype.getSlide = function() {
var result = this.el, done = false;
while (result && !done && this.scene.context.rmarkdown) {
switch(this.scene.context.rmarkdown) {
case "ioslides_presentation":
if (result.tagName === "SLIDE") return result;
break;
case "slidy_presentation":
if (result.tagName === "DIV" && result.classList.contains("slide"))
return result;
break;
default: return null;
}
result = result.parentElement;
}
return null;
};
/**
* Is this scene visible in the browser?
* @returns { boolean }
*/
rglwidgetClass.prototype.isInBrowserViewport = function() {
var rect = this.canvas.getBoundingClientRect(),
windHeight = (window.innerHeight || document.documentElement.clientHeight),
windWidth = (window.innerWidth || document.documentElement.clientWidth);
if (this.scene.context && this.scene.context.rmarkdown !== null) {
if (this.slide)
return (this.scene.context.rmarkdown === "ioslides_presentation" &&
this.slide.classList.contains("current")) ||
(this.scene.context.rmarkdown === "slidy_presentation" &&
!this.slide.classList.contains("hidden"));
}
return (
rect.top >= -windHeight &&
rect.left >= -windWidth &&
rect.bottom <= 2*windHeight &&
rect.right <= 2*windWidth);
};
/**
* Initialize WebGL
* @returns { Object } the WebGL context
*/
rglwidgetClass.prototype.initGL = function() {
var self = this;
if (this.gl) {
if (!this.drawing && this.gl.isContextLost())
this.restartCanvas();
else
return this.gl;
}
// if (!this.isInBrowserViewport()) return; Return what??? At this point we know this.gl is null.
this.canvas.addEventListener("webglcontextrestored",
this.onContextRestored, false);
this.canvas.addEventListener("webglcontextlost",
this.onContextLost, false);
this.gl = this.canvas.getContext("webgl", this.webGLoptions) ||
this.canvas.getContext("experimental-webgl", this.webGLoptions);
this.index_uint = this.gl.getExtension("OES_element_index_uint");
var save = this.startDrawing();
this.initSphereGL();
Object.keys(this.scene.objects).forEach(function(key){
self.initObj(parseInt(key, 10));
});
this.stopDrawing(save);
return this.gl;
};
/**
* Resize the display to match element
* @param { Object } el - DOM element to match
*/
rglwidgetClass.prototype.resize = function(el) {
this.canvas.width = el.width;
this.canvas.height = el.height;
};
/**
* Draw the whole scene
*/
rglwidgetClass.prototype.drawScene = function() {
var gl = this.gl || this.initGL(),
wasDrawing = this.startDrawing();
if (!wasDrawing) {
if (this.select.state !== "inactive")
this.selectionChanged();
gl.enable(gl.DEPTH_TEST);
gl.depthFunc(gl.LEQUAL);
gl.clearDepth(1.0);
gl.clearColor(1,1,1,1);
gl.depthMask(true); // Must be true before clearing depth buffer
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
this.drawSubscene(this.scene.rootSubscene, true);
this.drawSubscene(this.scene.rootSubscene, false);
}
this.stopDrawing(wasDrawing);
};
/**
* Change the displayed subset
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The subset control data.
*/
rglwidgetClass.prototype.subsetSetter = function(el, control) {
if (typeof control.subscenes === "undefined" ||
control.subscenes === null)
control.subscenes = this.scene.rootSubscene;
var value = Math.round(control.value),
subscenes = [].concat(control.subscenes),
fullset = [].concat(control.fullset),
i, j, entries, subsceneid,
adds = [], deletes = [],
ismissing = function(x) {
return fullset.indexOf(x) < 0;
},
tointeger = function(x) {
return parseInt(x, 10);
};
if (isNaN(value))
value = control.value = 0;
if (control.accumulate)
for (i=0; i <= value; i++)
adds = adds.concat(control.subsets[i]);
else
adds = adds.concat(control.subsets[value]);
deletes = fullset.filter(function(x) { return adds.indexOf(x) < 0; });
for (i = 0; i < subscenes.length; i++) {
subsceneid = subscenes[i];
if (typeof this.getObj(subsceneid) === "undefined")
this.alertOnce("typeof object is undefined");
for (j = 0; j < adds.length; j++)
this.addToSubscene(adds[j], subsceneid);
for (j = 0; j < deletes.length; j++)
this.delFromSubscene(deletes[j], subsceneid);
}
};
/**
* Change the requested property
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The property setter control data.
*/
rglwidgetClass.prototype.propertySetter = function(el, control)  {
var value = control.value,
values = [].concat(control.values),
svals = [].concat(control.param),
direct = values[0] === null,
entries = [].concat(control.entries),
ncol = entries.length,
nrow = values.length/ncol,
properties = this.repeatToLen(control.properties, ncol),
objids = this.repeatToLen(control.objids, ncol),
property, objid = objids[0],
obj = this.getObj(objid),
propvals, i, v1, v2, p, entry, gl, needsBinding,
newprop, newid,
getPropvals = function() {
if (property === "userMatrix")
return obj.par3d.userMatrix.getAsArray();
else if (property === "scale" || property === "FOV" || property === "zoom")
return [].concat(obj.par3d[property]);
else
return [].concat(obj[property]);
};
putPropvals = function(newvals) {
if (newvals.length == 1)
newvals = newvals[0];
if (property === "userMatrix")
obj.par3d.userMatrix.load(newvals);
else if (property === "scale" || property === "FOV" || property === "zoom")
obj.par3d[property] = newvals;
else
obj[property] = newvals;
};
if (direct && typeof value === "undefined")
return;
if (control.interp) {
values = values.slice(0, ncol).concat(values).
concat(values.slice(ncol*(nrow-1), ncol*nrow));
svals = [-Infinity].concat(svals).concat(Infinity);
for (i = 1; i < svals.length; i++) {
if (value <= svals[i]) {
if (svals[i] === Infinity)
p = 1;
else
p = (svals[i] - value)/(svals[i] - svals[i-1]);
break;
}
}
} else if (!direct) {
value = Math.round(value);
}
for (j=0; j<entries.length; j++) {
entry = entries[j];
newprop = properties[j];
newid = objids[j];
if (newprop !== property || newid != objid) {
if (typeof property !== "undefined")
putPropvals(propvals);
property = newprop;
objid = newid;
obj = this.getObj(objid);
propvals = getPropvals();
}
if (control.interp) {
v1 = values[ncol*(i-1) + j];
v2 = values[ncol*i + j];
this.setElement(propvals, entry, p*v1 + (1-p)*v2);
} else if (!direct) {
this.setElement(propvals, entry, values[ncol*value + j]);
} else {
this.setElement(propvals, entry, value[j]);
}
}
putPropvals(propvals);
needsBinding = [];
for (j=0; j < entries.length; j++) {
if (properties[j] === "values" &&
needsBinding.indexOf(objids[j]) === -1) {
needsBinding.push(objids[j]);
}
}
for (j=0; j < needsBinding.length; j++) {
gl = this.gl || this.initGL();
obj = this.getObj(needsBinding[j]);
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
}
};
/**
* Change the requested vertices
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The vertext setter control data.
*/
rglwidgetClass.prototype.vertexSetter = function(el, control)  {
var svals = [].concat(control.param),
j, k, p, a, propvals, stride, ofs, obj, entry,
attrib,
ofss    = {x:"vofs", y:"vofs", z:"vofs",
red:"cofs", green:"cofs", blue:"cofs",
alpha:"cofs", radii:"radofs",
nx:"nofs", ny:"nofs", nz:"nofs",
ox:"oofs", oy:"oofs", oz:"oofs",
ts:"tofs", tt:"tofs"},
pos     = {x:0, y:1, z:2,
red:0, green:1, blue:2,
alpha:3,radii:0,
nx:0, ny:1, nz:2,
ox:0, oy:1, oz:2,
ts:0, tt:1},
values = control.values,
direct = values === null,
ncol,
interp = control.interp,
vertices = [].concat(control.vertices),
attributes = [].concat(control.attributes),
value = control.value, newval, aliases, alias;
ncol = Math.max(vertices.length, attributes.length);
if (!ncol)
return;
vertices = this.repeatToLen(vertices, ncol);
attributes = this.repeatToLen(attributes, ncol);
if (direct)
interp = false;
/* JSON doesn't pass Infinity */
svals[0] = -Infinity;
svals[svals.length - 1] = Infinity;
for (j = 1; j < svals.length; j++) {
if (value <= svals[j]) {
if (interp) {
if (svals[j] === Infinity)
p = 1;
else
p = (svals[j] - value)/(svals[j] - svals[j-1]);
} else {
if (svals[j] - value > value - svals[j-1])
j = j - 1;
}
break;
}
}
obj = this.getObj(control.objid);
// First, make sure color attributes vary in original
if (typeof obj.vOffsets !== "undefined") {
varies = true;
for (k = 0; k < ncol; k++) {
attrib = attributes[k];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[attrib]];
if (ofs < 0) {
switch(attrib) {
case "alpha":
case "red":
case "green":
case "blue":
obj.colors = [obj.colors[0], obj.colors[0]];
break;
}
varies = false;
}
}
}
if (!varies)
this.initObj(control.objid);
}
propvals = obj.values;
aliases = obj.alias;
if (typeof aliases === "undefined")
aliases = [];
for (k=0; k<ncol; k++) {
if (direct) {
newval = value;
} else if (interp) {
newval = p*values[j-1][k] + (1-p)*values[j][k];
} else {
newval = values[j][k];
}      	
attrib = attributes[k];
vertex = vertices[k];
alias = aliases[vertex];
if (obj.type === "planes" || obj.type === "clipplanes") {
ofs = ["nx", "ny", "nz", "offset"].indexOf(attrib);
if (ofs >= 0) {
if (ofs < 3) {
if (obj.normals[vertex][ofs] != newval) {  // Assume no aliases here...
obj.normals[vertex][ofs] = newval;
obj.initialized = false;
}
} else {
if (obj.offsets[vertex][0] != newval) {
obj.offsets[vertex][0] = newval;
obj.initialized = false;
}
}
continue;
}
}
// Not a plane setting...
ofs = obj.vOffsets[ofss[attrib]];
if (ofs < 0)
this.alertOnce("Attribute '"+attrib+"' not found in object "+control.objid);
else {
stride = obj.vOffsets.stride;
ofs = ofs + pos[attrib];
entry = vertex*stride + ofs;
propvals[entry] = newval;
if (typeof alias !== "undefined")
for (a = 0; a < alias.length; a++)
propvals[alias[a]*stride + ofs] = newval;
}
}
if (typeof obj.buf !== "undefined") {
var gl = this.gl || this.initGL();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, propvals, gl.STATIC_DRAW);
}
};
/**
* Change the requested vertex properties by age
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The age setter control data.
*/
rglwidgetClass.prototype.ageSetter = function(el, control) {
var objids = [].concat(control.objids),
nobjs = objids.length,
time = control.value,
births = [].concat(control.births),
ages = [].concat(control.ages),
steps = births.length,
j = Array(steps),
p = Array(steps),
i, k, age, j0, propvals, stride, ofs, objid, obj,
attrib, dim, varies, alias, aliases, a, d,
attribs = ["colors", "alpha", "radii", "vertices",
"normals", "origins", "texcoords",
"x", "y", "z",
"red", "green", "blue"],
ofss    = ["cofs", "cofs", "radofs", "vofs",
"nofs", "oofs", "tofs",
"vofs", "vofs", "vofs",
"cofs", "cofs", "cofs"],
dims    = [3,1,1,3,
3,2,2,
1,1,1,
1,1,1],
pos     = [0,3,0,0,
0,0,0,
0,1,2,
0,1,2];
/* Infinity doesn't make it through JSON */
ages[0] = -Infinity;
ages[ages.length-1] = Infinity;
for (i = 0; i < steps; i++) {
if (births[i] !== null) {  // NA in R becomes null
age = time - births[i];
for (j0 = 1; age > ages[j0]; j0++);
if (ages[j0] == Infinity)
p[i] = 1;
else if (ages[j0] > ages[j0-1])
p[i] = (ages[j0] - age)/(ages[j0] - ages[j0-1]);
else
p[i] = 0;
j[i] = j0;
}
}
// First, make sure color attributes vary in original
for (l = 0; l < nobjs; l++) {
objid = objids[l];
obj = this.getObj(objid);
varies = true;
if (typeof obj.vOffsets === "undefined")
continue;
for (k = 0; k < attribs.length; k++) {
attrib = control[attribs[k]];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[k]];
if (ofs < 0) {
switch(attribs[k]) {
case "colors":
case "alpha":
case "red":
case "green":
case "blue":
obj.colors = [obj.colors[0], obj.colors[0]];
break;
}
varies = false;
}
}
}
if (!varies)
this.initObj(objid);
}
for (l = 0; l < nobjs; l++) {
objid = objids[l];
obj = this.getObj(objid);
if (typeof obj.vOffsets === "undefined")
continue;
aliases = obj.alias;
if (typeof aliases === "undefined")
aliases = [];
propvals = obj.values;
stride = obj.vOffsets.stride;
for (k = 0; k < attribs.length; k++) {
attrib = control[attribs[k]];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[k]];
if (ofs >= 0) {
dim = dims[k];
ofs = ofs + pos[k];
for (i = 0; i < steps; i++) {
alias = aliases[i];
if (births[i] !== null) {
for (d=0; d < dim; d++) {
propvals[i*stride + ofs + d] = p[i]*attrib[dim*(j[i]-1) + d] + (1-p[i])*attrib[dim*j[i] + d];
if (typeof alias !== "undefined")
for (a=0; a < alias.length; a++)
propvals[alias[a]*stride + ofs + d] = propvals[i*stride + ofs + d];
}
}
}
} else
this.alertOnce("\'"+attribs[k]+"\' property not found in object "+objid);
}
}
obj.values = propvals;
if (typeof obj.buf !== "undefined") {
gl = this.gl || this.initGL();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
}
}
};
/**
* Bridge to old style control
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The bridge control data.
*/
rglwidgetClass.prototype.oldBridge = function(el, control) {
var attrname, global = window[control.prefix + "rgl"];
if (global)
for (attrname in global)
this[attrname] = global[attrname];
window[control.prefix + "rgl"] = this;
};
/**
* Set up a player control
* @param { Object } el - The player control element
* @param { Object } control - The player data.
*/
rglwidgetClass.prototype.Player = function(el, control) {
var
self = this,
components = [].concat(control.components),
buttonLabels = [].concat(control.buttonLabels),
Tick = function() { /* "this" will be a timer */
var i,
nominal = this.value,
slider = this.Slider,
labels = this.outputLabels,
output = this.Output,
step;
if (typeof slider !== "undefined" && nominal != slider.value)
slider.value = nominal;
if (typeof output !== "undefined") {
step = Math.round((nominal - output.sliderMin)/output.sliderStep);
if (labels !== null) {
output.innerHTML = labels[step];
} else {
step = step*output.sliderStep + output.sliderMin;
output.innerHTML = step.toPrecision(output.outputPrecision);
}
}
for (i=0; i < this.actions.length; i++) {
this.actions[i].value = nominal;
}
self.applyControls(el, this.actions, false);
self.drawScene();
},
OnSliderInput = function() { /* "this" will be the slider */
this.rgltimer.value = Number(this.value);
this.rgltimer.Tick();
},
addSlider = function(min, max, step, value) {
var slider = document.createElement("input");
slider.type = "range";
slider.min = min;
slider.max = max;
slider.step = step;
slider.value = value;
slider.oninput = OnSliderInput;
slider.sliderActions = control.actions;
slider.sliderScene = this;
slider.className = "rgl-slider";
slider.id = el.id + "-slider";
el.rgltimer.Slider = slider;
slider.rgltimer = el.rgltimer;
el.appendChild(slider);
},
addLabel = function(labels, min, step, precision) {
var output = document.createElement("output");
output.sliderMin = min;
output.sliderStep = step;
output.outputPrecision = precision;
output.className = "rgl-label";
output.id = el.id + "-label";
el.rgltimer.Output = output;
el.rgltimer.outputLabels = labels;
el.appendChild(output);
},
addButton = function(which, label, active) {
var button = document.createElement("input"),
onclicks = {Reverse: function() { this.rgltimer.reverse();},
Play: function() { this.rgltimer.play();
this.value = this.rgltimer.enabled ? this.inactiveValue : this.activeValue; },
Slower: function() { this.rgltimer.slower(); },
Faster: function() { this.rgltimer.faster(); },
Reset: function() { this.rgltimer.reset(); },
Step:  function() { this.rgltimer.step(); }
};
button.rgltimer = el.rgltimer;
button.type = "button";
button.value = label;
button.activeValue = label;
button.inactiveValue = active;
if (which === "Play")
button.rgltimer.PlayButton = button;
button.onclick = onclicks[which];
button.className = "rgl-button";
button.id = el.id + "-" + which;
el.appendChild(button);
};
if (typeof control.reinit !== "undefined" && control.reinit !== null) {
control.actions.reinit = control.reinit;
}
el.rgltimer = new rgltimerClass(Tick, control.start, control.interval, control.stop,
control.step, control.value, control.rate, control.loop, control.actions);
for (var i=0; i < components.length; i++) {
switch(components[i]) {
case "Slider": addSlider(control.start, control.stop,
control.step, control.value);
break;
case "Label": addLabel(control.labels, control.start,
control.step, control.precision);
break;
default:
addButton(components[i], buttonLabels[i], control.pause);
}
}
el.rgltimer.Tick();
};
/**
* Apply all registered controls
* @param { Object } el - DOM element of the control
* @param { Object } x - List of actions to apply
* @param { boolean } [draw=true] - Whether to redraw after applying
*/
rglwidgetClass.prototype.applyControls = function(el, x, draw) {
var self = this, reinit = x.reinit, i, control, type;
for (i = 0; i < x.length; i++) {
control = x[i];
type = control.type;
self[type](el, control);
}
if (typeof reinit !== "undefined" && reinit !== null) {
reinit = [].concat(reinit);
for (i = 0; i < reinit.length; i++)
self.getObj(reinit[i]).initialized = false;
}
if (typeof draw === "undefined" || draw)
self.drawScene();
};
/**
* Handler for scene change
* @param { Object } message - What sort of scene change to do?
*/
rglwidgetClass.prototype.sceneChangeHandler = function(message) {
var self = document.getElementById(message.elementId).rglinstance,
objs = message.objects, mat = message.material,
root = message.rootSubscene,
initSubs = message.initSubscenes,
redraw = message.redrawScene,
skipRedraw = message.skipRedraw,
deletes, subs, allsubs = [], i,j;
if (typeof message.delete !== "undefined") {
deletes = [].concat(message.delete);
if (typeof message.delfromSubscenes !== "undefined")
subs = [].concat(message.delfromSubscenes);
else
subs = [];
for (i = 0; i < deletes.length; i++) {
for (j = 0; j < subs.length; j++) {
self.delFromSubscene(deletes[i], subs[j]);
}
delete self.scene.objects[deletes[i]];
}
}
if (typeof objs !== "undefined") {
Object.keys(objs).forEach(function(key){
key = parseInt(key, 10);
self.scene.objects[key] = objs[key];
self.initObj(key);
var obj = self.getObj(key),
subs = [].concat(obj.inSubscenes), k;
allsubs = allsubs.concat(subs);
for (k = 0; k < subs.length; k++)
self.addToSubscene(key, subs[k]);
});
}
if (typeof mat !== "undefined") {
self.scene.material = mat;
}
if (typeof root !== "undefined") {
self.scene.rootSubscene = root;
}
if (typeof initSubs !== "undefined")
allsubs = allsubs.concat(initSubs);
allsubs = self.unique(allsubs);
for (i = 0; i < allsubs.length; i++) {
self.initSubscene(allsubs[i]);
}
if (typeof skipRedraw !== "undefined") {
root = self.getObj(self.scene.rootSubscene);
root.par3d.skipRedraw = skipRedraw;
}
if (redraw)
self.drawScene();
};
/**
* Set mouse mode for a subscene
* @param { string } mode - name of mode
* @param { number } button - button number (1 to 3)
* @param { number } subscene - subscene id number
* @param { number } stayActive - if truthy, don't clear brush
*/
rglwidgetClass.prototype.setMouseMode = function(mode, button, subscene, stayActive) {
var sub = this.getObj(subscene),
which = ["left", "right", "middle"][button - 1];
if (!stayActive && sub.par3d.mouseMode[which] === "selecting")
this.clearBrush(null);
sub.par3d.mouseMode[which] = mode;
};
/**
* The class of an rgl timer object
* @class
*/
/**
* Construct an rgltimerClass object
* @constructor
* @param { function } Tick - action when timer fires
* @param { number } startTime - nominal start time in seconds
* @param { number } interval - seconds between updates
* @param { number } stopTime - nominal stop time in seconds
* @param { number } stepSize - nominal step size
* @param { number } value - current nominal time
* @param { number } rate - nominal units per second
* @param { string } loop - "none", "cycle" or "oscillate"
* @param { Object } actions - list of actions
*/
rgltimerClass = function(Tick, startTime, interval, stopTime, stepSize, value, rate, loop, actions) {
this.enabled = false;
this.timerId = 0;
/** nominal start time in seconds */
this.startTime = startTime;   
/** current nominal time */      
this.value = value;
/** seconds between updates */                 
this.interval = interval;
/** nominal stop time */           
this.stopTime = stopTime;
/** nominal step size */           
this.stepSize = stepSize;
/** nominal units per second */           
this.rate = rate;
/** "none", "cycle", or "oscillate" */                   
this.loop = loop;
/** real world start time */                   
this.realStart = undefined;
/** multiplier for fast-forward or reverse */         
this.multiplier = 1;                
this.actions = actions;
this.Tick = Tick;
};
/**
* Start playing timer object
*/
rgltimerClass.prototype.play = function() {
if (this.enabled) {
this.enabled = false;
window.clearInterval(this.timerId);
this.timerId = 0;
return;
}
var tick = function(self) {
var now = new Date();
self.value = self.multiplier*self.rate*(now - self.realStart)/1000 + self.startTime;
self.forceToRange();
if (typeof self.Tick !== "undefined") {
self.Tick(self.value);
}
};
this.realStart = new Date() - 1000*(this.value - this.startTime)/this.rate/this.multiplier;
this.timerId = window.setInterval(tick, 1000*this.interval, this);
this.enabled = true;
};
/**
* Force value into legal range
*/
rgltimerClass.prototype.forceToRange = function() {
if (this.value > this.stopTime + this.stepSize/2 || this.value < this.startTime - this.stepSize/2) {
if (!this.loop) {
this.reset();
} else {
var cycle = this.stopTime - this.startTime + this.stepSize,
newval = (this.value - this.startTime) % cycle + this.startTime;
if (newval < this.startTime) {
newval += cycle;
}
this.realStart += (this.value - newval)*1000/this.multiplier/this.rate;
this.value = newval;
}
}
};
/**
* Reset to start values
*/
rgltimerClass.prototype.reset = function() {
this.value = this.startTime;
this.newmultiplier(1);
if (typeof this.Tick !== "undefined") {
this.Tick(this.value);
}
if (this.enabled)
this.play();  /* really pause... */
if (typeof this.PlayButton !== "undefined")
this.PlayButton.value = "Play";
};
/**
* Increase the multiplier to play faster
*/
rgltimerClass.prototype.faster = function() {
this.newmultiplier(Math.SQRT2*this.multiplier);
};
/**
* Decrease the multiplier to play slower
*/
rgltimerClass.prototype.slower = function() {
this.newmultiplier(this.multiplier/Math.SQRT2);
};
/**
* Change sign of multiplier to reverse direction
*/
rgltimerClass.prototype.reverse = function() {
this.newmultiplier(-this.multiplier);
};
/**
* Set multiplier for play speed
* @param { number } newmult - new value
*/
rgltimerClass.prototype.newmultiplier = function(newmult) {
if (newmult != this.multiplier) {
this.realStart += 1000*(this.value - this.startTime)/this.rate*(1/this.multiplier - 1/newmult);
this.multiplier = newmult;
}
};
/**
* Take one step
*/
rgltimerClass.prototype.step = function() {
this.value += this.rate*this.multiplier;
this.forceToRange();
if (typeof this.Tick !== "undefined")
this.Tick(this.value);
};</script>
<div id="unnamed_chunk_3div" class="rglWebGL"></div>
<script type="text/javascript">
var unnamed_chunk_3div = document.getElementById("unnamed_chunk_3div"),
unnamed_chunk_3rgl = new rglwidgetClass();
unnamed_chunk_3div.width = 673;
unnamed_chunk_3div.height = 481;
unnamed_chunk_3rgl.initialize(unnamed_chunk_3div,
{"material":{"color":"#000000","alpha":1,"lit":true,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":false,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less","isTransparent":false},"rootSubscene":1,"objects":{"7":{"id":7,"type":"points","material":{"lit":false},"vertices":[[-3.041903,-5.688744,6.528015],[-1.034308,-4.574435,4.795337],[0.4177045,-2.200593,2.452976],[3.053637,-6.670028,7.403646],[5.728846,0.9114407,5.886458],[1.270159,-5.032366,5.285643],[3.94515,8.675275,9.582515],[-5.65379,3.542814,6.746619],[1.175283,1.488634,2.144136],[-1.725082,-3.531672,4.05569],[5.63594,-2.897137,6.41539],[1.680579,3.728849,4.210541],[8.12744,-9.73701,12.7226],[-5.651456,-3.215321,6.578544],[1.64935,8.409482,8.627848],[-0.1889722,-2.068238,2.305063],[0.9847068,-3.59719,3.861272],[1.296311,-6.403102,6.609096],[-5.089932,5.314185,7.426168],[-5.753172,8.90194,10.64629],[5.721361,9.524649,11.15585],[9.983981,6.752047,12.09421],[5.858741,-6.537451,8.835333],[-7.400498,7.71827,10.7396],[-8.83129,-4.904876,10.15133],[-8.578278,1.856538,8.833662],[6.432083,-4.143786,7.716389],[0.1177704,-7.41595,7.483995],[-9.981523,7.922033,12.78239],[4.191308,-5.047933,6.636919],[1.881472,-2.297105,3.13315],[-8.60346,-9.485167,12.84476],[-2.642886,-3.753031,4.697881],[-2.338262,2.819568,3.797029],[-5.579072,-7.00825,9.013412],[-3.63944,-4.245521,5.680666],[-2.468257,5.327007,5.955611],[-3.468616,8.663477,9.385474],[-1.524321,7.434886,7.655135],[-4.210275,0.9561087,4.431767],[-1.33887,-6.564769,6.774125],[-6.639915,-6.806081,9.560921],[-9.350463,-8.989547,13.00935],[2.690791,2.370103,3.722599],[-0.9226267,3.907552,4.137657],[-8.706333,-6.627686,10.98756],[6.52757,-9.695459,11.73078],[-1.030553,-4.609309,4.827812],[4.218264,-3.032192,5.290362],[-2.917331,-8.066046,8.635503],[-3.246051,5.374498,6.357836],[-4.95193,-6.700486,8.391551],[-4.875647,9.289361,10.5387],[0.9211302,-6.243089,6.389417],[6.027972,-3.41845,7.001589],[4.64023,9.228493,10.37771],[7.125714,-2.027264,7.475667],[-3.529205,6.673452,7.615133],[1.514978,4.820969,5.151397],[8.740385,-8.968707,12.56312],[1.696795,-2.678081,3.324339],[-4.390811,-8.741736,9.833472],[-3.001467,-4.344946,5.374697],[4.113747,-4.614683,6.262444],[-7.375525,9.542572,12.10203],[-4.785921,0.9573082,4.982115],[-6.989492,-5.24707,8.796861],[7.757634,8.367951,11.45441],[8.908677,-8.432468,12.30736],[3.200904,4.847893,5.894731],[-1.955935,-2.076082,3.022548],[4.324643,7.224643,8.479269],[-5.765095,5.771752,8.218846],[-4.65678,-4.080753,6.272013],[-8.449775,-4.514556,9.632233],[3.244244,6.620978,7.440597],[-1.371416,9.490365,9.640944],[0.7077767,-7.230311,7.333372],[3.995197,-9.948018,10.76683],[-0.3973687,-6.421848,6.511377],[-2.379451,5.591092,6.158092],[8.66735,1.109708,8.795135],[-4.772872,-6.160582,7.85704],[-3.817273,0.02893158,3.946189],[4.832943,0.3097725,4.945028],[-5.264491,-1.576612,5.585747],[-0.1962945,-5.440763,5.53538],[-7.834706,-7.423203,10.83912],[4.473332,-5.758923,7.360427],[-6.248059,8.160593,10.32635],[-1.804,-7.334364,7.618879],[4.512907,8.890119,10.02001],[-9.703561,-6.215767,11.56697],[-7.662439,9.079106,11.92238],[-3.927617,4.278402,5.893292],[-7.098456,-0.5332761,7.188355],[-1.573451,-4.533666,4.902028],[-7.94484,7.081223,10.68944],[-2.522551,-9.594135,9.970491],[2.297951,-7.867424,8.256933],[9.046039,6.997881,11.48047],[-6.347957,9.74299,11.67144],[8.336466,8.946384,12.26925],[-2.16044,-8.839568,9.154532],[-3.268503,5.297029,6.304096],[2.87926,2.458652,3.916007],[3.515402,-5.105279,6.278688],[5.007239,-5.90161,7.803938],[1.755195,-2.784528,3.440103],[0.9642107,4.82931,5.025131],[8.735256,5.014625,10.12182],[0.8725876,6.30115,6.439403],[-5.514135,4.856709,7.415747],[1.855446,2.519866,3.285179],[-5.775204,-5.120345,7.782732],[3.501835,8.35325,9.112608],[-0.4028516,9.189212,9.252237],[2.876619,1.594393,3.43759],[-3.579884,-7.429121,8.30707],[9.171999,-8.005317,12.21518],[-9.198751,-4.262588,10.18758],[7.821182,-2.290136,8.210701],[-8.251716,-7.402215,11.1303],[-3.768925,-9.533799,10.30039],[8.996481,-4.582543,10.14576],[7.338439,2.021595,7.677209],[-4.079342,-8.710998,9.670704],[-6.659851,-2.154688,7.070805],[-2.847126,-4.538592,5.450224],[7.875822,-9.548076,12.4175],[-2.999162,-6.379948,7.120303],[-5.063008,3.403615,6.182123],[-8.101959,4.376419,9.262547],[-6.623799,-9.202484,11.38246],[-2.423835,-0.9349811,2.783732],[9.482718,9.408211,13.39539],[0.1720216,-6.812796,6.887944],[7.101717,-0.2749321,7.177044],[9.079535,-2.627288,9.504767],[-3.395019,4.332125,5.594056],[-4.496328,-1.405118,4.815737],[-9.289365,7.138794,11.75817],[-2.362686,-6.400041,6.89513],[-4.148294,6.590561,7.851359],[-6.084896,1.643965,6.381895],[3.934465,6.393528,7.573454],[0.01481946,1.903846,2.150546],[0.3613038,6.088556,6.1807],[0.0797528,-5.184014,5.280186],[9.113288,0.2280227,9.170824],[-4.315263,-7.653095,8.842587],[-4.119615,-7.636724,8.734459],[5.947506,-1.503474,6.215566],[4.304824,-8.050396,9.183702],[2.871321,-0.7515653,3.131986],[-8.075539,-9.024336,12.15125],[-8.697647,-2.59446,9.131281],[-0.2368822,-4.731473,4.841792],[-1.979075,-4.02653,4.596704],[-2.595772,9.137438,9.551482],[4.219622,5.327444,6.86927],[-3.487123,6.053884,7.057587],[8.94285,-6.345466,11.01088],[2.862577,-8.254335,8.793656],[7.891471,0.7374632,7.98869],[-0.4485642,-5.098454,5.214925],[-4.893373,-6.438847,8.148855],[-6.581671,-9.636108,11.71209],[-9.552927,-0.9914612,9.656159],[-4.33479,-3.675765,5.770758],[-8.28151,7.073904,10.93725],[-4.224347,-2.192129,4.863182],[-3.244598,9.130409,9.741241],[-3.208265,2.21846,4.026727],[-7.010926,-6.166512,9.390365],[3.015615,4.883705,5.826192],[-9.420569,-4.863139,10.64882],[-1.230813,-8.168084,8.320606],[-0.03298241,-8.623958,8.681805],[-9.487806,-4.619227,10.5998],[-4.368209,-1.397872,4.694177],[5.465052,-4.451326,7.119066],[7.361863,-4.413455,8.641505],[3.832834,8.037712,8.960772],[-4.35842,9.840309,10.80868],[6.017657,0.8697866,6.161876],[4.106423,-6.989424,8.167911],[5.890417,-7.022818,9.220466],[-4.51641,6.103498,7.658371],[-0.5516427,-1.049928,1.551341],[8.081184,-5.971065,10.09748],[5.441123,-8.686115,10.29827],[-7.101994,3.687819,8.064634],[3.945943,4.397244,5.99218],[8.79126,-6.091478,10.74208],[2.33906,4.027888,4.763936],[-0.4214998,1.049843,1.509911],[6.843839,6.346877,9.387278],[-7.71771,-1.128193,7.863578],[-8.968542,-2.095448,9.264213],[4.912339,4.364938,6.647087],[1.511304,2.531263,3.113091],[-5.537863,-4.330377,7.100711],[-6.719273,-1.211207,6.900409],[-6.664672,0.2856675,6.745328],[1.367334,1.044973,1.99037],[-5.647286,6.869443,8.948804],[9.998017,-2.577188,10.37315],[-1.648835,2.701065,3.318796],[2.909037,-2.338545,3.864102],[5.932885,-9.779927,11.48243],[-3.630975,3.640532,5.238078],[8.672204,2.730684,9.146789],[9.545466,-9.000669,13.15781],[-7.553064,3.261241,8.287609],[2.112064,7.35958,7.721673],[-5.842088,0.9769291,6.007028],[1.406342,6.198298,6.434026],[-6.573358,4.240945,7.886358],[9.782089,-8.95244,13.29795],[4.893617,2.869255,5.760218],[-1.453643,-9.948372,10.10362],[2.517053,-0.6747959,2.791219],[-2.247485,-6.761582,7.19515],[3.212517,8.285376,8.942467],[5.815436,2.592171,6.445048],[3.928644,7.688426,8.691728],[4.926432,-7.568654,9.085938],[3.515446,4.771705,6.010618],[-7.852782,-6.149152,10.02388],[2.403997,-3.701366,4.525407],[1.576189,5.751432,6.046763],[5.552059,-1.358846,5.802743],[7.651042,-6.763002,10.26044],[-7.864532,-3.502975,8.667277],[4.481404,-2.777839,5.366505],[7.179633,2.106048,7.548679],[8.610368,-5.580025,10.30898],[-6.479118,4.756538,8.099606],[5.113476,-3.895552,6.50561],[-5.302091,-3.666935,6.523694],[3.34054,-8.8639,9.525121],[-8.469154,-2.081206,8.778269],[8.700352,9.25927,12.74481],[6.973804,7.850301,10.54804],[2.595432,-0.7370614,2.877417],[-6.337823,-1.769511,6.655762],[1.938265,-3.944125,4.506994],[1.191403,8.493879,8.635128],[-5.122084,2.782767,5.91435],[-8.044605,-3.808623,8.956633],[7.487816,-1.720374,7.747714],[2.908376,5.280206,6.110583],[9.27628,-9.411932,13.25269],[-0.0607456,-3.618773,3.754891],[-6.185519,-3.603159,7.22796],[-1.33394,7.23347,7.423105],[6.221944,6.866237,9.319753],[-4.32747,-7.172585,8.436407],[9.528131,-7.976361,12.46626],[-5.781956,5.891248,8.314916],[-8.504884,-9.186503,12.55886],[-2.177039,-8.971324,9.285696],[-5.918207,-8.664197,10.54009],[-5.110329,3.512228,6.281019],[2.930763,6.837099,7.505684],[-4.580698,-7.214484,8.60416],[0.4965666,7.802157,7.881639],[-3.304909,1.927518,3.954459],[3.891589,-0.4033292,4.03821],[-7.542045,-0.7222106,7.642253],[-0.9479923,2.802906,3.123295],[-9.598319,-3.554731,10.28415],[6.277579,-2.941359,7.004255],[9.260766,-8.1562,12.38085],[4.025662,-8.605541,9.553078],[-0.3855755,5.968174,6.063643],[-9.387498,0.7557794,9.470814],[0.6603565,9.759256,9.832556],[8.15541,-9.301682,12.41096],[-3.350342,-4.180826,5.450147],[8.575846,5.928659,10.4735],[-6.744457,-8.411697,10.82794],[-5.136079,7.59676,9.224428],[-7.314196,-9.131824,11.74256],[-9.044357,6.055414,10.93016],[2.239713,-3.213407,4.042561],[9.491448,5.013898,10.78085],[8.817914,0.4353626,8.885108],[-6.993684,-3.431875,7.854259],[4.518788,-2.354707,5.192697],[0.7050254,3.471202,3.680531],[2.210788,8.020265,8.379273],[4.985639,3.825325,6.363152],[-4.546945,-6.841832,8.275589],[4.777431,5.398863,7.278157],[7.143317,-2.43019,7.61136],[-4.58993,-3.731673,5.999403],[8.838252,-3.212069,9.456854],[-2.37239,7.73232,8.149663],[0.1647775,1.929151,2.179168],[-4.443655,-0.8199979,4.628008],[7.47997,-2.68896,8.011271],[9.848856,4.889863,11.04132],[7.136014,9.667871,12.0578],[3.537357,-4.299102,5.656428],[4.486667,-5.973381,7.537338],[-7.668743,8.414947,11.42895],[1.179935,-9.901007,10.02109],[-0.6659736,4.253392,4.419827],[8.69998,-1.017352,8.816158],[0.6480132,-4.972454,5.11324],[0.5026886,-6.480332,6.576275],[3.592055,-3.946205,5.429124],[0.02951845,4.867666,4.96941],[3.984587,3.587956,5.45439],[9.23109,1.867289,9.470997],[-6.547273,9.475042,11.56042],[-3.10706,7.420079,8.106256],[2.741327,0.8737422,3.04603],[3.468986,-1.028993,3.754023],[9.785803,-3.355498,10.39333],[-1.921161,4.548422,5.037758],[-2.130394,-2.094323,3.150361],[-0.258063,-2.448692,2.657572],[-3.07665,6.457794,7.222803],[5.906697,0.140908,5.992405],[-9.940585,-8.714071,13.25708],[8.09173,-1.896151,8.370872],[5.193403,-0.4298671,5.306244],[2.876321,-7.013688,7.646244],[-7.592597,-9.206544,11.97531],[3.9291,9.819221,10.62332],[-4.733774,1.303829,5.010847],[8.595806,-3.815913,9.457752],[9.984606,2.869265,10.43672],[3.31509,-4.724394,5.85745],[6.730862,-3.493743,7.649231],[2.209388,2.552856,3.521146],[5.527162,-0.778797,5.67063],[5.485619,-4.441967,7.129031],[-6.244771,5.592755,8.442516],[7.337988,3.714547,8.285163],[8.802517,-8.638209,12.37348],[5.081563,-2.895099,5.933286],[5.22298,2.580613,5.910929],[-8.182407,-3.961838,9.145925],[-4.204279,-2.154007,4.828634],[7.918489,3.945226,8.903216],[0.4553644,-2.511423,2.741277],[-6.180319,-0.4399933,6.276141],[6.018744,8.138547,10.17159],[5.616885,-0.1150081,5.706367],[0.6250942,-1.424698,1.849462],[-2.891747,5.9385,6.680418],[8.495856,-3.011796,9.069205],[-7.14734,0.558027,7.238499],[1.269702,1.386196,2.129245],[-3.569924,-1.848289,4.142527],[2.064136,-7.771055,8.102467],[-2.564502,-1.392696,3.084845],[-0.2942988,1.815555,2.093527],[7.139223,-3.765381,8.133057],[-5.792174,0.4365976,5.894056],[-3.287969,8.943099,9.580697],[-9.376264,-8.448777,12.66081],[-6.120036,-4.571502,7.70412],[1.308638,6.930924,7.12392],[3.888079,7.303708,8.334346],[9.93331,8.872841,13.35657],[4.465225,6.238637,7.736849],[2.07639,6.604617,6.995167],[6.458971,-0.6686956,6.570043],[-6.855697,-7.051288,9.885405],[9.812969,3.184412,10.36508],[9.249956,-1.951241,9.506263],[1.344536,-7.558774,7.742276],[-1.506272,-3.372873,3.826895],[-2.156861,-5.083169,5.611654],[9.27997,6.341208,11.28401],[-3.700741,0.08915697,3.834506],[-5.589389,-8.675991,10.3689],[8.945135,-2.457257,9.330249],[6.152593,8.193427,10.29498],[4.42568,1.918779,4.926292],[-6.841035,-5.74894,8.991667],[9.83373,9.569848,13.75806],[7.930473,2.045954,8.250959],[9.981849,1.937051,10.21712],[-0.5635027,-9.190079,9.261483],[-6.40193,-5.863482,8.738715],[-8.741631,-4.808119,10.02667],[-3.938452,-4.665986,6.187312],[3.222997,0.298709,3.387763],[9.827929,4.285499,10.76818],[3.58087,6.008193,7.065481],[-4.438308,8.625249,9.751589],[2.115459,-8.485417,8.802129],[6.481131,8.364964,10.6291],[1.407833,6.602992,6.825064],[6.162591,5.45037,8.287584],[3.230556,9.328724,9.922781],[9.172564,-6.421279,11.24139],[-6.997476,5.168924,8.756851],[8.156535,-5.475665,9.874815],[-3.147497,-7.65501,8.33702],[-0.7827289,-3.566025,3.785392],[-5.902806,-3.265833,6.819735],[-9.248297,1.037181,9.359847],[-4.788025,-6.877434,8.439448],[0.9194502,5.687163,5.847154],[-3.163739,-2.728009,4.295495],[4.246062,-1.368219,4.571769],[-3.840909,-8.296825,9.197277],[6.598755,9.207893,11.37228],[-3.798566,-2.826441,4.839201],[-6.170982,0.03028259,6.251554],[-5.110687,-5.309673,7.437187],[-1.977268,2.199503,3.122083],[-1.347537,9.688792,9.833034],[9.525033,-3.724202,10.27599],[0.1176193,-0.4161661,1.089508],[-4.259114,-2.36226,4.971954],[-0.3741134,0.2189756,1.089913],[-4.457707,-8.344141,9.51293],[5.149736,9.343604,10.71554],[6.102945,-5.238568,8.104846],[-3.320227,-8.776774,9.436932],[-2.376962,0.6968153,2.671235],[-6.660213,-0.0507873,6.735059],[8.124605,5.430187,9.823245],[9.180017,0.7250764,9.262745],[-6.51598,-6.067949,8.9598],[5.239398,3.254287,6.248333],[-2.575474,6.135159,6.728539],[0.327503,-7.575526,7.648258],[-6.371567,1.742287,6.680751],[-9.482172,-8.009757,12.45262],[-2.479106,3.567274,4.457736],[5.811526,4.317883,7.308758],[3.677737,9.434646,10.17538],[4.461787,4.406538,6.350206],[8.482294,9.872469,13.05431],[9.704349,-6.554043,11.75287],[-2.606091,-8.656819,9.095726],[3.696018,8.121447,8.978777],[-6.659905,-0.6786188,6.768668],[5.14638,7.365973,9.041172],[-2.071457,-8.210936,8.52704],[0.4747378,1.152311,1.597873],[5.46438,1.377041,5.723259],[9.771651,1.402305,9.922279],[-2.06874,4.045237,4.652271],[-4.074903,0.0921379,4.196824],[-0.2538389,6.898018,6.974746],[6.386501,-5.566225,8.530548],[6.650153,-3.655893,7.654417],[0.6512699,-0.2850249,1.226944],[8.878991,5.436614,10.45912],[9.991362,-9.927397,14.12022],[3.854037,-9.242485,10.06365],[-0.7051888,-1.1297,1.665387],[-5.708809,3.112792,6.578752],[8.30994,-4.484261,9.495457],[-6.258942,4.780906,7.939233],[8.721509,-8.213727,12.02206],[4.763298,9.071959,10.29512],[-6.833291,7.580168,10.2544],[-4.839773,1.621297,5.201154],[3.145813,1.545001,3.644608],[8.089378,-5.826228,10.01913],[-4.593126,7.384719,8.753906],[6.570402,-4.900798,8.257603],[-0.8645831,8.469284,8.571831],[-4.707701,3.852051,6.164474],[6.644974,-9.409233,11.56241],[3.119814,-9.38994,9.94506],[-8.936648,-7.037745,11.419],[8.945281,8.605769,12.453],[2.722086,-4.962667,5.747853],[2.499687,-9.795776,10.15902],[5.872604,1.863649,6.241848],[7.261394,-4.143202,8.419855],[3.226915,-4.497487,5.624977],[4.581976,3.055317,5.597273],[-1.220884,3.446172,3.790338],[-5.261484,0.3427791,5.366629],[-8.673301,-2.957429,9.218055],[-8.275361,1.425867,8.456636],[-3.850512,-6.168251,7.339875],[5.210942,1.015744,5.402375],[3.764964,6.849532,7.879786],[6.483525,8.703627,10.89905],[9.264662,-9.206747,13.09955],[4.868081,-8.711724,10.02957],[4.243781,-0.6278858,4.404987],[9.707273,-7.326916,12.20307],[5.16256,-4.053179,6.639298],[6.749945,-1.148111,6.919532],[-4.466007,5.279652,6.987127],[6.173522,-0.3672192,6.26476],[9.864346,5.003574,11.1059],[-9.648923,-9.275389,13.42142],[7.341928,7.068431,10.24044],[-3.712182,-6.179323,7.277659],[1.448943,4.955284,5.258733],[3.086207,2.411135,4.04206],[2.266205,-7.941073,8.318434],[6.973834,7.733959,10.46176],[0.2676625,-3.630218,3.774934],[-3.571945,5.506196,6.63905],[-4.100378,-4.859603,6.436524],[-8.31678,-6.671593,10.70883],[-4.672333,-8.22708,9.513966],[-9.521071,8.353239,12.70541],[4.346221,-2.995405,5.372345],[-8.370203,-1.889578,8.638913],[-4.163905,-9.671752,10.57738],[-7.702432,-6.681789,10.24567],[0.4029975,-5.696977,5.798099],[-4.896043,-5.017735,7.081589],[-0.2799949,-8.811684,8.872663],[-0.8567588,-7.850997,7.960666],[-2.519452,6.076905,6.654052],[2.201564,-0.9480549,2.597247],[1.987661,1.118784,2.490477],[-5.444849,6.719845,8.706474],[5.913722,5.172811,7.920232],[1.247939,1.911179,2.491979],[8.314281,9.432113,12.61317],[-0.5717853,7.614768,7.701405],[-6.994335,2.217281,7.405204],[6.644557,-5.732944,8.832711],[5.415682,-1.387766,5.679393],[8.780219,-3.282593,9.426964],[-1.415098,1.249111,2.136067],[-3.484509,-3.532863,5.061909],[5.839069,7.241543,9.355996],[-5.351105,-6.70214,8.634408],[-8.496284,-8.568311,12.10796],[-9.996485,3.074358,10.50625],[3.168969,-6.576552,7.368405],[-0.6096366,-7.786394,7.873981],[-2.318658,-3.568492,4.371534],[6.489139,-2.581017,7.054826],[-8.633159,-6.016428,10.57019],[-3.627624,-7.634729,8.511683],[-1.882741,-6.610703,6.945942],[3.753515,6.211136,7.325782],[-0.008680229,-2.114865,2.339386],[-4.800402,-7.233454,8.738805],[5.865456,8.880745,10.68977],[-2.618692,-2.199729,3.563195],[-5.08169,3.695873,6.362629],[8.107106,-9.597925,12.60339],[7.058955,-6.641448,9.743597],[9.6342,6.033559,11.41147],[-7.021075,-2.893579,7.659523],[9.683857,4.770857,10.8415],[1.733956,-5.074326,5.454849],[6.550436,-1.228368,6.739222],[6.782342,6.730544,9.607309],[-1.92065,5.611259,6.014575],[-7.722196,4.473583,8.98027],[6.377421,6.772868,9.356455],[-2.380948,6.603236,7.09025],[-3.897458,0.6243784,4.071857],[-2.21576,-6.526428,6.964471],[5.9151,9.713268,11.41648],[6.906823,-4.359046,8.228333],[-1.655639,7.454005,7.700866],[-0.6503811,5.488253,5.616397],[-7.166523,-9.708503,12.10843],[4.284157,2.552316,5.086091],[-7.928997,-4.178684,9.018336],[4.298151,2.501638,5.072701],[3.730556,1.062436,4.005723],[-9.915712,-2.881019,10.37408],[3.033545,8.635671,9.207454],[-9.927127,2.93804,10.40096],[-4.085099,0.5455044,4.240944],[-8.247934,-3.590074,9.050804],[-7.74084,6.450017,10.12538],[-4.853656,-5.747512,7.588931],[3.884557,-6.468047,7.610875],[9.324539,-6.358109,11.33016],[-8.045154,-7.532508,11.06631],[-2.137468,-5.235144,5.74243],[6.451472,2.431592,6.966644],[5.301823,0.4973912,5.418185],[-1.377372,2.560553,3.074669],[4.074917,-8.824548,9.771264],[4.941579,7.836434,9.318203],[4.563457,6.559665,8.053219],[-8.815145,1.95146,9.083775],[-9.547643,-8.365757,12.73355],[3.370833,2.486763,4.306565],[-5.890809,2.923168,6.651807],[-5.359407,-3.270872,6.357818],[3.740028,-2.109923,4.409034],[-1.900173,2.354209,3.186371],[-7.728576,-5.601779,9.597438],[-5.060655,-2.988601,5.961708],[8.653278,-7.631227,11.5808],[4.711415,-2.521313,5.4364],[7.794825,-6.265563,10.0507],[2.593345,5.492701,6.155908],[-8.062908,3.213484,8.737103],[-3.241966,4.468866,5.610803],[7.608766,-1.393468,7.799684],[9.547028,5.219147,10.92636],[-6.593624,-2.950086,7.292386],[0.5211073,1.628813,1.981057],[4.105025,-5.613443,7.025808],[-9.198268,-2.164328,9.502234],[6.899184,4.634,8.370944],[-9.441072,5.380912,10.91275],[3.237174,1.730742,3.804571],[1.97566,-9.29222,9.552413],[-7.110626,-1.181263,7.277113],[-9.972219,-4.023222,10.7996],[-2.150389,1.036545,2.588166],[9.210653,6.779508,11.48033],[2.880449,4.573123,5.496402],[5.511725,-9.803714,11.29123],[-6.438959,2.58029,7.00843],[0.4396359,-1.113428,1.559808],[-3.347202,3.255885,4.775411],[2.106855,1.573014,2.813043],[-0.5253502,-6.478419,6.576162],[0.006626537,-3.999259,4.122392],[-9.51368,-6.457925,11.54188],[6.089665,2.136668,6.530648],[-2.674759,-5.526684,6.220818],[-3.215332,2.31397,4.085684],[-9.205165,-3.621442,9.942328],[-2.959229,6.246099,6.983609],[8.89782,-9.810793,13.28243],[4.805547,-7.935825,9.331162],[8.557773,-6.684807,10.90514],[-8.639548,-6.840934,11.06527],[-6.87972,2.088519,7.258957],[-5.017636,0.523186,5.142994],[3.44954,-6.937848,7.812366],[3.747239,9.450079,10.21498],[6.50101,-8.485168,10.73598],[-2.699696,4.197338,5.089794],[8.811143,-4.978623,10.16971],[-0.2257219,-5.900246,5.988644],[2.171347,-9.325655,9.627179],[0.563679,3.934262,4.098311],[-6.016602,8.479056,10.4448],[-9.704993,-5.398738,11.15048],[1.86465,2.894588,3.585465],[3.093728,-0.6539504,3.316445],[4.541197,-7.624595,8.930673],[-1.970688,-6.526587,6.89057],[9.348808,-2.187417,9.653238],[7.979973,-1.59059,8.198168],[5.303439,-0.1274092,5.398397],[-0.1442079,-9.190939,9.246305],[3.945781,-5.511327,6.851563],[-0.5156551,6.619653,6.714589],[-7.239165,4.491143,8.577639],[7.22899,7.909896,10.76219],[-1.040086,-9.173748,9.286518],[-0.1158635,4.013148,4.137485],[8.383542,8.665289,12.09839],[9.657877,-4.631774,10.75769],[-0.6500534,4.95644,5.097928],[-4.213508,5.53698,7.029353],[6.379797,8.240031,10.469],[3.003027,-6.973892,7.658547],[-7.19857,8.06829,10.85895],[-3.558447,-1.626291,4.038239],[0.5741965,-4.617876,4.759672],[-6.607122,-0.3677726,6.692482],[0.3108141,-5.332093,5.43395],[-5.616878,4.506331,7.270236],[-3.416178,4.397916,5.657909],[-4.990015,-4.478836,6.779396],[-2.18823,7.758232,8.122716],[-1.591539,7.932154,8.151813],[-0.8714625,3.2324,3.493974],[1.4787,9.656899,9.820502],[-6.354913,-1.009002,6.511759],[-0.7315499,7.895685,7.992309],[7.202275,-6.842649,9.984719],[-9.614973,-7.082855,11.98393],[7.835856,3.018016,8.456303],[-5.68793,-4.310289,7.206326],[2.582724,-1.908109,3.363234],[-0.2722411,-4.853008,4.962439],[2.539633,4.276136,5.072975],[-7.876253,5.273706,9.531386],[-4.04956,-9.251941,10.14876],[-8.891197,-9.6392,13.15171],[4.785835,3.677992,6.118156],[5.458369,-4.245291,6.986866],[-1.363871,-1.038316,1.984501],[4.99707,-6.338666,8.133228],[-6.775901,-3.152407,7.539928],[-3.162222,1.276424,3.553717],[-2.198236,4.393112,5.013151],[7.866491,2.749083,8.392802],[-6.101921,9.083381,10.98823],[-1.846125,-9.999534,10.21758],[-6.458807,-9.877138,11.84373],[-1.147547,0.5980721,1.635406],[-5.398292,-5.461408,7.743935],[-4.618064,0.2597009,4.732226],[-6.333912,8.216282,10.42237],[3.412006,1.08476,3.717323],[0.01951983,0.7143508,1.229096],[-2.398876,-9.967411,10.30067],[-5.593959,4.258611,7.101278],[-2.304334,-1.279092,2.818871],[-1.991133,-1.682117,2.791797],[1.69195,3.728177,4.214499],[3.15646,-3.161589,4.578087],[-5.376165,-4.374596,7.002874],[-1.197415,-6.318246,6.507998],[-0.7435824,-6.081043,6.207415],[8.767601,-9.76186,13.15921],[2.488117,-5.995066,6.567461],[1.615795,9.754131,9.937498],[0.1692109,-0.4491381,1.109215],[-3.583409,6.059944,7.110819],[-9.155149,5.648639,10.80388],[-3.028139,-0.983506,3.337201],[1.318361,5.445743,5.691589],[9.814194,-8.132664,12.78509],[-0.1375082,9.462806,9.516491],[-9.967447,-7.00067,12.22127],[-3.548011,-0.2580783,3.695265],[-9.64862,-9.644037,13.67857],[3.855137,0.6912323,4.042262],[8.141286,3.109621,8.77213],[8.900586,6.432388,11.02706],[-4.508034,-5.986088,7.560134],[9.895302,1.763266,10.1008],[7.714651,-2.186776,8.080708],[1.598289,-2.493044,3.125667],[8.42073,-3.890986,9.329977],[-1.73958,2.57094,3.261268],[-1.107453,-4.396563,4.642868],[8.055679,9.937401,12.83144],[-5.614483,-9.447405,11.03521],[-4.124797,-8.763474,9.737167],[-9.190564,-2.899293,9.688776],[8.07986,8.668722,11.89247],[-3.644092,6.291048,7.338712],[0.3610228,-6.149076,6.24031],[-1.298302,-2.969804,3.39195],[5.510404,9.403978,10.94529],[-2.936406,-5.096062,5.96593],[-6.138706,1.963304,6.522137],[7.770466,-5.990461,9.862341],[-1.247762,-0.2406318,1.617038],[-4.645027,-9.478976,10.60317],[5.685143,7.116283,9.163097],[5.653201,7.899747,9.765485],[3.011887,1.094175,3.356886],[1.918545,0.2803718,2.18161],[5.492558,7.531768,9.375272],[-3.52991,-6.734656,7.669149],[5.409197,3.279938,6.404483],[7.830295,0.2067829,7.896599],[-1.579557,4.511926,4.8839],[-1.964247,5.57632,5.996133],[-3.038318,-4.046948,5.158407],[-5.879727,-1.362046,6.117709],[-3.067429,3.095804,4.471367],[4.20452,-1.41848,4.548634],[3.806179,-5.344582,6.637134],[6.415355,-5.133148,8.276834],[-3.374391,-2.726172,4.4518],[9.403295,0.1892082,9.458211],[-2.558141,-3.442218,4.403743],[-4.707647,3.467009,5.931449],[-5.288159,-4.685555,7.135758],[9.039901,0.295299,9.099835],[-9.604459,-5.01099,10.87914],[8.636831,7.620301,11.56131],[-0.2397155,6.698038,6.776516],[-6.815151,-6.500652,9.47126],[-1.830068,7.751893,8.027514],[-8.519813,5.09041,9.974943],[0.5539469,4.412366,4.558051],[-2.652797,4.529901,5.343906],[-5.434065,-7.937228,9.671021],[-7.002343,-1.431934,7.216872],[1.120567,-8.21558,8.351731],[7.947289,8.369146,11.58456],[6.226867,-1.99257,6.61394],[-9.384178,-9.056503,13.07987],[-9.330549,0.7833547,9.416623],[0.8363303,7.69295,7.802623],[-8.196482,3.421277,8.937978],[-6.859464,5.086468,8.59793],[-0.7466572,0.837224,1.502811],[-6.24236,-0.2529239,6.327008],[0.3246347,9.860938,9.916828],[-4.73045,-3.265831,5.834621],[-4.059841,-8.825868,9.766179],[-7.60701,5.864132,9.656844],[0.9435672,2.581935,2.925185],[4.627648,-5.615918,7.345315],[-3.565342,-5.155512,6.347517],[-0.235951,-9.500598,9.555994],[3.251614,1.414377,3.684217],[7.885126,8.810743,11.8661],[-8.273724,1.753649,8.516441],[4.546553,6.488283,7.985547],[-2.05517,-9.287406,9.564499],[7.380751,0.5589051,7.469127],[2.286138,8.545805,8.902652],[-7.643538,-5.896768,9.705439],[0.02049373,8.466708,8.525583],[5.743358,-2.5869,6.377947],[9.508097,4.16538,10.42853],[8.173876,9.777085,12.78294],[4.318057,-8.01525,9.159141],[6.750407,1.304072,6.94756],[0.1130702,7.46422,7.531757],[2.772575,-6.888364,7.492445],[9.711888,4.26615,10.65461],[6.38161,7.036725,9.551987],[-1.652751,1.832821,2.66286],[-2.373686,-0.9494892,2.745162],[-3.685802,-6.009079,7.119984],[-3.248025,-6.376012,7.225178],[-3.420874,-6.001142,6.979691],[9.074219,2.729853,9.528565],[-7.028851,8.631681,11.17634],[-7.003469,0.9171225,7.133701],[9.180106,-4.204924,10.14671],[-0.1882417,-2.888589,3.062578],[-0.8728055,0.7364391,1.517937],[1.973257,0.07500994,2.213452],[-1.814758,2.671747,3.381061],[-5.759304,-5.951791,8.342266],[4.71324,-5.12867,7.036895],[-8.910841,3.626144,9.672229],[-2.666356,-5.591552,6.274943],[-9.593958,6.012966,11.36661],[4.855868,0.7376733,5.012347],[-2.941201,9.606171,10.096],[2.771227,-0.8854892,3.076327],[2.985066,8.20958,8.792487],[-4.480999,7.468823,8.767136],[4.23523,-3.185581,5.393061],[-5.899967,6.425003,8.780107],[8.736698,-7.658174,11.66094],[-4.270422,5.71483,7.203873],[5.885511,2.774945,6.583279],[9.099333,6.770942,11.38611],[-9.226543,0.9493978,9.329012],[-5.037013,-6.986801,8.671037],[-1.981286,-0.4949469,2.273866],[5.407158,5.489035,7.769611],[7.918068,6.60015,10.35653],[4.008928,9.621238,10.4709],[6.303428,-2.358822,6.804208],[9.308292,-2.100394,9.59458],[8.769506,-5.41965,10.35745],[7.279232,-2.909902,7.902831],[4.151949,6.083782,7.433108],[9.69567,-9.46055,13.58337],[6.790306,-1.090174,6.949584],[4.606875,-5.780144,7.458777],[9.027784,2.285765,9.366195],[-6.427591,-6.664939,9.313181],[1.124015,-9.327994,9.448539],[1.900753,-0.8173243,2.298016],[6.192993,-4.948568,7.990087],[3.333558,-4.479685,5.672759],[0.06135894,-7.977505,8.040171],[2.917633,-9.112993,9.620771],[-9.090697,-2.175755,9.400781],[-2.792823,1.92535,3.5365],[3.154982,-1.804095,3.769439],[8.117728,0.2101777,8.18179],[-0.8287021,-9.003699,9.096887],[-1.983764,2.847507,3.611594],[-3.369469,3.040385,4.647285],[-1.059345,-6.825825,6.979548],[-5.393869,0.9008496,5.559258],[4.732203,0.1526559,4.839117],[-5.726658,-2.374327,6.279494],[5.452698,8.986813,10.55911],[-1.670591,3.244498,3.783867],[-3.985764,-5.511343,6.87468],[3.667432,-9.032956,9.800222],[-3.440557,4.381483,5.659932],[1.190985,5.97629,6.175313],[6.206937,-6.398638,8.970431],[9.081357,-5.820096,10.83257],[-9.307009,-2.439155,9.673153],[-0.09111921,-3.441918,3.585401],[5.517047,-8.299805,10.01621],[-0.0802413,7.929295,7.992506],[-3.543942,-1.563878,4.000655],[0.5349792,-2.208843,2.48298],[-3.894433,-8.353938,9.271186],[5.235967,-8.237469,9.81179],[2.879822,9.453516,9.932892],[-0.3568964,2.343184,2.572526],[9.95607,1.806061,10.16785],[-5.928872,-1.90463,6.30707],[3.388931,-3.041308,4.662018],[-2.474621,-3.866964,4.698634],[5.082086,-6.26346,8.12764],[-2.367739,-9.714494,10.04876],[-3.501747,6.835564,7.745138],[9.607226,-9.672073,13.66923],[6.708597,-4.585756,8.187455],[8.273551,1.108695,8.40719],[-1.944442,6.901646,7.239722],[6.452111,6.770212,9.405611],[2.709183,8.425428,8.9066],[-0.7923817,-8.696505,8.789599],[-3.65811,-7.247869,8.180058],[-7.349576,5.554318,9.26643],[-6.4676,-5.299846,8.421295],[7.884805,6.217567,10.091],[-4.07881,-1.546884,4.475438],[2.404434,7.92279,8.339779],[-5.605089,5.38428,7.836294],[7.898169,-3.179137,8.572514],[3.270693,0.5111493,3.458136],[3.18966,3.082465,4.547034],[8.760292,-0.3149024,8.822804],[-4.064149,-2.899484,5.091593],[9.257897,-9.026245,12.96849],[-6.744058,-9.765205,11.90973],[4.974947,3.208445,6.003684],[8.716031,-7.911117,11.81334],[3.39794,-4.279674,5.555323],[2.246479,4.616759,5.230787],[-9.588288,-6.107582,11.41218],[9.587842,8.648985,12.95113],[6.507744,7.653815,10.09612],[-5.370017,-2.663354,6.07705],[8.258739,-3.487273,9.020412],[9.171244,8.519852,12.55785],[-9.250583,1.960431,9.508762],[-8.067475,-7.527614,11.07922],[-2.632964,-5.416318,6.104835],[1.122984,-1.745171,2.303631],[-1.96018,1.237204,2.524477],[7.639277,-3.493596,8.459537],[9.989133,-5.761209,11.57473],[-9.452566,-3.74241,10.21551],[-8.568869,-8.911624,12.40333],[-5.287637,7.45386,9.193429],[-6.722613,0.1701316,6.798711],[-8.980945,3.401561,9.655463],[-2.079581,2.852481,3.668965],[-7.133776,2.699265,7.692645],[-2.73675,8.973417,9.434618],[-5.375832,-0.7364092,5.517415],[8.901663,1.017351,9.015243],[-2.197465,7.793251,8.158653],[-5.534708,-5.684674,7.996781],[-5.375805,3.592693,6.542685],[8.77229,-7.914689,11.85729],[-3.753534,5.091156,6.403818],[2.225373,-2.493938,3.488841],[-4.893599,7.084071,8.667835],[-6.465463,4.930681,8.192303],[-4.164371,-6.672914,7.929046],[7.532988,-1.103157,7.678728],[9.314239,3.362468,9.952951],[-8.80388,2.705448,9.264326],[-2.970509,0.9689353,3.280664],[3.915448,-7.601888,8.609264],[-8.65624,-2.292533,9.010338],[5.668028,-5.325031,7.841077],[-0.6850252,7.77432,7.868247],[6.750068,-6.743845,9.593897],[-5.686194,-4.086271,7.073217],[1.618654,-0.8071108,2.066753],[7.810856,-5.359083,9.52519],[-5.889349,-9.321992,11.07176],[5.521703,6.507422,8.592772],[-7.466011,-0.4967357,7.549044],[8.840009,-4.777817,10.09818],[4.235464,7.682645,8.829619],[7.999798,-2.678203,8.495265],[7.489807,4.617035,8.855181],[-5.704573,-9.509583,11.13438],[9.864993,-4.718752,10.98111],[-3.411546,-8.419056,9.138882],[-8.685945,-5.338408,10.24423],[-0.00544942,-5.719899,5.806657],[9.855411,9.527435,13.74413],[-1.381924,-7.261514,7.459176],[0.5713955,6.268704,6.373629],[1.177797,-7.121001,7.286691],[7.472037,9.811078,12.3729],[-1.81141,-1.440739,2.521297],[-7.488823,-1.293421,7.665208],[6.305003,2.229256,6.761852],[3.468614,-4.549651,5.807806],[-8.704941,-1.53771,8.896097],[4.989507,-6.508982,8.262084],[-9.723303,7.74386,12.47036],[6.592782,7.061751,9.712523],[-9.499548,-8.433242,12.74209],[7.42166,7.530411,10.62018],[1.964994,2.93893,3.674032],[4.131686,-5.467602,6.925713],[3.295332,7.661968,8.400295],[1.257933,-8.019286,8.178713],[-5.912465,-2.480287,6.48915],[-7.062162,-9.333306,11.74669],[8.272673,9.868433,12.916],[-5.324579,4.493381,7.03858],[-1.879703,-7.829208,8.113555],[-8.135071,-9.126273,12.26655],[-0.003457479,9.762694,9.813777],[2.305376,-6.843302,7.290099],[4.751412,9.623779,10.77929],[-3.090933,0.9387196,3.381577],[3.350428,9.0967,9.745528],[-8.449457,-3.87362,9.348703],[-5.626029,-7.75807,9.635344],[4.46301,3.285737,5.631566],[7.742088,5.108827,9.329525],[5.738068,1.388753,5.987825],[5.883861,7.41778,9.520677],[-7.81236,-0.4013245,7.886319],[5.901677,3.200914,6.787904],[1.075541,0.4974476,1.550562],[1.903574,-5.622466,6.019611],[4.570216,6.595976,8.086642],[4.626105,8.543122,9.766564],[4.630077,-2.791591,5.498236],[-5.925111,2.998719,6.715598],[2.815404,-1.812589,3.494564],[-7.91539,7.387867,10.87354],[-2.723855,8.178656,8.678122],[8.161918,-0.1731985,8.224773],[-9.870852,-7.164104,12.23757],[-0.6171274,8.892155,8.969462],[3.817678,6.292655,7.427798],[3.251376,-5.173341,6.191518],[-2.066576,-7.165576,7.524375],[5.696543,1.410834,5.953238],[-6.429852,-9.541462,11.54913],[-6.212714,-6.056158,8.733548],[9.427762,2.606971,9.832547],[-5.754764,7.469862,9.482412],[-3.947596,-6.655341,7.802376],[9.505708,3.204309,10.08098],[8.373031,-6.633307,10.72886],[-6.220325,9.296812,11.23046],[2.200725,6.536607,6.969248],[8.695606,6.098835,10.66815],[-1.796201,-2.207254,3.01634],[-6.211329,-3.729821,7.313834],[-1.574739,-0.517197,1.935793],[-1.572644,-8.880605,9.074048],[-5.206865,6.195558,8.154531],[-8.877953,5.409533,10.44419],[-1.04114,-1.245311,1.906507],[8.181485,-2.387106,8.581082],[-9.838553,6.3787,11.76796],[3.362804,-1.300751,3.741712],[4.996677,-9.470849,10.75471],[8.542764,-8.880493,12.36293],[-1.970145,-9.152161,9.415069],[2.557416,3.29369,4.288213],[8.162758,8.717281,11.98422],[-2.759927,-8.417048,8.914252],[-5.20458,6.487766,8.377276],[3.319734,7.991781,8.711441],[-0.6204314,8.566844,8.647297],[-6.480998,-0.3294943,6.565965],[-8.687576,8.218565,12.00078],[-1.335947,9.278155,9.427031],[-0.7707871,-8.356985,8.451823],[4.49949,5.403916,7.102655],[-0.7795951,1.783215,2.188064],[-2.852823,3.523244,4.642397],[9.409116,-7.916888,12.33728],[-6.608869,-1.630788,6.880161],[8.100712,-7.228788,10.90307],[9.295043,-3.088643,9.845687],[-2.760188,-0.8609962,3.059404],[-3.168357,-5.176823,6.151259],[7.93496,7.199534,10.7609],[8.949181,-3.569728,9.686629],[2.639462,-4.539268,5.345252],[-6.606663,1.949112,6.96039],[7.265266,0.4672716,7.348635],[4.436617,-4.954403,6.725302],[6.012594,1.716501,6.332271],[-6.264765,6.522858,9.099174],[5.191232,-8.67757,10.16116],[3.857165,2.011643,4.463678],[-9.798589,-9.291443,13.54043],[-2.116024,-1.80201,2.953777],[4.078426,-2.45347,4.863443],[-4.494982,1.810789,4.948112],[5.259814,4.504542,6.996895],[6.484181,2.080446,6.882795],[7.520535,-8.636292,11.49539],[0.444719,-6.175886,6.272108],[9.576776,2.208214,9.878808],[-6.16544,1.099843,6.342106],[-6.030856,1.933008,6.411532],[7.901597,3.515479,8.705965],[2.572615,-4.882704,5.608845],[5.780593,-9.922232,11.52675],[7.94468,9.948829,12.77095],[-8.909275,3.96881,9.804419],[-1.007452,-8.681378,8.796664],[6.777146,-0.06794503,6.850863],[4.532455,-7.402934,8.737653],[3.178219,-7.754284,8.439787],[-6.610613,-6.221748,9.132927],[-0.9663562,-8.458042,8.5716],[2.293383,-6.853894,7.296264],[-8.692527,2.556907,9.115799],[1.860486,-6.610416,6.93967],[7.821259,8.802216,11.81741],[3.022288,-4.371338,5.407663],[-3.210806,5.752988,6.663794],[3.383077,6.088171,7.036408],[-3.040464,-0.8164917,3.303192],[5.950975,1.066282,6.127892],[-9.485615,4.456438,10.5279],[-6.791014,-1.407247,7.007012],[4.380429,-1.794514,4.838227],[-6.257799,-1.064268,6.425941],[-2.367898,2.007349,3.261348],[6.539535,-2.019082,6.916806],[-5.988785,-8.150484,10.16346],[0.363221,9.150419,9.212063],[-2.240979,-1.605507,2.932514],[0.2511202,2.904192,3.081784],[-9.606741,-4.580312,10.68966],[9.373603,-2.85262,9.848953],[2.363446,6.968741,7.426253],[-8.263205,9.625187,12.72497],[3.141925,5.289403,6.232935],[-1.119785,9.248422,9.369484],[-2.32428,-7.629002,8.037659],[9.126457,-6.820256,11.43714],[8.5087,-7.842581,11.61482],[1.720424,-9.883305,10.08165],[-2.037954,7.393804,7.734442],[-1.499701,3.226796,3.696121],[3.514142,-8.813437,9.540748],[-1.73907,-3.30899,3.869596],[7.472314,2.574836,7.966508],[1.064879,-1.790243,2.310614],[-0.9414656,-0.05657198,1.374612],[8.378823,-1.623291,8.593006],[4.429738,-2.444871,5.157516],[8.164859,7.43662,11.0891],[2.50984,-2.977719,4.02071],[-5.904207,3.419781,6.895982],[-9.370181,7.358607,11.95614],[-5.453716,-0.7137897,5.590395],[-0.1486552,-9.346381,9.400901],[4.71305,-6.411644,8.020101],[5.09884,-3.738369,6.40106],[-6.937679,9.874072,12.10903],[-0.6793345,-7.2308,7.331164],[1.468221,-5.873801,6.136548],[6.601977,0.3127479,6.684602],[7.167895,7.989069,10.77979],[-1.5683,-5.819903,6.109897],[-1.184008,9.652411,9.776037],[4.113077,-9.313243,10.23005],[-8.496347,8.991637,12.41118],[-2.570241,-7.977763,8.441021],[9.20792,-1.4525,9.375262],[-3.259465,7.283296,8.041798],[8.160625,9.926207,12.88896],[-8.50174,-7.139536,11.14686],[-4.335377,4.12097,6.064477],[-7.064709,-8.5657,11.14815],[1.623756,-9.532433,9.721309],[-2.480522,9.775187,10.13446],[-0.7038638,-1.367937,1.83485],[-0.05145437,2.935561,3.101639],[2.619832,7.534458,8.039376],[-0.4654122,8.259791,8.333113],[8.820605,9.118629,12.72606],[-4.439463,1.902761,4.932477],[-9.60649,5.309828,11.02175],[-6.796242,5.787354,8.982337],[-5.292218,7.27561,9.052186],[-5.07464,-8.274911,9.758388],[2.973712,-3.746497,4.886635],[2.230644,-0.7734568,2.563982],[-6.010297,-9.361705,11.16983],[-7.365532,-0.4497307,7.446699],[-8.878706,0.7979292,8.970402],[0.07744554,-5.128317,5.22548],[9.415407,9.107209,13.13739],[-8.449056,-5.580894,10.17511],[-9.804804,-1.386515,9.952719],[-9.399762,6.104842,11.25276],[-9.004338,1.705922,9.218908],[-2.914934,8.57675,9.113588],[-1.327777,8.309918,8.474535],[-1.51345,0.5454924,1.894226],[3.188648,6.806147,7.582289],[-3.982998,7.786207,8.8028],[-9.04269,4.807896,10.2901],[-2.064922,-8.382483,8.690795],[7.103187,9.098195,11.58587],[-5.520591,-3.78564,6.76816],[9.649742,-5.729366,11.26691],[-8.850355,-8.21248,12.11502],[7.550882,2.737827,8.093919],[-8.327228,1.812205,8.580606],[-9.150446,1.844343,9.387878],[-8.131415,-7.07409,10.82417],[-6.962738,-8.815537,11.27801],[9.3894,3.022621,9.914488],[6.296584,-6.044283,8.785233],[-2.624699,3.153152,4.222726],[5.494321,8.930143,10.53257],[4.222104,-0.7573259,4.40451],[-2.864507,-3.711174,4.79356],[1.816288,-0.9339371,2.274014],[9.215417,2.490663,9.598297],[-8.940159,0.3440569,9.002489],[-6.428889,-0.9691069,6.577977],[4.733486,-8.174622,9.498964],[5.184879,-0.0189466,5.280467],[0.8592133,-5.006399,5.177092],[4.759046,2.282806,5.372124],[7.407982,-6.263682,9.752534],[6.483167,7.963304,10.31725],[1.266945,7.5522,7.72275],[-3.223343,6.045249,6.923509],[-8.500842,1.312464,8.659496],[-9.510337,-7.433927,12.11238],[-7.762035,9.529612,12.33137],[-3.939458,-5.178607,6.583107],[2.572414,1.094069,2.968889],[-5.229203,1.856552,5.638382],[-1.256112,-0.319253,1.636992],[5.780801,-4.52494,7.408963],[-1.478089,-2.26287,2.881897],[-0.05919762,-5.878773,5.963512],[0.4364631,0.06928477,1.093298],[-2.56362,0.04059999,2.752053],[1.128901,7.710281,7.856389],[1.409576,-1.320544,2.175027],[-5.241354,-3.339791,6.294918],[-3.572006,7.986213,8.805614],[4.076373,9.023294,9.951716],[7.568694,-5.990261,9.704039],[-5.427725,-6.362413,8.422618],[0.05472959,4.535043,4.64431],[2.47958,-6.240643,6.789251],[1.61984,-3.243638,3.760993],[9.511855,-3.893113,10.32626],[3.015509,-7.981602,8.59065],[-8.864993,-6.654546,11.12974],[1.503707,6.686837,6.926393],[-2.907784,4.052628,5.087141],[-2.299131,-8.970528,9.31431],[-3.333631,-6.537766,7.406449],[-9.501578,-4.781371,10.6837],[-5.376635,-3.660203,6.580676],[-5.662708,-3.261605,6.610924],[2.991921,7.301076,7.953446],[-2.698103,-0.1608225,2.881948],[-2.31315,3.113932,4.0059],[6.064968,2.171185,6.51904],[-2.717006,1.248646,3.152973],[6.674471,-9.193753,11.40498],[-5.659673,1.169478,5.865116],[-7.514844,-7.819026,10.89082],[-1.153888,2.922908,3.297704],[-1.97773,-0.03683682,2.216478],[5.787683,3.293335,6.733746],[0.6199435,-5.103593,5.237461],[4.351086,-9.075044,10.11377],[-2.903973,0.8470521,3.185994],[-2.570694,0.05032901,2.758804],[-9.709691,-8.681232,13.063],[4.258589,-5.923326,7.363516],[-2.93714,5.545444,6.354427],[8.270168,6.306636,10.44841],[-6.410083,-5.573445,8.552921],[8.447042,2.859041,8.973664],[-4.130958,-4.889364,6.478479],[-7.046788,1.014657,7.18935],[6.77867,1.200614,6.956424],[-3.285796,2.308695,4.138421],[3.915301,4.075672,5.739398],[5.389317,-8.195727,9.859751],[1.336507,-9.284917,9.433765],[-2.712457,-7.167543,7.72859],[2.439994,5.251873,5.876712],[-0.780812,1.474106,1.944905],[-5.978556,0.8224179,6.117148],[8.581075,2.349345,8.952891],[-9.17388,-2.383042,9.530948],[7.071238,-8.985312,11.47773],[6.328835,0.2006792,6.410494],[1.244294,4.740616,5.00217],[8.115897,4.705912,9.43469],[3.441625,5.118588,6.248578],[9.724601,-3.037157,10.2368],[6.417716,-3.339608,7.303428],[-5.739388,-7.077181,9.166628],[5.900565,0.2259179,5.988966],[1.833227,-0.3877175,2.123922],[-8.40111,3.045663,8.991924],[-9.746615,2.022884,10.00443],[1.707681,9.589739,9.791796],[-4.541277,7.75071,9.038623],[-6.593448,9.466991,11.58005],[-9.961249,0.4303156,10.02056],[-7.131291,-9.578871,11.98374],[5.513147,-8.100536,9.849542],[-6.330183,7.773947,10.07499],[-5.519099,5.409619,7.792588],[-7.596915,-7.660318,10.83483],[1.848001,-3.307423,3.918438],[-2.856755,2.368268,3.843143],[-3.86025,9.45337,10.26001],[9.284269,-9.669873,13.44262],[-6.845943,-0.2644165,6.923645],[5.952668,-3.602175,7.029219],[-2.436067,0.3849398,2.661316],[-2.413492,-8.405097,8.801739],[3.312071,-6.339179,7.221842],[-2.976202,7.152602,7.81137],[-6.823694,2.270509,7.260717],[1.20349,0.5668027,1.664227],[3.744715,-3.788441,5.419887],[-0.516924,3.486955,3.66416],[3.15053,1.797049,3.762342],[-5.768541,3.403114,6.771798],[-0.4776953,-3.958258,4.110474],[3.254518,-2.755924,4.380298],[-7.696437,6.984113,10.44093],[-7.122421,6.713025,9.838373],[-2.469656,-8.394485,8.807188],[2.728012,5.218865,5.973157],[-4.643526,8.923222,10.10872],[-3.877324,-3.67028,5.431813],[4.691253,3.222629,5.778685],[-9.615191,-2.34368,9.947097],[3.003509,-8.839321,9.389071],[-3.419664,-0.9618399,3.690425],[-3.574147,1.498411,4.00247],[-8.425005,-2.42517,8.823954],[-5.362371,7.489447,9.265357],[-1.450054,4.533877,4.86402],[7.035821,-3.751331,8.035873],[-7.047315,-3.979764,8.154947],[6.668497,-5.592751,8.760577],[-5.897414,-7.433497,9.541298],[-6.229954,5.712007,8.51113],[-2.874685,0.6963587,3.122295],[-6.970351,-4.091882,8.144279],[5.479525,0.2320946,5.57486],[-5.254241,-2.076076,5.737345],[-1.42186,8.886942,9.055353],[-5.97752,-7.350823,9.527084],[2.943729,-5.792251,6.573866],[7.398184,-9.034913,11.72019],[0.3650422,-1.960165,2.230584],[4.78452,-7.23552,8.731803],[-5.773848,-6.110262,8.465969],[-3.237452,0.5294455,3.429491],[-8.235762,-0.857227,8.34042],[-6.960609,-1.084925,7.115275],[-5.829744,-3.206735,6.728229],[9.50271,8.321171,12.67057],[-6.699147,-4.225137,7.98313],[-3.500146,2.310258,4.311417],[4.647257,6.48816,8.043209],[9.634215,0.5225286,9.700059],[7.883838,-8.245569,11.45183],[0.9985148,7.23731,7.373987],[-8.165297,5.130656,9.695138],[8.936748,9.703244,13.22945],[-4.993115,-6.506947,8.262661],[-1.468014,0.7796824,1.939837],[-6.450632,-2.707376,7.066862],[-0.2251106,0.559646,1.167852],[1.964695,5.817009,6.220742],[-7.20397,-6.444257,9.717285],[-0.9189651,1.430779,1.972721],[7.31983,-9.054502,11.68606],[-6.966276,6.096881,9.311336],[-6.87423,4.893336,8.497046],[-9.643953,8.767756,13.07208],[-5.75395,-2.344052,6.293053],[2.275598,7.839189,8.223821],[0.3471159,-6.377162,6.464417],[-9.472334,5.310347,10.90527],[6.704564,7.175653,9.87123],[-4.569597,2.528271,5.317271],[-5.245763,-5.730927,7.833362],[-8.490084,-8.418215,11.99783],[-3.702251,5.767375,6.925986],[4.451914,4.39353,6.334244],[-0.9171414,-8.818001,8.921787],[9.988764,5.867216,11.62754],[-6.139517,-8.345427,10.40864],[-4.852426,1.146938,5.08542],[-6.695445,-8.488696,10.85758],[1.863106,6.375288,6.716805],[-7.908834,6.119433,10.04973],[-5.168276,-4.509001,6.931246],[-8.166243,-6.039346,10.20594],[-4.562191,-2.822143,5.456929],[-8.50326,9.708751,12.9447],[-7.524321,-3.374022,8.306589],[5.990341,6.79079,9.110379],[0.6617281,-5.468819,5.598737],[9.741453,-8.815157,13.17585],[-7.239854,0.4968413,7.325458],[-5.575512,9.607954,11.15344],[5.762807,0.1593744,5.851098],[4.456146,-2.004892,4.987668],[2.681029,6.663498,7.251905],[5.994745,7.110694,9.354087],[1.844592,-8.083296,8.351179],[-4.01845,-7.24964,8.348965],[5.1214,1.832215,5.530438],[0.6415901,4.653295,4.802582],[-6.983442,-4.894094,8.586072],[5.709424,2.767353,6.423065],[7.426394,3.101728,8.109997],[9.197556,-2.512491,9.586847],[1.109242,-8.279277,8.412897],[3.108273,4.565097,5.612617],[6.818488,2.73529,7.414417],[-5.297457,-2.256478,5.844206],[1.126267,0.7037817,1.662464],[6.856763,-6.24603,9.328885],[-4.129322,9.575757,10.47599],[9.171185,1.767037,9.393246],[8.486213,-2.455209,8.890661],[-5.740573,-3.728183,6.917625],[-0.8691576,8.920292,9.018151],[-2.212044,4.2179,4.866602],[-9.138198,-9.324805,13.09422],[-1.90288,0.8802238,2.322874],[-0.1424326,-2.569984,2.76136],[-3.028126,6.604321,7.333935],[7.444137,-9.589713,12.18104],[4.359885,3.196768,5.497993],[1.594765,4.517179,4.893688],[-0.1869513,-5.010824,5.113053],[-9.057297,-6.362833,11.11397],[-0.8831979,-9.650416,9.742206],[3.117315,8.181191,8.811898],[7.713401,-2.15249,8.070301],[-0.7803711,3.8561,4.059371],[3.60577,-7.332383,8.231976],[2.165248,7.277583,7.658427],[-7.157082,6.967681,10.03855],[-3.900046,-9.214128,10.05537],[-1.712709,-0.2780989,2.002676],[-8.981532,2.870923,9.482095],[-9.492363,-3.388485,10.12851],[3.524731,-1.638453,4.013509],[1.827759,6.941776,7.247686],[-3.462868,3.244987,4.849886],[-3.974323,-5.087107,6.532526],[8.485037,-6.65832,10.83185],[8.36235,-5.661519,10.14799],[-7.290123,8.407651,11.17294],[4.487246,-2.226798,5.108229],[9.676881,1.130774,9.79391],[7.465437,-1.842117,7.754106],[-9.735509,4.827431,10.91257],[-4.189157,-0.1835597,4.310769],[5.994042,-6.47486,8.879885],[4.61795,-6.108637,7.722753],[4.213061,-6.141859,7.514807],[-8.915042,-5.561236,10.55487],[-8.932236,-3.678375,9.711605],[-8.668296,-2.471481,9.069044],[2.537138,8.381706,8.814197],[5.258278,1.358679,5.522273],[0.5611183,-4.862759,4.996126],[-8.406111,-7.371167,11.22483],[4.507922,-0.1913222,4.621468],[0.1436699,1.726958,2.000756],[-4.287156,-7.97215,9.106859],[-0.3479221,-9.045591,9.107347],[-4.788922,-0.5541776,4.923503],[-3.192727,-8.446579,9.085053],[-5.152153,9.493639,10.84776],[-9.157755,5.748074,10.8584],[7.393227,-3.337508,8.173051],[1.740336,8.484138,8.718335],[-8.709914,9.21586,12.71985],[-3.523593,5.145523,6.31602],[-2.584728,1.768667,3.287704],[2.692714,-3.724238,4.70326],[5.285431,-7.723823,9.412397],[-3.549942,-4.904076,6.136126],[-0.7855898,5.405947,5.553504],[7.297949,8.894229,11.54848],[4.691267,8.053204,9.373477],[-4.145994,-4.334078,6.080584],[-9.684891,3.725404,10.42477],[-2.705096,-9.38273,9.815964],[-3.428976,4.797397,5.981045],[2.509857,-6.245091,6.804451],[9.735221,-8.672724,13.07634],[1.347327,7.629707,7.812024],[-5.49386,2.96741,6.323608],[-2.563181,-0.09926005,2.753134],[5.754973,8.736856,10.50963],[8.592189,-9.303908,12.70387],[6.821752,6.367276,9.385015],[-2.408646,4.054587,4.820918],[-5.128446,6.404586,8.265572],[7.877079,3.977938,8.881012],[0.2965428,-6.636016,6.717488],[2.215573,-7.420852,7.808829],[8.678236,0.08833563,8.736108],[2.266876,8.293337,8.655528],[-8.378764,6.541925,10.6771],[3.31499,4.230096,5.466523],[9.804229,6.593243,11.85722],[6.61696,6.852157,9.577903],[9.537497,-3.864853,10.33929],[-7.079579,-6.102604,9.400118],[8.606421,1.763846,8.842038],[8.705057,5.797539,10.50664],[7.260916,-5.860518,9.384379],[-8.624317,-0.5296321,8.698238],[5.132399,0.6523948,5.269453],[4.317173,0.3672789,4.44667],[3.317909,-3.711586,5.077833],[-6.307526,2.275755,6.779671],[7.020151,6.480116,9.605958],[-5.634207,4.193324,7.094241],[-6.330136,4.603362,7.8906],[-2.184178,-1.912928,3.070819],[-3.560078,-1.854693,4.136912],[1.621059,-1.050693,2.175267],[-2.905499,7.01885,7.661996],[-9.148784,-4.873044,10.41378],[0.1969298,-6.755055,6.831512],[8.266631,7.837508,11.4352],[-6.909878,-7.750872,10.4318],[4.731191,-1.919491,5.202751],[4.870971,-8.030171,9.445106],[7.319163,-9.988521,12.42339],[-9.459288,-2.182966,9.759275],[-5.810571,-3.694773,6.958023],[4.422137,8.681581,9.794138],[-8.659067,-0.4776238,8.729695],[-0.02197551,2.611602,2.796596],[4.521853,-6.045698,7.615617],[-6.907262,4.784529,8.461796],[9.620697,-0.7951117,9.705154],[5.889467,3.632259,6.991361],[-8.995126,-7.08438,11.49351],[2.87243,-3.356023,4.52921],[-7.35588,2.786874,7.929416],[-7.905954,3.364197,8.649967],[8.218412,-1.476289,8.409621],[-6.601443,6.675242,9.441288],[3.019452,3.926336,5.05304],[6.758572,-0.9130701,6.892894],[-0.2999273,0.3041517,1.087412],[-9.96429,-2.997817,10.45342],[-8.582337,-2.643689,9.035795],[-1.250003,-7.663529,7.828933],[-3.36399,5.948492,6.90659],[5.896552,-0.9787188,6.060298],[1.493479,5.500616,5.786818],[-7.939598,-9.225152,12.21232],[-9.148705,1.264726,9.289689],[-9.007415,-1.522298,9.189718],[0.3650674,4.512632,4.636499],[8.787636,8.174843,12.0437],[2.110367,-0.4628544,2.380732],[-7.948184,-9.63125,12.52735],[-0.2728362,7.542438,7.613331],[4.787845,2.071367,5.311687],[-2.966353,-9.600279,10.09775],[-5.515339,-3.9374,6.84997],[-1.828264,5.825093,6.186619],[-2.609492,-2.340165,3.644972],[-9.466926,1.075727,9.580181],[-9.306424,-2.141585,9.601871],[3.457961,3.91825,5.320731],[-3.643557,1.202157,3.964932],[0.7261388,1.692997,2.096071],[3.382061,-0.2094166,3.533015],[-3.274327,9.963273,10.53508],[-8.157183,1.063578,8.286786],[1.092556,-0.3908903,1.531821],[-5.184419,9.936253,11.25199],[-8.694625,-7.847409,11.75493],[-8.775327,9.150921,12.71793],[4.507646,-8.98075,10.09816],[6.03896,-5.937541,8.527803],[9.882386,9.785092,13.94308],[-4.427474,-2.490897,5.177557],[3.608245,-7.157462,8.077666],[-6.965281,-8.257265,10.84885],[-4.318688,5.981488,7.445084],[9.122805,-9.001202,12.85485],[4.132092,9.289428,10.21605],[8.161431,-7.645637,11.22785],[0.6748958,3.048122,3.27819],[3.653655,9.284356,10.02738],[-3.790133,-0.0155499,3.919865],[-9.668244,-1.151604,9.787807],[-2.765187,1.280298,3.207089],[1.184197,-6.835649,7.009166],[-9.183529,4.109454,10.11063],[-6.500918,-2.794175,7.146282],[-9.208233,-4.284757,10.20543],[9.749644,5.378233,11.17949],[-7.883329,5.65545,9.753512],[1.147047,-2.077212,2.57498],[4.421711,9.180385,10.2387],[5.234329,2.160676,5.750367],[-2.018193,8.457796,8.752566],[6.945543,9.084182,11.4788],[8.882716,-5.39202,10.43918],[6.943231,6.423381,9.511481],[-9.051401,0.7070491,9.133881],[-2.885267,2.525401,3.962627],[9.070497,1.61535,9.267323],[-2.15434,-9.201559,9.50315],[-5.863227,-1.333914,6.095634],[-0.4234815,-6.863136,6.948523],[-0.2275611,9.693198,9.747301],[-8.205126,-1.953052,8.493439],[-1.252419,9.973078,10.10103],[1.767045,-7.006361,7.294624],[4.263052,-5.186948,6.788081],[-9.595199,-7.089537,11.97202],[9.02706,-6.001201,10.88587],[7.545714,-7.160617,10.45047],[9.676149,-3.436427,10.31683],[-3.7532,7.045309,8.045053],[-9.46703,-8.41836,12.70801],[7.809815,0.8347263,7.917701],[4.881371,9.219986,10.48026],[5.969259,-7.705784,9.798529],[-0.7320172,7.397609,7.500698],[-4.716808,5.049247,6.981631],[6.230476,8.299274,10.42578],[1.701212,9.886738,10.08175],[2.878656,-6.306434,7.004126],[-8.797114,-6.933459,11.24554],[-1.881626,1.459852,2.582961],[-5.947207,-9.258037,11.04901],[-0.6368671,-1.630587,2.01604],[-0.2918648,-4.303746,4.428026],[7.435112,-3.680739,8.356359],[9.44147,-0.8009455,9.528005],[9.795087,6.137315,11.60217],[4.500471,4.0172,6.114911],[7.209281,3.245232,7.969019],[8.288128,-1.933182,8.569145],[-3.874629,0.7026983,4.062824],[6.617379,-2.699894,7.216587],[1.744847,-7.515306,7.779738],[-1.306021,-9.668832,9.807752],[-5.989498,7.852294,9.926359],[1.413848,-0.8467727,1.92769],[8.4804,3.368994,9.179723],[2.230171,6.892572,7.313085],[3.368983,-3.628169,5.051105],[-5.73092,-0.1397433,5.81919],[7.632129,-6.910865,10.34454],[-5.881357,6.162466,8.577083],[2.290205,9.498583,9.821818],[5.122074,-2.182298,5.656683],[-6.887595,-3.661613,7.864246],[8.158918,2.599089,8.621091],[-0.2921467,3.274863,3.43658],[-8.342263,-1.435831,8.523788],[2.582807,-9.525176,9.919671],[0.3724781,-8.115988,8.185842],[8.532197,2.982601,9.09364],[4.313929,0.8972671,4.518303],[7.779123,-6.91765,10.45795],[7.348828,6.313142,9.739663],[-9.177999,-0.8653681,9.272784],[3.537408,9.666174,10.34157],[6.733428,-8.746271,11.08315],[-0.030704,6.798696,6.871914],[-8.979212,-9.036183,12.77806],[5.883031,-9.4387,11.16687],[7.464423,-4.614001,8.832135],[3.242755,-4.140954,5.35378],[-7.480505,0.8367196,7.59329],[3.230056,-7.111603,7.874526],[-1.377859,4.305697,4.630067],[-2.004286,-0.5537186,2.307329],[-0.449096,2.574614,2.798272],[-1.884219,6.071489,6.435313],[0.9480044,-1.367679,1.941457],[-6.833823,-8.861439,11.23504],[8.808703,4.424764,9.908168],[6.798666,-0.5070146,6.890494],[6.283653,8.427235,10.55948],[1.693074,2.464782,3.153039],[-5.16184,-7.663837,9.29403],[8.191246,-3.268884,8.875929],[-6.334036,-8.776405,10.86947],[3.187806,-2.477021,4.159055],[-4.520676,2.57581,5.298236],[5.036611,7.885295,9.409853],[8.897625,-8.78321,12.54243],[3.929451,9.781341,10.58845],[4.124517,-9.674155,10.56413],[6.647089,9.82032,11.90052],[1.340202,9.472386,9.618849],[-9.399371,-2.904678,9.888646],[3.324319,6.275362,7.17156],[8.722598,-1.178565,8.858484],[8.966121,7.368019,11.64813],[5.540499,7.787099,9.609165],[1.081987,8.181644,8.313243],[2.085647,-9.165516,9.452862],[5.566218,6.853775,8.885776],[-4.87812,1.951048,5.348144],[-9.299462,-7.041002,11.70708],[7.168901,0.9973215,7.306695],[3.478347,4.011228,5.40267],[4.529773,-3.796454,5.994323],[4.254022,-9.12751,10.11969],[7.323007,0.2468884,7.395092],[-2.78158,1.976785,3.555962],[1.723064,8.329702,8.564631],[0.7630793,9.267803,9.352778],[-3.652495,9.295551,10.03733],[6.827487,-0.6469295,6.930591],[-4.240455,-2.705044,5.128228],[-7.640895,-9.234953,12.02779],[-5.51957,-8.795327,10.43185],[-5.280343,6.504825,8.437699],[7.635977,-2.071634,7.974949],[-2.491942,-2.75112,3.844273],[8.111567,-9.814979,12.77229],[-5.8247,0.6120947,5.941531],[7.376803,6.38062,9.804567],[2.903178,-2.000878,3.664963],[6.332159,-5.13226,8.211963],[1.936882,9.684022,9.926319],[0.8890902,3.920027,4.142112],[-3.415543,-4.268127,5.557233],[-4.274508,0.1698053,4.393206],[7.980322,-2.901489,8.550097],[5.515811,5.653863,7.961805],[3.009576,1.676532,3.587243],[-1.606355,-0.6335289,1.995429],[6.499166,3.346721,7.378326],[-5.729475,-6.100859,8.42896],[-5.937152,-5.684358,8.280199],[-4.284411,-8.995752,10.01398],[9.715587,0.8149091,9.800853],[1.086373,7.67498,7.815722],[-6.512734,-1.965229,6.875887],[6.758175,6.393722,9.356956],[-4.335945,4.823946,6.56284],[9.374019,-1.639588,9.568724],[4.109079,-9.574764,10.46712],[-7.178827,6.059341,9.447284],[-9.793056,-9.076901,13.39007],[-1.74497,-3.668664,4.18378],[-2.533357,0.3574499,2.746938],[-3.862628,-8.986589,9.832532],[8.723919,3.519015,9.459928],[0.9209474,-4.988805,5.170718],[-8.251559,2.903401,8.804429],[6.070555,-9.018433,10.91713],[-6.83289,2.900652,7.490138],[4.858063,-1.209929,5.105361],[9.292111,-0.5787583,9.363669],[-1.966451,-9.110035,9.373348],[9.324001,2.572565,9.723945],[1.12903,3.87883,4.161734],[8.850939,-0.3142475,8.912792],[2.499379,2.07109,3.396514],[5.93015,-1.467184,6.190259],[-2.072576,-6.517341,6.911679],[5.818531,-5.515587,8.079418],[5.850018,8.622422,10.46751],[1.200441,4.627161,4.883818],[6.557717,-6.342932,9.178041],[-4.660152,3.547728,5.941666],[-6.089442,9.046528,10.95084],[-6.43328,5.21052,8.338861],[1.605479,-4.90378,5.255912],[-5.243912,-4.948011,7.278834],[8.230492,1.712941,8.466119],[8.918855,-5.560242,10.55757],[-5.216883,5.803392,7.867352],[7.493085,-1.429891,7.693563],[-4.751748,9.633662,10.78826],[-0.6516693,1.429783,1.862512],[-6.440693,5.422465,8.478541],[-9.978743,-9.708,13.95781],[4.931269,-5.124073,7.181472],[-9.18252,-9.579339,13.30723],[-3.565254,-1.814815,4.123662],[6.98989,-9.973835,12.22031],[-1.242719,-3.118166,3.502472],[-2.813835,7.157004,7.755022],[0.8432949,7.438224,7.552372],[-5.05622,-1.731807,5.437326],[-8.109459,2.957004,8.689488],[-1.176257,9.071175,9.201619],[9.295687,-5.831374,11.01884],[-4.378824,3.209398,5.520357],[-6.063727,-9.434914,11.25995],[-7.069687,0.6632364,7.170799],[9.567055,1.796476,9.785493],[-4.660811,-5.280585,7.113912],[-4.957922,-0.3099502,5.067254],[5.499045,-1.302452,5.738979],[-4.853376,9.963265,11.12753],[7.00768,-4.06804,8.164345],[-5.789689,6.803094,8.989026],[-2.349271,-9.695623,10.02618],[-4.73419,2.790305,5.585549],[-3.152541,-0.9710784,3.446957],[0.7728165,9.762372,9.843839],[9.663936,-7.053321,12.00587],[-6.933073,-7.866616,10.53333],[-7.339271,-4.277293,8.55337],[2.748998,-1.52041,3.296762],[-9.944197,6.868835,12.12716],[-2.447401,-8.547337,8.946884],[-9.80699,1.282377,9.940902],[-1.105541,5.891223,6.076901],[0.5164777,8.678631,8.751307],[-0.2445812,0.08158284,1.032703],[-4.254823,-0.9485197,4.472494],[6.722893,1.411044,6.941782],[-5.755992,0.9357589,5.916679],[4.457723,6.380291,7.847255],[6.221089,-6.368893,8.959059],[-9.543841,0.8297982,9.631899],[-2.99072,-2.634945,4.109421],[8.059358,-5.888727,10.03147],[-1.162493,-1.152232,1.91808],[9.434093,9.092658,13.14072],[-7.188223,6.706398,9.881616],[-5.689917,-8.926414,10.63278],[-0.07317816,-1.740349,2.008524],[0.3561026,6.918749,6.999707],[-1.804646,-9.823473,10.0378],[-2.691843,-7.740145,8.255656],[6.235535,2.593819,6.827137],[9.541096,6.586262,11.63664],[-9.845352,2.717472,10.26234],[-1.315269,-9.03635,9.186162],[-0.4380977,5.066693,5.182983],[-7.521572,-5.111415,9.148804],[7.484688,-6.974573,10.27936],[-5.470227,1.863917,5.864944],[-3.510298,-6.482588,7.439498],[-3.051395,5.632531,6.48355],[1.407881,3.808188,4.181438],[-7.813371,0.6502113,7.903894],[0.8627448,-2.125735,2.502614],[-6.101016,-6.676701,9.09949],[-3.417766,-9.157271,9.825311],[5.280532,-0.6020511,5.408002],[1.114885,-5.95552,6.140944],[9.527917,6.092309,11.3533],[-1.380333,-6.47958,6.700021],[-9.568471,-4.109879,10.46168],[-9.850197,-4.939722,11.06468],[7.892043,-5.512895,9.678655],[-6.937071,-6.424359,9.507647],[5.54477,-6.150361,8.340948],[-6.803483,-6.827603,9.690384],[-0.7230152,6.90731,7.016671],[0.2352304,5.563879,5.657923],[0.003069225,2.103029,2.328677],[0.9044126,-7.45222,7.573213],[7.173554,-0.9856173,7.309673],[-1.981353,-5.148424,5.606428],[8.311226,3.043813,8.907372],[4.96744,4.685901,6.901676],[-4.257229,-2.582718,5.078822],[1.608892,-8.236609,8.451642],[4.481431,-0.2090146,4.596402],[4.828547,-4.154133,6.447611],[8.156691,-7.79523,11.32684],[3.285307,-8.333434,9.013288],[9.720474,-2.42429,10.06801],[9.464022,6.864293,11.73398],[0.004337854,2.725754,2.903404],[9.590227,-1.501706,9.758462],[-7.627338,9.333088,12.09474],[-7.167334,3.024849,7.843494],[4.71315,8.752306,9.990829],[-5.644831,4.753794,7.447327],[9.825429,3.632952,10.52318],[-5.714248,0.1671683,5.803496],[-6.193954,7.727718,9.954029],[9.838564,1.560915,10.01168],[-6.988922,4.453662,8.347463],[4.384174,6.645096,8.023608],[7.535374,0.3730427,7.610586],[3.450114,-0.2299299,3.599466],[-1.216119,-4.362785,4.638194],[9.992735,3.3671,10.59208],[-5.024012,-8.669656,10.06994],[-3.303157,5.487618,6.482653],[4.364093,-3.59158,5.739752],[1.986302,-2.019601,3.004028],[-1.288948,-2.437577,2.933116],[-0.1256664,8.982393,9.03876],[0.405513,4.006213,4.148997],[-5.733899,7.320365,9.352291],[-9.080676,-8.762447,12.65856],[8.198731,4.237147,9.28292],[8.037115,-8.275572,11.5793],[0.5744638,-0.9545674,1.497066],[-5.852849,-5.086843,7.818684],[1.594472,6.186516,6.466477],[9.09143,-7.980768,12.13865],[-8.069131,-6.914532,10.67341],[-5.143257,2.231025,5.694784],[6.691443,-8.106286,10.55875],[7.252676,-5.586685,9.209362],[-2.347368,-7.647577,8.061983],[-1.342586,9.198783,9.349874],[-4.920785,-1.009701,5.121877],[-2.631611,-8.55188,9.003334],[4.45489,-4.313223,6.280919],[-4.331608,-9.114691,10.14103],[6.099317,-2.854352,6.80801],[2.168861,6.071543,6.524384],[6.930819,6.399509,9.486304],[9.592655,2.548592,9.975688],[-0.08237074,-0.7589402,1.258084],[-5.567073,4.523337,7.242436],[6.7435,4.892031,8.390874],[5.836103,-8.169882,10.08995],[-1.086303,-1.43025,2.055644],[-6.089606,3.06208,6.889095],[5.122951,-8.950181,10.361],[8.049097,-1.384505,8.228293],[-3.731539,-4.640787,6.038318],[7.503011,4.778867,8.951689],[1.891532,8.891033,9.144855],[-5.095928,1.989572,5.561194],[0.6509949,-9.724121,9.797056],[2.127204,4.254426,4.86057],[-7.783567,8.854142,11.8313],[-3.560956,-5.471616,6.604467],[-0.2545398,-9.938697,9.992121],[6.908202,-9.386826,11.69768],[2.913597,-8.833198,9.354915],[4.263649,9.976884,10.89573],[8.163691,-0.7209648,8.256248],[-6.713606,0.279117,6.793409],[-9.228199,6.881586,11.55491],[7.96655,-0.6702566,8.056994],[-4.64616,-5.617944,7.358539],[1.387929,1.496449,2.272819],[5.898706,6.857209,9.100332],[-7.837719,-6.674828,10.34327],[6.013824,6.215673,8.706358],[5.759248,3.063235,6.599421],[-2.80477,1.819968,3.489845],[-0.8029559,2.512514,2.820898],[-7.662939,3.009042,8.293067],[-7.164707,-4.353498,8.443102],[6.428756,5.533668,8.5411],[1.057776,-4.19933,4.444465],[-3.734309,1.056851,4.007742],[-6.440798,-2.222136,6.886346],[0.4504617,1.132124,1.576268],[-3.654546,-2.556334,4.570619],[-3.101261,-7.091191,7.804025],[-5.539345,6.677545,8.733496],[-2.465519,-4.193901,4.966648],[-4.666554,6.308002,7.909969],[1.445589,5.005058,5.304746],[-6.907197,-3.959167,8.023987],[3.030773,-4.559714,5.565661],[-8.144963,-1.143366,8.285391],[3.364002,6.566556,7.445546],[5.151148,6.073227,8.026108],[-6.973269,9.977259,12.21361],[-7.518181,2.657354,8.036453],[7.625678,5.104453,9.230732],[5.634518,-4.961909,7.574189],[4.941858,3.901241,6.37508],[-3.894489,7.218748,8.263012],[-3.975139,3.162305,5.177056],[-4.909255,-3.408205,6.059426],[3.937864,-7.748788,8.749313],[-6.51581,8.319539,10.61464],[-7.141914,6.813171,9.921],[-2.444913,1.387689,2.983836],[4.102855,2.034686,4.687575],[1.753845,-9.889277,10.09326],[-9.553993,-4.765972,10.72349],[-5.231761,-5.543099,7.687475],[0.07842987,5.81964,5.905452],[-1.718462,-1.412537,2.438929],[-7.203867,-0.896296,7.327964],[-4.128525,0.6044146,4.290691],[8.974292,-0.6554713,9.053594],[-2.911673,8.436675,8.980831],[-1.420046,8.557204,8.731683],[-9.50174,9.150206,13.22911],[-3.609398,0.5118783,3.780182],[3.559828,-5.260384,6.429931],[4.840874,-4.946115,6.992719],[-7.372895,-4.928119,8.924458],[-7.850753,-2.263053,8.231387],[1.825685,-4.411115,4.877608],[-3.949073,-7.594954,8.618498],[-2.181359,-2.996744,3.839115],[8.140444,-0.6305048,8.225835],[-5.630362,-1.702116,5.966421],[-8.38429,-2.673295,8.856795],[6.423937,1.807007,6.747758],[8.359583,7.118059,11.02494],[0.9696544,6.827619,6.968257],[-4.896592,-7.638462,9.128127],[9.083185,-3.216083,9.687489],[2.459254,-2.957395,3.974181],[-2.39681,-1.270382,2.891119],[4.541985,0.4873715,4.676233],[-9.460445,-0.9671503,9.562186],[7.590107,3.009293,8.225908],[2.638251,-3.191267,4.259642],[-9.139729,7.036569,11.57791],[6.935173,-7.476312,10.24655],[-8.131156,-2.484712,8.560928],[-8.20285,6.448049,10.48161],[-2.566108,0.5511991,2.808688],[8.33815,-9.181812,12.44309],[-7.653275,-0.697346,7.749769],[5.750815,-6.578878,8.795085],[-1.644157,2.502182,3.156607],[8.001544,-4.54412,9.256011],[9.592201,-3.788678,10.36168],[-5.459747,7.904014,9.658275],[-8.607169,8.760939,12.32223],[-6.242472,7.63276,9.910977],[3.471325,0.4339884,3.638467],[3.026366,-6.124321,6.904071],[7.894794,-1.07993,8.030817],[7.002524,-1.336339,7.198691],[-6.571298,8.696166,10.94556],[-4.171717,0.5407022,4.323839],[-7.931824,4.270163,9.06356],[2.265677,-1.447201,2.868394],[8.949616,3.969532,9.841383],[9.591448,-2.640757,9.998473],[1.899026,-1.049223,2.388968],[-4.657338,-7.414188,8.812547],[3.039388,0.2759309,3.211544],[0.7200609,-2.477669,2.767188],[5.521271,-1.360895,5.773774],[6.904535,-4.97416,8.568248],[-1.515025,1.422287,2.306122],[8.901109,-0.983125,9.010898],[1.096773,-4.102502,4.362732],[5.31308,-7.815383,9.503107],[-9.773732,3.479935,10.42285],[-1.573908,8.146759,8.357444],[-7.376445,7.544026,10.59831],[-5.656526,-9.459756,11.06722],[-6.659563,-7.868097,10.35648],[-8.211194,9.664722,12.72126],[6.065735,6.730831,9.115767],[-1.80612,6.07508,6.416281],[-1.93766,-3.934766,4.498546],[-1.043132,-8.535367,8.656824],[1.082779,-1.33746,1.990279],[-8.260476,3.73684,9.121372],[3.133033,-8.634405,9.239526],[7.918419,-5.564094,9.729363],[-0.2925778,-8.412028,8.476309],[9.66062,-1.66644,9.854167],[9.236466,-4.083993,10.14846],[2.641332,-4.880038,5.638387],[3.318348,-4.692002,5.833208],[-7.413627,-7.310647,10.4598],[4.094108,5.166577,6.667476],[2.441275,-4.066107,4.846963],[7.721624,1.030003,7.853941],[-4.33989,-8.356751,9.469421],[-6.351804,1.380281,6.576519],[-8.924212,4.922951,10.24095],[9.967822,8.887691,13.39211],[0.7842956,-8.444777,8.53987],[-5.610432,4.5746,7.307798],[0.632697,7.831387,7.920286],[9.380954,-6.912312,11.6954],[7.028822,1.140217,7.19058],[6.19697,8.67802,10.7103],[6.169463,7.095762,9.455798],[-0.1794305,-1.602341,1.897285],[1.632371,7.942024,8.169479],[4.584695,-7.443453,8.799115],[2.813052,-5.689362,6.425115],[3.173967,-3.632269,4.926199],[5.85272,6.894988,9.099187],[4.54765,2.13636,5.123003],[-7.795033,0.9679657,7.918301],[-4.262423,-6.057808,7.474309],[1.476844,-5.96887,6.229645],[6.600505,-8.905945,11.13025],[6.422215,0.5131214,6.519826],[-9.202706,5.420678,10.72723],[-5.424316,-4.608298,7.187462],[-8.290477,2.834357,8.81848],[-3.508288,6.525858,7.47629],[0.7399307,-1.676838,2.08789],[1.987607,5.903499,6.308873],[2.861574,4.889699,5.753066],[5.169264,-9.835417,11.15602],[-1.751613,-1.887963,2.762707],[-7.267958,4.759838,8.745243],[0.7858691,-6.917046,7.033002],[-1.523984,8.703543,8.892366],[7.512656,-1.580856,7.742034],[-2.081262,7.739766,8.076858],[-1.503352,-0.1419066,1.811133],[3.251731,-7.453951,8.193604],[-9.923244,-7.380457,12.40733],[6.246166,4.960128,8.038498],[4.892059,7.104119,8.68336],[-8.77663,-6.096045,10.73271],[-5.245758,-9.159153,10.60227],[7.422558,1.000876,7.556197],[3.213123,-8.24697,8.907114],[3.939478,3.001538,5.052595],[-9.516781,-9.371164,13.39357],[-2.523352,3.426829,4.371551],[-5.729075,-8.570105,10.35707],[5.348776,-5.536144,7.762621],[5.292478,2.401353,5.897187],[6.082316,8.890553,10.81834],[5.125439,3.887425,6.510161],[-3.77132,-7.313207,8.288899],[7.840828,6.23939,10.07018],[5.73671,8.823674,10.57199],[3.344836,8.766995,9.436532],[-8.570426,4.464202,9.715004],[-2.904012,5.164235,6.008545],[2.686827,6.999568,7.563927],[-8.936654,7.465593,11.68755],[8.849547,-3.689308,9.639786],[1.130079,-1.693978,2.268621],[7.602661,3.283154,8.341435],[-3.090892,4.873279,5.856831],[5.706783,1.888113,6.093631],[-6.285462,-1.66278,6.578135],[5.267685,-0.6284767,5.398471],[-0.4415352,-8.816381,8.883891],[-2.86185,6.643712,7.302678],[-9.146807,-6.444479,11.23367],[-9.585285,5.416756,11.05527],[-2.635964,6.69594,7.265254],[3.881209,-6.300143,7.466966],[-1.568742,5.421486,5.731794],[-4.357574,1.633952,4.760068],[-1.98475,2.256488,3.167171],[-6.210318,8.472669,10.55245],[4.095752,-1.295031,4.410475],[8.676723,-5.959231,10.57345],[8.079821,8.639658,11.87128],[-9.415824,9.416656,13.35407],[-0.8661263,5.349517,5.510672],[-7.702868,9.904681,12.58717],[0.7006514,-5.067911,5.21293],[-1.85542,4.551318,5.015684],[-2.379559,-1.351352,2.913495],[-6.887371,3.731286,7.896731],[5.155287,-2.94114,6.018911],[6.837196,3.92751,7.948118],[-8.317698,-4.149278,9.348829],[-7.89674,7.05005,10.63305],[6.95813,3.24177,7.741101],[0.7217754,9.695105,9.77323],[-4.6004,-7.038717,8.468011],[8.907845,-1.22187,9.046694],[-1.057764,-6.452589,6.614738],[-4.763662,0.9064254,4.95117],[5.253841,-2.097028,5.744596],[2.770725,-4.220033,5.146416],[7.682407,-7.484457,10.77202],[-8.700435,-6.033474,10.63487],[-0.9729668,-4.512551,4.723323],[-1.719562,-3.481554,4.009752],[-6.162062,8.920661,10.88803],[5.081709,6.195903,8.075456],[0.920469,-5.596223,5.758904],[3.355679,-5.736771,6.720946],[1.895612,-4.447153,4.93665],[7.032446,2.618459,7.570444],[-6.687431,-9.483097,11.64692],[-7.534453,4.344533,8.754596],[4.341943,-6.434105,7.826249],[0.9652938,8.653297,8.764209],[6.959234,7.108282,9.997931],[5.868962,3.708798,7.014264],[-6.015137,-8.088872,10.12974],[-3.977848,9.711598,10.54222],[-2.339062,2.581921,3.624573],[9.550542,-0.241664,9.605793],[2.983241,-8.792088,9.338122],[-8.929605,4.373181,9.993125],[0.06810029,-7.89494,7.958311],[-1.9054,8.059601,8.341926],[8.743557,2.420268,9.127294],[0.3870269,2.071447,2.332527],[7.717556,-0.983084,7.843923],[3.978788,8.016177,9.00499],[8.449105,8.914165,12.32273],[-7.755002,5.760336,9.711927],[-6.982034,-3.896673,8.058093],[6.329274,1.272274,6.53287],[9.046917,0.8706508,9.143563],[9.730608,-1.583254,9.909158],[5.149648,-3.549349,6.333778],[4.629429,-1.694898,5.030337],[1.033128,-4.118344,4.362122],[7.756213,1.936475,8.056599],[-6.316041,9.288143,11.27661],[9.484395,-2.358924,9.824371],[-7.664793,9.95261,12.60173],[-9.916412,-7.367805,12.39434],[-6.598246,0.04601313,6.673752],[9.580043,3.652634,10.30141],[5.942468,-2.85431,6.667834],[-7.307456,-8.399507,11.17813],[-4.165298,3.4255,5.484866],[-6.38595,-2.287893,6.856735],[1.641474,5.207185,5.550604],[-1.334195,-1.697057,2.379092],[3.236351,7.901875,8.597302],[7.572099,3.136566,8.256799],[-0.3554333,-5.286224,5.391706],[-0.2011173,-7.753931,7.820735],[2.305842,-0.4710882,2.557114],[-6.802968,-3.015779,7.508348],[2.831103,6.660454,7.305942],[3.834315,5.54005,6.811323],[6.021662,6.106284,8.634067],[9.554968,-8.989521,13.15709],[4.209889,-1.554219,4.597691],[4.471103,-5.756631,7.357279],[9.566424,-3.851477,10.361],[-6.714189,-5.285353,8.603213],[-1.960425,1.447983,2.634373],[-6.248015,5.433964,8.340603],[-9.72469,8.866408,13.19783],[-0.1106547,1.158469,1.534371],[5.554611,3.68553,6.740685],[3.05864,-2.891586,4.326262],[-6.538094,-5.412646,8.546544],[6.865629,-5.912524,9.115635],[1.420902,2.510075,3.052776],[-6.286458,2.637825,6.890405],[0.233871,-6.101477,6.187303],[1.083033,-0.9272058,1.741457],[-3.974171,-6.740233,7.888268],[-6.318327,-5.749727,8.601198],[-3.191135,-2.517413,4.185775],[0.4266985,5.520627,5.626667],[-1.521813,3.056114,3.557492],[-0.6417565,-0.8126413,1.439527],[4.367552,1.771806,4.818175],[2.593255,5.75035,6.386823],[1.868464,7.535292,7.827629],[3.450119,3.188926,4.803392],[0.8317133,-2.112475,2.480786],[6.534982,0.0442455,6.611198],[6.240709,-5.659387,8.483815],[-5.97239,8.756268,10.6462],[2.328843,5.210394,5.79411],[2.712455,3.620313,4.632934],[2.716785,-2.250359,3.666747],[-9.760301,0.2789031,9.815358],[9.66917,-6.137235,11.49602],[8.327071,0.3147519,8.392805],[-6.091465,-6.372541,8.87216],[4.106571,4.945037,6.505176],[5.536244,3.701572,6.734363],[-5.290969,-2.929846,6.130118],[5.180579,2.297593,5.754766],[0.6832566,0.08816669,1.214336],[-2.690414,9.533264,9.955976],[3.965437,6.270551,7.486288],[-7.413686,9.942178,12.44225],[5.37561,0.71338,5.514172],[-4.761615,5.871332,7.625321],[-1.961591,6.754912,7.104694],[-9.554995,0.8501869,9.644726],[4.492762,8.114991,9.329415],[-5.103177,-4.530147,6.896712],[7.859592,1.502844,8.064225],[-6.565852,-1.904069,6.909117],[5.588115,-5.661044,8.017135],[9.541453,2.329136,9.872396],[-7.996684,-7.369851,10.9207],[0.05165851,3.219379,3.371509],[3.851957,7.827186,8.780798],[-8.36538,0.9758527,8.481266],[-3.565815,4.11937,5.539336],[-9.905754,-1.005444,10.00674],[-9.818669,-7.35217,12.30694],[7.863722,-1.46149,8.06065],[3.921498,6.455255,7.618954],[-4.351451,4.40727,6.273688],[-4.320205,-1.097396,4.5682],[-1.466658,2.138135,2.778975],[8.439494,-2.970955,9.002869],[0.4265202,5.579271,5.684205],[-5.825909,-5.46436,8.049872],[-9.456888,-5.921628,11.20261],[5.461783,-3.96217,6.82128],[-4.82569,-4.642879,6.770791],[-3.788849,-7.504992,8.466421],[5.978655,-6.893345,9.179461],[0.116692,-0.7462668,1.253208],[5.378054,-7.652009,9.406206],[-2.360428,8.898854,9.260736],[9.965649,-8.621993,13.21563],[-4.170848,5.957799,7.341073],[-4.10597,7.112105,8.272909],[6.937434,5.077609,8.655062],[2.562068,-5.524598,6.171335],[-8.634955,-9.516273,12.88883],[2.099887,-1.132987,2.58712],[-3.012925,-9.312162,9.838399],[-1.709529,4.106616,4.559253],[-3.921282,-3.445087,5.31461],[2.085742,8.655048,8.958805],[1.359517,-4.407688,4.719746],[5.342777,9.518728,10.96136],[-4.940499,-5.061295,7.143195],[-0.247977,-1.805706,2.078958],[9.782456,5.382183,11.21001],[-2.331711,3.114439,4.017039],[8.966556,-5.507945,10.57055],[-3.464669,-4.609225,5.852255],[-7.351461,-6.866153,10.10881],[0.5360526,8.934681,9.006435],[4.384489,5.792258,7.333076],[6.544532,7.994675,10.38006],[7.459399,6.353863,9.849579],[4.98022,-4.943608,7.088149],[9.291154,-1.486595,9.46232],[2.139353,9.518391,9.806966],[-8.80283,-3.522274,9.533951],[8.395709,-0.4741697,8.468339],[9.731677,6.89488,11.9685],[-4.499157,-2.006398,5.026733],[-0.8001164,6.568672,6.692357],[6.836718,0.325233,6.917116],[1.652053,-1.482573,2.434605],[6.057894,-1.343887,6.28523],[-4.145637,-2.346744,4.867598],[-8.060881,-0.03882799,8.122765],[-5.71077,-4.473901,7.323161],[7.48904,-8.861048,11.64491],[5.740404,4.939793,7.638965],[-4.74698,1.924257,5.218868],[-2.823775,-7.717319,8.278328],[-8.628425,-5.341908,10.19734],[6.43164,-2.607637,7.01183],[8.430803,4.737818,9.722416],[-4.659626,7.431812,8.828587],[6.685815,-5.611357,8.785639],[0.1051155,-5.857368,5.943047],[-9.359874,0.5586219,9.429703],[4.074767,4.524537,6.170507],[-4.039064,-1.87703,4.564787],[3.148256,-2.764809,4.307631],[8.86189,-0.222809,8.920916],[-9.04057,-3.216115,9.647554],[-9.841917,0.6955872,9.917014],[4.467254,5.138409,6.881832],[-8.068205,-0.02004646,8.129965],[6.537325,8.001423,10.38072],[-6.597108,-6.702083,9.457259],[-6.373348,6.440958,9.116222],[-8.749536,-5.958732,10.63301],[5.149532,-0.5191855,5.27136],[-7.491467,-7.558377,10.68883],[0.6213144,5.101829,5.235904],[-4.063463,-7.085945,8.229359],[-2.493814,6.42942,6.968253],[7.437223,-3.112956,8.12421],[9.930758,-1.441631,10.08455],[1.246011,8.332426,8.484213],[2.162602,1.983187,3.099981],[-2.310179,2.136279,3.301608],[-2.569721,-7.256524,7.762771],[-8.60691,-2.811168,9.109422],[-5.856011,-3.09377,6.698081],[0.7385638,-3.430947,3.649229],[2.95355,3.311035,4.548232],[-9.52874,4.605996,10.63071],[-6.832178,3.556505,7.76707],[-2.355436,-5.667684,6.218579],[-8.485926,1.503675,8.675942],[2.195603,3.237589,4.037654],[-2.443788,-4.422638,5.150906],[5.173549,6.928321,8.704438],[-7.776254,-8.034103,11.22573],[8.296927,-0.5362085,8.374158],[8.904134,9.059843,12.74223],[-9.250561,9.672396,13.42118],[5.405952,0.1290019,5.499178],[-5.17282,2.242718,5.726068],[3.919319,-0.7608084,4.11581],[-2.278648,7.727407,8.118193],[9.304192,-2.816639,9.772483],[6.231879,-2.345387,6.733287],[-1.653624,-1.864072,2.685002],[5.264566,-3.480803,6.389964],[-9.706495,-7.054489,12.04084],[-2.679351,-4.590302,5.408308],[6.641344,-6.419457,9.290689],[-2.913373,-6.908608,7.564166],[3.870156,4.831626,6.270783],[-7.393209,4.330655,8.626361],[1.183056,-2.271774,2.749651],[9.08559,2.433326,9.458807],[5.684877,-6.011395,8.333948],[9.011618,-4.404053,10.07993],[-8.631409,-5.459059,10.2617],[7.876629,7.705898,11.06446],[8.418513,-9.029033,12.38527],[3.622995,6.075939,7.144447],[-6.705878,2.394955,7.190592],[4.718052,7.856261,9.218506],[6.456702,-0.2307421,6.537755],[2.685582,5.483787,6.187428],[6.417583,9.863062,11.80955],[-8.336783,-5.668727,10.13096],[-8.916729,6.150015,10.87799],[-0.3382897,-7.915481,7.985566],[5.394983,9.396948,10.88157],[-9.017365,8.893993,12.70496],[-5.133132,-7.645036,9.262592],[0.9121807,6.035666,6.185575],[2.698978,7.885425,8.394308],[9.213474,-3.386525,9.866948],[8.962705,7.412101,11.67344],[4.300461,5.755917,7.254278],[-9.070176,-1.521765,9.251155],[-4.831354,8.565942,9.885208],[8.493099,-4.161547,9.510584],[-5.254676,2.810923,6.042592],[-2.929318,-1.700692,3.531749],[-5.206558,2.557817,5.886482],[0.7493373,-7.718967,7.819461],[9.81639,-3.052752,10.32864],[2.225068,4.864552,5.441947],[0.1344084,-8.253077,8.314526],[2.339441,-0.8816288,2.692629],[-9.689688,-9.611759,13.68488],[3.323415,-8.120241,8.830821],[-7.273092,8.890711,11.53007],[-1.398099,4.331806,4.660388],[8.747486,-7.951686,11.86372],[-5.067562,2.287933,5.649321],[-5.519969,-0.3342389,5.619766],[-1.812067,-2.803627,3.484811],[-7.362425,-7.271195,10.39594],[-1.290146,-6.171833,6.384042],[-7.096188,-8.086367,10.80487],[-5.429584,2.564573,6.087481],[5.832578,-4.313067,7.322671],[-5.942492,4.803075,7.70602],[6.093147,-4.735665,7.781579],[-7.090342,-7.195148,10.15101],[8.974855,-9.789474,13.31848],[8.745869,6.83419,11.14434],[-1.178647,4.644712,4.895156],[-0.167893,8.34449,8.405873],[1.904197,5.89456,6.274696],[-6.109274,0.5557809,6.215475],[8.755855,4.906089,10.08636],[-9.880948,-1.798324,10.09292],[6.93046,-3.91654,8.023128],[-0.6973412,8.06389,8.155527],[5.670933,-4.990992,7.620334],[-2.2499,4.243398,4.905963],[9.680354,1.064406,9.789904],[-4.644562,5.451469,7.231215],[-1.301411,3.209804,3.605067],[7.548035,-8.489804,11.40393],[5.480314,1.018457,5.663135],[-2.775547,9.785541,10.22059],[4.002745,-9.348926,10.21882],[-0.6218412,0.227203,1.199295],[5.158777,-8.978812,10.40346],[6.243384,5.139518,8.148281],[-5.589895,-5.539083,7.93274],[-8.571853,-1.484201,8.756684],[8.012198,6.763923,10.53309],[-3.616441,-2.206655,4.352927],[-1.090332,-5.183389,5.390394],[-1.210217,1.705862,2.318317],[-3.317095,1.584062,3.809511],[1.391651,8.610571,8.779444],[9.417293,1.088147,9.532548],[-4.738892,7.317909,8.775471],[-3.088496,7.106491,7.812875],[-7.760248,-0.9514556,7.882051],[-8.694865,-4.343307,9.770618],[5.997252,9.162061,10.99592],[8.955977,9.724174,13.25779],[7.432267,-1.615341,7.67124],[-2.729711,7.706961,8.237024],[-4.181184,-9.063421,10.03135],[-4.154621,-5.545578,7.001022],[-1.364996,-8.247792,8.419579],[9.825677,4.709202,10.94169],[2.263019,-9.656537,9.968449],[-5.227311,4.568881,7.014232],[6.863929,-2.975228,7.547549],[5.31236,-2.96378,6.164832],[3.267597,9.184994,9.800067],[-1.843375,2.228842,3.060354],[-9.684142,6.058618,11.46689],[6.425085,8.608649,10.78845],[-1.680075,-0.01013124,1.955187],[7.469435,-7.385708,10.55183],[-0.278835,-2.467799,2.677271],[6.609545,4.503655,8.060334],[8.748904,-5.104685,10.17846],[0.1020295,-9.888182,9.939142],[-9.750089,-2.884565,10.2169],[-8.774037,9.647614,13.079],[-3.777452,-3.71931,5.394665],[-9.86069,1.210602,9.984927],[1.506357,-6.463502,6.711629],[-5.071323,8.040187,9.558395],[-5.555139,-8.524122,10.22351],[-5.0653,-1.573336,5.397467],[9.563263,0.9779875,9.665012],[0.7141263,7.939212,8.033746],[-5.553621,9.850652,11.35245],[-5.881377,-8.221259,10.15774],[9.005633,-0.3744994,9.06872],[-9.025229,-3.601507,9.768603],[6.65066,7.566576,10.12346],[0.8085372,5.029103,5.190916],[2.032837,-1.907108,2.961332],[8.789245,-0.4793704,8.858929],[-3.617323,-2.095429,4.298354],[-9.099251,-6.303992,11.11471],[0.6466826,0.6772971,1.370011],[-8.179686,-3.059276,8.790133],[-1.221609,0.6041986,1.69038],[-4.655384,4.092217,6.278443],[-6.971653,-0.1437718,7.044475],[7.96048,-4.929736,9.416556],[9.048763,6.648025,11.27281],[9.480151,-4.017835,10.34487],[-6.518622,-6.95547,9.584936],[-4.661863,7.60145,8.973016],[-7.492027,-8.984595,11.7411],[7.760184,-7.593283,10.90314],[-0.7040624,-3.416538,3.628834],[-0.1242139,-3.866629,3.995779],[8.908635,5.680697,10.61292],[3.989624,9.404256,10.26436],[2.744138,3.270545,4.384833],[2.820515,-8.217632,8.745558],[-1.45814,-2.891018,3.388828],[-2.992555,5.971439,6.753775],[-3.714221,8.616257,9.435853],[0.2560058,-1.414347,1.750976],[-4.498209,9.467051,10.52896],[-8.752057,-9.207722,12.74287],[1.334849,5.513638,5.760384],[8.461061,4.323266,9.554066],[-8.432881,9.241212,12.55044],[-3.217351,-9.346067,9.934804],[0.6229316,7.519641,7.611376],[1.628538,-1.296076,2.309101],[-9.465523,6.338847,11.43578],[7.065081,-6.036934,9.346654],[8.621999,-0.5058121,8.694521],[-4.247729,-9.543222,10.49363],[-0.9477508,0.03759045,1.378276],[8.890133,4.517643,10.02215],[2.0116,-2.082537,3.063249],[-7.505306,-9.225971,11.93517],[2.36611,2.101718,3.31899],[-3.33867,-2.447456,4.258727],[5.828794,-8.462157,10.3239],[-7.134597,5.976293,9.360478],[-8.185899,-5.263017,9.783061],[1.698678,6.850496,7.12845],[7.3473,-2.647092,7.873367],[8.634253,-5.799885,10.44935],[0.1925482,-4.227911,4.348828],[7.329391,-6.806287,10.05214],[-9.287465,-1.965493,9.545689],[-2.014692,-9.492306,9.755144],[-1.922591,2.764014,3.512283],[-1.29779,0.3498677,1.675311],[9.815413,9.297482,13.55675],[7.095695,4.883923,8.671885],[-5.638193,-1.126377,5.835919],[4.954532,6.950664,8.594132],[3.594721,7.875298,8.71449],[-0.1549462,9.639007,9.691979],[-0.8411581,-1.084186,1.697942],[2.096442,-3.98409,4.611729],[-2.665272,3.031994,4.158925],[-5.242934,-5.796334,7.879457],[-3.009352,1.360167,3.450544],[-2.700135,-1.152307,3.101377],[2.291648,-8.31062,8.678597],[3.056344,-1.038729,3.379378],[7.783367,9.582894,12.38599],[4.343237,-9.35814,10.36525],[-4.010896,8.03449,9.035502],[-5.968506,-0.4497551,6.068389],[-8.029069,3.757147,8.92088],[3.512609,2.072149,4.199074],[-1.011096,-0.5239977,1.515549],[-2.335566,9.126686,9.473715],[-5.474001,1.112651,5.67474],[1.702962,-3.95173,4.41772],[6.219316,-9.266407,11.20474],[6.044534,5.17817,8.021835],[-2.634613,7.970249,8.453759],[-6.672451,5.64056,8.794175],[4.519989,5.107052,6.892915],[-3.02192,-8.217034,8.812018],[1.33807,-2.131897,2.708398],[-4.765869,2.518663,5.482442],[3.194593,-3.837593,5.092401],[-8.492406,-8.415919,11.99786],[-0.1218517,-1.235991,1.594529],[-4.608742,1.736366,5.025482],[-2.874432,4.579149,5.49827],[-4.453241,7.003469,8.359422],[6.158055,5.206137,8.125608],[-5.731982,7.914561,9.823232],[5.960926,4.749156,7.686815],[-3.112435,-3.158636,4.545794],[0.2766018,-7.739667,7.808903],[-2.887481,8.668385,9.191216],[-6.764134,6.418889,9.378467],[-7.859723,-1.715155,8.106603],[4.126016,-7.35873,8.495583],[-3.542499,-4.532334,5.83878],[0.5339231,-0.8909119,1.441804],[-2.081234,1.796609,2.925634],[-2.842892,8.68109,9.189307],[-9.951938,4.897386,11.13667],[-0.3280921,0.220888,1.075377],[6.215904,7.040621,9.444988],[1.446348,5.096464,5.391277],[-3.708343,-7.751833,8.651169],[-0.002566022,1.479632,1.785866],[1.874298,8.665091,8.921703],[-6.506621,-7.233875,9.780852],[-1.275742,6.726372,6.91893],[-6.366201,-1.458828,6.607321],[-4.658998,4.66425,6.667945],[-9.433063,-5.543643,10.98702],[-2.501868,-0.2671383,2.707528],[-1.350676,-0.3601305,1.718727],[4.002132,7.692793,8.729039],[-4.401681,0.1279342,4.515658],[-0.2190793,-4.614385,4.726579],[3.480919,-4.216139,5.558114],[9.599639,-9.6509,13.64892],[-2.20225,-9.355509,9.663097],[3.807378,-5.496129,6.760441],[9.406713,-0.1573883,9.461027],[-8.592397,4.526348,9.763048],[-2.091454,-9.054244,9.346311],[-8.828965,-8.527289,12.31525],[8.764039,6.335141,10.86013],[-3.003377,8.59447,9.158886],[-7.862121,-1.402326,8.048569],[-6.490038,8.896433,11.05745],[-8.587867,6.218555,10.64997],[-5.865268,1.686724,6.184368],[5.694637,8.542978,10.31559],[-4.710813,-1.702457,5.107849],[3.718666,3.879225,5.465974],[-8.211199,-7.16979,10.94667],[-5.418053,-5.798284,7.998462],[5.477439,3.478985,6.565491],[-7.968633,-6.507578,10.33672],[-2.273689,-4.890684,5.485294],[-9.060187,-4.532425,10.17988],[-5.806609,3.062366,6.640391],[-6.570701,-2.363581,7.054121],[-0.6823145,-1.955002,2.299475],[8.201462,-9.440188,12.54516],[4.121156,-6.447497,7.717134],[-2.49813,-0.00123166,2.690846],[7.942626,-2.841238,8.494582],[1.612689,9.601096,9.786818],[3.402644,-1.511659,3.855269],[7.55603,8.795719,11.63865],[-5.553257,6.100834,8.310165],[3.398596,-5.177855,6.273805],[-9.544361,8.947153,13.12046],[4.647408,9.876268,10.9608],[5.288739,2.327195,5.864008],[0.719022,-6.080226,6.20372],[3.316679,-2.15897,4.081852],[-0.9522999,-0.89951,1.648027],[-5.447639,-5.487117,7.796489],[4.440183,-2.822706,5.355641],[-6.182523,-4.375196,7.63976],[-6.426242,4.640546,7.989447],[9.054606,-9.650908,13.27124],[-0.7563801,9.219437,9.304307],[-8.854386,-8.188512,12.10173],[-3.449108,-2.822445,4.567553],[5.776517,-9.365284,11.04883],[-6.185121,4.530274,7.731695],[-8.104848,8.766241,11.98063],[7.139179,1.205273,7.308937],[-6.296031,-8.536051,10.65383],[-3.497852,7.212914,8.078433],[8.575718,5.877668,10.44461],[0.7321961,3.73062,3.931111],[-3.17447,7.690309,8.379624],[-4.953062,3.430844,6.10766],[3.726999,1.441295,4.119205],[2.830591,-4.727957,5.60052],[1.253199,-8.023269,8.181891],[-3.300244,8.248303,8.940141],[2.598935,6.450672,7.026068],[6.118854,-0.8247897,6.25465],[2.535256,7.773664,8.237559],[1.750148,2.644239,3.324909],[-5.480698,7.616063,9.436232],[-2.477341,5.715315,6.308887],[-8.359075,-8.674232,12.08786],[2.058509,-1.386204,2.675634],[-0.2880091,-0.9446818,1.40548],[3.853998,-4.108811,5.721507],[3.345562,-2.07923,4.063986],[-3.488717,-7.060175,7.938338],[-9.906375,-3.078929,10.4219],[-2.838387,5.312733,6.105864],[9.849002,4.112174,10.71974],[-9.061871,-1.155482,9.189812],[9.111228,-0.08303156,9.166317],[-3.756477,1.371501,4.122152],[-1.507452,-1.886025,2.613331],[-0.58546,-1.205292,1.671973],[1.97936,-6.908501,7.255704],[2.201385,-1.234532,2.714805],[-8.938419,7.885831,11.96168],[-1.845297,4.45355,4.923336],[-4.036589,-6.875444,8.035284],[-3.328401,-6.001176,6.934866],[6.043948,1.617139,6.335964],[-5.187975,-0.5583677,5.312895],[9.902233,-8.881231,13.33906],[-5.784319,9.678064,11.31916],[2.685019,-8.34774,8.825764],[2.63303,-9.275112,9.693326],[6.277372,-1.153648,6.460364],[-2.016976,-5.513213,5.955142],[-2.421239,8.743479,9.127477],[-4.146562,2.954444,5.18871],[-9.037774,7.762835,11.95588],[-1.324682,8.197504,8.363843],[7.756094,4.20682,8.879997],[-4.411003,1.319431,4.711459],[-4.821663,-0.245895,4.930406],[3.719921,-2.870233,4.803753],[-7.615142,7.625823,10.82329],[4.9959,-2.83152,5.828938],[-2.855427,-4.817865,5.689049],[1.77477,3.524735,4.071065],[0.7634682,-1.379442,1.867015],[5.010292,-1.503087,5.325627],[0.5012153,0.01284967,1.118652],[8.257016,2.624504,8.721601],[1.435719,0.484439,1.815481],[-9.445323,6.765341,11.66122],[4.808167,0.07687366,4.911657],[6.356051,-8.039009,10.29685],[9.258588,6.240853,11.21025],[7.749983,-6.373483,10.08382],[-9.316336,-6.740025,11.54218],[-2.150416,1.224058,2.668822],[-0.6554601,-7.155839,7.255044],[-2.65005,-1.705379,3.306219],[2.650752,-3.973596,4.880158],[-1.006485,5.569822,5.747689],[-6.150459,0.3751894,6.242509],[-1.674111,-0.8397586,2.123168],[-3.943789,-9.676479,10.49703],[6.676079,1.451513,6.904847],[-8.09348,-2.54944,8.544241],[-5.580006,-2.607285,6.239744],[9.752867,1.389617,9.901992],[-8.06992,-1.837819,8.336737],[-7.41969,-9.952277,12.4539],[9.640874,8.730766,13.04503],[4.8515,-8.992544,10.26659],[-1.664575,-1.198587,2.281977],[0.1871968,1.574679,1.874742],[2.16701,3.788843,4.477864],[-0.8812708,-5.368802,5.531788],[-1.773844,9.317446,9.537364],[-8.699462,-5.168159,10.16811],[1.823223,-7.754296,8.028278],[-9.091809,-9.377187,13.09934],[-8.67146,-2.458014,9.06841],[-8.238868,7.921352,11.47287],[-5.036019,-1.109045,5.252758],[2.640137,-4.724244,5.503527],[-1.190008,-3.394661,3.73361],[-5.995418,7.623678,9.750154],[7.918768,4.151981,8.99699],[5.723397,-2.838668,6.466475],[-5.424822,-8.702673,10.30365],[-8.104388,-1.231211,8.258146],[4.144997,-1.315554,4.462251],[-0.763609,-7.416477,7.522448],[7.301703,-7.402964,10.44599],[-4.222275,-3.189601,5.385273],[-3.342825,6.512937,7.388695],[-6.236782,6.621729,9.151216],[-8.018916,-7.18769,10.81508],[1.040102,3.415309,3.70758],[-3.859879,-9.023204,9.864931],[8.19863,-2.979745,8.780457],[8.13558,1.004134,8.258084],[5.94966,5.234669,7.987504],[8.919057,-1.185498,9.052899],[7.272133,6.004961,9.483853],[-2.783334,-2.251608,3.71708],[-5.277263,-3.601145,6.466664],[3.63908,-2.274963,4.406627],[1.864789,-4.459811,4.936329],[-2.691373,6.117193,6.757481],[-5.012511,2.863964,5.858972],[1.545533,-1.426515,2.328866],[-6.161038,9.789037,11.60964],[-0.1473002,-9.674749,9.727408],[-9.10043,-2.526996,9.497554],[-6.200179,-1.494537,6.455685],[0.2768586,-6.791515,6.870323],[-7.889418,1.267584,8.05293],[-3.936018,-4.291722,5.908563],[6.210496,-2.032846,6.610804],[2.976514,2.645685,4.106006],[-3.909181,0.3753382,4.052478],[-0.3483818,-5.173019,5.280293],[3.75699,-5.004688,6.33734],[7.125709,-8.258592,10.95354],[7.008283,0.4115117,7.091218],[-6.948149,-1.918983,7.277312],[6.815001,-9.847961,12.01776],[-4.566868,-0.8811589,4.757387],[-2.187612,2.710229,3.62367],[-3.33573,2.843475,4.495825],[6.295224,9.358687,11.3232],[6.761526,2.545442,7.293662],[8.70917,7.545433,11.56647],[-7.607735,4.410984,8.850673],[8.811073,4.441425,9.917725],[-6.558119,-1.44778,6.790066],[-0.157008,-6.101405,6.184803],[-0.823268,-9.595239,9.682272],[3.031785,9.295607,9.828531],[5.298597,4.597368,7.085967],[0.9390346,3.085336,3.376549],[-7.215797,-2.252079,7.624932],[5.657963,-8.142364,9.965472],[-8.947058,-0.1616536,9.00422],[7.549945,-7.775246,10.88375],[-0.3071874,-6.258468,6.345296],[9.382704,-5.702458,11.02512],[-3.229882,-4.88838,5.943769],[7.024275,3.595139,7.953959],[-2.767973,2.992092,4.196939],[1.483588,8.040079,8.236741],[3.453118,4.525949,5.779986],[-8.9196,-3.160176,9.515566],[-5.792703,-6.936838,9.092587],[7.307197,-8.679286,11.38969],[9.564786,-6.247383,11.468],[-9.887691,0.3943029,9.945949],[8.665423,-0.02664633,8.722974],[-7.876329,2.243804,8.250528],[5.066976,-6.356473,8.190177],[-1.112815,8.163833,8.29979],[6.213939,-1.951916,6.589614],[-6.724212,0.2824129,6.804027],[-6.722588,-0.4982205,6.814793],[-1.460195,1.998938,2.669817],[3.238274,2.761164,4.37155],[-2.699486,-9.475277,9.902934],[8.953737,1.709341,9.170129],[3.679215,-5.06748,6.341607],[9.768892,2.79393,10.20967],[-1.92808,-3.840088,4.411776],[-4.570403,-1.164003,4.82115],[-7.929158,-8.693118,11.80855],[3.492252,-9.882214,10.52872],[-9.904985,0.45875,9.965901],[8.695321,-5.568851,10.37404],[3.841689,2.379962,4.628476],[5.798681,5.338654,7.945183],[1.129591,8.250124,8.386926],[-8.275185,-0.5251546,8.351913],[9.773652,-2.488147,10.13485],[-0.005044406,7.744898,7.809191],[1.490057,8.808585,8.989518],[-3.576774,6.918777,7.852566],[4.626541,-3.106527,5.661748],[-5.174949,8.846015,10.29719],[6.693248,-0.5275385,6.788068],[0.4719242,-4.972604,5.094065],[-7.243177,5.430484,9.107897],[-6.066568,9.1572,11.02985],[-7.072553,-6.128037,9.411367],[-4.507348,6.91926,8.318193],[6.87511,3.839505,7.937817],[5.224842,5.181442,7.426056],[-7.607229,2.33487,8.020072],[-6.814748,0.2450995,6.892087],[-9.334838,5.711603,10.98916],[-6.47277,-5.674035,8.66553],[8.364111,0.1291814,8.424669],[4.820153,-2.893068,5.709967],[8.507408,2.960811,9.063244],[-3.906602,9.89422,10.68443],[3.057073,-7.405908,8.074228],[-3.617354,-1.625133,4.089781],[9.408063,-6.921834,11.72277],[3.768478,-6.018914,7.171384],[2.13392,8.84993,9.158322],[3.323169,-6.733059,7.574796],[-3.851066,-5.083335,6.455308],[-0.6354638,9.982995,10.05306],[-2.691278,4.111516,5.014733],[-8.095251,-5.038843,9.58765],[1.146332,-6.596794,6.769916],[-7.868552,5.90678,9.889599],[2.39183,5.231067,5.838229],[4.609417,-2.658875,5.414457],[-8.647758,-9.548406,12.92114],[7.628523,3.240009,8.348175],[-2.48239,6.190083,6.74384],[6.336138,-9.344101,11.33397],[7.565905,4.233246,8.727159],[-3.508573,-7.707371,8.527229],[-4.697723,2.035825,5.216626],[8.59857,7.408839,11.39413],[4.80681,5.250265,7.188234],[-3.336297,8.888959,9.546961],[0.4867842,8.215775,8.290712],[-5.046217,8.070573,9.570708],[-6.626388,0.167133,6.703503],[4.465596,2.094603,5.032784],[7.138154,-3.55425,8.036537],[3.790457,0.7672266,3.994521],[7.991635,-0.7210651,8.086171],[7.747819,-1.794024,8.015436],[-2.753538,-8.350665,8.849609],[-3.225133,-5.777491,6.691853],[-1.229509,5.664421,5.881951],[6.108497,-6.421207,8.918837],[-8.343468,0.3801396,8.411775],[-9.227252,-9.006428,12.93282],[-0.397496,8.681831,8.748268],[4.404018,0.7060758,4.570987],[9.835362,-9.508087,13.71634],[-4.771787,5.149257,7.091177],[2.828335,-4.092783,5.07448],[1.008883,5.978785,6.145219],[4.263533,-2.58185,5.083666],[-8.485537,-4.598759,9.703243],[-9.48019,6.38921,11.47589],[6.712624,-1.365494,6.922708],[1.414343,-1.143638,2.075638],[-9.482386,-5.671178,11.09405],[-3.794306,-5.964787,7.139709],[-0.3386682,-2.933716,3.117914],[-3.829758,-7.527452,8.50468],[2.848239,-2.734752,4.073246],[-6.673223,7.120142,9.809604],[5.494308,-0.8725359,5.652322],[-0.8194822,7.531318,7.641485],[-3.192826,9.380002,9.958844],[8.104774,-0.05057115,8.16639],[-7.465996,0.3622955,7.541376],[-3.315508,-7.780128,8.516042],[8.902155,2.701534,9.356637],[5.451782,-7.083618,8.994419],[-9.95967,-7.174554,12.31541],[3.089015,-1.770146,3.69803],[-0.0772851,-5.052524,5.151114],[7.893746,-3.407162,8.655633],[1.33458,0.1459034,1.674034],[-8.996206,-9.785358,13.32985],[1.942913,-1.080992,2.43792],[6.074453,9.640243,11.43824],[0.9771609,8.476769,8.591302],[-5.720787,4.19653,7.165073],[3.699771,-3.064417,4.907032],[4.224403,-3.663934,5.680669],[-3.54878,-8.805574,9.546307],[-4.697133,8.407189,9.682142],[7.865429,-4.505397,9.119407],[-3.852514,7.852449,8.803568],[4.025621,0.1871599,4.152186],[7.304988,-3.252869,8.058784],[5.521527,1.889843,5.921044],[8.551996,-2.215318,8.890685],[2.73181,-4.204573,5.112848],[-5.109887,-8.29265,9.791781],[8.989524,2.96101,9.517306],[-6.209202,-4.516525,7.742945],[6.486707,6.613221,9.317299],[9.139816,-6.624868,11.33248],[-1.292964,5.224843,5.474554],[-0.1063254,-1.293281,1.638255],[2.259216,5.133118,5.696749],[1.303775,-0.7563099,1.808821],[2.842672,-8.042193,8.588227],[-7.415227,-0.3837342,7.492186],[5.003982,9.205214,10.52501],[-2.767699,2.86475,4.106939],[7.453221,5.480054,9.304918],[8.823992,8.389367,12.21656],[7.943764,-4.775422,9.322449],[6.796726,-5.774668,8.974535],[6.471035,4.971551,8.221352],[1.348951,-4.429848,4.737428],[4.810726,-7.823994,9.238937],[8.849734,9.221502,12.82006],[-1.6973,-8.576768,8.800101],[-4.549253,6.337482,7.865074],[-7.077578,0.8252182,7.195352],[-3.153342,-5.547021,6.458561],[0.3686684,-9.894412,9.951649],[1.689115,-1.891754,2.726141],[-4.762505,-8.633714,9.910726],[6.191169,-1.029447,6.35534],[5.478487,9.02305,10.60327],[-2.381128,-6.679209,7.161118],[0.8694938,1.342964,1.886683],[-1.479618,-2.535592,3.10137],[8.662861,-9.634588,12.99502],[-3.809545,-2.652849,4.748709],[4.292932,5.916048,7.377594],[-0.3986225,2.363383,2.597014],[-8.442391,5.360811,10.05049],[-1.901922,-8.864787,9.1215],[-0.3156001,7.236241,7.311825],[1.705611,-9.692336,9.891941],[-9.581105,8.171157,12.63192],[-9.992146,5.67351,11.53394],[5.47544,7.196645,9.09792],[7.156055,9.793629,12.17063],[5.616815,0.2763953,5.71183],[-7.906033,-0.5925217,7.991023],[3.280551,-3.067881,4.601511],[-2.436169,-3.429229,4.323718],[3.692118,-7.024133,7.998136],[6.384368,4.077352,7.641005],[-6.133347,-0.8226309,6.268545],[1.435963,-1.669516,2.418527],[-4.903309,-1.983239,5.382906],[-0.8198828,7.672681,7.78089],[-2.927693,-8.992625,9.509926],[3.293067,3.216767,4.710826],[6.500223,6.322356,9.122778],[-6.473027,-5.493163,8.548387],[-0.8160008,0.1336438,1.297582],[-2.711504,-8.749492,9.214437],[-9.991629,-4.805561,11.13221],[-7.331072,0.7562594,7.43751],[-1.695784,-7.786423,8.031442],[-8.593862,-8.500889,12.12929],[-2.587658,8.771084,9.199342],[-4.683978,-3.516619,5.941907],[5.394187,5.572039,7.819519],[-6.480431,0.5368044,6.579068],[1.141664,-4.216708,4.48152],[6.993794,4.528875,8.391893],[-6.397434,-9.912034,11.83958],[-9.476209,3.356431,10.10268],[4.301738,1.12315,4.557018],[-4.082791,-6.275594,7.553295],[1.525703,9.275139,9.452828],[-6.826807,-2.895879,7.482741],[-1.175551,3.201461,3.55405],[-3.677925,-8.575896,9.384729],[8.652877,9.630168,12.98508],[-3.274067,9.531051,10.12721],[-7.722453,-4.9268,9.214643],[3.781347,0.03422666,3.91149],[6.208839,3.547518,7.220427],[-3.652349,8.78375,9.565246],[-0.446119,0.5883811,1.243067],[7.720441,-4.937887,9.21889],[6.157681,-6.616125,9.093412],[4.341005,4.609493,6.410285],[5.130385,-1.593324,5.464388],[-7.642524,-4.509466,8.929919],[-4.883055,7.941676,9.376269],[0.7517663,7.243032,7.350283],[-1.604522,-5.258277,5.587841],[-2.188421,7.331568,7.716286],[9.119783,-4.371568,10.16273],[3.083964,-2.940983,4.377239],[3.437561,-3.895181,5.290488],[-8.163863,-4.848783,9.547741],[0.6851652,1.242979,1.736217],[1.009438,5.304302,5.491319],[9.957926,7.699821,12.62725],[1.001472,-3.677784,3.94069],[4.918794,-8.881055,10.20136],[-5.788497,8.994181,10.74253],[1.381,9.492838,9.644747],[7.120334,6.971377,10.01495],[4.899911,8.824444,10.14297],[0.05064856,7.516002,7.582404],[-4.266585,0.3669005,4.397541],[3.204297,9.932715,10.48458],[1.231031,-8.738987,8.881742],[-4.907528,3.52068,6.122011],[9.359031,-5.059325,10.68589],[6.2562,0.8751062,6.395769],[6.283878,4.647843,7.879694],[-2.570055,7.643559,8.125834],[-6.432033,5.068046,8.249615],[4.086849,-5.514896,6.9366],[2.950351,0.4111879,3.142236],[6.983173,-2.28801,7.416178],[8.102537,-4.798543,9.469801],[-7.845864,-3.315463,8.576122],[-7.096749,-5.40293,8.975271],[-7.437901,-8.455383,11.30557],[-4.561204,7.474215,8.812971],[-4.390179,-0.214652,4.507743],[9.229441,-8.587643,12.64635],[8.967897,-7.956067,12.03005],[0.264659,-2.943572,3.120042],[4.731253,-8.441511,9.728508],[1.193913,9.197316,9.328239],[-3.745926,6.616585,7.668844],[9.116789,-4.306604,10.13226],[-1.245416,-7.508409,7.67641],[-7.361577,5.392913,9.180213],[8.425856,-6.2451,10.53548],[-4.848618,0.822446,5.018517],[6.95961,-6.469832,9.554837],[6.434489,-5.877359,8.771887],[-3.37372,5.054373,6.158626],[-9.001662,-1.35842,9.158341],[-2.044367,-7.748425,8.075736],[-3.162754,4.331297,5.455561],[1.664819,8.809249,9.020781],[-6.46683,9.140206,11.24114],[7.956764,-7.237038,10.80207],[9.184902,3.552776,9.898719],[4.5175,-8.944674,10.0705],[-9.897186,7.741169,12.60476],[1.373851,-1.034538,1.989406],[-0.2811024,3.781803,3.921869],[7.023057,0.139555,7.095267],[0.1885484,-3.968369,4.096768],[-2.852903,9.584125,10.0496],[-5.129827,-8.878697,10.30274],[-5.098757,-5.021812,7.226058],[3.284099,9.989229,10.56267],[-9.069891,-9.509806,13.1795],[-6.721789,6.90901,9.691072],[-2.986058,7.610407,8.236191],[-2.463647,-4.868319,5.54708],[-4.547173,0.807107,4.725273],[-6.679747,-7.014564,9.737717],[-6.247161,9.64542,11.53521],[0.7707648,6.82713,6.942894],[-3.468883,2.358726,4.312394],[-6.761377,9.147274,11.41879],[7.079865,-3.720404,8.060143],[0.007705009,6.944863,7.016494],[-4.967925,-7.147223,8.761454],[-3.046752,-5.145854,6.06321],[-1.808607,-8.653306,8.896671],[5.466548,-5.95747,8.147061],[6.8967,1.710771,7.175737],[0.2790313,-5.683641,5.777684],[-8.680164,0.01532558,8.737591],[-0.06896258,-0.7183118,1.233178],[-8.829638,-5.288618,10.34079],[-9.513818,7.024251,11.86814],[-7.732433,-4.752549,9.131114],[-9.193524,1.798593,9.421031],[-2.032303,-2.823985,3.620103],[-0.8309576,1.60697,2.067086],[2.288458,-8.655633,9.008719],[5.747072,-9.592149,11.22667],[-0.8378792,5.315582,5.473341],[9.575475,6.467237,11.59805],[-8.789529,2.082773,9.088111],[5.219288,1.649848,5.564437],[-8.625629,-3.695544,9.437082],[6.040983,-1.133327,6.22719],[2.741824,6.415592,7.048221],[-8.485114,-3.696175,9.309074],[2.640933,5.277783,5.985776],[2.469792,9.504356,9.870799],[0.7398889,7.092432,7.200696],[-5.334754,-0.7022697,5.472913],[9.467709,-2.043286,9.737173],[-9.436759,2.404919,9.789589],[-6.294471,-7.161283,9.586675],[-6.481522,1.655963,6.764047],[-0.7341953,-5.560243,5.696959],[-0.5354129,-7.750636,7.8332],[7.79922,-8.954752,11.91702],[3.490864,1.368623,3.880626],[-8.162423,-7.03473,10.82186],[-7.297086,9.249571,11.82379],[2.439311,-4.893268,5.558265],[8.761016,7.196026,11.38148],[9.773023,8.033479,12.6905],[9.510495,-8.119328,12.54484],[-1.172419,-5.649761,5.856139],[-9.607784,-8.076482,12.59123],[3.292017,-2.411379,4.201443],[-0.7444665,-1.690541,2.100514],[-8.275818,-7.36447,11.12315],[8.012709,-8.94308,12.04916],[8.291071,2.01175,8.590052],[-6.172001,-5.438268,8.286637],[-4.105119,-7.221293,8.366545],[-0.3722795,-6.79686,6.880109],[-8.328864,2.214564,8.676074],[-1.544134,3.353922,3.82533],[1.986037,-5.987418,6.386981],[4.704031,0.1426945,4.811265],[1.148187,9.655889,9.775199],[0.7464197,3.047385,3.292977],[7.829722,-2.622285,8.317508],[-5.57087,0.4018295,5.674158],[-6.511107,-1.569089,6.771747],[0.7765609,7.785348,7.887629],[-6.186577,1.716758,6.497769],[5.077928,2.393133,5.701968],[-5.973591,4.570589,7.587758],[-6.099952,1.132228,6.284214],[-2.800604,5.364227,6.133377],[0.3003373,-6.732472,6.812957],[3.639672,-5.984831,7.075692],[7.724808,-8.724241,11.69551],[4.027124,7.783534,8.820495],[-0.7645735,-3.550792,3.767319],[-4.085245,0.4775698,4.232883],[0.9195137,-7.841184,7.957994],[5.101462,9.192356,10.56051],[0.5780753,5.777665,5.891993],[-2.561065,8.471248,8.906239],[-0.6121992,-0.8301218,1.436625],[7.901653,0.8812317,8.013282],[-4.573014,-4.305457,6.359985],[9.919281,9.452029,13.73801],[1.312574,-7.84256,8.014275],[-6.543033,3.614633,7.541675],[7.086831,-8.570868,11.16615],[0.6757597,-1.361085,1.819121],[6.717368,-1.815259,7.029808],[-4.229144,-1.43209,4.575646],[4.523955,9.938449,10.96535],[-7.598751,-5.265216,9.298575],[1.18102,-5.004004,5.23783],[0.763927,-1.085117,1.661645],[1.211819,-4.242982,4.524534],[4.850451,2.118425,5.38652],[1.962341,4.892951,5.365795],[-6.607779,9.919451,11.9607],[1.221171,9.297963,9.43098],[0.5271804,-7.09757,7.187031],[-9.678009,-2.703565,10.09817],[-9.09907,9.817823,13.42322],[6.683273,-9.351516,11.53763],[-5.314906,-6.976854,8.827498],[0.3049673,-6.528482,6.611663],[0.9843282,-9.253682,9.359463],[-7.665352,-4.591029,8.990838],[6.307841,3.232333,7.157991],[-9.233768,-1.716256,9.444998],[-5.466522,-0.1988414,5.560791],[-0.184316,-5.052709,5.154012],[-2.569852,8.400002,8.84105],[6.169427,-1.700746,6.47722],[-0.6498017,3.140564,3.359373],[6.435998,4.323072,7.817354],[6.803027,8.09551,10.6216],[6.352477,-9.034736,11.08965],[-1.895594,3.952019,4.495746],[-5.363084,-4.058688,6.799677],[-7.084423,-4.74158,8.583219],[1.712365,-8.505388,8.733488],[5.27159,-6.968347,8.794744],[-6.209102,-9.944238,11.76609],[6.63827,3.555931,7.596794],[-7.33247,-8.580679,11.33107],[9.697966,6.691155,11.82464],[-0.1269867,-5.709984,5.79828],[4.673099,-8.732858,9.954931],[-9.997308,7.459841,12.51381],[-9.341907,2.134025,9.634588],[5.259044,2.48706,5.902796],[-1.301864,6.548519,6.751144],[-0.9596233,-6.113801,6.268927],[-2.595451,1.381586,3.105663],[3.437685,0.4200046,3.60473],[0.1655611,-4.346909,4.463521],[-4.053598,-2.621016,4.929643],[4.945438,-0.5255925,5.07283],[-9.11295,5.326234,10.60258],[3.455571,-4.370979,5.660957],[-0.7705561,2.192904,2.530333],[3.451024,-6.308059,7.259557],[2.798734,-6.770386,7.393986],[9.738323,6.916812,11.98654],[4.443854,5.760645,7.3439],[-8.23856,8.830943,12.11856],[8.855104,1.932678,9.118559],[4.072318,-3.769139,5.638278],[-4.8446,-6.335761,8.038159],[3.934173,-0.7636639,4.130485],[-2.36902,6.872863,7.338154],[-9.166941,6.138456,11.07761],[-4.464229,-3.87862,5.997752],[0.8175203,-7.306541,7.41983],[5.720234,1.76786,6.070125],[-7.647174,-6.145559,9.861398],[2.924499,-7.046849,7.694853],[-0.08725037,0.8386303,1.308019],[-1.645242,-8.13881,8.363436],[7.236998,-7.617834,10.55488],[-7.849757,-5.640213,9.717546],[-9.239086,-6.287152,11.22003],[-7.492447,-8.612645,11.45925],[1.577802,-4.511866,4.883276],[-8.308249,-4.074402,9.307405],[6.281799,8.378735,10.5197],[-6.022283,-1.762617,6.354109],[3.423368,6.074414,7.044002],[-2.788911,-3.515841,4.597734],[1.612056,-8.943851,9.142823],[-2.370699,1.976176,3.244301],[-7.317078,-0.7115886,7.419298],[-4.053985,9.608095,10.47618],[9.33849,-2.678935,9.766478],[6.633677,3.807328,7.713716],[1.348844,-3.198411,3.61237],[-9.053053,-0.3889849,9.116418],[5.00709,0.2592806,5.112551],[6.957533,9.824132,12.07977],[-3.684072,5.841953,6.978596],[-0.2403985,-6.789228,6.866689],[6.263057,-1.620401,6.546112],[5.962756,2.202227,6.434614],[3.120473,-3.160677,4.552718],[-5.356354,5.635052,7.838644],[6.988067,2.059327,7.353497],[1.159495,-5.385854,5.599272],[-3.733474,7.839953,8.740921],[-8.07836,-9.467192,12.4855],[-1.187566,9.108995,9.240353],[-0.3347791,-0.3319912,1.105575],[-1.532578,9.994389,10.16054],[1.791123,-9.268029,9.492338],[0.2999002,-1.795169,2.076673],[-2.146001,2.822607,3.684078],[-0.7468649,5.325447,5.469752],[-4.151899,-8.718081,9.707893],[5.516128,-9.466512,11.00193],[1.082967,-8.782898,8.905735],[5.62564,8.740252,10.44221],[-7.96995,9.348073,12.32504],[-3.078927,-8.181225,8.798422],[-9.753015,5.781401,11.38182],[-5.694021,-7.643212,9.583348],[-1.213917,-5.71652,5.928929],[-8.153647,3.40878,8.893916],[0.3021806,1.682104,1.980097],[-4.732285,-0.1283279,4.83849],[9.423412,9.444643,13.37916],[0.3246038,-0.7324023,1.28132],[-7.070089,-0.2012226,7.143294],[1.244143,-6.826352,7.01049],[-7.451756,0.2967781,7.524409],[0.1533905,5.279818,5.375873],[-8.853349,-1.81026,9.091689],[6.714858,9.938375,12.0358],[1.911941,-3.718906,4.299509],[6.890465,-1.261867,7.076074],[-7.603061,8.142638,11.18522],[4.589641,2.88841,5.514319],[-5.827999,-9.055609,10.81525],[3.561408,4.197209,5.594657],[-0.2642165,-2.52885,2.732196],[9.186699,-0.1297752,9.241877],[-4.694407,3.880183,6.171975],[-6.886331,5.598645,8.931202],[6.524509,9.721381,11.75051],[-2.283896,-0.04342385,2.493605],[-1.150239,7.654446,7.804717],[6.439602,-9.009646,11.11945],[-2.681404,-2.405801,3.738691],[-0.0925548,6.323695,6.402943],[-1.396544,-8.747549,8.914592],[2.429743,-9.046222,9.420073],[7.490023,7.358797,10.54762],[1.558008,8.090616,8.299726],[-2.316818,-7.004644,7.445313],[-1.992821,6.97081,7.31871],[-8.845498,-7.167365,11.42865],[5.467544,1.315682,5.711834],[5.086253,-9.110887,10.48228],[1.89333,-1.11359,2.413458],[-8.265115,4.677238,9.549275],[-2.081071,9.787651,10.05629],[-3.15471,5.953996,6.811921],[8.770033,1.882806,9.025433],[-5.30592,0.933628,5.479457],[9.350909,-3.563342,10.05668],[1.247216,2.133069,2.66562],[5.019186,9.93304,11.17397],[4.732299,-4.283632,6.460971],[6.000791,-0.5384012,6.107321],[9.572927,-3.306265,10.17705],[5.745532,-4.120096,7.140471],[1.614259,-3.60828,4.07744],[2.263816,0.5671461,2.539],[-9.512125,-3.452347,10.16854],[-5.818967,3.038677,6.640326],[0.2256227,-5.920393,6.008491],[6.063024,-0.2007696,6.148217],[0.9010613,6.393887,6.534041],[-9.705169,-8.847168,13.17052],[8.509007,-1.111046,8.639307],[-0.4304795,6.827981,6.914235],[5.593645,-4.386387,7.178388],[2.423759,-7.813915,8.24208],[8.499418,6.178888,10.55551],[-1.129305,-4.586727,4.828394],[-5.961918,9.796649,11.51168],[-7.770928,8.702407,11.70979],[7.784481,5.324234,9.483966],[1.85839,6.88756,7.203617],[3.193527,0.4917679,3.382373],[-5.713337,-3.215566,6.631899],[-9.179932,4.847199,10.42912],[-7.297921,5.082871,8.949594],[-3.967735,6.736189,7.881571],[9.549225,-4.695562,10.68812],[-0.5834654,3.644828,3.824291],[7.575147,-2.321303,7.985693],[0.4854755,-8.296776,8.370913],[5.076977,-8.943995,10.33299],[0.1283744,8.624952,8.68368],[9.15706,-5.885106,10.93097],[-9.639448,-4.320027,10.61045],[-1.899052,-2.159673,3.044764],[0.7542663,-4.671247,4.836266],[-9.941556,-4.640815,11.01688],[-6.824981,-1.708938,7.106393],[1.539199,2.127271,2.8097],[-8.997989,9.704909,13.27212],[9.203363,7.57479,11.96158],[-9.45731,-8.032806,12.44856],[-5.893833,-2.751062,6.5807],[-3.483333,3.995053,5.393891],[-3.832588,5.008461,6.385406],[-9.673155,-0.8475894,9.761574],[-0.4823169,-4.889526,5.01399],[-7.898379,-5.849718,9.879454],[0.4478993,-3.592537,3.75592],[5.756533,0.6527089,5.87909],[-9.558059,-9.265389,13.3493],[-5.337088,7.261429,9.067131],[3.022364,2.362781,3.96452],[-2.803579,-1.145029,3.189224],[5.231949,3.637268,6.450039],[-7.605154,3.36605,8.376674],[-5.593075,0.965758,5.763261],[2.584531,-9.319409,9.722715],[8.798114,-8.216676,12.07976],[-0.02406261,9.64856,9.700273],[-6.234048,-5.630423,8.459611],[-8.994363,5.549534,10.61583],[-3.934919,-2.230843,4.632521],[-8.099369,-4.055617,9.113057],[-3.073812,-0.894674,3.353917],[5.006932,1.272436,5.261982],[-4.978136,-1.991595,5.454199],[-6.235122,-9.192406,11.15245],[1.468492,-5.148221,5.446159],[-8.09447,-1.517537,8.295985],[-2.659615,4.49363,5.316603],[-0.07068485,-5.309044,5.402865],[-5.053023,6.671087,8.428312],[8.817739,8.855864,12.5371],[-3.101916,2.781476,4.284681],[-7.93964,-8.167303,11.43428],[-0.9734198,9.840242,9.938708],[-9.483519,4.263031,10.4456],[9.651171,3.495205,10.31317],[-4.554325,2.838437,5.458809],[-9.797226,-5.483606,11.27189],[9.550074,-7.252088,12.03315],[-6.543791,-7.173251,9.760981],[9.731215,-6.716716,11.86637],[8.721338,-1.530789,8.910952],[3.13209,-5.897849,6.752378],[0.3623556,-1.335218,1.707076],[-8.048385,7.362542,10.9537],[-5.016875,-0.2769558,5.12306],[-5.434774,-7.253398,9.118583],[4.445865,1.895331,4.935382],[-9.350907,4.069662,10.24703],[9.854008,4.845221,11.02622],[9.244578,7.733581,12.09423],[0.5292944,-5.874404,5.982372],[3.187846,-8.898858,9.505369],[3.320816,-4.07546,5.351373],[9.593668,2.810463,10.04675],[7.93499,-7.807923,11.17711],[3.654922,8.270161,9.096924],[-2.659665,3.888083,4.815704],[-7.971523,-5.835155,9.929461],[-0.2112339,-7.139644,7.212429],[9.417537,5.670951,11.03855],[2.743155,6.495494,7.121541],[-0.8812202,6.617611,6.750506],[-4.057792,-5.683664,7.054765],[2.872254,7.472265,8.067502],[3.256895,-8.876364,9.507745],[2.592221,-6.163302,6.760614],[3.790209,-2.989434,4.929747],[0.8941025,-8.761087,8.863186],[-2.419285,7.437423,7.884681],[9.391314,8.886251,12.96774],[9.130524,1.489739,9.305148],[-8.796959,3.278603,9.441172],[-6.78279,8.60393,11.00154],[-1.603462,2.294063,2.972173],[-7.788879,-9.776871,12.54009],[8.747035,-5.921741,10.61026],[2.953391,-7.877439,8.472105],[-0.4669046,5.02558,5.145333],[-2.517196,-5.847114,6.443991],[-3.318007,-8.668057,9.335115],[8.285123,-1.783755,8.533759],[-2.359167,1.907107,3.194171],[-3.420571,-0.5741838,3.609708],[6.117384,-4.891331,7.896043],[-3.323358,-6.814673,7.647515],[-2.772804,6.798331,7.409841],[6.757722,-5.501614,8.771235],[1.985038,9.786057,10.0353],[-4.235989,-0.4532931,4.375966],[-1.379007,9.358368,9.512136],[7.865448,6.428217,10.20722],[-9.545401,1.731019,9.752492],[-3.885469,-1.711608,4.361935],[0.08586908,-5.063686,5.162199],[-2.206097,3.23831,4.043948],[0.5418423,3.972797,4.132397],[2.241916,5.372716,5.906967],[-0.3032296,-6.6926,6.773688],[-0.8791466,-1.77591,2.21963],[8.921905,0.6323236,9.000011],[4.798172,-8.797273,10.07047],[1.758845,9.439041,9.653446],[6.645394,7.756927,10.2631],[-4.502148,4.034481,6.12751],[6.597911,-8.129627,10.51776],[6.904244,3.074252,7.623622],[4.74366,8.657477,9.922411],[-4.986433,-1.591686,5.328975],[2.504702,3.565574,4.470666],[8.868769,6.908651,11.28648],[-2.612949,-7.628025,8.124917],[-1.770458,2.79946,3.459985],[3.176291,9.138727,9.726519],[5.820123,7.270753,9.366839],[-7.045428,-6.689942,9.766953],[5.595902,2.445398,6.188222],[-3.557668,-1.729897,4.080386],[-2.964814,-8.73433,9.277858],[8.511715,9.360169,12.69102],[1.49765,-2.233779,2.869272],[5.886727,3.404234,6.873309],[0.1698984,-5.785687,5.873929],[6.011628,-9.377918,11.18414],[9.787042,-5.612803,11.32651],[4.164275,3.74994,5.692384],[0.1928364,7.256134,7.327255],[6.299042,-0.1753928,6.380336],[2.708809,-1.630122,3.315862],[-3.33498,-7.812174,8.552903],[5.92335,-9.181023,10.97166],[7.153953,-3.932617,8.224629],[-2.555628,-9.798676,10.17572],[-8.616889,-8.320122,12.01978],[-3.388689,1.604545,3.880435],[0.1664538,-6.354692,6.435047],[4.092364,2.89539,5.111823],[4.581472,7.261191,8.643771],[-8.448915,6.270019,10.56869],[-9.050406,2.198597,9.36716],[7.100004,-4.617297,8.528158],[1.532559,-6.326555,6.585897],[0.3410364,6.253156,6.341787],[-3.034579,2.828206,4.267015],[2.510525,7.781698,8.23757],[0.6797224,-9.960518,10.03364],[6.942195,-3.740421,7.948888],[1.313348,-3.060961,3.477695],[7.712651,-0.9719962,7.837714],[-8.169613,8.937236,12.14976],[-7.42193,-9.612754,12.18565],[-6.911743,-7.098249,9.957777],[-3.42612,-1.973306,4.078264],[8.868762,2.612778,9.299545],[1.404593,9.097745,9.259689],[0.391002,5.885264,5.982409],[-3.756169,9.166971,9.957016],[-6.264457,-8.932402,10.95588],[-9.242496,-7.695589,12.06838],[-7.186439,7.193016,10.21687],[2.961105,2.994728,4.328572],[6.609887,-3.500645,7.546199],[4.264806,-2.959003,5.286234],[-4.640473,2.084695,5.184587],[4.858999,-0.2942559,4.969553],[2.290388,7.652032,8.049812],[-4.097424,9.383471,10.28778],[6.456778,7.940283,10.2829],[-1.944724,3.691883,4.290915],[-2.573021,8.341157,8.786088],[4.895772,7.386437,8.91785],[5.019363,7.118575,8.767447],[-6.335101,-9.684591,11.61571],[-3.062484,-8.666995,9.246384],[-9.49537,-6.591806,11.60233],[-0.3846353,-4.60069,4.723801],[1.142475,0.2174905,1.533803],[7.398667,9.363044,11.97526],[-6.941974,5.502146,8.914292],[-4.058337,4.865401,6.414221],[7.360017,-0.8161455,7.472345],[3.92758,-0.9051989,4.152742],[-3.012098,-1.90766,3.702958],[-4.877238,-9.336536,10.58104],[4.706364,2.817341,5.575596],[8.601871,-8.671347,12.25498],[-1.037801,1.899425,2.384291],[-3.150777,-2.664414,4.245762],[-0.1070798,9.052677,9.108372],[-3.413745,-9.096563,9.767349],[-4.246274,-6.023652,7.437421],[9.991524,-1.478807,10.14975],[-3.474917,5.728939,6.774643],[-2.42769,-8.389099,8.790373],[-7.912705,-6.337269,10.18685],[-3.54007,0.2667708,3.68826],[4.744718,-7.252232,8.723945],[0.6642042,-2.329581,2.620709],[-6.583406,2.233697,7.023578],[1.251512,-2.532171,2.99636],[6.695128,0.7445791,6.810223],[-9.863089,3.59788,10.54634],[-6.863056,9.171932,11.49895],[4.911162,6.69242,8.3611],[8.534121,-6.095148,10.5348],[-3.327712,-1.788103,3.90781],[1.797744,6.78762,7.092508],[-1.130452,9.377501,9.49818],[-4.917552,3.640016,6.199357],[5.255393,-3.457038,6.369479],[0.5862154,-0.9307338,1.486578],[-4.845455,6.719533,8.344493],[-1.665051,1.971965,2.767858],[1.668697,-6.764209,7.0384],[8.989677,-9.278499,12.95781],[-9.911304,8.471678,13.07682],[-9.185816,-3.801953,9.991699],[3.65679,6.145722,7.220942],[-3.864568,-4.718306,6.180396],[3.444493,4.613634,5.843813],[-0.5813892,-3.966104,4.131343],[-3.647323,-7.949504,8.80327],[-8.965711,-0.7316296,9.050925],[-7.944056,-2.004438,8.253835],[7.634384,1.499874,7.844325],[0.8717657,1.819062,2.251436],[1.332752,-7.18169,7.372442],[-2.394544,2.190699,3.396028],[-0.5636058,2.444658,2.700742],[2.439452,0.4832851,2.68039],[0.4899778,-8.655725,8.727065],[-3.847348,-7.299609,8.311821],[-5.427393,1.415544,5.6974],[5.742879,-7.221491,9.280657],[6.85184,-5.342892,8.746097],[-7.209309,2.649543,7.745594],[9.364814,-0.8147883,9.453234],[-2.11233,-2.692232,3.565116],[-6.783095,-9.762004,11.92925],[7.335895,-3.350999,8.12678],[7.280794,-4.636864,8.689676],[3.285004,9.469914,10.07326],[-0.1920985,-2.857993,3.033979],[9.973012,9.134368,13.56089],[-5.389814,6.13601,8.228045],[-0.1878814,6.142237,6.225944],[-1.76835,3.948017,4.440034],[7.774332,0.4561277,7.851642],[2.789953,-1.354859,3.258754],[-2.84064,-2.637895,4.003464],[-0.7772512,-3.396895,3.625329],[-9.319514,8.001427,12.3238],[5.554537,6.716581,8.772989],[-4.458842,8.110513,9.309226],[7.119703,-5.845469,9.26605],[-3.230946,7.751057,8.456825],[8.420546,-8.126653,11.74513],[0.6405681,-6.031434,6.147237],[6.174857,3.04193,6.955731],[-9.355068,1.001159,9.461481],[-7.499075,-3.906469,8.514495],[2.178648,2.056628,3.158517],[-1.123759,-8.675184,8.804638],[7.717688,0.9430032,7.83913],[-4.613144,-7.127548,8.548862],[-9.622925,6.954903,11.91517],[-5.8213,6.20454,8.566437],[-0.1227936,-7.073433,7.144826],[-6.642548,9.062164,11.28035],[7.850344,4.545984,9.126548],[-9.061531,-1.978853,9.328837],[7.64973,3.433202,8.444243],[-2.459876,2.20541,3.451786],[-0.7093611,5.221506,5.363517],[-4.494261,-7.562414,8.853728],[-7.474944,-8.238286,11.16889],[2.419461,-3.063378,4.02965],[1.423836,5.595123,5.859413],[5.855127,7.504602,9.570871],[5.090905,-9.420571,10.75474],[4.925076,-1.607487,5.2764],[-6.618388,-3.037699,7.350556],[-8.779198,-9.589725,13.03983],[7.780836,6.638401,10.27666],[-6.239818,2.143626,6.673114],[4.4191,-8.213828,9.380588],[-6.259062,3.395531,7.190653],[-4.803408,-7.567867,9.019165],[3.610719,-5.29274,6.484627],[5.953036,5.581099,8.22115],[8.361509,2.938143,8.918941],[-3.490412,3.505813,5.047147],[9.895079,5.529818,11.37943],[-1.349115,4.531281,4.832456],[3.557745,3.142116,4.850819],[3.330611,4.211916,5.461978],[4.550493,-7.957925,9.221473],[-0.7844535,5.214839,5.367487],[-6.271703,-5.109148,8.150929],[-1.272887,6.883322,7.071094],[5.594701,-5.622087,7.994282],[-2.620367,-4.740732,5.508254],[-4.122531,5.4164,6.879873],[-3.502323,2.652334,4.505679],[-9.312833,5.852195,11.04432],[-5.552279,-8.776051,10.43297],[-9.879997,-5.662083,11.43125],[-8.816473,-5.062438,10.2156],[-3.204977,7.870223,8.556418],[0.1417008,3.800217,3.932141],[0.06335647,3.322359,3.470171],[-4.098778,8.733688,9.699346],[3.975079,9.814874,10.6364],[-6.032164,3.494934,7.042838],[-1.838814,8.486774,8.741086],[-5.199818,9.895197,11.22288],[-5.039094,5.492867,7.520908],[-9.092926,1.972321,9.357957],[-3.910808,-1.131015,4.19209],[7.351281,-2.762437,7.916589],[-7.429469,0.8545226,7.545012],[5.32816,2.733668,6.071427],[8.336119,-2.076338,8.648819],[-4.949964,5.482555,7.453895],[-0.5404385,-4.246474,4.395977],[-5.546232,1.800582,5.916316],[-4.096731,5.247416,6.731908],[2.277301,-5.843974,6.351231],[-0.08136605,6.326296,6.40536],[-0.05548909,-0.272862,1.038043],[3.449883,-7.058935,7.920244],[-8.220737,-6.501094,10.52828],[9.720039,3.000892,10.22177],[-4.23629,-6.534389,7.851395],[4.321012,2.120915,4.916241],[-7.352894,6.86477,10.10891],[-0.7313371,2.796905,3.059008],[0.03783916,-7.310223,7.3784],[-1.94786,6.92364,7.261608],[3.261847,4.689376,5.799129],[-3.068762,-9.041254,9.600082],[2.135581,5.654165,6.126197],[-0.195546,-2.116759,2.349235],[-6.775344,-0.2905046,6.854902],[-9.096338,2.709102,9.543721],[-8.068222,-1.248489,8.225262],[-0.2166997,7.904926,7.970873],[-3.523556,0.4324432,3.68815],[-1.788757,-6.723916,7.029274],[1.688889,-6.234468,6.536125],[-1.888813,-1.475228,2.596904],[6.561693,-5.217651,8.44273],[-3.822665,5.732094,6.962017],[-7.527582,4.65997,8.909534],[-6.254177,-4.682565,7.87662],[6.097456,2.095912,6.524708],[7.425459,-5.158817,9.096748],[5.240229,-7.262443,9.011274],[-3.688795,-1.371752,4.060654],[6.094103,7.276927,9.5442],[5.572519,4.011515,6.938676],[-4.536556,2.54199,5.295475],[3.830873,-5.836015,7.05228],[-3.490408,-0.6514847,3.688818],[0.6689489,2.583772,2.850152],[-0.5533091,-4.075671,4.232877],[2.892832,-8.678652,9.202579],[-7.840051,9.622859,12.45254],[-7.77427,-7.468882,10.82698],[-6.73363,-6.20404,9.210423],[1.490287,-5.859417,6.128109],[4.832642,4.855544,6.923203],[2.096459,-0.6882451,2.422565],[-9.093729,2.450173,9.470969],[5.009709,-4.339595,6.70293],[-5.625989,2.443484,6.214689],[-5.524122,3.950332,6.864477],[-0.557893,-2.588303,2.830293],[-2.81859,-8.517387,9.0272],[-9.861232,-5.007495,11.1049],[-3.442376,0.4933105,3.618467],[-6.51197,5.732262,8.732959],[-2.567187,-8.674474,9.10148],[-9.838678,-4.573789,10.89583],[-8.304737,2.442508,8.714041],[8.686109,-4.089959,9.652784],[-4.071528,-1.877307,4.59365],[-5.66749,3.705119,6.844585],[2.739939,-2.204227,3.655939],[-7.402947,9.09726,11.77131],[-3.331964,9.658257,10.26567],[-2.715995,-8.826148,9.288569],[7.619859,-1.003143,7.75039],[2.637368,-8.233363,8.703102],[-8.255214,6.466034,10.53367],[-7.845716,-7.342427,10.79196],[4.256182,-3.063219,5.338389],[-0.4452148,-4.077673,4.222042],[-0.2529548,-8.993378,9.052339],[-0.09896415,6.56864,6.64506],[-3.694687,7.597183,8.506932],[8.682541,-6.771453,11.05618],[4.992578,-7.588313,9.138289],[-7.178981,-1.576249,7.417704],[4.78589,5.519712,7.373734],[-7.086592,-6.263816,9.510793],[8.873513,-5.174444,10.32057],[-2.848265,-4.973245,5.817713],[-7.047772,-1.007545,7.189314],[-0.6243341,7.820093,7.908454],[6.051969,-2.638283,6.67734],[7.438481,7.134314,10.35516],[9.50623,9.146602,13.22984],[9.97924,-4.014916,10.803],[-5.946162,-0.5519149,6.05487],[2.147895,9.975302,10.25281],[0.141526,-1.999357,2.239968],[0.4501656,-8.724562,8.793215],[-9.540568,6.827865,11.77464],[-8.984278,5.37029,10.51462],[5.436249,-3.372957,6.475311],[1.826314,2.444986,3.211445],[9.689237,4.694834,10.81308],[1.876511,-4.409919,4.895782],[7.092808,-2.311072,7.526551],[9.882794,-6.344429,11.78649],[-4.095236,-4.329861,6.043067],[-3.304761,-3.615643,4.999432],[-4.882585,-8.802597,10.1156],[7.636311,-6.466364,10.0562],[-5.914627,6.736696,9.020304],[1.154734,9.60603,9.726727],[-6.188161,-9.244102,11.16901],[-7.238706,1.375,7.43569],[-5.929667,0.5634254,6.039735],[2.143885,-8.820924,9.13263],[9.765995,-1.195538,9.889588],[-3.693328,-7.949773,8.822673],[-9.812684,4.73161,10.93969],[-6.768584,2.994825,7.468782],[4.632979,8.856627,10.04512],[4.497307,6.872924,8.274229],[3.718074,-0.5565043,3.890214],[0.05082548,-2.010869,2.24637],[7.919585,8.864242,11.92873],[-8.379777,6.701325,10.77629],[-1.497805,7.458305,7.672661],[5.45614,-3.320751,6.465048],[5.295333,-2.567695,5.969389],[5.270307,0.7327378,5.414151],[1.573168,-0.2496579,1.880741],[9.907186,4.390967,10.88269],[4.808665,4.791992,6.861957],[-4.971204,-9.200284,10.50515],[-6.289574,-8.584452,10.68885],[4.519662,7.105654,8.480429],[-2.112538,5.293636,5.786657],[-0.4803076,-5.043202,5.163776],[-1.850006,9.829858,10.05229],[-7.175415,-8.742315,11.35406],[-8.878262,3.980387,9.780951],[0.3798491,4.312988,4.443664],[-7.756751,-6.411272,10.11294],[1.832332,4.953068,5.374972],[-8.392366,-0.2724867,8.456125],[-5.938877,-2.669058,6.587422],[6.633727,-2.830252,7.281254],[-4.232026,-8.653081,9.68431],[-5.738872,4.356779,7.27435],[7.488274,-3.466636,8.312149],[5.672551,9.58562,11.18311],[-4.29986,1.006204,4.52783],[-9.891646,-8.813272,13.28602],[5.180959,8.622006,10.10848],[-0.9182124,5.507638,5.672494],[-6.730693,-0.9770041,6.874356],[1.943328,-6.803111,7.145547],[7.279415,8.154059,10.97627],[9.459602,2.842313,9.927881],[-0.8713915,-3.467912,3.712915],[-1.161102,-2.629477,3.043404],[8.588593,-4.692749,9.837979],[5.127622,5.849792,7.842995],[9.655338,6.721909,11.80718],[-3.034741,-0.6952713,3.270024],[-5.020337,-9.267941,10.58766],[-5.876033,-3.129963,6.732343],[2.182881,4.71429,5.29051],[8.415937,8.34199,11.89188],[1.718335,9.986978,10.18295],[-4.330139,3.054002,5.392313],[-4.404943,-7.833227,9.042288],[6.680939,6.717304,9.526653],[9.290026,-0.9768758,9.394619],[-2.294182,2.286829,3.390112],[-1.25054,-2.082019,2.626529],[-6.616131,2.891484,7.289298],[-5.362985,5.490157,7.73973],[3.664916,9.774693,10.48695],[6.047369,-7.060947,9.350276],[4.80045,-0.6764542,4.94994],[4.717824,6.514107,8.105026],[8.077683,9.202474,12.28554],[-7.436915,7.063358,10.30528],[6.764772,-4.243815,8.048113],[-2.312011,6.492326,6.963885],[3.466763,5.979393,6.983666],[4.704898,9.057532,10.25548],[0.0009140465,-2.245543,2.458143],[6.038418,6.449228,8.891291],[5.908764,6.322896,8.711631],[5.351612,3.37442,6.40519],[-2.423498,-1.105522,2.845263],[8.913103,-7.881099,11.93965],[4.329096,0.3101319,4.453903],[-4.257115,9.459668,10.42153],[-4.509526,0.5972896,4.657529],[-2.921969,-7.095994,7.73893],[2.448617,-2.469134,3.618335],[8.069372,8.819777,11.99597],[-3.830343,4.18327,5.759451],[7.19094,4.887028,8.751723],[-9.12351,7.836963,12.06882],[5.020235,6.85226,8.553142],[-6.438981,7.45382,9.9005],[-1.828428,2.24466,3.062947],[-8.200793,8.773265,12.05086],[3.611061,6.420071,7.433511],[3.391154,-6.302773,7.226677],[1.024138,3.735057,3.999939],[-6.199077,-5.500251,8.347533],[-6.50361,0.01360957,6.580056],[-8.59259,1.730311,8.821938],[-1.86308,-4.570709,5.036115],[8.346536,3.361562,9.053439],[5.894311,-9.162116,10.94017],[2.356128,7.656633,8.073126],[-7.344551,-2.259099,7.748932],[9.428171,8.666076,12.84489],[6.465826,-1.155618,6.643972],[-2.730127,-5.38705,6.121593],[8.898476,-7.796278,11.87286],[3.97973,-7.857392,8.864359],[-0.5253742,-6.519621,6.616757],[5.266454,-0.7857148,5.41783],[-2.245651,-2.443093,3.465783],[6.851363,-7.737328,10.38303],[9.907389,-3.301506,10.49077],[-5.17101,-5.980442,7.969004],[-7.363061,2.540963,7.853099],[7.536949,-8.093553,11.10456],[-6.090154,-3.405526,7.048942],[4.208602,-4.093874,5.955849],[9.239588,-2.31151,9.576694],[3.704305,-0.04474995,3.837171],[3.746574,-6.394134,7.478086],[7.091271,-0.8615729,7.213073],[-1.377693,2.850806,3.320412],[7.856383,-6.388833,10.17546],[7.206998,-8.714416,11.35262],[1.107793,1.538144,2.14315],[7.594915,-2.443495,8.040733],[-9.439507,3.776589,10.21601],[6.987806,4.10987,8.16826],[-1.100139,-6.552429,6.718976],[7.901615,9.072757,12.07271],[-1.047432,-1.258963,1.918881],[-6.637045,-7.177428,9.826792],[-4.274283,2.32119,4.965623],[8.637917,3.074792,9.223228],[-2.338374,9.54439,9.877417],[6.477477,9.641287,11.65814],[-0.1228338,0.2732408,1.04391],[2.039651,3.488674,4.163055],[8.574807,9.511909,12.84538],[-2.475487,-0.1553291,2.674352],[5.537508,-3.076454,6.413156],[-7.411102,1.918473,7.720425],[3.441199,7.849268,8.628607],[1.622842,-2.142086,2.867429],[8.733369,6.310763,10.82116],[3.50747,-6.237132,7.225244],[7.791735,-4.564802,9.085623],[1.006766,-5.694939,5.869064],[6.58714,-5.939056,8.925402],[-7.626366,-8.630653,11.56069],[-1.09759,-0.7745461,1.674702],[-7.263681,3.792223,8.254818],[5.252061,3.331281,6.299332],[0.299731,0.8084162,1.320369],[-6.907551,-2.877182,7.549334],[9.049226,4.672171,10.23316],[-1.755653,-3.012262,3.627126],[5.368212,-1.276995,5.607888],[-3.790822,-1.044921,4.057363],[-2.534649,-9.813009,10.18428],[-0.8320971,-1.256115,1.808372],[-7.367116,1.811843,7.652267],[5.910938,8.51913,10.41704],[-9.294811,0.9440248,9.395994],[1.397336,6.107323,6.344441],[5.0289,-2.886958,5.884247],[-1.481095,3.1196,3.595212],[-3.374809,-7.506688,8.29094],[2.097844,-7.741967,8.083254],[4.176259,-3.052999,5.26896],[4.399667,-9.204826,10.25114],[8.230752,-7.325108,11.06357],[-6.059896,-9.819654,11.58223],[-0.08060724,8.285413,8.345931],[7.63002,-5.777504,9.622721],[1.726767,-6.552061,6.849176],[-9.752542,3.579413,10.43668],[1.337365,-5.146305,5.410453],[5.459431,9.382244,10.901],[8.37096,9.954946,13.04507],[3.340443,3.444793,4.901546],[-9.363785,4.336634,10.36759],[2.395782,9.295922,9.651628],[8.823208,2.752994,9.296664],[-0.6886192,-1.502559,1.931808],[9.010512,4.962944,10.33538],[6.086441,8.263551,10.31169],[-3.601287,6.955719,7.896284],[-1.613504,7.272275,7.515942],[7.673546,-6.556836,10.14275],[6.126033,3.768699,7.261637],[-5.540016,4.993939,7.525371],[-8.064345,-4.729247,9.402098],[-8.790692,-3.463645,9.501216],[7.904502,4.737026,9.269336],[7.258043,0.6505223,7.355431],[-2.991899,6.975965,7.656078],[-3.612012,-4.467862,5.831674],[-3.342825,7.693839,8.448054],[-5.722999,-5.292297,7.858825],[-4.313746,-5.644767,7.174385],[-9.754655,-4.797706,10.91656],[4.753129,6.855307,8.401635],[-3.543037,-8.214815,9.002016],[0.1988629,-7.237789,7.30925],[7.66881,-5.22968,9.335962],[7.965008,-1.604094,8.186237],[-2.211025,7.419683,7.806429],[-2.410018,0.536045,2.663744],[-9.794659,-0.7047302,9.870764],[8.355259,-0.2084167,8.417469],[0.2095964,7.204988,7.277072],[0.1958178,-0.4305041,1.1062],[7.030733,-8.234029,10.87338],[-4.370622,6.987303,8.302093],[-6.548283,-8.296862,10.61687],[-5.439164,-4.701154,7.258467],[-8.803233,2.5812,9.228191],[-7.174647,-1.497169,7.397099],[3.830369,1.143819,4.120685],[-7.39859,-2.694618,7.93726],[7.388648,8.247266,11.11798],[7.533119,-6.555408,10.036],[8.33352,5.744198,10.17071],[-9.673656,9.621207,13.68018],[6.974899,6.813864,9.801936],[-9.47474,3.111123,10.02246],[6.185612,0.2417774,6.270586],[4.286698,-3.347763,5.530217],[-7.966925,-3.256915,8.664836],[-3.247722,-7.914569,8.613252],[-6.634429,-6.054555,9.037328],[4.694387,3.311934,5.831481],[0.3181615,7.082694,7.160012],[3.778829,1.983446,4.383333],[7.563684,8.095333,11.12402],[2.512037,8.625567,9.0394],[-9.027985,5.416164,10.57541],[4.634072,-0.1929671,4.744667],[0.237661,-4.320173,4.440763],[1.252222,-6.198136,6.401949],[3.445012,-3.257416,4.8455],[6.231672,9.99636,11.82205],[9.957458,8.451795,13.099],[6.66563,-7.513552,10.09376],[-2.48316,4.460933,5.2025],[-1.735369,6.345592,6.654175],[-9.512291,5.511035,11.0388],[-1.202129,5.488572,5.706973],[-5.842095,1.741488,6.177609],[-4.080179,-5.380664,6.826375],[-0.2171224,-7.257626,7.329412],[-5.521269,8.873476,10.49871],[8.522131,8.971009,12.41393],[1.814372,-3.893864,4.410683],[-1.166887,8.80185,8.934999],[-9.198527,-3.906793,10.0437],[2.567712,5.525129,6.174156],[-7.811181,-9.271245,12.16431],[-7.034051,1.465047,7.254257],[3.563511,6.845798,7.782259],[3.21813,3.677566,4.98807],[-2.570184,-8.177276,8.629813],[2.986188,9.331241,9.848318],[1.720468,-2.783075,3.421333],[-2.101199,-1.584903,2.815484],[-4.37086,-2.036557,4.92463],[-7.236173,-0.4550656,7.319104],[4.334809,-4.564884,6.374068],[8.28124,2.345389,8.66486],[-9.477623,8.484712,12.75992],[3.979211,-6.800369,7.942237],[6.604368,-9.990359,12.01769],[4.006103,2.186614,4.672273],[-0.7665653,-7.7161,7.8183],[6.860217,3.595786,7.809754],[-1.352559,3.3278,3.728763],[2.974854,0.2344117,3.147174],[-2.062276,-0.9550993,2.482981],[-1.539129,-1.807192,2.575823],[-6.782157,2.18609,7.195599],[8.49743,8.446326,12.02276],[3.538452,-3.436829,5.033134],[0.2872576,9.083922,9.143312],[-3.517957,-2.46332,4.409531],[5.980503,-5.553136,8.222149],[-5.996469,8.854427,10.74051],[-7.571315,0.2768982,7.642087],[-8.971998,-6.223009,10.96461],[-5.032347,-7.083431,8.746399],[-6.650301,-4.314278,7.989962],[-8.827759,-1.17574,8.961679],[1.012233,-1.338721,1.953661],[-0.4069908,-5.858973,5.957617],[-9.168438,-4.964858,10.47426],[-0.3618004,1.819035,2.10708],[-2.19035,-0.6046414,2.482584],[-4.321205,2.421154,5.053197],[8.664856,-1.635404,8.874361],[-1.486179,-0.05252942,1.792062],[-6.598129,4.658865,8.138816],[-8.626231,0.8010147,8.720865],[-4.006539,-2.399019,4.775735],[-6.460368,5.013347,8.238325],[3.115916,-4.339634,5.435196],[-2.317019,9.590232,9.916709],[5.79641,6.02563,8.420605],[4.088683,-8.206788,9.223269],[3.379387,-8.217355,8.941206],[-6.845811,-2.212047,7.26349],[9.010644,0.160674,9.067389],[-4.520114,-1.108102,4.760181],[-3.960474,-2.2276,4.652693],[-5.202273,4.951246,7.251102],[0.9739152,0.2286387,1.414492],[-2.613561,9.248532,9.662611],[-7.819681,7.10535,10.61289],[5.719466,8.959301,10.67621],[1.032918,-9.273552,9.384332],[2.673075,-7.191541,7.737157],[-5.098046,-4.597415,6.937312],[-3.266222,-5.875201,6.796042],[-7.416749,8.669814,11.45311],[1.999683,-4.949436,5.43099],[-4.360598,-1.251312,4.645492],[8.892723,5.431631,10.4682],[8.976098,3.03131,9.526761],[0.8242227,-7.296085,7.410277],[7.854018,1.31125,8.02527],[9.96048,0.761255,10.03946],[8.999065,-5.700532,10.6995],[9.961438,-5.296532,11.32623],[-2.231737,0.04797886,2.446007],[-5.31444,-0.9808491,5.495938],[-6.34016,9.736681,11.66193],[-0.1405704,-4.84863,4.952673],[-2.65626,5.21047,5.933356],[-5.752221,-5.460871,7.99432],[6.285916,-3.348216,7.191891],[1.958533,-9.935605,10.17605],[4.082737,3.940379,5.761539],[-0.1660168,-9.729791,9.782454],[7.83813,1.232879,7.997266],[4.175943,-5.310148,6.829068],[8.842064,-1.426971,9.012122],[5.377401,5.005413,7.414216],[1.427061,1.175881,2.102189],[5.354255,-4.58602,7.120367],[-7.088688,1.750088,7.369689],[-6.117498,4.855997,7.874293],[4.321116,-2.662959,5.173336],[-9.133011,-1.57402,9.32145],[-9.950322,-1.086755,10.05932],[-0.9085729,-0.4616627,1.427809],[-1.765643,-8.661813,8.896319],[2.495375,7.075105,7.56862],[4.304368,-8.708488,9.765517],[3.291614,-4.252625,5.469875],[-0.06350743,-5.132362,5.229261],[-5.756594,-5.998531,8.373813],[-3.446975,-1.202322,3.78513],[1.259475,-5.305825,5.544191],[3.248164,-7.686001,8.403879],[9.28395,-4.663032,10.43722],[6.422599,5.757902,8.683502],[9.968272,4.660773,11.0494],[-6.312305,-7.660166,9.976138],[-9.961235,8.549431,13.16507],[-8.831702,-5.049996,10.22259],[-8.098369,-6.473126,10.41561],[7.935949,8.021037,11.32768],[-1.054764,7.606422,7.744041],[7.319246,-0.6679088,7.417376],[-0.6228406,9.146373,9.221934],[8.213723,7.600577,11.23539],[4.806979,3.595523,6.085625],[5.305471,-8.058913,9.700212],[8.176721,6.692789,10.61377],[4.6984,4.270964,6.42776],[8.633535,-2.12911,8.948243],[9.910367,-7.401221,12.40941],[-5.992874,-9.716436,11.45965],[-3.516543,9.720257,10.38506],[9.964641,6.477415,11.9269],[-3.843204,-2.173472,4.527051],[-3.870059,6.384424,7.532478],[-6.396796,-5.749531,8.658875],[8.310374,-9.662182,12.78359],[0.003175242,-7.472034,7.538654],[6.388819,2.682924,7.001077],[-7.805254,-1.691582,8.048816],[2.018576,-3.048574,3.790574],[-2.626921,-9.608878,10.01156],[-1.326632,6.882082,7.07976],[8.595325,-8.313768,11.99993],[6.828905,-2.515126,7.345733],[-7.01678,-7.279289,10.15988],[-0.947649,2.422272,2.786654],[-0.1471628,3.325254,3.475482],[3.847013,4.484198,5.99229],[-7.254772,2.746905,7.821586],[-1.578596,-7.978799,8.194705],[-4.118902,-0.02146784,4.23861],[6.829428,3.29057,7.646498],[-2.414286,6.551425,7.053365],[-0.7022897,-3.593711,3.795784],[3.884145,4.122781,5.751861],[2.669055,2.522107,3.805901],[-2.391649,-6.242077,6.758957],[-8.533544,-6.480155,10.76168],[4.761231,-7.98952,9.354238],[-0.4931737,-4.955318,5.079213],[-7.708428,2.381491,8.12966],[-1.329093,-8.366303,8.530037],[7.919778,-8.208394,11.44992],[-1.685743,-0.1149085,1.963399],[0.1251292,9.953254,10.00414],[-0.4837623,-3.924105,4.078312],[3.546972,-7.121738,8.018738],[9.136032,2.180532,9.44573],[0.5201251,-8.735956,8.808374],[5.674231,8.875386,10.58156],[-4.301173,8.83085,9.873398],[-0.8669236,-9.371727,9.464714],[-1.78464,-9.987376,10.19473],[7.143572,-9.672669,12.06612],[8.567078,3.208896,9.202817],[-9.730177,-6.735032,11.8759],[3.775859,0.05804045,3.906467],[8.746965,7.849939,11.79538],[-5.368512,2.414229,5.970714],[-1.038313,8.572183,8.69255],[-1.849439,-3.773161,4.319394],[3.583935,-7.010862,7.93705],[2.345286,6.449784,6.935422],[2.945334,-2.586705,4.045496],[-9.28803,4.758542,10.48386],[8.347279,8.847899,12.20502],[-8.825763,-8.560483,12.33596],[-0.511731,9.758428,9.82287],[-3.826276,7.594543,8.562563],[-2.297189,3.976526,4.699982],[-2.379452,-6.231309,6.744702],[-4.746875,2.527446,5.469992],[-6.675619,-6.247814,9.197776],[-6.456285,6.233543,9.029987],[-1.124454,1.377658,2.040181],[-9.371878,1.302352,9.514631],[-6.300137,3.942327,7.498911],[3.36234,-5.924472,6.885107],[0.02877293,8.209517,8.270247],[1.879692,-2.314492,3.144855],[-2.019574,-9.314815,9.583551],[4.050168,9.377419,10.26352],[0.06567475,-5.29344,5.387469],[1.142947,6.533425,6.707605],[6.454525,-8.207373,10.48913],[1.978135,8.494253,8.778687],[-6.413297,6.861104,9.444846],[-3.777156,0.9050322,4.010734],[0.1367473,-1.598056,1.890101],[4.929979,1.355541,5.209816],[-5.935758,3.706128,7.068848],[-6.840974,-9.309042,11.59557],[7.55067,-1.696861,7.803329],[9.236504,2.459469,9.610514],[-4.360991,-6.862153,8.19191],[-6.749981,9.751724,11.90203],[2.887198,-1.655367,3.475076],[4.555217,1.696464,4.96266],[6.938653,4.30997,8.229261],[-0.6749889,6.712042,6.819613],[-4.206535,2.998337,5.261651],[-0.7753767,-2.29432,2.620136],[-1.097536,1.321962,1.988006],[5.450697,7.635018,9.434172],[1.204,8.813084,8.95098],[0.2611371,-3.57574,3.722111],[3.34717,5.233631,6.292411],[5.807158,-7.103186,9.229211],[-2.057605,-2.960222,3.74121],[0.9492258,8.585353,8.695361],[1.487402,5.987214,6.249727],[9.780889,-2.428483,10.12735],[-3.856296,8.626837,9.502281],[-3.303499,2.641207,4.346157],[-3.849544,-5.088033,6.458101],[4.510579,6.078043,7.634653],[4.733953,4.385271,6.530001],[4.093221,-7.403226,8.518346],[-4.765079,7.905339,9.284416],[6.881997,8.986346,11.36294],[3.546852,0.5798545,3.730468],[0.6422952,-5.587445,5.71245],[9.756817,-8.002405,12.65836],[-3.078542,-7.112895,7.814775],[-9.759953,7.928243,12.61403],[8.917852,5.695872,10.62878],[-9.141186,6.993979,11.55323],[5.649253,-2.744061,6.359554],[6.009523,-9.114164,10.96277],[-5.701095,0.0003382005,5.788133],[1.162601,-5.123074,5.347666],[-3.417685,-1.639606,3.920316],[-5.920545,-1.206758,6.124469],[-6.055091,-6.870179,9.212137],[1.59142,-1.282266,2.275263],[1.511324,0.5587375,1.896388],[2.25783,2.509166,3.52047],[-3.014926,7.432176,8.082513],[-1.363776,8.54944,8.715092],[9.053772,-3.016051,9.595173],[-0.5845368,4.075113,4.236535],[4.873462,-3.898479,6.320504],[-5.198618,-2.989377,6.079639],[6.555161,8.182425,10.53196],[6.738505,5.376369,8.678294],[-1.587449,-2.443388,3.080607],[-9.210304,-0.8920859,9.307283],[4.244548,-6.98519,8.234627],[9.634701,-0.8139758,9.720597],[6.92641,-7.002208,9.899801],[-3.480829,-6.72266,7.63612],[-8.131789,0.2806815,8.197851],[-3.936451,1.292557,4.2622],[-3.306552,1.339463,3.705057],[0.5845514,-5.390221,5.513273],[1.986386,4.623194,5.130269],[-0.3905732,-1.537562,1.875272],[-0.4433313,-8.900284,8.967251],[-9.896967,8.224089,12.9068],[-5.539051,-2.244497,6.059608],[-6.192106,-6.598693,9.104116],[7.799154,1.498992,8.00461],[0.7315909,9.500031,9.580491],[-2.69218,-9.468432,9.894395],[0.2602515,-4.05347,4.183103],[5.091085,-1.426486,5.380893],[-2.066404,5.749497,6.190859],[-4.727398,3.31195,5.858097],[5.071301,7.835826,9.387134],[6.385142,-2.269845,6.849981],[-6.931363,7.424224,10.20602],[7.994401,7.560511,11.04861],[-4.704292,-5.908688,7.618593],[4.958176,5.433969,7.423713],[1.035075,8.064593,8.192012],[9.435819,-6.943178,11.75765],[1.484686,-7.530641,7.740468],[-7.137606,1.135475,7.296213],[1.133024,1.750322,2.312438],[7.621129,6.338502,9.962842],[8.642584,-2.936575,9.182468],[-5.191862,-5.779867,7.833408],[-8.512844,-7.440871,11.35055],[3.024447,-9.49743,10.01741],[-9.909455,-5.079688,11.18036],[5.817839,-3.282792,6.754552],[-1.720743,-8.72197,8.946156],[-0.1134006,-1.58053,1.873749],[3.926838,1.501001,4.321233],[-4.81981,-5.807367,7.612889],[3.21061,-6.803831,7.589475],[4.84507,-1.634395,5.210178],[6.407008,9.788041,11.74119],[-4.412595,-0.002801558,4.524489],[-7.183946,0.01592841,7.253229],[-1.874779,7.037568,7.351337],[9.419631,8.117165,12.47469],[-1.328788,4.756385,5.038737],[-8.366064,-8.125026,11.705],[1.858195,4.742432,5.190717],[-0.2389697,-4.995475,5.100184],[5.537696,5.293663,7.725862],[-5.994616,3.943197,7.2446],[-3.669861,-2.530679,4.568613],[-7.666212,8.158738,11.23992],[9.434498,-8.508505,12.7438],[9.4807,3.053211,10.01028],[5.149949,6.711606,8.518663],[0.8383255,-6.732842,6.858131],[-3.834799,-8.891646,9.734838],[-3.210676,-2.783138,4.365123],[5.958632,-7.724224,9.806576],[6.09935,-4.630532,7.722946],[4.031431,1.79892,4.526428],[-8.418694,9.856373,13.00086],[8.432146,-9.852248,13.00646],[-3.17346,9.039262,9.632191],[7.702759,-0.1876267,7.769665],[2.725812,-5.132061,5.896448],[9.259423,1.211481,9.39173],[2.226624,0.5350823,2.498834],[6.708278,7.881425,10.39797],[-7.579759,-9.813981,12.44054],[6.090719,-0.7107639,6.213055],[9.842042,-6.904179,12.06373],[6.639154,5.115967,8.44106],[1.943455,-9.899028,10.13744],[-7.16399,9.496553,11.93764],[-6.431551,6.76273,9.386126],[-4.810117,-6.545177,8.183921],[-3.343397,5.969669,6.914858],[3.148456,-4.370961,5.478876],[8.59661,6.060768,10.56573],[-3.821444,8.991271,9.820712],[-0.2248116,0.4333252,1.112794],[-8.312072,-2.366546,8.700062],[2.515895,-0.7276482,2.803426],[-1.785973,5.690784,6.047704],[-8.623596,4.222074,9.653617],[0.8864586,-9.720638,9.812065],[-2.833818,-9.506208,9.969881],[5.229098,9.661953,11.03163],[-1.203379,-1.720033,2.325217],[5.600996,1.508218,5.886075],[1.548139,0.5629035,1.927069],[2.492042,-0.261804,2.697928],[-8.79223,1.579307,8.988744],[7.561719,9.694558,12.33548],[5.772565,-1.269017,5.994407],[-0.796079,3.052852,3.309629],[-9.654892,-5.881423,11.34937],[-6.114015,-4.458129,7.632567],[-7.078677,-2.212794,7.48359],[2.98707,2.182968,3.832484],[-4.574138,1.844025,5.032213],[-5.818757,8.280632,10.1699],[9.239558,-9.492204,13.28425],[0.1273028,4.950349,5.051946],[1.439345,0.1029962,1.755654],[-9.031464,-2.038116,9.312424],[-1.606379,-1.149148,2.213819],[3.142701,-9.549229,10.10269],[-8.505375,-4.104911,9.496931],[-7.185192,-3.811416,8.194747],[-3.714039,-2.151504,4.407159],[-6.428707,-5.342975,8.418768],[-6.815545,-5.210689,8.637299],[3.961151,7.913762,8.906085],[0.2867004,0.6015444,1.201687],[2.743747,-7.070983,7.650291],[6.172889,2.129431,6.605984],[1.86723,1.465592,2.575754],[5.668011,6.405434,8.611384],[7.047507,9.039029,11.50528],[-4.166893,-7.385161,8.538361],[-3.214375,0.7783158,3.455139],[-3.6436,4.677432,6.012836],[-8.715343,5.527618,10.36879],[2.664623,9.147576,9.580103],[4.784992,1.801154,5.209636],[-0.4562527,7.237981,7.320965],[-8.938846,-8.782934,12.57151],[5.684935,0.2493311,5.777599],[8.634586,0.6135774,8.713929],[3.724366,-5.156191,6.438727],[2.584017,8.706952,9.137185],[0.7252603,4.985926,5.136678],[9.829041,-8.184972,12.8298],[8.173757,-2.371691,8.569435],[-0.2424591,-5.357921,5.455832],[-0.637948,-7.898927,7.987492],[5.171413,-3.71479,6.4454],[-6.583494,-5.908904,8.90267],[1.222323,-8.732979,8.874626],[6.176297,4.035224,7.445111],[9.043328,-6.020391,10.90994],[4.343039,-1.314122,4.646386],[-7.406632,-0.2129447,7.476867],[-5.930724,-6.058285,8.536762],[6.479351,-5.116382,8.31621],[-2.201322,0.8423689,2.560352],[0.3672808,-3.526938,3.684316],[-4.28589,-0.1351238,4.40308],[-1.449916,-4.212266,4.565681],[-1.971353,-5.75412,6.1641],[-2.354739,-9.764953,10.09451],[-7.132049,-6.765631,9.881289],[6.375224,-2.129933,6.795594],[7.556833,4.175402,8.691358],[2.60592,-3.249162,4.283442],[5.326232,-1.714092,5.683912],[-4.390086,-5.343,6.987167],[2.037257,-6.180394,6.583896],[8.631883,3.818434,9.491567],[5.822118,-5.575985,8.123341],[-1.215822,-0.2257919,1.590347],[-6.489497,-3.996876,7.68691],[0.4339459,-9.077522,9.142741],[-5.691982,-2.500682,6.296989],[4.318553,6.510588,7.876399],[-0.9427246,1.41735,1.974236],[3.144178,-2.922055,4.407296],[-5.373308,7.819366,9.540174],[3.378784,-9.442305,10.07836],[8.230807,-1.441455,8.415698],[-1.70247,-1.889684,2.733004],[-6.675728,1.263037,6.867358],[-2.193156,0.159548,2.415655],[3.134127,-1.415264,3.581302],[-1.589675,8.230441,8.441992],[7.216843,-1.173036,7.379623],[6.115482,6.648402,9.088474],[-4.718324,-8.224687,9.534572],[9.796232,-2.069963,10.06235],[7.644711,-9.458405,12.20258],[-4.778167,-9.182371,10.39937],[3.114856,-6.393063,7.181475],[-7.412675,-8.695173,11.46969],[2.659851,6.391689,6.994891],[5.090868,5.302477,7.418436],[-8.303514,-4.045014,9.290343],[-5.781781,6.563475,8.803874],[-2.282662,8.432644,8.793181],[-8.718794,-6.687311,11.03347],[-1.797253,-6.033904,6.374803],[4.670062,-4.464409,6.537616],[3.211351,8.692034,9.320098],[-6.210711,6.632183,9.141049],[2.930666,8.983961,9.50265],[-7.116538,7.187198,10.16371],[-8.724978,4.838112,10.02659],[9.900552,9.94648,14.06959],[-6.598035,9.708784,11.78111],[-6.978902,0.5033729,7.06813],[-7.714168,-1.409159,7.905322],[-5.779727,8.725163,10.5135],[8.561168,-7.696928,11.55579],[9.05477,1.715871,9.27001],[-3.672422,-0.7436646,3.878108],[-1.860377,-1.166564,2.412856],[-9.331932,-1.991064,9.594233],[5.232635,-7.602649,9.283359],[4.402773,0.8718396,4.598317],[9.289529,7.433735,11.93967],[-8.841322,2.582144,9.264796],[-2.423716,0.3212551,2.641515],[-7.504774,6.786659,10.16761],[4.908409,-8.590297,9.944128],[5.678101,-7.476893,9.441649],[8.971113,0.361403,9.033907],[4.612993,-6.309518,7.879703],[-2.801332,-8.532589,9.03618],[1.859999,-9.764791,9.990534],[-0.1223905,-2.152766,2.376842],[2.302488,9.828669,10.14417],[0.5685661,-6.249756,6.35474],[-1.978482,7.366503,7.692838],[-1.869934,8.447145,8.709243],[3.526625,-7.297998,8.166876],[8.810927,7.269131,11.46616],[-3.277145,0.2397907,3.434702],[-4.59306,-4.094045,6.233571],[-5.124409,-8.153001,9.681478],[9.226227,-6.7167,11.45588],[8.329406,-3.852328,9.231437],[9.345326,-0.8438118,9.43648],[-6.216916,2.564122,6.79888],[0.4501454,-7.007302,7.092596],[6.901176,8.498323,10.99308],[6.968523,7.471736,10.26582],[-5.121204,7.937611,9.499073],[-5.862463,4.403567,7.399992],[1.943803,4.026189,4.581328],[-4.539941,9.038058,10.16354],[6.092238,-2.613747,6.704256],[-0.5523768,6.722203,6.818587],[-0.8807149,6.290572,6.43016],[5.17356,5.816909,7.848703],[-3.61358,0.146569,3.752259],[-5.954983,-4.063736,7.278446],[2.460107,-6.288762,6.826468],[2.765574,4.440061,5.325649],[-9.213046,0.4470236,9.277933],[0.7065157,3.178584,3.406253],[4.651875,3.240661,5.756895],[-5.737328,9.083543,10.79017],[-5.973075,-4.505481,7.54831],[3.403608,6.682033,7.565323],[5.225115,2.009542,5.686834],[-7.29562,6.935678,10.11582],[8.341672,-7.670549,11.37633],[-3.155721,-5.422273,6.352922],[2.368886,1.363965,2.910674],[-3.64349,8.323347,9.140739],[-4.535423,-1.0259,4.756315],[1.713539,-8.055218,8.295948],[6.425728,-5.428614,8.471118],[-1.126094,8.788054,8.916163],[-6.74098,-5.139996,8.535829],[5.632917,0.866414,5.786227],[1.026796,-1.406616,2.008203],[0.5926866,7.976244,8.060506],[9.796929,-2.180968,10.08645],[-7.692049,-8.227604,11.30757],[8.352587,3.062636,8.952399],[8.817795,3.112678,9.404375],[3.367173,1.638844,3.876038],[6.772399,-5.985216,9.093306],[9.480343,8.715034,12.91622],[1.259418,-3.196919,3.578607],[-7.784199,-3.230407,8.487007],[3.548122,-0.1848098,3.690979],[5.666448,0.4495046,5.771541],[6.273116,-2.600337,6.863944],[-3.247134,2.639091,4.302172],[0.6157043,-0.9285924,1.497122],[-0.397537,4.462512,4.590429],[-2.363607,7.24121,7.682563],[-5.024988,-8.88337,10.25499],[0.02714735,-6.807304,6.880416],[5.724142,2.73471,6.422183],[-1.556867,0.1652712,1.857727],[2.437657,5.699063,6.278654],[-1.453669,8.222518,8.409694],[7.891091,7.27129,10.77687],[4.973039,3.644826,6.246268],[-2.148363,-1.290172,2.698149],[-9.824419,-0.8524811,9.911909],[-5.071531,0.1047136,5.170241],[8.824073,-2.311228,9.176385],[-9.950633,-2.528361,10.31541],[-1.570319,8.983089,9.173974],[9.715645,-4.918021,10.9353],[-0.4915651,4.475212,4.611849],[-9.966254,-3.040321,10.46756],[9.072865,-3.916114,9.932412],[6.96347,5.874078,9.164863],[-6.153059,-5.597968,8.378387],[6.768522,-1.299626,6.964332],[9.994459,-5.592236,11.49619],[3.145407,0.855683,3.40966],[3.645443,-6.468386,7.491947],[7.3457,-5.992078,9.532277],[-9.367178,-0.6231905,9.440995],[4.726673,3.897429,6.207366],[2.685268,-2.809078,4.012678],[9.635968,1.522252,9.806586],[0.2299369,5.514014,5.608674],[-9.454101,-4.84646,10.6709],[0.4608471,-5.578248,5.685879],[-1.186673,-5.862148,6.064073],[-3.431759,2.963291,4.643066],[-6.378344,0.6249974,6.486439],[-2.741036,5.544598,6.265449],[-6.133494,0.3199404,6.22271],[7.957261,6.265694,10.17727],[5.798728,3.169422,6.683598],[-3.768661,7.400171,8.364529],[3.408745,1.011114,3.693494],[-1.139079,5.473227,5.679235],[-6.1916,6.327544,8.909193],[-0.7524023,0.6898185,1.428971],[-0.2387759,-3.716099,3.855698],[3.492849,4.148553,5.51457],[5.868623,-4.040486,7.194878],[-4.964653,2.2621,5.546609],[-2.064678,-1.370818,2.672459],[4.545645,-1.863695,5.013607],[1.614517,2.21834,2.920222],[-1.98102,-9.795766,10.04398],[3.999891,-2.120296,4.636247],[-0.2504907,-6.307937,6.39162],[0.5715886,-3.276291,3.472866],[-2.232058,2.856742,3.760726],[-2.052534,4.344109,4.907563],[3.162111,-2.354582,4.067309],[-4.224827,-4.067536,5.949287],[-5.236049,-8.583732,10.10429],[-1.377895,0.6234703,1.813094],[-0.03061666,-9.565313,9.617492],[-2.097532,8.389739,8.705593],[-3.39277,8.748192,9.436193],[-1.579357,-0.8387722,2.04888],[-4.559032,6.773366,8.225769],[6.824577,4.727161,8.361872],[7.636345,7.220038,10.55664],[-7.5246,-7.620345,10.75589],[3.663876,3.856132,5.412369],[-9.530212,-7.413086,12.11523],[-2.85618,0.7955877,3.129013],[-7.926993,-1.767873,8.183068],[1.395448,-1.604666,2.349942],[-8.58397,5.97853,10.50844],[8.204805,-8.482701,11.84378],[-2.528067,0.1748075,2.724277],[0.1843074,0.5097941,1.137479],[3.592715,5.228517,6.422226],[-0.8022833,-3.228786,3.474006],[-3.161693,4.916305,5.930122],[6.057115,8.054401,10.12729],[4.605575,8.162604,9.425467],[-8.472075,-2.825166,8.986525],[7.976295,2.750408,8.496236],[7.995648,-4.052494,9.019595],[-8.750119,5.161572,10.20815],[0.8302956,7.378136,7.491748],[-5.487214,-2.759067,6.222698],[-2.988625,1.270672,3.398012],[6.740742,6.875594,9.680465],[3.909247,5.975696,7.210489],[5.770918,-0.7618713,5.906263],[-0.6957083,6.174659,6.293682],[6.895833,5.249075,8.723835],[-0.3838115,-0.6928713,1.275689],[-0.7511274,-6.397105,6.518216],[3.012791,-1.551948,3.533476],[4.677028,0.09618674,4.783706],[-8.901204,-1.24233,9.042943],[-3.030371,-7.915993,8.534992],[-0.9613995,-1.37038,1.949931],[4.934332,1.557121,5.269939],[-8.174003,5.386798,9.840321],[8.881634,4.569489,10.03811],[-7.979372,6.954469,10.63179],[5.003873,-0.1734396,5.105764],[3.199261,-0.9488475,3.483616],[2.167207,1.98123,3.101944],[7.853066,0.9784937,7.976722],[3.572113,-8.587242,9.354182],[-8.462676,9.99012,13.13086],[0.4816113,2.502291,2.737409],[-2.530841,5.836735,6.439925],[7.056058,-9.389963,11.7881],[-1.731091,-9.9827,10.18091],[-1.212852,3.021595,3.406031],[2.751721,-3.254798,4.377862],[-2.199744,0.0502214,2.416898],[5.794235,5.050643,7.751268],[4.580724,4.234765,6.317933],[9.077262,5.069532,10.44494],[-7.234488,7.747513,10.64715],[-1.554882,-7.907702,8.120925],[-5.018642,0.1946119,5.121],[-1.631149,6.671094,6.940039],[6.545647,-3.618124,7.545616],[-2.776318,-8.624308,9.115187],[-0.5863408,9.134809,9.20807],[0.919656,4.392591,4.597893],[9.871752,5.913122,11.55061],[2.794753,-5.047592,5.855667],[-3.159183,-7.935737,8.599789],[-2.079922,1.286963,2.642414],[-4.510835,1.42228,4.834306],[-8.268965,-1.657295,8.492491],[-9.368095,-4.447832,10.41847],[-9.941953,6.770451,12.06986],[6.451815,-0.3535357,6.538418],[3.131323,-5.428794,6.346416],[-9.36201,-4.089815,10.26517],[-4.94766,5.702868,7.615907],[2.407736,9.415155,9.76946],[1.078312,4.063225,4.321175],[-9.115693,8.124866,12.25191],[-5.166039,3.129997,6.122486],[0.446642,-4.485253,4.617032],[-2.786937,8.314138,8.82564],[7.313338,0.1356409,7.382636],[-5.648239,4.250359,7.139199],[-8.738784,9.357939,12.84279],[7.146266,6.306078,9.583096],[0.220617,-9.172583,9.229568],[-0.03352026,-2.108375,2.333745],[2.691955,3.511127,4.535927],[9.384059,1.785999,9.604705],[-6.156694,5.077703,8.042882],[6.981015,-4.238972,8.22821],[4.585124,-0.07844605,4.693562],[-9.805724,-5.060467,11.07974],[-3.262731,-3.058494,4.582553],[5.667136,-5.269591,7.802885],[-8.675536,-6.925065,11.14547],[1.756378,3.281664,3.854112],[2.769966,-4.278598,5.194142],[-3.536299,-0.0502138,3.675313],[-1.234987,1.445216,2.147986],[9.114923,-1.239418,9.252998],[-8.814472,7.052705,11.33294],[6.703195,8.07312,10.54078],[6.624823,0.5549707,6.722817],[-7.23088,1.530312,7.458383],[-3.773198,-3.300552,5.111816],[-8.484586,7.244505,11.20139],[0.9218155,0.1220218,1.365516],[-8.632879,7.231766,11.30597],[-3.8937,-1.567511,4.314857],[0.7757208,1.317376,1.826807],[8.974139,4.182344,9.95124],[8.732759,-9.264638,12.77085],[-8.449502,1.61081,8.659607],[8.78466,0.5525085,8.858641],[-9.39983,2.131582,9.690225],[-7.868934,-8.226571,11.42789],[-1.735859,7.818471,8.07104],[2.627481,-9.197844,9.617899],[-6.687459,-4.552818,8.151703],[-1.675011,-1.921468,2.738193],[3.055664,6.630907,7.369261],[-8.150497,5.363699,9.808153],[-3.491322,-8.374535,9.128098],[-0.4435148,-3.081098,3.269536],[4.929679,-4.024912,6.442178],[6.137059,2.033488,6.542061],[5.574124,0.6252104,5.697521],[5.001685,8.240156,9.69108],[-2.001446,7.894809,8.205717],[7.789741,-2.567505,8.262696],[6.353422,-4.066763,7.609502],[-2.142015,-4.900218,5.440622],[-1.289749,6.903225,7.093516],[5.031083,-5.002819,7.165193],[2.564011,7.359951,7.857674],[0.06994602,5.516962,5.607295],[5.13666,-7.235562,8.929648],[-6.05153,5.860044,8.482991],[-1.300622,-3.144835,3.547055],[4.870419,-1.725187,5.262818],[-3.288009,0.04176957,3.436968],[-6.159395,9.229239,11.14078],[-5.711587,8.442889,10.24229],[4.84835,-6.177685,7.916457],[-3.634527,-6.022688,7.105108],[-2.845341,0.3525958,3.036493],[6.906506,-7.791764,10.45999],[1.120914,3.787187,4.074216],[1.04476,5.511952,5.698521],[7.180264,1.073161,7.328565],[5.877617,-3.672277,7.002285],[2.860896,-2.604263,3.995862],[-5.713141,-6.098545,8.416189],[-7.830185,0.6976722,7.924553],[-7.618977,-4.151505,8.73406],[5.226352,3.264827,6.242904],[1.715426,-0.2553042,2.001966],[8.872216,-1.597548,9.070192],[-3.067755,-4.950715,5.909374],[-1.941772,6.735886,7.081147],[-0.2626123,-9.659324,9.714499],[6.035718,6.456255,8.894556],[-4.963308,4.437322,6.732329],[5.801674,1.333385,6.036335],[-0.5648732,0.9085909,1.464452],[2.564757,-0.1222894,2.755528],[5.314218,-5.021285,7.379309],[-7.895869,-6.990839,10.59323],[-6.660478,5.588873,8.751998],[-9.261439,-8.30958,12.48292],[-0.1152086,6.36035,6.439513],[-7.834631,0.8368785,7.942406],[-1.781197,-5.205369,5.591826],[-9.801784,5.088519,11.08909],[3.198986,-9.589025,10.1579],[-1.204233,3.041464,3.420626],[7.677082,0.6278552,7.767354],[8.025846,9.753731,12.67081],[0.4239451,-0.4594648,1.179338],[-0.3165215,-5.900624,5.993125],[4.667445,-2.119202,5.222649],[5.186003,1.029087,5.380859],[-7.693999,-6.290645,9.988485],[-4.583825,2.670511,5.398433],[0.3795631,7.32406,7.401752],[-4.913333,5.655811,7.558376],[-0.1229653,-9.061055,9.116899],[-7.420374,-6.911229,10.18956],[-2.103309,-6.048307,6.481198],[-3.908605,2.018751,4.51138],[6.865092,-1.542077,7.106863],[6.098489,-9.54079,11.36742],[0.02303421,-2.53566,2.725821],[6.656744,-6.327019,9.23815],[3.894964,-8.596059,9.49015],[-8.904092,-1.74284,9.127997],[-5.448986,-5.452226,7.772915],[-1.082812,-1.581103,2.161566],[9.925782,-7.682561,12.59138],[-0.02044907,1.927344,2.171422],[7.741354,6.244067,9.995846],[-0.5604116,-1.233433,1.68387],[-4.893585,-0.6596279,5.038083],[8.263647,-8.941373,12.21622],[-9.371228,6.87865,11.66772],[-2.887908,5.712253,6.478414],[8.395977,-3.263701,9.063343],[5.558252,-6.340277,8.490776],[-0.361561,8.190516,8.259254],[-9.340508,1.933191,9.590742],[9.625546,2.643094,10.0318],[-2.610256,-8.511839,8.959065],[-2.067486,-1.576757,2.785797],[9.0611,-7.923067,12.07802],[-8.36918,-0.262283,8.432792],[-8.729993,3.661694,9.519494],[-9.963963,-7.091597,12.27075],[0.9667211,-0.2322506,1.410138],[-7.892538,4.037686,8.921606],[9.423474,-3.793584,10.2075],[5.5897,-8.841227,10.50771],[2.02856,9.640434,9.902173],[1.408355,6.508661,6.733953],[-4.168839,-3.392972,5.46731],[0.3661405,4.988121,5.10053],[-0.3243754,1.854155,2.131458],[2.787305,4.820713,5.657592],[8.691643,1.541047,8.883663],[-9.599483,-7.879524,12.45941],[-1.928683,-1.132439,2.449946],[-7.620471,1.319894,7.798314],[1.760347,2.004556,2.849046],[-2.977897,-5.412057,6.257654],[-7.436267,3.551436,8.30125],[-6.587712,-3.878249,7.709654],[-3.114157,5.741285,6.607597],[5.770242,4.143603,7.173921],[4.90867,-4.12878,6.491676],[0.1061024,-5.758353,5.845502],[-1.15713,-0.8682531,1.75864],[8.455318,6.024734,10.43024],[2.87105,3.220447,4.428793],[-1.491341,-5.536266,5.820167],[0.05917628,3.15267,3.307995],[-2.342754,-8.000742,8.396449],[0.6526024,-2.251445,2.548508],[-0.5315793,-2.08757,2.374979],[-9.483057,-1.670519,9.680858],[2.560334,-3.862617,4.740794],[9.344026,-3.675591,10.09063],[-1.620749,-7.15423,7.403367],[-7.54264,9.813334,12.41745],[0.6581783,-0.999714,1.559688],[-3.36078,-4.642632,5.817978],[-1.547312,-6.890494,7.132537],[8.59678,-6.860675,11.04416],[-4.328981,4.72014,6.482268],[-7.059129,-8.230525,10.88912],[-3.707352,7.228424,8.18502],[2.223587,-9.408253,9.71903],[0.9695383,6.595372,6.740841],[9.972059,1.442284,10.12532],[-4.45276,3.444577,5.717708],[0.9263209,-6.223936,6.371456],[-4.770729,2.517374,5.486075],[-4.860175,8.948103,10.23181],[-7.119716,9.231068,11.70055],[-8.271671,-1.72569,8.508734],[4.312801,0.5269353,4.458465],[6.273713,4.736818,7.924451],[-0.7954014,6.875453,6.993177],[2.944426,-2.759865,4.157703],[-1.625123,5.270832,5.605594],[2.010818,-9.987833,10.2372],[-6.566801,8.392813,10.70337],[5.338665,-3.76864,6.610899],[-3.707789,-5.991301,7.116417],[6.727559,9.525448,11.70445],[5.886319,0.4115959,5.984828],[-2.896024,5.736445,6.503365],[6.719172,-8.248734,10.68592],[2.655728,6.225114,6.841413],[-2.520917,8.344558,8.774204],[-8.155748,8.656614,11.93538],[4.750352,-1.773141,5.168159],[7.103421,0.8974991,7.229391],[1.579468,7.742548,7.965034],[9.391182,-8.091545,12.43653],[-4.504276,4.089758,6.1656],[0.8510399,2.62388,2.934113],[7.882614,-1.804811,8.148187],[0.2948385,-8.459833,8.523831],[1.859791,-5.400042,5.798213],[-0.07765443,6.769044,6.842952],[-4.361794,-1.410521,4.691995],[0.8175639,0.626614,1.435638],[-6.932276,7.058878,9.944054],[5.409409,9.917595,11.34109],[1.595567,9.875236,10.05316],[-4.085325,-4.574531,6.214195],[-3.074847,1.709825,3.657621],[0.6487843,5.463596,5.59212],[1.315113,-5.677998,5.913475],[5.064538,1.347727,5.335345],[-4.432533,-3.711971,5.867373],[3.958381,-6.54045,7.71014],[2.567716,-9.019423,9.430968],[-8.280606,2.08127,8.596518],[6.030274,3.853254,7.225771],[5.519918,7.67179,9.503991],[9.794312,-8.575317,13.05621],[1.572587,-4.751002,5.103435],[-9.367412,-2.469072,9.738826],[-5.567148,-4.551611,7.260186],[-5.529717,6.835659,8.848955],[-2.598328,-9.113197,9.528991],[-9.057766,-2.931741,9.572785],[2.310779,9.310882,9.645322],[6.568303,-1.19531,6.750658],[-4.979407,-4.265796,6.632609],[-6.056858,6.028311,8.60384],[-3.534643,-9.023995,9.743008],[7.727399,-9.153692,12.02093],[-3.039054,6.900986,7.60654],[-0.6871678,-4.013977,4.193353],[-2.730926,-0.6901178,2.989017],[-1.276487,-2.072218,2.631256],[-7.862317,1.740177,8.114447],[-0.04788458,-4.815825,4.918787],[-3.218704,-8.877201,9.495512],[-5.738143,-2.598269,6.377875],[4.946555,-8.913199,10.24273],[6.840295,-9.689208,11.90254],[1.514121,-4.667672,5.007966],[9.331351,6.844305,11.61545],[-6.68945,-8.129333,10.5752],[4.622447,-9.699542,10.79111],[8.073965,-5.566976,9.857999],[7.422898,-0.08794449,7.49047],[1.59819,-7.253628,7.49462],[8.495255,-5.438328,10.13631],[0.9893582,-5.891506,6.057117],[-1.302723,-3.134779,3.538916],[-2.2684,1.672571,2.990508],[7.293083,-7.491862,10.50319],[4.559189,3.938592,6.107267],[-7.476426,-4.836211,8.960239],[0.07257926,3.489533,3.630717],[2.578427,8.671405,9.101733],[-1.700708,0.2090265,1.983961],[-3.563574,-1.275161,3.914728],[8.83339,0.5629784,8.907621],[0.8228799,1.545611,2.016444],[6.160541,0.4300254,6.255972],[2.652571,4.527829,5.342038],[-4.908752,7.580643,9.086363],[-0.2475922,0.7442766,1.270925],[8.98802,3.47213,9.687115],[-7.484377,-3.622884,8.375034],[-0.9712468,-3.349067,3.627612],[-5.128825,7.334453,9.005501],[-6.478066,7.813004,10.19845],[2.763827,-0.02792751,2.939306],[8.919458,-6.318466,10.97633],[-1.582768,8.81579,9.012398],[0.5572333,-6.405056,6.506555],[9.146947,4.755019,10.35745],[-2.271315,1.142491,2.732061],[4.451296,-9.975582,10.96933],[-5.91888,4.544013,7.528691],[8.205305,-2.244459,8.565315],[9.367803,-9.146473,13.13064],[0.9444498,6.05405,6.208342],[8.708511,-1.896828,8.968618],[-6.577792,5.737506,8.785575],[-3.668125,3.065928,4.884163],[-6.853799,8.109402,10.66475],[-9.179476,4.13652,10.11798],[2.987286,-7.978736,8.578118],[-1.871645,-4.631244,5.094259],[-8.128564,2.601316,8.593043],[-2.342243,-5.08891,5.690616],[-2.13654,7.546615,7.906718],[-8.120636,5.304657,9.751108],[7.735485,1.75302,7.994423],[2.810771,5.446122,6.209725],[8.534062,-0.2220986,8.595322],[-5.088343,-6.111478,8.015073],[2.031003,7.844086,8.16423],[-7.822056,2.019354,8.140168],[-2.871026,2.710015,4.072711],[-5.293822,3.317883,6.327156],[-4.53763,-2.406813,5.232861],[-0.7144206,8.494853,8.583294],[-0.6928185,-9.255297,9.334908],[-9.031346,6.358394,11.09028],[-4.548903,8.057827,9.307045],[-7.572444,4.761558,9.000797],[-3.907696,2.849155,4.938397],[6.811929,-1.571467,7.062003],[8.026845,7.65428,11.13635],[-5.014406,1.80473,5.422298],[-4.452302,-3.895666,5.999934],[-9.745532,-6.849249,11.95356],[-8.252456,-8.149432,11.64115],[-1.287291,-1.569174,2.262614],[2.035094,-6.387805,6.778323],[-7.998936,6.158072,10.1442],[-5.582253,-6.905292,8.935581],[-4.205225,-0.7847086,4.39314],[-4.977543,1.90756,5.423534],[6.833866,7.364063,10.0961],[8.840868,9.035715,12.68089],[8.119843,2.896391,8.678763],[1.024766,-6.154596,6.318955],[7.523291,2.491617,7.987995],[9.884971,-0.7096169,9.960733],[-4.473674,-9.859607,10.87316],[-1.211977,4.070703,4.363429],[-7.130369,-0.9394756,7.261183],[-8.756721,1.229472,8.898975],[2.015134,0.2345587,2.261809],[5.792131,-0.09246592,5.878549],[-8.942019,5.242732,10.41374],[-5.283311,-6.989999,8.818926],[-3.068269,-5.448164,6.3322],[9.422103,-5.866987,11.1444],[2.217899,6.027002,6.499526],[1.475196,-0.2644929,1.80171],[0.5483621,-0.569298,1.274677],[-6.723847,3.127157,7.482596],[-3.220378,-6.362167,7.200556],[3.885213,1.779958,4.388979],[8.008543,9.819411,12.71053],[-5.400459,7.780539,9.523746],[-3.80077,-5.28319,6.584675],[1.515672,4.31562,4.682076],[0.04704114,-0.4570099,1.100487],[-9.3807,4.534233,10.46694],[-4.317955,2.627459,5.152502],[-5.386005,2.18635,5.898235],[-3.078362,1.657046,3.636222],[-1.55085,3.22317,3.714022],[9.329842,1.534254,9.507886],[9.618481,3.020506,10.13107],[-2.672339,-8.200042,8.682285],[-4.10764,-3.578767,5.538978],[-6.229263,5.48944,8.362874],[-2.213726,-1.306887,2.758357],[-6.686526,-8.625918,10.95975],[-0.8037808,7.062741,7.178327],[7.366685,7.079023,10.26551],[2.67157,-9.112498,9.548555],[0.4371137,-6.308893,6.402593],[-7.375981,-8.662299,11.42106],[4.227077,-3.703878,5.708493],[-4.415784,7.307591,8.596513],[-9.612047,5.260453,11.0029],[0.9410959,2.172843,2.570391],[-9.349108,1.818529,9.576683],[7.367763,7.844082,10.80803],[-5.048435,-9.108312,10.46174],[-8.680584,-4.560285,9.856405],[-1.689643,-7.92039,8.160114],[1.941117,6.196791,6.570248],[3.106349,-2.68605,4.226614],[-4.227938,0.587585,4.384143],[-5.182187,-6.347508,8.255054],[-4.037848,-9.009819,9.923762],[9.963615,-1.0582,10.06943],[-5.58919,-2.017205,6.025626],[-4.398933,9.518809,10.53368],[-8.203061,2.292809,8.575965],[-2.950682,-1.014693,3.276603],[-8.859926,-0.6380008,8.938978],[5.805517,-5.400095,7.991561],[-0.4993446,-8.088068,8.164936],[-7.288155,3.507306,8.149748],[-0.3446994,-9.071862,9.133318],[-4.813023,-4.186533,6.456953],[6.631924,-8.602888,10.90835],[6.793623,-5.574646,8.844771],[-6.295321,5.316723,8.300518],[3.581292,3.971518,5.44046],[-5.594051,-9.496297,11.06675],[-7.024911,-8.619894,11.16476],[-4.089489,-0.3160424,4.221825],[0.9346508,-5.175265,5.353218],[1.88957,-6.630456,6.966593],[4.958542,7.824816,9.317451],[7.695049,8.875834,11.78958],[-5.708388,-9.720697,11.31714],[-7.658229,-5.881218,9.707584],[-3.470033,5.236098,6.360649],[3.423062,-7.583859,8.380469],[-8.505283,-6.043657,10.48168],[9.874221,-6.004384,11.59969],[-4.314255,2.801472,5.240328],[-6.125445,2.006508,6.522818],[-4.701215,8.121949,9.437557],[-9.440939,-5.981909,11.22117],[6.973144,-6.553458,9.621463],[1.958218,-6.354204,6.723877],[-2.080184,2.418175,3.342863],[1.132305,9.261869,9.38426],[0.3114181,-5.076701,5.183616],[5.000511,-4.823949,7.019657],[3.958211,3.3041,5.252096],[4.586076,-4.007062,6.1716],[3.646121,-3.674597,5.272273],[-9.04536,9.094843,12.86603],[3.210931,7.584123,8.296325],[-6.204061,2.103299,6.626782],[-8.139106,-3.34698,8.857049],[0.1923764,-5.829826,5.918098],[-0.5821831,8.26426,8.344875],[6.312841,9.974486,11.84662],[8.134214,2.507849,8.570575],[4.805101,1.150662,5.041132],[-6.5411,7.301963,9.854169],[-9.252667,-1.949382,9.508519],[1.740101,-0.2142676,2.018381],[0.9483807,-8.587557,8.697446],[-1.740935,6.754721,7.04678],[-0.3919378,-5.768196,5.867342],[-1.119636,-1.340234,2.012414],[-4.225517,-5.248886,6.812181],[0.05004711,1.070789,1.465979],[0.2205344,8.511806,8.573183],[2.833916,-3.228281,4.410542],[-1.900849,4.026293,4.563361],[2.40012,4.037052,4.801913],[2.084653,-0.8789852,2.473539],[-0.3727897,-1.107105,1.537743],[7.574677,0.6735095,7.670029],[-2.082255,6.426575,6.829103],[8.700549,-7.504076,11.53303],[0.3297594,3.711161,3.857649],[3.57172,-8.592532,9.358889],[9.486297,-7.634013,12.21753],[-1.455543,6.328349,6.57013],[6.775012,-8.803459,11.15355],[1.237194,-6.189594,6.390753],[-0.8071821,6.896319,7.015038],[-3.590311,-6.186021,7.221994],[-9.25328,2.32729,9.593719],[-2.66138,0.4080593,2.872187],[3.661843,-3.306708,5.034224],[-5.030243,-8.949168,10.3146],[-4.717939,7.70677,9.091384],[-0.800963,-5.171043,5.327404],[-1.514609,7.327822,7.549239],[-5.392818,-3.050046,6.275768],[-9.834284,1.345805,9.976188],[-2.257185,5.88502,6.381876],[8.582619,-7.583958,11.49686],[-3.778784,8.526208,9.379521],[-6.605089,-0.4704065,6.696901],[-1.53807,2.799272,3.346876],[0.3782699,6.548408,6.635114],[7.344037,-9.033677,11.68513],[-5.475315,-1.467997,5.756222],[-8.443909,-5.244874,9.99041],[-9.097788,3.67546,9.862999],[-8.605014,3.831026,9.472224],[6.22117,-3.094599,7.019936],[-7.454628,-2.541923,7.939323],[-0.1863018,0.7802246,1.281975],[-6.003625,-7.085527,9.340675],[8.657167,-7.444147,11.46132],[9.15439,-3.104486,9.71806],[-5.208347,7.554413,9.230169],[4.407073,1.329419,4.710588],[-7.283471,0.2506438,7.35607],[-6.674237,1.936529,7.021081],[0.513122,-9.802138,9.866367],[-2.481669,8.504194,8.915155],[8.048215,2.029599,8.360206],[5.976627,4.436513,7.510175],[0.1913482,0.165432,1.031495],[1.736068,-3.002151,3.609272],[-3.159579,-1.220147,3.531529],[-3.618634,3.80265,5.343656],[-3.839973,7.819675,8.768848],[4.49885,-8.24042,9.441619],[1.557806,-9.492188,9.671007],[-6.694396,3.114524,7.450852],[-5.840158,4.604533,7.503944],[-9.250802,-3.492455,9.93854],[3.199217,5.167014,6.158979],[-2.865358,-6.597273,7.261837],[0.07649726,5.792635,5.878816],[-8.201728,6.090166,10.26443],[5.913901,-2.800275,6.619348],[-5.040035,1.4115,5.328629],[7.689444,-4.93036,9.188906],[2.89882,2.154769,3.747824],[-1.376317,2.320995,2.877719],[1.097776,-9.577212,9.691651],[-4.198784,5.162857,6.729404],[-5.333177,-3.173281,6.285896],[0.4731708,-1.678766,2.010509],[-3.467787,-2.215908,4.235067],[-3.398185,2.353945,4.253083],[9.247824,9.027206,12.96197],[5.12629,7.297447,8.973939],[1.827306,1.292038,2.451205],[6.19278,-3.053747,6.976811],[0.6183995,2.701298,2.946087],[-9.857362,-5.289218,11.23136],[-9.698234,-9.671102,13.73266],[-6.664303,2.340588,7.133813],[-3.959518,8.326156,9.273761],[6.376746,-3.131662,7.174273],[3.058139,-0.3306349,3.23443],[-3.477865,-9.830137,10.47507],[-5.777637,3.487473,6.822284],[-2.272362,3.205062,4.05414],[-4.502324,8.11825,9.336857],[0.8084341,4.404276,4.588161],[-3.20337,1.255677,3.583058],[6.671663,3.509382,7.604397],[-2.843734,-7.924769,8.478725],[-1.746291,6.936928,7.222915],[-5.173424,1.445329,5.463817],[-9.246421,-2.147739,9.545107],[6.034188,-0.393842,6.129155],[8.285094,6.476793,10.56369],[-9.945242,3.268205,10.51613],[6.602489,-9.642103,11.72873],[1.081629,-3.75785,4.036256],[1.593793,-4.964276,5.30888],[-4.380206,0.3940746,4.510155],[-1.460837,-6.576234,6.810353],[-9.279702,-4.899477,10.54124],[-9.552313,4.663843,10.67699],[4.920699,-3.801898,6.29823],[1.330257,-3.425948,3.808767],[9.245344,-0.8129395,9.334734],[7.989747,-5.170611,9.569288],[-5.38701,-4.040065,6.807496],[8.188102,4.037775,9.184152],[9.458018,0.1028924,9.511293],[5.518744,7.790508,9.599403],[-0.5603711,0.2490118,1.17304],[-7.244575,0.4245221,7.325577],[5.485299,-9.575283,11.08037],[8.648648,4.587071,9.840749],[9.384186,1.277118,9.523338],[-3.406436,8.834037,9.520716],[-6.227909,-8.416624,10.51791],[-4.587909,1.575971,4.953038],[7.546268,6.464905,9.987049],[5.326536,1.533739,5.632437],[-5.521974,-0.4231442,5.627721],[-7.884479,-9.739418,12.57065],[3.065085,6.434003,7.196606],[-6.585783,7.794355,10.25302],[-2.756704,-6.531453,7.15956],[-8.055808,-5.821294,9.98917],[-0.4318567,-9.627213,9.688639],[5.503541,-3.025461,6.359432],[7.399652,-5.180482,9.088027],[1.502369,-8.386772,8.578756],[8.741837,9.959373,13.28942],[4.989625,-6.40066,8.17709],[9.907209,-5.484608,11.3681],[-1.398442,5.695308,5.949132],[-5.847764,4.550245,7.476702],[3.574004,2.360997,4.398614],[-7.103703,8.621171,11.21549],[1.33427,-4.727232,5.012683],[6.317934,-8.280746,10.46361],[3.974149,4.737729,6.264178],[5.433543,-8.299885,9.97053],[7.789483,9.804729,12.56219],[6.897159,9.86858,12.08138],[1.368209,8.553389,8.719659],[9.591915,2.331086,9.921633],[-9.112831,4.995437,10.44022],[-8.81927,-4.187082,9.813826],[-3.899167,9.232754,10.0721],[-8.131642,8.784856,12.01238],[-8.234715,-5.553664,9.982671],[-3.63055,-8.058772,8.895206],[-0.5004715,0.3757097,1.179673],[1.288012,-9.501261,9.640173],[-7.020773,-9.883094,12.16416],[-6.672575,-0.2375668,6.751274],[7.073048,-6.087205,9.385205],[-8.147676,1.985496,8.445521],[-2.34111,-0.4088805,2.578367],[-8.655838,6.977238,11.16268],[4.115394,7.500473,8.613568],[-8.026224,-4.965421,9.490821],[3.412475,5.097459,6.215229],[3.463235,0.4871455,3.637486],[4.070653,-9.131362,10.04749],[3.569018,2.724166,4.599888],[3.777491,1.124452,4.066181],[-8.363262,-3.870089,9.269398],[-8.030842,-5.095548,9.563421],[-5.180842,-5.806614,7.845884],[4.726517,-4.58236,6.658678],[2.861328,5.767391,6.515367],[-6.186208,3.37994,7.119913],[-6.204551,3.676137,7.280827],[2.654691,2.364103,3.692745],[0.00846419,-2.976504,3.140008],[7.124101,-3.932832,8.198779],[5.507616,8.216025,9.941675],[-3.408413,5.998132,6.971002],[-0.7528598,2.766899,3.036862],[2.378039,-3.539226,4.379634],[1.215792,0.5692243,1.673967],[-6.472194,5.067191,8.280442],[6.871382,5.323272,8.749464],[3.833009,4.241039,5.803307],[-1.657754,4.503138,4.901673],[-4.142649,9.351177,10.27648],[4.245843,-9.610056,10.55369],[-7.776928,1.081365,7.915174],[0.6313806,-8.498402,8.580296],[7.914124,5.182645,9.51279],[1.999021,-1.984016,2.988713],[-8.608974,-8.834377,12.37581],[-0.1842621,-9.104429,9.161036],[-4.463251,-2.947528,5.441372],[2.617846,6.323191,6.916348],[-4.388184,8.522636,9.638023],[2.353544,-0.8377615,2.690914],[-1.778016,-7.021984,7.312291],[9.007226,-7.36319,11.67676],[4.927439,4.919424,7.03423],[-5.102855,4.182783,6.67344],[8.886512,-2.715543,9.345816],[-6.661032,9.249573,11.4422],[-3.094651,-9.657151,10.19007],[-3.241611,-1.4479,3.688422],[0.8944147,2.710047,3.02396],[2.345039,4.716369,5.361282],[-0.6578109,-2.479692,2.753468],[6.721481,1.631987,6.988683],[4.470157,-3.212457,5.594835],[8.277605,-2.32863,8.656862],[-4.952847,3.876831,6.368713],[-0.1358875,7.167522,7.238221],[0.3642903,8.061838,8.131785],[1.641988,2.535964,3.182332],[2.860137,-2.187816,3.737234],[5.73823,-7.24674,9.297447],[-6.782026,-6.292163,9.305224],[5.100291,-4.968967,7.190522],[-1.883052,6.345166,6.693804],[0.1683278,0.08346415,1.017497],[-0.2377991,-8.33717,8.400294],[5.627065,6.29052,8.499088],[9.353195,-9.818936,13.59757],[-1.621771,-3.198608,3.723068],[7.192931,0.9689987,7.326474],[8.63508,5.086285,10.07149],[2.102231,8.65371,8.961366],[9.894201,-2.238774,10.19349],[8.190392,-7.405055,11.08681],[-6.345327,4.307199,7.734025],[4.208673,-8.746495,9.757771],[6.62853,-5.426552,8.624667],[9.473272,-9.005304,13.10871],[-6.394259,3.020326,7.142053],[-6.431546,-4.909522,8.152802],[-9.964492,-8.711798,13.27353],[4.491734,-7.062324,8.42924],[-7.067148,-3.036106,7.75645],[-5.554258,3.732126,6.765985],[-8.709225,3.220101,9.339146],[9.855475,7.142904,12.21276],[-8.672705,6.125474,10.66477],[-0.3241166,-8.519677,8.584285],[-8.106476,2.352405,8.499928],[5.463415,-4.317992,7.035194],[-1.520319,-8.007737,8.211895],[3.646886,-7.204578,8.13669],[4.52784,-3.519991,5.821656],[3.401963,9.963077,10.57527],[9.547862,4.640257,10.66272],[9.855782,8.786025,13.24125],[3.43933,9.114237,9.792768],[-4.008897,-7.795266,8.822552],[-9.374822,-7.859648,12.27442],[-8.837395,1.359964,8.997169],[4.081876,-5.138642,6.638325],[9.103245,3.191905,9.698316],[4.696304,-1.960799,5.186522],[-4.296092,-1.033506,4.530401],[-7.506926,1.670453,7.755279],[-3.653022,-9.183108,9.933481],[4.386205,8.122967,9.285547],[7.104578,-5.202581,8.862385],[-6.335342,1.511901,6.589568],[8.844131,3.67276,9.628489],[-0.0836802,6.65635,6.731567],[-4.216298,0.9566904,4.437614],[0.9624963,-9.949294,10.04564],[1.696162,-5.460142,5.804318],[-6.894196,2.252523,7.321461],[-1.098641,-4.780403,5.005923],[-0.92011,7.827037,7.944125],[-3.273946,8.528016,9.189439],[-0.5921013,-4.662289,4.804947],[-3.734113,-3.244865,5.047054],[-2.821681,8.847816,9.340543],[5.973403,1.718313,6.295565],[-5.413475,5.122678,7.51981],[0.2748265,-9.278308,9.336087],[-2.416136,-8.074759,8.487606],[-2.438126,-9.655097,10.00826],[-0.860733,-3.591719,3.826396],[6.814395,3.712067,7.824028],[-0.1512861,9.715086,9.767588],[2.079297,9.722048,9.992081],[0.1602747,6.15919,6.241899],[7.150445,-4.351374,8.429907],[-4.912046,7.253247,8.816903],[3.313258,-0.4364766,3.488294],[-8.007559,-5.401042,9.710421],[-8.773735,5.343706,10.32151],[1.812785,9.764929,9.981985],[-6.330198,6.353443,9.024282],[4.117129,-8.961232,9.912337],[4.135772,2.781761,5.083582],[-3.709543,8.093154,8.958786],[-6.903333,-9.40776,11.71162],[4.920719,3.360218,6.041899],[5.517754,-6.235108,8.385832],[7.483998,7.84889,10.89106],[-1.422322,5.242317,5.523123],[-1.297191,3.041083,3.454112],[-0.6932096,5.432197,5.566803],[-2.660855,-2.715499,3.931168],[-6.526903,-4.02987,7.735652],[-8.994722,-2.477776,9.383198],[1.344755,1.114255,2.012444],[-8.310623,8.257821,11.75832],[4.012216,5.260782,6.691316],[1.860358,4.230139,4.728109],[-5.680151,6.141972,8.425433],[-4.631479,-6.192647,7.797402],[8.768841,4.147613,9.75168],[7.27931,-7.327649,10.37703],[7.475484,-4.87677,8.981411],[2.175358,6.071966,6.52694],[9.343933,-8.546568,12.70248],[5.521414,-4.566196,7.234373],[2.110999,-0.1277481,2.339367],[6.133238,-0.6194916,6.245028],[7.636235,1.936895,7.941262],[-8.193365,-6.21764,10.33394],[-4.197696,9.970716,10.86443],[-3.293712,-6.870378,7.684441],[6.501111,-2.273584,6.959427],[-5.169088,8.631076,10.11014],[0.6799436,-2.828709,3.076349],[-4.625515,9.848477,10.92648],[-4.039587,9.976832,10.80997],[-2.381156,0.7002983,2.675878],[-4.320964,2.103901,4.908884],[-8.828771,-7.127019,11.39042],[-8.424398,6.392712,10.62249],[-9.468686,1.256209,9.603857],[8.409105,-5.301296,9.990835],[9.542516,7.661332,12.27826],[-9.376445,9.593121,13.4516],[-9.22235,-4.24204,10.20033],[7.007424,9.567676,11.90145],[-0.5157887,6.579267,6.674788],[-8.498514,7.320752,11.26136],[-1.186122,-4.447841,4.710644],[-6.154821,6.352644,8.901568],[-2.483309,0.5262555,2.728327],[-6.133219,3.440567,7.103089],[-1.850801,9.342095,9.576023],[-1.910659,8.499621,8.768932],[8.674284,8.792012,12.39123],[1.21802,1.568901,2.223741],[2.443402,-6.278926,6.811397],[-1.049228,-4.955646,5.163265],[0.4471389,-7.240188,7.322586],[-8.004623,-8.935256,12.03797],[5.096827,8.700178,10.13266],[-3.356565,7.372504,8.162128],[6.5272,-2.342362,7.006497],[7.527561,-2.266067,7.924597],[8.803021,5.676377,10.52209],[-3.289887,-8.364706,9.043875],[-7.199174,-5.939231,9.386297],[5.4246,5.061669,7.48644],[-3.763376,5.353172,6.619627],[-7.779036,-1.10065,7.919901],[-9.037905,2.002374,9.31092],[-6.209653,4.973032,8.018157],[-6.868597,3.27031,7.672845],[-1.084399,-8.011011,8.145687],[-0.19871,8.661882,8.72168],[-7.474057,-6.421443,9.904365],[-5.824407,-1.944304,6.221257],[-4.237559,9.614012,10.55396],[2.97298,6.17842,6.929032],[-5.288781,-4.640368,7.106632],[0.09748384,9.271699,9.325979],[-2.733221,0.342806,2.930531],[3.660259,-2.800046,4.715693],[2.85258,-7.243559,7.848972],[-4.827399,4.741879,6.840263],[-8.312787,-1.768331,8.55742],[3.003031,-4.953661,5.878517],[-5.238241,3.793713,6.544572],[8.819159,0.5777015,8.894453],[-2.218428,2.262863,3.322946],[-9.891909,2.924544,10.36353],[-7.605673,5.091452,9.207016],[-8.677004,3.638836,9.46211],[-7.400281,-3.116132,8.091627],[-9.340113,3.056729,9.878325],[0.7392623,9.973129,10.05036],[5.011908,-6.588912,8.338643],[-6.113444,-5.701231,8.418921],[-1.37468,-9.605685,9.754945],[-5.036319,7.219014,8.858819],[7.015351,1.579307,7.260121],[5.337308,3.877641,6.672553],[6.078102,0.1742336,6.162279],[1.857093,3.106973,3.755273],[9.646472,-8.18922,12.69322],[-7.806624,7.608266,10.94665],[3.953963,-1.293243,4.278586],[4.878132,-8.193668,9.588138],[-4.241272,-8.245793,9.326387],[4.583147,-6.013002,7.626364],[-5.163519,-1.259569,5.408183],[5.424809,8.846239,10.42519],[6.77919,-9.802999,11.96061],[-5.816376,-8.00568,9.94591],[3.37724,0.2355142,3.530045],[-1.272229,9.507941,9.644662],[-7.036113,-4.112913,8.211147],[2.528866,-7.632853,8.102814],[-5.978845,-7.470517,9.620562],[-2.004041,4.917244,5.403283],[-6.527682,2.639809,7.111907],[-6.273365,4.752095,7.933318],[-1.736331,-3.436656,3.978121],[9.446387,7.11002,11.86535],[2.355834,-7.306991,7.742226],[-7.407826,-9.341113,11.96379],[-4.105353,-3.655481,5.58717],[7.408773,3.199331,8.131766],[-4.096444,5.550698,6.970732],[1.229261,-6.558821,6.747534],[-1.404189,5.2798,5.554101],[-7.562,4.633255,8.924735],[-7.32924,4.482044,8.649074],[4.327305,-3.061116,5.394071],[5.433433,9.411162,10.91294],[3.35107,7.270765,8.068067],[1.658204,7.71437,7.953687],[2.350625,-7.385228,7.814539],[6.32415,-5.017839,8.134715],[5.729139,-9.89349,11.47624],[-9.578288,-0.5986705,9.648938],[-4.395489,7.866179,9.066261],[8.887831,8.516375,12.34999],[-6.822834,-0.6655062,6.927767],[-8.697207,1.200575,8.836447],[8.101689,0.754006,8.19792],[-5.055212,3.602833,6.287732],[7.599844,9.451752,12.16936],[-1.059859,-3.5209,3.810517],[1.470723,4.609729,4.940913],[-7.235186,9.021031,11.60719],[-8.737629,2.930336,9.270007],[8.50043,6.655081,10.84193],[-9.985733,5.154525,11.28202],[-2.324172,7.092688,7.530471],[2.539665,7.234468,7.732234],[5.943471,2.960564,6.714893],[-4.998644,-2.594707,5.720047],[-8.557516,4.090441,9.537442],[-5.356759,-4.27804,6.92795],[3.506959,-5.752959,6.811409],[2.303844,3.468144,4.282023],[6.788889,4.353317,8.126523],[-5.470558,4.186804,6.961058],[-9.826497,-3.168422,10.37299],[-0.6345074,-0.6013285,1.32823],[-1.662757,1.622342,2.529181],[3.746644,-0.4937233,3.909105],[-6.482711,-0.7381079,6.600784],[-0.9160557,3.810504,4.044639],[5.727883,8.777951,10.52906],[-2.557914,-4.135497,4.964399],[-7.242374,-6.882777,10.04115],[6.706203,-5.216043,8.554546],[-5.296378,-8.201425,9.814019],[1.614036,6.729176,6.991918],[4.228331,-9.966308,10.87226],[2.584314,0.7999793,2.884206],[1.098211,-5.781582,5.969318],[-7.94208,6.718913,10.45086],[9.699252,-7.688447,12.41723],[-2.472052,8.95604,9.344608],[-2.566027,-5.532864,6.180378],[-7.843612,9.764907,12.56486],[-5.953922,-8.764696,10.64279],[-7.635367,2.173394,8.001404],[5.715504,-5.938817,8.302803],[3.022935,-7.281115,7.946872],[-2.195314,-5.32821,5.848866],[-6.696933,-4.284548,8.012881],[-6.22815,-5.397714,8.302118],[2.953192,7.704,8.311014],[3.170845,-6.042594,6.896898],[4.82937,1.535101,5.165206],[-6.632435,5.815989,8.877777],[3.752564,-8.694406,9.522312],[6.510014,-6.787874,9.458093],[8.186159,-6.36848,10.41973],[-9.778513,1.346002,9.921242],[-3.487041,-7.269691,8.124522],[6.533079,9.02444,11.18578],[-3.277369,6.754363,7.573808],[-1.585491,-1.478712,2.387545],[3.667328,7.29227,8.223534],[-3.456783,7.962462,8.737858],[-5.026861,-3.042894,5.960582],[-8.215939,-7.917788,11.45395],[-9.427879,-8.75769,12.90667],[5.036247,7.801427,9.339489],[2.685361,-6.427624,7.037437],[-3.024619,5.015242,5.941462],[8.214381,-1.504491,8.410681],[2.52175,-3.294428,4.267608],[-1.028031,-8.491248,8.611512],[0.4200528,7.577144,7.654382],[-4.506796,-3.651003,5.885664],[3.507754,-5.383912,6.503141],[1.14214,-7.265862,7.422751],[-6.23141,-7.945045,10.14664],[-8.28712,-2.834662,8.815421],[-2.693166,-7.707397,8.225395],[-4.678989,-6.301057,7.91178],[-4.621791,-6.385323,7.945647],[9.883827,8.125211,12.8339],[0.175563,-3.808017,3.941043],[9.229416,-2.438194,9.598276],[2.489372,-0.8037568,2.800535],[-0.8205624,-5.489109,5.639471],[9.596234,0.1865652,9.650001],[2.065513,9.06792,9.353797],[-5.297056,0.699661,5.435837],[-6.003742,6.078227,8.60173],[-7.575471,-9.330177,12.05985],[-8.079981,0.5243261,8.158493],[-4.66063,1.377823,4.961841],[-1.742158,-5.782901,6.121851],[-1.7212,-6.557142,6.852638],[-2.866317,-4.316349,5.276991],[-2.947071,5.304501,6.150037],[6.001,-9.7129,11.46091],[0.3601148,6.074141,6.16643],[9.771941,-4.692154,10.8861],[-5.913303,7.668077,9.734811],[-4.647981,-3.48841,5.896841],[-0.8620473,2.128571,2.504783],[-7.634768,1.733547,7.89271],[1.300079,-2.671776,3.135058],[6.810694,0.3985775,6.895247],[2.860226,-0.9004071,3.160953],[-9.223289,-0.3708122,9.284749],[2.855296,-4.329863,5.282085],[-3.104126,5.606508,6.486026],[-6.328829,-5.23476,8.273862],[5.052829,4.11543,6.593014],[7.551884,0.6950386,7.649446],[-0.3361078,-5.098318,5.206325],[6.510005,-7.51157,9.990189],[-0.5843338,3.991308,4.155958],[8.556939,4.188618,9.579442],[5.527631,4.296954,7.072377],[-6.362906,5.712953,8.609553],[-0.04065396,3.533752,3.672745],[-6.249732,-2.306056,6.736248],[-1.376093,6.638126,6.852617],[4.829661,-5.226271,7.186066],[2.653857,4.368217,5.208097],[-8.712836,5.090327,10.14026],[5.254259,8.438421,9.990705],[-7.454636,-8.122212,11.06987],[9.444942,-6.833362,11.7005],[0.4633803,-0.1555532,1.113067],[7.822089,4.155527,8.913669],[-4.341193,7.507366,8.729633],[7.153741,3.069317,7.848358],[-5.025961,-0.900012,5.202914],[-7.083825,1.885503,7.398357],[2.71598,-7.053759,7.624438],[1.125293,-8.533052,8.664829],[-4.097926,3.829888,5.697459],[5.153508,-5.663762,7.722489],[9.027788,-6.824181,11.36092],[6.511296,-6.288052,9.106952],[5.097913,-2.640175,5.827456],[-8.679426,9.772799,13.10878],[-8.267003,-9.416638,12.57046],[1.516385,-9.359751,9.534378],[2.678054,3.352383,4.405728],[3.856893,-5.029574,6.41656],[-0.9105905,-1.371008,1.925834],[-4.811333,-5.347414,7.262491],[-5.365844,9.008847,10.53336],[6.744812,-3.000938,7.449706],[0.9295557,0.7049063,1.536544],[8.848479,-9.129352,12.75306],[5.762689,4.940704,7.656314],[-6.756377,-8.254164,10.71354],[3.949971,-6.261746,7.470725],[-5.73892,8.180806,10.04295],[-4.480615,-6.12644,7.655663],[-1.640548,1.789111,2.625322],[7.8557,1.276933,8.021382],[-2.967768,8.063034,8.649865],[8.881578,3.32091,9.534719],[9.012138,0.5138011,9.081995],[2.733267,-9.012237,9.470541],[-7.119005,4.453341,8.456506],[8.670476,-9.234056,12.7061],[-8.626374,-8.543777,12.18238],[-0.1266658,4.551849,4.662122],[-2.712544,0.9101889,3.030897],[1.982368,8.471878,8.757997],[-9.130811,-2.352025,9.481758],[9.435747,6.764711,11.6531],[7.465111,-4.810232,8.93679],[2.978899,7.785901,8.396075],[2.308105,-3.50225,4.311973],[-7.500843,6.320781,9.859763],[4.012511,4.153818,5.861267],[-5.885612,-4.282483,7.347115],[-4.847686,-4.755383,6.863944],[-1.242047,-7.948472,8.106842],[4.30695,-1.006376,4.534602],[6.221229,-3.092966,7.019268],[-0.9588438,-4.719742,4.918877],[4.396189,-9.739716,10.73259],[2.474442,-3.967891,4.781947],[7.926678,0.4565772,8.002542],[-7.582646,6.028792,9.738729],[8.742036,-5.678376,10.47221],[-0.8497955,0.1298555,1.318717],[-4.61554,-9.15615,10.30234],[-0.04802005,-8.787118,8.843966],[-4.769739,7.729351,9.137465],[3.526048,-5.351829,6.486531],[3.774125,-3.670161,5.358554],[-6.791216,-9.304497,11.56262],[-0.3528547,-8.474853,8.540939],[-3.200744,-4.202963,5.37677],[-5.111928,1.781026,5.504894],[8.154201,7.637369,11.21697],[-4.762461,3.893652,6.2323],[7.594437,-6.039321,9.754428],[1.734547,1.476875,2.487933],[-6.351333,-8.567682,10.71189],[1.016546,-4.496115,4.716823],[-0.8780234,-4.051511,4.264465],[7.209317,0.09051945,7.278904],[-1.768653,8.024791,8.278007],[7.548349,-0.7616227,7.652297],[9.37438,7.219262,11.8742],[4.994645,-0.8333626,5.161489],[3.048025,6.721196,7.447478],[-6.640135,8.003376,10.44727],[0.4618702,-6.983134,7.069476],[8.019113,9.899201,12.7789],[6.478812,4.825059,8.139791],[-6.137762,-7.88561,10.04266],[-6.404296,9.114817,11.18458],[0.3082063,-5.21245,5.316449],[-6.811054,-6.385899,9.389896],[-2.667733,2.832849,4.01769],[-6.17162,1.795247,6.504752],[3.297254,4.605373,5.751638],[-4.450625,7.798906,9.034987],[6.384295,-3.07068,7.154599],[-1.736159,-5.807827,6.143703],[-9.513812,1.13933,9.633831],[8.790201,1.701384,9.009015],[-1.848478,-3.955659,4.479298],[2.771293,-0.07399832,2.947124],[-4.046422,-5.522366,6.918819],[-6.003879,6.486986,8.895367],[-0.3133254,9.241992,9.301214],[-5.444848,-9.289721,10.81412],[5.033078,3.130363,6.010911],[-2.52609,8.729855,9.142838],[-8.615582,1.57312,8.814929],[6.994019,-7.619315,10.39087],[4.720222,8.544467,9.812665],[-4.65344,8.673052,9.893247],[3.47357,-9.514304,10.1778],[7.65667,7.613498,10.84389],[8.7818,0.7629068,8.871417],[5.211605,8.286202,9.839816],[0.5741524,7.707094,7.792878],[3.041765,4.235093,5.30927],[9.733126,-8.778148,13.14495],[-3.172102,-6.62955,7.417086],[-3.585435,-1.200684,3.911136],[9.319433,2.984244,9.836541],[-3.978042,-0.3739226,4.118815],[6.70028,4.199738,7.970668],[8.481762,3.883112,9.381836],[0.9457899,-3.957003,4.189557],[-5.914163,-9.869108,11.54888],[9.938349,-5.185156,11.25418],[-2.505458,0.0580698,2.698276],[-1.921153,-4.55277,5.041681],[-9.590449,2.764656,10.03095],[-0.2582407,-4.309671,4.431699],[-7.280403,-0.9532853,7.410332],[-0.518421,-4.985213,5.110881],[4.585272,-6.970812,8.403388],[7.681074,-1.959543,7.989913],[8.328671,0.9565075,8.442846],[-1.704406,-7.396946,7.656358],[-8.901085,-3.698011,9.690439],[-5.670837,-1.433594,5.934103],[-8.260428,2.595154,8.716048],[-9.404058,-8.026395,12.40401],[-0.07616752,5.155213,5.251859],[6.418629,-8.617049,10.79131],[-1.321802,9.798922,9.93811],[3.892201,4.120944,5.755989],[8.069643,7.748674,11.23215],[-4.156975,4.420312,6.149765],[7.140532,4.984724,8.765539],[3.324951,7.085734,7.890687],[-1.543996,-8.189936,8.393984],[8.899469,-2.217864,9.226022],[-8.608215,6.265348,10.69374],[-9.483435,4.900776,10.72162],[-1.598838,4.653937,5.021495],[2.767783,7.695614,8.23912],[-4.943755,0.7661147,5.10173],[7.747374,-7.632473,10.92138],[-8.817636,2.353732,9.181001],[7.324375,-9.891635,12.34872],[0.9977325,5.799835,5.969385],[-4.369578,-1.798315,4.829819],[9.340625,-8.486835,12.65992],[7.087702,3.186739,7.83523],[-4.889749,5.190933,7.201072],[-9.815722,-8.026232,12.71884],[-0.225533,3.598889,3.742041],[5.235047,-0.5369509,5.356681],[-9.27564,-7.786441,12.1518],[-9.706514,-8.311066,12.81757],[2.013272,5.192687,5.65838],[1.072282,-7.981505,8.115061],[4.311873,-0.05827542,4.426697],[8.598612,-0.3647695,8.664248],[8.673545,-2.605576,9.111499],[-0.7689646,9.789647,9.870587],[-0.178725,-2.133846,2.363312],[-4.351927,5.546595,7.120673],[-4.274554,0.6669027,4.440335],[8.089087,-8.6682,11.89836],[4.806416,-3.858615,6.244241],[-9.11662,-9.429894,13.1543],[-7.916379,-9.074356,12.08358],[-0.1070311,1.895266,2.145574],[6.279397,9.520947,11.44899],[7.370306,-8.617706,11.3836],[2.786037,6.506184,7.147897],[6.540406,-6.873806,9.540761],[-0.1309329,-5.147186,5.245061],[6.774303,9.086212,11.37763],[0.5188863,-7.642477,7.725069],[-2.906109,-4.373397,5.345285],[0.8935785,-1.982793,2.393731],[6.086888,3.611365,7.147878],[4.562513,-2.79034,5.44082],[1.430404,2.75789,3.263742],[4.350205,0.7609032,4.528052],[-2.373484,9.629734,9.968209],[1.343691,2.337634,2.875767],[-6.779666,0.09667151,6.853701],[-2.513795,2.378929,3.602564],[-5.542003,-0.02929618,5.631577],[7.728151,3.468508,8.529646],[0.008541658,-9.269371,9.32316],[1.68286,0.962463,2.181365],[0.4229776,-9.827638,9.887435],[-5.29739,-3.713336,6.546083],[6.972241,3.466967,7.850605],[-3.050101,0.7290826,3.291607],[3.931128,6.104918,7.329651],[9.961454,4.77858,11.09348],[6.509902,-2.89937,7.196192],[-5.904479,4.440134,7.455043],[9.843227,8.289136,12.90732],[4.748306,1.980513,5.241073],[8.919592,-9.707728,13.22116],[3.76708,-2.384317,4.569011],[-4.564688,0.02764288,4.673023],[1.20863,-0.5613021,1.666087],[-7.018301,5.095172,8.730254],[-4.429016,-5.080137,6.813514],[9.094856,2.965081,9.618114],[-1.921833,-8.718095,8.98324],[-4.996882,-6.135042,7.975436],[-8.384532,5.8051,10.24693],[8.089399,8.134583,11.51563],[2.501767,2.164517,3.456006],[7.609864,4.317207,8.806151],[5.714136,2.287381,6.235661],[4.434064,-5.597634,7.210717],[-8.737966,7.970977,11.86965],[1.774915,9.445719,9.662915],[-2.916885,6.122565,6.855218],[4.253862,-2.637668,5.104178],[5.357363,7.898254,9.596028],[3.397299,0.6603549,3.602459],[5.529535,7.687595,9.522336],[8.930492,-7.229106,11.53316],[-4.3671,-4.140095,6.10016],[5.508986,-3.62624,6.670722],[2.380994,-6.316705,6.824214],[7.988701,0.8453712,8.095306],[-4.116508,3.902705,5.759926],[0.7929236,8.930547,9.021275],[2.837278,7.881638,8.436253],[2.297065,-2.520831,3.554026],[-0.1972847,-8.534089,8.594743],[5.619083,8.636233,10.35174],[-8.418441,7.475679,11.30292],[7.322587,9.53005,12.05994],[8.554487,-7.799382,11.61936],[-8.7601,-3.913902,9.646657],[4.757094,0.6563044,4.905169],[2.412698,-8.433375,8.828529],[0.06849057,-5.372314,5.465021],[1.325732,-7.086471,7.278437],[-2.758644,-6.454037,7.089761],[-0.738654,9.836713,9.914965],[-9.615586,-4.923341,10.84891],[-2.425974,2.019445,3.311119],[9.72805,0.7636763,9.809086],[-7.754606,3.921351,8.747051],[3.712618,-6.628274,7.662737],[-5.12988,5.487108,7.577864],[8.336836,8.632369,12.04245],[-7.926429,6.572209,10.34515],[-8.406075,9.532824,12.74899],[-3.080389,4.831837,5.816824],[0.386309,6.745194,6.829852],[-8.867058,6.639728,11.12253],[-8.535024,-9.468293,12.78652],[3.586074,0.2255294,3.729717],[-1.308353,9.477373,9.619375],[-2.239104,5.420539,5.94944],[2.722982,-7.1358,7.702874],[6.04682,-7.668456,9.816784],[-0.09225792,9.006059,9.061877],[-9.720293,-3.254287,10.29925],[-8.613657,-9.097463,12.56817],[-0.5259464,-8.394044,8.469746],[-6.943306,2.036621,7.30461],[-7.239037,7.385971,10.3902],[-2.24103,-9.768764,10.07229],[-9.445622,4.854628,10.66711],[0.6811336,0.3635109,1.263362],[6.12316,-2.148799,6.565853],[4.218693,-2.479451,4.994502],[-6.700991,-0.748849,6.816455],[-6.710789,4.876767,8.35569],[-3.499238,3.743524,5.220981],[-4.232322,7.146689,8.365866],[-0.8796306,-4.681028,4.866804],[5.61755,-7.425023,9.364179],[9.018975,4.259015,10.02403],[-1.573974,-6.32147,6.59078],[-9.562936,-2.469244,9.92708],[6.490421,-8.973615,11.11986],[-4.667984,-1.712872,5.071884],[4.201763,-8.720414,9.731415],[-0.5117332,-0.1817929,1.137945],[-2.281942,-1.269995,2.796453],[-4.285253,5.141473,6.767432],[-8.945993,7.867944,11.95556],[9.914631,7.015563,12.1868],[0.8510534,-1.829189,2.251716],[-9.816579,3.302569,10.40539],[9.536884,-6.549827,11.6126],[-9.640514,8.866539,13.13602],[3.529039,-5.030427,6.225698],[-1.025225,0.5439321,1.531975],[1.80765,3.104873,3.729321],[4.133254,9.614061,10.51256],[-4.882473,-9.331009,10.57858],[-6.652891,-7.755595,10.26695],[6.480993,-6.496277,9.230649],[-0.496778,4.433576,4.572022],[3.992181,-9.018338,9.913018],[-1.836475,-5.866832,6.228351],[6.925438,-1.582324,7.173942],[-0.9155752,-6.797914,6.931804],[-3.774951,-2.70752,4.751939],[2.781934,6.791641,7.407128],[-3.183555,0.6631463,3.402173],[8.736612,-9.68168,13.07912],[-2.718114,-0.3984241,2.923506],[-4.051464,0.4496708,4.197209],[6.415598,1.271961,6.616478],[-2.380747,7.2704,7.715353],[5.449648,9.015038,10.58157],[7.686534,1.748236,7.946014],[9.042533,8.943469,12.75747],[-6.532913,-3.223767,7.353341],[-7.083174,-4.068359,8.229392],[-5.9045,2.622423,6.537601],[0.4420698,-3.385568,3.557738],[-0.08360122,2.324774,2.532106],[6.982682,-6.520504,9.605978],[-0.5552057,9.227955,9.29857],[-5.692018,-0.3693993,5.790986],[-1.630938,4.236986,4.648871],[8.466703,3.677413,9.284849],[-7.725689,-6.379553,10.06901],[6.721283,-2.046502,7.096747],[3.708242,-9.61476,10.35349],[9.647709,-9.304895,13.44096],[-7.42188,-9.826145,12.35465],[5.797779,-9.084334,10.82309],[0.3312819,1.217724,1.610155],[7.198078,-1.491186,7.418622],[-5.623948,0.182703,5.715083],[-9.472306,-5.447819,10.97285],[6.241331,7.408358,9.73848],[4.733553,-4.11825,6.353464],[0.5441359,6.243438,6.346385],[-1.268018,3.850943,4.175839],[-8.260001,9.783139,12.8428],[2.35613,2.974897,3.924456],[-8.311852,-7.454811,11.20987],[-6.795971,7.802593,10.39546],[-8.157111,0.1687097,8.21991],[-4.371194,-0.6546242,4.531652],[-7.272781,-1.971629,7.601359],[5.472095,-1.802168,5.847361],[-0.9563451,2.883364,3.198184],[-0.2231359,-1.119438,1.517541],[-1.809216,-8.575486,8.821124],[-5.635319,-2.525741,6.255892],[8.175455,1.087373,8.307855],[-9.798319,-8.028528,12.70686],[-8.469098,1.430329,8.647049],[8.250371,-9.245125,12.43145],[5.656221,-5.795274,8.159537],[1.629996,-5.637427,5.952939],[7.583842,3.821485,8.55093],[-2.549341,1.930063,3.350266],[-4.422685,8.96384,10.04543],[-9.218658,-6.26207,11.18915],[-5.019064,-1.304018,5.281238],[1.921373,7.531694,7.83697],[9.759215,6.267434,11.64143],[0.5014624,-7.000504,7.089324],[4.988159,3.023966,5.918285],[-8.875,9.796667,13.25671],[-0.3532771,2.080671,2.335379],[-0.6190489,-0.7729082,1.407341],[-2.888308,-5.814084,6.568553],[-1.645858,5.028585,5.384748],[7.00354,0.7855525,7.118052],[-0.3237876,-4.532771,4.653047],[-5.292446,-4.137267,6.791683],[1.91639,2.640853,3.412719],[8.231642,8.906877,12.16932],[-7.397657,-8.296646,11.16063],[9.817116,6.358974,11.73935],[2.039426,-9.406781,9.677127],[-9.398036,1.444491,9.560839],[-6.738308,1.102808,6.900796],[2.209942,-5.23971,5.773942],[-0.2415967,-0.2630748,1.061874],[6.92971,-4.775293,8.474922],[8.826113,-6.960929,11.28516],[-4.835144,6.964148,8.53686],[-1.221025,5.507308,5.728992],[2.929848,-6.134605,6.871491],[8.161365,6.524289,10.49639],[4.945938,-0.765834,5.103803],[-5.212814,4.67788,7.075027],[4.317467,-9.252325,10.25895],[-3.126554,1.746284,3.718177],[2.14129,-4.961522,5.495618],[-6.597948,6.896719,9.596752],[0.822346,9.968937,10.05266],[-6.406776,-9.870483,11.80988],[2.099644,-2.945622,3.753025],[6.547447,-3.071651,7.300966],[-3.87085,-9.452687,10.26337],[8.872138,3.275737,9.510272],[1.920628,6.6042,6.950127],[3.088315,2.746783,4.252353],[-8.494967,-3.473881,9.232135],[-6.246355,8.406932,10.5211],[3.931929,-4.315178,5.922907],[-3.165942,7.030084,7.774655],[7.217946,7.315255,10.32529],[3.589253,-0.2991308,3.737943],[-8.570781,-2.231697,8.912842],[-1.274755,3.903455,4.226342],[-5.964087,4.391073,7.473411],[-7.845904,-0.7299004,7.942983],[1.024719,-7.396264,7.533576],[-7.454608,-0.2362343,7.525091],[-1.767584,-8.497818,8.73712],[3.434939,-5.737662,6.761625],[-3.414565,-3.706123,5.137568],[1.418871,-1.683507,2.418138],[1.751529,5.164093,5.543979],[-6.761735,-6.287154,9.287054],[-2.597889,4.734739,5.49243],[-5.509061,-1.511925,5.799626],[9.839247,-4.603515,10.90886],[-7.809974,-8.466151,11.56164],[-1.109558,-3.221675,3.5511],[-1.394479,4.363064,4.68838],[-7.972812,-7.28617,10.84684],[2.33563,-4.706736,5.348694],[8.471682,-0.9432539,8.582488],[-9.764753,9.03308,13.33968],[7.617754,0.1650439,7.684882],[-2.816867,1.282319,3.25255],[6.717492,1.944764,7.064475],[-5.94455,-3.437358,6.939244],[-2.209268,3.631077,4.366416],[7.674986,2.102244,8.020277],[-0.7320896,1.629951,2.047607],[5.566275,2.773414,6.298829],[3.752674,-1.642559,4.216701],[6.355454,-6.441626,9.104194],[-2.703317,2.753822,3.986409],[-8.461588,7.608008,11.4228],[-4.949049,-5.252151,7.285477],[-3.267402,2.830409,4.437018],[-1.765477,4.828528,5.237518],[2.320202,8.377325,8.750023],[-5.393049,-7.835573,9.564579],[0.7215912,7.913606,8.00911],[2.868868,9.791408,10.25193],[-2.165875,-0.6681592,2.477388],[0.280549,-7.765381,7.834529],[5.521624,5.142292,7.611274],[7.441223,-7.0288,10.28474],[1.396245,6.806418,7.019746],[-7.24221,-1.53345,7.470012],[-6.317887,1.530392,6.577066],[2.256967,-3.612664,4.375528],[-6.605617,2.942278,7.30008],[1.368516,2.078359,2.681868],[0.222781,-1.880353,2.141345],[3.418133,-3.938334,5.309813],[-9.742201,-1.973329,9.99022],[-7.380557,-0.007511387,7.447998],[-2.510121,3.229632,4.210847],[1.645018,-2.025418,2.794352],[-0.9620118,6.062205,6.218987],[-7.196429,7.901275,10.734],[-1.658927,9.484046,9.679832],[9.359571,-6.749017,11.58235],[-6.174525,-5.415334,8.273488],[5.876275,-4.347985,7.378047],[-4.443119,-3.842043,5.958406],[8.806245,4.577332,9.975065],[-9.998322,-8.189841,12.96302],[3.670466,6.390255,7.436913],[6.010252,-7.739885,9.850328],[-1.013413,-6.738624,6.887384],[1.47961,9.349658,9.518684],[-6.357447,0.3907436,6.447465],[-3.478887,-4.81395,6.02302],[7.11959,-7.525304,10.40763],[-7.686966,2.503084,8.14585],[0.8672729,1.379579,1.91191],[3.966667,-4.116413,5.803388],[-8.752355,-6.099651,10.71492],[9.705988,-4.489369,10.74061],[9.345048,-9.938581,13.67865],[6.325536,5.847777,8.672307],[-9.467117,-2.886945,9.947903],[0.550096,8.114392,8.194264],[-1.707835,7.55361,7.808567],[8.699325,2.255898,9.042528],[-0.443489,5.047714,5.164891],[-6.014341,-6.877272,9.190711],[-1.868284,-3.825175,4.372922],[-5.479156,-7.823905,9.603887],[3.293312,-7.579367,8.324224],[5.176101,8.492621,9.995831],[-1.088759,0.9390795,1.751361],[-0.4976019,7.390312,7.474244],[-4.986975,4.434894,6.7482],[-7.340069,-3.008807,7.995595],[-6.200695,5.339779,8.243898],[1.887079,4.416033,4.905346],[-9.834669,8.284545,12.89785],[8.00197,-0.4207813,8.075183],[-0.9788389,-8.19805,8.31662],[6.242282,-3.260337,7.113079],[-7.630072,4.633099,8.982405],[9.341702,6.547118,11.45129],[-8.921689,-3.456297,9.619903],[-9.970397,-9.732785,13.96911],[3.204634,-0.7150405,3.43234],[8.882448,-5.76064,10.63404],[-1.608783,-9.725355,9.908113],[-4.105629,8.825674,9.785127],[-9.784104,7.76367,12.53009],[6.189972,5.452042,8.309062],[5.159957,4.275333,6.775222],[5.335212,-5.345921,7.618619],[-2.142636,5.338643,5.838835],[-8.125752,1.054561,8.254693],[-7.319963,6.186022,9.635803],[2.862478,2.712537,4.068371],[-3.785566,3.611607,5.326745],[-6.445819,8.58559,10.78244],[0.7414798,4.000434,4.189662],[5.835517,8.341513,10.22908],[-3.802618,-9.4677,10.25169],[5.53199,-4.757628,7.364641],[-5.671235,-6.614383,8.770003],[2.979932,0.09383439,3.144646],[8.995069,7.873106,11.99571],[-4.044803,-4.981507,6.494293],[2.591845,-1.869556,3.348567],[6.07908,-2.206645,6.544043],[8.639713,7.967947,11.79546],[-0.7113718,2.960058,3.204371],[4.988286,8.573673,9.969497],[-4.301503,-7.557499,8.753212],[7.970456,-6.633963,10.41814],[3.144678,8.330637,8.960385],[-4.955563,7.297562,8.877614],[6.771559,0.4401241,6.859134],[0.2522663,7.438441,7.509597],[-5.741791,9.877705,11.46897],[-0.002834536,7.033933,7.104662],[5.292902,3.186608,6.258537],[-7.451675,7.161532,10.3834],[-6.796552,8.698745,11.08428],[-3.521904,-4.797263,6.034694],[-3.546559,-6.791486,7.72673],[-9.960214,-9.814773,14.01912],[-8.79664,-7.164227,11.3889],[-1.187238,3.865725,4.165737],[8.571535,8.423882,12.05956],[8.433756,-3.330383,9.122483],[1.280959,4.176704,4.481709],[4.859322,5.570313,7.459316],[-3.527965,9.010957,9.728508],[-3.513196,-8.335714,9.100916],[1.323209,-4.507433,4.802898],[6.291962,-6.677093,9.228887],[4.85391,-4.003003,6.370595],[8.097837,-1.426059,8.283031],[1.341087,5.405328,5.658276],[2.768852,0.3351713,2.962917],[-5.96591,-1.048702,6.139369],[9.898243,9.204383,13.55345],[9.064104,6.470388,11.18141],[-6.962765,6.457912,9.54907],[5.385591,-9.138425,10.65436],[-2.879481,-0.9181861,3.18347],[-6.333363,-1.869528,6.678819],[-3.519088,-8.846932,9.573515],[-6.49287,-6.636992,9.33847],[6.871681,3.178057,7.636756],[8.953896,1.719174,9.172122],[8.583587,7.548365,11.47414],[-1.834816,-1.007321,2.319751],[-3.467229,-7.182373,8.03792],[0.937993,-2.727779,3.052967],[-7.853124,5.494349,9.636359],[7.529009,-6.300434,9.868204],[0.4841012,-8.753901,8.824123],[-2.520921,0.7145643,2.804575],[-3.904924,7.536099,8.546416],[2.54444,-9.460559,9.847657],[-8.409608,-2.035116,8.709949],[6.222199,-1.144862,6.405191],[4.041777,1.261301,4.350499],[-8.574822,1.523805,8.766389],[2.530471,8.909701,9.315904],[-6.876236,-5.032562,8.579586],[2.266021,3.242768,4.08049],[-2.8887,-0.2234807,3.06505],[-8.580026,6.004681,10.52013],[-6.366804,8.897782,10.98666],[-6.307605,-7.48838,9.841835],[-2.533159,-7.415536,7.899815],[-4.17792,-2.795904,5.125631],[7.470123,-6.048381,9.663625],[8.491374,-2.163338,8.819493],[4.554606,-5.359749,7.10432],[-2.017992,5.530516,5.971508],[8.559461,-1.929476,8.831038],[5.414207,7.047172,8.942945],[8.284659,3.08646,8.897293],[2.302403,-0.822262,2.641434],[7.864484,2.507157,8.314803],[-0.8254211,-0.8056305,1.526552],[-9.931387,0.4388292,9.991247],[1.151211,-9.726581,9.845387],[-6.518139,-9.838907,11.84442],[-7.782867,3.737033,8.691285],[9.007822,0.5141098,9.077729],[5.715435,-6.499272,8.712447],[-9.568385,-4.646553,10.68384],[4.772173,-5.041345,7.013472],[-2.406569,0.2801538,2.62108],[-3.004109,-4.448497,5.460202],[-5.228158,2.356396,5.821189],[-0.8205496,-0.4607213,1.373159],[-5.745471,-9.670102,11.29253],[-0.5490607,-0.9196024,1.465311],[5.345141,-0.8636051,5.506028],[2.396254,4.997114,5.631446],[8.421698,8.006932,11.66344],[5.967124,-0.9656389,6.12691],[-5.409371,-3.648533,6.600992],[-8.68363,-4.08937,9.650304],[7.137286,8.129291,10.86399],[1.9712,1.366222,2.598498],[-2.578277,8.371614,8.816544],[-4.72879,-7.233075,8.699358],[-2.239137,0.01216827,2.452322],[8.055722,-6.173094,10.19812],[4.124795,9.563983,10.46345],[-1.23172,-8.347272,8.49671],[-0.2258907,-6.558723,6.638364],[-1.274325,9.18642,9.328142],[-0.6958172,2.813554,3.065982],[6.11798,6.792273,9.195904],[0.4199962,-4.095762,4.23694],[0.1417869,0.1486078,1.020876],[-0.7565235,-6.574657,6.693164],[9.174877,9.56471,13.29143],[-3.266567,-8.444529,9.109365],[8.591202,7.280721,11.30565],[-7.766782,9.029067,11.95186],[-5.801124,-7.201901,9.301635],[4.569283,-5.351764,7.107723],[3.064027,8.408289,9.004864],[-3.694114,9.42176,10.16937],[-2.080818,0.5668491,2.377209],[-5.16078,2.177069,5.689753],[-6.589136,0.4477324,6.679608],[-3.127911,-8.528039,9.13845],[6.72633,0.2068954,6.803405],[-9.025735,-7.433285,11.73532],[-5.992648,-5.940362,8.497043],[3.900394,2.075798,4.530123],[-4.484804,3.791527,5.957277],[6.87816,-7.937977,10.55086],[-8.687666,-4.535642,9.851274],[-4.642164,-1.216177,4.901915],[-0.3067089,-1.225388,1.611101],[8.144987,8.89753,12.104],[-0.4832915,2.515373,2.749668],[-2.640578,-9.576648,9.98423],[-1.668362,0.4760434,2.002511],[7.676082,-1.016946,7.807458],[9.945144,3.958249,10.75052],[-9.405484,4.90723,10.6557],[3.290848,4.205225,5.432642],[-6.399308,3.395197,7.3129],[3.164304,-5.593063,6.503475],[8.867822,-5.596076,10.53349],[1.406776,9.221146,9.381287],[7.260642,5.892334,9.404069],[7.271656,-2.523253,7.761687],[-6.643249,-2.441034,7.147825],[-6.578782,-6.067419,9.005218],[-2.061043,0.4447858,2.33361],[5.685795,4.593679,7.377679],[-1.990419,7.370585,7.699824],[3.070413,9.440065,9.977088],[6.792292,2.566762,7.329631],[8.354246,0.9878359,8.471673],[-5.455597,2.748026,6.189927],[-8.793015,-9.925571,13.2979],[-7.197706,-7.34424,10.33174],[3.53539,-5.099667,6.285347],[6.122721,-9.816615,11.61265],[-0.9693071,-8.163097,8.281045],[-3.307282,-3.000848,4.576374],[-8.529248,-4.069667,9.503172],[-9.003881,0.7178803,9.087642],[-8.770717,-8.065888,11.95759],[0.6704646,8.58518,8.669189],[0.05736151,6.087897,6.169747],[6.258059,-3.320621,7.154706],[1.693938,4.917612,5.296446],[-5.811575,0.1034274,5.89789],[-6.842583,-5.485727,8.8269],[-5.27635,-0.2817936,5.377665],[-7.793952,-9.187656,12.08961],[-1.203214,-1.553657,2.204898],[2.958686,-9.569278,10.06603],[-2.59116,-1.496461,3.154918],[-7.520459,-0.1013831,7.587331],[-1.677687,-3.469742,3.981676],[7.318183,6.005535,9.519572],[-2.09769,-1.101953,2.571887],[-9.4067,-2.138538,9.69842],[1.802655,-7.127885,7.419994],[9.49691,-7.296138,12.01769],[-9.876997,3.608595,10.563],[0.7077064,0.114424,1.230423],[-5.410233,7.970649,9.685137],[-6.787271,9.734004,11.90873],[6.980133,-2.177514,7.379961],[4.781888,8.557879,9.854124],[-6.413872,8.138393,10.41015],[-6.456966,-5.041533,8.252846],[-9.498833,-6.782908,11.71476],[6.444081,4.73064,8.056373],[-7.915591,-6.771502,10.46469],[5.640233,-1.542073,5.932134],[-5.978745,6.308198,8.748644],[-4.64818,8.366838,9.623385],[2.872935,5.838068,6.583069],[2.443709,8.345262,8.753005],[-4.299799,-9.224187,10.22614],[-4.526229,-4.795987,6.66995],[-6.525562,-9.486687,11.55769],[8.119173,4.447168,9.311191],[9.963889,-7.431947,12.47048],[6.885765,3.729117,7.894306],[-7.745547,2.848518,8.313095],[-7.346038,7.798335,10.76003],[-9.565331,-4.590172,10.6567],[-2.080207,7.853196,8.18535],[0.2762056,5.840543,5.931967],[3.553886,-6.193332,7.210234],[-4.710434,1.213523,4.965967],[9.096131,-2.248778,9.423195],[5.603229,8.463868,10.19967],[8.752498,-9.342013,12.84054],[6.805602,-8.234557,10.72959],[-2.354595,-9.895359,10.22068],[-6.195251,-9.395123,11.29821],[7.051396,-0.05801253,7.122187],[-1.51895,-7.029256,7.260693],[0.5710357,-1.165546,1.638469],[8.38316,-1.309111,8.543486],[-7.584769,0.003224481,7.650407],[9.627203,2.901395,10.10451],[8.701589,-1.372108,8.865683],[2.349599,1.105853,2.782719],[-1.844727,7.30598,7.601339],[4.51886,-4.733383,6.620047],[-6.152534,8.048641,10.18009],[9.609738,-1.191238,9.73479],[-2.510166,-9.905671,10.26758],[-6.032952,4.628039,7.669111],[-9.032045,0.6683925,9.111783],[-3.314516,2.49281,4.26616],[9.188288,6.56379,11.33613],[-8.598759,3.626505,9.385637],[-6.74834,-6.286643,9.276958],[8.863714,7.005636,11.34215],[7.068816,-7.068326,10.04636],[4.458429,6.11583,7.634197],[-9.900784,5.023381,11.14719],[2.112659,-9.1468,9.440724],[2.36831,-3.426184,4.283413],[-3.644206,6.476838,7.498644],[0.02813184,7.893926,7.957063],[2.212456,-2.410221,3.421129],[9.157416,8.498555,12.5333],[-8.138524,6.500087,10.46359],[0.4170667,-3.210383,3.388289],[5.338561,5.993906,8.088704],[-7.296755,-5.61996,9.264264],[-0.03984099,3.084425,3.242725],[7.906857,5.409262,9.632159],[5.696153,1.465044,5.965946],[-5.54676,-4.15445,7.001857],[-5.793713,2.831575,6.525712],[3.170079,4.859415,5.887556],[5.224588,8.51926,10.04361],[-2.23266,6.549489,6.991465],[-6.722463,4.312824,8.049345],[7.85836,-2.704596,8.370702],[4.342119,2.407593,5.064632],[9.256673,7.891332,12.20488],[-1.041666,8.17993,8.306402],[4.711588,0.3116349,4.826612],[6.007087,7.804272,9.899078],[-1.523333,-6.341169,6.5978],[-2.321968,1.540341,2.960437],[3.005346,9.485739,10.00057],[-7.468617,-7.321367,10.50632],[4.700648,0.27274,4.813572],[8.644441,5.807283,10.46188],[-5.224864,9.661751,11.02944],[1.240047,9.206403,9.343209],[2.696657,-5.490264,6.19798],[7.721027,-6.455823,10.11395],[9.555413,8.194212,12.62739],[-1.966033,4.338076,4.86664],[-6.780257,-5.631152,8.870274],[1.432619,7.431477,7.634085],[9.205764,-6.19849,11.14304],[2.320294,9.456417,9.788135],[1.695852,-4.86161,5.245109],[-7.060384,-3.076795,7.766318],[5.654979,1.144943,5.85574],[-5.860897,7.481451,9.556266],[-9.090671,-6.211082,11.05522],[0.9908937,7.421219,7.553566],[-3.273881,3.159265,4.658246],[-5.444359,-0.2433726,5.540783],[8.641319,-7.963031,11.79331],[6.721455,-4.73029,8.27971],[-9.714141,-6.128237,11.52909],[9.137045,9.32002,13.09001],[0.3965994,7.540221,7.616575],[-3.225683,-0.05041387,3.37751],[-3.879953,4.478992,6.009609],[-0.9086257,1.449317,1.981444],[-4.28878,0.008847909,4.403829],[3.574009,8.891496,9.634949],[6.323309,0.8350607,6.456126],[-4.902241,8.94894,10.25258],[7.557505,5.451907,9.372255],[-8.173478,5.778943,10.05992],[6.604413,7.981186,10.40757],[9.706312,-6.623049,11.7931],[7.680892,5.182306,9.319463],[-9.95693,-3.108499,10.4787],[5.753266,-2.907979,6.523528],[-3.399902,5.199487,6.292376],[-8.099285,-6.52687,10.44981],[0.6323549,-7.040657,7.139379],[6.007474,7.809853,9.903714],[-6.27541,-6.315988,8.959491],[-9.576747,-5.913208,11.29956],[-1.350085,-5.739144,5.980009],[4.715127,-6.869254,8.391607],[-4.700253,-6.650467,8.204944],[-1.052815,-0.05567932,1.453107],[-4.202512,9.994757,10.88835],[-2.460445,-9.521001,9.884495],[3.277577,3.238438,4.714869],[6.502559,-3.107726,7.276072],[-9.319098,7.485343,11.99483],[5.054916,-9.039697,10.4052],[8.651809,4.774879,9.932435],[7.919752,7.032094,10.63827],[-9.342142,-3.313599,9.962708],[-6.119645,2.802037,6.804518],[-1.285865,5.749808,5.976097],[7.072503,4.836334,8.626147],[4.502183,5.781779,7.395852],[5.135805,-7.734806,9.338293],[8.010149,-7.023561,10.70014],[1.985417,8.246566,8.540944],[5.811714,-3.973017,7.110617],[0.2955076,-6.817666,6.896948],[-6.923053,-1.116154,7.083393],[1.011084,-6.1286,6.291425],[-4.780539,2.177586,5.347469],[-3.246174,-9.668044,10.24738],[0.8294585,-1.384034,1.898302],[-1.287022,-4.059087,4.374084],[7.659674,1.756609,7.921886],[-3.774156,8.564518,9.412503],[-7.166386,-5.144437,8.878193],[2.89283,-0.4807405,3.098318],[-0.5939379,7.674165,7.761802],[2.779091,7.315125,7.888879],[-2.934112,7.74526,8.342546],[-7.141005,-0.1173069,7.211637],[5.711208,-7.698069,9.637332],[7.562585,-1.776493,7.832536],[3.515161,-0.08226947,3.65556],[0.04887685,-4.180752,4.298962],[-9.418813,1.211273,9.548885],[-5.144225,7.987236,9.552958],[-8.493486,8.921155,12.35825],[-9.794226,9.685065,13.81041],[-6.622498,-2.170884,7.040612],[1.84421,-2.980332,3.644652],[5.551663,1.629288,5.871588],[4.866793,-6.271788,8.001312],[-3.883588,-1.720004,4.363561],[1.633496,-2.305029,2.99691],[-4.791341,-5.337177,7.241713],[-0.153621,-5.659955,5.749668],[4.91619,3.88261,6.343783],[-8.752515,-6.508084,10.9527],[6.670581,3.683131,7.685188],[6.959377,6.153641,9.343459],[7.711674,-6.639959,10.22541],[8.664337,-8.999319,12.5323],[7.190095,6.67357,9.86073],[-7.714577,-3.813399,8.663527],[-8.843811,7.645998,11.73347],[-1.236038,2.33171,2.822174],[6.86543,-9.782062,11.99262],[-4.81969,-4.016256,6.352931],[-2.716782,0.9658375,3.051843],[6.279785,3.226591,7.13068],[5.988545,7.062302,9.313366],[0.3497157,3.259281,3.427129],[-0.400224,-2.957847,3.147862],[-2.771131,5.727732,6.440969],[-9.060739,-9.464213,13.14033],[-2.686362,4.666284,5.476381],[-5.426566,-4.923005,7.394836],[7.001914,2.398277,7.468502],[-1.27111,-9.327424,9.466602],[-6.433523,7.228132,9.728109],[8.146307,4.436872,9.32996],[-9.416797,0.4648163,9.481146],[2.115339,-8.598804,8.911459],[0.1275717,-5.450495,5.542939],[-9.390732,6.101947,11.24365],[9.173678,1.481514,9.346189],[-3.514993,5.934194,6.969206],[8.872467,-3.105265,9.453218],[4.635201,0.3042365,4.751594],[0.8521704,-2.50765,2.83099],[-1.264532,-0.5772447,1.712382],[9.594765,0.3381583,9.65266],[-8.372036,-8.990547,12.32562],[6.153417,-4.174575,7.502774],[0.05528133,6.724392,6.798567],[-7.295497,-9.276687,11.84404],[7.213897,-7.255341,10.28009],[-9.690878,-7.606582,12.36015],[-6.816602,-5.490444,8.809712],[9.828436,-1.108944,9.941222],[2.461898,8.214951,8.634023],[-8.031164,-4.560596,9.289705],[-2.749908,4.296051,5.197888],[-9.660964,3.596449,10.35706],[1.745205,2.79693,3.445077],[5.656064,-0.5055121,5.765987],[-5.517062,-7.139977,9.078394],[-2.987266,-7.487069,8.122806],[-3.665091,-7.776266,8.654664],[3.908836,-3.311644,5.219768],[-8.000709,-7.430827,10.96488],[-5.536983,4.382638,7.13202],[5.073999,-9.928991,11.1951],[5.979218,-7.368114,9.541497],[-9.741603,-4.551137,10.79869],[-5.692732,8.883421,10.59822],[-9.520691,4.359598,10.51901],[-4.118342,6.827909,8.036236],[2.721877,-5.980291,6.646239],[-3.570181,3.029938,4.788185],[6.733442,-7.02594,9.782795],[9.138887,6.516537,11.26874],[-2.102552,9.265126,9.553183],[-9.850201,2.459501,10.20175],[-1.372385,-1.455762,2.236668],[5.112638,-2.416369,5.74264],[-0.9729397,8.068336,8.188081],[-8.845244,4.332002,9.899726],[-4.36812,-3.785616,5.86612],[-8.291666,0.3297997,8.358259],[5.203055,-4.293825,6.81973],[7.585103,8.379109,11.34651],[9.685826,6.426669,11.66693],[5.68926,-4.481298,7.310931],[-8.260763,-0.7904697,8.358531],[5.930641,-5.287305,8.008002],[2.115996,8.640783,8.952127],[5.371409,8.283354,9.923003],[2.569029,9.425808,9.820681],[-1.364214,0.643064,1.809589],[8.484279,-2.065276,8.789104],[1.919727,-8.023396,8.310248],[-4.85901,1.329636,5.135943],[5.782844,3.782283,6.981902],[6.844118,-4.294352,8.141462],[9.725384,7.933613,12.59068],[6.427474,0.776082,6.550933],[4.73378,-2.455393,5.425645],[-1.155095,0.9355804,1.791523],[-4.330986,6.633783,7.985268],[0.3120895,4.209913,4.338291],[-7.870994,3.749007,8.775397],[-0.5163726,2.609716,2.842052],[-1.974608,9.29828,9.558089],[-7.22308,-3.424016,8.055853],[8.58983,-3.005651,9.155278],[-6.137481,-8.10546,10.21602],[0.1509345,4.450158,4.563627],[-6.348761,-6.90885,9.436047],[-9.152072,5.700014,10.82823],[-0.7172685,3.248672,3.473952],[2.923665,4.441242,5.410402],[-2.979861,7.959851,8.557966],[-3.394727,-8.39138,9.107109],[2.347251,-1.456962,2.938082],[5.05723,6.286762,8.130126],[8.105049,9.156861,12.26947],[-6.956933,8.350556,10.9147],[4.714888,-5.170071,7.068224],[-9.30477,1.944593,9.558252],[-4.266386,-9.944015,10.86671],[-3.496035,-6.184555,7.174328],[-5.242507,5.198597,7.450456],[4.374572,0.4114265,4.506235],[-0.5373961,3.77066,3.937852],[5.224508,1.003154,5.413114],[-0.2638269,7.755358,7.824013],[6.042891,2.299481,6.542488],[6.961604,-2.098565,7.339475],[-5.976563,9.71031,11.44594],[1.531032,4.352366,4.720927],[-0.4361767,3.635112,3.795298],[6.651188,-8.053465,10.49269],[-2.094604,6.966311,7.342809],[8.144788,-7.383139,11.03849],[-0.4258833,4.337403,4.471515],[-7.691981,-6.565997,10.16262],[5.349708,2.027448,5.807746],[-9.986711,6.215974,11.80562],[1.390699,-3.992712,4.344628],[-9.976563,7.714045,12.65062],[6.307154,9.084847,11.10471],[-3.561745,-7.074172,7.983103],[3.837214,-4.057848,5.673654],[5.75317,2.948879,6.541778],[5.303426,9.384875,10.826],[-8.555077,2.325507,8.921733],[-5.507264,-3.052977,6.375784],[0.6302136,-0.628988,1.338953],[6.883108,-0.2492657,6.959835],[-7.661137,-0.3497174,7.734037],[6.809366,2.122405,7.202226],[-8.270583,-8.308651,11.76589],[-1.14907,0.266005,1.546325],[-1.331891,-2.732104,3.199738],[2.674235,5.148112,5.886815],[-5.435168,6.1555,8.272318],[4.759968,6.239854,7.911578],[9.939418,-7.696709,12.61076],[-5.412753,-3.040012,6.28805],[-8.404799,-0.4041634,8.473723],[-7.941461,9.65003,12.53754],[2.161283,3.524569,4.253672],[-2.584022,2.257523,3.574014],[3.719075,0.4108611,3.873025],[9.397729,-5.395422,10.88246],[5.493707,-3.001473,6.339531],[7.834795,-5.247403,9.482576],[-4.657112,-1.908442,5.131359],[-8.869919,3.889807,9.73684],[-2.898595,4.3419,5.315445],[8.609422,0.5180905,8.682774],[-2.119409,4.600142,5.162673],[1.427867,-9.888586,10.04106],[2.845344,-0.3770888,3.039437],[5.380958,-2.146908,5.879109],[0.2785438,5.524589,5.621269],[8.36851,8.07515,11.67219],[-5.47829,-3.0914,6.369334],[7.561787,-3.172675,8.261144],[9.905645,-2.393294,10.23961],[-8.689283,-1.221714,8.831547],[-2.009772,5.321614,5.775704],[3.303096,-2.229831,4.108842],[-3.00502,7.089302,7.764557],[3.832423,1.785506,4.344594],[-1.22665,-0.7005389,1.73073],[3.061147,-4.515763,5.546416],[0.2950814,9.956347,10.01079],[2.342495,-5.496477,6.057932],[9.704756,3.116942,10.24195],[5.292198,1.735932,5.658695],[2.588292,8.449986,8.893904],[-1.937122,4.032635,4.584167],[-1.629061,1.446221,2.396955],[4.388939,1.226719,4.665579],[-5.912254,-0.5252813,6.019192],[9.38886,8.074286,12.42356],[-3.721497,-3.53158,5.227007],[6.687713,4.050403,7.882339],[9.886271,2.925684,10.35847],[-7.33911,7.526112,10.55959],[-8.687601,-1.454512,8.865101],[5.457984,-0.2595952,5.554906],[-3.38138,7.87986,8.63284],[9.314254,-7.614414,12.07206],[9.717654,5.841908,11.38247],[-0.8098171,4.376632,4.561876],[7.515048,-1.390892,7.707822],[7.611284,-7.851686,10.98092],[3.060411,1.466578,3.537933],[-3.747635,-1.522596,4.166901],[-9.771288,9.061381,13.36363],[-9.713408,3.138533,10.25674],[-2.506113,-1.881171,3.289287],[-7.033851,-4.093378,8.199439],[8.790252,3.079035,9.367442],[9.582214,-1.590485,9.764654],[-4.51204,-0.6480803,4.666746],[2.913148,-8.22659,8.78426],[3.073523,2.457875,4.060504],[6.849484,3.938343,7.964043],[-2.189016,-8.700482,9.027191],[-2.726192,-8.054372,8.561835],[-7.967922,7.574862,11.03931],[-5.646032,8.481114,10.23753],[0.1006383,2.149233,2.372621],[-9.330975,8.562786,12.70387],[-3.647635,-7.956758,8.809951],[4.140971,6.679741,7.922536],[0.9725227,-1.204681,1.84311],[-2.875521,-7.934724,8.498734],[3.659951,6.148551,7.224951],[6.285611,0.6529648,6.398067],[4.303754,5.399015,6.976508],[-7.389679,-5.060134,9.011787],[-7.489662,-9.477468,12.12095],[8.218142,7.660279,11.27908],[3.308093,-5.689506,6.656873],[6.800603,2.850459,7.441325],[-8.042593,-7.117435,10.78616],[9.796212,-0.2718988,9.850873],[-9.643706,4.984173,10.90152],[-5.091618,4.836809,7.093609],[-6.061507,1.343542,6.288638],[8.900368,8.189474,12.13606],[6.420856,-5.281527,8.373882],[9.673691,-8.02802,12.61069],[-9.236411,8.962557,12.90886],[6.111218,-2.311999,6.610017],[4.618452,-1.548904,4.972846],[5.262911,8.928064,10.41194],[-0.712313,-2.685911,2.95322],[-7.49828,-6.241155,9.806948],[5.803868,2.012737,6.223825],[8.520808,6.568493,10.80506],[-3.79365,6.086175,7.241085],[-0.2100788,-0.906535,1.365994],[7.461927,3.229228,8.191963],[-4.014753,3.97233,5.735647],[-5.488386,-9.542683,11.05374],[6.355106,6.217136,8.946517],[-8.47031,7.459631,11.33103],[3.687391,-3.181161,4.971583],[1.00514,-3.314373,3.60491],[0.4590615,3.746172,3.904426],[6.354249,-5.435336,8.421363],[-2.72963,-8.211126,8.710539],[3.145156,-0.9194884,3.425999],[-3.422776,-6.361998,7.293175],[6.83184,9.609143,11.83257],[-5.112109,-1.922883,5.552579],[-5.737499,2.370407,6.287903],[-8.785292,-4.460123,9.903234],[1.917969,-2.711665,3.468679],[-4.034719,-9.761673,10.60986],[-0.704213,7.331947,7.43326],[-3.407119,-0.3346226,3.566571],[-5.843982,7.806048,9.802373],[5.677042,4.569586,7.355944],[9.42035,4.769665,10.60626],[6.611796,3.854523,7.718367],[4.871259,-1.757893,5.274405],[-8.599042,4.399346,9.710704],[8.596562,4.167639,9.605733],[0.7838565,-2.431814,2.743748],[-7.997217,-8.364579,11.61558],[-1.952592,1.948595,2.934219],[-7.393856,-9.232882,11.87077],[-3.896001,-6.842731,7.937366],[-9.529223,5.078297,10.84413],[2.598686,-0.988889,2.954839],[-2.177667,-8.75069,9.072861],[2.466465,0.4048691,2.692093],[-2.003557,-4.448223,4.980053],[4.968081,-8.890176,10.23314],[7.65771,-4.287101,8.832879],[-6.867488,8.425926,10.91598],[0.400919,-4.90643,5.023325],[-8.545948,2.604914,8.989928],[4.081016,7.779835,8.841975],[-0.8879083,-0.1065813,1.341544],[-9.426552,5.551732,10.98552],[-2.155964,2.346442,3.339756],[1.036439,2.721967,3.079498],[1.53487,-7.005282,7.240843],[-4.863965,-8.549818,9.887241],[2.899323,2.800843,4.153408],[4.84918,2.371097,5.489686],[-3.035087,-0.9401883,3.331022],[5.219366,2.129823,5.725201],[-3.275616,2.877387,4.473144],[1.346292,9.667265,9.811652],[9.499433,-8.969755,13.10327],[-1.215968,2.830204,3.238616],[-7.64644,3.633471,8.524679],[5.66618,0.6697079,5.79259],[5.546146,-9.375782,10.93915],[7.142503,-1.368971,7.340942],[8.766945,-9.774076,13.16783],[-9.867897,-7.895108,12.67707],[-2.833221,2.83088,4.128078],[-6.134891,9.111846,11.03008],[-0.06419513,9.07541,9.130564],[-9.295587,-2.840542,9.771214],[0.495623,2.029551,2.316186],[-6.160142,8.467386,10.51874],[4.948281,3.558643,6.176522],[2.179674,-5.345382,5.858676],[-5.873667,-2.838712,6.599867],[8.115447,-1.417761,8.298826],[-8.273927,-5.048388,9.743925],[2.562944,3.539474,4.482919],[-8.172512,1.985257,8.469427],[6.368373,0.436836,6.461192],[-0.8734065,-0.2534992,1.351703],[7.623536,-7.200106,10.53375],[-1.199864,6.344592,6.534028],[-2.314285,8.35978,8.731658],[-1.854185,-3.868324,4.404763],[-5.225898,1.703422,5.58674],[8.162394,7.273822,10.97876],[-4.53942,-4.512261,6.478181],[-6.27558,-6.751106,9.27148],[-6.430559,-1.377687,6.652075],[8.490353,6.328948,10.63681],[-9.783928,0.1828912,9.8366],[-0.4069764,-3.467886,3.63206],[0.9609663,6.996747,7.132876],[5.680158,4.494083,7.311702],[7.81978,7.879951,11.14642],[3.696755,0.2708317,3.839186],[6.753907,-6.788193,9.627815],[1.644007,-1.456293,2.413203],[-0.4227698,9.492949,9.554832],[1.54039,-2.564134,3.153979],[0.6216,-3.085399,3.302434],[0.4085452,-5.1043,5.217354],[-8.58848,3.80134,9.44522],[-3.201164,5.323492,6.291821],[3.609313,6.981125,7.922327],[-9.920141,-8.02449,12.7985],[-1.795236,5.646761,6.009058],[-7.231731,-2.202168,7.62545],[-7.576827,-7.799958,10.92006],[-4.150321,-1.925343,4.683173],[3.054161,1.029894,3.374697],[2.167949,4.965298,5.509463],[8.492682,-7.338418,11.26845],[4.036742,8.873575,9.799776],[1.238819,-2.10774,2.641447],[5.53414,1.267627,5.764858],[5.209553,2.438741,5.838399],[-5.804586,-4.038043,7.141358],[3.302181,7.341842,8.112154],[0.4445978,-6.78814,6.875792],[-8.180478,-4.373916,9.330132],[9.765917,-2.499104,10.13009],[3.551634,-2.781557,4.620732],[-0.9520357,-5.446328,5.618617],[2.539263,0.6036174,2.795033],[-2.575618,-4.461639,5.24786],[3.921132,-3.903223,5.622315],[0.08766359,3.615979,3.752731],[6.45107,-0.2302451,6.532176],[-8.230626,-1.789837,8.482141],[1.486085,3.349394,3.798274],[1.458294,0.6839184,1.895881],[-8.351471,-5.954769,10.30564],[-1.242331,6.992116,7.171686],[-8.512174,2.379429,8.894874],[-5.539298,-2.208616,6.046636],[9.883962,1.407377,10.03361],[-2.948647,9.470695,9.969381],[-2.004489,8.043962,8.350048],[4.317366,4.297235,6.172996],[0.7603397,5.014472,5.169434],[-1.03705,2.728283,3.085288],[-9.839784,2.771929,10.27156],[-5.109217,5.539526,7.602003],[1.258254,2.065607,2.617238],[-7.590696,-7.505877,10.72179],[-4.392595,0.8461297,4.583757],[4.731454,-4.357172,6.509348],[4.531736,-3.610459,5.879799],[9.422297,-2.922859,9.915784],[6.549312,-6.505973,9.285536],[-8.341627,2.976743,8.913122],[3.28631,-6.167225,7.059355],[1.212044,7.902057,8.05677],[4.718584,-3.008519,5.684736],[2.576645,5.792448,6.418065],[3.622828,-2.033818,4.273324],[9.028528,0.02504543,9.083775],[9.457362,-6.049624,11.27119],[-9.7353,-8.349508,12.86431],[4.581546,1.633982,4.96593],[-3.218571,-4.515484,5.634607],[-8.826765,-5.317999,10.3534],[7.906235,7.640223,11.04],[-8.593467,-4.036452,9.54676],[-6.583915,-8.07842,10.46942],[0.01355066,8.782938,8.839693],[7.725895,-0.0001989026,7.790344],[2.526086,5.634169,6.254995],[-9.814339,-5.97842,11.53528],[-3.017777,-7.206866,7.876922],[-4.998536,-7.955591,9.44864],[-4.89253,1.426206,5.193353],[7.91209,-1.596901,8.133343],[-7.134211,5.966797,9.354124],[-7.623719,-9.315269,12.07871],[7.411042,3.109991,8.099111],[-9.956186,0.9988199,10.05601],[3.685178,-6.114183,7.20859],[-7.470538,2.709942,8.00954],[-8.710026,1.586316,8.909598],[-3.242844,2.730262,4.355498],[8.37985,-9.908996,13.01576],[-7.993836,-0.0338022,8.056212],[2.334594,8.231142,8.614059],[-1.287272,-3.404085,3.774237],[-6.653809,-0.7597569,6.771292],[-2.318077,0.1299019,2.527915],[-3.566018,-3.487973,5.087479],[5.030667,-8.7639,10.15448],[-0.4099193,2.558668,2.777556],[-1.607957,0.4748652,1.952184],[1.13855,8.227063,8.365457],[2.278853,8.162086,8.533043],[-8.028608,-3.276806,8.729033],[-4.942138,-8.991689,10.30899],[4.583262,0.4444767,4.712097],[-6.995543,-2.365533,7.452072],[-8.738464,0.4322799,8.806113],[3.36752,-0.4778393,3.545211],[-4.933935,-5.145176,7.198371],[1.661434,6.284175,6.576566],[2.247256,2.871178,3.780717],[8.341209,-1.446343,8.524534],[8.640924,-8.922264,12.46083],[9.148058,9.873708,13.4973],[7.048779,-1.982616,7.390268],[-0.4077353,-0.6749324,1.273492],[-3.441042,7.156018,8.003085],[7.109136,-5.202234,8.865836],[-8.51216,5.200642,10.02515],[0.110187,-8.898687,8.955378],[8.19046,-2.312141,8.569109],[-0.7714875,8.794306,8.884538],[5.291195,2.369772,5.883244],[2.290035,-3.823735,4.567846],[8.365217,-3.792887,9.239202],[3.914841,-0.9055041,4.140763],[0.7836291,-4.072394,4.265967],[-5.321082,5.440382,7.675393],[7.224546,9.662083,12.10578],[-6.269326,2.559697,6.84518],[5.010461,-3.619702,6.261546],[-8.167981,-3.235171,8.842072],[-6.41914,-6.093175,8.906859],[-1.836858,-3.862442,4.392323],[-6.693031,-6.203371,9.180331],[6.487657,8.567205,10.7929],[1.227096,1.886894,2.462952],[-5.739096,2.43817,6.315211],[0.05224429,-9.298409,9.352172],[6.49795,3.007214,7.22957],[-2.603885,6.166858,6.768334],[-6.383596,-3.87799,7.535855],[-8.229796,-1.752917,8.473621],[-1.121638,6.225115,6.403915],[-4.676883,-5.206608,7.069795],[8.607748,4.204974,9.631986],[1.449065,-3.403363,3.831798],[5.60028,2.50017,6.214015],[9.470739,-7.789732,12.30345],[6.275753,-1.104111,6.450127],[-3.176027,-3.635947,4.93024],[7.188738,9.777526,12.17694],[-8.702328,-3.398479,9.395752],[-0.4473366,-0.2466115,1.12291],[9.868558,5.265122,11.22987],[-7.603662,5.264144,9.301983],[1.539964,7.210263,7.440388],[4.820233,9.385746,10.59844],[-1.083886,8.240823,8.371737],[3.130792,8.966063,9.549459],[7.783566,3.700001,8.676053],[-3.753003,-9.901081,10.63562],[5.966269,1.42789,6.215724],[4.992762,-5.363596,7.395663],[-8.597426,-9.813262,13.08495],[-8.804231,-5.534162,10.44708],[-9.456139,-8.532581,12.7759],[-4.848384,-9.449927,10.66808],[3.437904,-1.69978,3.963387],[-7.806082,3.604158,8.655915],[-7.816789,1.291956,7.985696],[2.09898,5.968791,6.405636],[9.902542,5.849402,11.54452],[-0.7812575,3.470397,3.695135],[6.131123,-1.289393,6.344541],[-2.280894,0.4838836,2.537049],[1.392283,0.3589211,1.751364],[5.262171,7.133982,8.920995],[-4.395662,-7.302435,8.581806],[3.679153,5.143423,6.402419],[-6.922775,-3.049582,7.630516],[-7.011956,6.381824,9.533897],[2.749953,-1.380092,3.235258],[-3.169292,-2.804842,4.348741],[9.347296,-1.025422,9.456396],[5.317281,-8.352743,9.951974],[7.877711,-5.106984,9.441378],[-6.817086,1.664488,7.088242],[1.221074,-8.525931,8.670786],[5.847062,-1.533443,6.126956],[-7.332918,-1.228738,7.502099],[-9.506182,-7.040963,11.87193],[0.8802206,-6.052578,6.197458],[-7.718095,5.147407,9.330852],[6.239284,2.275912,6.716282],[2.088439,-3.296391,4.028371],[-4.885619,1.261187,5.143915],[2.29446,-5.161702,5.736525],[-4.315888,-0.427763,4.450828],[5.178983,7.132309,8.870833],[6.416711,-8.093504,10.37685],[-4.668981,7.766009,9.116485],[-6.567901,-6.958319,9.620578],[-5.68936,6.735595,8.873391],[-7.544181,-5.604098,9.450956],[2.374298,-9.859146,10.19019],[3.701605,7.812054,8.702302],[0.4328016,8.734986,8.802688],[-2.707819,9.823438,10.23876],[-7.015729,1.687685,7.28483],[7.128049,-4.114085,8.290645],[-4.506529,2.17575,5.103204],[-2.176206,6.114533,6.56684],[-8.853323,0.4589798,8.921434],[8.761862,-5.620487,10.45754],[-7.432215,5.14282,9.093207],[-7.975033,6.203743,10.1532],[-2.271354,7.525535,7.924186],[-7.913134,-6.978258,10.59782],[8.853261,-2.216847,9.181211],[-4.262769,1.346763,4.580936],[7.881732,3.646243,8.741669],[3.767196,-5.796587,6.985141],[9.663768,-3.231267,10.23863],[-0.838113,5.440815,5.595078],[0.01758915,9.283054,9.336777],[7.18789,4.472132,8.52442],[-4.496712,-5.753802,7.370662],[0.4737675,3.855541,4.011191],[-6.395357,-3.903155,7.558784],[-2.071427,0.2525704,2.314002],[-0.583937,-1.922532,2.244351],[6.40242,-7.378621,9.820134],[-9.901577,4.32034,10.84927],[-9.15307,8.737449,12.69337],[-9.562348,7.469228,12.17489],[-8.502807,5.833885,10.36011],[-2.103853,1.51075,2.77643],[0.6278409,6.974176,7.073423],[-3.76777,8.062303,8.955267],[1.628041,-0.4235212,1.957009],[1.834822,-1.520167,2.584082],[0.8803871,8.174051,8.28192],[-2.864612,3.478539,4.615868],[8.502811,-2.216656,8.843719],[-4.292207,5.367061,6.944665],[3.665708,2.218157,4.399731],[-3.282311,-3.815122,5.131152],[-7.334166,8.857468,11.54317],[-2.108232,7.612644,7.962223],[-7.679197,-2.714946,8.206156],[6.280956,2.978266,7.022854],[8.575798,4.261388,9.628278],[-1.176049,7.865938,8.015988],[-6.910703,9.303584,11.63248],[6.801386,5.664369,8.907521],[-0.1826481,-4.544296,4.656607],[-0.6958897,-8.236432,8.326048],[-6.372269,5.869689,8.721184],[1.164694,3.561532,3.878275],[-6.570951,-2.050509,6.955717],[3.26231,5.26378,6.272961],[9.028268,-1.615678,9.226052],[-4.494044,-4.23007,6.252193],[8.239147,-9.330477,12.48765],[-7.753506,-0.6346967,7.84345],[9.594437,-1.063586,9.704866],[8.27592,-2.296895,8.646767],[-1.527382,-7.628835,7.844235],[3.135595,6.379183,7.178156],[1.281318,1.168977,2.00207],[3.354242,5.211828,6.278064],[6.797619,1.093817,6.957303],[-3.8277,-6.451684,7.568058],[-4.386517,9.210752,10.25083],[-5.361558,-2.414589,5.964608],[4.077808,1.411443,4.429525],[-0.3566394,8.445326,8.511799],[0.4379539,-3.192174,3.373689],[2.614228,1.382634,3.121837],[-4.522262,2.579446,5.301357],[4.852912,9.44384,10.66475],[-8.778483,1.609621,8.980682],[0.7209939,2.097156,2.432672],[-1.352192,-9.621225,9.767108],[2.43314,0.1913206,2.63757],[5.420979,-7.458278,9.274316],[-4.339304,-0.6493874,4.500141],[-6.930082,4.164631,8.146791],[-7.389629,-6.509341,9.898391],[-6.570888,-5.396851,8.561691],[2.457463,6.661322,7.170239],[6.415688,0.9102389,6.556644],[6.140507,-4.608752,7.742508],[7.654753,7.662013,10.87666],[8.578066,-8.720964,12.27349],[-1.410535,9.227992,9.38858],[-6.033105,-5.634171,8.315181],[8.878057,1.694841,9.093535],[2.216245,-6.124629,6.5896],[-8.759223,2.822046,9.256778],[-8.151454,-1.765532,8.400197],[-4.614935,-4.416363,6.465438],[0.580506,-2.147643,2.439131],[-7.554883,9.255985,11.98956],[8.039597,6.923586,10.65698],[1.294822,6.451134,6.655351],[0.5737064,1.917551,2.237441],[-0.2189846,7.62392,7.692341],[7.815151,4.466852,9.057006],[-2.230784,0.767038,2.562176],[1.735363,-1.973717,2.811947],[-0.8918424,2.172843,2.552769],[6.036754,-8.545058,10.51001],[5.601614,2.373123,6.165208],[-7.324876,1.060416,7.468486],[-1.284686,7.063678,7.24886],[-6.125092,6.264773,8.818398],[-3.673465,5.753618,6.899164],[-1.094023,-9.221831,9.340185],[7.110207,3.60023,8.032228],[-4.654252,-7.49536,8.879329],[-5.158544,9.769765,11.09319],[-6.350033,-2.80598,7.014018],[-1.342552,-6.484852,6.697443],[4.64857,9.549127,10.66747],[2.224705,5.961815,6.441471],[-8.019834,0.1616334,8.083555],[-5.288442,-7.935313,9.588368],[0.9381549,4.47694,4.682214],[2.146832,6.086014,6.530579],[-6.820617,8.808073,11.18494],[1.073914,5.406758,5.60235],[7.458989,3.638862,8.359296],[-3.144139,-1.922775,3.818727],[5.676257,2.776208,6.397439],[-1.884646,-6.292674,6.644519],[2.678828,-8.324955,8.802329],[-6.694368,4.006064,7.865311],[-8.935085,0.57959,9.009532],[-2.037726,2.866172,3.656127],[-1.378841,-3.47156,3.866902],[7.195556,-1.368278,7.392443],[9.189773,-1.452697,9.35747],[-2.250911,3.380278,4.182449],[9.742982,0.1139229,9.794829],[-2.002393,-1.46029,2.672456],[5.452972,-3.386751,6.496537],[-7.601783,3.563102,8.45475],[0.8048021,8.852252,8.944836],[-6.482598,-3.680233,7.521183],[9.397771,-7.175802,11.86635],[-9.992327,9.915764,14.11272],[2.761565,2.253851,3.702173],[1.579102,-2.120915,2.826985],[-6.466414,7.706761,10.10983],[6.286707,-0.967658,6.43887],[-9.807929,3.922649,10.6105],[0.1588592,-9.758213,9.810605],[0.9871626,7.163414,7.299931],[-5.67415,8.944867,10.63986],[-8.878024,2.234915,9.20946],[1.425817,9.028227,9.194663],[5.05016,-6.314926,8.14754],[-0.01506646,9.522015,9.574392],[7.66149,-9.332746,12.11605],[9.517308,-2.339504,9.85152],[3.339201,-7.324954,8.112041],[0.3698737,-9.019782,9.082582],[8.68538,-6.043775,10.62841],[-7.102218,-1.712586,7.373903],[-3.639417,8.757776,9.536456],[-3.060217,3.024961,4.417614],[6.136308,-3.579839,7.174226],[-6.271366,-5.309749,8.27789],[-0.02308518,7.861613,7.924992],[3.804981,2.848308,4.857029],[-4.973302,4.514073,6.790478],[-6.193156,9.112843,11.06341],[-3.326665,7.388991,8.164796],[9.344018,-9.423883,13.30865],[-1.189527,4.666159,4.918131],[6.26242,-5.895755,8.658974],[-7.232396,5.33037,9.039933],[-9.910393,-5.096123,11.18867],[8.586522,-2.835258,9.09764],[7.568013,9.056513,11.84463],[9.585442,6.698571,11.73676],[4.386148,-7.281508,8.559127],[0.1045197,-0.3170305,1.054245],[6.643882,-1.7954,6.954468],[6.332513,-8.384063,10.5543],[-7.816001,-1.40168,8.00341],[-3.013778,-5.9563,6.749843],[4.089204,5.651335,7.046926],[4.808768,-5.773274,7.579904],[-0.1417471,-6.633094,6.709548],[9.291702,9.527193,13.34553],[5.056002,-4.849937,7.077079],[-4.577136,8.676988,9.861049],[-1.037334,3.789022,4.053733],[-1.030892,5.166046,5.361975],[-8.145175,-1.940997,8.432755],[6.990674,-8.781493,11.26872],[-6.316634,5.520186,8.448214],[8.731665,-3.539948,9.474873],[-6.39118,0.858331,6.525635],[0.9387794,-8.334926,8.447028],[9.591467,-6.420211,11.58513],[5.207167,6.152712,8.12222],[-9.841431,-0.8555321,9.929033],[1.319028,8.397737,8.559311],[-9.530882,8.586939,12.86753],[8.269636,-0.8072497,8.368902],[3.162071,-7.132278,7.865627],[9.175628,7.813777,12.09327],[4.769326,7.845448,9.235666],[7.822712,-3.488331,8.623415],[0.4605006,6.694979,6.784895],[8.513551,-1.144994,8.648211],[1.738642,-7.085684,7.364088],[-3.055666,5.637854,6.490184],[2.717718,-7.909828,8.423264],[5.122303,-9.203934,10.58066],[2.243475,-7.662331,8.046397],[-8.278835,-7.205778,11.021],[-4.389565,-9.346236,10.37403],[-5.575451,0.1534223,5.666497],[7.921429,-9.041698,12.06239],[2.471972,-3.363847,4.292564],[9.804796,8.126839,12.77417],[-5.766611,-2.140164,6.231702],[2.815241,4.238675,5.185744],[1.195815,-2.386564,2.850555],[2.328872,4.043753,4.772377],[3.174418,-3.627871,4.923248],[1.194175,-5.215586,5.443196],[7.99209,-9.440923,12.40986],[7.33506,-1.372923,7.529145],[-3.875699,9.138737,9.976851],[-1.093602,-0.02147418,1.482035],[0.5499411,-8.830357,8.903798],[7.357243,3.140731,8.061836],[-9.647991,5.660707,11.23064],[2.167501,-0.9435672,2.566784],[1.973622,-9.091444,9.356791],[-0.5766166,6.523121,6.62447],[-3.645449,3.816537,5.371709],[3.028321,-3.199838,4.517708],[-4.259498,-2.37993,4.980702],[1.573673,-3.574972,4.031981],[9.928385,1.432314,10.08089],[5.668653,-4.023506,7.022978],[-3.015524,-6.736376,7.447963],[8.741213,-0.08338013,8.798622],[-4.664971,-3.820917,6.112394],[-8.833728,8.225812,12.11193],[-9.139448,0.8278351,9.231187],[8.898272,-0.04355759,8.954392],[9.644762,9.09332,13.29323],[-3.540021,9.289449,9.991277],[5.201117,-7.782029,9.413373],[-8.826631,-8.85052,12.53958],[5.789567,-2.985605,6.590365],[2.089862,6.618273,7.012065],[3.846438,-1.821747,4.371939],[8.172541,6.143625,10.273],[8.927029,4.608509,10.09605],[-4.411568,-4.609733,6.458449],[-5.735736,-8.023814,9.91364],[8.031073,6.148355,10.16368],[5.306867,-1.68213,5.656183],[2.992514,1.44649,3.470947],[-9.735484,9.176992,13.41629],[-8.850456,3.223975,9.472307],[5.359424,-4.959246,7.370044],[-6.850361,-7.549759,10.24335],[1.390629,-6.773109,6.986334],[6.527538,-1.252648,6.721449],[7.724956,-0.6646653,7.817719],[-9.097095,5.674687,10.76844],[-7.881074,-0.4224228,7.955487],[-8.019501,1.667505,8.251846],[0.6658719,5.052963,5.193826],[-3.719379,7.9034,8.791901],[9.648835,-7.592917,12.31878],[-8.421021,3.082421,9.023022],[9.322601,-7.927401,12.27822],[2.746874,4.97211,5.767772],[-6.321476,-6.113978,8.85109],[-8.943552,-7.790077,11.90262],[3.735515,-7.85208,8.75267],[-2.271551,4.570876,5.201235],[-5.290196,-9.966021,11.3273],[-8.550871,3.290435,9.216527],[-4.704498,8.763965,9.996968],[4.555178,7.043856,8.447813],[-1.393412,-6.027716,6.266973],[6.197787,-2.924438,6.92567],[-7.470422,2.070308,7.816226],[-9.002561,-5.003296,10.3479],[4.021904,-8.849232,9.771623],[-4.371559,-0.5781825,4.521595],[-1.394048,7.005796,7.212805],[5.151854,5.657102,7.716502],[0.783473,-5.586363,5.728987],[-3.893329,6.327384,7.496252],[-4.106097,-0.03957271,4.226298],[-0.6583478,-4.874676,5.019551],[-4.296177,-3.784023,5.81171],[2.435812,1.556926,3.058954],[-2.280282,-3.362389,4.183939],[5.392343,-6.210691,8.285532],[3.902527,2.041545,4.516373],[-6.583261,4.775305,8.194076],[-6.536058,-6.967498,9.605524],[-8.436217,6.876896,10.92984],[-3.190295,6.07799,6.936854],[-1.69324,9.154835,9.363656],[-6.37684,-8.837935,10.94409],[3.026399,1.806178,3.663519],[-4.392001,-7.701091,8.921685],[9.019797,-6.558826,11.1971],[-7.64614,-7.250592,10.58464],[-8.776855,-9.227051,12.77387],[8.761072,2.539946,9.176476],[4.766675,9.24598,10.45033],[7.529822,-4.864341,9.019979],[3.063493,-8.498067,9.088572],[2.558691,-5.823748,6.439172],[-5.859529,-0.5744239,5.971939],[-0.108496,8.145699,8.207568],[4.801093,-2.700113,5.598312],[-0.3307176,-0.9465841,1.41612],[9.975103,-5.641694,11.50354],[-6.956714,-1.798866,7.254777],[3.258018,-4.504288,5.6483],[-8.936925,-8.976037,12.70582],[-8.590981,-6.434268,10.77983],[7.084657,-6.528081,9.685464],[-5.187905,2.153354,5.705374],[-0.977889,-7.181745,7.316674],[-6.881467,4.318819,8.185767],[8.238012,-7.040447,10.88268],[2.254118,9.661015,9.97077],[-3.607501,-0.7381164,3.81561],[8.605948,4.560483,9.79083],[-5.889764,4.122805,7.258571],[9.474532,-1.491558,9.64321],[-8.784048,-2.765752,9.263308],[-0.931716,-4.695389,4.890274],[-5.544023,5.933953,8.182175],[-7.793084,1.441166,7.988061],[-5.456663,-8.520458,10.16727],[-4.4685,7.284647,8.604278],[-0.4957681,0.5219213,1.232148],[8.088213,-6.67859,10.53673],[7.223837,-9.917446,12.31014],[-4.832296,-5.756984,7.582476],[-4.832038,9.176151,10.41875],[6.330236,-2.287717,6.804818],[-2.143137,-1.943817,3.061284],[7.444915,0.04331499,7.5119],[-4.550204,2.669474,5.3694],[6.337031,3.287456,7.208698],[-6.49394,1.94681,6.852834],[-3.683126,-4.98513,6.278291],[3.347065,0.6817288,3.559157],[8.116897,7.779674,11.28749],[-6.062768,9.108914,10.9877],[5.748271,6.935249,9.063128],[-5.349812,1.882278,5.758772],[-2.015656,-8.373559,8.670604],[2.711586,0.1076271,2.892106],[-6.189147,-3.397052,7.130603],[6.975004,-5.159878,8.733557],[-4.958733,8.991615,10.31689],[7.526658,6.973432,10.30919],[8.567953,-0.9580713,8.679154],[9.975401,1.598837,10.15209],[6.041345,-9.076909,10.94934],[0.720969,1.088291,1.644437],[-9.08722,-1.188121,9.218958],[-5.056805,-6.935944,8.641678],[-9.950916,1.643418,10.13517],[-1.629544,2.149483,2.87675],[-1.049801,-1.286897,1.938604],[-3.32688,2.700846,4.400307],[-5.421143,6.374725,8.427687],[4.811043,2.814184,5.662665],[-8.4691,-5.725873,10.27187],[-6.151237,1.877115,6.508554],[5.530796,-4.031139,6.916631],[7.170907,0.4671847,7.255354],[7.624505,-3.574853,8.480132],[6.535028,-6.724841,9.430275],[-8.015165,0.3090455,8.083217],[-9.763792,-8.852089,13.21708],[-4.966646,9.179585,10.48486],[8.264064,-9.493368,12.62612],[0.792177,4.757202,4.925293],[-7.544542,0.2121989,7.613484],[7.418346,9.358835,11.98414],[0.6238882,-3.55512,3.745412],[4.611703,-0.7680504,4.780974],[2.755252,8.299909,8.802267],[2.670033,4.995256,5.751667],[4.142551,-3.068704,5.251445],[-0.6547215,4.372291,4.532724],[5.292856,2.901232,6.118127],[-2.095939,3.986862,4.613895],[4.869792,7.657408,9.129664],[-6.415383,7.651375,10.03497],[-4.856083,6.020186,7.798986],[-9.141344,4.07512,10.05837],[0.01039044,-4.947556,5.047616],[-1.786279,8.403962,8.649704],[8.07268,8.399835,11.69296],[6.13905,3.425942,7.101058],[5.234118,-4.162678,6.761944],[4.071911,1.868322,4.590325],[3.754162,-2.850402,4.81856],[-1.249953,-4.171374,4.467969],[-9.665374,4.591088,10.74698],[-1.538786,-9.451509,9.628027],[3.459909,-6.121176,7.102096],[-8.718925,8.848378,12.46248],[-1.721624,-1.390285,2.428349],[-9.694456,-7.900251,12.54577],[-1.322114,-4.182253,4.498803],[1.961344,-7.368572,7.690431],[7.04318,-4.892656,8.633914],[-3.476261,1.772403,4.028127],[5.383033,4.757623,7.253415],[6.691832,3.184972,7.478279],[7.325391,-6.526647,9.861971],[1.287592,3.93273,4.25726],[-3.995226,-8.459115,9.408424],[-1.07775,-1.074826,1.821207],[-5.844141,2.550693,6.454457],[-8.540077,-2.51236,8.957951],[8.675397,-7.796515,11.70676],[7.784962,8.286751,11.41385],[-0.7455349,-0.6947383,1.427755],[-5.850648,6.104971,8.514738],[-3.595155,-2.165032,4.314221],[-3.448289,4.770738,5.970816],[-1.969774,0.1866018,2.216941],[5.405728,-2.552995,6.061327],[6.135682,2.741385,6.794246],[-5.114353,6.120713,8.03864],[-3.072687,-3.195878,4.544782],[0.9492246,0.9237468,1.659619],[-9.078024,1.475537,9.251364],[-3.824661,1.910584,4.390714],[-9.181222,5.48505,10.74154],[-7.568254,2.270387,7.964491],[1.718228,4.593277,5.005047],[-6.650422,7.433678,10.02435],[-5.867073,6.514933,8.824222],[6.994081,7.994789,10.66929],[-8.978614,-7.681035,11.85807],[4.02502,-5.457002,6.854171],[-1.260802,0.7790422,1.787884],[-7.328275,-2.503732,7.808475],[0.7039286,-6.548667,6.661873],[0.2769133,5.23121,5.333127],[8.422505,-8.56652,12.05503],[-7.559659,-9.154775,11.91463],[2.533987,2.124338,3.454547],[1.38788,7.290708,7.4887],[7.434774,-8.968016,11.69193],[5.004758,1.641106,5.361048],[-6.868747,1.726777,7.152723],[4.383638,-3.968881,5.997358],[7.916744,-9.508169,12.4129],[-7.566618,-5.755548,9.559291],[-4.09993,5.438019,6.88342],[9.56525,5.999496,11.33525],[-7.727228,9.550424,12.32561],[8.463212,-4.984123,9.87256],[-5.784758,-4.220578,7.230263],[-1.402976,5.879892,6.12711],[1.783133,-2.683408,3.373461],[3.910944,-9.39973,10.22988],[3.278206,5.724522,6.67209],[9.024782,-9.307799,13.00315],[-5.218582,-6.007776,8.02041],[-9.922637,0.9718275,10.02014],[1.338774,1.172606,2.041401],[9.999704,-0.6919439,10.07337],[-7.05016,1.017414,7.193045],[4.560174,-6.900802,8.331641],[-1.84303,-2.004161,2.90059],[8.999015,5.252398,10.46757],[-3.200114,6.322969,7.156862],[-6.804102,-1.431819,7.024664],[-6.464664,-5.204144,8.359127],[7.274015,-0.9632553,7.405347],[-2.695135,-8.431828,8.908395],[-1.406049,-5.427601,5.695246],[9.552304,-4.591383,10.64553],[4.168467,-2.778877,5.108647],[7.764115,-6.202511,9.987623],[6.151637,-4.665113,7.78498],[-7.912083,0.8083904,8.015894],[3.668775,7.283665,8.216549],[-8.057689,-2.170293,8.404553],[-0.6293378,0.4412126,1.261243],[-7.865199,-3.892821,8.832633],[5.529031,-7.950016,9.735139],[7.529834,-6.324677,9.884328],[1.222427,-8.840517,8.980483],[-2.99599,-1.362308,3.439744],[2.805335,-4.665936,5.535419],[5.860879,-6.081991,8.505323],[9.51487,0.1619051,9.568645],[3.147053,3.428732,4.760268],[7.807281,7.948546,11.18629],[6.140275,-7.819922,9.992705],[3.439135,4.902809,6.071671],[-2.432884,-2.605927,3.702672],[6.162295,-6.15423,8.766323],[4.279854,-4.92226,6.598924],[1.969283,-2.876886,3.62692],[0.4123502,5.364861,5.472821],[-8.826879,-4.462129,9.941046],[4.154095,-5.021147,6.593059],[3.627757,8.607985,9.394575],[0.7112548,3.729844,3.926528],[-1.993246,-2.53299,3.374769],[-5.544368,6.877139,8.890166],[3.851819,-0.129206,3.981608],[-0.8989518,-4.910514,5.091293],[-8.890851,-7.50208,11.67598],[4.23878,-6.034218,7.44171],[3.805643,-2.373648,4.595336],[2.823805,-1.868224,3.530458],[-2.574398,8.30294,8.75022],[-9.212195,9.744407,13.44686],[-9.554458,8.467089,12.80544],[-5.332816,0.0009258324,5.425765],[-0.9352102,9.504898,9.603005],[8.596414,-2.744575,9.079153],[-3.264433,-5.331437,6.330936],[-3.293745,-7.557235,8.304249],[-8.672374,4.857015,9.990028],[-1.314236,0.8300928,1.848316],[-7.240139,-4.592278,8.631838],[-9.172147,4.696818,10.35318],[8.774025,9.635779,13.07026],[2.92818,-2.39205,3.911027],[0.8265426,-2.424136,2.749474],[7.053493,-9.82773,12.13821],[2.448555,-1.017013,2.833679],[-5.855382,0.8942896,6.0071],[-5.791171,0.3031885,5.884691],[7.696047,-2.334458,8.104248],[-5.895906,1.339982,6.128397],[5.642587,2.923238,6.433048],[-1.476298,-9.705963,9.868393],[1.199568,-9.519218,9.646475],[-6.060663,-6.947519,9.2736],[7.246321,3.176785,7.975032],[2.540748,-9.233099,9.62837],[-9.383009,5.4703,10.90711],[3.183561,-9.36678,9.943421],[8.713938,0.4855293,8.784557],[-3.936769,-2.939817,5.014048],[2.288825,-4.950212,5.544666],[6.049231,5.801808,8.441218],[-3.310434,3.704368,5.067674],[-8.692095,9.634013,13.0141],[8.133362,9.612084,12.63106],[0.2367634,6.693104,6.771535],[-3.289825,0.8122867,3.533095],[2.159907,7.402987,7.776208],[3.632519,-6.466856,7.484345],[-7.370052,-3.641711,8.281288],[8.692006,1.778056,8.928183],[3.560568,-7.012141,7.927659],[-0.08190703,-8.172306,8.233668],[-3.67495,0.5479887,3.847798],[-0.9451037,4.791742,4.98538],[5.55002,4.012274,6.92106],[1.400135,-6.18532,6.420168],[-1.488609,-8.145972,8.341032],[-2.17092,-3.973686,4.637141],[-1.560124,5.96696,6.248088],[-4.682688,8.050237,9.366637],[-5.294898,4.340503,6.919242],[-2.170773,-6.510557,6.935388],[4.201369,1.870418,4.706375],[-5.858255,0.4665318,5.961276],[-4.277745,-8.112552,9.225649],[-5.506805,4.040856,6.903145],[-2.374068,1.03564,2.776463],[3.265788,-5.774648,6.709093],[-8.118911,-5.426979,9.816762],[-0.7242826,7.706254,7.804546],[7.914749,3.395502,8.670218],[8.893354,2.763813,9.366451],[7.144872,-3.62748,8.075134],[5.121372,-7.210237,8.900334],[-8.18878,-7.486156,11.13996],[-7.578259,-9.826997,12.44989],[-6.090753,-3.248941,6.975162],[-3.956844,5.246691,6.647133],[8.056703,8.81746,11.98574],[6.687154,-8.964364,11.22844],[5.290472,1.380838,5.558399],[3.720553,-7.417311,8.358171],[-1.965638,9.277436,9.535961],[-4.075217,-9.418243,10.31071],[-7.758079,0.8201176,7.865137],[-5.283161,-8.855816,10.36037],[6.283272,6.495422,9.092305],[-8.565834,-6.997786,11.10597],[9.785643,8.870648,13.24565],[3.527998,-9.785699,10.4502],[5.018051,-9.675211,10.94489],[9.153428,-7.888123,12.12467],[-4.045822,5.620289,6.99688],[-8.835875,0.1141629,8.893014],[-6.659587,-1.741175,6.955702],[-5.539852,-0.6598487,5.667924],[-1.904704,0.4626071,2.200433],[3.239472,-3.512659,4.8819],[-7.717862,-3.63658,8.590117],[6.68314,0.725951,6.796423],[1.892064,8.29648,8.56805],[-9.954054,-2.807493,10.39063],[2.279047,1.640579,2.980865],[-4.346849,-9.43797,10.43889],[3.346283,8.527424,9.21491],[1.97928,-8.949726,9.220366],[-5.320354,-6.095983,8.15274],[-3.523422,-9.8565,10.51499],[-3.336872,-5.108668,6.1833],[-4.722311,-5.293189,7.163662],[0.0736405,-1.856202,2.109717],[2.243383,5.463639,5.990335],[-0.711205,0.7465994,1.436392],[4.693914,-0.9250808,4.887597],[-6.137969,-0.1914783,6.221842],[-8.098494,-5.737525,9.975209],[-3.131614,3.321217,4.67306],[9.953545,7.441852,12.46813],[6.453738,4.913034,8.172432],[8.816946,9.869316,13.27185],[-5.342802,5.281577,7.578957],[8.221967,-8.2325,11.67796],[4.191467,-2.613358,5.039647],[-3.5353,-9.21453,9.919975],[1.803554,8.237839,8.492043],[-9.466272,-8.853731,12.99996],[-1.522986,3.313057,3.780983],[-8.295181,5.374096,9.934331],[-0.5882798,5.95962,6.071502],[-3.389609,-9.942723,10.55212],[-2.939474,-8.926507,9.451087],[-1.414424,-3.136452,3.583005],[-7.412049,4.524833,8.741428],[-2.014258,-0.3502819,2.275947],[6.236513,6.371331,8.971508],[-1.581704,6.24323,6.517646],[0.5131085,0.612433,1.279982],[-0.7568983,-6.707337,6.823582],[-1.393305,8.373328,8.547159],[3.101409,-8.435834,9.043343],[-3.08002,1.424913,3.537923],[8.464237,-6.593049,10.77551],[0.9276742,1.240577,1.843803],[-4.879653,7.17759,8.736636],[6.035867,-0.02667818,6.118203],[7.412883,-2.26095,7.814265],[-4.928728,7.198143,8.78098],[-8.565647,7.523532,11.44438],[6.303196,7.135712,9.573331],[-4.296556,6.112131,7.537807],[8.503992,-0.9757835,8.618007],[7.671831,-0.779765,7.775926],[8.785624,3.955292,9.686666],[2.497396,-5.616624,6.227636],[0.7682014,-1.883587,2.266723],[-3.402829,9.304497,9.957556],[2.87362,-4.376971,5.330626],[3.828948,-7.480005,8.462347],[5.714533,5.295773,7.855004],[-6.797708,9.181025,11.46735],[8.640022,-1.152795,8.773764],[0.03164914,-8.938588,8.994407],[8.499027,8.124404,11.79998],[3.160216,1.088999,3.488966],[3.928357,-8.372392,9.302094],[7.342614,5.456699,9.202692],[9.45303,2.720441,9.887395],[-3.230067,-2.176985,4.021517],[-1.342089,-8.329427,8.495914],[1.230018,4.268458,4.553316],[-8.531075,-7.460285,11.37696],[-8.198013,9.449995,12.55029],[5.159562,-0.2912545,5.26364],[-8.352279,-9.535994,12.71596],[-8.431263,6.84596,10.90657],[9.131682,-0.1572132,9.187619],[-8.313376,1.304107,8.474251],[0.4557128,-2.061645,2.336248],[1.700167,-3.839171,4.316226],[3.884583,-5.498198,6.805892],[5.773306,3.497218,6.823606],[7.41402,-4.240247,8.599267],[-2.942019,2.405438,3.929582],[-3.302036,-4.220436,5.451195],[-4.294317,-0.3947437,4.426848],[-1.241212,7.682295,7.845907],[6.907396,-4.056017,8.072386],[0.9935642,0.5196707,1.502407],[-3.637551,-3.440467,5.105741],[6.775263,-3.066444,7.503817],[1.265105,2.632757,3.087377],[-8.00841,1.299985,8.174631],[-0.2172656,7.461051,7.530901],[-9.778321,-4.549777,10.83125],[-0.6144368,6.911464,7.010411],[7.844446,-1.856752,8.122984],[-0.1265412,1.148571,1.528145],[1.0589,-9.837396,9.944628],[6.392903,-3.821232,7.514721],[3.944023,9.161726,10.0246],[-7.657671,-2.919956,8.256274],[-9.008518,-0.7284731,9.093079],[-3.868994,8.203053,9.124648],[1.618934,-6.773863,7.036062],[-0.9805614,2.878367,3.201015],[7.956067,4.648369,9.268567],[5.056073,-8.352242,9.81447],[5.104334,-6.208344,8.099244],[3.617593,-0.6865909,3.815545],[9.172071,7.334178,11.78631],[-9.912702,0.6206676,9.982328],[0.7140785,8.723207,8.809328],[-4.645825,-7.40603,8.799601],[7.955473,-3.644958,8.807683],[-8.705577,-0.5884283,8.782557],[-6.85583,1.86845,7.175898],[7.758985,1.574982,7.980126],[9.500363,5.356575,10.95216],[-7.086548,9.562607,11.94415],[6.524494,9.888138,11.88883],[-7.78012,8.737201,11.74176],[7.2475,-5.788234,9.328982],[3.910177,-7.044117,8.118441],[3.375908,9.453892,10.08825],[-8.426269,5.493369,10.10837],[-4.907058,-3.54528,6.135815],[-6.25844,7.458004,9.787231],[2.795188,1.345079,3.259189],[-3.092826,-2.883555,4.345165],[-1.881412,-2.403254,3.211751],[3.913477,-9.071575,9.930196],[-8.622274,0.2618967,8.68402],[-0.7237727,-4.273581,4.448296],[-2.745336,6.618085,7.234357],[1.608967,-6.53871,6.807607],[0.8165967,1.655607,2.099491],[2.062727,-6.638676,7.023309],[4.187735,-1.109694,4.446184],[-5.060585,3.411557,6.184517],[-4.555074,1.884863,5.03005],[-5.793154,1.701741,6.120176],[-5.690922,-6.000469,8.330199],[-1.156797,-4.580519,4.82901],[-0.7730715,-7.014665,7.127634],[-0.5024052,-8.783413,8.854421],[-9.239203,8.027947,12.28051],[1.125546,-0.3600941,1.54807],[-0.3860157,4.944257,5.059119],[1.757559,1.945758,2.806241],[1.284689,-4.376941,4.669908],[-8.001415,-2.577013,8.465439],[1.573816,1.489577,2.386575],[3.392912,-3.719109,5.132604],[-3.465774,0.5096141,3.642979],[1.421127,-2.903837,3.384061],[8.839343,6.99048,11.31374],[-4.118101,-1.985243,4.679737],[-8.086476,0.6165342,8.171365],[-1.621389,6.415367,6.692222],[6.521732,-8.601357,10.8405],[-2.816975,-8.856125,9.346994],[-7.756667,4.020779,8.793893],[-0.1591686,-8.545244,8.60503],[-4.190453,0.8239931,4.386213],[5.078227,-6.090568,7.992709],[-6.396691,-2.76756,7.041097],[9.088747,5.826882,10.84241],[-3.593327,0.7480958,3.804162],[-1.961426,4.122888,4.673906],[3.413041,-9.669045,10.30239],[-4.711914,-3.124594,5.741535],[-5.072612,5.135477,7.287284],[-3.018935,-5.932097,6.730806],[-2.230787,-3.689917,4.426274],[5.274719,7.409838,9.15032],[-7.073236,1.381904,7.27601],[2.044195,-8.662006,8.955953],[-7.896027,-2.239056,8.268047],[9.178702,1.314698,9.326146],[-3.869946,-6.842029,7.924005],[-5.868618,-8.027988,9.994462],[1.466221,0.5784342,1.866652],[-3.897302,-0.04512352,4.023804],[1.378318,-3.488064,3.881539],[1.096284,-8.56199,8.689621],[9.181385,-0.08216676,9.236048],[9.646998,9.204898,13.37142],[-3.861032,-1.322104,4.201848],[7.402659,-0.3434136,7.477787],[-5.103405,-2.59179,5.810518],[5.346535,4.985034,7.378075],[7.395342,-6.546007,9.926796],[-2.752658,9.401275,9.846883],[0.9455354,6.609884,6.751637],[8.459929,-9.250417,12.5754],[-5.667725,1.091076,5.857778],[-7.681514,5.589438,9.552354],[3.97705,-4.89619,6.386674],[0.6886155,7.557859,7.654765],[-8.172928,-7.013574,10.81605],[-4.647655,-8.440685,9.687407],[1.972289,0.5841413,2.28717],[-4.61665,7.949775,9.247291],[8.160316,-3.868616,9.086085],[5.690356,-4.399115,7.261706],[3.956909,8.911331,9.801477],[-0.1949912,4.317881,4.436453],[5.862401,7.821668,9.825794],[3.105625,-1.675085,3.667536],[-4.032102,-9.683037,10.53656],[-1.313643,9.163334,9.310873],[-3.825372,-7.171228,8.189016],[9.458878,3.039787,9.985523],[7.832903,-5.233397,9.473269],[6.392731,-1.636837,6.674297],[7.161544,-3.373244,7.979128],[4.956673,6.113239,7.933493],[8.604821,7.649539,11.55675],[-2.944657,-9.340287,9.844388],[3.42394,6.831837,7.706968],[8.656447,-7.75092,11.66237],[9.708483,3.337308,10.31466],[-4.535065,-3.024656,5.542144],[9.336043,0.3869878,9.397418],[2.933005,3.685217,4.814909],[2.71901,-1.106613,3.101227],[2.132254,4.725662,5.279999],[-6.66116,-5.696234,8.821459],[7.975079,5.153825,9.547974],[7.036418,-8.186328,10.84099],[-1.742679,8.4203,8.656695],[-0.4990797,-2.765613,2.9829],[7.775609,2.484919,8.224045],[-4.103766,5.366756,6.829566],[2.275336,-3.08777,3.963771],[5.349056,-6.786456,8.698758],[1.755314,7.147473,7.427482],[4.533585,-6.465954,7.960022],[4.773693,2.207315,5.353539],[-1.824898,-2.552813,3.293494],[-9.11576,-5.125728,10.50572],[4.807636,-8.137838,9.504619],[5.00536,-1.630387,5.358338],[-1.539052,-3.082522,3.587565],[5.812063,3.763015,6.995738],[-4.002409,-1.400848,4.356793],[3.371857,-4.558335,5.757416],[5.732923,0.06749072,5.819876],[4.851419,-6.608312,8.258696],[3.692251,-0.7638254,3.900788],[9.116572,1.798292,9.345895],[4.769442,8.844292,10.09797],[4.547272,-9.749542,10.80422],[9.863743,8.947307,13.35469],[4.352216,-8.853247,9.915733],[-9.889881,-0.5347189,9.954681],[9.309204,-9.932529,13.64978],[5.673664,-1.703194,6.007607],[-0.5956883,-4.122004,4.283195],[3.836388,6.140979,7.309548],[8.18159,6.457518,10.47081],[8.332293,8.397147,11.87178],[-4.209202,0.5073656,4.356008],[-7.775099,-6.186798,9.986423],[-0.9909362,8.003054,8.125936],[-8.425021,-1.822727,8.677749],[-9.978869,8.070058,12.87259],[-7.481596,0.1370275,7.549375],[-2.495724,6.290946,6.841392],[-1.090874,8.683802,8.808996],[8.785177,-9.453193,12.94381],[-2.183738,2.365319,3.370971],[9.343595,-4.953338,10.62254],[-8.337256,9.915389,12.99326],[6.92354,-7.58918,10.32139],[-2.601277,6.513917,7.085038],[-2.105528,-9.924904,10.19495],[6.467261,0.3042453,6.551185],[8.813169,-4.258473,9.839032],[8.997934,3.259588,9.622252],[3.9247,4.507986,6.060133],[-5.553974,1.132167,5.75573],[9.076445,-2.747355,9.535712],[-3.265733,-8.678535,9.326413],[8.652298,5.265437,10.17777],[-2.305889,9.93751,10.25043],[5.63966,9.888097,11.42717],[-3.467799,-7.871281,8.659255],[-2.520455,1.654067,3.17626],[2.402801,-0.9343385,2.76522],[-5.986443,2.140195,6.435677],[5.106097,-2.257247,5.67163],[4.336025,-1.898591,4.83795],[-9.067716,-5.574539,10.69107],[1.102418,-7.011524,7.16776],[-0.1045234,8.041994,8.104604],[-2.607809,-7.580831,8.078965],[1.627296,0.8914442,2.107787],[9.53703,2.528933,9.917179],[-7.954943,-0.9258319,8.070829],[0.9461079,5.620181,5.786324],[7.831127,-3.684137,8.712027],[-8.347431,-3.432151,9.080709],[-9.555387,-5.241915,10.94455],[4.215145,4.907137,6.545796],[5.462389,8.148331,9.860679],[9.101638,-7.973267,12.14137],[-2.201456,-3.678858,4.402318],[-3.395812,-8.970531,9.643753],[3.906039,5.375623,6.719707],[3.611562,2.168442,4.32961],[1.232855,-7.173386,7.346931],[3.978066,7.842439,8.850359],[2.780921,-2.756804,4.041471],[4.105873,-7.507646,8.615273],[-0.5040063,-3.702314,3.867964],[4.220964,3.276284,5.436044],[-1.579259,-5.149135,5.477924],[-0.5051661,-0.6325328,1.286581],[3.260526,7.675269,8.398856],[6.053866,9.304554,11.14558],[4.734195,2.156251,5.29736],[3.267867,3.622482,4.980093],[-0.6172408,0.9578243,1.516052],[-4.697692,7.119747,8.588312],[-2.688927,-8.328742,8.808988],[-3.196027,-4.617124,5.70372],[2.427577,-5.541341,6.13185],[-2.150837,9.529118,9.819887],[-1.320893,-8.228212,8.393345],[2.067204,5.898046,6.329319],[5.339082,5.582166,7.788863],[4.512345,-9.911207,10.93587],[7.45048,-1.079528,7.594408],[1.417976,-3.187948,3.629555],[-6.291715,2.334022,6.784787],[-6.080217,-5.88175,8.518452],[-9.147601,-4.219443,10.12335],[3.750812,-1.620354,4.206439],[0.3241295,-3.177555,3.346926],[-6.179671,4.504398,7.712194],[-6.08424,7.664049,9.836444],[8.316434,-9.808362,12.89833],[1.923817,-9.146397,9.399875],[-4.468796,1.134504,4.717758],[2.356532,8.779983,9.145563],[0.8285567,-6.005584,6.144391],[-6.992067,7.235173,10.11122],[5.15147,1.684553,5.511384],[9.296889,-4.551106,10.39927],[-1.847455,-6.277057,6.619255],[8.059102,0.06481592,8.121165],[-3.889967,0.2191482,4.022421],[5.410938,-3.943845,6.769946],[2.6709,-8.405908,8.876541],[9.638663,8.062624,12.60594],[7.054151,2.757605,7.639727],[-0.09691031,-6.591172,6.667304],[-9.656896,2.050215,9.922651],[-4.325148,-7.947258,9.103066],[-0.9131834,-6.018161,6.168644],[-6.497045,-0.3163068,6.581158],[5.372506,2.394362,5.966305],[5.930819,8.248613,10.20854],[5.689508,-4.188765,7.135562],[-2.567661,7.567753,8.053804],[5.127347,-4.922411,7.177731],[8.42156,-2.324631,8.793553],[2.029898,2.409332,3.305354],[-8.415689,-4.036847,9.387223],[-5.967998,-7.034938,9.279405],[8.781093,-2.130947,9.091124],[4.463771,6.290485,7.777883],[0.08966027,9.374501,9.428113],[1.150754,-0.3150168,1.55675],[7.970949,-7.072564,10.70314],[5.182117,5.436411,7.576867],[-6.135428,7.959688,10.09951],[2.610349,-9.97078,10.35521],[7.076854,1.20478,7.24799],[-6.481087,-1.162195,6.659968],[1.895936,7.711941,8.004287],[-6.667254,4.935257,8.35518],[-5.296018,2.445245,5.918364],[-0.5569239,-8.435905,8.513205],[0.1368878,8.119997,8.182487],[-3.482261,-3.570922,5.087005],[3.966849,-3.264455,5.23379],[-2.449784,2.049441,3.346886],[4.240481,-7.701833,8.848724],[0.6043187,7.837348,7.923965],[-8.278341,8.984646,12.25785],[-6.161838,4.471349,7.67862],[-5.487699,5.878633,8.103898],[2.370713,8.315944,8.704896],[-3.844634,8.736897,9.597634],[6.007084,5.301047,8.073794],[7.791266,-8.603666,11.65019],[6.70668,-6.936698,9.700378],[1.358993,-4.007942,4.348616],[5.277054,8.271118,9.861982],[-4.599102,-5.383015,7.150426],[7.969011,-4.946503,9.432551],[-7.164047,-3.292398,7.947544],[-1.916175,-6.108802,6.479907],[9.444728,-3.534904,10.13402],[-6.236247,-0.7823165,6.36418],[-1.491895,-7.745686,7.951188],[2.211836,-3.187731,4.006725],[1.084832,9.602946,9.715629],[-7.796807,6.173684,9.995227],[-7.631226,-7.459112,10.71793],[9.66765,3.641966,10.37918],[-7.335021,3.068607,8.013669],[-1.090342,1.170644,1.886599],[-7.582307,6.215689,9.855261],[3.180058,-0.6471342,3.395813],[2.933642,4.250311,5.260361],[-2.732233,-1.409439,3.232896],[2.301126,3.853234,4.598107],[-3.981205,6.310432,7.528051],[0.8219969,0.3833369,1.350047],[9.932079,-1.64114,10.1163],[8.292449,-6.574603,10.62968],[8.287403,4.382566,9.42804],[-9.848344,0.9773747,9.947118],[9.762671,-3.059087,10.27948],[9.74319,3.064892,10.26271],[7.743554,-9.675999,12.43332],[-3.020437,3.808021,4.962264],[-9.522521,3.212657,10.09948],[-4.564081,-7.022442,8.434781],[-1.069704,-1.786735,2.310127],[5.018824,-2.62077,5.749524],[0.1771266,1.236532,1.60012],[-8.939567,0.03542064,8.995394],[8.677719,-1.215152,8.819263],[7.328279,-7.001496,10.18453],[-4.281662,-7.557171,8.743196],[-7.810601,-1.636136,8.042539],[-6.456265,5.897349,8.801255],[-7.107834,-9.0732,11.56911],[4.959996,-1.091613,5.176212],[-5.300412,6.672563,8.580062],[-9.575088,8.472,12.82408],[5.361254,8.754495,10.31427],[-4.373386,-2.996463,5.394932],[-3.531559,9.871614,10.53189],[5.410559,-7.411525,9.230647],[-4.617996,2.156578,5.193912],[-8.544888,-1.239986,8.692105],[-1.830522,8.194637,8.455938],[-9.675075,-1.651282,9.865789],[3.906553,1.44596,4.283918],[-9.606283,-2.653954,10.01619],[3.768788,9.509815,10.27815],[-4.01157,6.817894,7.973479],[-8.408295,-8.283887,11.84577],[-9.262033,5.830087,10.98977],[2.604995,-7.574852,8.072447],[-5.482109,-5.821367,8.05865],[5.420883,-8.258027,9.928796],[7.201836,5.135153,8.901474],[3.350782,-6.56102,7.434697],[-0.1416261,7.293512,7.363109],[-7.454583,-5.13459,9.106855],[-0.3547255,8.396651,8.463427],[-5.425512,5.441562,7.748986],[9.317666,9.670343,13.46605],[7.863175,-2.718521,8.37973],[2.88214,2.472476,3.926814],[3.309585,-2.386222,4.200882],[8.763013,8.741624,12.41799],[-7.972623,-1.486765,8.171486],[-2.923015,6.109588,6.846246],[-0.03019246,-3.511232,3.650981],[1.793492,-7.185418,7.473074],[-5.713451,2.561077,6.340555],[-5.431937,-2.668948,6.134266],[9.065042,-3.073666,9.624055],[6.120029,9.066875,10.98467],[-9.213751,-3.279983,9.831149],[-4.052538,5.73931,7.096672],[-1.886525,2.995689,3.67874],[0.9600481,-6.303245,6.453882],[-3.295546,-5.545703,6.52805],[8.648563,0.4917242,8.720059],[9.667921,-4.706967,10.79927],[-6.913043,-0.8779425,7.039954],[3.22879,8.885417,9.506614],[-3.586608,-9.519876,10.22212],[1.552099,3.108614,3.61559],[-4.85323,-1.697058,5.237733],[-0.3494792,-1.275768,1.658228],[-9.065651,7.923227,12.08154],[6.161022,-2.578762,6.753385],[1.517751,-1.669029,2.467636],[5.015759,8.398524,9.833262],[9.134568,1.402121,9.295498],[-5.784668,-1.840547,6.152236],[-6.487428,4.764674,8.111032],[9.721638,-1.810348,9.939196],[7.995433,2.74288,8.511777],[6.96268,-0.8029315,7.079803],[-7.067807,-4.134576,8.249158],[0.08312746,0.3273996,1.05551],[6.319187,6.652804,9.229948],[2.744331,-9.593898,10.02867],[6.340886,-5.302591,8.326122],[-6.533237,-9.609505,11.66301],[-7.847186,3.975194,8.853276],[9.276086,-9.42109,13.25906],[-2.567999,-2.066199,3.444386],[1.507167,5.727664,6.006471],[-4.377569,-2.875891,5.332341],[7.277967,-1.319549,7.463914],[-2.462607,-7.705773,8.15128],[-7.871677,-5.854197,9.860778],[3.859045,-6.647035,7.750825],[1.396679,7.714634,7.903562],[-2.547669,-2.898422,3.986411],[-5.288947,5.996678,8.058108],[-6.080444,7.658011,9.829391],[0.05636749,7.471814,7.538646],[-8.501192,9.417291,12.72618],[-0.7431832,0.2273937,1.266503],[-3.544834,1.713762,4.062367],[-6.611368,-9.969915,12.00456],[3.658901,-6.477952,7.506758],[-1.937926,3.586168,4.197161],[7.396028,-4.337619,8.632275],[-6.790079,-8.161118,10.66344],[4.299419,-9.880243,10.82147],[4.858232,-8.468619,9.814271],[-0.6557209,4.207897,4.374513],[6.667487,8.118443,10.55294],[0.5836998,1.159444,1.638601],[-6.488178,-6.974902,9.578398],[2.494824,-2.234003,3.494985],[3.418528,-5.970246,6.951991],[-5.434537,5.846418,8.044551],[-2.984448,1.861628,3.656855],[8.471598,0.356605,8.537865],[1.022493,-9.687049,9.792058],[3.30419,-9.289217,9.909955],[4.736517,8.961623,10.18554],[4.489333,5.908398,7.487541],[7.227281,-6.272296,9.621605],[0.1291775,-6.084775,6.167753],[-3.429725,-1.406945,3.839597],[5.127396,3.732782,6.42058],[-1.607248,-3.821559,4.264688],[1.977532,5.87404,6.278134],[-3.763825,7.213028,8.197205],[8.385225,9.06651,12.39006],[-3.180408,8.625345,9.247247],[-6.858607,-2.488442,7.364295],[0.0793945,2.436035,2.634496],[9.708916,-2.579351,10.09535],[-0.3738451,3.868665,4.013269],[1.51049,-7.787713,7.995627],[-2.040331,1.895448,2.958999],[-9.041863,6.453693,11.15372],[3.129939,-4.861234,5.867547],[0.9626269,4.778871,4.976369],[-2.750983,-2.521036,3.863099],[-2.193401,8.803934,9.127994],[5.204125,7.932351,9.53966],[1.154904,8.249766,8.39002],[4.943279,-4.522977,6.774462],[2.000416,7.748052,8.064364],[0.7200015,-9.387012,9.467544],[-3.897451,-4.355172,5.929388],[3.905285,5.09567,6.497469],[-8.482476,0.3613054,8.548857],[8.510964,6.877924,10.98828],[-0.743115,0.9819769,1.586348],[8.531843,-5.302313,10.09489],[6.217713,4.946216,8.007809],[1.633975,4.232572,4.645917],[-4.717134,-6.524675,8.113122],[-6.248425,-6.603447,9.145947],[4.59359,6.371948,7.918509],[-1.281706,-3.075392,3.478622],[5.410211,4.201266,6.9225],[-0.1705726,7.193582,7.264758],[-1.385369,9.24435,9.400918],[-4.40051,-1.268039,4.687474],[-9.958926,-4.981251,11.18003],[8.818188,-1.761084,9.047753],[-7.385609,5.140891,9.054059],[-3.32,-1.631031,3.831796],[-7.729207,-4.905153,9.208755],[-2.461737,2.965813,3.981984],[9.031911,-3.402581,9.703245],[-8.864042,-3.734098,9.670301],[-9.782412,5.585088,11.30879],[8.7902,8.92574,12.56728],[-9.818149,5.329083,11.21585],[8.925853,3.06943,9.491693],[-8.117188,8.771082,11.99252],[-9.397514,-4.039348,10.27763],[6.94065,-9.412771,11.73767],[-5.480756,-6.987105,8.936349],[3.096016,8.265113,8.882421],[8.752796,4.648306,9.960833],[-3.05083,-1.48056,3.535481],[-8.737484,-6.986384,11.23179],[-2.19229,2.750065,3.656363],[-1.080054,-6.289668,6.459601],[-2.671841,-2.975795,4.122388],[4.840734,6.740364,8.358541],[-2.715479,2.035814,3.53813],[9.194938,-0.0387273,9.249236],[-9.348006,-6.445932,11.39892],[2.019125,7.586775,7.914292],[2.274336,-1.356365,2.830606],[-1.594979,8.413983,8.622011],[-7.005465,-4.985167,8.656121],[0.9863861,-5.400437,5.580114],[-0.2337728,0.7345783,1.262638],[-8.779316,2.177715,9.100486],[1.630014,-6.452364,6.729781],[-3.666707,-4.291912,5.732822],[7.657179,-9.014653,11.86998],[9.951585,-8.479381,13.11236],[-1.216702,6.454689,6.644048],[-3.417977,1.30143,3.791607],[1.323982,3.030725,3.455173],[-3.007133,-1.294123,3.423099],[2.046302,6.532696,6.918343],[-2.579396,3.102377,4.156684],[9.671131,-5.627823,11.23402],[-0.4459549,-3.59103,3.754247],[-5.913123,-0.6204882,6.029099],[5.300447,-6.734713,8.628505],[6.883778,-6.678185,9.64285],[-6.016191,-1.448033,6.268282],[-1.129391,4.797773,5.02933],[-1.562144,-2.601199,3.194766],[-8.516914,5.644536,10.26638],[8.637904,-4.944825,10.00323],[1.84358,-5.70316,6.07658],[-0.3087837,4.514617,4.63434],[0.9342974,-9.003939,9.107349],[-2.630512,0.168294,2.819205],[-0.7690254,8.927218,9.015909],[-1.017821,9.267973,9.377168],[1.557907,9.513385,9.69183],[-4.86746,5.224019,7.209891],[4.000099,0.09461228,4.124287],[-8.753507,-1.648803,8.963394],[-6.55331,3.447573,7.472057],[-0.429565,1.317621,1.708991],[-8.199811,3.487483,8.966573],[-1.661065,-0.3228753,1.96555],[-8.716265,3.264412,9.361071],[2.647088,1.676738,3.289152],[9.638319,7.559939,12.29024],[-0.8764858,-7.066384,7.190411],[9.429744,0.7529694,9.512466],[9.847383,0.9298601,9.941609],[6.25801,-3.090954,7.051006],[2.462018,0.5678892,2.717358],[-3.838681,9.902414,10.66739],[-5.681934,8.464965,10.24402],[7.963779,-2.552094,8.42229],[-7.913896,-2.666478,8.410699],[6.1714,1.13241,6.353623],[2.744889,-2.826753,4.065089],[-7.316671,-9.798336,12.26952],[-8.941626,9.025511,12.74412],[8.116389,-6.092865,10.19798],[1.540646,-1.394418,2.306077],[8.462674,-1.353776,8.628416],[-5.331528,-4.361621,6.960526],[-5.127351,-7.991109,9.547123],[-8.982678,-7.297693,11.61658],[5.724152,0.6507374,5.847167],[-6.261768,-4.235706,7.625677],[-4.329756,9.172103,10.19187],[-6.218832,0.7122765,6.338865],[2.125359,4.520913,5.094684],[9.453441,-5.808335,11.14021],[-7.354702,7.084458,10.26066],[-2.505737,7.667425,8.12823],[2.32421,9.505895,9.836869],[3.741202,-4.286661,5.776855],[3.16271,2.165846,3.961518],[-9.434541,-8.305854,12.60943],[6.499024,4.276935,7.844074],[2.17498,-2.170666,3.231459],[2.405977,9.225812,9.586675],[-9.409488,-6.607799,11.54129],[-4.317581,9.747323,10.70756],[-2.293085,-5.737193,6.258883],[3.757021,2.519901,4.633046],[2.855282,-6.063993,6.776772],[8.343461,-6.098679,10.38303],[7.955764,-9.258489,12.24801],[-4.024394,-0.3841014,4.164526],[-7.083542,6.918124,9.951735],[4.833675,8.394437,9.73812],[5.05583,-2.940102,5.933432],[-1.082645,6.879809,7.035901],[-5.895067,3.754334,7.06023],[1.185799,-4.713655,4.962325],[-4.070023,-6.788394,7.977931],[8.315164,5.564547,10.05515],[-3.799644,-5.545042,6.795938],[-0.5556321,2.717601,2.948573],[0.697813,-9.862013,9.937114],[9.569399,4.85715,10.778],[6.351316,-7.730144,10.05457],[-0.2588238,7.303119,7.375807],[9.396074,7.172751,11.86316],[-8.984544,7.629799,11.82945],[-0.9248022,-6.441074,6.583517],[-1.947302,-5.151423,5.597245],[9.733056,-0.04013788,9.784374],[-5.630575,0.2609905,5.724639],[-8.08064,-0.9857324,8.201732],[5.64075,-1.339162,5.883147],[9.004729,5.28183,10.48727],[7.214299,-4.113894,8.364821],[9.563101,8.317739,12.71368],[3.910252,-0.5536227,4.073888],[-9.105172,-2.409622,9.47156],[4.256504,-6.586244,7.905468],[3.356564,-7.598704,8.367008],[-5.277678,3.32234,6.315998],[5.948592,-9.703931,11.42594],[1.619027,6.061372,6.353069],[-9.965819,-3.541204,10.62345],[-3.924944,8.475065,9.393185],[-4.330619,9.054907,10.0869],[-6.263229,-3.795506,7.391475],[-9.186654,2.210999,9.501743],[8.500083,-6.186438,10.56047],[-7.469251,1.400294,7.66489],[-9.455531,-4.467549,10.50552],[-4.28089,-7.482431,8.678294],[0.9238409,-3.029005,3.320896],[-7.879797,6.760355,10.43042],[1.082923,-9.571944,9.684773],[-7.324514,-1.56565,7.556439],[0.6990712,-5.407317,5.543264],[-1.723192,2.579258,3.259135],[-1.119988,-2.531333,2.943131],[7.441556,-1.44042,7.645362],[0.2856756,3.33345,3.491919],[-6.832502,-7.796116,10.41453],[8.321278,-8.687914,12.0716],[5.257184,8.648355,10.17015],[7.909396,-8.182192,11.42396],[6.925184,-3.29527,7.734144],[-6.423004,2.243521,6.876654],[7.252946,-5.434628,9.118137],[3.481032,-1.259807,3.83467],[7.387582,5.195416,9.086733],[-0.7358399,-5.006735,5.158377],[8.238169,3.628763,9.057336],[-7.849524,-7.095903,10.62859],[0.4084393,-2.939633,3.131815],[2.475915,-2.223273,3.474637],[2.484237,-4.228388,5.005067],[2.672413,0.7210938,2.943088],[5.559503,-8.497533,10.20373],[7.05869,7.164809,10.1074],[5.343874,1.415699,5.617935],[4.948767,1.932361,5.405951],[-8.889015,-9.877321,13.32577],[0.3097768,-5.704572,5.799836],[2.855456,0.6082838,3.086039],[-7.777306,5.970829,9.855825],[6.339019,2.823967,7.011273],[4.38928,2.375654,5.090138],[-0.8851867,-8.24357,8.351048],[-4.935272,-7.361081,8.918656],[-0.3454079,3.016747,3.196884],[-8.562421,-2.557631,8.992025],[5.967784,-3.541412,7.011138],[-5.973212,4.921917,7.804135],[-3.101011,-1.373932,3.536093],[-0.8852762,2.235829,2.604351],[-6.651644,0.2874828,6.732534],[-5.439477,-2.336606,6.003969],[-9.93841,9.883267,14.05172],[7.775953,-5.284043,9.454447],[3.118805,-5.941008,6.78399],[-6.627821,-3.398061,7.514974],[6.410026,-3.302818,7.279906],[9.532999,5.566576,11.08444],[4.654861,0.1477374,4.763356],[-1.353864,-7.756852,7.937361],[8.615736,8.607346,12.21955],[-8.804266,5.553329,10.45727],[5.190955,0.4593525,5.306318],[-0.8879007,7.361665,7.482144],[-6.94941,2.677927,7.514359],[4.618527,5.386382,7.165465],[4.805251,7.843675,9.252766],[0.4806498,2.432758,2.673824],[3.177994,-9.376887,9.951163],[-4.120412,5.589307,7.015565],[5.537644,8.3939,10.1056],[7.030504,-0.5521873,7.122703],[-3.002629,-0.6751144,3.235979],[3.615765,6.796436,7.763073],[1.177503,3.677153,3.98848],[-0.9363456,0.1661804,1.379985],[-2.195665,-8.031215,8.385783],[-9.185344,4.001422,10.06886],[-2.214412,5.284571,5.816383],[-9.868307,-3.693821,10.58432],[1.210772,-0.5270593,1.65643],[6.906159,-3.512191,7.812203],[-5.051567,0.9651259,5.239255],[-6.390922,-4.823406,8.069022],[8.361105,-4.991146,9.78875],[9.867305,-8.722146,13.20756],[-8.018476,-7.034538,10.71358],[-8.334758,-2.553138,8.774207],[-3.279745,8.879598,9.518613],[8.280025,9.32918,12.51369],[8.996838,2.416065,9.369123],[9.808517,-6.212262,11.65329],[-8.899082,7.206081,11.4944],[4.579411,-2.002954,5.097335],[-2.257466,-8.840346,9.178664],[5.661117,-4.171127,7.102573],[-9.706571,3.015545,10.21328],[9.247133,-2.036185,9.521319],[2.182299,9.573478,9.869848],[-9.361732,5.658043,10.98433],[9.290874,0.9936484,9.397216],[-9.304005,-6.980574,11.67446],[-8.858253,2.869965,9.365113],[6.853761,-5.708704,8.975709],[-9.56776,-8.77217,13.01895],[-4.722635,6.142035,7.812035],[-3.116193,-8.241952,8.867944],[-1.397599,-2.523278,3.052902],[2.346619,1.314062,2.869386],[-8.473837,-4.406491,9.603285],[1.900305,0.414188,2.186941],[5.424891,2.682921,6.134126],[9.547637,-2.248614,9.859697],[5.486035,6.287983,8.404482],[0.4940209,2.062742,2.344987],[-7.433337,-0.9041066,7.554595],[0.3947026,0.8626363,1.37838],[-3.759305,8.543415,9.387349],[8.295027,-0.5868743,8.375672],[8.455276,8.332938,11.91342],[4.394027,-7.110199,8.417981],[6.650087,5.329773,8.5808],[7.677896,7.868862,11.03943],[6.472447,9.353907,11.41876],[9.21455,-6.979062,11.60238],[-1.798314,8.037925,8.297117],[-1.914138,6.981771,7.308149],[9.024523,9.476378,13.12417],[-8.071795,8.306405,11.62541],[-1.077911,1.995729,2.478876],[-7.155457,-8.362831,11.05158],[-1.594652,-7.542224,7.773549],[6.378462,3.086607,7.15625],[8.818546,7.570209,11.66511],[5.489308,9.402327,10.93326],[-5.106024,-5.434169,7.523408],[8.750553,4.735037,9.999639],[8.284713,1.36708,8.456084],[-9.045237,-7.759229,11.95918],[9.025335,0.9987243,9.135324],[3.472428,5.418213,6.512664],[-2.210507,-9.815309,10.11072],[-8.794622,-8.910691,12.55969],[-6.100082,-7.784616,9.940385],[6.977111,0.06723807,7.04873],[9.907228,-8.108665,12.84148],[-9.012413,-8.793814,12.6315],[-2.468261,-2.407471,3.590018],[-2.934753,6.199063,6.931173],[2.613584,-5.414074,6.094507],[4.666449,1.776078,5.09217],[9.263367,7.876267,12.20023],[6.877477,-9.865987,12.06803],[4.711473,-2.280931,5.329223],[7.243971,-5.148477,8.943262],[-7.087565,0.719841,7.193869],[8.529179,-2.523646,8.950736],[-7.297935,4.273604,8.516075],[0.7457405,5.071621,5.222784],[-8.039132,7.635794,11.13252],[-2.51892,-2.633455,3.778895],[5.692453,3.141192,6.578078],[1.057502,-4.921965,5.132646],[-5.773711,-0.6523147,5.895867],[7.179937,-0.68894,7.281904],[-0.1919281,-2.592232,2.78505],[0.1538497,-0.2567036,1.043823],[0.03407677,-5.107383,5.204471],[-5.008235,-2.108597,5.52527],[-1.714579,-7.740038,7.990492],[2.250948,-5.439478,5.971154],[6.544402,-9.318598,11.4309],[-7.483016,-3.676786,8.397279],[9.136044,-0.8862554,9.233241],[7.460351,-3.955342,8.503033],[-9.107127,7.357878,11.75067],[-5.327539,-9.462308,10.90495],[-7.006706,-8.42065,11.00006],[-0.7100271,-3.661831,3.861754],[0.3741361,-5.086576,5.197425],[0.958651,-5.518213,5.689437],[5.427458,-2.096516,5.903615],[-4.576718,0.4855451,4.709788],[-2.346684,-0.09556466,2.552657],[8.365328,1.285214,8.522351],[5.09964,7.815167,9.385262],[2.974661,-0.2812568,3.150827],[-5.771164,-8.110549,10.00437],[0.7960012,1.23444,1.776924],[7.917594,-2.997341,8.524807],[8.270347,-1.308018,8.432648],[-5.954382,2.308191,6.463932],[1.307313,-1.271762,2.080011],[4.988663,3.867118,6.390724],[3.855012,-1.397443,4.220659],[4.124998,0.6837839,4.299206],[-5.439649,4.295053,7.002661],[-9.270196,-1.217458,9.403124],[9.105769,3.48705,9.801763],[9.65147,-1.255855,9.784072],[-3.063665,4.793519,5.776146],[-9.83349,-6.446624,11.8007],[0.4257514,-3.380692,3.551104],[-0.3831196,-5.636263,5.737094],[-7.355731,-5.95379,9.516007],[2.531679,7.629434,8.100473],[-0.4568526,5.181347,5.296704],[7.337766,1.873867,7.638991],[-7.102428,1.117481,7.259011],[-7.883779,6.161119,10.05551],[-8.987316,-3.924709,9.857748],[-0.515686,3.337099,3.52167],[8.071041,-4.655997,9.371233],[4.276924,4.736269,6.459437],[2.275834,8.343608,8.706045],[1.366447,7.73754,7.92065],[-4.559934,4.329689,6.36704],[-8.447812,-2.664625,8.914356],[6.323609,1.35744,6.544515],[-4.161082,4.542749,6.241087],[7.365647,6.737235,10.0321],[-1.817082,-2.847834,3.523059],[9.727248,2.873117,10.19187],[-5.960362,-9.897388,11.59673],[4.036988,-2.478403,4.841462],[-9.250761,-4.856565,10.49585],[-1.132484,2.118851,2.602316],[-8.344337,2.478607,8.761932],[5.077918,-5.682603,7.686171],[-8.413991,-7.507917,11.32096],[-0.7363116,-1.909726,2.277983],[-3.009999,-4.644807,5.62444],[-5.230474,-8.534769,10.05983],[6.851517,3.140395,7.602984],[-7.53136,-2.349385,7.952421],[-0.4390952,-5.659335,5.763756],[8.499512,3.481005,9.238998],[-3.196815,-7.186887,7.929122],[-3.130087,8.023958,8.670718],[3.970209,-3.335944,5.281201],[-8.186887,1.03436,8.312341],[-4.695768,-3.577466,5.987361],[-6.971942,-8.455154,11.00444],[1.778772,-4.321733,4.779268],[-2.702448,-0.6580647,2.955719],[-4.201855,8.575105,9.601459],[-2.191721,3.359065,4.133638],[4.384438,6.808509,8.1596],[-6.455182,-9.214353,11.29485],[-8.050876,2.949149,8.632154],[2.574652,-5.308305,5.98389],[0.8364726,2.764832,3.056793],[-8.078938,4.322987,9.217237],[3.041436,-8.48842,9.072134],[3.700261,7.359215,8.297588],[3.840542,-9.143258,9.967394],[-9.288502,-6.977803,11.66045],[5.543,3.40953,6.584053],[7.562519,9.385516,12.09461],[1.995383,-7.548348,7.871411],[0.1607402,8.672624,8.731565],[-7.317855,8.158545,11.00513],[-9.871163,-6.564285,11.89663],[-0.5437045,8.12111,8.20049],[-6.555251,-9.907906,11.92216],[-5.179623,6.540651,8.402893],[1.506051,-4.499685,4.849263],[-6.442364,-9.336809,11.38772],[-3.121041,-1.875009,3.775785],[6.176118,-0.7271423,6.298664],[4.145356,-3.882703,5.767093],[5.262579,-8.86446,10.35729],[1.039092,-3.858583,4.119269],[-2.30193,-3.201142,4.067701],[-9.220319,-4.424801,10.27585],[5.756808,1.359505,5.999091],[-8.91712,-5.040965,10.29205],[2.84891,-4.409703,5.344321],[-3.010178,4.648291,5.627413],[4.263401,-1.404338,4.598777],[-3.723753,4.29454,5.771431],[1.067071,8.691072,8.81325],[-3.89783,-0.8290271,4.108572],[-2.487268,-2.863611,3.922598],[-6.660068,5.747478,8.853812],[5.073587,-0.1774408,5.174242],[-9.616675,6.026206,11.39279],[8.820393,5.628447,10.51089],[2.569104,8.063481,8.521739],[-7.297857,2.993166,7.95096],[-3.803301,4.360104,5.871593],[8.732265,-3.495281,9.458828],[7.745537,1.674449,7.98731],[-9.834803,-0.8376865,9.920941],[-1.061328,-5.705637,5.889033],[-9.651464,-1.91869,9.891012],[-3.171285,9.903316,10.44666],[0.0563279,6.719602,6.793838],[9.395549,-2.880148,9.877833],[-9.020248,-8.094055,12.16053],[9.686136,5.317018,11.09468],[5.865831,-6.701387,8.961951],[-0.058093,5.389501,5.481796],[2.435114,9.618417,9.972147],[-9.568503,-2.100797,9.847315],[-8.739391,5.114414,10.17517],[-1.405107,-0.5090637,1.798185],[-4.989345,-8.920851,10.27011],[-1.671831,2.586332,3.237921],[-9.600349,-9.539162,13.57064],[4.902878,8.550792,9.907283],[7.549505,8.80075,11.63822],[6.146396,0.8615125,6.286524],[-7.67172,-1.63717,7.907946],[1.296679,5.000852,5.262119],[-8.552695,7.066598,11.13936],[-2.205918,-2.941158,3.81005],[8.983545,-6.291642,11.01312],[-7.524231,-2.989878,8.158029],[-8.577415,9.709098,12.99379],[-7.745353,-3.269469,8.466399],[-4.024871,-5.874642,7.191037],[-3.883561,4.448117,5.988972],[0.2720359,-3.040547,3.212309],[-0.2212987,9.267292,9.323715],[-3.274174,-2.113061,4.023088],[-6.614052,-8.68004,10.9585],[-5.657049,5.461384,7.926469],[-3.313859,-4.48102,5.662261],[1.382255,4.425842,4.74328],[-1.179213,5.114813,5.343394],[5.542331,-4.534362,7.230344],[-7.995508,-8.392656,11.63464],[8.871778,-2.000868,9.149422],[-7.521232,-1.386672,7.713092],[1.171833,6.72613,6.900291],[-8.170115,4.82053,9.538778],[3.359535,8.777493,9.4515],[-2.30741,-8.876643,9.225992],[3.150689,-5.484211,6.403391],[-4.27342,8.376751,9.456853],[8.452131,3.715903,9.286897],[-4.573811,6.525752,8.031511],[4.730932,-4.457972,6.576871],[3.958866,9.707193,10.53101],[0.1436155,2.801457,2.978051],[4.475457,-1.030222,4.700115],[8.353098,5.265168,9.924527],[4.056514,4.313087,6.004833],[-9.540756,-3.528484,10.22136],[4.457321,8.237097,9.418995],[-5.725309,9.117539,10.81243],[-8.418668,0.1974885,8.480151],[3.153896,0.3782662,3.330187],[-7.570054,-4.379375,8.802536],[-1.707923,-6.775105,7.058261],[-0.1217515,4.211831,4.330628],[6.060525,-7.78697,9.918007],[0.6945292,5.137134,5.279443],[4.200402,8.38931,9.435248],[-9.819388,-7.282794,12.26619],[7.560737,-2.054918,7.898571],[-5.710106,3.833735,6.950024],[-2.280266,3.984998,4.698917],[1.779105,9.735011,9.94664],[9.709916,-9.917535,13.91546],[9.813367,-4.770545,10.9572],[-4.358712,-8.01297,9.176386],[2.477756,8.376047,8.791896],[6.396301,-9.453664,11.45794],[-6.741276,-7.237201,9.940919],[8.94155,-8.540244,12.40512],[5.702064,-0.6252913,5.822759],[7.848375,-1.457853,8.045018],[5.084004,8.064917,9.585926],[-7.603934,9.789922,12.43633],[-0.4134541,1.000094,1.473476],[-6.675208,-6.295908,9.230215],[-5.096969,-9.818861,11.10807],[7.001001,-8.881346,11.35308],[5.881155,9.524705,11.23868],[-2.980885,0.6924743,3.219503],[1.042603,-8.80103,8.918809],[0.3174655,5.015146,5.123716],[8.959971,1.257269,9.102846],[-5.469697,6.930364,8.885242],[9.154168,7.294065,11.74743],[5.526065,-4.839615,7.413452],[8.117255,-8.715508,11.95198],[9.461662,8.866683,13.00543],[-2.414317,6.040425,6.581464],[0.05368252,9.735518,9.786889],[-9.349236,5.23558,10.76195],[-3.267764,1.324418,3.665019],[-3.894006,4.664481,6.157976],[1.267164,-6.683411,6.875586],[-0.2298904,-0.6323171,1.205269],[7.076639,8.347242,10.98887],[-1.885097,7.76726,8.055056],[-6.502326,9.749555,11.76155],[3.678898,-2.589121,4.608453],[-4.798928,-5.144135,7.105761],[3.03733,5.459948,6.327433],[9.670519,5.460569,11.15064],[3.877644,-4.350502,5.912951],[-6.437234,8.539247,10.74042],[1.069013,1.821746,2.336996],[-9.760849,0.1468956,9.81304],[3.611541,-9.904011,10.58927],[1.989929,9.131886,9.399529],[6.723576,9.037912,11.30886],[7.610737,-0.08751834,7.676651],[-2.140759,-3.664471,4.360183],[-2.159015,9.968084,10.24812],[-5.676472,9.424576,11.0474],[6.678045,-8.747154,11.05029],[-8.498635,2.365088,8.878088],[-4.073749,0.2831024,4.204233],[5.991421,8.64347,10.56441],[4.383621,6.577181,7.967147],[1.802249,6.228503,6.560667],[8.890569,3.798294,9.719529],[3.995889,-5.661216,7.001178],[9.778047,-5.933062,11.48092],[-5.48572,3.96854,6.844153],[-6.530749,-1.483261,6.771318],[6.606818,7.501803,10.04625],[5.417221,0.3155325,5.517776],[8.675779,8.12309,11.92702],[7.085883,4.390355,8.395532],[-7.256907,7.879449,10.75864],[-8.343482,-2.234293,8.695157],[6.697455,-1.52475,6.941237],[0.006046193,-5.435939,5.527157],[-0.6282765,-0.243369,1.205803],[-4.192049,9.323976,10.2718],[-4.070932,3.466961,5.439881],[-3.068352,7.229818,7.917389],[1.888747,-9.863667,10.09254],[-2.247804,-1.520722,2.892269],[-5.839566,7.98635,9.943958],[6.335563,9.79397,11.70731],[4.901665,-0.9386669,5.089932],[8.54783,-1.597589,8.753153],[2.812118,3.535459,4.626821],[-0.09436949,9.8448,9.895908],[-2.320908,1.299834,2.841862],[-7.114805,5.868664,9.276943],[4.727992,-6.032749,7.729681],[-0.8061469,-8.341537,8.439853],[-7.19957,-3.346104,8.001888],[9.983344,-2.989781,10.46929],[2.772044,-5.211648,5.987112],[2.986814,2.100066,3.785674],[-8.301259,-5.242431,9.868839],[9.180977,8.360816,12.45767],[2.217172,-9.364913,9.675611],[-7.042462,-9.893325,12.185],[4.494464,-3.126139,5.565335],[-9.062602,5.468362,10.63173],[9.814893,-7.994506,12.6982],[-3.322364,1.395748,3.739815],[-3.028989,2.933724,4.333764],[1.817012,-1.796409,2.743833],[8.227825,-4.385048,9.376873],[-2.050683,-2.450463,3.348144],[8.524322,-5.665154,10.28387],[-8.127131,-2.524321,8.56869],[0.9374768,-5.218096,5.395126],[-1.725521,-5.143566,5.516674],[-8.373802,5.612187,10.13002],[5.335165,-6.947776,8.816778],[0.3416752,-6.444231,6.530303],[5.866871,-0.3149542,5.959813],[2.089319,0.3375451,2.340767],[4.749691,-8.752966,10.00869],[1.776799,-8.748232,8.982681],[-1.621284,-4.996727,5.347507],[-1.683876,-7.576716,7.825731],[2.516871,4.978923,5.667832],[-1.767169,-4.030402,4.512984],[1.014109,-2.238785,2.653408],[-0.6826468,9.000643,9.081717],[-0.1900559,7.803774,7.86988],[6.085577,8.89895,10.82708],[4.495818,0.8640206,4.686033],[5.196871,-8.063067,9.644715],[-9.730252,-7.007998,12.03286],[1.798787,-7.877105,8.141524],[9.463124,-9.206002,13.24014],[1.467239,4.475827,4.815165],[-9.050193,5.327938,10.54955],[6.665246,2.402804,7.155346],[-8.013946,-3.918139,8.976365],[-6.486549,-7.495129,9.962543],[7.926046,4.590718,9.213951],[1.415715,-3.383583,3.801695],[1.627648,-3.870811,4.316528],[-7.144736,3.842362,8.1738],[7.02283,-5.275131,8.840088],[-6.722955,-5.197782,8.556581],[-2.561205,2.655481,3.82248],[0.2445347,-0.6587573,1.222194],[-7.60045,7.513473,10.73402],[2.932808,9.970504,10.4409],[-0.4193568,2.539734,2.761541],[2.763345,-8.258009,8.765317],[-7.818409,-1.436448,8.011923],[2.585891,-7.482489,7.979629],[7.946949,-2.699323,8.452238],[-9.423691,9.352599,13.31454],[-3.776657,7.28552,8.266918],[6.987598,-0.461293,7.073848],[5.485225,-6.808428,8.800135],[-7.51793,6.473206,9.971042],[3.347569,-3.142716,4.699243],[0.3795152,7.255478,7.333893],[8.353857,1.519195,8.549555],[7.884916,0.934549,8.00283],[0.7779847,-2.294533,2.621096],[9.332724,7.178868,11.81676],[5.264858,2.272345,5.820848],[8.682915,-6.8624,11.1124],[-7.559869,2.969754,8.183584],[-1.540396,-5.715488,6.003301],[8.905748,3.171482,9.506349],[-8.752276,-3.099894,9.33872],[5.259305,5.168972,7.441678],[-1.317375,6.53665,6.742645],[-7.07714,1.005969,7.217886],[0.04448572,-4.641075,4.747795],[5.548146,-3.282238,6.52342],[4.574085,6.998615,8.420384],[9.535464,2.508167,9.910398],[-2.542828,3.08561,4.121525],[-5.376565,-3.420582,6.450414],[0.5573341,7.275316,7.364839],[-3.597444,1.242685,3.935209],[-3.034036,-2.06164,3.80207],[-7.089075,-1.116634,7.245817],[-9.931342,8.164815,12.89557],[-1.103872,-8.039215,8.176033],[5.801492,-1.722915,6.133983],[3.585934,-9.767879,10.45325],[-0.2613812,3.058623,3.228544],[9.765979,5.928327,11.46819],[5.067154,-8.258693,9.740743],[9.833155,3.788495,10.58507],[2.705681,5.870229,6.540665],[-3.033623,-4.640895,5.633895],[0.8303958,-2.665977,2.965972],[1.43105,-9.579439,9.737226],[7.571125,3.08919,8.238023],[9.394603,7.092769,11.8138],[-1.108578,-1.271363,1.960946],[-7.270456,9.294601,11.84268],[8.079309,5.896771,10.05222],[0.1225288,-6.99917,7.071308],[4.202919,4.619568,6.324946],[-7.107224,9.033723,11.5378],[-6.447438,9.033813,11.14357],[-3.442085,7.897871,8.673195],[-2.892744,-6.92124,7.567796],[9.609443,-1.787272,9.82526],[2.531256,-8.320899,8.754691],[9.44817,9.898406,13.72029],[-0.4416661,9.36996,9.433516],[-3.212806,7.925436,8.610148],[-9.24003,-9.400044,13.21889],[6.185919,1.330051,6.405828],[0.6978033,-9.095151,9.17653],[1.534233,6.784938,7.027749],[7.258026,-8.493993,11.21726],[2.171502,6.701437,7.115102],[-0.9774226,5.637635,5.808466],[-6.110554,0.1837512,6.194565],[3.456857,-7.246337,8.090691],[-7.192708,-4.507564,8.547116],[8.175223,4.494832,9.382846],[-6.481455,9.852799,11.83583],[5.347162,-9.051738,10.56059],[1.128879,8.839611,8.967335],[6.334618,5.997565,8.780556],[-9.757228,6.237847,11.62387],[-7.047923,6.658036,9.746931],[0.3554588,-7.729477,7.801998],[-1.569239,-4.758078,5.108994],[-0.2862947,-6.050424,6.139185],[9.085858,1.666216,9.291346],[-1.328022,2.174167,2.736905],[-4.711282,7.376837,8.809875],[5.814189,8.94837,10.71812],[-6.467264,0.0189399,6.544147],[-2.98152,6.630057,7.33806],[-9.888857,6.256556,11.74453],[6.452325,6.837381,9.454221],[-4.870245,-8.893534,10.18893],[-3.822709,0.1115255,3.952916],[-4.019165,-8.500865,9.45613],[6.358894,-1.31986,6.570964],[-3.931115,-2.903332,4.988286],[9.09477,7.889347,12.08125],[5.410104,-6.150962,8.252488],[6.449808,6.77903,9.410381],[-9.103885,-4.740682,10.31285],[-1.787893,0.545781,2.120009],[8.988338,-6.252124,10.99451],[6.263295,-3.97896,7.487388],[-9.200784,-8.153522,12.33428],[-1.672017,-1.71004,2.592273],[3.543728,7.130653,8.025225],[-8.614675,4.021636,9.559611],[9.000957,2.153265,9.308801],[-8.109016,-1.139103,8.249467],[-3.146516,0.1120077,3.303499],[2.56628,1.581393,3.17594],[-3.455392,4.47798,5.74387],[-6.648096,8.058856,10.49487],[-9.35427,2.422297,9.714417],[-0.8878292,7.141387,7.265511],[8.781408,0.05091817,8.83831],[-0.6522209,-3.614348,3.806429],[4.352376,2.473576,5.105072],[-7.616106,7.600242,10.80596],[-0.2082577,7.686616,7.754189],[-7.147081,2.936642,7.791318],[-8.581828,1.149056,8.715968],[-5.158416,-6.69232,8.508608],[7.23813,3.773354,8.223669],[-1.14139,-8.561301,8.694748],[-3.914879,4.446796,6.00835],[2.866461,-9.521122,9.993416],[2.607465,9.941417,10.32621],[8.893923,-9.2972,12.90503],[0.3165582,1.885522,2.157638],[2.995021,-5.631243,6.456086],[-6.495008,0.2526748,6.576395],[5.975975,8.62524,10.54073],[-3.961747,-7.484834,8.527495],[8.43512,-1.002244,8.553113],[4.967566,-6.585263,8.309176],[1.45342,-5.02632,5.326943],[1.161562,0.4391521,1.594391],[-0.1787593,-6.628766,6.706153],[-2.047884,-3.052995,3.809804],[3.49123,4.908309,6.105751],[3.438717,-1.447785,3.862752],[-1.829091,-9.598876,9.822627],[-8.947632,6.020453,10.83079],[6.126363,1.458726,6.376536],[-4.202673,2.493553,4.988012],[-4.593033,-2.922624,5.535132],[-8.605408,-9.469215,12.83429],[9.337125,-6.337077,11.32874],[-4.753941,-4.067977,6.336276],[-5.301907,5.276918,7.546926],[-6.424314,-9.432957,11.45655],[-6.760554,-4.480027,8.171641],[7.629502,-3.776235,8.571421],[-6.235689,-1.586765,6.511654],[-1.136132,7.316988,7.471888],[2.742484,9.874319,10.29677],[9.720802,-9.777529,13.82368],[4.059032,5.529339,6.931763],[-5.001196,2.907128,5.87055],[-6.735356,-6.129835,9.161872],[6.354763,0.4406797,6.44804],[9.768438,0.1920953,9.821369],[-4.126097,6.183279,7.500507],[-2.401377,1.689492,3.101772],[4.724698,3.630991,6.042091],[6.708174,-4.079123,7.91447],[-5.457417,-2.731731,6.184315],[2.710279,-9.62024,10.04463],[8.724179,-6.580532,10.97336],[-8.453034,-0.6089199,8.533731],[3.95682,1.271329,4.274658],[5.32605,1.258763,5.563388],[-0.3532307,1.662414,1.9719],[-3.154928,-2.157859,3.95094],[6.873587,-7.360641,10.12053],[1.568869,5.943872,6.22824],[5.849385,9.656276,11.33397],[3.220246,-1.576554,3.722299],[1.001647,-6.673859,6.822293],[-2.500759,-8.906477,9.30479],[-8.969245,1.72965,9.189073],[8.322105,3.98756,9.282137],[-7.114598,-0.009739595,7.184539],[-8.638371,3.267454,9.289657],[-4.973714,4.143791,6.550483],[-3.057006,-7.038534,7.73862],[5.167496,6.676783,8.501908],[7.149391,-2.799904,7.742949],[3.479694,-4.998685,6.172124],[-2.57984,7.905544,8.375751],[-9.345629,-4.396735,10.37651],[-2.919431,-5.876326,6.63734],[-9.787966,-4.307847,10.74066],[-6.643061,-2.349177,7.116804],[-5.91692,-6.800507,9.069555],[-0.719407,9.700235,9.778145],[0.08967771,1.262665,1.613185],[9.308549,-9.968415,13.67547],[-8.082753,-4.721406,9.413957],[-8.245977,-1.416735,8.426344],[4.769582,-3.125557,5.789475],[3.891599,-4.107701,5.746108],[0.7644521,9.936826,10.01623],[-2.136371,-8.575363,8.893871],[4.598984,4.788814,6.714417],[-9.146389,9.111232,12.94878],[1.217703,2.65502,3.087383],[-3.284445,-7.923078,8.634973],[-8.471087,0.7460657,8.562472],[7.285286,-1.61646,7.529166],[8.389912,-6.863554,10.88572],[9.025734,0.5374122,9.09685],[-1.180281,8.293049,8.436096],[-1.376158,-4.843579,5.133621],[-9.705334,4.568925,10.77351],[-7.028446,4.149483,8.222971],[8.344283,2.479819,8.762223],[2.105952,4.307898,4.898267],[1.336579,1.37788,2.164485],[-0.9777235,5.332594,5.512939],[-1.434183,-9.695643,9.852024],[-8.548654,-7.203107,11.22338],[5.407984,6.711101,8.676702],[-5.230673,-2.134977,5.737427],[-1.551727,6.830709,7.075765],[2.156929,1.652314,2.895252],[-5.377623,6.035294,8.145159],[6.155353,-8.575088,10.60285],[1.799313,0.8816051,2.239365],[-0.08841647,-9.264,9.318235],[5.686598,-5.369621,7.88481],[8.336019,-9.838454,12.93385],[5.923166,0.5627087,6.033286],[-3.086261,-8.329372,8.938873],[1.952016,1.502451,2.658519],[7.033173,7.717665,10.48942],[7.365803,9.641955,12.17466],[7.177139,-2.878791,7.797356],[9.015049,5.749077,10.73885],[1.114204,6.516705,6.686471],[2.88548,-7.152742,7.777385],[7.909256,-2.233466,8.279173],[-8.236672,-5.266643,9.827527],[0.05108798,9.262949,9.316912],[8.053933,5.695292,9.914746],[3.461235,-1.086674,3.763112],[-2.362585,1.039229,2.767997],[-6.676546,-7.085462,9.786728],[9.728827,1.869912,9.957242],[6.592977,6.372203,9.223465],[-3.103691,5.101085,6.054252],[-3.525172,6.865457,7.782116],[4.359023,9.0679,10.11078],[-1.630761,-7.875064,8.104074],[5.995005,-5.937586,8.496765],[6.186158,-6.858808,9.290415],[-7.556148,9.927172,12.51575],[8.841129,-8.880692,12.57109],[-9.45821,8.99941,13.09378],[7.053361,8.899969,11.39997],[-6.520841,8.818426,11.01299],[6.645033,-9.071957,11.28968],[9.19484,-8.112863,12.30299],[2.432688,0.8865439,2.775595],[-7.227872,2.353487,7.666879],[8.720655,-8.925817,12.51879],[7.566294,-3.470777,8.384217],[2.886936,-6.156673,6.873065],[6.367738,-7.604116,9.968484],[-4.769611,0.3506742,4.885915],[4.303854,-1.584648,4.694068],[9.965402,4.13878,10.83692],[-8.512636,7.053324,11.1002],[-2.203094,-4.816769,5.390258],[-3.588014,-2.250613,4.351908],[8.801136,2.202357,9.127452],[-2.558908,-6.282612,6.857056],[-9.308473,9.584642,13.39825],[-5.991918,-2.5914,6.604425],[0.5623642,6.222707,6.327586],[-0.6571087,-8.876643,8.956929],[-9.622786,-0.683592,9.698727],[7.648746,0.4170209,7.725104],[9.36125,-4.955989,10.63931],[2.750953,-8.432299,8.925884],[-2.785273,-3.030056,4.235444],[4.209073,3.331489,5.460321],[8.442007,6.878038,10.93503],[-1.007618,-7.805453,7.933498],[-9.179793,-1.923817,9.432374],[9.988093,0.6577474,10.05956],[-0.2984064,-9.596757,9.653331],[6.122307,-2.987323,6.885255],[0.2787552,-4.312842,4.436024],[7.276969,7.039121,10.17367],[-9.319702,-6.869609,11.62103],[-1.397509,0.5471867,1.803454],[-1.59139,1.929206,2.693392],[3.586832,8.910706,9.657434],[-6.80502,-2.046342,7.176058],[4.229972,2.642745,5.086921],[3.822251,-9.495728,10.28487],[7.708927,9.746727,12.46701],[-2.193439,8.969528,9.287821],[-1.089389,-5.339833,5.540811],[-3.71719,7.513505,8.442172],[4.330385,2.46795,5.083602],[8.706561,-1.413793,8.877106],[2.725078,-6.346901,6.979198],[2.886431,-9.144958,9.641666],[-5.701901,-1.516788,5.984339],[9.125038,-0.1299488,9.180589],[-6.35394,2.689199,6.971681],[-2.022666,5.041821,5.523689],[-7.977414,9.369103,12.34582],[7.975768,6.842147,10.55594],[3.228816,-1.648467,3.760678],[6.227774,7.339038,9.67712],[-8.176476,-8.50629,11.8411],[-5.391586,-7.096285,8.96808],[-5.930111,3.837183,7.133736],[9.538124,-1.851256,9.767444],[-6.45814,0.7773647,6.581176],[-5.491472,-5.770138,8.028123],[-2.364752,-7.059816,7.512194],[6.733721,7.232028,9.932031],[9.976801,4.176486,10.86184],[4.379723,3.549239,5.725301],[-7.962668,3.026568,8.576958],[6.291411,-8.900658,10.94548],[-0.6730929,7.737399,7.830734],[8.074352,5.507142,9.824652],[1.724149,-2.291896,3.037347],[8.594369,-3.557268,9.35507],[-6.462662,-6.583423,9.27941],[1.541058,-0.0315937,1.837351],[-0.3924099,-6.33608,6.4265],[-2.9456,-1.561783,3.480765],[3.558359,-2.937077,4.721053],[-5.768657,1.546993,6.055624],[2.444483,1.632054,3.104689],[3.778471,2.283903,4.526926],[2.535997,6.272407,6.839179],[-9.940829,3.06462,10.45045],[8.284396,-9.410103,12.57701],[3.343011,5.572033,6.574441],[1.397865,-1.863913,2.535389],[-3.337565,5.737277,6.712353],[8.09921,-4.957433,9.548473],[-7.24733,1.242639,7.420778],[9.024527,5.593688,10.66449],[7.877152,4.904559,9.332964],[2.986919,1.031407,3.314436],[-9.677384,5.677188,11.2642],[9.230555,-0.4082595,9.293536],[5.113084,5.734038,7.747439],[-7.471643,-7.778997,10.83228],[7.739139,-4.310154,8.91469],[8.200529,-0.03130321,8.261334],[-5.986668,3.295812,6.906704],[9.900335,-3.551423,10.56547],[4.935203,5.698133,7.604271],[6.679972,-6.38053,9.291566],[7.635226,-0.9659994,7.760788],[7.315536,-1.581691,7.551081],[2.87922,4.31873,5.285956],[8.166503,2.875467,8.715508],[-9.76383,9.195138,13.44927],[1.727226,3.843928,4.331177],[-5.012064,4.65922,6.91586],[5.136998,-4.923748,7.185544],[0.2537997,-2.018875,2.267216],[3.602673,1.834312,4.164607],[-4.750686,-7.978419,9.339389],[8.19693,-8.433284,11.80296],[4.231031,-5.941548,7.362311],[-9.704782,-9.70581,13.76174],[8.215248,2.755201,8.722467],[1.861456,4.764777,5.212304],[-7.958154,-8.126219,11.41786],[-3.395,9.26473,9.917725],[8.95293,5.626122,10.62112],[3.61915,1.516023,4.049268],[-1.974789,-9.494608,9.749224],[4.02235,-5.576089,6.947811],[-8.254399,8.762672,12.07972],[9.188443,-6.129087,11.09023],[-2.467096,-7.478632,7.938293],[-6.018831,-6.90888,9.217319],[-0.9246269,-8.153808,8.266771],[-8.681619,-5.875215,10.53037],[5.007256,4.028661,6.504054],[8.868887,-1.117399,8.994761],[-6.903319,9.552498,11.82819],[4.445877,7.22493,8.541981],[7.497507,1.72818,7.758815],[1.589667,8.083004,8.298312],[-3.373686,6.959839,7.798789],[2.497607,-9.614494,9.983814],[0.3086573,-0.3875857,1.116016],[-6.66591,3.472776,7.582514],[7.535618,3.081184,8.202392],[-6.416052,1.917197,6.770626],[9.819991,9.529582,13.72025],[9.186749,-5.752191,10.88504],[0.08668057,4.448658,4.560491],[-4.397523,-5.128588,6.829394],[3.137731,3.527668,4.82595],[7.030819,7.043518,10.00218],[-2.232719,-0.0212764,2.446526],[-9.795912,-2.826422,10.24444],[1.392345,0.03213672,1.714543],[2.7033,2.463945,3.791946],[-5.16085,-3.483878,6.306488],[4.433296,0.5097912,4.573183],[-1.137217,3.851399,4.138422],[5.03279,-2.616686,5.759863],[-8.867289,6.931862,11.29954],[-2.764022,4.76533,5.598946],[2.802314,-6.402376,7.059985],[-2.121371,7.877507,8.219205],[-1.590977,9.286454,9.474673],[-8.816046,-0.6699641,8.897838],[-1.511165,7.416255,7.634427],[2.674253,5.386782,6.096642],[-2.450363,0.8564859,2.781699],[8.622622,9.274506,12.70299],[-8.027534,6.630534,10.4597],[9.735451,-8.171624,12.74968],[-0.7665257,1.042218,1.63517],[-7.717442,-3.626319,8.585401],[4.574441,-2.108612,5.135344],[8.7827,-0.3181025,8.845168],[4.919716,-6.695331,8.368457],[9.190126,1.021974,9.300691],[5.341414,8.733173,10.28586],[4.44009,3.365635,5.660556],[5.662347,-3.5732,6.76978],[1.651683,8.187069,8.411668],[-9.590072,-6.125689,11.42338],[4.541033,-6.895466,8.316756],[-2.195767,5.466592,5.975368],[-9.308893,7.49321,11.99182],[3.285446,-5.056606,6.112563],[-1.483943,-2.245931,2.871636],[-2.51213,-1.271853,2.988044],[1.16326,8.189335,8.331769],[-5.721913,0.4359732,5.824977],[-2.771897,3.870532,4.86461],[9.995577,8.279236,13.01758],[-1.880692,-4.626969,5.093706],[8.364713,6.612001,10.7092],[-5.335392,-8.135152,9.779933],[6.793189,7.190413,9.942307],[9.12758,1.118447,9.250061],[-2.168602,2.384898,3.374993],[1.555179,-2.796114,3.352139],[6.265238,6.2058,8.874973],[5.782384,-2.465317,6.365042],[4.70989,2.887762,5.614466],[-8.6951,6.518927,10.91335],[1.430928,-2.350515,2.927879],[1.996061,0.4968781,2.28717],[-8.098332,8.560811,11.82668],[4.086466,-0.5491541,4.242732],[9.133516,-7.440553,11.82298],[-6.092925,5.946148,8.572071],[-9.833797,2.244525,10.13615],[-0.1375821,-8.464081,8.52406],[9.551702,-5.100273,10.87418],[9.633879,7.942067,12.5255],[-2.577386,5.476626,6.134847],[-3.698343,-7.947155,8.822416],[6.758987,5.999448,9.092706],[6.157216,7.992589,10.13868],[-2.724678,-2.892081,4.097317],[-0.208334,-4.237328,4.35871],[-4.247097,1.053309,4.488574],[-4.107233,-7.94145,8.996443],[3.821431,-1.193496,4.126472],[8.021711,-4.04451,9.039132],[6.607568,-4.007018,7.792056],[-8.146136,0.5319188,8.224504],[3.762984,-4.321376,5.81673],[2.536395,3.067427,4.10395],[-1.669106,6.22882,6.525651],[-0.06719616,-2.348359,2.553293],[4.93654,7.260211,8.836294],[-7.985876,-1.190074,8.135754],[8.974156,-5.601364,10.62595],[3.222809,-8.075059,8.751747],[-4.947002,5.72608,7.632878],[9.239606,9.281075,13.13426],[1.102046,-5.040099,5.255198],[3.736293,-6.597278,7.647481],[0.9232691,5.910995,6.065665],[4.290192,0.8101659,4.479075],[5.860132,1.041424,6.035371],[-9.987302,-7.913805,12.7818],[7.493257,0.1628384,7.561442],[-8.959002,-9.892291,13.38361],[4.858397,-9.849054,11.0276],[-4.517028,-2.604101,5.308944],[0.2754277,8.509908,8.572887],[1.933881,-7.178392,7.50128],[-4.499557,4.019481,6.115737],[8.273193,-1.567188,8.479493],[7.50668,-8.972349,11.74109],[1.729105,-7.248565,7.518743],[-5.081962,-3.962711,6.521458],[-9.400064,1.025065,9.508521],[-6.964202,-6.613914,9.656292],[2.643056,4.718027,5.499593],[1.401744,-9.659814,9.812078],[0.3279737,4.702696,4.819016],[5.917629,-0.03103761,6.001608],[-3.005771,-6.84682,7.544111],[9.658606,7.529014,12.28718],[0.9521365,-4.4302,4.640391],[2.281769,9.334179,9.660919],[8.873232,-9.168705,12.79841],[-5.429903,-8.828409,10.41272],[-7.624604,-3.280442,8.360376],[-4.926692,1.922196,5.382112],[-3.413877,9.663634,10.29759],[1.722756,-3.406496,3.94615],[9.752666,6.841311,11.95483],[2.368335,-4.466183,5.153232],[3.056464,2.812296,4.272117],[-8.36513,1.766829,8.607966],[-3.109254,5.397729,6.308957],[7.545311,4.411655,8.79741],[-5.879078,6.736372,8.996792],[2.104112,7.277143,7.640949],[5.607522,7.028739,9.04696],[7.212952,-4.548128,8.585577],[-3.598277,-8.186029,8.997705],[-8.967277,2.428614,9.343994],[-4.014578,-9.444023,10.3105],[1.005391,8.127069,8.249853],[-4.016331,-6.88186,8.030624],[3.200202,8.972046,9.578043],[-5.649019,8.783291,10.49084],[1.166136,-9.565235,9.687807],[9.387866,-1.83944,9.618502],[2.103669,9.518847,9.799687],[2.327504,-3.712446,4.494389],[3.508416,-4.948055,6.147539],[-7.675206,5.53871,9.517673],[-5.920835,2.441941,6.482234],[2.458244,-3.195405,4.153742],[4.273759,-0.1761405,4.392726],[0.4734253,3.540241,3.709102],[7.449901,7.0696,10.31893],[-3.378277,-4.854597,5.998322],[-2.602253,0.01921602,2.787847],[5.212995,0.7491653,5.36065],[6.749562,1.988581,7.107112],[-9.035769,-8.103188,12.17813],[-6.138606,0.5411665,6.243023],[-5.960494,7.580759,9.695123],[-6.355009,1.782804,6.675667],[4.925193,4.638552,6.839129],[1.001271,6.501187,6.653419],[-3.332393,-9.905554,10.4988],[-8.99926,-6.484553,11.13715],[2.969697,-0.4663791,3.16806],[-4.348559,-8.389857,9.502614],[-0.8416639,-9.91759,10.00335],[4.142294,7.84552,8.928089],[-0.07629377,-3.213423,3.36629],[4.850558,8.469481,9.81122],[-1.200345,-5.942368,6.144312],[4.976585,-7.066247,8.700475],[7.769999,-7.284159,10.69728],[-6.881454,5.5069,8.870195],[0.831877,-2.660263,2.961253],[2.342125,0.6807148,2.636081],[5.510836,-3.155994,6.428811],[-9.445552,-0.5687637,9.515353],[-4.661745,6.335663,7.929218],[-5.721874,-2.923104,6.502644],[-5.456796,0.6304021,5.583371],[9.70949,5.253099,11.08464],[2.766297,0.3725084,2.964989],[-2.697643,7.53075,8.061604],[9.037578,-1.044074,9.152481],[-1.105024,3.466515,3.773302],[4.981349,-0.3324722,5.091598],[-1.258563,7.078665,7.258889],[-1.326885,-4.352457,4.65881],[7.080709,4.938965,8.690789],[-4.172138,-3.929942,5.818176],[5.865673,1.400075,6.1128],[5.779924,-7.636144,9.62903],[-7.033457,-9.869432,12.16039],[7.106021,0.2346319,7.179874],[-1.802717,-2.165224,2.989646],[-4.059247,-2.699739,4.976553],[4.066726,-1.749202,4.538498],[1.512066,-0.00311173,1.81283],[7.164791,-4.746936,8.652609],[-1.814752,-0.8564341,2.242053],[4.325398,-8.249094,9.36785],[-3.541405,-3.075261,4.795704],[0.3424194,7.053362,7.132123],[-9.223133,2.958226,9.737416],[8.395632,1.758363,8.635883],[4.52785,4.677428,6.586331],[-0.2945758,0.7672083,1.294366],[0.5397591,0.836157,1.41085],[-0.0756669,-5.620868,5.70963],[-4.499197,-9.694325,10.73419],[-2.728481,-6.284762,6.924077],[9.595323,4.621696,10.69721],[-2.703941,-8.741934,9.205036],[-5.876688,-1.279203,6.09687],[7.911356,-5.286425,9.567437],[-4.932025,1.046085,5.139958],[7.606833,-3.059488,8.259805],[-3.174063,0.2862914,3.340156],[-1.745291,-2.417857,3.145166],[2.43953,-8.260708,8.671252],[6.550898,-4.358201,7.931468],[-0.1887187,-6.087053,6.171534],[-5.813763,-8.49759,10.34451],[5.140632,4.287502,6.768218],[-0.823933,5.903151,6.043679],[4.662083,-7.369172,8.777227],[-4.819428,-0.2963587,4.930995],[9.863588,-6.758821,11.99883],[0.219812,9.871028,9.923986],[-6.208371,-1.612859,6.491932],[-4.005644,5.337673,6.748032],[3.466263,6.335796,7.290905],[6.925158,7.291933,10.10594],[-9.521317,4.004838,10.37758],[-4.058208,2.781599,5.020592],[4.024561,0.05919294,4.14736],[0.9981067,6.914491,7.057365],[-2.990633,6.932555,7.616049],[-3.742201,-9.852781,10.58685],[5.384836,-9.57519,11.0309],[1.336231,-3.552461,3.924983],[-1.152076,7.141569,7.302691],[-8.123281,-1.050738,8.251772],[1.340369,-9.391962,9.539682],[1.302532,-8.588085,8.743671],[-6.195855,3.711069,7.291135],[-0.4900052,-8.054438,8.131056],[-0.384175,4.239369,4.372624],[4.671275,7.024289,8.494789],[8.028484,8.067039,11.42513],[-9.545441,2.180571,9.842273],[3.586485,-3.863996,5.365942],[-2.693748,5.407343,6.123368],[-9.756125,-5.685719,11.3362],[-8.323656,-6.730685,10.75106],[6.715686,3.911841,7.836002],[-0.5331475,-0.1605588,1.144563],[5.426686,-0.4125124,5.533452],[-3.091128,-3.183656,4.548707],[1.226619,8.797573,8.938786],[-1.643916,-7.290673,7.540317],[-8.73745,3.160568,9.345171],[4.609836,9.810455,10.88557],[-5.145938,-7.250679,8.947235],[3.539937,-9.170094,9.880374],[-1.181191,0.0220716,1.547804],[-7.842568,-2.282121,8.228848],[-9.340626,-6.973419,11.6994],[-7.59757,-8.456885,11.41236],[8.824499,0.8088707,8.917738],[-9.292551,-5.981956,11.09663],[2.679304,-6.723249,7.306213],[-9.773339,-3.998205,10.60678],[-1.231197,4.324911,4.606593],[2.647211,-9.205343,9.630476],[-3.647274,-4.101408,5.578903],[-7.99578,-2.528956,8.445597],[3.346361,2.129011,4.090333],[-7.206783,0.3089704,7.282388],[5.662699,-5.16792,7.731336],[-6.622648,-2.856936,7.28159],[8.087882,-2.338845,8.478445],[7.950572,-8.582885,11.74213],[-5.678049,0.04600915,5.765618],[-6.245456,-4.856081,7.974161],[8.224795,-0.877148,8.331666],[-2.677323,-7.949429,8.447573],[-8.083122,4.858543,9.483792],[-4.492181,3.573102,5.826384],[9.660138,-2.292943,9.97877],[-8.395582,4.852795,9.748611],[3.511057,-1.045152,3.79735],[-5.656837,9.748077,11.31481],[-4.425598,2.004305,4.960156],[-7.981494,-3.708503,8.85761],[8.10352,7.087535,10.81204],[2.535029,-8.444074,8.872923],[1.205052,-3.708745,4.025784],[-4.685283,7.180394,8.631913],[-8.130489,3.515763,8.914339],[5.918495,-6.840697,9.100754],[0.5885875,1.728997,2.082274],[7.349698,-7.711635,10.69988],[3.457119,1.87545,4.058199],[0.7777522,-0.7449881,1.469662],[-3.812182,-3.85495,5.513019],[-6.263912,-0.8475741,6.399607],[9.803756,-3.740854,10.54076],[7.468621,-1.014493,7.603256],[-1.970203,-5.905822,6.305587],[4.240612,6.131842,7.522119],[6.221618,8.014585,10.1952],[-1.50288,-8.030087,8.230489],[-0.4661069,-4.893843,5.016667],[3.192487,-0.9891909,3.48862],[-7.773095,-2.618037,8.262877],[-1.550666,-3.357833,3.831398],[1.418966,-2.742363,3.245616],[-9.234354,2.845675,9.714482],[5.805252,4.190757,7.229342],[-6.437168,0.8579642,6.570634],[-7.715916,-5.376927,9.457627],[-9.270813,-7.393644,11.90017],[1.991181,-4.085368,4.653497],[-7.930744,-2.391779,8.343699],[3.486256,-8.814675,9.531657],[-4.281052,-9.027506,10.04108],[7.867065,5.990905,9.938896],[0.1555087,-1.632546,1.920779],[-1.381177,2.201705,2.784808],[4.546048,0.04054535,4.654911],[7.846531,5.429543,9.594164],[-2.628924,-5.448638,6.131795],[1.400033,3.004203,3.461983],[-8.32434,-5.573797,10.06786],[-7.490795,3.721234,8.423752],[-6.838972,7.89924,10.49617],[-9.674415,-6.652012,11.78319],[7.311444,-2.915282,7.934487],[9.077174,-7.231893,11.64884],[7.693282,2.890139,8.278858],[-6.52078,-1.5308,6.772291],[-5.212349,7.728436,9.375357],[3.317288,-2.135614,4.070042],[-1.628779,9.846558,10.03033],[-0.05581499,-9.821731,9.872664],[3.464185,-6.748479,7.65131],[1.035163,-5.880105,6.053692],[-1.79543,-6.76516,7.070428],[1.410998,4.414864,4.741512],[-4.471932,9.654901,10.68715],[3.524183,8.008802,8.806859]],"colors":[[0,0,0,1]],"centers":[[-3.041903,-5.688744,6.528015],[-1.034308,-4.574435,4.795337],[0.4177045,-2.200593,2.452976],[3.053637,-6.670028,7.403646],[5.728846,0.9114407,5.886458],[1.270159,-5.032366,5.285643],[3.94515,8.675275,9.582515],[-5.65379,3.542814,6.746619],[1.175283,1.488634,2.144136],[-1.725082,-3.531672,4.05569],[5.63594,-2.897137,6.41539],[1.680579,3.728849,4.210541],[8.12744,-9.73701,12.7226],[-5.651456,-3.215321,6.578544],[1.64935,8.409482,8.627848],[-0.1889722,-2.068238,2.305063],[0.9847068,-3.59719,3.861272],[1.296311,-6.403102,6.609096],[-5.089932,5.314185,7.426168],[-5.753172,8.90194,10.64629],[5.721361,9.524649,11.15585],[9.983981,6.752047,12.09421],[5.858741,-6.537451,8.835333],[-7.400498,7.71827,10.7396],[-8.83129,-4.904876,10.15133],[-8.578278,1.856538,8.833662],[6.432083,-4.143786,7.716389],[0.1177704,-7.41595,7.483995],[-9.981523,7.922033,12.78239],[4.191308,-5.047933,6.636919],[1.881472,-2.297105,3.13315],[-8.60346,-9.485167,12.84476],[-2.642886,-3.753031,4.697881],[-2.338262,2.819568,3.797029],[-5.579072,-7.00825,9.013412],[-3.63944,-4.245521,5.680666],[-2.468257,5.327007,5.955611],[-3.468616,8.663477,9.385474],[-1.524321,7.434886,7.655135],[-4.210275,0.9561087,4.431767],[-1.33887,-6.564769,6.774125],[-6.639915,-6.806081,9.560921],[-9.350463,-8.989547,13.00935],[2.690791,2.370103,3.722599],[-0.9226267,3.907552,4.137657],[-8.706333,-6.627686,10.98756],[6.52757,-9.695459,11.73078],[-1.030553,-4.609309,4.827812],[4.218264,-3.032192,5.290362],[-2.917331,-8.066046,8.635503],[-3.246051,5.374498,6.357836],[-4.95193,-6.700486,8.391551],[-4.875647,9.289361,10.5387],[0.9211302,-6.243089,6.389417],[6.027972,-3.41845,7.001589],[4.64023,9.228493,10.37771],[7.125714,-2.027264,7.475667],[-3.529205,6.673452,7.615133],[1.514978,4.820969,5.151397],[8.740385,-8.968707,12.56312],[1.696795,-2.678081,3.324339],[-4.390811,-8.741736,9.833472],[-3.001467,-4.344946,5.374697],[4.113747,-4.614683,6.262444],[-7.375525,9.542572,12.10203],[-4.785921,0.9573082,4.982115],[-6.989492,-5.24707,8.796861],[7.757634,8.367951,11.45441],[8.908677,-8.432468,12.30736],[3.200904,4.847893,5.894731],[-1.955935,-2.076082,3.022548],[4.324643,7.224643,8.479269],[-5.765095,5.771752,8.218846],[-4.65678,-4.080753,6.272013],[-8.449775,-4.514556,9.632233],[3.244244,6.620978,7.440597],[-1.371416,9.490365,9.640944],[0.7077767,-7.230311,7.333372],[3.995197,-9.948018,10.76683],[-0.3973687,-6.421848,6.511377],[-2.379451,5.591092,6.158092],[8.66735,1.109708,8.795135],[-4.772872,-6.160582,7.85704],[-3.817273,0.02893158,3.946189],[4.832943,0.3097725,4.945028],[-5.264491,-1.576612,5.585747],[-0.1962945,-5.440763,5.53538],[-7.834706,-7.423203,10.83912],[4.473332,-5.758923,7.360427],[-6.248059,8.160593,10.32635],[-1.804,-7.334364,7.618879],[4.512907,8.890119,10.02001],[-9.703561,-6.215767,11.56697],[-7.662439,9.079106,11.92238],[-3.927617,4.278402,5.893292],[-7.098456,-0.5332761,7.188355],[-1.573451,-4.533666,4.902028],[-7.94484,7.081223,10.68944],[-2.522551,-9.594135,9.970491],[2.297951,-7.867424,8.256933],[9.046039,6.997881,11.48047],[-6.347957,9.74299,11.67144],[8.336466,8.946384,12.26925],[-2.16044,-8.839568,9.154532],[-3.268503,5.297029,6.304096],[2.87926,2.458652,3.916007],[3.515402,-5.105279,6.278688],[5.007239,-5.90161,7.803938],[1.755195,-2.784528,3.440103],[0.9642107,4.82931,5.025131],[8.735256,5.014625,10.12182],[0.8725876,6.30115,6.439403],[-5.514135,4.856709,7.415747],[1.855446,2.519866,3.285179],[-5.775204,-5.120345,7.782732],[3.501835,8.35325,9.112608],[-0.4028516,9.189212,9.252237],[2.876619,1.594393,3.43759],[-3.579884,-7.429121,8.30707],[9.171999,-8.005317,12.21518],[-9.198751,-4.262588,10.18758],[7.821182,-2.290136,8.210701],[-8.251716,-7.402215,11.1303],[-3.768925,-9.533799,10.30039],[8.996481,-4.582543,10.14576],[7.338439,2.021595,7.677209],[-4.079342,-8.710998,9.670704],[-6.659851,-2.154688,7.070805],[-2.847126,-4.538592,5.450224],[7.875822,-9.548076,12.4175],[-2.999162,-6.379948,7.120303],[-5.063008,3.403615,6.182123],[-8.101959,4.376419,9.262547],[-6.623799,-9.202484,11.38246],[-2.423835,-0.9349811,2.783732],[9.482718,9.408211,13.39539],[0.1720216,-6.812796,6.887944],[7.101717,-0.2749321,7.177044],[9.079535,-2.627288,9.504767],[-3.395019,4.332125,5.594056],[-4.496328,-1.405118,4.815737],[-9.289365,7.138794,11.75817],[-2.362686,-6.400041,6.89513],[-4.148294,6.590561,7.851359],[-6.084896,1.643965,6.381895],[3.934465,6.393528,7.573454],[0.01481946,1.903846,2.150546],[0.3613038,6.088556,6.1807],[0.0797528,-5.184014,5.280186],[9.113288,0.2280227,9.170824],[-4.315263,-7.653095,8.842587],[-4.119615,-7.636724,8.734459],[5.947506,-1.503474,6.215566],[4.304824,-8.050396,9.183702],[2.871321,-0.7515653,3.131986],[-8.075539,-9.024336,12.15125],[-8.697647,-2.59446,9.131281],[-0.2368822,-4.731473,4.841792],[-1.979075,-4.02653,4.596704],[-2.595772,9.137438,9.551482],[4.219622,5.327444,6.86927],[-3.487123,6.053884,7.057587],[8.94285,-6.345466,11.01088],[2.862577,-8.254335,8.793656],[7.891471,0.7374632,7.98869],[-0.4485642,-5.098454,5.214925],[-4.893373,-6.438847,8.148855],[-6.581671,-9.636108,11.71209],[-9.552927,-0.9914612,9.656159],[-4.33479,-3.675765,5.770758],[-8.28151,7.073904,10.93725],[-4.224347,-2.192129,4.863182],[-3.244598,9.130409,9.741241],[-3.208265,2.21846,4.026727],[-7.010926,-6.166512,9.390365],[3.015615,4.883705,5.826192],[-9.420569,-4.863139,10.64882],[-1.230813,-8.168084,8.320606],[-0.03298241,-8.623958,8.681805],[-9.487806,-4.619227,10.5998],[-4.368209,-1.397872,4.694177],[5.465052,-4.451326,7.119066],[7.361863,-4.413455,8.641505],[3.832834,8.037712,8.960772],[-4.35842,9.840309,10.80868],[6.017657,0.8697866,6.161876],[4.106423,-6.989424,8.167911],[5.890417,-7.022818,9.220466],[-4.51641,6.103498,7.658371],[-0.5516427,-1.049928,1.551341],[8.081184,-5.971065,10.09748],[5.441123,-8.686115,10.29827],[-7.101994,3.687819,8.064634],[3.945943,4.397244,5.99218],[8.79126,-6.091478,10.74208],[2.33906,4.027888,4.763936],[-0.4214998,1.049843,1.509911],[6.843839,6.346877,9.387278],[-7.71771,-1.128193,7.863578],[-8.968542,-2.095448,9.264213],[4.912339,4.364938,6.647087],[1.511304,2.531263,3.113091],[-5.537863,-4.330377,7.100711],[-6.719273,-1.211207,6.900409],[-6.664672,0.2856675,6.745328],[1.367334,1.044973,1.99037],[-5.647286,6.869443,8.948804],[9.998017,-2.577188,10.37315],[-1.648835,2.701065,3.318796],[2.909037,-2.338545,3.864102],[5.932885,-9.779927,11.48243],[-3.630975,3.640532,5.238078],[8.672204,2.730684,9.146789],[9.545466,-9.000669,13.15781],[-7.553064,3.261241,8.287609],[2.112064,7.35958,7.721673],[-5.842088,0.9769291,6.007028],[1.406342,6.198298,6.434026],[-6.573358,4.240945,7.886358],[9.782089,-8.95244,13.29795],[4.893617,2.869255,5.760218],[-1.453643,-9.948372,10.10362],[2.517053,-0.6747959,2.791219],[-2.247485,-6.761582,7.19515],[3.212517,8.285376,8.942467],[5.815436,2.592171,6.445048],[3.928644,7.688426,8.691728],[4.926432,-7.568654,9.085938],[3.515446,4.771705,6.010618],[-7.852782,-6.149152,10.02388],[2.403997,-3.701366,4.525407],[1.576189,5.751432,6.046763],[5.552059,-1.358846,5.802743],[7.651042,-6.763002,10.26044],[-7.864532,-3.502975,8.667277],[4.481404,-2.777839,5.366505],[7.179633,2.106048,7.548679],[8.610368,-5.580025,10.30898],[-6.479118,4.756538,8.099606],[5.113476,-3.895552,6.50561],[-5.302091,-3.666935,6.523694],[3.34054,-8.8639,9.525121],[-8.469154,-2.081206,8.778269],[8.700352,9.25927,12.74481],[6.973804,7.850301,10.54804],[2.595432,-0.7370614,2.877417],[-6.337823,-1.769511,6.655762],[1.938265,-3.944125,4.506994],[1.191403,8.493879,8.635128],[-5.122084,2.782767,5.91435],[-8.044605,-3.808623,8.956633],[7.487816,-1.720374,7.747714],[2.908376,5.280206,6.110583],[9.27628,-9.411932,13.25269],[-0.0607456,-3.618773,3.754891],[-6.185519,-3.603159,7.22796],[-1.33394,7.23347,7.423105],[6.221944,6.866237,9.319753],[-4.32747,-7.172585,8.436407],[9.528131,-7.976361,12.46626],[-5.781956,5.891248,8.314916],[-8.504884,-9.186503,12.55886],[-2.177039,-8.971324,9.285696],[-5.918207,-8.664197,10.54009],[-5.110329,3.512228,6.281019],[2.930763,6.837099,7.505684],[-4.580698,-7.214484,8.60416],[0.4965666,7.802157,7.881639],[-3.304909,1.927518,3.954459],[3.891589,-0.4033292,4.03821],[-7.542045,-0.7222106,7.642253],[-0.9479923,2.802906,3.123295],[-9.598319,-3.554731,10.28415],[6.277579,-2.941359,7.004255],[9.260766,-8.1562,12.38085],[4.025662,-8.605541,9.553078],[-0.3855755,5.968174,6.063643],[-9.387498,0.7557794,9.470814],[0.6603565,9.759256,9.832556],[8.15541,-9.301682,12.41096],[-3.350342,-4.180826,5.450147],[8.575846,5.928659,10.4735],[-6.744457,-8.411697,10.82794],[-5.136079,7.59676,9.224428],[-7.314196,-9.131824,11.74256],[-9.044357,6.055414,10.93016],[2.239713,-3.213407,4.042561],[9.491448,5.013898,10.78085],[8.817914,0.4353626,8.885108],[-6.993684,-3.431875,7.854259],[4.518788,-2.354707,5.192697],[0.7050254,3.471202,3.680531],[2.210788,8.020265,8.379273],[4.985639,3.825325,6.363152],[-4.546945,-6.841832,8.275589],[4.777431,5.398863,7.278157],[7.143317,-2.43019,7.61136],[-4.58993,-3.731673,5.999403],[8.838252,-3.212069,9.456854],[-2.37239,7.73232,8.149663],[0.1647775,1.929151,2.179168],[-4.443655,-0.8199979,4.628008],[7.47997,-2.68896,8.011271],[9.848856,4.889863,11.04132],[7.136014,9.667871,12.0578],[3.537357,-4.299102,5.656428],[4.486667,-5.973381,7.537338],[-7.668743,8.414947,11.42895],[1.179935,-9.901007,10.02109],[-0.6659736,4.253392,4.419827],[8.69998,-1.017352,8.816158],[0.6480132,-4.972454,5.11324],[0.5026886,-6.480332,6.576275],[3.592055,-3.946205,5.429124],[0.02951845,4.867666,4.96941],[3.984587,3.587956,5.45439],[9.23109,1.867289,9.470997],[-6.547273,9.475042,11.56042],[-3.10706,7.420079,8.106256],[2.741327,0.8737422,3.04603],[3.468986,-1.028993,3.754023],[9.785803,-3.355498,10.39333],[-1.921161,4.548422,5.037758],[-2.130394,-2.094323,3.150361],[-0.258063,-2.448692,2.657572],[-3.07665,6.457794,7.222803],[5.906697,0.140908,5.992405],[-9.940585,-8.714071,13.25708],[8.09173,-1.896151,8.370872],[5.193403,-0.4298671,5.306244],[2.876321,-7.013688,7.646244],[-7.592597,-9.206544,11.97531],[3.9291,9.819221,10.62332],[-4.733774,1.303829,5.010847],[8.595806,-3.815913,9.457752],[9.984606,2.869265,10.43672],[3.31509,-4.724394,5.85745],[6.730862,-3.493743,7.649231],[2.209388,2.552856,3.521146],[5.527162,-0.778797,5.67063],[5.485619,-4.441967,7.129031],[-6.244771,5.592755,8.442516],[7.337988,3.714547,8.285163],[8.802517,-8.638209,12.37348],[5.081563,-2.895099,5.933286],[5.22298,2.580613,5.910929],[-8.182407,-3.961838,9.145925],[-4.204279,-2.154007,4.828634],[7.918489,3.945226,8.903216],[0.4553644,-2.511423,2.741277],[-6.180319,-0.4399933,6.276141],[6.018744,8.138547,10.17159],[5.616885,-0.1150081,5.706367],[0.6250942,-1.424698,1.849462],[-2.891747,5.9385,6.680418],[8.495856,-3.011796,9.069205],[-7.14734,0.558027,7.238499],[1.269702,1.386196,2.129245],[-3.569924,-1.848289,4.142527],[2.064136,-7.771055,8.102467],[-2.564502,-1.392696,3.084845],[-0.2942988,1.815555,2.093527],[7.139223,-3.765381,8.133057],[-5.792174,0.4365976,5.894056],[-3.287969,8.943099,9.580697],[-9.376264,-8.448777,12.66081],[-6.120036,-4.571502,7.70412],[1.308638,6.930924,7.12392],[3.888079,7.303708,8.334346],[9.93331,8.872841,13.35657],[4.465225,6.238637,7.736849],[2.07639,6.604617,6.995167],[6.458971,-0.6686956,6.570043],[-6.855697,-7.051288,9.885405],[9.812969,3.184412,10.36508],[9.249956,-1.951241,9.506263],[1.344536,-7.558774,7.742276],[-1.506272,-3.372873,3.826895],[-2.156861,-5.083169,5.611654],[9.27997,6.341208,11.28401],[-3.700741,0.08915697,3.834506],[-5.589389,-8.675991,10.3689],[8.945135,-2.457257,9.330249],[6.152593,8.193427,10.29498],[4.42568,1.918779,4.926292],[-6.841035,-5.74894,8.991667],[9.83373,9.569848,13.75806],[7.930473,2.045954,8.250959],[9.981849,1.937051,10.21712],[-0.5635027,-9.190079,9.261483],[-6.40193,-5.863482,8.738715],[-8.741631,-4.808119,10.02667],[-3.938452,-4.665986,6.187312],[3.222997,0.298709,3.387763],[9.827929,4.285499,10.76818],[3.58087,6.008193,7.065481],[-4.438308,8.625249,9.751589],[2.115459,-8.485417,8.802129],[6.481131,8.364964,10.6291],[1.407833,6.602992,6.825064],[6.162591,5.45037,8.287584],[3.230556,9.328724,9.922781],[9.172564,-6.421279,11.24139],[-6.997476,5.168924,8.756851],[8.156535,-5.475665,9.874815],[-3.147497,-7.65501,8.33702],[-0.7827289,-3.566025,3.785392],[-5.902806,-3.265833,6.819735],[-9.248297,1.037181,9.359847],[-4.788025,-6.877434,8.439448],[0.9194502,5.687163,5.847154],[-3.163739,-2.728009,4.295495],[4.246062,-1.368219,4.571769],[-3.840909,-8.296825,9.197277],[6.598755,9.207893,11.37228],[-3.798566,-2.826441,4.839201],[-6.170982,0.03028259,6.251554],[-5.110687,-5.309673,7.437187],[-1.977268,2.199503,3.122083],[-1.347537,9.688792,9.833034],[9.525033,-3.724202,10.27599],[0.1176193,-0.4161661,1.089508],[-4.259114,-2.36226,4.971954],[-0.3741134,0.2189756,1.089913],[-4.457707,-8.344141,9.51293],[5.149736,9.343604,10.71554],[6.102945,-5.238568,8.104846],[-3.320227,-8.776774,9.436932],[-2.376962,0.6968153,2.671235],[-6.660213,-0.0507873,6.735059],[8.124605,5.430187,9.823245],[9.180017,0.7250764,9.262745],[-6.51598,-6.067949,8.9598],[5.239398,3.254287,6.248333],[-2.575474,6.135159,6.728539],[0.327503,-7.575526,7.648258],[-6.371567,1.742287,6.680751],[-9.482172,-8.009757,12.45262],[-2.479106,3.567274,4.457736],[5.811526,4.317883,7.308758],[3.677737,9.434646,10.17538],[4.461787,4.406538,6.350206],[8.482294,9.872469,13.05431],[9.704349,-6.554043,11.75287],[-2.606091,-8.656819,9.095726],[3.696018,8.121447,8.978777],[-6.659905,-0.6786188,6.768668],[5.14638,7.365973,9.041172],[-2.071457,-8.210936,8.52704],[0.4747378,1.152311,1.597873],[5.46438,1.377041,5.723259],[9.771651,1.402305,9.922279],[-2.06874,4.045237,4.652271],[-4.074903,0.0921379,4.196824],[-0.2538389,6.898018,6.974746],[6.386501,-5.566225,8.530548],[6.650153,-3.655893,7.654417],[0.6512699,-0.2850249,1.226944],[8.878991,5.436614,10.45912],[9.991362,-9.927397,14.12022],[3.854037,-9.242485,10.06365],[-0.7051888,-1.1297,1.665387],[-5.708809,3.112792,6.578752],[8.30994,-4.484261,9.495457],[-6.258942,4.780906,7.939233],[8.721509,-8.213727,12.02206],[4.763298,9.071959,10.29512],[-6.833291,7.580168,10.2544],[-4.839773,1.621297,5.201154],[3.145813,1.545001,3.644608],[8.089378,-5.826228,10.01913],[-4.593126,7.384719,8.753906],[6.570402,-4.900798,8.257603],[-0.8645831,8.469284,8.571831],[-4.707701,3.852051,6.164474],[6.644974,-9.409233,11.56241],[3.119814,-9.38994,9.94506],[-8.936648,-7.037745,11.419],[8.945281,8.605769,12.453],[2.722086,-4.962667,5.747853],[2.499687,-9.795776,10.15902],[5.872604,1.863649,6.241848],[7.261394,-4.143202,8.419855],[3.226915,-4.497487,5.624977],[4.581976,3.055317,5.597273],[-1.220884,3.446172,3.790338],[-5.261484,0.3427791,5.366629],[-8.673301,-2.957429,9.218055],[-8.275361,1.425867,8.456636],[-3.850512,-6.168251,7.339875],[5.210942,1.015744,5.402375],[3.764964,6.849532,7.879786],[6.483525,8.703627,10.89905],[9.264662,-9.206747,13.09955],[4.868081,-8.711724,10.02957],[4.243781,-0.6278858,4.404987],[9.707273,-7.326916,12.20307],[5.16256,-4.053179,6.639298],[6.749945,-1.148111,6.919532],[-4.466007,5.279652,6.987127],[6.173522,-0.3672192,6.26476],[9.864346,5.003574,11.1059],[-9.648923,-9.275389,13.42142],[7.341928,7.068431,10.24044],[-3.712182,-6.179323,7.277659],[1.448943,4.955284,5.258733],[3.086207,2.411135,4.04206],[2.266205,-7.941073,8.318434],[6.973834,7.733959,10.46176],[0.2676625,-3.630218,3.774934],[-3.571945,5.506196,6.63905],[-4.100378,-4.859603,6.436524],[-8.31678,-6.671593,10.70883],[-4.672333,-8.22708,9.513966],[-9.521071,8.353239,12.70541],[4.346221,-2.995405,5.372345],[-8.370203,-1.889578,8.638913],[-4.163905,-9.671752,10.57738],[-7.702432,-6.681789,10.24567],[0.4029975,-5.696977,5.798099],[-4.896043,-5.017735,7.081589],[-0.2799949,-8.811684,8.872663],[-0.8567588,-7.850997,7.960666],[-2.519452,6.076905,6.654052],[2.201564,-0.9480549,2.597247],[1.987661,1.118784,2.490477],[-5.444849,6.719845,8.706474],[5.913722,5.172811,7.920232],[1.247939,1.911179,2.491979],[8.314281,9.432113,12.61317],[-0.5717853,7.614768,7.701405],[-6.994335,2.217281,7.405204],[6.644557,-5.732944,8.832711],[5.415682,-1.387766,5.679393],[8.780219,-3.282593,9.426964],[-1.415098,1.249111,2.136067],[-3.484509,-3.532863,5.061909],[5.839069,7.241543,9.355996],[-5.351105,-6.70214,8.634408],[-8.496284,-8.568311,12.10796],[-9.996485,3.074358,10.50625],[3.168969,-6.576552,7.368405],[-0.6096366,-7.786394,7.873981],[-2.318658,-3.568492,4.371534],[6.489139,-2.581017,7.054826],[-8.633159,-6.016428,10.57019],[-3.627624,-7.634729,8.511683],[-1.882741,-6.610703,6.945942],[3.753515,6.211136,7.325782],[-0.008680229,-2.114865,2.339386],[-4.800402,-7.233454,8.738805],[5.865456,8.880745,10.68977],[-2.618692,-2.199729,3.563195],[-5.08169,3.695873,6.362629],[8.107106,-9.597925,12.60339],[7.058955,-6.641448,9.743597],[9.6342,6.033559,11.41147],[-7.021075,-2.893579,7.659523],[9.683857,4.770857,10.8415],[1.733956,-5.074326,5.454849],[6.550436,-1.228368,6.739222],[6.782342,6.730544,9.607309],[-1.92065,5.611259,6.014575],[-7.722196,4.473583,8.98027],[6.377421,6.772868,9.356455],[-2.380948,6.603236,7.09025],[-3.897458,0.6243784,4.071857],[-2.21576,-6.526428,6.964471],[5.9151,9.713268,11.41648],[6.906823,-4.359046,8.228333],[-1.655639,7.454005,7.700866],[-0.6503811,5.488253,5.616397],[-7.166523,-9.708503,12.10843],[4.284157,2.552316,5.086091],[-7.928997,-4.178684,9.018336],[4.298151,2.501638,5.072701],[3.730556,1.062436,4.005723],[-9.915712,-2.881019,10.37408],[3.033545,8.635671,9.207454],[-9.927127,2.93804,10.40096],[-4.085099,0.5455044,4.240944],[-8.247934,-3.590074,9.050804],[-7.74084,6.450017,10.12538],[-4.853656,-5.747512,7.588931],[3.884557,-6.468047,7.610875],[9.324539,-6.358109,11.33016],[-8.045154,-7.532508,11.06631],[-2.137468,-5.235144,5.74243],[6.451472,2.431592,6.966644],[5.301823,0.4973912,5.418185],[-1.377372,2.560553,3.074669],[4.074917,-8.824548,9.771264],[4.941579,7.836434,9.318203],[4.563457,6.559665,8.053219],[-8.815145,1.95146,9.083775],[-9.547643,-8.365757,12.73355],[3.370833,2.486763,4.306565],[-5.890809,2.923168,6.651807],[-5.359407,-3.270872,6.357818],[3.740028,-2.109923,4.409034],[-1.900173,2.354209,3.186371],[-7.728576,-5.601779,9.597438],[-5.060655,-2.988601,5.961708],[8.653278,-7.631227,11.5808],[4.711415,-2.521313,5.4364],[7.794825,-6.265563,10.0507],[2.593345,5.492701,6.155908],[-8.062908,3.213484,8.737103],[-3.241966,4.468866,5.610803],[7.608766,-1.393468,7.799684],[9.547028,5.219147,10.92636],[-6.593624,-2.950086,7.292386],[0.5211073,1.628813,1.981057],[4.105025,-5.613443,7.025808],[-9.198268,-2.164328,9.502234],[6.899184,4.634,8.370944],[-9.441072,5.380912,10.91275],[3.237174,1.730742,3.804571],[1.97566,-9.29222,9.552413],[-7.110626,-1.181263,7.277113],[-9.972219,-4.023222,10.7996],[-2.150389,1.036545,2.588166],[9.210653,6.779508,11.48033],[2.880449,4.573123,5.496402],[5.511725,-9.803714,11.29123],[-6.438959,2.58029,7.00843],[0.4396359,-1.113428,1.559808],[-3.347202,3.255885,4.775411],[2.106855,1.573014,2.813043],[-0.5253502,-6.478419,6.576162],[0.006626537,-3.999259,4.122392],[-9.51368,-6.457925,11.54188],[6.089665,2.136668,6.530648],[-2.674759,-5.526684,6.220818],[-3.215332,2.31397,4.085684],[-9.205165,-3.621442,9.942328],[-2.959229,6.246099,6.983609],[8.89782,-9.810793,13.28243],[4.805547,-7.935825,9.331162],[8.557773,-6.684807,10.90514],[-8.639548,-6.840934,11.06527],[-6.87972,2.088519,7.258957],[-5.017636,0.523186,5.142994],[3.44954,-6.937848,7.812366],[3.747239,9.450079,10.21498],[6.50101,-8.485168,10.73598],[-2.699696,4.197338,5.089794],[8.811143,-4.978623,10.16971],[-0.2257219,-5.900246,5.988644],[2.171347,-9.325655,9.627179],[0.563679,3.934262,4.098311],[-6.016602,8.479056,10.4448],[-9.704993,-5.398738,11.15048],[1.86465,2.894588,3.585465],[3.093728,-0.6539504,3.316445],[4.541197,-7.624595,8.930673],[-1.970688,-6.526587,6.89057],[9.348808,-2.187417,9.653238],[7.979973,-1.59059,8.198168],[5.303439,-0.1274092,5.398397],[-0.1442079,-9.190939,9.246305],[3.945781,-5.511327,6.851563],[-0.5156551,6.619653,6.714589],[-7.239165,4.491143,8.577639],[7.22899,7.909896,10.76219],[-1.040086,-9.173748,9.286518],[-0.1158635,4.013148,4.137485],[8.383542,8.665289,12.09839],[9.657877,-4.631774,10.75769],[-0.6500534,4.95644,5.097928],[-4.213508,5.53698,7.029353],[6.379797,8.240031,10.469],[3.003027,-6.973892,7.658547],[-7.19857,8.06829,10.85895],[-3.558447,-1.626291,4.038239],[0.5741965,-4.617876,4.759672],[-6.607122,-0.3677726,6.692482],[0.3108141,-5.332093,5.43395],[-5.616878,4.506331,7.270236],[-3.416178,4.397916,5.657909],[-4.990015,-4.478836,6.779396],[-2.18823,7.758232,8.122716],[-1.591539,7.932154,8.151813],[-0.8714625,3.2324,3.493974],[1.4787,9.656899,9.820502],[-6.354913,-1.009002,6.511759],[-0.7315499,7.895685,7.992309],[7.202275,-6.842649,9.984719],[-9.614973,-7.082855,11.98393],[7.835856,3.018016,8.456303],[-5.68793,-4.310289,7.206326],[2.582724,-1.908109,3.363234],[-0.2722411,-4.853008,4.962439],[2.539633,4.276136,5.072975],[-7.876253,5.273706,9.531386],[-4.04956,-9.251941,10.14876],[-8.891197,-9.6392,13.15171],[4.785835,3.677992,6.118156],[5.458369,-4.245291,6.986866],[-1.363871,-1.038316,1.984501],[4.99707,-6.338666,8.133228],[-6.775901,-3.152407,7.539928],[-3.162222,1.276424,3.553717],[-2.198236,4.393112,5.013151],[7.866491,2.749083,8.392802],[-6.101921,9.083381,10.98823],[-1.846125,-9.999534,10.21758],[-6.458807,-9.877138,11.84373],[-1.147547,0.5980721,1.635406],[-5.398292,-5.461408,7.743935],[-4.618064,0.2597009,4.732226],[-6.333912,8.216282,10.42237],[3.412006,1.08476,3.717323],[0.01951983,0.7143508,1.229096],[-2.398876,-9.967411,10.30067],[-5.593959,4.258611,7.101278],[-2.304334,-1.279092,2.818871],[-1.991133,-1.682117,2.791797],[1.69195,3.728177,4.214499],[3.15646,-3.161589,4.578087],[-5.376165,-4.374596,7.002874],[-1.197415,-6.318246,6.507998],[-0.7435824,-6.081043,6.207415],[8.767601,-9.76186,13.15921],[2.488117,-5.995066,6.567461],[1.615795,9.754131,9.937498],[0.1692109,-0.4491381,1.109215],[-3.583409,6.059944,7.110819],[-9.155149,5.648639,10.80388],[-3.028139,-0.983506,3.337201],[1.318361,5.445743,5.691589],[9.814194,-8.132664,12.78509],[-0.1375082,9.462806,9.516491],[-9.967447,-7.00067,12.22127],[-3.548011,-0.2580783,3.695265],[-9.64862,-9.644037,13.67857],[3.855137,0.6912323,4.042262],[8.141286,3.109621,8.77213],[8.900586,6.432388,11.02706],[-4.508034,-5.986088,7.560134],[9.895302,1.763266,10.1008],[7.714651,-2.186776,8.080708],[1.598289,-2.493044,3.125667],[8.42073,-3.890986,9.329977],[-1.73958,2.57094,3.261268],[-1.107453,-4.396563,4.642868],[8.055679,9.937401,12.83144],[-5.614483,-9.447405,11.03521],[-4.124797,-8.763474,9.737167],[-9.190564,-2.899293,9.688776],[8.07986,8.668722,11.89247],[-3.644092,6.291048,7.338712],[0.3610228,-6.149076,6.24031],[-1.298302,-2.969804,3.39195],[5.510404,9.403978,10.94529],[-2.936406,-5.096062,5.96593],[-6.138706,1.963304,6.522137],[7.770466,-5.990461,9.862341],[-1.247762,-0.2406318,1.617038],[-4.645027,-9.478976,10.60317],[5.685143,7.116283,9.163097],[5.653201,7.899747,9.765485],[3.011887,1.094175,3.356886],[1.918545,0.2803718,2.18161],[5.492558,7.531768,9.375272],[-3.52991,-6.734656,7.669149],[5.409197,3.279938,6.404483],[7.830295,0.2067829,7.896599],[-1.579557,4.511926,4.8839],[-1.964247,5.57632,5.996133],[-3.038318,-4.046948,5.158407],[-5.879727,-1.362046,6.117709],[-3.067429,3.095804,4.471367],[4.20452,-1.41848,4.548634],[3.806179,-5.344582,6.637134],[6.415355,-5.133148,8.276834],[-3.374391,-2.726172,4.4518],[9.403295,0.1892082,9.458211],[-2.558141,-3.442218,4.403743],[-4.707647,3.467009,5.931449],[-5.288159,-4.685555,7.135758],[9.039901,0.295299,9.099835],[-9.604459,-5.01099,10.87914],[8.636831,7.620301,11.56131],[-0.2397155,6.698038,6.776516],[-6.815151,-6.500652,9.47126],[-1.830068,7.751893,8.027514],[-8.519813,5.09041,9.974943],[0.5539469,4.412366,4.558051],[-2.652797,4.529901,5.343906],[-5.434065,-7.937228,9.671021],[-7.002343,-1.431934,7.216872],[1.120567,-8.21558,8.351731],[7.947289,8.369146,11.58456],[6.226867,-1.99257,6.61394],[-9.384178,-9.056503,13.07987],[-9.330549,0.7833547,9.416623],[0.8363303,7.69295,7.802623],[-8.196482,3.421277,8.937978],[-6.859464,5.086468,8.59793],[-0.7466572,0.837224,1.502811],[-6.24236,-0.2529239,6.327008],[0.3246347,9.860938,9.916828],[-4.73045,-3.265831,5.834621],[-4.059841,-8.825868,9.766179],[-7.60701,5.864132,9.656844],[0.9435672,2.581935,2.925185],[4.627648,-5.615918,7.345315],[-3.565342,-5.155512,6.347517],[-0.235951,-9.500598,9.555994],[3.251614,1.414377,3.684217],[7.885126,8.810743,11.8661],[-8.273724,1.753649,8.516441],[4.546553,6.488283,7.985547],[-2.05517,-9.287406,9.564499],[7.380751,0.5589051,7.469127],[2.286138,8.545805,8.902652],[-7.643538,-5.896768,9.705439],[0.02049373,8.466708,8.525583],[5.743358,-2.5869,6.377947],[9.508097,4.16538,10.42853],[8.173876,9.777085,12.78294],[4.318057,-8.01525,9.159141],[6.750407,1.304072,6.94756],[0.1130702,7.46422,7.531757],[2.772575,-6.888364,7.492445],[9.711888,4.26615,10.65461],[6.38161,7.036725,9.551987],[-1.652751,1.832821,2.66286],[-2.373686,-0.9494892,2.745162],[-3.685802,-6.009079,7.119984],[-3.248025,-6.376012,7.225178],[-3.420874,-6.001142,6.979691],[9.074219,2.729853,9.528565],[-7.028851,8.631681,11.17634],[-7.003469,0.9171225,7.133701],[9.180106,-4.204924,10.14671],[-0.1882417,-2.888589,3.062578],[-0.8728055,0.7364391,1.517937],[1.973257,0.07500994,2.213452],[-1.814758,2.671747,3.381061],[-5.759304,-5.951791,8.342266],[4.71324,-5.12867,7.036895],[-8.910841,3.626144,9.672229],[-2.666356,-5.591552,6.274943],[-9.593958,6.012966,11.36661],[4.855868,0.7376733,5.012347],[-2.941201,9.606171,10.096],[2.771227,-0.8854892,3.076327],[2.985066,8.20958,8.792487],[-4.480999,7.468823,8.767136],[4.23523,-3.185581,5.393061],[-5.899967,6.425003,8.780107],[8.736698,-7.658174,11.66094],[-4.270422,5.71483,7.203873],[5.885511,2.774945,6.583279],[9.099333,6.770942,11.38611],[-9.226543,0.9493978,9.329012],[-5.037013,-6.986801,8.671037],[-1.981286,-0.4949469,2.273866],[5.407158,5.489035,7.769611],[7.918068,6.60015,10.35653],[4.008928,9.621238,10.4709],[6.303428,-2.358822,6.804208],[9.308292,-2.100394,9.59458],[8.769506,-5.41965,10.35745],[7.279232,-2.909902,7.902831],[4.151949,6.083782,7.433108],[9.69567,-9.46055,13.58337],[6.790306,-1.090174,6.949584],[4.606875,-5.780144,7.458777],[9.027784,2.285765,9.366195],[-6.427591,-6.664939,9.313181],[1.124015,-9.327994,9.448539],[1.900753,-0.8173243,2.298016],[6.192993,-4.948568,7.990087],[3.333558,-4.479685,5.672759],[0.06135894,-7.977505,8.040171],[2.917633,-9.112993,9.620771],[-9.090697,-2.175755,9.400781],[-2.792823,1.92535,3.5365],[3.154982,-1.804095,3.769439],[8.117728,0.2101777,8.18179],[-0.8287021,-9.003699,9.096887],[-1.983764,2.847507,3.611594],[-3.369469,3.040385,4.647285],[-1.059345,-6.825825,6.979548],[-5.393869,0.9008496,5.559258],[4.732203,0.1526559,4.839117],[-5.726658,-2.374327,6.279494],[5.452698,8.986813,10.55911],[-1.670591,3.244498,3.783867],[-3.985764,-5.511343,6.87468],[3.667432,-9.032956,9.800222],[-3.440557,4.381483,5.659932],[1.190985,5.97629,6.175313],[6.206937,-6.398638,8.970431],[9.081357,-5.820096,10.83257],[-9.307009,-2.439155,9.673153],[-0.09111921,-3.441918,3.585401],[5.517047,-8.299805,10.01621],[-0.0802413,7.929295,7.992506],[-3.543942,-1.563878,4.000655],[0.5349792,-2.208843,2.48298],[-3.894433,-8.353938,9.271186],[5.235967,-8.237469,9.81179],[2.879822,9.453516,9.932892],[-0.3568964,2.343184,2.572526],[9.95607,1.806061,10.16785],[-5.928872,-1.90463,6.30707],[3.388931,-3.041308,4.662018],[-2.474621,-3.866964,4.698634],[5.082086,-6.26346,8.12764],[-2.367739,-9.714494,10.04876],[-3.501747,6.835564,7.745138],[9.607226,-9.672073,13.66923],[6.708597,-4.585756,8.187455],[8.273551,1.108695,8.40719],[-1.944442,6.901646,7.239722],[6.452111,6.770212,9.405611],[2.709183,8.425428,8.9066],[-0.7923817,-8.696505,8.789599],[-3.65811,-7.247869,8.180058],[-7.349576,5.554318,9.26643],[-6.4676,-5.299846,8.421295],[7.884805,6.217567,10.091],[-4.07881,-1.546884,4.475438],[2.404434,7.92279,8.339779],[-5.605089,5.38428,7.836294],[7.898169,-3.179137,8.572514],[3.270693,0.5111493,3.458136],[3.18966,3.082465,4.547034],[8.760292,-0.3149024,8.822804],[-4.064149,-2.899484,5.091593],[9.257897,-9.026245,12.96849],[-6.744058,-9.765205,11.90973],[4.974947,3.208445,6.003684],[8.716031,-7.911117,11.81334],[3.39794,-4.279674,5.555323],[2.246479,4.616759,5.230787],[-9.588288,-6.107582,11.41218],[9.587842,8.648985,12.95113],[6.507744,7.653815,10.09612],[-5.370017,-2.663354,6.07705],[8.258739,-3.487273,9.020412],[9.171244,8.519852,12.55785],[-9.250583,1.960431,9.508762],[-8.067475,-7.527614,11.07922],[-2.632964,-5.416318,6.104835],[1.122984,-1.745171,2.303631],[-1.96018,1.237204,2.524477],[7.639277,-3.493596,8.459537],[9.989133,-5.761209,11.57473],[-9.452566,-3.74241,10.21551],[-8.568869,-8.911624,12.40333],[-5.287637,7.45386,9.193429],[-6.722613,0.1701316,6.798711],[-8.980945,3.401561,9.655463],[-2.079581,2.852481,3.668965],[-7.133776,2.699265,7.692645],[-2.73675,8.973417,9.434618],[-5.375832,-0.7364092,5.517415],[8.901663,1.017351,9.015243],[-2.197465,7.793251,8.158653],[-5.534708,-5.684674,7.996781],[-5.375805,3.592693,6.542685],[8.77229,-7.914689,11.85729],[-3.753534,5.091156,6.403818],[2.225373,-2.493938,3.488841],[-4.893599,7.084071,8.667835],[-6.465463,4.930681,8.192303],[-4.164371,-6.672914,7.929046],[7.532988,-1.103157,7.678728],[9.314239,3.362468,9.952951],[-8.80388,2.705448,9.264326],[-2.970509,0.9689353,3.280664],[3.915448,-7.601888,8.609264],[-8.65624,-2.292533,9.010338],[5.668028,-5.325031,7.841077],[-0.6850252,7.77432,7.868247],[6.750068,-6.743845,9.593897],[-5.686194,-4.086271,7.073217],[1.618654,-0.8071108,2.066753],[7.810856,-5.359083,9.52519],[-5.889349,-9.321992,11.07176],[5.521703,6.507422,8.592772],[-7.466011,-0.4967357,7.549044],[8.840009,-4.777817,10.09818],[4.235464,7.682645,8.829619],[7.999798,-2.678203,8.495265],[7.489807,4.617035,8.855181],[-5.704573,-9.509583,11.13438],[9.864993,-4.718752,10.98111],[-3.411546,-8.419056,9.138882],[-8.685945,-5.338408,10.24423],[-0.00544942,-5.719899,5.806657],[9.855411,9.527435,13.74413],[-1.381924,-7.261514,7.459176],[0.5713955,6.268704,6.373629],[1.177797,-7.121001,7.286691],[7.472037,9.811078,12.3729],[-1.81141,-1.440739,2.521297],[-7.488823,-1.293421,7.665208],[6.305003,2.229256,6.761852],[3.468614,-4.549651,5.807806],[-8.704941,-1.53771,8.896097],[4.989507,-6.508982,8.262084],[-9.723303,7.74386,12.47036],[6.592782,7.061751,9.712523],[-9.499548,-8.433242,12.74209],[7.42166,7.530411,10.62018],[1.964994,2.93893,3.674032],[4.131686,-5.467602,6.925713],[3.295332,7.661968,8.400295],[1.257933,-8.019286,8.178713],[-5.912465,-2.480287,6.48915],[-7.062162,-9.333306,11.74669],[8.272673,9.868433,12.916],[-5.324579,4.493381,7.03858],[-1.879703,-7.829208,8.113555],[-8.135071,-9.126273,12.26655],[-0.003457479,9.762694,9.813777],[2.305376,-6.843302,7.290099],[4.751412,9.623779,10.77929],[-3.090933,0.9387196,3.381577],[3.350428,9.0967,9.745528],[-8.449457,-3.87362,9.348703],[-5.626029,-7.75807,9.635344],[4.46301,3.285737,5.631566],[7.742088,5.108827,9.329525],[5.738068,1.388753,5.987825],[5.883861,7.41778,9.520677],[-7.81236,-0.4013245,7.886319],[5.901677,3.200914,6.787904],[1.075541,0.4974476,1.550562],[1.903574,-5.622466,6.019611],[4.570216,6.595976,8.086642],[4.626105,8.543122,9.766564],[4.630077,-2.791591,5.498236],[-5.925111,2.998719,6.715598],[2.815404,-1.812589,3.494564],[-7.91539,7.387867,10.87354],[-2.723855,8.178656,8.678122],[8.161918,-0.1731985,8.224773],[-9.870852,-7.164104,12.23757],[-0.6171274,8.892155,8.969462],[3.817678,6.292655,7.427798],[3.251376,-5.173341,6.191518],[-2.066576,-7.165576,7.524375],[5.696543,1.410834,5.953238],[-6.429852,-9.541462,11.54913],[-6.212714,-6.056158,8.733548],[9.427762,2.606971,9.832547],[-5.754764,7.469862,9.482412],[-3.947596,-6.655341,7.802376],[9.505708,3.204309,10.08098],[8.373031,-6.633307,10.72886],[-6.220325,9.296812,11.23046],[2.200725,6.536607,6.969248],[8.695606,6.098835,10.66815],[-1.796201,-2.207254,3.01634],[-6.211329,-3.729821,7.313834],[-1.574739,-0.517197,1.935793],[-1.572644,-8.880605,9.074048],[-5.206865,6.195558,8.154531],[-8.877953,5.409533,10.44419],[-1.04114,-1.245311,1.906507],[8.181485,-2.387106,8.581082],[-9.838553,6.3787,11.76796],[3.362804,-1.300751,3.741712],[4.996677,-9.470849,10.75471],[8.542764,-8.880493,12.36293],[-1.970145,-9.152161,9.415069],[2.557416,3.29369,4.288213],[8.162758,8.717281,11.98422],[-2.759927,-8.417048,8.914252],[-5.20458,6.487766,8.377276],[3.319734,7.991781,8.711441],[-0.6204314,8.566844,8.647297],[-6.480998,-0.3294943,6.565965],[-8.687576,8.218565,12.00078],[-1.335947,9.278155,9.427031],[-0.7707871,-8.356985,8.451823],[4.49949,5.403916,7.102655],[-0.7795951,1.783215,2.188064],[-2.852823,3.523244,4.642397],[9.409116,-7.916888,12.33728],[-6.608869,-1.630788,6.880161],[8.100712,-7.228788,10.90307],[9.295043,-3.088643,9.845687],[-2.760188,-0.8609962,3.059404],[-3.168357,-5.176823,6.151259],[7.93496,7.199534,10.7609],[8.949181,-3.569728,9.686629],[2.639462,-4.539268,5.345252],[-6.606663,1.949112,6.96039],[7.265266,0.4672716,7.348635],[4.436617,-4.954403,6.725302],[6.012594,1.716501,6.332271],[-6.264765,6.522858,9.099174],[5.191232,-8.67757,10.16116],[3.857165,2.011643,4.463678],[-9.798589,-9.291443,13.54043],[-2.116024,-1.80201,2.953777],[4.078426,-2.45347,4.863443],[-4.494982,1.810789,4.948112],[5.259814,4.504542,6.996895],[6.484181,2.080446,6.882795],[7.520535,-8.636292,11.49539],[0.444719,-6.175886,6.272108],[9.576776,2.208214,9.878808],[-6.16544,1.099843,6.342106],[-6.030856,1.933008,6.411532],[7.901597,3.515479,8.705965],[2.572615,-4.882704,5.608845],[5.780593,-9.922232,11.52675],[7.94468,9.948829,12.77095],[-8.909275,3.96881,9.804419],[-1.007452,-8.681378,8.796664],[6.777146,-0.06794503,6.850863],[4.532455,-7.402934,8.737653],[3.178219,-7.754284,8.439787],[-6.610613,-6.221748,9.132927],[-0.9663562,-8.458042,8.5716],[2.293383,-6.853894,7.296264],[-8.692527,2.556907,9.115799],[1.860486,-6.610416,6.93967],[7.821259,8.802216,11.81741],[3.022288,-4.371338,5.407663],[-3.210806,5.752988,6.663794],[3.383077,6.088171,7.036408],[-3.040464,-0.8164917,3.303192],[5.950975,1.066282,6.127892],[-9.485615,4.456438,10.5279],[-6.791014,-1.407247,7.007012],[4.380429,-1.794514,4.838227],[-6.257799,-1.064268,6.425941],[-2.367898,2.007349,3.261348],[6.539535,-2.019082,6.916806],[-5.988785,-8.150484,10.16346],[0.363221,9.150419,9.212063],[-2.240979,-1.605507,2.932514],[0.2511202,2.904192,3.081784],[-9.606741,-4.580312,10.68966],[9.373603,-2.85262,9.848953],[2.363446,6.968741,7.426253],[-8.263205,9.625187,12.72497],[3.141925,5.289403,6.232935],[-1.119785,9.248422,9.369484],[-2.32428,-7.629002,8.037659],[9.126457,-6.820256,11.43714],[8.5087,-7.842581,11.61482],[1.720424,-9.883305,10.08165],[-2.037954,7.393804,7.734442],[-1.499701,3.226796,3.696121],[3.514142,-8.813437,9.540748],[-1.73907,-3.30899,3.869596],[7.472314,2.574836,7.966508],[1.064879,-1.790243,2.310614],[-0.9414656,-0.05657198,1.374612],[8.378823,-1.623291,8.593006],[4.429738,-2.444871,5.157516],[8.164859,7.43662,11.0891],[2.50984,-2.977719,4.02071],[-5.904207,3.419781,6.895982],[-9.370181,7.358607,11.95614],[-5.453716,-0.7137897,5.590395],[-0.1486552,-9.346381,9.400901],[4.71305,-6.411644,8.020101],[5.09884,-3.738369,6.40106],[-6.937679,9.874072,12.10903],[-0.6793345,-7.2308,7.331164],[1.468221,-5.873801,6.136548],[6.601977,0.3127479,6.684602],[7.167895,7.989069,10.77979],[-1.5683,-5.819903,6.109897],[-1.184008,9.652411,9.776037],[4.113077,-9.313243,10.23005],[-8.496347,8.991637,12.41118],[-2.570241,-7.977763,8.441021],[9.20792,-1.4525,9.375262],[-3.259465,7.283296,8.041798],[8.160625,9.926207,12.88896],[-8.50174,-7.139536,11.14686],[-4.335377,4.12097,6.064477],[-7.064709,-8.5657,11.14815],[1.623756,-9.532433,9.721309],[-2.480522,9.775187,10.13446],[-0.7038638,-1.367937,1.83485],[-0.05145437,2.935561,3.101639],[2.619832,7.534458,8.039376],[-0.4654122,8.259791,8.333113],[8.820605,9.118629,12.72606],[-4.439463,1.902761,4.932477],[-9.60649,5.309828,11.02175],[-6.796242,5.787354,8.982337],[-5.292218,7.27561,9.052186],[-5.07464,-8.274911,9.758388],[2.973712,-3.746497,4.886635],[2.230644,-0.7734568,2.563982],[-6.010297,-9.361705,11.16983],[-7.365532,-0.4497307,7.446699],[-8.878706,0.7979292,8.970402],[0.07744554,-5.128317,5.22548],[9.415407,9.107209,13.13739],[-8.449056,-5.580894,10.17511],[-9.804804,-1.386515,9.952719],[-9.399762,6.104842,11.25276],[-9.004338,1.705922,9.218908],[-2.914934,8.57675,9.113588],[-1.327777,8.309918,8.474535],[-1.51345,0.5454924,1.894226],[3.188648,6.806147,7.582289],[-3.982998,7.786207,8.8028],[-9.04269,4.807896,10.2901],[-2.064922,-8.382483,8.690795],[7.103187,9.098195,11.58587],[-5.520591,-3.78564,6.76816],[9.649742,-5.729366,11.26691],[-8.850355,-8.21248,12.11502],[7.550882,2.737827,8.093919],[-8.327228,1.812205,8.580606],[-9.150446,1.844343,9.387878],[-8.131415,-7.07409,10.82417],[-6.962738,-8.815537,11.27801],[9.3894,3.022621,9.914488],[6.296584,-6.044283,8.785233],[-2.624699,3.153152,4.222726],[5.494321,8.930143,10.53257],[4.222104,-0.7573259,4.40451],[-2.864507,-3.711174,4.79356],[1.816288,-0.9339371,2.274014],[9.215417,2.490663,9.598297],[-8.940159,0.3440569,9.002489],[-6.428889,-0.9691069,6.577977],[4.733486,-8.174622,9.498964],[5.184879,-0.0189466,5.280467],[0.8592133,-5.006399,5.177092],[4.759046,2.282806,5.372124],[7.407982,-6.263682,9.752534],[6.483167,7.963304,10.31725],[1.266945,7.5522,7.72275],[-3.223343,6.045249,6.923509],[-8.500842,1.312464,8.659496],[-9.510337,-7.433927,12.11238],[-7.762035,9.529612,12.33137],[-3.939458,-5.178607,6.583107],[2.572414,1.094069,2.968889],[-5.229203,1.856552,5.638382],[-1.256112,-0.319253,1.636992],[5.780801,-4.52494,7.408963],[-1.478089,-2.26287,2.881897],[-0.05919762,-5.878773,5.963512],[0.4364631,0.06928477,1.093298],[-2.56362,0.04059999,2.752053],[1.128901,7.710281,7.856389],[1.409576,-1.320544,2.175027],[-5.241354,-3.339791,6.294918],[-3.572006,7.986213,8.805614],[4.076373,9.023294,9.951716],[7.568694,-5.990261,9.704039],[-5.427725,-6.362413,8.422618],[0.05472959,4.535043,4.64431],[2.47958,-6.240643,6.789251],[1.61984,-3.243638,3.760993],[9.511855,-3.893113,10.32626],[3.015509,-7.981602,8.59065],[-8.864993,-6.654546,11.12974],[1.503707,6.686837,6.926393],[-2.907784,4.052628,5.087141],[-2.299131,-8.970528,9.31431],[-3.333631,-6.537766,7.406449],[-9.501578,-4.781371,10.6837],[-5.376635,-3.660203,6.580676],[-5.662708,-3.261605,6.610924],[2.991921,7.301076,7.953446],[-2.698103,-0.1608225,2.881948],[-2.31315,3.113932,4.0059],[6.064968,2.171185,6.51904],[-2.717006,1.248646,3.152973],[6.674471,-9.193753,11.40498],[-5.659673,1.169478,5.865116],[-7.514844,-7.819026,10.89082],[-1.153888,2.922908,3.297704],[-1.97773,-0.03683682,2.216478],[5.787683,3.293335,6.733746],[0.6199435,-5.103593,5.237461],[4.351086,-9.075044,10.11377],[-2.903973,0.8470521,3.185994],[-2.570694,0.05032901,2.758804],[-9.709691,-8.681232,13.063],[4.258589,-5.923326,7.363516],[-2.93714,5.545444,6.354427],[8.270168,6.306636,10.44841],[-6.410083,-5.573445,8.552921],[8.447042,2.859041,8.973664],[-4.130958,-4.889364,6.478479],[-7.046788,1.014657,7.18935],[6.77867,1.200614,6.956424],[-3.285796,2.308695,4.138421],[3.915301,4.075672,5.739398],[5.389317,-8.195727,9.859751],[1.336507,-9.284917,9.433765],[-2.712457,-7.167543,7.72859],[2.439994,5.251873,5.876712],[-0.780812,1.474106,1.944905],[-5.978556,0.8224179,6.117148],[8.581075,2.349345,8.952891],[-9.17388,-2.383042,9.530948],[7.071238,-8.985312,11.47773],[6.328835,0.2006792,6.410494],[1.244294,4.740616,5.00217],[8.115897,4.705912,9.43469],[3.441625,5.118588,6.248578],[9.724601,-3.037157,10.2368],[6.417716,-3.339608,7.303428],[-5.739388,-7.077181,9.166628],[5.900565,0.2259179,5.988966],[1.833227,-0.3877175,2.123922],[-8.40111,3.045663,8.991924],[-9.746615,2.022884,10.00443],[1.707681,9.589739,9.791796],[-4.541277,7.75071,9.038623],[-6.593448,9.466991,11.58005],[-9.961249,0.4303156,10.02056],[-7.131291,-9.578871,11.98374],[5.513147,-8.100536,9.849542],[-6.330183,7.773947,10.07499],[-5.519099,5.409619,7.792588],[-7.596915,-7.660318,10.83483],[1.848001,-3.307423,3.918438],[-2.856755,2.368268,3.843143],[-3.86025,9.45337,10.26001],[9.284269,-9.669873,13.44262],[-6.845943,-0.2644165,6.923645],[5.952668,-3.602175,7.029219],[-2.436067,0.3849398,2.661316],[-2.413492,-8.405097,8.801739],[3.312071,-6.339179,7.221842],[-2.976202,7.152602,7.81137],[-6.823694,2.270509,7.260717],[1.20349,0.5668027,1.664227],[3.744715,-3.788441,5.419887],[-0.516924,3.486955,3.66416],[3.15053,1.797049,3.762342],[-5.768541,3.403114,6.771798],[-0.4776953,-3.958258,4.110474],[3.254518,-2.755924,4.380298],[-7.696437,6.984113,10.44093],[-7.122421,6.713025,9.838373],[-2.469656,-8.394485,8.807188],[2.728012,5.218865,5.973157],[-4.643526,8.923222,10.10872],[-3.877324,-3.67028,5.431813],[4.691253,3.222629,5.778685],[-9.615191,-2.34368,9.947097],[3.003509,-8.839321,9.389071],[-3.419664,-0.9618399,3.690425],[-3.574147,1.498411,4.00247],[-8.425005,-2.42517,8.823954],[-5.362371,7.489447,9.265357],[-1.450054,4.533877,4.86402],[7.035821,-3.751331,8.035873],[-7.047315,-3.979764,8.154947],[6.668497,-5.592751,8.760577],[-5.897414,-7.433497,9.541298],[-6.229954,5.712007,8.51113],[-2.874685,0.6963587,3.122295],[-6.970351,-4.091882,8.144279],[5.479525,0.2320946,5.57486],[-5.254241,-2.076076,5.737345],[-1.42186,8.886942,9.055353],[-5.97752,-7.350823,9.527084],[2.943729,-5.792251,6.573866],[7.398184,-9.034913,11.72019],[0.3650422,-1.960165,2.230584],[4.78452,-7.23552,8.731803],[-5.773848,-6.110262,8.465969],[-3.237452,0.5294455,3.429491],[-8.235762,-0.857227,8.34042],[-6.960609,-1.084925,7.115275],[-5.829744,-3.206735,6.728229],[9.50271,8.321171,12.67057],[-6.699147,-4.225137,7.98313],[-3.500146,2.310258,4.311417],[4.647257,6.48816,8.043209],[9.634215,0.5225286,9.700059],[7.883838,-8.245569,11.45183],[0.9985148,7.23731,7.373987],[-8.165297,5.130656,9.695138],[8.936748,9.703244,13.22945],[-4.993115,-6.506947,8.262661],[-1.468014,0.7796824,1.939837],[-6.450632,-2.707376,7.066862],[-0.2251106,0.559646,1.167852],[1.964695,5.817009,6.220742],[-7.20397,-6.444257,9.717285],[-0.9189651,1.430779,1.972721],[7.31983,-9.054502,11.68606],[-6.966276,6.096881,9.311336],[-6.87423,4.893336,8.497046],[-9.643953,8.767756,13.07208],[-5.75395,-2.344052,6.293053],[2.275598,7.839189,8.223821],[0.3471159,-6.377162,6.464417],[-9.472334,5.310347,10.90527],[6.704564,7.175653,9.87123],[-4.569597,2.528271,5.317271],[-5.245763,-5.730927,7.833362],[-8.490084,-8.418215,11.99783],[-3.702251,5.767375,6.925986],[4.451914,4.39353,6.334244],[-0.9171414,-8.818001,8.921787],[9.988764,5.867216,11.62754],[-6.139517,-8.345427,10.40864],[-4.852426,1.146938,5.08542],[-6.695445,-8.488696,10.85758],[1.863106,6.375288,6.716805],[-7.908834,6.119433,10.04973],[-5.168276,-4.509001,6.931246],[-8.166243,-6.039346,10.20594],[-4.562191,-2.822143,5.456929],[-8.50326,9.708751,12.9447],[-7.524321,-3.374022,8.306589],[5.990341,6.79079,9.110379],[0.6617281,-5.468819,5.598737],[9.741453,-8.815157,13.17585],[-7.239854,0.4968413,7.325458],[-5.575512,9.607954,11.15344],[5.762807,0.1593744,5.851098],[4.456146,-2.004892,4.987668],[2.681029,6.663498,7.251905],[5.994745,7.110694,9.354087],[1.844592,-8.083296,8.351179],[-4.01845,-7.24964,8.348965],[5.1214,1.832215,5.530438],[0.6415901,4.653295,4.802582],[-6.983442,-4.894094,8.586072],[5.709424,2.767353,6.423065],[7.426394,3.101728,8.109997],[9.197556,-2.512491,9.586847],[1.109242,-8.279277,8.412897],[3.108273,4.565097,5.612617],[6.818488,2.73529,7.414417],[-5.297457,-2.256478,5.844206],[1.126267,0.7037817,1.662464],[6.856763,-6.24603,9.328885],[-4.129322,9.575757,10.47599],[9.171185,1.767037,9.393246],[8.486213,-2.455209,8.890661],[-5.740573,-3.728183,6.917625],[-0.8691576,8.920292,9.018151],[-2.212044,4.2179,4.866602],[-9.138198,-9.324805,13.09422],[-1.90288,0.8802238,2.322874],[-0.1424326,-2.569984,2.76136],[-3.028126,6.604321,7.333935],[7.444137,-9.589713,12.18104],[4.359885,3.196768,5.497993],[1.594765,4.517179,4.893688],[-0.1869513,-5.010824,5.113053],[-9.057297,-6.362833,11.11397],[-0.8831979,-9.650416,9.742206],[3.117315,8.181191,8.811898],[7.713401,-2.15249,8.070301],[-0.7803711,3.8561,4.059371],[3.60577,-7.332383,8.231976],[2.165248,7.277583,7.658427],[-7.157082,6.967681,10.03855],[-3.900046,-9.214128,10.05537],[-1.712709,-0.2780989,2.002676],[-8.981532,2.870923,9.482095],[-9.492363,-3.388485,10.12851],[3.524731,-1.638453,4.013509],[1.827759,6.941776,7.247686],[-3.462868,3.244987,4.849886],[-3.974323,-5.087107,6.532526],[8.485037,-6.65832,10.83185],[8.36235,-5.661519,10.14799],[-7.290123,8.407651,11.17294],[4.487246,-2.226798,5.108229],[9.676881,1.130774,9.79391],[7.465437,-1.842117,7.754106],[-9.735509,4.827431,10.91257],[-4.189157,-0.1835597,4.310769],[5.994042,-6.47486,8.879885],[4.61795,-6.108637,7.722753],[4.213061,-6.141859,7.514807],[-8.915042,-5.561236,10.55487],[-8.932236,-3.678375,9.711605],[-8.668296,-2.471481,9.069044],[2.537138,8.381706,8.814197],[5.258278,1.358679,5.522273],[0.5611183,-4.862759,4.996126],[-8.406111,-7.371167,11.22483],[4.507922,-0.1913222,4.621468],[0.1436699,1.726958,2.000756],[-4.287156,-7.97215,9.106859],[-0.3479221,-9.045591,9.107347],[-4.788922,-0.5541776,4.923503],[-3.192727,-8.446579,9.085053],[-5.152153,9.493639,10.84776],[-9.157755,5.748074,10.8584],[7.393227,-3.337508,8.173051],[1.740336,8.484138,8.718335],[-8.709914,9.21586,12.71985],[-3.523593,5.145523,6.31602],[-2.584728,1.768667,3.287704],[2.692714,-3.724238,4.70326],[5.285431,-7.723823,9.412397],[-3.549942,-4.904076,6.136126],[-0.7855898,5.405947,5.553504],[7.297949,8.894229,11.54848],[4.691267,8.053204,9.373477],[-4.145994,-4.334078,6.080584],[-9.684891,3.725404,10.42477],[-2.705096,-9.38273,9.815964],[-3.428976,4.797397,5.981045],[2.509857,-6.245091,6.804451],[9.735221,-8.672724,13.07634],[1.347327,7.629707,7.812024],[-5.49386,2.96741,6.323608],[-2.563181,-0.09926005,2.753134],[5.754973,8.736856,10.50963],[8.592189,-9.303908,12.70387],[6.821752,6.367276,9.385015],[-2.408646,4.054587,4.820918],[-5.128446,6.404586,8.265572],[7.877079,3.977938,8.881012],[0.2965428,-6.636016,6.717488],[2.215573,-7.420852,7.808829],[8.678236,0.08833563,8.736108],[2.266876,8.293337,8.655528],[-8.378764,6.541925,10.6771],[3.31499,4.230096,5.466523],[9.804229,6.593243,11.85722],[6.61696,6.852157,9.577903],[9.537497,-3.864853,10.33929],[-7.079579,-6.102604,9.400118],[8.606421,1.763846,8.842038],[8.705057,5.797539,10.50664],[7.260916,-5.860518,9.384379],[-8.624317,-0.5296321,8.698238],[5.132399,0.6523948,5.269453],[4.317173,0.3672789,4.44667],[3.317909,-3.711586,5.077833],[-6.307526,2.275755,6.779671],[7.020151,6.480116,9.605958],[-5.634207,4.193324,7.094241],[-6.330136,4.603362,7.8906],[-2.184178,-1.912928,3.070819],[-3.560078,-1.854693,4.136912],[1.621059,-1.050693,2.175267],[-2.905499,7.01885,7.661996],[-9.148784,-4.873044,10.41378],[0.1969298,-6.755055,6.831512],[8.266631,7.837508,11.4352],[-6.909878,-7.750872,10.4318],[4.731191,-1.919491,5.202751],[4.870971,-8.030171,9.445106],[7.319163,-9.988521,12.42339],[-9.459288,-2.182966,9.759275],[-5.810571,-3.694773,6.958023],[4.422137,8.681581,9.794138],[-8.659067,-0.4776238,8.729695],[-0.02197551,2.611602,2.796596],[4.521853,-6.045698,7.615617],[-6.907262,4.784529,8.461796],[9.620697,-0.7951117,9.705154],[5.889467,3.632259,6.991361],[-8.995126,-7.08438,11.49351],[2.87243,-3.356023,4.52921],[-7.35588,2.786874,7.929416],[-7.905954,3.364197,8.649967],[8.218412,-1.476289,8.409621],[-6.601443,6.675242,9.441288],[3.019452,3.926336,5.05304],[6.758572,-0.9130701,6.892894],[-0.2999273,0.3041517,1.087412],[-9.96429,-2.997817,10.45342],[-8.582337,-2.643689,9.035795],[-1.250003,-7.663529,7.828933],[-3.36399,5.948492,6.90659],[5.896552,-0.9787188,6.060298],[1.493479,5.500616,5.786818],[-7.939598,-9.225152,12.21232],[-9.148705,1.264726,9.289689],[-9.007415,-1.522298,9.189718],[0.3650674,4.512632,4.636499],[8.787636,8.174843,12.0437],[2.110367,-0.4628544,2.380732],[-7.948184,-9.63125,12.52735],[-0.2728362,7.542438,7.613331],[4.787845,2.071367,5.311687],[-2.966353,-9.600279,10.09775],[-5.515339,-3.9374,6.84997],[-1.828264,5.825093,6.186619],[-2.609492,-2.340165,3.644972],[-9.466926,1.075727,9.580181],[-9.306424,-2.141585,9.601871],[3.457961,3.91825,5.320731],[-3.643557,1.202157,3.964932],[0.7261388,1.692997,2.096071],[3.382061,-0.2094166,3.533015],[-3.274327,9.963273,10.53508],[-8.157183,1.063578,8.286786],[1.092556,-0.3908903,1.531821],[-5.184419,9.936253,11.25199],[-8.694625,-7.847409,11.75493],[-8.775327,9.150921,12.71793],[4.507646,-8.98075,10.09816],[6.03896,-5.937541,8.527803],[9.882386,9.785092,13.94308],[-4.427474,-2.490897,5.177557],[3.608245,-7.157462,8.077666],[-6.965281,-8.257265,10.84885],[-4.318688,5.981488,7.445084],[9.122805,-9.001202,12.85485],[4.132092,9.289428,10.21605],[8.161431,-7.645637,11.22785],[0.6748958,3.048122,3.27819],[3.653655,9.284356,10.02738],[-3.790133,-0.0155499,3.919865],[-9.668244,-1.151604,9.787807],[-2.765187,1.280298,3.207089],[1.184197,-6.835649,7.009166],[-9.183529,4.109454,10.11063],[-6.500918,-2.794175,7.146282],[-9.208233,-4.284757,10.20543],[9.749644,5.378233,11.17949],[-7.883329,5.65545,9.753512],[1.147047,-2.077212,2.57498],[4.421711,9.180385,10.2387],[5.234329,2.160676,5.750367],[-2.018193,8.457796,8.752566],[6.945543,9.084182,11.4788],[8.882716,-5.39202,10.43918],[6.943231,6.423381,9.511481],[-9.051401,0.7070491,9.133881],[-2.885267,2.525401,3.962627],[9.070497,1.61535,9.267323],[-2.15434,-9.201559,9.50315],[-5.863227,-1.333914,6.095634],[-0.4234815,-6.863136,6.948523],[-0.2275611,9.693198,9.747301],[-8.205126,-1.953052,8.493439],[-1.252419,9.973078,10.10103],[1.767045,-7.006361,7.294624],[4.263052,-5.186948,6.788081],[-9.595199,-7.089537,11.97202],[9.02706,-6.001201,10.88587],[7.545714,-7.160617,10.45047],[9.676149,-3.436427,10.31683],[-3.7532,7.045309,8.045053],[-9.46703,-8.41836,12.70801],[7.809815,0.8347263,7.917701],[4.881371,9.219986,10.48026],[5.969259,-7.705784,9.798529],[-0.7320172,7.397609,7.500698],[-4.716808,5.049247,6.981631],[6.230476,8.299274,10.42578],[1.701212,9.886738,10.08175],[2.878656,-6.306434,7.004126],[-8.797114,-6.933459,11.24554],[-1.881626,1.459852,2.582961],[-5.947207,-9.258037,11.04901],[-0.6368671,-1.630587,2.01604],[-0.2918648,-4.303746,4.428026],[7.435112,-3.680739,8.356359],[9.44147,-0.8009455,9.528005],[9.795087,6.137315,11.60217],[4.500471,4.0172,6.114911],[7.209281,3.245232,7.969019],[8.288128,-1.933182,8.569145],[-3.874629,0.7026983,4.062824],[6.617379,-2.699894,7.216587],[1.744847,-7.515306,7.779738],[-1.306021,-9.668832,9.807752],[-5.989498,7.852294,9.926359],[1.413848,-0.8467727,1.92769],[8.4804,3.368994,9.179723],[2.230171,6.892572,7.313085],[3.368983,-3.628169,5.051105],[-5.73092,-0.1397433,5.81919],[7.632129,-6.910865,10.34454],[-5.881357,6.162466,8.577083],[2.290205,9.498583,9.821818],[5.122074,-2.182298,5.656683],[-6.887595,-3.661613,7.864246],[8.158918,2.599089,8.621091],[-0.2921467,3.274863,3.43658],[-8.342263,-1.435831,8.523788],[2.582807,-9.525176,9.919671],[0.3724781,-8.115988,8.185842],[8.532197,2.982601,9.09364],[4.313929,0.8972671,4.518303],[7.779123,-6.91765,10.45795],[7.348828,6.313142,9.739663],[-9.177999,-0.8653681,9.272784],[3.537408,9.666174,10.34157],[6.733428,-8.746271,11.08315],[-0.030704,6.798696,6.871914],[-8.979212,-9.036183,12.77806],[5.883031,-9.4387,11.16687],[7.464423,-4.614001,8.832135],[3.242755,-4.140954,5.35378],[-7.480505,0.8367196,7.59329],[3.230056,-7.111603,7.874526],[-1.377859,4.305697,4.630067],[-2.004286,-0.5537186,2.307329],[-0.449096,2.574614,2.798272],[-1.884219,6.071489,6.435313],[0.9480044,-1.367679,1.941457],[-6.833823,-8.861439,11.23504],[8.808703,4.424764,9.908168],[6.798666,-0.5070146,6.890494],[6.283653,8.427235,10.55948],[1.693074,2.464782,3.153039],[-5.16184,-7.663837,9.29403],[8.191246,-3.268884,8.875929],[-6.334036,-8.776405,10.86947],[3.187806,-2.477021,4.159055],[-4.520676,2.57581,5.298236],[5.036611,7.885295,9.409853],[8.897625,-8.78321,12.54243],[3.929451,9.781341,10.58845],[4.124517,-9.674155,10.56413],[6.647089,9.82032,11.90052],[1.340202,9.472386,9.618849],[-9.399371,-2.904678,9.888646],[3.324319,6.275362,7.17156],[8.722598,-1.178565,8.858484],[8.966121,7.368019,11.64813],[5.540499,7.787099,9.609165],[1.081987,8.181644,8.313243],[2.085647,-9.165516,9.452862],[5.566218,6.853775,8.885776],[-4.87812,1.951048,5.348144],[-9.299462,-7.041002,11.70708],[7.168901,0.9973215,7.306695],[3.478347,4.011228,5.40267],[4.529773,-3.796454,5.994323],[4.254022,-9.12751,10.11969],[7.323007,0.2468884,7.395092],[-2.78158,1.976785,3.555962],[1.723064,8.329702,8.564631],[0.7630793,9.267803,9.352778],[-3.652495,9.295551,10.03733],[6.827487,-0.6469295,6.930591],[-4.240455,-2.705044,5.128228],[-7.640895,-9.234953,12.02779],[-5.51957,-8.795327,10.43185],[-5.280343,6.504825,8.437699],[7.635977,-2.071634,7.974949],[-2.491942,-2.75112,3.844273],[8.111567,-9.814979,12.77229],[-5.8247,0.6120947,5.941531],[7.376803,6.38062,9.804567],[2.903178,-2.000878,3.664963],[6.332159,-5.13226,8.211963],[1.936882,9.684022,9.926319],[0.8890902,3.920027,4.142112],[-3.415543,-4.268127,5.557233],[-4.274508,0.1698053,4.393206],[7.980322,-2.901489,8.550097],[5.515811,5.653863,7.961805],[3.009576,1.676532,3.587243],[-1.606355,-0.6335289,1.995429],[6.499166,3.346721,7.378326],[-5.729475,-6.100859,8.42896],[-5.937152,-5.684358,8.280199],[-4.284411,-8.995752,10.01398],[9.715587,0.8149091,9.800853],[1.086373,7.67498,7.815722],[-6.512734,-1.965229,6.875887],[6.758175,6.393722,9.356956],[-4.335945,4.823946,6.56284],[9.374019,-1.639588,9.568724],[4.109079,-9.574764,10.46712],[-7.178827,6.059341,9.447284],[-9.793056,-9.076901,13.39007],[-1.74497,-3.668664,4.18378],[-2.533357,0.3574499,2.746938],[-3.862628,-8.986589,9.832532],[8.723919,3.519015,9.459928],[0.9209474,-4.988805,5.170718],[-8.251559,2.903401,8.804429],[6.070555,-9.018433,10.91713],[-6.83289,2.900652,7.490138],[4.858063,-1.209929,5.105361],[9.292111,-0.5787583,9.363669],[-1.966451,-9.110035,9.373348],[9.324001,2.572565,9.723945],[1.12903,3.87883,4.161734],[8.850939,-0.3142475,8.912792],[2.499379,2.07109,3.396514],[5.93015,-1.467184,6.190259],[-2.072576,-6.517341,6.911679],[5.818531,-5.515587,8.079418],[5.850018,8.622422,10.46751],[1.200441,4.627161,4.883818],[6.557717,-6.342932,9.178041],[-4.660152,3.547728,5.941666],[-6.089442,9.046528,10.95084],[-6.43328,5.21052,8.338861],[1.605479,-4.90378,5.255912],[-5.243912,-4.948011,7.278834],[8.230492,1.712941,8.466119],[8.918855,-5.560242,10.55757],[-5.216883,5.803392,7.867352],[7.493085,-1.429891,7.693563],[-4.751748,9.633662,10.78826],[-0.6516693,1.429783,1.862512],[-6.440693,5.422465,8.478541],[-9.978743,-9.708,13.95781],[4.931269,-5.124073,7.181472],[-9.18252,-9.579339,13.30723],[-3.565254,-1.814815,4.123662],[6.98989,-9.973835,12.22031],[-1.242719,-3.118166,3.502472],[-2.813835,7.157004,7.755022],[0.8432949,7.438224,7.552372],[-5.05622,-1.731807,5.437326],[-8.109459,2.957004,8.689488],[-1.176257,9.071175,9.201619],[9.295687,-5.831374,11.01884],[-4.378824,3.209398,5.520357],[-6.063727,-9.434914,11.25995],[-7.069687,0.6632364,7.170799],[9.567055,1.796476,9.785493],[-4.660811,-5.280585,7.113912],[-4.957922,-0.3099502,5.067254],[5.499045,-1.302452,5.738979],[-4.853376,9.963265,11.12753],[7.00768,-4.06804,8.164345],[-5.789689,6.803094,8.989026],[-2.349271,-9.695623,10.02618],[-4.73419,2.790305,5.585549],[-3.152541,-0.9710784,3.446957],[0.7728165,9.762372,9.843839],[9.663936,-7.053321,12.00587],[-6.933073,-7.866616,10.53333],[-7.339271,-4.277293,8.55337],[2.748998,-1.52041,3.296762],[-9.944197,6.868835,12.12716],[-2.447401,-8.547337,8.946884],[-9.80699,1.282377,9.940902],[-1.105541,5.891223,6.076901],[0.5164777,8.678631,8.751307],[-0.2445812,0.08158284,1.032703],[-4.254823,-0.9485197,4.472494],[6.722893,1.411044,6.941782],[-5.755992,0.9357589,5.916679],[4.457723,6.380291,7.847255],[6.221089,-6.368893,8.959059],[-9.543841,0.8297982,9.631899],[-2.99072,-2.634945,4.109421],[8.059358,-5.888727,10.03147],[-1.162493,-1.152232,1.91808],[9.434093,9.092658,13.14072],[-7.188223,6.706398,9.881616],[-5.689917,-8.926414,10.63278],[-0.07317816,-1.740349,2.008524],[0.3561026,6.918749,6.999707],[-1.804646,-9.823473,10.0378],[-2.691843,-7.740145,8.255656],[6.235535,2.593819,6.827137],[9.541096,6.586262,11.63664],[-9.845352,2.717472,10.26234],[-1.315269,-9.03635,9.186162],[-0.4380977,5.066693,5.182983],[-7.521572,-5.111415,9.148804],[7.484688,-6.974573,10.27936],[-5.470227,1.863917,5.864944],[-3.510298,-6.482588,7.439498],[-3.051395,5.632531,6.48355],[1.407881,3.808188,4.181438],[-7.813371,0.6502113,7.903894],[0.8627448,-2.125735,2.502614],[-6.101016,-6.676701,9.09949],[-3.417766,-9.157271,9.825311],[5.280532,-0.6020511,5.408002],[1.114885,-5.95552,6.140944],[9.527917,6.092309,11.3533],[-1.380333,-6.47958,6.700021],[-9.568471,-4.109879,10.46168],[-9.850197,-4.939722,11.06468],[7.892043,-5.512895,9.678655],[-6.937071,-6.424359,9.507647],[5.54477,-6.150361,8.340948],[-6.803483,-6.827603,9.690384],[-0.7230152,6.90731,7.016671],[0.2352304,5.563879,5.657923],[0.003069225,2.103029,2.328677],[0.9044126,-7.45222,7.573213],[7.173554,-0.9856173,7.309673],[-1.981353,-5.148424,5.606428],[8.311226,3.043813,8.907372],[4.96744,4.685901,6.901676],[-4.257229,-2.582718,5.078822],[1.608892,-8.236609,8.451642],[4.481431,-0.2090146,4.596402],[4.828547,-4.154133,6.447611],[8.156691,-7.79523,11.32684],[3.285307,-8.333434,9.013288],[9.720474,-2.42429,10.06801],[9.464022,6.864293,11.73398],[0.004337854,2.725754,2.903404],[9.590227,-1.501706,9.758462],[-7.627338,9.333088,12.09474],[-7.167334,3.024849,7.843494],[4.71315,8.752306,9.990829],[-5.644831,4.753794,7.447327],[9.825429,3.632952,10.52318],[-5.714248,0.1671683,5.803496],[-6.193954,7.727718,9.954029],[9.838564,1.560915,10.01168],[-6.988922,4.453662,8.347463],[4.384174,6.645096,8.023608],[7.535374,0.3730427,7.610586],[3.450114,-0.2299299,3.599466],[-1.216119,-4.362785,4.638194],[9.992735,3.3671,10.59208],[-5.024012,-8.669656,10.06994],[-3.303157,5.487618,6.482653],[4.364093,-3.59158,5.739752],[1.986302,-2.019601,3.004028],[-1.288948,-2.437577,2.933116],[-0.1256664,8.982393,9.03876],[0.405513,4.006213,4.148997],[-5.733899,7.320365,9.352291],[-9.080676,-8.762447,12.65856],[8.198731,4.237147,9.28292],[8.037115,-8.275572,11.5793],[0.5744638,-0.9545674,1.497066],[-5.852849,-5.086843,7.818684],[1.594472,6.186516,6.466477],[9.09143,-7.980768,12.13865],[-8.069131,-6.914532,10.67341],[-5.143257,2.231025,5.694784],[6.691443,-8.106286,10.55875],[7.252676,-5.586685,9.209362],[-2.347368,-7.647577,8.061983],[-1.342586,9.198783,9.349874],[-4.920785,-1.009701,5.121877],[-2.631611,-8.55188,9.003334],[4.45489,-4.313223,6.280919],[-4.331608,-9.114691,10.14103],[6.099317,-2.854352,6.80801],[2.168861,6.071543,6.524384],[6.930819,6.399509,9.486304],[9.592655,2.548592,9.975688],[-0.08237074,-0.7589402,1.258084],[-5.567073,4.523337,7.242436],[6.7435,4.892031,8.390874],[5.836103,-8.169882,10.08995],[-1.086303,-1.43025,2.055644],[-6.089606,3.06208,6.889095],[5.122951,-8.950181,10.361],[8.049097,-1.384505,8.228293],[-3.731539,-4.640787,6.038318],[7.503011,4.778867,8.951689],[1.891532,8.891033,9.144855],[-5.095928,1.989572,5.561194],[0.6509949,-9.724121,9.797056],[2.127204,4.254426,4.86057],[-7.783567,8.854142,11.8313],[-3.560956,-5.471616,6.604467],[-0.2545398,-9.938697,9.992121],[6.908202,-9.386826,11.69768],[2.913597,-8.833198,9.354915],[4.263649,9.976884,10.89573],[8.163691,-0.7209648,8.256248],[-6.713606,0.279117,6.793409],[-9.228199,6.881586,11.55491],[7.96655,-0.6702566,8.056994],[-4.64616,-5.617944,7.358539],[1.387929,1.496449,2.272819],[5.898706,6.857209,9.100332],[-7.837719,-6.674828,10.34327],[6.013824,6.215673,8.706358],[5.759248,3.063235,6.599421],[-2.80477,1.819968,3.489845],[-0.8029559,2.512514,2.820898],[-7.662939,3.009042,8.293067],[-7.164707,-4.353498,8.443102],[6.428756,5.533668,8.5411],[1.057776,-4.19933,4.444465],[-3.734309,1.056851,4.007742],[-6.440798,-2.222136,6.886346],[0.4504617,1.132124,1.576268],[-3.654546,-2.556334,4.570619],[-3.101261,-7.091191,7.804025],[-5.539345,6.677545,8.733496],[-2.465519,-4.193901,4.966648],[-4.666554,6.308002,7.909969],[1.445589,5.005058,5.304746],[-6.907197,-3.959167,8.023987],[3.030773,-4.559714,5.565661],[-8.144963,-1.143366,8.285391],[3.364002,6.566556,7.445546],[5.151148,6.073227,8.026108],[-6.973269,9.977259,12.21361],[-7.518181,2.657354,8.036453],[7.625678,5.104453,9.230732],[5.634518,-4.961909,7.574189],[4.941858,3.901241,6.37508],[-3.894489,7.218748,8.263012],[-3.975139,3.162305,5.177056],[-4.909255,-3.408205,6.059426],[3.937864,-7.748788,8.749313],[-6.51581,8.319539,10.61464],[-7.141914,6.813171,9.921],[-2.444913,1.387689,2.983836],[4.102855,2.034686,4.687575],[1.753845,-9.889277,10.09326],[-9.553993,-4.765972,10.72349],[-5.231761,-5.543099,7.687475],[0.07842987,5.81964,5.905452],[-1.718462,-1.412537,2.438929],[-7.203867,-0.896296,7.327964],[-4.128525,0.6044146,4.290691],[8.974292,-0.6554713,9.053594],[-2.911673,8.436675,8.980831],[-1.420046,8.557204,8.731683],[-9.50174,9.150206,13.22911],[-3.609398,0.5118783,3.780182],[3.559828,-5.260384,6.429931],[4.840874,-4.946115,6.992719],[-7.372895,-4.928119,8.924458],[-7.850753,-2.263053,8.231387],[1.825685,-4.411115,4.877608],[-3.949073,-7.594954,8.618498],[-2.181359,-2.996744,3.839115],[8.140444,-0.6305048,8.225835],[-5.630362,-1.702116,5.966421],[-8.38429,-2.673295,8.856795],[6.423937,1.807007,6.747758],[8.359583,7.118059,11.02494],[0.9696544,6.827619,6.968257],[-4.896592,-7.638462,9.128127],[9.083185,-3.216083,9.687489],[2.459254,-2.957395,3.974181],[-2.39681,-1.270382,2.891119],[4.541985,0.4873715,4.676233],[-9.460445,-0.9671503,9.562186],[7.590107,3.009293,8.225908],[2.638251,-3.191267,4.259642],[-9.139729,7.036569,11.57791],[6.935173,-7.476312,10.24655],[-8.131156,-2.484712,8.560928],[-8.20285,6.448049,10.48161],[-2.566108,0.5511991,2.808688],[8.33815,-9.181812,12.44309],[-7.653275,-0.697346,7.749769],[5.750815,-6.578878,8.795085],[-1.644157,2.502182,3.156607],[8.001544,-4.54412,9.256011],[9.592201,-3.788678,10.36168],[-5.459747,7.904014,9.658275],[-8.607169,8.760939,12.32223],[-6.242472,7.63276,9.910977],[3.471325,0.4339884,3.638467],[3.026366,-6.124321,6.904071],[7.894794,-1.07993,8.030817],[7.002524,-1.336339,7.198691],[-6.571298,8.696166,10.94556],[-4.171717,0.5407022,4.323839],[-7.931824,4.270163,9.06356],[2.265677,-1.447201,2.868394],[8.949616,3.969532,9.841383],[9.591448,-2.640757,9.998473],[1.899026,-1.049223,2.388968],[-4.657338,-7.414188,8.812547],[3.039388,0.2759309,3.211544],[0.7200609,-2.477669,2.767188],[5.521271,-1.360895,5.773774],[6.904535,-4.97416,8.568248],[-1.515025,1.422287,2.306122],[8.901109,-0.983125,9.010898],[1.096773,-4.102502,4.362732],[5.31308,-7.815383,9.503107],[-9.773732,3.479935,10.42285],[-1.573908,8.146759,8.357444],[-7.376445,7.544026,10.59831],[-5.656526,-9.459756,11.06722],[-6.659563,-7.868097,10.35648],[-8.211194,9.664722,12.72126],[6.065735,6.730831,9.115767],[-1.80612,6.07508,6.416281],[-1.93766,-3.934766,4.498546],[-1.043132,-8.535367,8.656824],[1.082779,-1.33746,1.990279],[-8.260476,3.73684,9.121372],[3.133033,-8.634405,9.239526],[7.918419,-5.564094,9.729363],[-0.2925778,-8.412028,8.476309],[9.66062,-1.66644,9.854167],[9.236466,-4.083993,10.14846],[2.641332,-4.880038,5.638387],[3.318348,-4.692002,5.833208],[-7.413627,-7.310647,10.4598],[4.094108,5.166577,6.667476],[2.441275,-4.066107,4.846963],[7.721624,1.030003,7.853941],[-4.33989,-8.356751,9.469421],[-6.351804,1.380281,6.576519],[-8.924212,4.922951,10.24095],[9.967822,8.887691,13.39211],[0.7842956,-8.444777,8.53987],[-5.610432,4.5746,7.307798],[0.632697,7.831387,7.920286],[9.380954,-6.912312,11.6954],[7.028822,1.140217,7.19058],[6.19697,8.67802,10.7103],[6.169463,7.095762,9.455798],[-0.1794305,-1.602341,1.897285],[1.632371,7.942024,8.169479],[4.584695,-7.443453,8.799115],[2.813052,-5.689362,6.425115],[3.173967,-3.632269,4.926199],[5.85272,6.894988,9.099187],[4.54765,2.13636,5.123003],[-7.795033,0.9679657,7.918301],[-4.262423,-6.057808,7.474309],[1.476844,-5.96887,6.229645],[6.600505,-8.905945,11.13025],[6.422215,0.5131214,6.519826],[-9.202706,5.420678,10.72723],[-5.424316,-4.608298,7.187462],[-8.290477,2.834357,8.81848],[-3.508288,6.525858,7.47629],[0.7399307,-1.676838,2.08789],[1.987607,5.903499,6.308873],[2.861574,4.889699,5.753066],[5.169264,-9.835417,11.15602],[-1.751613,-1.887963,2.762707],[-7.267958,4.759838,8.745243],[0.7858691,-6.917046,7.033002],[-1.523984,8.703543,8.892366],[7.512656,-1.580856,7.742034],[-2.081262,7.739766,8.076858],[-1.503352,-0.1419066,1.811133],[3.251731,-7.453951,8.193604],[-9.923244,-7.380457,12.40733],[6.246166,4.960128,8.038498],[4.892059,7.104119,8.68336],[-8.77663,-6.096045,10.73271],[-5.245758,-9.159153,10.60227],[7.422558,1.000876,7.556197],[3.213123,-8.24697,8.907114],[3.939478,3.001538,5.052595],[-9.516781,-9.371164,13.39357],[-2.523352,3.426829,4.371551],[-5.729075,-8.570105,10.35707],[5.348776,-5.536144,7.762621],[5.292478,2.401353,5.897187],[6.082316,8.890553,10.81834],[5.125439,3.887425,6.510161],[-3.77132,-7.313207,8.288899],[7.840828,6.23939,10.07018],[5.73671,8.823674,10.57199],[3.344836,8.766995,9.436532],[-8.570426,4.464202,9.715004],[-2.904012,5.164235,6.008545],[2.686827,6.999568,7.563927],[-8.936654,7.465593,11.68755],[8.849547,-3.689308,9.639786],[1.130079,-1.693978,2.268621],[7.602661,3.283154,8.341435],[-3.090892,4.873279,5.856831],[5.706783,1.888113,6.093631],[-6.285462,-1.66278,6.578135],[5.267685,-0.6284767,5.398471],[-0.4415352,-8.816381,8.883891],[-2.86185,6.643712,7.302678],[-9.146807,-6.444479,11.23367],[-9.585285,5.416756,11.05527],[-2.635964,6.69594,7.265254],[3.881209,-6.300143,7.466966],[-1.568742,5.421486,5.731794],[-4.357574,1.633952,4.760068],[-1.98475,2.256488,3.167171],[-6.210318,8.472669,10.55245],[4.095752,-1.295031,4.410475],[8.676723,-5.959231,10.57345],[8.079821,8.639658,11.87128],[-9.415824,9.416656,13.35407],[-0.8661263,5.349517,5.510672],[-7.702868,9.904681,12.58717],[0.7006514,-5.067911,5.21293],[-1.85542,4.551318,5.015684],[-2.379559,-1.351352,2.913495],[-6.887371,3.731286,7.896731],[5.155287,-2.94114,6.018911],[6.837196,3.92751,7.948118],[-8.317698,-4.149278,9.348829],[-7.89674,7.05005,10.63305],[6.95813,3.24177,7.741101],[0.7217754,9.695105,9.77323],[-4.6004,-7.038717,8.468011],[8.907845,-1.22187,9.046694],[-1.057764,-6.452589,6.614738],[-4.763662,0.9064254,4.95117],[5.253841,-2.097028,5.744596],[2.770725,-4.220033,5.146416],[7.682407,-7.484457,10.77202],[-8.700435,-6.033474,10.63487],[-0.9729668,-4.512551,4.723323],[-1.719562,-3.481554,4.009752],[-6.162062,8.920661,10.88803],[5.081709,6.195903,8.075456],[0.920469,-5.596223,5.758904],[3.355679,-5.736771,6.720946],[1.895612,-4.447153,4.93665],[7.032446,2.618459,7.570444],[-6.687431,-9.483097,11.64692],[-7.534453,4.344533,8.754596],[4.341943,-6.434105,7.826249],[0.9652938,8.653297,8.764209],[6.959234,7.108282,9.997931],[5.868962,3.708798,7.014264],[-6.015137,-8.088872,10.12974],[-3.977848,9.711598,10.54222],[-2.339062,2.581921,3.624573],[9.550542,-0.241664,9.605793],[2.983241,-8.792088,9.338122],[-8.929605,4.373181,9.993125],[0.06810029,-7.89494,7.958311],[-1.9054,8.059601,8.341926],[8.743557,2.420268,9.127294],[0.3870269,2.071447,2.332527],[7.717556,-0.983084,7.843923],[3.978788,8.016177,9.00499],[8.449105,8.914165,12.32273],[-7.755002,5.760336,9.711927],[-6.982034,-3.896673,8.058093],[6.329274,1.272274,6.53287],[9.046917,0.8706508,9.143563],[9.730608,-1.583254,9.909158],[5.149648,-3.549349,6.333778],[4.629429,-1.694898,5.030337],[1.033128,-4.118344,4.362122],[7.756213,1.936475,8.056599],[-6.316041,9.288143,11.27661],[9.484395,-2.358924,9.824371],[-7.664793,9.95261,12.60173],[-9.916412,-7.367805,12.39434],[-6.598246,0.04601313,6.673752],[9.580043,3.652634,10.30141],[5.942468,-2.85431,6.667834],[-7.307456,-8.399507,11.17813],[-4.165298,3.4255,5.484866],[-6.38595,-2.287893,6.856735],[1.641474,5.207185,5.550604],[-1.334195,-1.697057,2.379092],[3.236351,7.901875,8.597302],[7.572099,3.136566,8.256799],[-0.3554333,-5.286224,5.391706],[-0.2011173,-7.753931,7.820735],[2.305842,-0.4710882,2.557114],[-6.802968,-3.015779,7.508348],[2.831103,6.660454,7.305942],[3.834315,5.54005,6.811323],[6.021662,6.106284,8.634067],[9.554968,-8.989521,13.15709],[4.209889,-1.554219,4.597691],[4.471103,-5.756631,7.357279],[9.566424,-3.851477,10.361],[-6.714189,-5.285353,8.603213],[-1.960425,1.447983,2.634373],[-6.248015,5.433964,8.340603],[-9.72469,8.866408,13.19783],[-0.1106547,1.158469,1.534371],[5.554611,3.68553,6.740685],[3.05864,-2.891586,4.326262],[-6.538094,-5.412646,8.546544],[6.865629,-5.912524,9.115635],[1.420902,2.510075,3.052776],[-6.286458,2.637825,6.890405],[0.233871,-6.101477,6.187303],[1.083033,-0.9272058,1.741457],[-3.974171,-6.740233,7.888268],[-6.318327,-5.749727,8.601198],[-3.191135,-2.517413,4.185775],[0.4266985,5.520627,5.626667],[-1.521813,3.056114,3.557492],[-0.6417565,-0.8126413,1.439527],[4.367552,1.771806,4.818175],[2.593255,5.75035,6.386823],[1.868464,7.535292,7.827629],[3.450119,3.188926,4.803392],[0.8317133,-2.112475,2.480786],[6.534982,0.0442455,6.611198],[6.240709,-5.659387,8.483815],[-5.97239,8.756268,10.6462],[2.328843,5.210394,5.79411],[2.712455,3.620313,4.632934],[2.716785,-2.250359,3.666747],[-9.760301,0.2789031,9.815358],[9.66917,-6.137235,11.49602],[8.327071,0.3147519,8.392805],[-6.091465,-6.372541,8.87216],[4.106571,4.945037,6.505176],[5.536244,3.701572,6.734363],[-5.290969,-2.929846,6.130118],[5.180579,2.297593,5.754766],[0.6832566,0.08816669,1.214336],[-2.690414,9.533264,9.955976],[3.965437,6.270551,7.486288],[-7.413686,9.942178,12.44225],[5.37561,0.71338,5.514172],[-4.761615,5.871332,7.625321],[-1.961591,6.754912,7.104694],[-9.554995,0.8501869,9.644726],[4.492762,8.114991,9.329415],[-5.103177,-4.530147,6.896712],[7.859592,1.502844,8.064225],[-6.565852,-1.904069,6.909117],[5.588115,-5.661044,8.017135],[9.541453,2.329136,9.872396],[-7.996684,-7.369851,10.9207],[0.05165851,3.219379,3.371509],[3.851957,7.827186,8.780798],[-8.36538,0.9758527,8.481266],[-3.565815,4.11937,5.539336],[-9.905754,-1.005444,10.00674],[-9.818669,-7.35217,12.30694],[7.863722,-1.46149,8.06065],[3.921498,6.455255,7.618954],[-4.351451,4.40727,6.273688],[-4.320205,-1.097396,4.5682],[-1.466658,2.138135,2.778975],[8.439494,-2.970955,9.002869],[0.4265202,5.579271,5.684205],[-5.825909,-5.46436,8.049872],[-9.456888,-5.921628,11.20261],[5.461783,-3.96217,6.82128],[-4.82569,-4.642879,6.770791],[-3.788849,-7.504992,8.466421],[5.978655,-6.893345,9.179461],[0.116692,-0.7462668,1.253208],[5.378054,-7.652009,9.406206],[-2.360428,8.898854,9.260736],[9.965649,-8.621993,13.21563],[-4.170848,5.957799,7.341073],[-4.10597,7.112105,8.272909],[6.937434,5.077609,8.655062],[2.562068,-5.524598,6.171335],[-8.634955,-9.516273,12.88883],[2.099887,-1.132987,2.58712],[-3.012925,-9.312162,9.838399],[-1.709529,4.106616,4.559253],[-3.921282,-3.445087,5.31461],[2.085742,8.655048,8.958805],[1.359517,-4.407688,4.719746],[5.342777,9.518728,10.96136],[-4.940499,-5.061295,7.143195],[-0.247977,-1.805706,2.078958],[9.782456,5.382183,11.21001],[-2.331711,3.114439,4.017039],[8.966556,-5.507945,10.57055],[-3.464669,-4.609225,5.852255],[-7.351461,-6.866153,10.10881],[0.5360526,8.934681,9.006435],[4.384489,5.792258,7.333076],[6.544532,7.994675,10.38006],[7.459399,6.353863,9.849579],[4.98022,-4.943608,7.088149],[9.291154,-1.486595,9.46232],[2.139353,9.518391,9.806966],[-8.80283,-3.522274,9.533951],[8.395709,-0.4741697,8.468339],[9.731677,6.89488,11.9685],[-4.499157,-2.006398,5.026733],[-0.8001164,6.568672,6.692357],[6.836718,0.325233,6.917116],[1.652053,-1.482573,2.434605],[6.057894,-1.343887,6.28523],[-4.145637,-2.346744,4.867598],[-8.060881,-0.03882799,8.122765],[-5.71077,-4.473901,7.323161],[7.48904,-8.861048,11.64491],[5.740404,4.939793,7.638965],[-4.74698,1.924257,5.218868],[-2.823775,-7.717319,8.278328],[-8.628425,-5.341908,10.19734],[6.43164,-2.607637,7.01183],[8.430803,4.737818,9.722416],[-4.659626,7.431812,8.828587],[6.685815,-5.611357,8.785639],[0.1051155,-5.857368,5.943047],[-9.359874,0.5586219,9.429703],[4.074767,4.524537,6.170507],[-4.039064,-1.87703,4.564787],[3.148256,-2.764809,4.307631],[8.86189,-0.222809,8.920916],[-9.04057,-3.216115,9.647554],[-9.841917,0.6955872,9.917014],[4.467254,5.138409,6.881832],[-8.068205,-0.02004646,8.129965],[6.537325,8.001423,10.38072],[-6.597108,-6.702083,9.457259],[-6.373348,6.440958,9.116222],[-8.749536,-5.958732,10.63301],[5.149532,-0.5191855,5.27136],[-7.491467,-7.558377,10.68883],[0.6213144,5.101829,5.235904],[-4.063463,-7.085945,8.229359],[-2.493814,6.42942,6.968253],[7.437223,-3.112956,8.12421],[9.930758,-1.441631,10.08455],[1.246011,8.332426,8.484213],[2.162602,1.983187,3.099981],[-2.310179,2.136279,3.301608],[-2.569721,-7.256524,7.762771],[-8.60691,-2.811168,9.109422],[-5.856011,-3.09377,6.698081],[0.7385638,-3.430947,3.649229],[2.95355,3.311035,4.548232],[-9.52874,4.605996,10.63071],[-6.832178,3.556505,7.76707],[-2.355436,-5.667684,6.218579],[-8.485926,1.503675,8.675942],[2.195603,3.237589,4.037654],[-2.443788,-4.422638,5.150906],[5.173549,6.928321,8.704438],[-7.776254,-8.034103,11.22573],[8.296927,-0.5362085,8.374158],[8.904134,9.059843,12.74223],[-9.250561,9.672396,13.42118],[5.405952,0.1290019,5.499178],[-5.17282,2.242718,5.726068],[3.919319,-0.7608084,4.11581],[-2.278648,7.727407,8.118193],[9.304192,-2.816639,9.772483],[6.231879,-2.345387,6.733287],[-1.653624,-1.864072,2.685002],[5.264566,-3.480803,6.389964],[-9.706495,-7.054489,12.04084],[-2.679351,-4.590302,5.408308],[6.641344,-6.419457,9.290689],[-2.913373,-6.908608,7.564166],[3.870156,4.831626,6.270783],[-7.393209,4.330655,8.626361],[1.183056,-2.271774,2.749651],[9.08559,2.433326,9.458807],[5.684877,-6.011395,8.333948],[9.011618,-4.404053,10.07993],[-8.631409,-5.459059,10.2617],[7.876629,7.705898,11.06446],[8.418513,-9.029033,12.38527],[3.622995,6.075939,7.144447],[-6.705878,2.394955,7.190592],[4.718052,7.856261,9.218506],[6.456702,-0.2307421,6.537755],[2.685582,5.483787,6.187428],[6.417583,9.863062,11.80955],[-8.336783,-5.668727,10.13096],[-8.916729,6.150015,10.87799],[-0.3382897,-7.915481,7.985566],[5.394983,9.396948,10.88157],[-9.017365,8.893993,12.70496],[-5.133132,-7.645036,9.262592],[0.9121807,6.035666,6.185575],[2.698978,7.885425,8.394308],[9.213474,-3.386525,9.866948],[8.962705,7.412101,11.67344],[4.300461,5.755917,7.254278],[-9.070176,-1.521765,9.251155],[-4.831354,8.565942,9.885208],[8.493099,-4.161547,9.510584],[-5.254676,2.810923,6.042592],[-2.929318,-1.700692,3.531749],[-5.206558,2.557817,5.886482],[0.7493373,-7.718967,7.819461],[9.81639,-3.052752,10.32864],[2.225068,4.864552,5.441947],[0.1344084,-8.253077,8.314526],[2.339441,-0.8816288,2.692629],[-9.689688,-9.611759,13.68488],[3.323415,-8.120241,8.830821],[-7.273092,8.890711,11.53007],[-1.398099,4.331806,4.660388],[8.747486,-7.951686,11.86372],[-5.067562,2.287933,5.649321],[-5.519969,-0.3342389,5.619766],[-1.812067,-2.803627,3.484811],[-7.362425,-7.271195,10.39594],[-1.290146,-6.171833,6.384042],[-7.096188,-8.086367,10.80487],[-5.429584,2.564573,6.087481],[5.832578,-4.313067,7.322671],[-5.942492,4.803075,7.70602],[6.093147,-4.735665,7.781579],[-7.090342,-7.195148,10.15101],[8.974855,-9.789474,13.31848],[8.745869,6.83419,11.14434],[-1.178647,4.644712,4.895156],[-0.167893,8.34449,8.405873],[1.904197,5.89456,6.274696],[-6.109274,0.5557809,6.215475],[8.755855,4.906089,10.08636],[-9.880948,-1.798324,10.09292],[6.93046,-3.91654,8.023128],[-0.6973412,8.06389,8.155527],[5.670933,-4.990992,7.620334],[-2.2499,4.243398,4.905963],[9.680354,1.064406,9.789904],[-4.644562,5.451469,7.231215],[-1.301411,3.209804,3.605067],[7.548035,-8.489804,11.40393],[5.480314,1.018457,5.663135],[-2.775547,9.785541,10.22059],[4.002745,-9.348926,10.21882],[-0.6218412,0.227203,1.199295],[5.158777,-8.978812,10.40346],[6.243384,5.139518,8.148281],[-5.589895,-5.539083,7.93274],[-8.571853,-1.484201,8.756684],[8.012198,6.763923,10.53309],[-3.616441,-2.206655,4.352927],[-1.090332,-5.183389,5.390394],[-1.210217,1.705862,2.318317],[-3.317095,1.584062,3.809511],[1.391651,8.610571,8.779444],[9.417293,1.088147,9.532548],[-4.738892,7.317909,8.775471],[-3.088496,7.106491,7.812875],[-7.760248,-0.9514556,7.882051],[-8.694865,-4.343307,9.770618],[5.997252,9.162061,10.99592],[8.955977,9.724174,13.25779],[7.432267,-1.615341,7.67124],[-2.729711,7.706961,8.237024],[-4.181184,-9.063421,10.03135],[-4.154621,-5.545578,7.001022],[-1.364996,-8.247792,8.419579],[9.825677,4.709202,10.94169],[2.263019,-9.656537,9.968449],[-5.227311,4.568881,7.014232],[6.863929,-2.975228,7.547549],[5.31236,-2.96378,6.164832],[3.267597,9.184994,9.800067],[-1.843375,2.228842,3.060354],[-9.684142,6.058618,11.46689],[6.425085,8.608649,10.78845],[-1.680075,-0.01013124,1.955187],[7.469435,-7.385708,10.55183],[-0.278835,-2.467799,2.677271],[6.609545,4.503655,8.060334],[8.748904,-5.104685,10.17846],[0.1020295,-9.888182,9.939142],[-9.750089,-2.884565,10.2169],[-8.774037,9.647614,13.079],[-3.777452,-3.71931,5.394665],[-9.86069,1.210602,9.984927],[1.506357,-6.463502,6.711629],[-5.071323,8.040187,9.558395],[-5.555139,-8.524122,10.22351],[-5.0653,-1.573336,5.397467],[9.563263,0.9779875,9.665012],[0.7141263,7.939212,8.033746],[-5.553621,9.850652,11.35245],[-5.881377,-8.221259,10.15774],[9.005633,-0.3744994,9.06872],[-9.025229,-3.601507,9.768603],[6.65066,7.566576,10.12346],[0.8085372,5.029103,5.190916],[2.032837,-1.907108,2.961332],[8.789245,-0.4793704,8.858929],[-3.617323,-2.095429,4.298354],[-9.099251,-6.303992,11.11471],[0.6466826,0.6772971,1.370011],[-8.179686,-3.059276,8.790133],[-1.221609,0.6041986,1.69038],[-4.655384,4.092217,6.278443],[-6.971653,-0.1437718,7.044475],[7.96048,-4.929736,9.416556],[9.048763,6.648025,11.27281],[9.480151,-4.017835,10.34487],[-6.518622,-6.95547,9.584936],[-4.661863,7.60145,8.973016],[-7.492027,-8.984595,11.7411],[7.760184,-7.593283,10.90314],[-0.7040624,-3.416538,3.628834],[-0.1242139,-3.866629,3.995779],[8.908635,5.680697,10.61292],[3.989624,9.404256,10.26436],[2.744138,3.270545,4.384833],[2.820515,-8.217632,8.745558],[-1.45814,-2.891018,3.388828],[-2.992555,5.971439,6.753775],[-3.714221,8.616257,9.435853],[0.2560058,-1.414347,1.750976],[-4.498209,9.467051,10.52896],[-8.752057,-9.207722,12.74287],[1.334849,5.513638,5.760384],[8.461061,4.323266,9.554066],[-8.432881,9.241212,12.55044],[-3.217351,-9.346067,9.934804],[0.6229316,7.519641,7.611376],[1.628538,-1.296076,2.309101],[-9.465523,6.338847,11.43578],[7.065081,-6.036934,9.346654],[8.621999,-0.5058121,8.694521],[-4.247729,-9.543222,10.49363],[-0.9477508,0.03759045,1.378276],[8.890133,4.517643,10.02215],[2.0116,-2.082537,3.063249],[-7.505306,-9.225971,11.93517],[2.36611,2.101718,3.31899],[-3.33867,-2.447456,4.258727],[5.828794,-8.462157,10.3239],[-7.134597,5.976293,9.360478],[-8.185899,-5.263017,9.783061],[1.698678,6.850496,7.12845],[7.3473,-2.647092,7.873367],[8.634253,-5.799885,10.44935],[0.1925482,-4.227911,4.348828],[7.329391,-6.806287,10.05214],[-9.287465,-1.965493,9.545689],[-2.014692,-9.492306,9.755144],[-1.922591,2.764014,3.512283],[-1.29779,0.3498677,1.675311],[9.815413,9.297482,13.55675],[7.095695,4.883923,8.671885],[-5.638193,-1.126377,5.835919],[4.954532,6.950664,8.594132],[3.594721,7.875298,8.71449],[-0.1549462,9.639007,9.691979],[-0.8411581,-1.084186,1.697942],[2.096442,-3.98409,4.611729],[-2.665272,3.031994,4.158925],[-5.242934,-5.796334,7.879457],[-3.009352,1.360167,3.450544],[-2.700135,-1.152307,3.101377],[2.291648,-8.31062,8.678597],[3.056344,-1.038729,3.379378],[7.783367,9.582894,12.38599],[4.343237,-9.35814,10.36525],[-4.010896,8.03449,9.035502],[-5.968506,-0.4497551,6.068389],[-8.029069,3.757147,8.92088],[3.512609,2.072149,4.199074],[-1.011096,-0.5239977,1.515549],[-2.335566,9.126686,9.473715],[-5.474001,1.112651,5.67474],[1.702962,-3.95173,4.41772],[6.219316,-9.266407,11.20474],[6.044534,5.17817,8.021835],[-2.634613,7.970249,8.453759],[-6.672451,5.64056,8.794175],[4.519989,5.107052,6.892915],[-3.02192,-8.217034,8.812018],[1.33807,-2.131897,2.708398],[-4.765869,2.518663,5.482442],[3.194593,-3.837593,5.092401],[-8.492406,-8.415919,11.99786],[-0.1218517,-1.235991,1.594529],[-4.608742,1.736366,5.025482],[-2.874432,4.579149,5.49827],[-4.453241,7.003469,8.359422],[6.158055,5.206137,8.125608],[-5.731982,7.914561,9.823232],[5.960926,4.749156,7.686815],[-3.112435,-3.158636,4.545794],[0.2766018,-7.739667,7.808903],[-2.887481,8.668385,9.191216],[-6.764134,6.418889,9.378467],[-7.859723,-1.715155,8.106603],[4.126016,-7.35873,8.495583],[-3.542499,-4.532334,5.83878],[0.5339231,-0.8909119,1.441804],[-2.081234,1.796609,2.925634],[-2.842892,8.68109,9.189307],[-9.951938,4.897386,11.13667],[-0.3280921,0.220888,1.075377],[6.215904,7.040621,9.444988],[1.446348,5.096464,5.391277],[-3.708343,-7.751833,8.651169],[-0.002566022,1.479632,1.785866],[1.874298,8.665091,8.921703],[-6.506621,-7.233875,9.780852],[-1.275742,6.726372,6.91893],[-6.366201,-1.458828,6.607321],[-4.658998,4.66425,6.667945],[-9.433063,-5.543643,10.98702],[-2.501868,-0.2671383,2.707528],[-1.350676,-0.3601305,1.718727],[4.002132,7.692793,8.729039],[-4.401681,0.1279342,4.515658],[-0.2190793,-4.614385,4.726579],[3.480919,-4.216139,5.558114],[9.599639,-9.6509,13.64892],[-2.20225,-9.355509,9.663097],[3.807378,-5.496129,6.760441],[9.406713,-0.1573883,9.461027],[-8.592397,4.526348,9.763048],[-2.091454,-9.054244,9.346311],[-8.828965,-8.527289,12.31525],[8.764039,6.335141,10.86013],[-3.003377,8.59447,9.158886],[-7.862121,-1.402326,8.048569],[-6.490038,8.896433,11.05745],[-8.587867,6.218555,10.64997],[-5.865268,1.686724,6.184368],[5.694637,8.542978,10.31559],[-4.710813,-1.702457,5.107849],[3.718666,3.879225,5.465974],[-8.211199,-7.16979,10.94667],[-5.418053,-5.798284,7.998462],[5.477439,3.478985,6.565491],[-7.968633,-6.507578,10.33672],[-2.273689,-4.890684,5.485294],[-9.060187,-4.532425,10.17988],[-5.806609,3.062366,6.640391],[-6.570701,-2.363581,7.054121],[-0.6823145,-1.955002,2.299475],[8.201462,-9.440188,12.54516],[4.121156,-6.447497,7.717134],[-2.49813,-0.00123166,2.690846],[7.942626,-2.841238,8.494582],[1.612689,9.601096,9.786818],[3.402644,-1.511659,3.855269],[7.55603,8.795719,11.63865],[-5.553257,6.100834,8.310165],[3.398596,-5.177855,6.273805],[-9.544361,8.947153,13.12046],[4.647408,9.876268,10.9608],[5.288739,2.327195,5.864008],[0.719022,-6.080226,6.20372],[3.316679,-2.15897,4.081852],[-0.9522999,-0.89951,1.648027],[-5.447639,-5.487117,7.796489],[4.440183,-2.822706,5.355641],[-6.182523,-4.375196,7.63976],[-6.426242,4.640546,7.989447],[9.054606,-9.650908,13.27124],[-0.7563801,9.219437,9.304307],[-8.854386,-8.188512,12.10173],[-3.449108,-2.822445,4.567553],[5.776517,-9.365284,11.04883],[-6.185121,4.530274,7.731695],[-8.104848,8.766241,11.98063],[7.139179,1.205273,7.308937],[-6.296031,-8.536051,10.65383],[-3.497852,7.212914,8.078433],[8.575718,5.877668,10.44461],[0.7321961,3.73062,3.931111],[-3.17447,7.690309,8.379624],[-4.953062,3.430844,6.10766],[3.726999,1.441295,4.119205],[2.830591,-4.727957,5.60052],[1.253199,-8.023269,8.181891],[-3.300244,8.248303,8.940141],[2.598935,6.450672,7.026068],[6.118854,-0.8247897,6.25465],[2.535256,7.773664,8.237559],[1.750148,2.644239,3.324909],[-5.480698,7.616063,9.436232],[-2.477341,5.715315,6.308887],[-8.359075,-8.674232,12.08786],[2.058509,-1.386204,2.675634],[-0.2880091,-0.9446818,1.40548],[3.853998,-4.108811,5.721507],[3.345562,-2.07923,4.063986],[-3.488717,-7.060175,7.938338],[-9.906375,-3.078929,10.4219],[-2.838387,5.312733,6.105864],[9.849002,4.112174,10.71974],[-9.061871,-1.155482,9.189812],[9.111228,-0.08303156,9.166317],[-3.756477,1.371501,4.122152],[-1.507452,-1.886025,2.613331],[-0.58546,-1.205292,1.671973],[1.97936,-6.908501,7.255704],[2.201385,-1.234532,2.714805],[-8.938419,7.885831,11.96168],[-1.845297,4.45355,4.923336],[-4.036589,-6.875444,8.035284],[-3.328401,-6.001176,6.934866],[6.043948,1.617139,6.335964],[-5.187975,-0.5583677,5.312895],[9.902233,-8.881231,13.33906],[-5.784319,9.678064,11.31916],[2.685019,-8.34774,8.825764],[2.63303,-9.275112,9.693326],[6.277372,-1.153648,6.460364],[-2.016976,-5.513213,5.955142],[-2.421239,8.743479,9.127477],[-4.146562,2.954444,5.18871],[-9.037774,7.762835,11.95588],[-1.324682,8.197504,8.363843],[7.756094,4.20682,8.879997],[-4.411003,1.319431,4.711459],[-4.821663,-0.245895,4.930406],[3.719921,-2.870233,4.803753],[-7.615142,7.625823,10.82329],[4.9959,-2.83152,5.828938],[-2.855427,-4.817865,5.689049],[1.77477,3.524735,4.071065],[0.7634682,-1.379442,1.867015],[5.010292,-1.503087,5.325627],[0.5012153,0.01284967,1.118652],[8.257016,2.624504,8.721601],[1.435719,0.484439,1.815481],[-9.445323,6.765341,11.66122],[4.808167,0.07687366,4.911657],[6.356051,-8.039009,10.29685],[9.258588,6.240853,11.21025],[7.749983,-6.373483,10.08382],[-9.316336,-6.740025,11.54218],[-2.150416,1.224058,2.668822],[-0.6554601,-7.155839,7.255044],[-2.65005,-1.705379,3.306219],[2.650752,-3.973596,4.880158],[-1.006485,5.569822,5.747689],[-6.150459,0.3751894,6.242509],[-1.674111,-0.8397586,2.123168],[-3.943789,-9.676479,10.49703],[6.676079,1.451513,6.904847],[-8.09348,-2.54944,8.544241],[-5.580006,-2.607285,6.239744],[9.752867,1.389617,9.901992],[-8.06992,-1.837819,8.336737],[-7.41969,-9.952277,12.4539],[9.640874,8.730766,13.04503],[4.8515,-8.992544,10.26659],[-1.664575,-1.198587,2.281977],[0.1871968,1.574679,1.874742],[2.16701,3.788843,4.477864],[-0.8812708,-5.368802,5.531788],[-1.773844,9.317446,9.537364],[-8.699462,-5.168159,10.16811],[1.823223,-7.754296,8.028278],[-9.091809,-9.377187,13.09934],[-8.67146,-2.458014,9.06841],[-8.238868,7.921352,11.47287],[-5.036019,-1.109045,5.252758],[2.640137,-4.724244,5.503527],[-1.190008,-3.394661,3.73361],[-5.995418,7.623678,9.750154],[7.918768,4.151981,8.99699],[5.723397,-2.838668,6.466475],[-5.424822,-8.702673,10.30365],[-8.104388,-1.231211,8.258146],[4.144997,-1.315554,4.462251],[-0.763609,-7.416477,7.522448],[7.301703,-7.402964,10.44599],[-4.222275,-3.189601,5.385273],[-3.342825,6.512937,7.388695],[-6.236782,6.621729,9.151216],[-8.018916,-7.18769,10.81508],[1.040102,3.415309,3.70758],[-3.859879,-9.023204,9.864931],[8.19863,-2.979745,8.780457],[8.13558,1.004134,8.258084],[5.94966,5.234669,7.987504],[8.919057,-1.185498,9.052899],[7.272133,6.004961,9.483853],[-2.783334,-2.251608,3.71708],[-5.277263,-3.601145,6.466664],[3.63908,-2.274963,4.406627],[1.864789,-4.459811,4.936329],[-2.691373,6.117193,6.757481],[-5.012511,2.863964,5.858972],[1.545533,-1.426515,2.328866],[-6.161038,9.789037,11.60964],[-0.1473002,-9.674749,9.727408],[-9.10043,-2.526996,9.497554],[-6.200179,-1.494537,6.455685],[0.2768586,-6.791515,6.870323],[-7.889418,1.267584,8.05293],[-3.936018,-4.291722,5.908563],[6.210496,-2.032846,6.610804],[2.976514,2.645685,4.106006],[-3.909181,0.3753382,4.052478],[-0.3483818,-5.173019,5.280293],[3.75699,-5.004688,6.33734],[7.125709,-8.258592,10.95354],[7.008283,0.4115117,7.091218],[-6.948149,-1.918983,7.277312],[6.815001,-9.847961,12.01776],[-4.566868,-0.8811589,4.757387],[-2.187612,2.710229,3.62367],[-3.33573,2.843475,4.495825],[6.295224,9.358687,11.3232],[6.761526,2.545442,7.293662],[8.70917,7.545433,11.56647],[-7.607735,4.410984,8.850673],[8.811073,4.441425,9.917725],[-6.558119,-1.44778,6.790066],[-0.157008,-6.101405,6.184803],[-0.823268,-9.595239,9.682272],[3.031785,9.295607,9.828531],[5.298597,4.597368,7.085967],[0.9390346,3.085336,3.376549],[-7.215797,-2.252079,7.624932],[5.657963,-8.142364,9.965472],[-8.947058,-0.1616536,9.00422],[7.549945,-7.775246,10.88375],[-0.3071874,-6.258468,6.345296],[9.382704,-5.702458,11.02512],[-3.229882,-4.88838,5.943769],[7.024275,3.595139,7.953959],[-2.767973,2.992092,4.196939],[1.483588,8.040079,8.236741],[3.453118,4.525949,5.779986],[-8.9196,-3.160176,9.515566],[-5.792703,-6.936838,9.092587],[7.307197,-8.679286,11.38969],[9.564786,-6.247383,11.468],[-9.887691,0.3943029,9.945949],[8.665423,-0.02664633,8.722974],[-7.876329,2.243804,8.250528],[5.066976,-6.356473,8.190177],[-1.112815,8.163833,8.29979],[6.213939,-1.951916,6.589614],[-6.724212,0.2824129,6.804027],[-6.722588,-0.4982205,6.814793],[-1.460195,1.998938,2.669817],[3.238274,2.761164,4.37155],[-2.699486,-9.475277,9.902934],[8.953737,1.709341,9.170129],[3.679215,-5.06748,6.341607],[9.768892,2.79393,10.20967],[-1.92808,-3.840088,4.411776],[-4.570403,-1.164003,4.82115],[-7.929158,-8.693118,11.80855],[3.492252,-9.882214,10.52872],[-9.904985,0.45875,9.965901],[8.695321,-5.568851,10.37404],[3.841689,2.379962,4.628476],[5.798681,5.338654,7.945183],[1.129591,8.250124,8.386926],[-8.275185,-0.5251546,8.351913],[9.773652,-2.488147,10.13485],[-0.005044406,7.744898,7.809191],[1.490057,8.808585,8.989518],[-3.576774,6.918777,7.852566],[4.626541,-3.106527,5.661748],[-5.174949,8.846015,10.29719],[6.693248,-0.5275385,6.788068],[0.4719242,-4.972604,5.094065],[-7.243177,5.430484,9.107897],[-6.066568,9.1572,11.02985],[-7.072553,-6.128037,9.411367],[-4.507348,6.91926,8.318193],[6.87511,3.839505,7.937817],[5.224842,5.181442,7.426056],[-7.607229,2.33487,8.020072],[-6.814748,0.2450995,6.892087],[-9.334838,5.711603,10.98916],[-6.47277,-5.674035,8.66553],[8.364111,0.1291814,8.424669],[4.820153,-2.893068,5.709967],[8.507408,2.960811,9.063244],[-3.906602,9.89422,10.68443],[3.057073,-7.405908,8.074228],[-3.617354,-1.625133,4.089781],[9.408063,-6.921834,11.72277],[3.768478,-6.018914,7.171384],[2.13392,8.84993,9.158322],[3.323169,-6.733059,7.574796],[-3.851066,-5.083335,6.455308],[-0.6354638,9.982995,10.05306],[-2.691278,4.111516,5.014733],[-8.095251,-5.038843,9.58765],[1.146332,-6.596794,6.769916],[-7.868552,5.90678,9.889599],[2.39183,5.231067,5.838229],[4.609417,-2.658875,5.414457],[-8.647758,-9.548406,12.92114],[7.628523,3.240009,8.348175],[-2.48239,6.190083,6.74384],[6.336138,-9.344101,11.33397],[7.565905,4.233246,8.727159],[-3.508573,-7.707371,8.527229],[-4.697723,2.035825,5.216626],[8.59857,7.408839,11.39413],[4.80681,5.250265,7.188234],[-3.336297,8.888959,9.546961],[0.4867842,8.215775,8.290712],[-5.046217,8.070573,9.570708],[-6.626388,0.167133,6.703503],[4.465596,2.094603,5.032784],[7.138154,-3.55425,8.036537],[3.790457,0.7672266,3.994521],[7.991635,-0.7210651,8.086171],[7.747819,-1.794024,8.015436],[-2.753538,-8.350665,8.849609],[-3.225133,-5.777491,6.691853],[-1.229509,5.664421,5.881951],[6.108497,-6.421207,8.918837],[-8.343468,0.3801396,8.411775],[-9.227252,-9.006428,12.93282],[-0.397496,8.681831,8.748268],[4.404018,0.7060758,4.570987],[9.835362,-9.508087,13.71634],[-4.771787,5.149257,7.091177],[2.828335,-4.092783,5.07448],[1.008883,5.978785,6.145219],[4.263533,-2.58185,5.083666],[-8.485537,-4.598759,9.703243],[-9.48019,6.38921,11.47589],[6.712624,-1.365494,6.922708],[1.414343,-1.143638,2.075638],[-9.482386,-5.671178,11.09405],[-3.794306,-5.964787,7.139709],[-0.3386682,-2.933716,3.117914],[-3.829758,-7.527452,8.50468],[2.848239,-2.734752,4.073246],[-6.673223,7.120142,9.809604],[5.494308,-0.8725359,5.652322],[-0.8194822,7.531318,7.641485],[-3.192826,9.380002,9.958844],[8.104774,-0.05057115,8.16639],[-7.465996,0.3622955,7.541376],[-3.315508,-7.780128,8.516042],[8.902155,2.701534,9.356637],[5.451782,-7.083618,8.994419],[-9.95967,-7.174554,12.31541],[3.089015,-1.770146,3.69803],[-0.0772851,-5.052524,5.151114],[7.893746,-3.407162,8.655633],[1.33458,0.1459034,1.674034],[-8.996206,-9.785358,13.32985],[1.942913,-1.080992,2.43792],[6.074453,9.640243,11.43824],[0.9771609,8.476769,8.591302],[-5.720787,4.19653,7.165073],[3.699771,-3.064417,4.907032],[4.224403,-3.663934,5.680669],[-3.54878,-8.805574,9.546307],[-4.697133,8.407189,9.682142],[7.865429,-4.505397,9.119407],[-3.852514,7.852449,8.803568],[4.025621,0.1871599,4.152186],[7.304988,-3.252869,8.058784],[5.521527,1.889843,5.921044],[8.551996,-2.215318,8.890685],[2.73181,-4.204573,5.112848],[-5.109887,-8.29265,9.791781],[8.989524,2.96101,9.517306],[-6.209202,-4.516525,7.742945],[6.486707,6.613221,9.317299],[9.139816,-6.624868,11.33248],[-1.292964,5.224843,5.474554],[-0.1063254,-1.293281,1.638255],[2.259216,5.133118,5.696749],[1.303775,-0.7563099,1.808821],[2.842672,-8.042193,8.588227],[-7.415227,-0.3837342,7.492186],[5.003982,9.205214,10.52501],[-2.767699,2.86475,4.106939],[7.453221,5.480054,9.304918],[8.823992,8.389367,12.21656],[7.943764,-4.775422,9.322449],[6.796726,-5.774668,8.974535],[6.471035,4.971551,8.221352],[1.348951,-4.429848,4.737428],[4.810726,-7.823994,9.238937],[8.849734,9.221502,12.82006],[-1.6973,-8.576768,8.800101],[-4.549253,6.337482,7.865074],[-7.077578,0.8252182,7.195352],[-3.153342,-5.547021,6.458561],[0.3686684,-9.894412,9.951649],[1.689115,-1.891754,2.726141],[-4.762505,-8.633714,9.910726],[6.191169,-1.029447,6.35534],[5.478487,9.02305,10.60327],[-2.381128,-6.679209,7.161118],[0.8694938,1.342964,1.886683],[-1.479618,-2.535592,3.10137],[8.662861,-9.634588,12.99502],[-3.809545,-2.652849,4.748709],[4.292932,5.916048,7.377594],[-0.3986225,2.363383,2.597014],[-8.442391,5.360811,10.05049],[-1.901922,-8.864787,9.1215],[-0.3156001,7.236241,7.311825],[1.705611,-9.692336,9.891941],[-9.581105,8.171157,12.63192],[-9.992146,5.67351,11.53394],[5.47544,7.196645,9.09792],[7.156055,9.793629,12.17063],[5.616815,0.2763953,5.71183],[-7.906033,-0.5925217,7.991023],[3.280551,-3.067881,4.601511],[-2.436169,-3.429229,4.323718],[3.692118,-7.024133,7.998136],[6.384368,4.077352,7.641005],[-6.133347,-0.8226309,6.268545],[1.435963,-1.669516,2.418527],[-4.903309,-1.983239,5.382906],[-0.8198828,7.672681,7.78089],[-2.927693,-8.992625,9.509926],[3.293067,3.216767,4.710826],[6.500223,6.322356,9.122778],[-6.473027,-5.493163,8.548387],[-0.8160008,0.1336438,1.297582],[-2.711504,-8.749492,9.214437],[-9.991629,-4.805561,11.13221],[-7.331072,0.7562594,7.43751],[-1.695784,-7.786423,8.031442],[-8.593862,-8.500889,12.12929],[-2.587658,8.771084,9.199342],[-4.683978,-3.516619,5.941907],[5.394187,5.572039,7.819519],[-6.480431,0.5368044,6.579068],[1.141664,-4.216708,4.48152],[6.993794,4.528875,8.391893],[-6.397434,-9.912034,11.83958],[-9.476209,3.356431,10.10268],[4.301738,1.12315,4.557018],[-4.082791,-6.275594,7.553295],[1.525703,9.275139,9.452828],[-6.826807,-2.895879,7.482741],[-1.175551,3.201461,3.55405],[-3.677925,-8.575896,9.384729],[8.652877,9.630168,12.98508],[-3.274067,9.531051,10.12721],[-7.722453,-4.9268,9.214643],[3.781347,0.03422666,3.91149],[6.208839,3.547518,7.220427],[-3.652349,8.78375,9.565246],[-0.446119,0.5883811,1.243067],[7.720441,-4.937887,9.21889],[6.157681,-6.616125,9.093412],[4.341005,4.609493,6.410285],[5.130385,-1.593324,5.464388],[-7.642524,-4.509466,8.929919],[-4.883055,7.941676,9.376269],[0.7517663,7.243032,7.350283],[-1.604522,-5.258277,5.587841],[-2.188421,7.331568,7.716286],[9.119783,-4.371568,10.16273],[3.083964,-2.940983,4.377239],[3.437561,-3.895181,5.290488],[-8.163863,-4.848783,9.547741],[0.6851652,1.242979,1.736217],[1.009438,5.304302,5.491319],[9.957926,7.699821,12.62725],[1.001472,-3.677784,3.94069],[4.918794,-8.881055,10.20136],[-5.788497,8.994181,10.74253],[1.381,9.492838,9.644747],[7.120334,6.971377,10.01495],[4.899911,8.824444,10.14297],[0.05064856,7.516002,7.582404],[-4.266585,0.3669005,4.397541],[3.204297,9.932715,10.48458],[1.231031,-8.738987,8.881742],[-4.907528,3.52068,6.122011],[9.359031,-5.059325,10.68589],[6.2562,0.8751062,6.395769],[6.283878,4.647843,7.879694],[-2.570055,7.643559,8.125834],[-6.432033,5.068046,8.249615],[4.086849,-5.514896,6.9366],[2.950351,0.4111879,3.142236],[6.983173,-2.28801,7.416178],[8.102537,-4.798543,9.469801],[-7.845864,-3.315463,8.576122],[-7.096749,-5.40293,8.975271],[-7.437901,-8.455383,11.30557],[-4.561204,7.474215,8.812971],[-4.390179,-0.214652,4.507743],[9.229441,-8.587643,12.64635],[8.967897,-7.956067,12.03005],[0.264659,-2.943572,3.120042],[4.731253,-8.441511,9.728508],[1.193913,9.197316,9.328239],[-3.745926,6.616585,7.668844],[9.116789,-4.306604,10.13226],[-1.245416,-7.508409,7.67641],[-7.361577,5.392913,9.180213],[8.425856,-6.2451,10.53548],[-4.848618,0.822446,5.018517],[6.95961,-6.469832,9.554837],[6.434489,-5.877359,8.771887],[-3.37372,5.054373,6.158626],[-9.001662,-1.35842,9.158341],[-2.044367,-7.748425,8.075736],[-3.162754,4.331297,5.455561],[1.664819,8.809249,9.020781],[-6.46683,9.140206,11.24114],[7.956764,-7.237038,10.80207],[9.184902,3.552776,9.898719],[4.5175,-8.944674,10.0705],[-9.897186,7.741169,12.60476],[1.373851,-1.034538,1.989406],[-0.2811024,3.781803,3.921869],[7.023057,0.139555,7.095267],[0.1885484,-3.968369,4.096768],[-2.852903,9.584125,10.0496],[-5.129827,-8.878697,10.30274],[-5.098757,-5.021812,7.226058],[3.284099,9.989229,10.56267],[-9.069891,-9.509806,13.1795],[-6.721789,6.90901,9.691072],[-2.986058,7.610407,8.236191],[-2.463647,-4.868319,5.54708],[-4.547173,0.807107,4.725273],[-6.679747,-7.014564,9.737717],[-6.247161,9.64542,11.53521],[0.7707648,6.82713,6.942894],[-3.468883,2.358726,4.312394],[-6.761377,9.147274,11.41879],[7.079865,-3.720404,8.060143],[0.007705009,6.944863,7.016494],[-4.967925,-7.147223,8.761454],[-3.046752,-5.145854,6.06321],[-1.808607,-8.653306,8.896671],[5.466548,-5.95747,8.147061],[6.8967,1.710771,7.175737],[0.2790313,-5.683641,5.777684],[-8.680164,0.01532558,8.737591],[-0.06896258,-0.7183118,1.233178],[-8.829638,-5.288618,10.34079],[-9.513818,7.024251,11.86814],[-7.732433,-4.752549,9.131114],[-9.193524,1.798593,9.421031],[-2.032303,-2.823985,3.620103],[-0.8309576,1.60697,2.067086],[2.288458,-8.655633,9.008719],[5.747072,-9.592149,11.22667],[-0.8378792,5.315582,5.473341],[9.575475,6.467237,11.59805],[-8.789529,2.082773,9.088111],[5.219288,1.649848,5.564437],[-8.625629,-3.695544,9.437082],[6.040983,-1.133327,6.22719],[2.741824,6.415592,7.048221],[-8.485114,-3.696175,9.309074],[2.640933,5.277783,5.985776],[2.469792,9.504356,9.870799],[0.7398889,7.092432,7.200696],[-5.334754,-0.7022697,5.472913],[9.467709,-2.043286,9.737173],[-9.436759,2.404919,9.789589],[-6.294471,-7.161283,9.586675],[-6.481522,1.655963,6.764047],[-0.7341953,-5.560243,5.696959],[-0.5354129,-7.750636,7.8332],[7.79922,-8.954752,11.91702],[3.490864,1.368623,3.880626],[-8.162423,-7.03473,10.82186],[-7.297086,9.249571,11.82379],[2.439311,-4.893268,5.558265],[8.761016,7.196026,11.38148],[9.773023,8.033479,12.6905],[9.510495,-8.119328,12.54484],[-1.172419,-5.649761,5.856139],[-9.607784,-8.076482,12.59123],[3.292017,-2.411379,4.201443],[-0.7444665,-1.690541,2.100514],[-8.275818,-7.36447,11.12315],[8.012709,-8.94308,12.04916],[8.291071,2.01175,8.590052],[-6.172001,-5.438268,8.286637],[-4.105119,-7.221293,8.366545],[-0.3722795,-6.79686,6.880109],[-8.328864,2.214564,8.676074],[-1.544134,3.353922,3.82533],[1.986037,-5.987418,6.386981],[4.704031,0.1426945,4.811265],[1.148187,9.655889,9.775199],[0.7464197,3.047385,3.292977],[7.829722,-2.622285,8.317508],[-5.57087,0.4018295,5.674158],[-6.511107,-1.569089,6.771747],[0.7765609,7.785348,7.887629],[-6.186577,1.716758,6.497769],[5.077928,2.393133,5.701968],[-5.973591,4.570589,7.587758],[-6.099952,1.132228,6.284214],[-2.800604,5.364227,6.133377],[0.3003373,-6.732472,6.812957],[3.639672,-5.984831,7.075692],[7.724808,-8.724241,11.69551],[4.027124,7.783534,8.820495],[-0.7645735,-3.550792,3.767319],[-4.085245,0.4775698,4.232883],[0.9195137,-7.841184,7.957994],[5.101462,9.192356,10.56051],[0.5780753,5.777665,5.891993],[-2.561065,8.471248,8.906239],[-0.6121992,-0.8301218,1.436625],[7.901653,0.8812317,8.013282],[-4.573014,-4.305457,6.359985],[9.919281,9.452029,13.73801],[1.312574,-7.84256,8.014275],[-6.543033,3.614633,7.541675],[7.086831,-8.570868,11.16615],[0.6757597,-1.361085,1.819121],[6.717368,-1.815259,7.029808],[-4.229144,-1.43209,4.575646],[4.523955,9.938449,10.96535],[-7.598751,-5.265216,9.298575],[1.18102,-5.004004,5.23783],[0.763927,-1.085117,1.661645],[1.211819,-4.242982,4.524534],[4.850451,2.118425,5.38652],[1.962341,4.892951,5.365795],[-6.607779,9.919451,11.9607],[1.221171,9.297963,9.43098],[0.5271804,-7.09757,7.187031],[-9.678009,-2.703565,10.09817],[-9.09907,9.817823,13.42322],[6.683273,-9.351516,11.53763],[-5.314906,-6.976854,8.827498],[0.3049673,-6.528482,6.611663],[0.9843282,-9.253682,9.359463],[-7.665352,-4.591029,8.990838],[6.307841,3.232333,7.157991],[-9.233768,-1.716256,9.444998],[-5.466522,-0.1988414,5.560791],[-0.184316,-5.052709,5.154012],[-2.569852,8.400002,8.84105],[6.169427,-1.700746,6.47722],[-0.6498017,3.140564,3.359373],[6.435998,4.323072,7.817354],[6.803027,8.09551,10.6216],[6.352477,-9.034736,11.08965],[-1.895594,3.952019,4.495746],[-5.363084,-4.058688,6.799677],[-7.084423,-4.74158,8.583219],[1.712365,-8.505388,8.733488],[5.27159,-6.968347,8.794744],[-6.209102,-9.944238,11.76609],[6.63827,3.555931,7.596794],[-7.33247,-8.580679,11.33107],[9.697966,6.691155,11.82464],[-0.1269867,-5.709984,5.79828],[4.673099,-8.732858,9.954931],[-9.997308,7.459841,12.51381],[-9.341907,2.134025,9.634588],[5.259044,2.48706,5.902796],[-1.301864,6.548519,6.751144],[-0.9596233,-6.113801,6.268927],[-2.595451,1.381586,3.105663],[3.437685,0.4200046,3.60473],[0.1655611,-4.346909,4.463521],[-4.053598,-2.621016,4.929643],[4.945438,-0.5255925,5.07283],[-9.11295,5.326234,10.60258],[3.455571,-4.370979,5.660957],[-0.7705561,2.192904,2.530333],[3.451024,-6.308059,7.259557],[2.798734,-6.770386,7.393986],[9.738323,6.916812,11.98654],[4.443854,5.760645,7.3439],[-8.23856,8.830943,12.11856],[8.855104,1.932678,9.118559],[4.072318,-3.769139,5.638278],[-4.8446,-6.335761,8.038159],[3.934173,-0.7636639,4.130485],[-2.36902,6.872863,7.338154],[-9.166941,6.138456,11.07761],[-4.464229,-3.87862,5.997752],[0.8175203,-7.306541,7.41983],[5.720234,1.76786,6.070125],[-7.647174,-6.145559,9.861398],[2.924499,-7.046849,7.694853],[-0.08725037,0.8386303,1.308019],[-1.645242,-8.13881,8.363436],[7.236998,-7.617834,10.55488],[-7.849757,-5.640213,9.717546],[-9.239086,-6.287152,11.22003],[-7.492447,-8.612645,11.45925],[1.577802,-4.511866,4.883276],[-8.308249,-4.074402,9.307405],[6.281799,8.378735,10.5197],[-6.022283,-1.762617,6.354109],[3.423368,6.074414,7.044002],[-2.788911,-3.515841,4.597734],[1.612056,-8.943851,9.142823],[-2.370699,1.976176,3.244301],[-7.317078,-0.7115886,7.419298],[-4.053985,9.608095,10.47618],[9.33849,-2.678935,9.766478],[6.633677,3.807328,7.713716],[1.348844,-3.198411,3.61237],[-9.053053,-0.3889849,9.116418],[5.00709,0.2592806,5.112551],[6.957533,9.824132,12.07977],[-3.684072,5.841953,6.978596],[-0.2403985,-6.789228,6.866689],[6.263057,-1.620401,6.546112],[5.962756,2.202227,6.434614],[3.120473,-3.160677,4.552718],[-5.356354,5.635052,7.838644],[6.988067,2.059327,7.353497],[1.159495,-5.385854,5.599272],[-3.733474,7.839953,8.740921],[-8.07836,-9.467192,12.4855],[-1.187566,9.108995,9.240353],[-0.3347791,-0.3319912,1.105575],[-1.532578,9.994389,10.16054],[1.791123,-9.268029,9.492338],[0.2999002,-1.795169,2.076673],[-2.146001,2.822607,3.684078],[-0.7468649,5.325447,5.469752],[-4.151899,-8.718081,9.707893],[5.516128,-9.466512,11.00193],[1.082967,-8.782898,8.905735],[5.62564,8.740252,10.44221],[-7.96995,9.348073,12.32504],[-3.078927,-8.181225,8.798422],[-9.753015,5.781401,11.38182],[-5.694021,-7.643212,9.583348],[-1.213917,-5.71652,5.928929],[-8.153647,3.40878,8.893916],[0.3021806,1.682104,1.980097],[-4.732285,-0.1283279,4.83849],[9.423412,9.444643,13.37916],[0.3246038,-0.7324023,1.28132],[-7.070089,-0.2012226,7.143294],[1.244143,-6.826352,7.01049],[-7.451756,0.2967781,7.524409],[0.1533905,5.279818,5.375873],[-8.853349,-1.81026,9.091689],[6.714858,9.938375,12.0358],[1.911941,-3.718906,4.299509],[6.890465,-1.261867,7.076074],[-7.603061,8.142638,11.18522],[4.589641,2.88841,5.514319],[-5.827999,-9.055609,10.81525],[3.561408,4.197209,5.594657],[-0.2642165,-2.52885,2.732196],[9.186699,-0.1297752,9.241877],[-4.694407,3.880183,6.171975],[-6.886331,5.598645,8.931202],[6.524509,9.721381,11.75051],[-2.283896,-0.04342385,2.493605],[-1.150239,7.654446,7.804717],[6.439602,-9.009646,11.11945],[-2.681404,-2.405801,3.738691],[-0.0925548,6.323695,6.402943],[-1.396544,-8.747549,8.914592],[2.429743,-9.046222,9.420073],[7.490023,7.358797,10.54762],[1.558008,8.090616,8.299726],[-2.316818,-7.004644,7.445313],[-1.992821,6.97081,7.31871],[-8.845498,-7.167365,11.42865],[5.467544,1.315682,5.711834],[5.086253,-9.110887,10.48228],[1.89333,-1.11359,2.413458],[-8.265115,4.677238,9.549275],[-2.081071,9.787651,10.05629],[-3.15471,5.953996,6.811921],[8.770033,1.882806,9.025433],[-5.30592,0.933628,5.479457],[9.350909,-3.563342,10.05668],[1.247216,2.133069,2.66562],[5.019186,9.93304,11.17397],[4.732299,-4.283632,6.460971],[6.000791,-0.5384012,6.107321],[9.572927,-3.306265,10.17705],[5.745532,-4.120096,7.140471],[1.614259,-3.60828,4.07744],[2.263816,0.5671461,2.539],[-9.512125,-3.452347,10.16854],[-5.818967,3.038677,6.640326],[0.2256227,-5.920393,6.008491],[6.063024,-0.2007696,6.148217],[0.9010613,6.393887,6.534041],[-9.705169,-8.847168,13.17052],[8.509007,-1.111046,8.639307],[-0.4304795,6.827981,6.914235],[5.593645,-4.386387,7.178388],[2.423759,-7.813915,8.24208],[8.499418,6.178888,10.55551],[-1.129305,-4.586727,4.828394],[-5.961918,9.796649,11.51168],[-7.770928,8.702407,11.70979],[7.784481,5.324234,9.483966],[1.85839,6.88756,7.203617],[3.193527,0.4917679,3.382373],[-5.713337,-3.215566,6.631899],[-9.179932,4.847199,10.42912],[-7.297921,5.082871,8.949594],[-3.967735,6.736189,7.881571],[9.549225,-4.695562,10.68812],[-0.5834654,3.644828,3.824291],[7.575147,-2.321303,7.985693],[0.4854755,-8.296776,8.370913],[5.076977,-8.943995,10.33299],[0.1283744,8.624952,8.68368],[9.15706,-5.885106,10.93097],[-9.639448,-4.320027,10.61045],[-1.899052,-2.159673,3.044764],[0.7542663,-4.671247,4.836266],[-9.941556,-4.640815,11.01688],[-6.824981,-1.708938,7.106393],[1.539199,2.127271,2.8097],[-8.997989,9.704909,13.27212],[9.203363,7.57479,11.96158],[-9.45731,-8.032806,12.44856],[-5.893833,-2.751062,6.5807],[-3.483333,3.995053,5.393891],[-3.832588,5.008461,6.385406],[-9.673155,-0.8475894,9.761574],[-0.4823169,-4.889526,5.01399],[-7.898379,-5.849718,9.879454],[0.4478993,-3.592537,3.75592],[5.756533,0.6527089,5.87909],[-9.558059,-9.265389,13.3493],[-5.337088,7.261429,9.067131],[3.022364,2.362781,3.96452],[-2.803579,-1.145029,3.189224],[5.231949,3.637268,6.450039],[-7.605154,3.36605,8.376674],[-5.593075,0.965758,5.763261],[2.584531,-9.319409,9.722715],[8.798114,-8.216676,12.07976],[-0.02406261,9.64856,9.700273],[-6.234048,-5.630423,8.459611],[-8.994363,5.549534,10.61583],[-3.934919,-2.230843,4.632521],[-8.099369,-4.055617,9.113057],[-3.073812,-0.894674,3.353917],[5.006932,1.272436,5.261982],[-4.978136,-1.991595,5.454199],[-6.235122,-9.192406,11.15245],[1.468492,-5.148221,5.446159],[-8.09447,-1.517537,8.295985],[-2.659615,4.49363,5.316603],[-0.07068485,-5.309044,5.402865],[-5.053023,6.671087,8.428312],[8.817739,8.855864,12.5371],[-3.101916,2.781476,4.284681],[-7.93964,-8.167303,11.43428],[-0.9734198,9.840242,9.938708],[-9.483519,4.263031,10.4456],[9.651171,3.495205,10.31317],[-4.554325,2.838437,5.458809],[-9.797226,-5.483606,11.27189],[9.550074,-7.252088,12.03315],[-6.543791,-7.173251,9.760981],[9.731215,-6.716716,11.86637],[8.721338,-1.530789,8.910952],[3.13209,-5.897849,6.752378],[0.3623556,-1.335218,1.707076],[-8.048385,7.362542,10.9537],[-5.016875,-0.2769558,5.12306],[-5.434774,-7.253398,9.118583],[4.445865,1.895331,4.935382],[-9.350907,4.069662,10.24703],[9.854008,4.845221,11.02622],[9.244578,7.733581,12.09423],[0.5292944,-5.874404,5.982372],[3.187846,-8.898858,9.505369],[3.320816,-4.07546,5.351373],[9.593668,2.810463,10.04675],[7.93499,-7.807923,11.17711],[3.654922,8.270161,9.096924],[-2.659665,3.888083,4.815704],[-7.971523,-5.835155,9.929461],[-0.2112339,-7.139644,7.212429],[9.417537,5.670951,11.03855],[2.743155,6.495494,7.121541],[-0.8812202,6.617611,6.750506],[-4.057792,-5.683664,7.054765],[2.872254,7.472265,8.067502],[3.256895,-8.876364,9.507745],[2.592221,-6.163302,6.760614],[3.790209,-2.989434,4.929747],[0.8941025,-8.761087,8.863186],[-2.419285,7.437423,7.884681],[9.391314,8.886251,12.96774],[9.130524,1.489739,9.305148],[-8.796959,3.278603,9.441172],[-6.78279,8.60393,11.00154],[-1.603462,2.294063,2.972173],[-7.788879,-9.776871,12.54009],[8.747035,-5.921741,10.61026],[2.953391,-7.877439,8.472105],[-0.4669046,5.02558,5.145333],[-2.517196,-5.847114,6.443991],[-3.318007,-8.668057,9.335115],[8.285123,-1.783755,8.533759],[-2.359167,1.907107,3.194171],[-3.420571,-0.5741838,3.609708],[6.117384,-4.891331,7.896043],[-3.323358,-6.814673,7.647515],[-2.772804,6.798331,7.409841],[6.757722,-5.501614,8.771235],[1.985038,9.786057,10.0353],[-4.235989,-0.4532931,4.375966],[-1.379007,9.358368,9.512136],[7.865448,6.428217,10.20722],[-9.545401,1.731019,9.752492],[-3.885469,-1.711608,4.361935],[0.08586908,-5.063686,5.162199],[-2.206097,3.23831,4.043948],[0.5418423,3.972797,4.132397],[2.241916,5.372716,5.906967],[-0.3032296,-6.6926,6.773688],[-0.8791466,-1.77591,2.21963],[8.921905,0.6323236,9.000011],[4.798172,-8.797273,10.07047],[1.758845,9.439041,9.653446],[6.645394,7.756927,10.2631],[-4.502148,4.034481,6.12751],[6.597911,-8.129627,10.51776],[6.904244,3.074252,7.623622],[4.74366,8.657477,9.922411],[-4.986433,-1.591686,5.328975],[2.504702,3.565574,4.470666],[8.868769,6.908651,11.28648],[-2.612949,-7.628025,8.124917],[-1.770458,2.79946,3.459985],[3.176291,9.138727,9.726519],[5.820123,7.270753,9.366839],[-7.045428,-6.689942,9.766953],[5.595902,2.445398,6.188222],[-3.557668,-1.729897,4.080386],[-2.964814,-8.73433,9.277858],[8.511715,9.360169,12.69102],[1.49765,-2.233779,2.869272],[5.886727,3.404234,6.873309],[0.1698984,-5.785687,5.873929],[6.011628,-9.377918,11.18414],[9.787042,-5.612803,11.32651],[4.164275,3.74994,5.692384],[0.1928364,7.256134,7.327255],[6.299042,-0.1753928,6.380336],[2.708809,-1.630122,3.315862],[-3.33498,-7.812174,8.552903],[5.92335,-9.181023,10.97166],[7.153953,-3.932617,8.224629],[-2.555628,-9.798676,10.17572],[-8.616889,-8.320122,12.01978],[-3.388689,1.604545,3.880435],[0.1664538,-6.354692,6.435047],[4.092364,2.89539,5.111823],[4.581472,7.261191,8.643771],[-8.448915,6.270019,10.56869],[-9.050406,2.198597,9.36716],[7.100004,-4.617297,8.528158],[1.532559,-6.326555,6.585897],[0.3410364,6.253156,6.341787],[-3.034579,2.828206,4.267015],[2.510525,7.781698,8.23757],[0.6797224,-9.960518,10.03364],[6.942195,-3.740421,7.948888],[1.313348,-3.060961,3.477695],[7.712651,-0.9719962,7.837714],[-8.169613,8.937236,12.14976],[-7.42193,-9.612754,12.18565],[-6.911743,-7.098249,9.957777],[-3.42612,-1.973306,4.078264],[8.868762,2.612778,9.299545],[1.404593,9.097745,9.259689],[0.391002,5.885264,5.982409],[-3.756169,9.166971,9.957016],[-6.264457,-8.932402,10.95588],[-9.242496,-7.695589,12.06838],[-7.186439,7.193016,10.21687],[2.961105,2.994728,4.328572],[6.609887,-3.500645,7.546199],[4.264806,-2.959003,5.286234],[-4.640473,2.084695,5.184587],[4.858999,-0.2942559,4.969553],[2.290388,7.652032,8.049812],[-4.097424,9.383471,10.28778],[6.456778,7.940283,10.2829],[-1.944724,3.691883,4.290915],[-2.573021,8.341157,8.786088],[4.895772,7.386437,8.91785],[5.019363,7.118575,8.767447],[-6.335101,-9.684591,11.61571],[-3.062484,-8.666995,9.246384],[-9.49537,-6.591806,11.60233],[-0.3846353,-4.60069,4.723801],[1.142475,0.2174905,1.533803],[7.398667,9.363044,11.97526],[-6.941974,5.502146,8.914292],[-4.058337,4.865401,6.414221],[7.360017,-0.8161455,7.472345],[3.92758,-0.9051989,4.152742],[-3.012098,-1.90766,3.702958],[-4.877238,-9.336536,10.58104],[4.706364,2.817341,5.575596],[8.601871,-8.671347,12.25498],[-1.037801,1.899425,2.384291],[-3.150777,-2.664414,4.245762],[-0.1070798,9.052677,9.108372],[-3.413745,-9.096563,9.767349],[-4.246274,-6.023652,7.437421],[9.991524,-1.478807,10.14975],[-3.474917,5.728939,6.774643],[-2.42769,-8.389099,8.790373],[-7.912705,-6.337269,10.18685],[-3.54007,0.2667708,3.68826],[4.744718,-7.252232,8.723945],[0.6642042,-2.329581,2.620709],[-6.583406,2.233697,7.023578],[1.251512,-2.532171,2.99636],[6.695128,0.7445791,6.810223],[-9.863089,3.59788,10.54634],[-6.863056,9.171932,11.49895],[4.911162,6.69242,8.3611],[8.534121,-6.095148,10.5348],[-3.327712,-1.788103,3.90781],[1.797744,6.78762,7.092508],[-1.130452,9.377501,9.49818],[-4.917552,3.640016,6.199357],[5.255393,-3.457038,6.369479],[0.5862154,-0.9307338,1.486578],[-4.845455,6.719533,8.344493],[-1.665051,1.971965,2.767858],[1.668697,-6.764209,7.0384],[8.989677,-9.278499,12.95781],[-9.911304,8.471678,13.07682],[-9.185816,-3.801953,9.991699],[3.65679,6.145722,7.220942],[-3.864568,-4.718306,6.180396],[3.444493,4.613634,5.843813],[-0.5813892,-3.966104,4.131343],[-3.647323,-7.949504,8.80327],[-8.965711,-0.7316296,9.050925],[-7.944056,-2.004438,8.253835],[7.634384,1.499874,7.844325],[0.8717657,1.819062,2.251436],[1.332752,-7.18169,7.372442],[-2.394544,2.190699,3.396028],[-0.5636058,2.444658,2.700742],[2.439452,0.4832851,2.68039],[0.4899778,-8.655725,8.727065],[-3.847348,-7.299609,8.311821],[-5.427393,1.415544,5.6974],[5.742879,-7.221491,9.280657],[6.85184,-5.342892,8.746097],[-7.209309,2.649543,7.745594],[9.364814,-0.8147883,9.453234],[-2.11233,-2.692232,3.565116],[-6.783095,-9.762004,11.92925],[7.335895,-3.350999,8.12678],[7.280794,-4.636864,8.689676],[3.285004,9.469914,10.07326],[-0.1920985,-2.857993,3.033979],[9.973012,9.134368,13.56089],[-5.389814,6.13601,8.228045],[-0.1878814,6.142237,6.225944],[-1.76835,3.948017,4.440034],[7.774332,0.4561277,7.851642],[2.789953,-1.354859,3.258754],[-2.84064,-2.637895,4.003464],[-0.7772512,-3.396895,3.625329],[-9.319514,8.001427,12.3238],[5.554537,6.716581,8.772989],[-4.458842,8.110513,9.309226],[7.119703,-5.845469,9.26605],[-3.230946,7.751057,8.456825],[8.420546,-8.126653,11.74513],[0.6405681,-6.031434,6.147237],[6.174857,3.04193,6.955731],[-9.355068,1.001159,9.461481],[-7.499075,-3.906469,8.514495],[2.178648,2.056628,3.158517],[-1.123759,-8.675184,8.804638],[7.717688,0.9430032,7.83913],[-4.613144,-7.127548,8.548862],[-9.622925,6.954903,11.91517],[-5.8213,6.20454,8.566437],[-0.1227936,-7.073433,7.144826],[-6.642548,9.062164,11.28035],[7.850344,4.545984,9.126548],[-9.061531,-1.978853,9.328837],[7.64973,3.433202,8.444243],[-2.459876,2.20541,3.451786],[-0.7093611,5.221506,5.363517],[-4.494261,-7.562414,8.853728],[-7.474944,-8.238286,11.16889],[2.419461,-3.063378,4.02965],[1.423836,5.595123,5.859413],[5.855127,7.504602,9.570871],[5.090905,-9.420571,10.75474],[4.925076,-1.607487,5.2764],[-6.618388,-3.037699,7.350556],[-8.779198,-9.589725,13.03983],[7.780836,6.638401,10.27666],[-6.239818,2.143626,6.673114],[4.4191,-8.213828,9.380588],[-6.259062,3.395531,7.190653],[-4.803408,-7.567867,9.019165],[3.610719,-5.29274,6.484627],[5.953036,5.581099,8.22115],[8.361509,2.938143,8.918941],[-3.490412,3.505813,5.047147],[9.895079,5.529818,11.37943],[-1.349115,4.531281,4.832456],[3.557745,3.142116,4.850819],[3.330611,4.211916,5.461978],[4.550493,-7.957925,9.221473],[-0.7844535,5.214839,5.367487],[-6.271703,-5.109148,8.150929],[-1.272887,6.883322,7.071094],[5.594701,-5.622087,7.994282],[-2.620367,-4.740732,5.508254],[-4.122531,5.4164,6.879873],[-3.502323,2.652334,4.505679],[-9.312833,5.852195,11.04432],[-5.552279,-8.776051,10.43297],[-9.879997,-5.662083,11.43125],[-8.816473,-5.062438,10.2156],[-3.204977,7.870223,8.556418],[0.1417008,3.800217,3.932141],[0.06335647,3.322359,3.470171],[-4.098778,8.733688,9.699346],[3.975079,9.814874,10.6364],[-6.032164,3.494934,7.042838],[-1.838814,8.486774,8.741086],[-5.199818,9.895197,11.22288],[-5.039094,5.492867,7.520908],[-9.092926,1.972321,9.357957],[-3.910808,-1.131015,4.19209],[7.351281,-2.762437,7.916589],[-7.429469,0.8545226,7.545012],[5.32816,2.733668,6.071427],[8.336119,-2.076338,8.648819],[-4.949964,5.482555,7.453895],[-0.5404385,-4.246474,4.395977],[-5.546232,1.800582,5.916316],[-4.096731,5.247416,6.731908],[2.277301,-5.843974,6.351231],[-0.08136605,6.326296,6.40536],[-0.05548909,-0.272862,1.038043],[3.449883,-7.058935,7.920244],[-8.220737,-6.501094,10.52828],[9.720039,3.000892,10.22177],[-4.23629,-6.534389,7.851395],[4.321012,2.120915,4.916241],[-7.352894,6.86477,10.10891],[-0.7313371,2.796905,3.059008],[0.03783916,-7.310223,7.3784],[-1.94786,6.92364,7.261608],[3.261847,4.689376,5.799129],[-3.068762,-9.041254,9.600082],[2.135581,5.654165,6.126197],[-0.195546,-2.116759,2.349235],[-6.775344,-0.2905046,6.854902],[-9.096338,2.709102,9.543721],[-8.068222,-1.248489,8.225262],[-0.2166997,7.904926,7.970873],[-3.523556,0.4324432,3.68815],[-1.788757,-6.723916,7.029274],[1.688889,-6.234468,6.536125],[-1.888813,-1.475228,2.596904],[6.561693,-5.217651,8.44273],[-3.822665,5.732094,6.962017],[-7.527582,4.65997,8.909534],[-6.254177,-4.682565,7.87662],[6.097456,2.095912,6.524708],[7.425459,-5.158817,9.096748],[5.240229,-7.262443,9.011274],[-3.688795,-1.371752,4.060654],[6.094103,7.276927,9.5442],[5.572519,4.011515,6.938676],[-4.536556,2.54199,5.295475],[3.830873,-5.836015,7.05228],[-3.490408,-0.6514847,3.688818],[0.6689489,2.583772,2.850152],[-0.5533091,-4.075671,4.232877],[2.892832,-8.678652,9.202579],[-7.840051,9.622859,12.45254],[-7.77427,-7.468882,10.82698],[-6.73363,-6.20404,9.210423],[1.490287,-5.859417,6.128109],[4.832642,4.855544,6.923203],[2.096459,-0.6882451,2.422565],[-9.093729,2.450173,9.470969],[5.009709,-4.339595,6.70293],[-5.625989,2.443484,6.214689],[-5.524122,3.950332,6.864477],[-0.557893,-2.588303,2.830293],[-2.81859,-8.517387,9.0272],[-9.861232,-5.007495,11.1049],[-3.442376,0.4933105,3.618467],[-6.51197,5.732262,8.732959],[-2.567187,-8.674474,9.10148],[-9.838678,-4.573789,10.89583],[-8.304737,2.442508,8.714041],[8.686109,-4.089959,9.652784],[-4.071528,-1.877307,4.59365],[-5.66749,3.705119,6.844585],[2.739939,-2.204227,3.655939],[-7.402947,9.09726,11.77131],[-3.331964,9.658257,10.26567],[-2.715995,-8.826148,9.288569],[7.619859,-1.003143,7.75039],[2.637368,-8.233363,8.703102],[-8.255214,6.466034,10.53367],[-7.845716,-7.342427,10.79196],[4.256182,-3.063219,5.338389],[-0.4452148,-4.077673,4.222042],[-0.2529548,-8.993378,9.052339],[-0.09896415,6.56864,6.64506],[-3.694687,7.597183,8.506932],[8.682541,-6.771453,11.05618],[4.992578,-7.588313,9.138289],[-7.178981,-1.576249,7.417704],[4.78589,5.519712,7.373734],[-7.086592,-6.263816,9.510793],[8.873513,-5.174444,10.32057],[-2.848265,-4.973245,5.817713],[-7.047772,-1.007545,7.189314],[-0.6243341,7.820093,7.908454],[6.051969,-2.638283,6.67734],[7.438481,7.134314,10.35516],[9.50623,9.146602,13.22984],[9.97924,-4.014916,10.803],[-5.946162,-0.5519149,6.05487],[2.147895,9.975302,10.25281],[0.141526,-1.999357,2.239968],[0.4501656,-8.724562,8.793215],[-9.540568,6.827865,11.77464],[-8.984278,5.37029,10.51462],[5.436249,-3.372957,6.475311],[1.826314,2.444986,3.211445],[9.689237,4.694834,10.81308],[1.876511,-4.409919,4.895782],[7.092808,-2.311072,7.526551],[9.882794,-6.344429,11.78649],[-4.095236,-4.329861,6.043067],[-3.304761,-3.615643,4.999432],[-4.882585,-8.802597,10.1156],[7.636311,-6.466364,10.0562],[-5.914627,6.736696,9.020304],[1.154734,9.60603,9.726727],[-6.188161,-9.244102,11.16901],[-7.238706,1.375,7.43569],[-5.929667,0.5634254,6.039735],[2.143885,-8.820924,9.13263],[9.765995,-1.195538,9.889588],[-3.693328,-7.949773,8.822673],[-9.812684,4.73161,10.93969],[-6.768584,2.994825,7.468782],[4.632979,8.856627,10.04512],[4.497307,6.872924,8.274229],[3.718074,-0.5565043,3.890214],[0.05082548,-2.010869,2.24637],[7.919585,8.864242,11.92873],[-8.379777,6.701325,10.77629],[-1.497805,7.458305,7.672661],[5.45614,-3.320751,6.465048],[5.295333,-2.567695,5.969389],[5.270307,0.7327378,5.414151],[1.573168,-0.2496579,1.880741],[9.907186,4.390967,10.88269],[4.808665,4.791992,6.861957],[-4.971204,-9.200284,10.50515],[-6.289574,-8.584452,10.68885],[4.519662,7.105654,8.480429],[-2.112538,5.293636,5.786657],[-0.4803076,-5.043202,5.163776],[-1.850006,9.829858,10.05229],[-7.175415,-8.742315,11.35406],[-8.878262,3.980387,9.780951],[0.3798491,4.312988,4.443664],[-7.756751,-6.411272,10.11294],[1.832332,4.953068,5.374972],[-8.392366,-0.2724867,8.456125],[-5.938877,-2.669058,6.587422],[6.633727,-2.830252,7.281254],[-4.232026,-8.653081,9.68431],[-5.738872,4.356779,7.27435],[7.488274,-3.466636,8.312149],[5.672551,9.58562,11.18311],[-4.29986,1.006204,4.52783],[-9.891646,-8.813272,13.28602],[5.180959,8.622006,10.10848],[-0.9182124,5.507638,5.672494],[-6.730693,-0.9770041,6.874356],[1.943328,-6.803111,7.145547],[7.279415,8.154059,10.97627],[9.459602,2.842313,9.927881],[-0.8713915,-3.467912,3.712915],[-1.161102,-2.629477,3.043404],[8.588593,-4.692749,9.837979],[5.127622,5.849792,7.842995],[9.655338,6.721909,11.80718],[-3.034741,-0.6952713,3.270024],[-5.020337,-9.267941,10.58766],[-5.876033,-3.129963,6.732343],[2.182881,4.71429,5.29051],[8.415937,8.34199,11.89188],[1.718335,9.986978,10.18295],[-4.330139,3.054002,5.392313],[-4.404943,-7.833227,9.042288],[6.680939,6.717304,9.526653],[9.290026,-0.9768758,9.394619],[-2.294182,2.286829,3.390112],[-1.25054,-2.082019,2.626529],[-6.616131,2.891484,7.289298],[-5.362985,5.490157,7.73973],[3.664916,9.774693,10.48695],[6.047369,-7.060947,9.350276],[4.80045,-0.6764542,4.94994],[4.717824,6.514107,8.105026],[8.077683,9.202474,12.28554],[-7.436915,7.063358,10.30528],[6.764772,-4.243815,8.048113],[-2.312011,6.492326,6.963885],[3.466763,5.979393,6.983666],[4.704898,9.057532,10.25548],[0.0009140465,-2.245543,2.458143],[6.038418,6.449228,8.891291],[5.908764,6.322896,8.711631],[5.351612,3.37442,6.40519],[-2.423498,-1.105522,2.845263],[8.913103,-7.881099,11.93965],[4.329096,0.3101319,4.453903],[-4.257115,9.459668,10.42153],[-4.509526,0.5972896,4.657529],[-2.921969,-7.095994,7.73893],[2.448617,-2.469134,3.618335],[8.069372,8.819777,11.99597],[-3.830343,4.18327,5.759451],[7.19094,4.887028,8.751723],[-9.12351,7.836963,12.06882],[5.020235,6.85226,8.553142],[-6.438981,7.45382,9.9005],[-1.828428,2.24466,3.062947],[-8.200793,8.773265,12.05086],[3.611061,6.420071,7.433511],[3.391154,-6.302773,7.226677],[1.024138,3.735057,3.999939],[-6.199077,-5.500251,8.347533],[-6.50361,0.01360957,6.580056],[-8.59259,1.730311,8.821938],[-1.86308,-4.570709,5.036115],[8.346536,3.361562,9.053439],[5.894311,-9.162116,10.94017],[2.356128,7.656633,8.073126],[-7.344551,-2.259099,7.748932],[9.428171,8.666076,12.84489],[6.465826,-1.155618,6.643972],[-2.730127,-5.38705,6.121593],[8.898476,-7.796278,11.87286],[3.97973,-7.857392,8.864359],[-0.5253742,-6.519621,6.616757],[5.266454,-0.7857148,5.41783],[-2.245651,-2.443093,3.465783],[6.851363,-7.737328,10.38303],[9.907389,-3.301506,10.49077],[-5.17101,-5.980442,7.969004],[-7.363061,2.540963,7.853099],[7.536949,-8.093553,11.10456],[-6.090154,-3.405526,7.048942],[4.208602,-4.093874,5.955849],[9.239588,-2.31151,9.576694],[3.704305,-0.04474995,3.837171],[3.746574,-6.394134,7.478086],[7.091271,-0.8615729,7.213073],[-1.377693,2.850806,3.320412],[7.856383,-6.388833,10.17546],[7.206998,-8.714416,11.35262],[1.107793,1.538144,2.14315],[7.594915,-2.443495,8.040733],[-9.439507,3.776589,10.21601],[6.987806,4.10987,8.16826],[-1.100139,-6.552429,6.718976],[7.901615,9.072757,12.07271],[-1.047432,-1.258963,1.918881],[-6.637045,-7.177428,9.826792],[-4.274283,2.32119,4.965623],[8.637917,3.074792,9.223228],[-2.338374,9.54439,9.877417],[6.477477,9.641287,11.65814],[-0.1228338,0.2732408,1.04391],[2.039651,3.488674,4.163055],[8.574807,9.511909,12.84538],[-2.475487,-0.1553291,2.674352],[5.537508,-3.076454,6.413156],[-7.411102,1.918473,7.720425],[3.441199,7.849268,8.628607],[1.622842,-2.142086,2.867429],[8.733369,6.310763,10.82116],[3.50747,-6.237132,7.225244],[7.791735,-4.564802,9.085623],[1.006766,-5.694939,5.869064],[6.58714,-5.939056,8.925402],[-7.626366,-8.630653,11.56069],[-1.09759,-0.7745461,1.674702],[-7.263681,3.792223,8.254818],[5.252061,3.331281,6.299332],[0.299731,0.8084162,1.320369],[-6.907551,-2.877182,7.549334],[9.049226,4.672171,10.23316],[-1.755653,-3.012262,3.627126],[5.368212,-1.276995,5.607888],[-3.790822,-1.044921,4.057363],[-2.534649,-9.813009,10.18428],[-0.8320971,-1.256115,1.808372],[-7.367116,1.811843,7.652267],[5.910938,8.51913,10.41704],[-9.294811,0.9440248,9.395994],[1.397336,6.107323,6.344441],[5.0289,-2.886958,5.884247],[-1.481095,3.1196,3.595212],[-3.374809,-7.506688,8.29094],[2.097844,-7.741967,8.083254],[4.176259,-3.052999,5.26896],[4.399667,-9.204826,10.25114],[8.230752,-7.325108,11.06357],[-6.059896,-9.819654,11.58223],[-0.08060724,8.285413,8.345931],[7.63002,-5.777504,9.622721],[1.726767,-6.552061,6.849176],[-9.752542,3.579413,10.43668],[1.337365,-5.146305,5.410453],[5.459431,9.382244,10.901],[8.37096,9.954946,13.04507],[3.340443,3.444793,4.901546],[-9.363785,4.336634,10.36759],[2.395782,9.295922,9.651628],[8.823208,2.752994,9.296664],[-0.6886192,-1.502559,1.931808],[9.010512,4.962944,10.33538],[6.086441,8.263551,10.31169],[-3.601287,6.955719,7.896284],[-1.613504,7.272275,7.515942],[7.673546,-6.556836,10.14275],[6.126033,3.768699,7.261637],[-5.540016,4.993939,7.525371],[-8.064345,-4.729247,9.402098],[-8.790692,-3.463645,9.501216],[7.904502,4.737026,9.269336],[7.258043,0.6505223,7.355431],[-2.991899,6.975965,7.656078],[-3.612012,-4.467862,5.831674],[-3.342825,7.693839,8.448054],[-5.722999,-5.292297,7.858825],[-4.313746,-5.644767,7.174385],[-9.754655,-4.797706,10.91656],[4.753129,6.855307,8.401635],[-3.543037,-8.214815,9.002016],[0.1988629,-7.237789,7.30925],[7.66881,-5.22968,9.335962],[7.965008,-1.604094,8.186237],[-2.211025,7.419683,7.806429],[-2.410018,0.536045,2.663744],[-9.794659,-0.7047302,9.870764],[8.355259,-0.2084167,8.417469],[0.2095964,7.204988,7.277072],[0.1958178,-0.4305041,1.1062],[7.030733,-8.234029,10.87338],[-4.370622,6.987303,8.302093],[-6.548283,-8.296862,10.61687],[-5.439164,-4.701154,7.258467],[-8.803233,2.5812,9.228191],[-7.174647,-1.497169,7.397099],[3.830369,1.143819,4.120685],[-7.39859,-2.694618,7.93726],[7.388648,8.247266,11.11798],[7.533119,-6.555408,10.036],[8.33352,5.744198,10.17071],[-9.673656,9.621207,13.68018],[6.974899,6.813864,9.801936],[-9.47474,3.111123,10.02246],[6.185612,0.2417774,6.270586],[4.286698,-3.347763,5.530217],[-7.966925,-3.256915,8.664836],[-3.247722,-7.914569,8.613252],[-6.634429,-6.054555,9.037328],[4.694387,3.311934,5.831481],[0.3181615,7.082694,7.160012],[3.778829,1.983446,4.383333],[7.563684,8.095333,11.12402],[2.512037,8.625567,9.0394],[-9.027985,5.416164,10.57541],[4.634072,-0.1929671,4.744667],[0.237661,-4.320173,4.440763],[1.252222,-6.198136,6.401949],[3.445012,-3.257416,4.8455],[6.231672,9.99636,11.82205],[9.957458,8.451795,13.099],[6.66563,-7.513552,10.09376],[-2.48316,4.460933,5.2025],[-1.735369,6.345592,6.654175],[-9.512291,5.511035,11.0388],[-1.202129,5.488572,5.706973],[-5.842095,1.741488,6.177609],[-4.080179,-5.380664,6.826375],[-0.2171224,-7.257626,7.329412],[-5.521269,8.873476,10.49871],[8.522131,8.971009,12.41393],[1.814372,-3.893864,4.410683],[-1.166887,8.80185,8.934999],[-9.198527,-3.906793,10.0437],[2.567712,5.525129,6.174156],[-7.811181,-9.271245,12.16431],[-7.034051,1.465047,7.254257],[3.563511,6.845798,7.782259],[3.21813,3.677566,4.98807],[-2.570184,-8.177276,8.629813],[2.986188,9.331241,9.848318],[1.720468,-2.783075,3.421333],[-2.101199,-1.584903,2.815484],[-4.37086,-2.036557,4.92463],[-7.236173,-0.4550656,7.319104],[4.334809,-4.564884,6.374068],[8.28124,2.345389,8.66486],[-9.477623,8.484712,12.75992],[3.979211,-6.800369,7.942237],[6.604368,-9.990359,12.01769],[4.006103,2.186614,4.672273],[-0.7665653,-7.7161,7.8183],[6.860217,3.595786,7.809754],[-1.352559,3.3278,3.728763],[2.974854,0.2344117,3.147174],[-2.062276,-0.9550993,2.482981],[-1.539129,-1.807192,2.575823],[-6.782157,2.18609,7.195599],[8.49743,8.446326,12.02276],[3.538452,-3.436829,5.033134],[0.2872576,9.083922,9.143312],[-3.517957,-2.46332,4.409531],[5.980503,-5.553136,8.222149],[-5.996469,8.854427,10.74051],[-7.571315,0.2768982,7.642087],[-8.971998,-6.223009,10.96461],[-5.032347,-7.083431,8.746399],[-6.650301,-4.314278,7.989962],[-8.827759,-1.17574,8.961679],[1.012233,-1.338721,1.953661],[-0.4069908,-5.858973,5.957617],[-9.168438,-4.964858,10.47426],[-0.3618004,1.819035,2.10708],[-2.19035,-0.6046414,2.482584],[-4.321205,2.421154,5.053197],[8.664856,-1.635404,8.874361],[-1.486179,-0.05252942,1.792062],[-6.598129,4.658865,8.138816],[-8.626231,0.8010147,8.720865],[-4.006539,-2.399019,4.775735],[-6.460368,5.013347,8.238325],[3.115916,-4.339634,5.435196],[-2.317019,9.590232,9.916709],[5.79641,6.02563,8.420605],[4.088683,-8.206788,9.223269],[3.379387,-8.217355,8.941206],[-6.845811,-2.212047,7.26349],[9.010644,0.160674,9.067389],[-4.520114,-1.108102,4.760181],[-3.960474,-2.2276,4.652693],[-5.202273,4.951246,7.251102],[0.9739152,0.2286387,1.414492],[-2.613561,9.248532,9.662611],[-7.819681,7.10535,10.61289],[5.719466,8.959301,10.67621],[1.032918,-9.273552,9.384332],[2.673075,-7.191541,7.737157],[-5.098046,-4.597415,6.937312],[-3.266222,-5.875201,6.796042],[-7.416749,8.669814,11.45311],[1.999683,-4.949436,5.43099],[-4.360598,-1.251312,4.645492],[8.892723,5.431631,10.4682],[8.976098,3.03131,9.526761],[0.8242227,-7.296085,7.410277],[7.854018,1.31125,8.02527],[9.96048,0.761255,10.03946],[8.999065,-5.700532,10.6995],[9.961438,-5.296532,11.32623],[-2.231737,0.04797886,2.446007],[-5.31444,-0.9808491,5.495938],[-6.34016,9.736681,11.66193],[-0.1405704,-4.84863,4.952673],[-2.65626,5.21047,5.933356],[-5.752221,-5.460871,7.99432],[6.285916,-3.348216,7.191891],[1.958533,-9.935605,10.17605],[4.082737,3.940379,5.761539],[-0.1660168,-9.729791,9.782454],[7.83813,1.232879,7.997266],[4.175943,-5.310148,6.829068],[8.842064,-1.426971,9.012122],[5.377401,5.005413,7.414216],[1.427061,1.175881,2.102189],[5.354255,-4.58602,7.120367],[-7.088688,1.750088,7.369689],[-6.117498,4.855997,7.874293],[4.321116,-2.662959,5.173336],[-9.133011,-1.57402,9.32145],[-9.950322,-1.086755,10.05932],[-0.9085729,-0.4616627,1.427809],[-1.765643,-8.661813,8.896319],[2.495375,7.075105,7.56862],[4.304368,-8.708488,9.765517],[3.291614,-4.252625,5.469875],[-0.06350743,-5.132362,5.229261],[-5.756594,-5.998531,8.373813],[-3.446975,-1.202322,3.78513],[1.259475,-5.305825,5.544191],[3.248164,-7.686001,8.403879],[9.28395,-4.663032,10.43722],[6.422599,5.757902,8.683502],[9.968272,4.660773,11.0494],[-6.312305,-7.660166,9.976138],[-9.961235,8.549431,13.16507],[-8.831702,-5.049996,10.22259],[-8.098369,-6.473126,10.41561],[7.935949,8.021037,11.32768],[-1.054764,7.606422,7.744041],[7.319246,-0.6679088,7.417376],[-0.6228406,9.146373,9.221934],[8.213723,7.600577,11.23539],[4.806979,3.595523,6.085625],[5.305471,-8.058913,9.700212],[8.176721,6.692789,10.61377],[4.6984,4.270964,6.42776],[8.633535,-2.12911,8.948243],[9.910367,-7.401221,12.40941],[-5.992874,-9.716436,11.45965],[-3.516543,9.720257,10.38506],[9.964641,6.477415,11.9269],[-3.843204,-2.173472,4.527051],[-3.870059,6.384424,7.532478],[-6.396796,-5.749531,8.658875],[8.310374,-9.662182,12.78359],[0.003175242,-7.472034,7.538654],[6.388819,2.682924,7.001077],[-7.805254,-1.691582,8.048816],[2.018576,-3.048574,3.790574],[-2.626921,-9.608878,10.01156],[-1.326632,6.882082,7.07976],[8.595325,-8.313768,11.99993],[6.828905,-2.515126,7.345733],[-7.01678,-7.279289,10.15988],[-0.947649,2.422272,2.786654],[-0.1471628,3.325254,3.475482],[3.847013,4.484198,5.99229],[-7.254772,2.746905,7.821586],[-1.578596,-7.978799,8.194705],[-4.118902,-0.02146784,4.23861],[6.829428,3.29057,7.646498],[-2.414286,6.551425,7.053365],[-0.7022897,-3.593711,3.795784],[3.884145,4.122781,5.751861],[2.669055,2.522107,3.805901],[-2.391649,-6.242077,6.758957],[-8.533544,-6.480155,10.76168],[4.761231,-7.98952,9.354238],[-0.4931737,-4.955318,5.079213],[-7.708428,2.381491,8.12966],[-1.329093,-8.366303,8.530037],[7.919778,-8.208394,11.44992],[-1.685743,-0.1149085,1.963399],[0.1251292,9.953254,10.00414],[-0.4837623,-3.924105,4.078312],[3.546972,-7.121738,8.018738],[9.136032,2.180532,9.44573],[0.5201251,-8.735956,8.808374],[5.674231,8.875386,10.58156],[-4.301173,8.83085,9.873398],[-0.8669236,-9.371727,9.464714],[-1.78464,-9.987376,10.19473],[7.143572,-9.672669,12.06612],[8.567078,3.208896,9.202817],[-9.730177,-6.735032,11.8759],[3.775859,0.05804045,3.906467],[8.746965,7.849939,11.79538],[-5.368512,2.414229,5.970714],[-1.038313,8.572183,8.69255],[-1.849439,-3.773161,4.319394],[3.583935,-7.010862,7.93705],[2.345286,6.449784,6.935422],[2.945334,-2.586705,4.045496],[-9.28803,4.758542,10.48386],[8.347279,8.847899,12.20502],[-8.825763,-8.560483,12.33596],[-0.511731,9.758428,9.82287],[-3.826276,7.594543,8.562563],[-2.297189,3.976526,4.699982],[-2.379452,-6.231309,6.744702],[-4.746875,2.527446,5.469992],[-6.675619,-6.247814,9.197776],[-6.456285,6.233543,9.029987],[-1.124454,1.377658,2.040181],[-9.371878,1.302352,9.514631],[-6.300137,3.942327,7.498911],[3.36234,-5.924472,6.885107],[0.02877293,8.209517,8.270247],[1.879692,-2.314492,3.144855],[-2.019574,-9.314815,9.583551],[4.050168,9.377419,10.26352],[0.06567475,-5.29344,5.387469],[1.142947,6.533425,6.707605],[6.454525,-8.207373,10.48913],[1.978135,8.494253,8.778687],[-6.413297,6.861104,9.444846],[-3.777156,0.9050322,4.010734],[0.1367473,-1.598056,1.890101],[4.929979,1.355541,5.209816],[-5.935758,3.706128,7.068848],[-6.840974,-9.309042,11.59557],[7.55067,-1.696861,7.803329],[9.236504,2.459469,9.610514],[-4.360991,-6.862153,8.19191],[-6.749981,9.751724,11.90203],[2.887198,-1.655367,3.475076],[4.555217,1.696464,4.96266],[6.938653,4.30997,8.229261],[-0.6749889,6.712042,6.819613],[-4.206535,2.998337,5.261651],[-0.7753767,-2.29432,2.620136],[-1.097536,1.321962,1.988006],[5.450697,7.635018,9.434172],[1.204,8.813084,8.95098],[0.2611371,-3.57574,3.722111],[3.34717,5.233631,6.292411],[5.807158,-7.103186,9.229211],[-2.057605,-2.960222,3.74121],[0.9492258,8.585353,8.695361],[1.487402,5.987214,6.249727],[9.780889,-2.428483,10.12735],[-3.856296,8.626837,9.502281],[-3.303499,2.641207,4.346157],[-3.849544,-5.088033,6.458101],[4.510579,6.078043,7.634653],[4.733953,4.385271,6.530001],[4.093221,-7.403226,8.518346],[-4.765079,7.905339,9.284416],[6.881997,8.986346,11.36294],[3.546852,0.5798545,3.730468],[0.6422952,-5.587445,5.71245],[9.756817,-8.002405,12.65836],[-3.078542,-7.112895,7.814775],[-9.759953,7.928243,12.61403],[8.917852,5.695872,10.62878],[-9.141186,6.993979,11.55323],[5.649253,-2.744061,6.359554],[6.009523,-9.114164,10.96277],[-5.701095,0.0003382005,5.788133],[1.162601,-5.123074,5.347666],[-3.417685,-1.639606,3.920316],[-5.920545,-1.206758,6.124469],[-6.055091,-6.870179,9.212137],[1.59142,-1.282266,2.275263],[1.511324,0.5587375,1.896388],[2.25783,2.509166,3.52047],[-3.014926,7.432176,8.082513],[-1.363776,8.54944,8.715092],[9.053772,-3.016051,9.595173],[-0.5845368,4.075113,4.236535],[4.873462,-3.898479,6.320504],[-5.198618,-2.989377,6.079639],[6.555161,8.182425,10.53196],[6.738505,5.376369,8.678294],[-1.587449,-2.443388,3.080607],[-9.210304,-0.8920859,9.307283],[4.244548,-6.98519,8.234627],[9.634701,-0.8139758,9.720597],[6.92641,-7.002208,9.899801],[-3.480829,-6.72266,7.63612],[-8.131789,0.2806815,8.197851],[-3.936451,1.292557,4.2622],[-3.306552,1.339463,3.705057],[0.5845514,-5.390221,5.513273],[1.986386,4.623194,5.130269],[-0.3905732,-1.537562,1.875272],[-0.4433313,-8.900284,8.967251],[-9.896967,8.224089,12.9068],[-5.539051,-2.244497,6.059608],[-6.192106,-6.598693,9.104116],[7.799154,1.498992,8.00461],[0.7315909,9.500031,9.580491],[-2.69218,-9.468432,9.894395],[0.2602515,-4.05347,4.183103],[5.091085,-1.426486,5.380893],[-2.066404,5.749497,6.190859],[-4.727398,3.31195,5.858097],[5.071301,7.835826,9.387134],[6.385142,-2.269845,6.849981],[-6.931363,7.424224,10.20602],[7.994401,7.560511,11.04861],[-4.704292,-5.908688,7.618593],[4.958176,5.433969,7.423713],[1.035075,8.064593,8.192012],[9.435819,-6.943178,11.75765],[1.484686,-7.530641,7.740468],[-7.137606,1.135475,7.296213],[1.133024,1.750322,2.312438],[7.621129,6.338502,9.962842],[8.642584,-2.936575,9.182468],[-5.191862,-5.779867,7.833408],[-8.512844,-7.440871,11.35055],[3.024447,-9.49743,10.01741],[-9.909455,-5.079688,11.18036],[5.817839,-3.282792,6.754552],[-1.720743,-8.72197,8.946156],[-0.1134006,-1.58053,1.873749],[3.926838,1.501001,4.321233],[-4.81981,-5.807367,7.612889],[3.21061,-6.803831,7.589475],[4.84507,-1.634395,5.210178],[6.407008,9.788041,11.74119],[-4.412595,-0.002801558,4.524489],[-7.183946,0.01592841,7.253229],[-1.874779,7.037568,7.351337],[9.419631,8.117165,12.47469],[-1.328788,4.756385,5.038737],[-8.366064,-8.125026,11.705],[1.858195,4.742432,5.190717],[-0.2389697,-4.995475,5.100184],[5.537696,5.293663,7.725862],[-5.994616,3.943197,7.2446],[-3.669861,-2.530679,4.568613],[-7.666212,8.158738,11.23992],[9.434498,-8.508505,12.7438],[9.4807,3.053211,10.01028],[5.149949,6.711606,8.518663],[0.8383255,-6.732842,6.858131],[-3.834799,-8.891646,9.734838],[-3.210676,-2.783138,4.365123],[5.958632,-7.724224,9.806576],[6.09935,-4.630532,7.722946],[4.031431,1.79892,4.526428],[-8.418694,9.856373,13.00086],[8.432146,-9.852248,13.00646],[-3.17346,9.039262,9.632191],[7.702759,-0.1876267,7.769665],[2.725812,-5.132061,5.896448],[9.259423,1.211481,9.39173],[2.226624,0.5350823,2.498834],[6.708278,7.881425,10.39797],[-7.579759,-9.813981,12.44054],[6.090719,-0.7107639,6.213055],[9.842042,-6.904179,12.06373],[6.639154,5.115967,8.44106],[1.943455,-9.899028,10.13744],[-7.16399,9.496553,11.93764],[-6.431551,6.76273,9.386126],[-4.810117,-6.545177,8.183921],[-3.343397,5.969669,6.914858],[3.148456,-4.370961,5.478876],[8.59661,6.060768,10.56573],[-3.821444,8.991271,9.820712],[-0.2248116,0.4333252,1.112794],[-8.312072,-2.366546,8.700062],[2.515895,-0.7276482,2.803426],[-1.785973,5.690784,6.047704],[-8.623596,4.222074,9.653617],[0.8864586,-9.720638,9.812065],[-2.833818,-9.506208,9.969881],[5.229098,9.661953,11.03163],[-1.203379,-1.720033,2.325217],[5.600996,1.508218,5.886075],[1.548139,0.5629035,1.927069],[2.492042,-0.261804,2.697928],[-8.79223,1.579307,8.988744],[7.561719,9.694558,12.33548],[5.772565,-1.269017,5.994407],[-0.796079,3.052852,3.309629],[-9.654892,-5.881423,11.34937],[-6.114015,-4.458129,7.632567],[-7.078677,-2.212794,7.48359],[2.98707,2.182968,3.832484],[-4.574138,1.844025,5.032213],[-5.818757,8.280632,10.1699],[9.239558,-9.492204,13.28425],[0.1273028,4.950349,5.051946],[1.439345,0.1029962,1.755654],[-9.031464,-2.038116,9.312424],[-1.606379,-1.149148,2.213819],[3.142701,-9.549229,10.10269],[-8.505375,-4.104911,9.496931],[-7.185192,-3.811416,8.194747],[-3.714039,-2.151504,4.407159],[-6.428707,-5.342975,8.418768],[-6.815545,-5.210689,8.637299],[3.961151,7.913762,8.906085],[0.2867004,0.6015444,1.201687],[2.743747,-7.070983,7.650291],[6.172889,2.129431,6.605984],[1.86723,1.465592,2.575754],[5.668011,6.405434,8.611384],[7.047507,9.039029,11.50528],[-4.166893,-7.385161,8.538361],[-3.214375,0.7783158,3.455139],[-3.6436,4.677432,6.012836],[-8.715343,5.527618,10.36879],[2.664623,9.147576,9.580103],[4.784992,1.801154,5.209636],[-0.4562527,7.237981,7.320965],[-8.938846,-8.782934,12.57151],[5.684935,0.2493311,5.777599],[8.634586,0.6135774,8.713929],[3.724366,-5.156191,6.438727],[2.584017,8.706952,9.137185],[0.7252603,4.985926,5.136678],[9.829041,-8.184972,12.8298],[8.173757,-2.371691,8.569435],[-0.2424591,-5.357921,5.455832],[-0.637948,-7.898927,7.987492],[5.171413,-3.71479,6.4454],[-6.583494,-5.908904,8.90267],[1.222323,-8.732979,8.874626],[6.176297,4.035224,7.445111],[9.043328,-6.020391,10.90994],[4.343039,-1.314122,4.646386],[-7.406632,-0.2129447,7.476867],[-5.930724,-6.058285,8.536762],[6.479351,-5.116382,8.31621],[-2.201322,0.8423689,2.560352],[0.3672808,-3.526938,3.684316],[-4.28589,-0.1351238,4.40308],[-1.449916,-4.212266,4.565681],[-1.971353,-5.75412,6.1641],[-2.354739,-9.764953,10.09451],[-7.132049,-6.765631,9.881289],[6.375224,-2.129933,6.795594],[7.556833,4.175402,8.691358],[2.60592,-3.249162,4.283442],[5.326232,-1.714092,5.683912],[-4.390086,-5.343,6.987167],[2.037257,-6.180394,6.583896],[8.631883,3.818434,9.491567],[5.822118,-5.575985,8.123341],[-1.215822,-0.2257919,1.590347],[-6.489497,-3.996876,7.68691],[0.4339459,-9.077522,9.142741],[-5.691982,-2.500682,6.296989],[4.318553,6.510588,7.876399],[-0.9427246,1.41735,1.974236],[3.144178,-2.922055,4.407296],[-5.373308,7.819366,9.540174],[3.378784,-9.442305,10.07836],[8.230807,-1.441455,8.415698],[-1.70247,-1.889684,2.733004],[-6.675728,1.263037,6.867358],[-2.193156,0.159548,2.415655],[3.134127,-1.415264,3.581302],[-1.589675,8.230441,8.441992],[7.216843,-1.173036,7.379623],[6.115482,6.648402,9.088474],[-4.718324,-8.224687,9.534572],[9.796232,-2.069963,10.06235],[7.644711,-9.458405,12.20258],[-4.778167,-9.182371,10.39937],[3.114856,-6.393063,7.181475],[-7.412675,-8.695173,11.46969],[2.659851,6.391689,6.994891],[5.090868,5.302477,7.418436],[-8.303514,-4.045014,9.290343],[-5.781781,6.563475,8.803874],[-2.282662,8.432644,8.793181],[-8.718794,-6.687311,11.03347],[-1.797253,-6.033904,6.374803],[4.670062,-4.464409,6.537616],[3.211351,8.692034,9.320098],[-6.210711,6.632183,9.141049],[2.930666,8.983961,9.50265],[-7.116538,7.187198,10.16371],[-8.724978,4.838112,10.02659],[9.900552,9.94648,14.06959],[-6.598035,9.708784,11.78111],[-6.978902,0.5033729,7.06813],[-7.714168,-1.409159,7.905322],[-5.779727,8.725163,10.5135],[8.561168,-7.696928,11.55579],[9.05477,1.715871,9.27001],[-3.672422,-0.7436646,3.878108],[-1.860377,-1.166564,2.412856],[-9.331932,-1.991064,9.594233],[5.232635,-7.602649,9.283359],[4.402773,0.8718396,4.598317],[9.289529,7.433735,11.93967],[-8.841322,2.582144,9.264796],[-2.423716,0.3212551,2.641515],[-7.504774,6.786659,10.16761],[4.908409,-8.590297,9.944128],[5.678101,-7.476893,9.441649],[8.971113,0.361403,9.033907],[4.612993,-6.309518,7.879703],[-2.801332,-8.532589,9.03618],[1.859999,-9.764791,9.990534],[-0.1223905,-2.152766,2.376842],[2.302488,9.828669,10.14417],[0.5685661,-6.249756,6.35474],[-1.978482,7.366503,7.692838],[-1.869934,8.447145,8.709243],[3.526625,-7.297998,8.166876],[8.810927,7.269131,11.46616],[-3.277145,0.2397907,3.434702],[-4.59306,-4.094045,6.233571],[-5.124409,-8.153001,9.681478],[9.226227,-6.7167,11.45588],[8.329406,-3.852328,9.231437],[9.345326,-0.8438118,9.43648],[-6.216916,2.564122,6.79888],[0.4501454,-7.007302,7.092596],[6.901176,8.498323,10.99308],[6.968523,7.471736,10.26582],[-5.121204,7.937611,9.499073],[-5.862463,4.403567,7.399992],[1.943803,4.026189,4.581328],[-4.539941,9.038058,10.16354],[6.092238,-2.613747,6.704256],[-0.5523768,6.722203,6.818587],[-0.8807149,6.290572,6.43016],[5.17356,5.816909,7.848703],[-3.61358,0.146569,3.752259],[-5.954983,-4.063736,7.278446],[2.460107,-6.288762,6.826468],[2.765574,4.440061,5.325649],[-9.213046,0.4470236,9.277933],[0.7065157,3.178584,3.406253],[4.651875,3.240661,5.756895],[-5.737328,9.083543,10.79017],[-5.973075,-4.505481,7.54831],[3.403608,6.682033,7.565323],[5.225115,2.009542,5.686834],[-7.29562,6.935678,10.11582],[8.341672,-7.670549,11.37633],[-3.155721,-5.422273,6.352922],[2.368886,1.363965,2.910674],[-3.64349,8.323347,9.140739],[-4.535423,-1.0259,4.756315],[1.713539,-8.055218,8.295948],[6.425728,-5.428614,8.471118],[-1.126094,8.788054,8.916163],[-6.74098,-5.139996,8.535829],[5.632917,0.866414,5.786227],[1.026796,-1.406616,2.008203],[0.5926866,7.976244,8.060506],[9.796929,-2.180968,10.08645],[-7.692049,-8.227604,11.30757],[8.352587,3.062636,8.952399],[8.817795,3.112678,9.404375],[3.367173,1.638844,3.876038],[6.772399,-5.985216,9.093306],[9.480343,8.715034,12.91622],[1.259418,-3.196919,3.578607],[-7.784199,-3.230407,8.487007],[3.548122,-0.1848098,3.690979],[5.666448,0.4495046,5.771541],[6.273116,-2.600337,6.863944],[-3.247134,2.639091,4.302172],[0.6157043,-0.9285924,1.497122],[-0.397537,4.462512,4.590429],[-2.363607,7.24121,7.682563],[-5.024988,-8.88337,10.25499],[0.02714735,-6.807304,6.880416],[5.724142,2.73471,6.422183],[-1.556867,0.1652712,1.857727],[2.437657,5.699063,6.278654],[-1.453669,8.222518,8.409694],[7.891091,7.27129,10.77687],[4.973039,3.644826,6.246268],[-2.148363,-1.290172,2.698149],[-9.824419,-0.8524811,9.911909],[-5.071531,0.1047136,5.170241],[8.824073,-2.311228,9.176385],[-9.950633,-2.528361,10.31541],[-1.570319,8.983089,9.173974],[9.715645,-4.918021,10.9353],[-0.4915651,4.475212,4.611849],[-9.966254,-3.040321,10.46756],[9.072865,-3.916114,9.932412],[6.96347,5.874078,9.164863],[-6.153059,-5.597968,8.378387],[6.768522,-1.299626,6.964332],[9.994459,-5.592236,11.49619],[3.145407,0.855683,3.40966],[3.645443,-6.468386,7.491947],[7.3457,-5.992078,9.532277],[-9.367178,-0.6231905,9.440995],[4.726673,3.897429,6.207366],[2.685268,-2.809078,4.012678],[9.635968,1.522252,9.806586],[0.2299369,5.514014,5.608674],[-9.454101,-4.84646,10.6709],[0.4608471,-5.578248,5.685879],[-1.186673,-5.862148,6.064073],[-3.431759,2.963291,4.643066],[-6.378344,0.6249974,6.486439],[-2.741036,5.544598,6.265449],[-6.133494,0.3199404,6.22271],[7.957261,6.265694,10.17727],[5.798728,3.169422,6.683598],[-3.768661,7.400171,8.364529],[3.408745,1.011114,3.693494],[-1.139079,5.473227,5.679235],[-6.1916,6.327544,8.909193],[-0.7524023,0.6898185,1.428971],[-0.2387759,-3.716099,3.855698],[3.492849,4.148553,5.51457],[5.868623,-4.040486,7.194878],[-4.964653,2.2621,5.546609],[-2.064678,-1.370818,2.672459],[4.545645,-1.863695,5.013607],[1.614517,2.21834,2.920222],[-1.98102,-9.795766,10.04398],[3.999891,-2.120296,4.636247],[-0.2504907,-6.307937,6.39162],[0.5715886,-3.276291,3.472866],[-2.232058,2.856742,3.760726],[-2.052534,4.344109,4.907563],[3.162111,-2.354582,4.067309],[-4.224827,-4.067536,5.949287],[-5.236049,-8.583732,10.10429],[-1.377895,0.6234703,1.813094],[-0.03061666,-9.565313,9.617492],[-2.097532,8.389739,8.705593],[-3.39277,8.748192,9.436193],[-1.579357,-0.8387722,2.04888],[-4.559032,6.773366,8.225769],[6.824577,4.727161,8.361872],[7.636345,7.220038,10.55664],[-7.5246,-7.620345,10.75589],[3.663876,3.856132,5.412369],[-9.530212,-7.413086,12.11523],[-2.85618,0.7955877,3.129013],[-7.926993,-1.767873,8.183068],[1.395448,-1.604666,2.349942],[-8.58397,5.97853,10.50844],[8.204805,-8.482701,11.84378],[-2.528067,0.1748075,2.724277],[0.1843074,0.5097941,1.137479],[3.592715,5.228517,6.422226],[-0.8022833,-3.228786,3.474006],[-3.161693,4.916305,5.930122],[6.057115,8.054401,10.12729],[4.605575,8.162604,9.425467],[-8.472075,-2.825166,8.986525],[7.976295,2.750408,8.496236],[7.995648,-4.052494,9.019595],[-8.750119,5.161572,10.20815],[0.8302956,7.378136,7.491748],[-5.487214,-2.759067,6.222698],[-2.988625,1.270672,3.398012],[6.740742,6.875594,9.680465],[3.909247,5.975696,7.210489],[5.770918,-0.7618713,5.906263],[-0.6957083,6.174659,6.293682],[6.895833,5.249075,8.723835],[-0.3838115,-0.6928713,1.275689],[-0.7511274,-6.397105,6.518216],[3.012791,-1.551948,3.533476],[4.677028,0.09618674,4.783706],[-8.901204,-1.24233,9.042943],[-3.030371,-7.915993,8.534992],[-0.9613995,-1.37038,1.949931],[4.934332,1.557121,5.269939],[-8.174003,5.386798,9.840321],[8.881634,4.569489,10.03811],[-7.979372,6.954469,10.63179],[5.003873,-0.1734396,5.105764],[3.199261,-0.9488475,3.483616],[2.167207,1.98123,3.101944],[7.853066,0.9784937,7.976722],[3.572113,-8.587242,9.354182],[-8.462676,9.99012,13.13086],[0.4816113,2.502291,2.737409],[-2.530841,5.836735,6.439925],[7.056058,-9.389963,11.7881],[-1.731091,-9.9827,10.18091],[-1.212852,3.021595,3.406031],[2.751721,-3.254798,4.377862],[-2.199744,0.0502214,2.416898],[5.794235,5.050643,7.751268],[4.580724,4.234765,6.317933],[9.077262,5.069532,10.44494],[-7.234488,7.747513,10.64715],[-1.554882,-7.907702,8.120925],[-5.018642,0.1946119,5.121],[-1.631149,6.671094,6.940039],[6.545647,-3.618124,7.545616],[-2.776318,-8.624308,9.115187],[-0.5863408,9.134809,9.20807],[0.919656,4.392591,4.597893],[9.871752,5.913122,11.55061],[2.794753,-5.047592,5.855667],[-3.159183,-7.935737,8.599789],[-2.079922,1.286963,2.642414],[-4.510835,1.42228,4.834306],[-8.268965,-1.657295,8.492491],[-9.368095,-4.447832,10.41847],[-9.941953,6.770451,12.06986],[6.451815,-0.3535357,6.538418],[3.131323,-5.428794,6.346416],[-9.36201,-4.089815,10.26517],[-4.94766,5.702868,7.615907],[2.407736,9.415155,9.76946],[1.078312,4.063225,4.321175],[-9.115693,8.124866,12.25191],[-5.166039,3.129997,6.122486],[0.446642,-4.485253,4.617032],[-2.786937,8.314138,8.82564],[7.313338,0.1356409,7.382636],[-5.648239,4.250359,7.139199],[-8.738784,9.357939,12.84279],[7.146266,6.306078,9.583096],[0.220617,-9.172583,9.229568],[-0.03352026,-2.108375,2.333745],[2.691955,3.511127,4.535927],[9.384059,1.785999,9.604705],[-6.156694,5.077703,8.042882],[6.981015,-4.238972,8.22821],[4.585124,-0.07844605,4.693562],[-9.805724,-5.060467,11.07974],[-3.262731,-3.058494,4.582553],[5.667136,-5.269591,7.802885],[-8.675536,-6.925065,11.14547],[1.756378,3.281664,3.854112],[2.769966,-4.278598,5.194142],[-3.536299,-0.0502138,3.675313],[-1.234987,1.445216,2.147986],[9.114923,-1.239418,9.252998],[-8.814472,7.052705,11.33294],[6.703195,8.07312,10.54078],[6.624823,0.5549707,6.722817],[-7.23088,1.530312,7.458383],[-3.773198,-3.300552,5.111816],[-8.484586,7.244505,11.20139],[0.9218155,0.1220218,1.365516],[-8.632879,7.231766,11.30597],[-3.8937,-1.567511,4.314857],[0.7757208,1.317376,1.826807],[8.974139,4.182344,9.95124],[8.732759,-9.264638,12.77085],[-8.449502,1.61081,8.659607],[8.78466,0.5525085,8.858641],[-9.39983,2.131582,9.690225],[-7.868934,-8.226571,11.42789],[-1.735859,7.818471,8.07104],[2.627481,-9.197844,9.617899],[-6.687459,-4.552818,8.151703],[-1.675011,-1.921468,2.738193],[3.055664,6.630907,7.369261],[-8.150497,5.363699,9.808153],[-3.491322,-8.374535,9.128098],[-0.4435148,-3.081098,3.269536],[4.929679,-4.024912,6.442178],[6.137059,2.033488,6.542061],[5.574124,0.6252104,5.697521],[5.001685,8.240156,9.69108],[-2.001446,7.894809,8.205717],[7.789741,-2.567505,8.262696],[6.353422,-4.066763,7.609502],[-2.142015,-4.900218,5.440622],[-1.289749,6.903225,7.093516],[5.031083,-5.002819,7.165193],[2.564011,7.359951,7.857674],[0.06994602,5.516962,5.607295],[5.13666,-7.235562,8.929648],[-6.05153,5.860044,8.482991],[-1.300622,-3.144835,3.547055],[4.870419,-1.725187,5.262818],[-3.288009,0.04176957,3.436968],[-6.159395,9.229239,11.14078],[-5.711587,8.442889,10.24229],[4.84835,-6.177685,7.916457],[-3.634527,-6.022688,7.105108],[-2.845341,0.3525958,3.036493],[6.906506,-7.791764,10.45999],[1.120914,3.787187,4.074216],[1.04476,5.511952,5.698521],[7.180264,1.073161,7.328565],[5.877617,-3.672277,7.002285],[2.860896,-2.604263,3.995862],[-5.713141,-6.098545,8.416189],[-7.830185,0.6976722,7.924553],[-7.618977,-4.151505,8.73406],[5.226352,3.264827,6.242904],[1.715426,-0.2553042,2.001966],[8.872216,-1.597548,9.070192],[-3.067755,-4.950715,5.909374],[-1.941772,6.735886,7.081147],[-0.2626123,-9.659324,9.714499],[6.035718,6.456255,8.894556],[-4.963308,4.437322,6.732329],[5.801674,1.333385,6.036335],[-0.5648732,0.9085909,1.464452],[2.564757,-0.1222894,2.755528],[5.314218,-5.021285,7.379309],[-7.895869,-6.990839,10.59323],[-6.660478,5.588873,8.751998],[-9.261439,-8.30958,12.48292],[-0.1152086,6.36035,6.439513],[-7.834631,0.8368785,7.942406],[-1.781197,-5.205369,5.591826],[-9.801784,5.088519,11.08909],[3.198986,-9.589025,10.1579],[-1.204233,3.041464,3.420626],[7.677082,0.6278552,7.767354],[8.025846,9.753731,12.67081],[0.4239451,-0.4594648,1.179338],[-0.3165215,-5.900624,5.993125],[4.667445,-2.119202,5.222649],[5.186003,1.029087,5.380859],[-7.693999,-6.290645,9.988485],[-4.583825,2.670511,5.398433],[0.3795631,7.32406,7.401752],[-4.913333,5.655811,7.558376],[-0.1229653,-9.061055,9.116899],[-7.420374,-6.911229,10.18956],[-2.103309,-6.048307,6.481198],[-3.908605,2.018751,4.51138],[6.865092,-1.542077,7.106863],[6.098489,-9.54079,11.36742],[0.02303421,-2.53566,2.725821],[6.656744,-6.327019,9.23815],[3.894964,-8.596059,9.49015],[-8.904092,-1.74284,9.127997],[-5.448986,-5.452226,7.772915],[-1.082812,-1.581103,2.161566],[9.925782,-7.682561,12.59138],[-0.02044907,1.927344,2.171422],[7.741354,6.244067,9.995846],[-0.5604116,-1.233433,1.68387],[-4.893585,-0.6596279,5.038083],[8.263647,-8.941373,12.21622],[-9.371228,6.87865,11.66772],[-2.887908,5.712253,6.478414],[8.395977,-3.263701,9.063343],[5.558252,-6.340277,8.490776],[-0.361561,8.190516,8.259254],[-9.340508,1.933191,9.590742],[9.625546,2.643094,10.0318],[-2.610256,-8.511839,8.959065],[-2.067486,-1.576757,2.785797],[9.0611,-7.923067,12.07802],[-8.36918,-0.262283,8.432792],[-8.729993,3.661694,9.519494],[-9.963963,-7.091597,12.27075],[0.9667211,-0.2322506,1.410138],[-7.892538,4.037686,8.921606],[9.423474,-3.793584,10.2075],[5.5897,-8.841227,10.50771],[2.02856,9.640434,9.902173],[1.408355,6.508661,6.733953],[-4.168839,-3.392972,5.46731],[0.3661405,4.988121,5.10053],[-0.3243754,1.854155,2.131458],[2.787305,4.820713,5.657592],[8.691643,1.541047,8.883663],[-9.599483,-7.879524,12.45941],[-1.928683,-1.132439,2.449946],[-7.620471,1.319894,7.798314],[1.760347,2.004556,2.849046],[-2.977897,-5.412057,6.257654],[-7.436267,3.551436,8.30125],[-6.587712,-3.878249,7.709654],[-3.114157,5.741285,6.607597],[5.770242,4.143603,7.173921],[4.90867,-4.12878,6.491676],[0.1061024,-5.758353,5.845502],[-1.15713,-0.8682531,1.75864],[8.455318,6.024734,10.43024],[2.87105,3.220447,4.428793],[-1.491341,-5.536266,5.820167],[0.05917628,3.15267,3.307995],[-2.342754,-8.000742,8.396449],[0.6526024,-2.251445,2.548508],[-0.5315793,-2.08757,2.374979],[-9.483057,-1.670519,9.680858],[2.560334,-3.862617,4.740794],[9.344026,-3.675591,10.09063],[-1.620749,-7.15423,7.403367],[-7.54264,9.813334,12.41745],[0.6581783,-0.999714,1.559688],[-3.36078,-4.642632,5.817978],[-1.547312,-6.890494,7.132537],[8.59678,-6.860675,11.04416],[-4.328981,4.72014,6.482268],[-7.059129,-8.230525,10.88912],[-3.707352,7.228424,8.18502],[2.223587,-9.408253,9.71903],[0.9695383,6.595372,6.740841],[9.972059,1.442284,10.12532],[-4.45276,3.444577,5.717708],[0.9263209,-6.223936,6.371456],[-4.770729,2.517374,5.486075],[-4.860175,8.948103,10.23181],[-7.119716,9.231068,11.70055],[-8.271671,-1.72569,8.508734],[4.312801,0.5269353,4.458465],[6.273713,4.736818,7.924451],[-0.7954014,6.875453,6.993177],[2.944426,-2.759865,4.157703],[-1.625123,5.270832,5.605594],[2.010818,-9.987833,10.2372],[-6.566801,8.392813,10.70337],[5.338665,-3.76864,6.610899],[-3.707789,-5.991301,7.116417],[6.727559,9.525448,11.70445],[5.886319,0.4115959,5.984828],[-2.896024,5.736445,6.503365],[6.719172,-8.248734,10.68592],[2.655728,6.225114,6.841413],[-2.520917,8.344558,8.774204],[-8.155748,8.656614,11.93538],[4.750352,-1.773141,5.168159],[7.103421,0.8974991,7.229391],[1.579468,7.742548,7.965034],[9.391182,-8.091545,12.43653],[-4.504276,4.089758,6.1656],[0.8510399,2.62388,2.934113],[7.882614,-1.804811,8.148187],[0.2948385,-8.459833,8.523831],[1.859791,-5.400042,5.798213],[-0.07765443,6.769044,6.842952],[-4.361794,-1.410521,4.691995],[0.8175639,0.626614,1.435638],[-6.932276,7.058878,9.944054],[5.409409,9.917595,11.34109],[1.595567,9.875236,10.05316],[-4.085325,-4.574531,6.214195],[-3.074847,1.709825,3.657621],[0.6487843,5.463596,5.59212],[1.315113,-5.677998,5.913475],[5.064538,1.347727,5.335345],[-4.432533,-3.711971,5.867373],[3.958381,-6.54045,7.71014],[2.567716,-9.019423,9.430968],[-8.280606,2.08127,8.596518],[6.030274,3.853254,7.225771],[5.519918,7.67179,9.503991],[9.794312,-8.575317,13.05621],[1.572587,-4.751002,5.103435],[-9.367412,-2.469072,9.738826],[-5.567148,-4.551611,7.260186],[-5.529717,6.835659,8.848955],[-2.598328,-9.113197,9.528991],[-9.057766,-2.931741,9.572785],[2.310779,9.310882,9.645322],[6.568303,-1.19531,6.750658],[-4.979407,-4.265796,6.632609],[-6.056858,6.028311,8.60384],[-3.534643,-9.023995,9.743008],[7.727399,-9.153692,12.02093],[-3.039054,6.900986,7.60654],[-0.6871678,-4.013977,4.193353],[-2.730926,-0.6901178,2.989017],[-1.276487,-2.072218,2.631256],[-7.862317,1.740177,8.114447],[-0.04788458,-4.815825,4.918787],[-3.218704,-8.877201,9.495512],[-5.738143,-2.598269,6.377875],[4.946555,-8.913199,10.24273],[6.840295,-9.689208,11.90254],[1.514121,-4.667672,5.007966],[9.331351,6.844305,11.61545],[-6.68945,-8.129333,10.5752],[4.622447,-9.699542,10.79111],[8.073965,-5.566976,9.857999],[7.422898,-0.08794449,7.49047],[1.59819,-7.253628,7.49462],[8.495255,-5.438328,10.13631],[0.9893582,-5.891506,6.057117],[-1.302723,-3.134779,3.538916],[-2.2684,1.672571,2.990508],[7.293083,-7.491862,10.50319],[4.559189,3.938592,6.107267],[-7.476426,-4.836211,8.960239],[0.07257926,3.489533,3.630717],[2.578427,8.671405,9.101733],[-1.700708,0.2090265,1.983961],[-3.563574,-1.275161,3.914728],[8.83339,0.5629784,8.907621],[0.8228799,1.545611,2.016444],[6.160541,0.4300254,6.255972],[2.652571,4.527829,5.342038],[-4.908752,7.580643,9.086363],[-0.2475922,0.7442766,1.270925],[8.98802,3.47213,9.687115],[-7.484377,-3.622884,8.375034],[-0.9712468,-3.349067,3.627612],[-5.128825,7.334453,9.005501],[-6.478066,7.813004,10.19845],[2.763827,-0.02792751,2.939306],[8.919458,-6.318466,10.97633],[-1.582768,8.81579,9.012398],[0.5572333,-6.405056,6.506555],[9.146947,4.755019,10.35745],[-2.271315,1.142491,2.732061],[4.451296,-9.975582,10.96933],[-5.91888,4.544013,7.528691],[8.205305,-2.244459,8.565315],[9.367803,-9.146473,13.13064],[0.9444498,6.05405,6.208342],[8.708511,-1.896828,8.968618],[-6.577792,5.737506,8.785575],[-3.668125,3.065928,4.884163],[-6.853799,8.109402,10.66475],[-9.179476,4.13652,10.11798],[2.987286,-7.978736,8.578118],[-1.871645,-4.631244,5.094259],[-8.128564,2.601316,8.593043],[-2.342243,-5.08891,5.690616],[-2.13654,7.546615,7.906718],[-8.120636,5.304657,9.751108],[7.735485,1.75302,7.994423],[2.810771,5.446122,6.209725],[8.534062,-0.2220986,8.595322],[-5.088343,-6.111478,8.015073],[2.031003,7.844086,8.16423],[-7.822056,2.019354,8.140168],[-2.871026,2.710015,4.072711],[-5.293822,3.317883,6.327156],[-4.53763,-2.406813,5.232861],[-0.7144206,8.494853,8.583294],[-0.6928185,-9.255297,9.334908],[-9.031346,6.358394,11.09028],[-4.548903,8.057827,9.307045],[-7.572444,4.761558,9.000797],[-3.907696,2.849155,4.938397],[6.811929,-1.571467,7.062003],[8.026845,7.65428,11.13635],[-5.014406,1.80473,5.422298],[-4.452302,-3.895666,5.999934],[-9.745532,-6.849249,11.95356],[-8.252456,-8.149432,11.64115],[-1.287291,-1.569174,2.262614],[2.035094,-6.387805,6.778323],[-7.998936,6.158072,10.1442],[-5.582253,-6.905292,8.935581],[-4.205225,-0.7847086,4.39314],[-4.977543,1.90756,5.423534],[6.833866,7.364063,10.0961],[8.840868,9.035715,12.68089],[8.119843,2.896391,8.678763],[1.024766,-6.154596,6.318955],[7.523291,2.491617,7.987995],[9.884971,-0.7096169,9.960733],[-4.473674,-9.859607,10.87316],[-1.211977,4.070703,4.363429],[-7.130369,-0.9394756,7.261183],[-8.756721,1.229472,8.898975],[2.015134,0.2345587,2.261809],[5.792131,-0.09246592,5.878549],[-8.942019,5.242732,10.41374],[-5.283311,-6.989999,8.818926],[-3.068269,-5.448164,6.3322],[9.422103,-5.866987,11.1444],[2.217899,6.027002,6.499526],[1.475196,-0.2644929,1.80171],[0.5483621,-0.569298,1.274677],[-6.723847,3.127157,7.482596],[-3.220378,-6.362167,7.200556],[3.885213,1.779958,4.388979],[8.008543,9.819411,12.71053],[-5.400459,7.780539,9.523746],[-3.80077,-5.28319,6.584675],[1.515672,4.31562,4.682076],[0.04704114,-0.4570099,1.100487],[-9.3807,4.534233,10.46694],[-4.317955,2.627459,5.152502],[-5.386005,2.18635,5.898235],[-3.078362,1.657046,3.636222],[-1.55085,3.22317,3.714022],[9.329842,1.534254,9.507886],[9.618481,3.020506,10.13107],[-2.672339,-8.200042,8.682285],[-4.10764,-3.578767,5.538978],[-6.229263,5.48944,8.362874],[-2.213726,-1.306887,2.758357],[-6.686526,-8.625918,10.95975],[-0.8037808,7.062741,7.178327],[7.366685,7.079023,10.26551],[2.67157,-9.112498,9.548555],[0.4371137,-6.308893,6.402593],[-7.375981,-8.662299,11.42106],[4.227077,-3.703878,5.708493],[-4.415784,7.307591,8.596513],[-9.612047,5.260453,11.0029],[0.9410959,2.172843,2.570391],[-9.349108,1.818529,9.576683],[7.367763,7.844082,10.80803],[-5.048435,-9.108312,10.46174],[-8.680584,-4.560285,9.856405],[-1.689643,-7.92039,8.160114],[1.941117,6.196791,6.570248],[3.106349,-2.68605,4.226614],[-4.227938,0.587585,4.384143],[-5.182187,-6.347508,8.255054],[-4.037848,-9.009819,9.923762],[9.963615,-1.0582,10.06943],[-5.58919,-2.017205,6.025626],[-4.398933,9.518809,10.53368],[-8.203061,2.292809,8.575965],[-2.950682,-1.014693,3.276603],[-8.859926,-0.6380008,8.938978],[5.805517,-5.400095,7.991561],[-0.4993446,-8.088068,8.164936],[-7.288155,3.507306,8.149748],[-0.3446994,-9.071862,9.133318],[-4.813023,-4.186533,6.456953],[6.631924,-8.602888,10.90835],[6.793623,-5.574646,8.844771],[-6.295321,5.316723,8.300518],[3.581292,3.971518,5.44046],[-5.594051,-9.496297,11.06675],[-7.024911,-8.619894,11.16476],[-4.089489,-0.3160424,4.221825],[0.9346508,-5.175265,5.353218],[1.88957,-6.630456,6.966593],[4.958542,7.824816,9.317451],[7.695049,8.875834,11.78958],[-5.708388,-9.720697,11.31714],[-7.658229,-5.881218,9.707584],[-3.470033,5.236098,6.360649],[3.423062,-7.583859,8.380469],[-8.505283,-6.043657,10.48168],[9.874221,-6.004384,11.59969],[-4.314255,2.801472,5.240328],[-6.125445,2.006508,6.522818],[-4.701215,8.121949,9.437557],[-9.440939,-5.981909,11.22117],[6.973144,-6.553458,9.621463],[1.958218,-6.354204,6.723877],[-2.080184,2.418175,3.342863],[1.132305,9.261869,9.38426],[0.3114181,-5.076701,5.183616],[5.000511,-4.823949,7.019657],[3.958211,3.3041,5.252096],[4.586076,-4.007062,6.1716],[3.646121,-3.674597,5.272273],[-9.04536,9.094843,12.86603],[3.210931,7.584123,8.296325],[-6.204061,2.103299,6.626782],[-8.139106,-3.34698,8.857049],[0.1923764,-5.829826,5.918098],[-0.5821831,8.26426,8.344875],[6.312841,9.974486,11.84662],[8.134214,2.507849,8.570575],[4.805101,1.150662,5.041132],[-6.5411,7.301963,9.854169],[-9.252667,-1.949382,9.508519],[1.740101,-0.2142676,2.018381],[0.9483807,-8.587557,8.697446],[-1.740935,6.754721,7.04678],[-0.3919378,-5.768196,5.867342],[-1.119636,-1.340234,2.012414],[-4.225517,-5.248886,6.812181],[0.05004711,1.070789,1.465979],[0.2205344,8.511806,8.573183],[2.833916,-3.228281,4.410542],[-1.900849,4.026293,4.563361],[2.40012,4.037052,4.801913],[2.084653,-0.8789852,2.473539],[-0.3727897,-1.107105,1.537743],[7.574677,0.6735095,7.670029],[-2.082255,6.426575,6.829103],[8.700549,-7.504076,11.53303],[0.3297594,3.711161,3.857649],[3.57172,-8.592532,9.358889],[9.486297,-7.634013,12.21753],[-1.455543,6.328349,6.57013],[6.775012,-8.803459,11.15355],[1.237194,-6.189594,6.390753],[-0.8071821,6.896319,7.015038],[-3.590311,-6.186021,7.221994],[-9.25328,2.32729,9.593719],[-2.66138,0.4080593,2.872187],[3.661843,-3.306708,5.034224],[-5.030243,-8.949168,10.3146],[-4.717939,7.70677,9.091384],[-0.800963,-5.171043,5.327404],[-1.514609,7.327822,7.549239],[-5.392818,-3.050046,6.275768],[-9.834284,1.345805,9.976188],[-2.257185,5.88502,6.381876],[8.582619,-7.583958,11.49686],[-3.778784,8.526208,9.379521],[-6.605089,-0.4704065,6.696901],[-1.53807,2.799272,3.346876],[0.3782699,6.548408,6.635114],[7.344037,-9.033677,11.68513],[-5.475315,-1.467997,5.756222],[-8.443909,-5.244874,9.99041],[-9.097788,3.67546,9.862999],[-8.605014,3.831026,9.472224],[6.22117,-3.094599,7.019936],[-7.454628,-2.541923,7.939323],[-0.1863018,0.7802246,1.281975],[-6.003625,-7.085527,9.340675],[8.657167,-7.444147,11.46132],[9.15439,-3.104486,9.71806],[-5.208347,7.554413,9.230169],[4.407073,1.329419,4.710588],[-7.283471,0.2506438,7.35607],[-6.674237,1.936529,7.021081],[0.513122,-9.802138,9.866367],[-2.481669,8.504194,8.915155],[8.048215,2.029599,8.360206],[5.976627,4.436513,7.510175],[0.1913482,0.165432,1.031495],[1.736068,-3.002151,3.609272],[-3.159579,-1.220147,3.531529],[-3.618634,3.80265,5.343656],[-3.839973,7.819675,8.768848],[4.49885,-8.24042,9.441619],[1.557806,-9.492188,9.671007],[-6.694396,3.114524,7.450852],[-5.840158,4.604533,7.503944],[-9.250802,-3.492455,9.93854],[3.199217,5.167014,6.158979],[-2.865358,-6.597273,7.261837],[0.07649726,5.792635,5.878816],[-8.201728,6.090166,10.26443],[5.913901,-2.800275,6.619348],[-5.040035,1.4115,5.328629],[7.689444,-4.93036,9.188906],[2.89882,2.154769,3.747824],[-1.376317,2.320995,2.877719],[1.097776,-9.577212,9.691651],[-4.198784,5.162857,6.729404],[-5.333177,-3.173281,6.285896],[0.4731708,-1.678766,2.010509],[-3.467787,-2.215908,4.235067],[-3.398185,2.353945,4.253083],[9.247824,9.027206,12.96197],[5.12629,7.297447,8.973939],[1.827306,1.292038,2.451205],[6.19278,-3.053747,6.976811],[0.6183995,2.701298,2.946087],[-9.857362,-5.289218,11.23136],[-9.698234,-9.671102,13.73266],[-6.664303,2.340588,7.133813],[-3.959518,8.326156,9.273761],[6.376746,-3.131662,7.174273],[3.058139,-0.3306349,3.23443],[-3.477865,-9.830137,10.47507],[-5.777637,3.487473,6.822284],[-2.272362,3.205062,4.05414],[-4.502324,8.11825,9.336857],[0.8084341,4.404276,4.588161],[-3.20337,1.255677,3.583058],[6.671663,3.509382,7.604397],[-2.843734,-7.924769,8.478725],[-1.746291,6.936928,7.222915],[-5.173424,1.445329,5.463817],[-9.246421,-2.147739,9.545107],[6.034188,-0.393842,6.129155],[8.285094,6.476793,10.56369],[-9.945242,3.268205,10.51613],[6.602489,-9.642103,11.72873],[1.081629,-3.75785,4.036256],[1.593793,-4.964276,5.30888],[-4.380206,0.3940746,4.510155],[-1.460837,-6.576234,6.810353],[-9.279702,-4.899477,10.54124],[-9.552313,4.663843,10.67699],[4.920699,-3.801898,6.29823],[1.330257,-3.425948,3.808767],[9.245344,-0.8129395,9.334734],[7.989747,-5.170611,9.569288],[-5.38701,-4.040065,6.807496],[8.188102,4.037775,9.184152],[9.458018,0.1028924,9.511293],[5.518744,7.790508,9.599403],[-0.5603711,0.2490118,1.17304],[-7.244575,0.4245221,7.325577],[5.485299,-9.575283,11.08037],[8.648648,4.587071,9.840749],[9.384186,1.277118,9.523338],[-3.406436,8.834037,9.520716],[-6.227909,-8.416624,10.51791],[-4.587909,1.575971,4.953038],[7.546268,6.464905,9.987049],[5.326536,1.533739,5.632437],[-5.521974,-0.4231442,5.627721],[-7.884479,-9.739418,12.57065],[3.065085,6.434003,7.196606],[-6.585783,7.794355,10.25302],[-2.756704,-6.531453,7.15956],[-8.055808,-5.821294,9.98917],[-0.4318567,-9.627213,9.688639],[5.503541,-3.025461,6.359432],[7.399652,-5.180482,9.088027],[1.502369,-8.386772,8.578756],[8.741837,9.959373,13.28942],[4.989625,-6.40066,8.17709],[9.907209,-5.484608,11.3681],[-1.398442,5.695308,5.949132],[-5.847764,4.550245,7.476702],[3.574004,2.360997,4.398614],[-7.103703,8.621171,11.21549],[1.33427,-4.727232,5.012683],[6.317934,-8.280746,10.46361],[3.974149,4.737729,6.264178],[5.433543,-8.299885,9.97053],[7.789483,9.804729,12.56219],[6.897159,9.86858,12.08138],[1.368209,8.553389,8.719659],[9.591915,2.331086,9.921633],[-9.112831,4.995437,10.44022],[-8.81927,-4.187082,9.813826],[-3.899167,9.232754,10.0721],[-8.131642,8.784856,12.01238],[-8.234715,-5.553664,9.982671],[-3.63055,-8.058772,8.895206],[-0.5004715,0.3757097,1.179673],[1.288012,-9.501261,9.640173],[-7.020773,-9.883094,12.16416],[-6.672575,-0.2375668,6.751274],[7.073048,-6.087205,9.385205],[-8.147676,1.985496,8.445521],[-2.34111,-0.4088805,2.578367],[-8.655838,6.977238,11.16268],[4.115394,7.500473,8.613568],[-8.026224,-4.965421,9.490821],[3.412475,5.097459,6.215229],[3.463235,0.4871455,3.637486],[4.070653,-9.131362,10.04749],[3.569018,2.724166,4.599888],[3.777491,1.124452,4.066181],[-8.363262,-3.870089,9.269398],[-8.030842,-5.095548,9.563421],[-5.180842,-5.806614,7.845884],[4.726517,-4.58236,6.658678],[2.861328,5.767391,6.515367],[-6.186208,3.37994,7.119913],[-6.204551,3.676137,7.280827],[2.654691,2.364103,3.692745],[0.00846419,-2.976504,3.140008],[7.124101,-3.932832,8.198779],[5.507616,8.216025,9.941675],[-3.408413,5.998132,6.971002],[-0.7528598,2.766899,3.036862],[2.378039,-3.539226,4.379634],[1.215792,0.5692243,1.673967],[-6.472194,5.067191,8.280442],[6.871382,5.323272,8.749464],[3.833009,4.241039,5.803307],[-1.657754,4.503138,4.901673],[-4.142649,9.351177,10.27648],[4.245843,-9.610056,10.55369],[-7.776928,1.081365,7.915174],[0.6313806,-8.498402,8.580296],[7.914124,5.182645,9.51279],[1.999021,-1.984016,2.988713],[-8.608974,-8.834377,12.37581],[-0.1842621,-9.104429,9.161036],[-4.463251,-2.947528,5.441372],[2.617846,6.323191,6.916348],[-4.388184,8.522636,9.638023],[2.353544,-0.8377615,2.690914],[-1.778016,-7.021984,7.312291],[9.007226,-7.36319,11.67676],[4.927439,4.919424,7.03423],[-5.102855,4.182783,6.67344],[8.886512,-2.715543,9.345816],[-6.661032,9.249573,11.4422],[-3.094651,-9.657151,10.19007],[-3.241611,-1.4479,3.688422],[0.8944147,2.710047,3.02396],[2.345039,4.716369,5.361282],[-0.6578109,-2.479692,2.753468],[6.721481,1.631987,6.988683],[4.470157,-3.212457,5.594835],[8.277605,-2.32863,8.656862],[-4.952847,3.876831,6.368713],[-0.1358875,7.167522,7.238221],[0.3642903,8.061838,8.131785],[1.641988,2.535964,3.182332],[2.860137,-2.187816,3.737234],[5.73823,-7.24674,9.297447],[-6.782026,-6.292163,9.305224],[5.100291,-4.968967,7.190522],[-1.883052,6.345166,6.693804],[0.1683278,0.08346415,1.017497],[-0.2377991,-8.33717,8.400294],[5.627065,6.29052,8.499088],[9.353195,-9.818936,13.59757],[-1.621771,-3.198608,3.723068],[7.192931,0.9689987,7.326474],[8.63508,5.086285,10.07149],[2.102231,8.65371,8.961366],[9.894201,-2.238774,10.19349],[8.190392,-7.405055,11.08681],[-6.345327,4.307199,7.734025],[4.208673,-8.746495,9.757771],[6.62853,-5.426552,8.624667],[9.473272,-9.005304,13.10871],[-6.394259,3.020326,7.142053],[-6.431546,-4.909522,8.152802],[-9.964492,-8.711798,13.27353],[4.491734,-7.062324,8.42924],[-7.067148,-3.036106,7.75645],[-5.554258,3.732126,6.765985],[-8.709225,3.220101,9.339146],[9.855475,7.142904,12.21276],[-8.672705,6.125474,10.66477],[-0.3241166,-8.519677,8.584285],[-8.106476,2.352405,8.499928],[5.463415,-4.317992,7.035194],[-1.520319,-8.007737,8.211895],[3.646886,-7.204578,8.13669],[4.52784,-3.519991,5.821656],[3.401963,9.963077,10.57527],[9.547862,4.640257,10.66272],[9.855782,8.786025,13.24125],[3.43933,9.114237,9.792768],[-4.008897,-7.795266,8.822552],[-9.374822,-7.859648,12.27442],[-8.837395,1.359964,8.997169],[4.081876,-5.138642,6.638325],[9.103245,3.191905,9.698316],[4.696304,-1.960799,5.186522],[-4.296092,-1.033506,4.530401],[-7.506926,1.670453,7.755279],[-3.653022,-9.183108,9.933481],[4.386205,8.122967,9.285547],[7.104578,-5.202581,8.862385],[-6.335342,1.511901,6.589568],[8.844131,3.67276,9.628489],[-0.0836802,6.65635,6.731567],[-4.216298,0.9566904,4.437614],[0.9624963,-9.949294,10.04564],[1.696162,-5.460142,5.804318],[-6.894196,2.252523,7.321461],[-1.098641,-4.780403,5.005923],[-0.92011,7.827037,7.944125],[-3.273946,8.528016,9.189439],[-0.5921013,-4.662289,4.804947],[-3.734113,-3.244865,5.047054],[-2.821681,8.847816,9.340543],[5.973403,1.718313,6.295565],[-5.413475,5.122678,7.51981],[0.2748265,-9.278308,9.336087],[-2.416136,-8.074759,8.487606],[-2.438126,-9.655097,10.00826],[-0.860733,-3.591719,3.826396],[6.814395,3.712067,7.824028],[-0.1512861,9.715086,9.767588],[2.079297,9.722048,9.992081],[0.1602747,6.15919,6.241899],[7.150445,-4.351374,8.429907],[-4.912046,7.253247,8.816903],[3.313258,-0.4364766,3.488294],[-8.007559,-5.401042,9.710421],[-8.773735,5.343706,10.32151],[1.812785,9.764929,9.981985],[-6.330198,6.353443,9.024282],[4.117129,-8.961232,9.912337],[4.135772,2.781761,5.083582],[-3.709543,8.093154,8.958786],[-6.903333,-9.40776,11.71162],[4.920719,3.360218,6.041899],[5.517754,-6.235108,8.385832],[7.483998,7.84889,10.89106],[-1.422322,5.242317,5.523123],[-1.297191,3.041083,3.454112],[-0.6932096,5.432197,5.566803],[-2.660855,-2.715499,3.931168],[-6.526903,-4.02987,7.735652],[-8.994722,-2.477776,9.383198],[1.344755,1.114255,2.012444],[-8.310623,8.257821,11.75832],[4.012216,5.260782,6.691316],[1.860358,4.230139,4.728109],[-5.680151,6.141972,8.425433],[-4.631479,-6.192647,7.797402],[8.768841,4.147613,9.75168],[7.27931,-7.327649,10.37703],[7.475484,-4.87677,8.981411],[2.175358,6.071966,6.52694],[9.343933,-8.546568,12.70248],[5.521414,-4.566196,7.234373],[2.110999,-0.1277481,2.339367],[6.133238,-0.6194916,6.245028],[7.636235,1.936895,7.941262],[-8.193365,-6.21764,10.33394],[-4.197696,9.970716,10.86443],[-3.293712,-6.870378,7.684441],[6.501111,-2.273584,6.959427],[-5.169088,8.631076,10.11014],[0.6799436,-2.828709,3.076349],[-4.625515,9.848477,10.92648],[-4.039587,9.976832,10.80997],[-2.381156,0.7002983,2.675878],[-4.320964,2.103901,4.908884],[-8.828771,-7.127019,11.39042],[-8.424398,6.392712,10.62249],[-9.468686,1.256209,9.603857],[8.409105,-5.301296,9.990835],[9.542516,7.661332,12.27826],[-9.376445,9.593121,13.4516],[-9.22235,-4.24204,10.20033],[7.007424,9.567676,11.90145],[-0.5157887,6.579267,6.674788],[-8.498514,7.320752,11.26136],[-1.186122,-4.447841,4.710644],[-6.154821,6.352644,8.901568],[-2.483309,0.5262555,2.728327],[-6.133219,3.440567,7.103089],[-1.850801,9.342095,9.576023],[-1.910659,8.499621,8.768932],[8.674284,8.792012,12.39123],[1.21802,1.568901,2.223741],[2.443402,-6.278926,6.811397],[-1.049228,-4.955646,5.163265],[0.4471389,-7.240188,7.322586],[-8.004623,-8.935256,12.03797],[5.096827,8.700178,10.13266],[-3.356565,7.372504,8.162128],[6.5272,-2.342362,7.006497],[7.527561,-2.266067,7.924597],[8.803021,5.676377,10.52209],[-3.289887,-8.364706,9.043875],[-7.199174,-5.939231,9.386297],[5.4246,5.061669,7.48644],[-3.763376,5.353172,6.619627],[-7.779036,-1.10065,7.919901],[-9.037905,2.002374,9.31092],[-6.209653,4.973032,8.018157],[-6.868597,3.27031,7.672845],[-1.084399,-8.011011,8.145687],[-0.19871,8.661882,8.72168],[-7.474057,-6.421443,9.904365],[-5.824407,-1.944304,6.221257],[-4.237559,9.614012,10.55396],[2.97298,6.17842,6.929032],[-5.288781,-4.640368,7.106632],[0.09748384,9.271699,9.325979],[-2.733221,0.342806,2.930531],[3.660259,-2.800046,4.715693],[2.85258,-7.243559,7.848972],[-4.827399,4.741879,6.840263],[-8.312787,-1.768331,8.55742],[3.003031,-4.953661,5.878517],[-5.238241,3.793713,6.544572],[8.819159,0.5777015,8.894453],[-2.218428,2.262863,3.322946],[-9.891909,2.924544,10.36353],[-7.605673,5.091452,9.207016],[-8.677004,3.638836,9.46211],[-7.400281,-3.116132,8.091627],[-9.340113,3.056729,9.878325],[0.7392623,9.973129,10.05036],[5.011908,-6.588912,8.338643],[-6.113444,-5.701231,8.418921],[-1.37468,-9.605685,9.754945],[-5.036319,7.219014,8.858819],[7.015351,1.579307,7.260121],[5.337308,3.877641,6.672553],[6.078102,0.1742336,6.162279],[1.857093,3.106973,3.755273],[9.646472,-8.18922,12.69322],[-7.806624,7.608266,10.94665],[3.953963,-1.293243,4.278586],[4.878132,-8.193668,9.588138],[-4.241272,-8.245793,9.326387],[4.583147,-6.013002,7.626364],[-5.163519,-1.259569,5.408183],[5.424809,8.846239,10.42519],[6.77919,-9.802999,11.96061],[-5.816376,-8.00568,9.94591],[3.37724,0.2355142,3.530045],[-1.272229,9.507941,9.644662],[-7.036113,-4.112913,8.211147],[2.528866,-7.632853,8.102814],[-5.978845,-7.470517,9.620562],[-2.004041,4.917244,5.403283],[-6.527682,2.639809,7.111907],[-6.273365,4.752095,7.933318],[-1.736331,-3.436656,3.978121],[9.446387,7.11002,11.86535],[2.355834,-7.306991,7.742226],[-7.407826,-9.341113,11.96379],[-4.105353,-3.655481,5.58717],[7.408773,3.199331,8.131766],[-4.096444,5.550698,6.970732],[1.229261,-6.558821,6.747534],[-1.404189,5.2798,5.554101],[-7.562,4.633255,8.924735],[-7.32924,4.482044,8.649074],[4.327305,-3.061116,5.394071],[5.433433,9.411162,10.91294],[3.35107,7.270765,8.068067],[1.658204,7.71437,7.953687],[2.350625,-7.385228,7.814539],[6.32415,-5.017839,8.134715],[5.729139,-9.89349,11.47624],[-9.578288,-0.5986705,9.648938],[-4.395489,7.866179,9.066261],[8.887831,8.516375,12.34999],[-6.822834,-0.6655062,6.927767],[-8.697207,1.200575,8.836447],[8.101689,0.754006,8.19792],[-5.055212,3.602833,6.287732],[7.599844,9.451752,12.16936],[-1.059859,-3.5209,3.810517],[1.470723,4.609729,4.940913],[-7.235186,9.021031,11.60719],[-8.737629,2.930336,9.270007],[8.50043,6.655081,10.84193],[-9.985733,5.154525,11.28202],[-2.324172,7.092688,7.530471],[2.539665,7.234468,7.732234],[5.943471,2.960564,6.714893],[-4.998644,-2.594707,5.720047],[-8.557516,4.090441,9.537442],[-5.356759,-4.27804,6.92795],[3.506959,-5.752959,6.811409],[2.303844,3.468144,4.282023],[6.788889,4.353317,8.126523],[-5.470558,4.186804,6.961058],[-9.826497,-3.168422,10.37299],[-0.6345074,-0.6013285,1.32823],[-1.662757,1.622342,2.529181],[3.746644,-0.4937233,3.909105],[-6.482711,-0.7381079,6.600784],[-0.9160557,3.810504,4.044639],[5.727883,8.777951,10.52906],[-2.557914,-4.135497,4.964399],[-7.242374,-6.882777,10.04115],[6.706203,-5.216043,8.554546],[-5.296378,-8.201425,9.814019],[1.614036,6.729176,6.991918],[4.228331,-9.966308,10.87226],[2.584314,0.7999793,2.884206],[1.098211,-5.781582,5.969318],[-7.94208,6.718913,10.45086],[9.699252,-7.688447,12.41723],[-2.472052,8.95604,9.344608],[-2.566027,-5.532864,6.180378],[-7.843612,9.764907,12.56486],[-5.953922,-8.764696,10.64279],[-7.635367,2.173394,8.001404],[5.715504,-5.938817,8.302803],[3.022935,-7.281115,7.946872],[-2.195314,-5.32821,5.848866],[-6.696933,-4.284548,8.012881],[-6.22815,-5.397714,8.302118],[2.953192,7.704,8.311014],[3.170845,-6.042594,6.896898],[4.82937,1.535101,5.165206],[-6.632435,5.815989,8.877777],[3.752564,-8.694406,9.522312],[6.510014,-6.787874,9.458093],[8.186159,-6.36848,10.41973],[-9.778513,1.346002,9.921242],[-3.487041,-7.269691,8.124522],[6.533079,9.02444,11.18578],[-3.277369,6.754363,7.573808],[-1.585491,-1.478712,2.387545],[3.667328,7.29227,8.223534],[-3.456783,7.962462,8.737858],[-5.026861,-3.042894,5.960582],[-8.215939,-7.917788,11.45395],[-9.427879,-8.75769,12.90667],[5.036247,7.801427,9.339489],[2.685361,-6.427624,7.037437],[-3.024619,5.015242,5.941462],[8.214381,-1.504491,8.410681],[2.52175,-3.294428,4.267608],[-1.028031,-8.491248,8.611512],[0.4200528,7.577144,7.654382],[-4.506796,-3.651003,5.885664],[3.507754,-5.383912,6.503141],[1.14214,-7.265862,7.422751],[-6.23141,-7.945045,10.14664],[-8.28712,-2.834662,8.815421],[-2.693166,-7.707397,8.225395],[-4.678989,-6.301057,7.91178],[-4.621791,-6.385323,7.945647],[9.883827,8.125211,12.8339],[0.175563,-3.808017,3.941043],[9.229416,-2.438194,9.598276],[2.489372,-0.8037568,2.800535],[-0.8205624,-5.489109,5.639471],[9.596234,0.1865652,9.650001],[2.065513,9.06792,9.353797],[-5.297056,0.699661,5.435837],[-6.003742,6.078227,8.60173],[-7.575471,-9.330177,12.05985],[-8.079981,0.5243261,8.158493],[-4.66063,1.377823,4.961841],[-1.742158,-5.782901,6.121851],[-1.7212,-6.557142,6.852638],[-2.866317,-4.316349,5.276991],[-2.947071,5.304501,6.150037],[6.001,-9.7129,11.46091],[0.3601148,6.074141,6.16643],[9.771941,-4.692154,10.8861],[-5.913303,7.668077,9.734811],[-4.647981,-3.48841,5.896841],[-0.8620473,2.128571,2.504783],[-7.634768,1.733547,7.89271],[1.300079,-2.671776,3.135058],[6.810694,0.3985775,6.895247],[2.860226,-0.9004071,3.160953],[-9.223289,-0.3708122,9.284749],[2.855296,-4.329863,5.282085],[-3.104126,5.606508,6.486026],[-6.328829,-5.23476,8.273862],[5.052829,4.11543,6.593014],[7.551884,0.6950386,7.649446],[-0.3361078,-5.098318,5.206325],[6.510005,-7.51157,9.990189],[-0.5843338,3.991308,4.155958],[8.556939,4.188618,9.579442],[5.527631,4.296954,7.072377],[-6.362906,5.712953,8.609553],[-0.04065396,3.533752,3.672745],[-6.249732,-2.306056,6.736248],[-1.376093,6.638126,6.852617],[4.829661,-5.226271,7.186066],[2.653857,4.368217,5.208097],[-8.712836,5.090327,10.14026],[5.254259,8.438421,9.990705],[-7.454636,-8.122212,11.06987],[9.444942,-6.833362,11.7005],[0.4633803,-0.1555532,1.113067],[7.822089,4.155527,8.913669],[-4.341193,7.507366,8.729633],[7.153741,3.069317,7.848358],[-5.025961,-0.900012,5.202914],[-7.083825,1.885503,7.398357],[2.71598,-7.053759,7.624438],[1.125293,-8.533052,8.664829],[-4.097926,3.829888,5.697459],[5.153508,-5.663762,7.722489],[9.027788,-6.824181,11.36092],[6.511296,-6.288052,9.106952],[5.097913,-2.640175,5.827456],[-8.679426,9.772799,13.10878],[-8.267003,-9.416638,12.57046],[1.516385,-9.359751,9.534378],[2.678054,3.352383,4.405728],[3.856893,-5.029574,6.41656],[-0.9105905,-1.371008,1.925834],[-4.811333,-5.347414,7.262491],[-5.365844,9.008847,10.53336],[6.744812,-3.000938,7.449706],[0.9295557,0.7049063,1.536544],[8.848479,-9.129352,12.75306],[5.762689,4.940704,7.656314],[-6.756377,-8.254164,10.71354],[3.949971,-6.261746,7.470725],[-5.73892,8.180806,10.04295],[-4.480615,-6.12644,7.655663],[-1.640548,1.789111,2.625322],[7.8557,1.276933,8.021382],[-2.967768,8.063034,8.649865],[8.881578,3.32091,9.534719],[9.012138,0.5138011,9.081995],[2.733267,-9.012237,9.470541],[-7.119005,4.453341,8.456506],[8.670476,-9.234056,12.7061],[-8.626374,-8.543777,12.18238],[-0.1266658,4.551849,4.662122],[-2.712544,0.9101889,3.030897],[1.982368,8.471878,8.757997],[-9.130811,-2.352025,9.481758],[9.435747,6.764711,11.6531],[7.465111,-4.810232,8.93679],[2.978899,7.785901,8.396075],[2.308105,-3.50225,4.311973],[-7.500843,6.320781,9.859763],[4.012511,4.153818,5.861267],[-5.885612,-4.282483,7.347115],[-4.847686,-4.755383,6.863944],[-1.242047,-7.948472,8.106842],[4.30695,-1.006376,4.534602],[6.221229,-3.092966,7.019268],[-0.9588438,-4.719742,4.918877],[4.396189,-9.739716,10.73259],[2.474442,-3.967891,4.781947],[7.926678,0.4565772,8.002542],[-7.582646,6.028792,9.738729],[8.742036,-5.678376,10.47221],[-0.8497955,0.1298555,1.318717],[-4.61554,-9.15615,10.30234],[-0.04802005,-8.787118,8.843966],[-4.769739,7.729351,9.137465],[3.526048,-5.351829,6.486531],[3.774125,-3.670161,5.358554],[-6.791216,-9.304497,11.56262],[-0.3528547,-8.474853,8.540939],[-3.200744,-4.202963,5.37677],[-5.111928,1.781026,5.504894],[8.154201,7.637369,11.21697],[-4.762461,3.893652,6.2323],[7.594437,-6.039321,9.754428],[1.734547,1.476875,2.487933],[-6.351333,-8.567682,10.71189],[1.016546,-4.496115,4.716823],[-0.8780234,-4.051511,4.264465],[7.209317,0.09051945,7.278904],[-1.768653,8.024791,8.278007],[7.548349,-0.7616227,7.652297],[9.37438,7.219262,11.8742],[4.994645,-0.8333626,5.161489],[3.048025,6.721196,7.447478],[-6.640135,8.003376,10.44727],[0.4618702,-6.983134,7.069476],[8.019113,9.899201,12.7789],[6.478812,4.825059,8.139791],[-6.137762,-7.88561,10.04266],[-6.404296,9.114817,11.18458],[0.3082063,-5.21245,5.316449],[-6.811054,-6.385899,9.389896],[-2.667733,2.832849,4.01769],[-6.17162,1.795247,6.504752],[3.297254,4.605373,5.751638],[-4.450625,7.798906,9.034987],[6.384295,-3.07068,7.154599],[-1.736159,-5.807827,6.143703],[-9.513812,1.13933,9.633831],[8.790201,1.701384,9.009015],[-1.848478,-3.955659,4.479298],[2.771293,-0.07399832,2.947124],[-4.046422,-5.522366,6.918819],[-6.003879,6.486986,8.895367],[-0.3133254,9.241992,9.301214],[-5.444848,-9.289721,10.81412],[5.033078,3.130363,6.010911],[-2.52609,8.729855,9.142838],[-8.615582,1.57312,8.814929],[6.994019,-7.619315,10.39087],[4.720222,8.544467,9.812665],[-4.65344,8.673052,9.893247],[3.47357,-9.514304,10.1778],[7.65667,7.613498,10.84389],[8.7818,0.7629068,8.871417],[5.211605,8.286202,9.839816],[0.5741524,7.707094,7.792878],[3.041765,4.235093,5.30927],[9.733126,-8.778148,13.14495],[-3.172102,-6.62955,7.417086],[-3.585435,-1.200684,3.911136],[9.319433,2.984244,9.836541],[-3.978042,-0.3739226,4.118815],[6.70028,4.199738,7.970668],[8.481762,3.883112,9.381836],[0.9457899,-3.957003,4.189557],[-5.914163,-9.869108,11.54888],[9.938349,-5.185156,11.25418],[-2.505458,0.0580698,2.698276],[-1.921153,-4.55277,5.041681],[-9.590449,2.764656,10.03095],[-0.2582407,-4.309671,4.431699],[-7.280403,-0.9532853,7.410332],[-0.518421,-4.985213,5.110881],[4.585272,-6.970812,8.403388],[7.681074,-1.959543,7.989913],[8.328671,0.9565075,8.442846],[-1.704406,-7.396946,7.656358],[-8.901085,-3.698011,9.690439],[-5.670837,-1.433594,5.934103],[-8.260428,2.595154,8.716048],[-9.404058,-8.026395,12.40401],[-0.07616752,5.155213,5.251859],[6.418629,-8.617049,10.79131],[-1.321802,9.798922,9.93811],[3.892201,4.120944,5.755989],[8.069643,7.748674,11.23215],[-4.156975,4.420312,6.149765],[7.140532,4.984724,8.765539],[3.324951,7.085734,7.890687],[-1.543996,-8.189936,8.393984],[8.899469,-2.217864,9.226022],[-8.608215,6.265348,10.69374],[-9.483435,4.900776,10.72162],[-1.598838,4.653937,5.021495],[2.767783,7.695614,8.23912],[-4.943755,0.7661147,5.10173],[7.747374,-7.632473,10.92138],[-8.817636,2.353732,9.181001],[7.324375,-9.891635,12.34872],[0.9977325,5.799835,5.969385],[-4.369578,-1.798315,4.829819],[9.340625,-8.486835,12.65992],[7.087702,3.186739,7.83523],[-4.889749,5.190933,7.201072],[-9.815722,-8.026232,12.71884],[-0.225533,3.598889,3.742041],[5.235047,-0.5369509,5.356681],[-9.27564,-7.786441,12.1518],[-9.706514,-8.311066,12.81757],[2.013272,5.192687,5.65838],[1.072282,-7.981505,8.115061],[4.311873,-0.05827542,4.426697],[8.598612,-0.3647695,8.664248],[8.673545,-2.605576,9.111499],[-0.7689646,9.789647,9.870587],[-0.178725,-2.133846,2.363312],[-4.351927,5.546595,7.120673],[-4.274554,0.6669027,4.440335],[8.089087,-8.6682,11.89836],[4.806416,-3.858615,6.244241],[-9.11662,-9.429894,13.1543],[-7.916379,-9.074356,12.08358],[-0.1070311,1.895266,2.145574],[6.279397,9.520947,11.44899],[7.370306,-8.617706,11.3836],[2.786037,6.506184,7.147897],[6.540406,-6.873806,9.540761],[-0.1309329,-5.147186,5.245061],[6.774303,9.086212,11.37763],[0.5188863,-7.642477,7.725069],[-2.906109,-4.373397,5.345285],[0.8935785,-1.982793,2.393731],[6.086888,3.611365,7.147878],[4.562513,-2.79034,5.44082],[1.430404,2.75789,3.263742],[4.350205,0.7609032,4.528052],[-2.373484,9.629734,9.968209],[1.343691,2.337634,2.875767],[-6.779666,0.09667151,6.853701],[-2.513795,2.378929,3.602564],[-5.542003,-0.02929618,5.631577],[7.728151,3.468508,8.529646],[0.008541658,-9.269371,9.32316],[1.68286,0.962463,2.181365],[0.4229776,-9.827638,9.887435],[-5.29739,-3.713336,6.546083],[6.972241,3.466967,7.850605],[-3.050101,0.7290826,3.291607],[3.931128,6.104918,7.329651],[9.961454,4.77858,11.09348],[6.509902,-2.89937,7.196192],[-5.904479,4.440134,7.455043],[9.843227,8.289136,12.90732],[4.748306,1.980513,5.241073],[8.919592,-9.707728,13.22116],[3.76708,-2.384317,4.569011],[-4.564688,0.02764288,4.673023],[1.20863,-0.5613021,1.666087],[-7.018301,5.095172,8.730254],[-4.429016,-5.080137,6.813514],[9.094856,2.965081,9.618114],[-1.921833,-8.718095,8.98324],[-4.996882,-6.135042,7.975436],[-8.384532,5.8051,10.24693],[8.089399,8.134583,11.51563],[2.501767,2.164517,3.456006],[7.609864,4.317207,8.806151],[5.714136,2.287381,6.235661],[4.434064,-5.597634,7.210717],[-8.737966,7.970977,11.86965],[1.774915,9.445719,9.662915],[-2.916885,6.122565,6.855218],[4.253862,-2.637668,5.104178],[5.357363,7.898254,9.596028],[3.397299,0.6603549,3.602459],[5.529535,7.687595,9.522336],[8.930492,-7.229106,11.53316],[-4.3671,-4.140095,6.10016],[5.508986,-3.62624,6.670722],[2.380994,-6.316705,6.824214],[7.988701,0.8453712,8.095306],[-4.116508,3.902705,5.759926],[0.7929236,8.930547,9.021275],[2.837278,7.881638,8.436253],[2.297065,-2.520831,3.554026],[-0.1972847,-8.534089,8.594743],[5.619083,8.636233,10.35174],[-8.418441,7.475679,11.30292],[7.322587,9.53005,12.05994],[8.554487,-7.799382,11.61936],[-8.7601,-3.913902,9.646657],[4.757094,0.6563044,4.905169],[2.412698,-8.433375,8.828529],[0.06849057,-5.372314,5.465021],[1.325732,-7.086471,7.278437],[-2.758644,-6.454037,7.089761],[-0.738654,9.836713,9.914965],[-9.615586,-4.923341,10.84891],[-2.425974,2.019445,3.311119],[9.72805,0.7636763,9.809086],[-7.754606,3.921351,8.747051],[3.712618,-6.628274,7.662737],[-5.12988,5.487108,7.577864],[8.336836,8.632369,12.04245],[-7.926429,6.572209,10.34515],[-8.406075,9.532824,12.74899],[-3.080389,4.831837,5.816824],[0.386309,6.745194,6.829852],[-8.867058,6.639728,11.12253],[-8.535024,-9.468293,12.78652],[3.586074,0.2255294,3.729717],[-1.308353,9.477373,9.619375],[-2.239104,5.420539,5.94944],[2.722982,-7.1358,7.702874],[6.04682,-7.668456,9.816784],[-0.09225792,9.006059,9.061877],[-9.720293,-3.254287,10.29925],[-8.613657,-9.097463,12.56817],[-0.5259464,-8.394044,8.469746],[-6.943306,2.036621,7.30461],[-7.239037,7.385971,10.3902],[-2.24103,-9.768764,10.07229],[-9.445622,4.854628,10.66711],[0.6811336,0.3635109,1.263362],[6.12316,-2.148799,6.565853],[4.218693,-2.479451,4.994502],[-6.700991,-0.748849,6.816455],[-6.710789,4.876767,8.35569],[-3.499238,3.743524,5.220981],[-4.232322,7.146689,8.365866],[-0.8796306,-4.681028,4.866804],[5.61755,-7.425023,9.364179],[9.018975,4.259015,10.02403],[-1.573974,-6.32147,6.59078],[-9.562936,-2.469244,9.92708],[6.490421,-8.973615,11.11986],[-4.667984,-1.712872,5.071884],[4.201763,-8.720414,9.731415],[-0.5117332,-0.1817929,1.137945],[-2.281942,-1.269995,2.796453],[-4.285253,5.141473,6.767432],[-8.945993,7.867944,11.95556],[9.914631,7.015563,12.1868],[0.8510534,-1.829189,2.251716],[-9.816579,3.302569,10.40539],[9.536884,-6.549827,11.6126],[-9.640514,8.866539,13.13602],[3.529039,-5.030427,6.225698],[-1.025225,0.5439321,1.531975],[1.80765,3.104873,3.729321],[4.133254,9.614061,10.51256],[-4.882473,-9.331009,10.57858],[-6.652891,-7.755595,10.26695],[6.480993,-6.496277,9.230649],[-0.496778,4.433576,4.572022],[3.992181,-9.018338,9.913018],[-1.836475,-5.866832,6.228351],[6.925438,-1.582324,7.173942],[-0.9155752,-6.797914,6.931804],[-3.774951,-2.70752,4.751939],[2.781934,6.791641,7.407128],[-3.183555,0.6631463,3.402173],[8.736612,-9.68168,13.07912],[-2.718114,-0.3984241,2.923506],[-4.051464,0.4496708,4.197209],[6.415598,1.271961,6.616478],[-2.380747,7.2704,7.715353],[5.449648,9.015038,10.58157],[7.686534,1.748236,7.946014],[9.042533,8.943469,12.75747],[-6.532913,-3.223767,7.353341],[-7.083174,-4.068359,8.229392],[-5.9045,2.622423,6.537601],[0.4420698,-3.385568,3.557738],[-0.08360122,2.324774,2.532106],[6.982682,-6.520504,9.605978],[-0.5552057,9.227955,9.29857],[-5.692018,-0.3693993,5.790986],[-1.630938,4.236986,4.648871],[8.466703,3.677413,9.284849],[-7.725689,-6.379553,10.06901],[6.721283,-2.046502,7.096747],[3.708242,-9.61476,10.35349],[9.647709,-9.304895,13.44096],[-7.42188,-9.826145,12.35465],[5.797779,-9.084334,10.82309],[0.3312819,1.217724,1.610155],[7.198078,-1.491186,7.418622],[-5.623948,0.182703,5.715083],[-9.472306,-5.447819,10.97285],[6.241331,7.408358,9.73848],[4.733553,-4.11825,6.353464],[0.5441359,6.243438,6.346385],[-1.268018,3.850943,4.175839],[-8.260001,9.783139,12.8428],[2.35613,2.974897,3.924456],[-8.311852,-7.454811,11.20987],[-6.795971,7.802593,10.39546],[-8.157111,0.1687097,8.21991],[-4.371194,-0.6546242,4.531652],[-7.272781,-1.971629,7.601359],[5.472095,-1.802168,5.847361],[-0.9563451,2.883364,3.198184],[-0.2231359,-1.119438,1.517541],[-1.809216,-8.575486,8.821124],[-5.635319,-2.525741,6.255892],[8.175455,1.087373,8.307855],[-9.798319,-8.028528,12.70686],[-8.469098,1.430329,8.647049],[8.250371,-9.245125,12.43145],[5.656221,-5.795274,8.159537],[1.629996,-5.637427,5.952939],[7.583842,3.821485,8.55093],[-2.549341,1.930063,3.350266],[-4.422685,8.96384,10.04543],[-9.218658,-6.26207,11.18915],[-5.019064,-1.304018,5.281238],[1.921373,7.531694,7.83697],[9.759215,6.267434,11.64143],[0.5014624,-7.000504,7.089324],[4.988159,3.023966,5.918285],[-8.875,9.796667,13.25671],[-0.3532771,2.080671,2.335379],[-0.6190489,-0.7729082,1.407341],[-2.888308,-5.814084,6.568553],[-1.645858,5.028585,5.384748],[7.00354,0.7855525,7.118052],[-0.3237876,-4.532771,4.653047],[-5.292446,-4.137267,6.791683],[1.91639,2.640853,3.412719],[8.231642,8.906877,12.16932],[-7.397657,-8.296646,11.16063],[9.817116,6.358974,11.73935],[2.039426,-9.406781,9.677127],[-9.398036,1.444491,9.560839],[-6.738308,1.102808,6.900796],[2.209942,-5.23971,5.773942],[-0.2415967,-0.2630748,1.061874],[6.92971,-4.775293,8.474922],[8.826113,-6.960929,11.28516],[-4.835144,6.964148,8.53686],[-1.221025,5.507308,5.728992],[2.929848,-6.134605,6.871491],[8.161365,6.524289,10.49639],[4.945938,-0.765834,5.103803],[-5.212814,4.67788,7.075027],[4.317467,-9.252325,10.25895],[-3.126554,1.746284,3.718177],[2.14129,-4.961522,5.495618],[-6.597948,6.896719,9.596752],[0.822346,9.968937,10.05266],[-6.406776,-9.870483,11.80988],[2.099644,-2.945622,3.753025],[6.547447,-3.071651,7.300966],[-3.87085,-9.452687,10.26337],[8.872138,3.275737,9.510272],[1.920628,6.6042,6.950127],[3.088315,2.746783,4.252353],[-8.494967,-3.473881,9.232135],[-6.246355,8.406932,10.5211],[3.931929,-4.315178,5.922907],[-3.165942,7.030084,7.774655],[7.217946,7.315255,10.32529],[3.589253,-0.2991308,3.737943],[-8.570781,-2.231697,8.912842],[-1.274755,3.903455,4.226342],[-5.964087,4.391073,7.473411],[-7.845904,-0.7299004,7.942983],[1.024719,-7.396264,7.533576],[-7.454608,-0.2362343,7.525091],[-1.767584,-8.497818,8.73712],[3.434939,-5.737662,6.761625],[-3.414565,-3.706123,5.137568],[1.418871,-1.683507,2.418138],[1.751529,5.164093,5.543979],[-6.761735,-6.287154,9.287054],[-2.597889,4.734739,5.49243],[-5.509061,-1.511925,5.799626],[9.839247,-4.603515,10.90886],[-7.809974,-8.466151,11.56164],[-1.109558,-3.221675,3.5511],[-1.394479,4.363064,4.68838],[-7.972812,-7.28617,10.84684],[2.33563,-4.706736,5.348694],[8.471682,-0.9432539,8.582488],[-9.764753,9.03308,13.33968],[7.617754,0.1650439,7.684882],[-2.816867,1.282319,3.25255],[6.717492,1.944764,7.064475],[-5.94455,-3.437358,6.939244],[-2.209268,3.631077,4.366416],[7.674986,2.102244,8.020277],[-0.7320896,1.629951,2.047607],[5.566275,2.773414,6.298829],[3.752674,-1.642559,4.216701],[6.355454,-6.441626,9.104194],[-2.703317,2.753822,3.986409],[-8.461588,7.608008,11.4228],[-4.949049,-5.252151,7.285477],[-3.267402,2.830409,4.437018],[-1.765477,4.828528,5.237518],[2.320202,8.377325,8.750023],[-5.393049,-7.835573,9.564579],[0.7215912,7.913606,8.00911],[2.868868,9.791408,10.25193],[-2.165875,-0.6681592,2.477388],[0.280549,-7.765381,7.834529],[5.521624,5.142292,7.611274],[7.441223,-7.0288,10.28474],[1.396245,6.806418,7.019746],[-7.24221,-1.53345,7.470012],[-6.317887,1.530392,6.577066],[2.256967,-3.612664,4.375528],[-6.605617,2.942278,7.30008],[1.368516,2.078359,2.681868],[0.222781,-1.880353,2.141345],[3.418133,-3.938334,5.309813],[-9.742201,-1.973329,9.99022],[-7.380557,-0.007511387,7.447998],[-2.510121,3.229632,4.210847],[1.645018,-2.025418,2.794352],[-0.9620118,6.062205,6.218987],[-7.196429,7.901275,10.734],[-1.658927,9.484046,9.679832],[9.359571,-6.749017,11.58235],[-6.174525,-5.415334,8.273488],[5.876275,-4.347985,7.378047],[-4.443119,-3.842043,5.958406],[8.806245,4.577332,9.975065],[-9.998322,-8.189841,12.96302],[3.670466,6.390255,7.436913],[6.010252,-7.739885,9.850328],[-1.013413,-6.738624,6.887384],[1.47961,9.349658,9.518684],[-6.357447,0.3907436,6.447465],[-3.478887,-4.81395,6.02302],[7.11959,-7.525304,10.40763],[-7.686966,2.503084,8.14585],[0.8672729,1.379579,1.91191],[3.966667,-4.116413,5.803388],[-8.752355,-6.099651,10.71492],[9.705988,-4.489369,10.74061],[9.345048,-9.938581,13.67865],[6.325536,5.847777,8.672307],[-9.467117,-2.886945,9.947903],[0.550096,8.114392,8.194264],[-1.707835,7.55361,7.808567],[8.699325,2.255898,9.042528],[-0.443489,5.047714,5.164891],[-6.014341,-6.877272,9.190711],[-1.868284,-3.825175,4.372922],[-5.479156,-7.823905,9.603887],[3.293312,-7.579367,8.324224],[5.176101,8.492621,9.995831],[-1.088759,0.9390795,1.751361],[-0.4976019,7.390312,7.474244],[-4.986975,4.434894,6.7482],[-7.340069,-3.008807,7.995595],[-6.200695,5.339779,8.243898],[1.887079,4.416033,4.905346],[-9.834669,8.284545,12.89785],[8.00197,-0.4207813,8.075183],[-0.9788389,-8.19805,8.31662],[6.242282,-3.260337,7.113079],[-7.630072,4.633099,8.982405],[9.341702,6.547118,11.45129],[-8.921689,-3.456297,9.619903],[-9.970397,-9.732785,13.96911],[3.204634,-0.7150405,3.43234],[8.882448,-5.76064,10.63404],[-1.608783,-9.725355,9.908113],[-4.105629,8.825674,9.785127],[-9.784104,7.76367,12.53009],[6.189972,5.452042,8.309062],[5.159957,4.275333,6.775222],[5.335212,-5.345921,7.618619],[-2.142636,5.338643,5.838835],[-8.125752,1.054561,8.254693],[-7.319963,6.186022,9.635803],[2.862478,2.712537,4.068371],[-3.785566,3.611607,5.326745],[-6.445819,8.58559,10.78244],[0.7414798,4.000434,4.189662],[5.835517,8.341513,10.22908],[-3.802618,-9.4677,10.25169],[5.53199,-4.757628,7.364641],[-5.671235,-6.614383,8.770003],[2.979932,0.09383439,3.144646],[8.995069,7.873106,11.99571],[-4.044803,-4.981507,6.494293],[2.591845,-1.869556,3.348567],[6.07908,-2.206645,6.544043],[8.639713,7.967947,11.79546],[-0.7113718,2.960058,3.204371],[4.988286,8.573673,9.969497],[-4.301503,-7.557499,8.753212],[7.970456,-6.633963,10.41814],[3.144678,8.330637,8.960385],[-4.955563,7.297562,8.877614],[6.771559,0.4401241,6.859134],[0.2522663,7.438441,7.509597],[-5.741791,9.877705,11.46897],[-0.002834536,7.033933,7.104662],[5.292902,3.186608,6.258537],[-7.451675,7.161532,10.3834],[-6.796552,8.698745,11.08428],[-3.521904,-4.797263,6.034694],[-3.546559,-6.791486,7.72673],[-9.960214,-9.814773,14.01912],[-8.79664,-7.164227,11.3889],[-1.187238,3.865725,4.165737],[8.571535,8.423882,12.05956],[8.433756,-3.330383,9.122483],[1.280959,4.176704,4.481709],[4.859322,5.570313,7.459316],[-3.527965,9.010957,9.728508],[-3.513196,-8.335714,9.100916],[1.323209,-4.507433,4.802898],[6.291962,-6.677093,9.228887],[4.85391,-4.003003,6.370595],[8.097837,-1.426059,8.283031],[1.341087,5.405328,5.658276],[2.768852,0.3351713,2.962917],[-5.96591,-1.048702,6.139369],[9.898243,9.204383,13.55345],[9.064104,6.470388,11.18141],[-6.962765,6.457912,9.54907],[5.385591,-9.138425,10.65436],[-2.879481,-0.9181861,3.18347],[-6.333363,-1.869528,6.678819],[-3.519088,-8.846932,9.573515],[-6.49287,-6.636992,9.33847],[6.871681,3.178057,7.636756],[8.953896,1.719174,9.172122],[8.583587,7.548365,11.47414],[-1.834816,-1.007321,2.319751],[-3.467229,-7.182373,8.03792],[0.937993,-2.727779,3.052967],[-7.853124,5.494349,9.636359],[7.529009,-6.300434,9.868204],[0.4841012,-8.753901,8.824123],[-2.520921,0.7145643,2.804575],[-3.904924,7.536099,8.546416],[2.54444,-9.460559,9.847657],[-8.409608,-2.035116,8.709949],[6.222199,-1.144862,6.405191],[4.041777,1.261301,4.350499],[-8.574822,1.523805,8.766389],[2.530471,8.909701,9.315904],[-6.876236,-5.032562,8.579586],[2.266021,3.242768,4.08049],[-2.8887,-0.2234807,3.06505],[-8.580026,6.004681,10.52013],[-6.366804,8.897782,10.98666],[-6.307605,-7.48838,9.841835],[-2.533159,-7.415536,7.899815],[-4.17792,-2.795904,5.125631],[7.470123,-6.048381,9.663625],[8.491374,-2.163338,8.819493],[4.554606,-5.359749,7.10432],[-2.017992,5.530516,5.971508],[8.559461,-1.929476,8.831038],[5.414207,7.047172,8.942945],[8.284659,3.08646,8.897293],[2.302403,-0.822262,2.641434],[7.864484,2.507157,8.314803],[-0.8254211,-0.8056305,1.526552],[-9.931387,0.4388292,9.991247],[1.151211,-9.726581,9.845387],[-6.518139,-9.838907,11.84442],[-7.782867,3.737033,8.691285],[9.007822,0.5141098,9.077729],[5.715435,-6.499272,8.712447],[-9.568385,-4.646553,10.68384],[4.772173,-5.041345,7.013472],[-2.406569,0.2801538,2.62108],[-3.004109,-4.448497,5.460202],[-5.228158,2.356396,5.821189],[-0.8205496,-0.4607213,1.373159],[-5.745471,-9.670102,11.29253],[-0.5490607,-0.9196024,1.465311],[5.345141,-0.8636051,5.506028],[2.396254,4.997114,5.631446],[8.421698,8.006932,11.66344],[5.967124,-0.9656389,6.12691],[-5.409371,-3.648533,6.600992],[-8.68363,-4.08937,9.650304],[7.137286,8.129291,10.86399],[1.9712,1.366222,2.598498],[-2.578277,8.371614,8.816544],[-4.72879,-7.233075,8.699358],[-2.239137,0.01216827,2.452322],[8.055722,-6.173094,10.19812],[4.124795,9.563983,10.46345],[-1.23172,-8.347272,8.49671],[-0.2258907,-6.558723,6.638364],[-1.274325,9.18642,9.328142],[-0.6958172,2.813554,3.065982],[6.11798,6.792273,9.195904],[0.4199962,-4.095762,4.23694],[0.1417869,0.1486078,1.020876],[-0.7565235,-6.574657,6.693164],[9.174877,9.56471,13.29143],[-3.266567,-8.444529,9.109365],[8.591202,7.280721,11.30565],[-7.766782,9.029067,11.95186],[-5.801124,-7.201901,9.301635],[4.569283,-5.351764,7.107723],[3.064027,8.408289,9.004864],[-3.694114,9.42176,10.16937],[-2.080818,0.5668491,2.377209],[-5.16078,2.177069,5.689753],[-6.589136,0.4477324,6.679608],[-3.127911,-8.528039,9.13845],[6.72633,0.2068954,6.803405],[-9.025735,-7.433285,11.73532],[-5.992648,-5.940362,8.497043],[3.900394,2.075798,4.530123],[-4.484804,3.791527,5.957277],[6.87816,-7.937977,10.55086],[-8.687666,-4.535642,9.851274],[-4.642164,-1.216177,4.901915],[-0.3067089,-1.225388,1.611101],[8.144987,8.89753,12.104],[-0.4832915,2.515373,2.749668],[-2.640578,-9.576648,9.98423],[-1.668362,0.4760434,2.002511],[7.676082,-1.016946,7.807458],[9.945144,3.958249,10.75052],[-9.405484,4.90723,10.6557],[3.290848,4.205225,5.432642],[-6.399308,3.395197,7.3129],[3.164304,-5.593063,6.503475],[8.867822,-5.596076,10.53349],[1.406776,9.221146,9.381287],[7.260642,5.892334,9.404069],[7.271656,-2.523253,7.761687],[-6.643249,-2.441034,7.147825],[-6.578782,-6.067419,9.005218],[-2.061043,0.4447858,2.33361],[5.685795,4.593679,7.377679],[-1.990419,7.370585,7.699824],[3.070413,9.440065,9.977088],[6.792292,2.566762,7.329631],[8.354246,0.9878359,8.471673],[-5.455597,2.748026,6.189927],[-8.793015,-9.925571,13.2979],[-7.197706,-7.34424,10.33174],[3.53539,-5.099667,6.285347],[6.122721,-9.816615,11.61265],[-0.9693071,-8.163097,8.281045],[-3.307282,-3.000848,4.576374],[-8.529248,-4.069667,9.503172],[-9.003881,0.7178803,9.087642],[-8.770717,-8.065888,11.95759],[0.6704646,8.58518,8.669189],[0.05736151,6.087897,6.169747],[6.258059,-3.320621,7.154706],[1.693938,4.917612,5.296446],[-5.811575,0.1034274,5.89789],[-6.842583,-5.485727,8.8269],[-5.27635,-0.2817936,5.377665],[-7.793952,-9.187656,12.08961],[-1.203214,-1.553657,2.204898],[2.958686,-9.569278,10.06603],[-2.59116,-1.496461,3.154918],[-7.520459,-0.1013831,7.587331],[-1.677687,-3.469742,3.981676],[7.318183,6.005535,9.519572],[-2.09769,-1.101953,2.571887],[-9.4067,-2.138538,9.69842],[1.802655,-7.127885,7.419994],[9.49691,-7.296138,12.01769],[-9.876997,3.608595,10.563],[0.7077064,0.114424,1.230423],[-5.410233,7.970649,9.685137],[-6.787271,9.734004,11.90873],[6.980133,-2.177514,7.379961],[4.781888,8.557879,9.854124],[-6.413872,8.138393,10.41015],[-6.456966,-5.041533,8.252846],[-9.498833,-6.782908,11.71476],[6.444081,4.73064,8.056373],[-7.915591,-6.771502,10.46469],[5.640233,-1.542073,5.932134],[-5.978745,6.308198,8.748644],[-4.64818,8.366838,9.623385],[2.872935,5.838068,6.583069],[2.443709,8.345262,8.753005],[-4.299799,-9.224187,10.22614],[-4.526229,-4.795987,6.66995],[-6.525562,-9.486687,11.55769],[8.119173,4.447168,9.311191],[9.963889,-7.431947,12.47048],[6.885765,3.729117,7.894306],[-7.745547,2.848518,8.313095],[-7.346038,7.798335,10.76003],[-9.565331,-4.590172,10.6567],[-2.080207,7.853196,8.18535],[0.2762056,5.840543,5.931967],[3.553886,-6.193332,7.210234],[-4.710434,1.213523,4.965967],[9.096131,-2.248778,9.423195],[5.603229,8.463868,10.19967],[8.752498,-9.342013,12.84054],[6.805602,-8.234557,10.72959],[-2.354595,-9.895359,10.22068],[-6.195251,-9.395123,11.29821],[7.051396,-0.05801253,7.122187],[-1.51895,-7.029256,7.260693],[0.5710357,-1.165546,1.638469],[8.38316,-1.309111,8.543486],[-7.584769,0.003224481,7.650407],[9.627203,2.901395,10.10451],[8.701589,-1.372108,8.865683],[2.349599,1.105853,2.782719],[-1.844727,7.30598,7.601339],[4.51886,-4.733383,6.620047],[-6.152534,8.048641,10.18009],[9.609738,-1.191238,9.73479],[-2.510166,-9.905671,10.26758],[-6.032952,4.628039,7.669111],[-9.032045,0.6683925,9.111783],[-3.314516,2.49281,4.26616],[9.188288,6.56379,11.33613],[-8.598759,3.626505,9.385637],[-6.74834,-6.286643,9.276958],[8.863714,7.005636,11.34215],[7.068816,-7.068326,10.04636],[4.458429,6.11583,7.634197],[-9.900784,5.023381,11.14719],[2.112659,-9.1468,9.440724],[2.36831,-3.426184,4.283413],[-3.644206,6.476838,7.498644],[0.02813184,7.893926,7.957063],[2.212456,-2.410221,3.421129],[9.157416,8.498555,12.5333],[-8.138524,6.500087,10.46359],[0.4170667,-3.210383,3.388289],[5.338561,5.993906,8.088704],[-7.296755,-5.61996,9.264264],[-0.03984099,3.084425,3.242725],[7.906857,5.409262,9.632159],[5.696153,1.465044,5.965946],[-5.54676,-4.15445,7.001857],[-5.793713,2.831575,6.525712],[3.170079,4.859415,5.887556],[5.224588,8.51926,10.04361],[-2.23266,6.549489,6.991465],[-6.722463,4.312824,8.049345],[7.85836,-2.704596,8.370702],[4.342119,2.407593,5.064632],[9.256673,7.891332,12.20488],[-1.041666,8.17993,8.306402],[4.711588,0.3116349,4.826612],[6.007087,7.804272,9.899078],[-1.523333,-6.341169,6.5978],[-2.321968,1.540341,2.960437],[3.005346,9.485739,10.00057],[-7.468617,-7.321367,10.50632],[4.700648,0.27274,4.813572],[8.644441,5.807283,10.46188],[-5.224864,9.661751,11.02944],[1.240047,9.206403,9.343209],[2.696657,-5.490264,6.19798],[7.721027,-6.455823,10.11395],[9.555413,8.194212,12.62739],[-1.966033,4.338076,4.86664],[-6.780257,-5.631152,8.870274],[1.432619,7.431477,7.634085],[9.205764,-6.19849,11.14304],[2.320294,9.456417,9.788135],[1.695852,-4.86161,5.245109],[-7.060384,-3.076795,7.766318],[5.654979,1.144943,5.85574],[-5.860897,7.481451,9.556266],[-9.090671,-6.211082,11.05522],[0.9908937,7.421219,7.553566],[-3.273881,3.159265,4.658246],[-5.444359,-0.2433726,5.540783],[8.641319,-7.963031,11.79331],[6.721455,-4.73029,8.27971],[-9.714141,-6.128237,11.52909],[9.137045,9.32002,13.09001],[0.3965994,7.540221,7.616575],[-3.225683,-0.05041387,3.37751],[-3.879953,4.478992,6.009609],[-0.9086257,1.449317,1.981444],[-4.28878,0.008847909,4.403829],[3.574009,8.891496,9.634949],[6.323309,0.8350607,6.456126],[-4.902241,8.94894,10.25258],[7.557505,5.451907,9.372255],[-8.173478,5.778943,10.05992],[6.604413,7.981186,10.40757],[9.706312,-6.623049,11.7931],[7.680892,5.182306,9.319463],[-9.95693,-3.108499,10.4787],[5.753266,-2.907979,6.523528],[-3.399902,5.199487,6.292376],[-8.099285,-6.52687,10.44981],[0.6323549,-7.040657,7.139379],[6.007474,7.809853,9.903714],[-6.27541,-6.315988,8.959491],[-9.576747,-5.913208,11.29956],[-1.350085,-5.739144,5.980009],[4.715127,-6.869254,8.391607],[-4.700253,-6.650467,8.204944],[-1.052815,-0.05567932,1.453107],[-4.202512,9.994757,10.88835],[-2.460445,-9.521001,9.884495],[3.277577,3.238438,4.714869],[6.502559,-3.107726,7.276072],[-9.319098,7.485343,11.99483],[5.054916,-9.039697,10.4052],[8.651809,4.774879,9.932435],[7.919752,7.032094,10.63827],[-9.342142,-3.313599,9.962708],[-6.119645,2.802037,6.804518],[-1.285865,5.749808,5.976097],[7.072503,4.836334,8.626147],[4.502183,5.781779,7.395852],[5.135805,-7.734806,9.338293],[8.010149,-7.023561,10.70014],[1.985417,8.246566,8.540944],[5.811714,-3.973017,7.110617],[0.2955076,-6.817666,6.896948],[-6.923053,-1.116154,7.083393],[1.011084,-6.1286,6.291425],[-4.780539,2.177586,5.347469],[-3.246174,-9.668044,10.24738],[0.8294585,-1.384034,1.898302],[-1.287022,-4.059087,4.374084],[7.659674,1.756609,7.921886],[-3.774156,8.564518,9.412503],[-7.166386,-5.144437,8.878193],[2.89283,-0.4807405,3.098318],[-0.5939379,7.674165,7.761802],[2.779091,7.315125,7.888879],[-2.934112,7.74526,8.342546],[-7.141005,-0.1173069,7.211637],[5.711208,-7.698069,9.637332],[7.562585,-1.776493,7.832536],[3.515161,-0.08226947,3.65556],[0.04887685,-4.180752,4.298962],[-9.418813,1.211273,9.548885],[-5.144225,7.987236,9.552958],[-8.493486,8.921155,12.35825],[-9.794226,9.685065,13.81041],[-6.622498,-2.170884,7.040612],[1.84421,-2.980332,3.644652],[5.551663,1.629288,5.871588],[4.866793,-6.271788,8.001312],[-3.883588,-1.720004,4.363561],[1.633496,-2.305029,2.99691],[-4.791341,-5.337177,7.241713],[-0.153621,-5.659955,5.749668],[4.91619,3.88261,6.343783],[-8.752515,-6.508084,10.9527],[6.670581,3.683131,7.685188],[6.959377,6.153641,9.343459],[7.711674,-6.639959,10.22541],[8.664337,-8.999319,12.5323],[7.190095,6.67357,9.86073],[-7.714577,-3.813399,8.663527],[-8.843811,7.645998,11.73347],[-1.236038,2.33171,2.822174],[6.86543,-9.782062,11.99262],[-4.81969,-4.016256,6.352931],[-2.716782,0.9658375,3.051843],[6.279785,3.226591,7.13068],[5.988545,7.062302,9.313366],[0.3497157,3.259281,3.427129],[-0.400224,-2.957847,3.147862],[-2.771131,5.727732,6.440969],[-9.060739,-9.464213,13.14033],[-2.686362,4.666284,5.476381],[-5.426566,-4.923005,7.394836],[7.001914,2.398277,7.468502],[-1.27111,-9.327424,9.466602],[-6.433523,7.228132,9.728109],[8.146307,4.436872,9.32996],[-9.416797,0.4648163,9.481146],[2.115339,-8.598804,8.911459],[0.1275717,-5.450495,5.542939],[-9.390732,6.101947,11.24365],[9.173678,1.481514,9.346189],[-3.514993,5.934194,6.969206],[8.872467,-3.105265,9.453218],[4.635201,0.3042365,4.751594],[0.8521704,-2.50765,2.83099],[-1.264532,-0.5772447,1.712382],[9.594765,0.3381583,9.65266],[-8.372036,-8.990547,12.32562],[6.153417,-4.174575,7.502774],[0.05528133,6.724392,6.798567],[-7.295497,-9.276687,11.84404],[7.213897,-7.255341,10.28009],[-9.690878,-7.606582,12.36015],[-6.816602,-5.490444,8.809712],[9.828436,-1.108944,9.941222],[2.461898,8.214951,8.634023],[-8.031164,-4.560596,9.289705],[-2.749908,4.296051,5.197888],[-9.660964,3.596449,10.35706],[1.745205,2.79693,3.445077],[5.656064,-0.5055121,5.765987],[-5.517062,-7.139977,9.078394],[-2.987266,-7.487069,8.122806],[-3.665091,-7.776266,8.654664],[3.908836,-3.311644,5.219768],[-8.000709,-7.430827,10.96488],[-5.536983,4.382638,7.13202],[5.073999,-9.928991,11.1951],[5.979218,-7.368114,9.541497],[-9.741603,-4.551137,10.79869],[-5.692732,8.883421,10.59822],[-9.520691,4.359598,10.51901],[-4.118342,6.827909,8.036236],[2.721877,-5.980291,6.646239],[-3.570181,3.029938,4.788185],[6.733442,-7.02594,9.782795],[9.138887,6.516537,11.26874],[-2.102552,9.265126,9.553183],[-9.850201,2.459501,10.20175],[-1.372385,-1.455762,2.236668],[5.112638,-2.416369,5.74264],[-0.9729397,8.068336,8.188081],[-8.845244,4.332002,9.899726],[-4.36812,-3.785616,5.86612],[-8.291666,0.3297997,8.358259],[5.203055,-4.293825,6.81973],[7.585103,8.379109,11.34651],[9.685826,6.426669,11.66693],[5.68926,-4.481298,7.310931],[-8.260763,-0.7904697,8.358531],[5.930641,-5.287305,8.008002],[2.115996,8.640783,8.952127],[5.371409,8.283354,9.923003],[2.569029,9.425808,9.820681],[-1.364214,0.643064,1.809589],[8.484279,-2.065276,8.789104],[1.919727,-8.023396,8.310248],[-4.85901,1.329636,5.135943],[5.782844,3.782283,6.981902],[6.844118,-4.294352,8.141462],[9.725384,7.933613,12.59068],[6.427474,0.776082,6.550933],[4.73378,-2.455393,5.425645],[-1.155095,0.9355804,1.791523],[-4.330986,6.633783,7.985268],[0.3120895,4.209913,4.338291],[-7.870994,3.749007,8.775397],[-0.5163726,2.609716,2.842052],[-1.974608,9.29828,9.558089],[-7.22308,-3.424016,8.055853],[8.58983,-3.005651,9.155278],[-6.137481,-8.10546,10.21602],[0.1509345,4.450158,4.563627],[-6.348761,-6.90885,9.436047],[-9.152072,5.700014,10.82823],[-0.7172685,3.248672,3.473952],[2.923665,4.441242,5.410402],[-2.979861,7.959851,8.557966],[-3.394727,-8.39138,9.107109],[2.347251,-1.456962,2.938082],[5.05723,6.286762,8.130126],[8.105049,9.156861,12.26947],[-6.956933,8.350556,10.9147],[4.714888,-5.170071,7.068224],[-9.30477,1.944593,9.558252],[-4.266386,-9.944015,10.86671],[-3.496035,-6.184555,7.174328],[-5.242507,5.198597,7.450456],[4.374572,0.4114265,4.506235],[-0.5373961,3.77066,3.937852],[5.224508,1.003154,5.413114],[-0.2638269,7.755358,7.824013],[6.042891,2.299481,6.542488],[6.961604,-2.098565,7.339475],[-5.976563,9.71031,11.44594],[1.531032,4.352366,4.720927],[-0.4361767,3.635112,3.795298],[6.651188,-8.053465,10.49269],[-2.094604,6.966311,7.342809],[8.144788,-7.383139,11.03849],[-0.4258833,4.337403,4.471515],[-7.691981,-6.565997,10.16262],[5.349708,2.027448,5.807746],[-9.986711,6.215974,11.80562],[1.390699,-3.992712,4.344628],[-9.976563,7.714045,12.65062],[6.307154,9.084847,11.10471],[-3.561745,-7.074172,7.983103],[3.837214,-4.057848,5.673654],[5.75317,2.948879,6.541778],[5.303426,9.384875,10.826],[-8.555077,2.325507,8.921733],[-5.507264,-3.052977,6.375784],[0.6302136,-0.628988,1.338953],[6.883108,-0.2492657,6.959835],[-7.661137,-0.3497174,7.734037],[6.809366,2.122405,7.202226],[-8.270583,-8.308651,11.76589],[-1.14907,0.266005,1.546325],[-1.331891,-2.732104,3.199738],[2.674235,5.148112,5.886815],[-5.435168,6.1555,8.272318],[4.759968,6.239854,7.911578],[9.939418,-7.696709,12.61076],[-5.412753,-3.040012,6.28805],[-8.404799,-0.4041634,8.473723],[-7.941461,9.65003,12.53754],[2.161283,3.524569,4.253672],[-2.584022,2.257523,3.574014],[3.719075,0.4108611,3.873025],[9.397729,-5.395422,10.88246],[5.493707,-3.001473,6.339531],[7.834795,-5.247403,9.482576],[-4.657112,-1.908442,5.131359],[-8.869919,3.889807,9.73684],[-2.898595,4.3419,5.315445],[8.609422,0.5180905,8.682774],[-2.119409,4.600142,5.162673],[1.427867,-9.888586,10.04106],[2.845344,-0.3770888,3.039437],[5.380958,-2.146908,5.879109],[0.2785438,5.524589,5.621269],[8.36851,8.07515,11.67219],[-5.47829,-3.0914,6.369334],[7.561787,-3.172675,8.261144],[9.905645,-2.393294,10.23961],[-8.689283,-1.221714,8.831547],[-2.009772,5.321614,5.775704],[3.303096,-2.229831,4.108842],[-3.00502,7.089302,7.764557],[3.832423,1.785506,4.344594],[-1.22665,-0.7005389,1.73073],[3.061147,-4.515763,5.546416],[0.2950814,9.956347,10.01079],[2.342495,-5.496477,6.057932],[9.704756,3.116942,10.24195],[5.292198,1.735932,5.658695],[2.588292,8.449986,8.893904],[-1.937122,4.032635,4.584167],[-1.629061,1.446221,2.396955],[4.388939,1.226719,4.665579],[-5.912254,-0.5252813,6.019192],[9.38886,8.074286,12.42356],[-3.721497,-3.53158,5.227007],[6.687713,4.050403,7.882339],[9.886271,2.925684,10.35847],[-7.33911,7.526112,10.55959],[-8.687601,-1.454512,8.865101],[5.457984,-0.2595952,5.554906],[-3.38138,7.87986,8.63284],[9.314254,-7.614414,12.07206],[9.717654,5.841908,11.38247],[-0.8098171,4.376632,4.561876],[7.515048,-1.390892,7.707822],[7.611284,-7.851686,10.98092],[3.060411,1.466578,3.537933],[-3.747635,-1.522596,4.166901],[-9.771288,9.061381,13.36363],[-9.713408,3.138533,10.25674],[-2.506113,-1.881171,3.289287],[-7.033851,-4.093378,8.199439],[8.790252,3.079035,9.367442],[9.582214,-1.590485,9.764654],[-4.51204,-0.6480803,4.666746],[2.913148,-8.22659,8.78426],[3.073523,2.457875,4.060504],[6.849484,3.938343,7.964043],[-2.189016,-8.700482,9.027191],[-2.726192,-8.054372,8.561835],[-7.967922,7.574862,11.03931],[-5.646032,8.481114,10.23753],[0.1006383,2.149233,2.372621],[-9.330975,8.562786,12.70387],[-3.647635,-7.956758,8.809951],[4.140971,6.679741,7.922536],[0.9725227,-1.204681,1.84311],[-2.875521,-7.934724,8.498734],[3.659951,6.148551,7.224951],[6.285611,0.6529648,6.398067],[4.303754,5.399015,6.976508],[-7.389679,-5.060134,9.011787],[-7.489662,-9.477468,12.12095],[8.218142,7.660279,11.27908],[3.308093,-5.689506,6.656873],[6.800603,2.850459,7.441325],[-8.042593,-7.117435,10.78616],[9.796212,-0.2718988,9.850873],[-9.643706,4.984173,10.90152],[-5.091618,4.836809,7.093609],[-6.061507,1.343542,6.288638],[8.900368,8.189474,12.13606],[6.420856,-5.281527,8.373882],[9.673691,-8.02802,12.61069],[-9.236411,8.962557,12.90886],[6.111218,-2.311999,6.610017],[4.618452,-1.548904,4.972846],[5.262911,8.928064,10.41194],[-0.712313,-2.685911,2.95322],[-7.49828,-6.241155,9.806948],[5.803868,2.012737,6.223825],[8.520808,6.568493,10.80506],[-3.79365,6.086175,7.241085],[-0.2100788,-0.906535,1.365994],[7.461927,3.229228,8.191963],[-4.014753,3.97233,5.735647],[-5.488386,-9.542683,11.05374],[6.355106,6.217136,8.946517],[-8.47031,7.459631,11.33103],[3.687391,-3.181161,4.971583],[1.00514,-3.314373,3.60491],[0.4590615,3.746172,3.904426],[6.354249,-5.435336,8.421363],[-2.72963,-8.211126,8.710539],[3.145156,-0.9194884,3.425999],[-3.422776,-6.361998,7.293175],[6.83184,9.609143,11.83257],[-5.112109,-1.922883,5.552579],[-5.737499,2.370407,6.287903],[-8.785292,-4.460123,9.903234],[1.917969,-2.711665,3.468679],[-4.034719,-9.761673,10.60986],[-0.704213,7.331947,7.43326],[-3.407119,-0.3346226,3.566571],[-5.843982,7.806048,9.802373],[5.677042,4.569586,7.355944],[9.42035,4.769665,10.60626],[6.611796,3.854523,7.718367],[4.871259,-1.757893,5.274405],[-8.599042,4.399346,9.710704],[8.596562,4.167639,9.605733],[0.7838565,-2.431814,2.743748],[-7.997217,-8.364579,11.61558],[-1.952592,1.948595,2.934219],[-7.393856,-9.232882,11.87077],[-3.896001,-6.842731,7.937366],[-9.529223,5.078297,10.84413],[2.598686,-0.988889,2.954839],[-2.177667,-8.75069,9.072861],[2.466465,0.4048691,2.692093],[-2.003557,-4.448223,4.980053],[4.968081,-8.890176,10.23314],[7.65771,-4.287101,8.832879],[-6.867488,8.425926,10.91598],[0.400919,-4.90643,5.023325],[-8.545948,2.604914,8.989928],[4.081016,7.779835,8.841975],[-0.8879083,-0.1065813,1.341544],[-9.426552,5.551732,10.98552],[-2.155964,2.346442,3.339756],[1.036439,2.721967,3.079498],[1.53487,-7.005282,7.240843],[-4.863965,-8.549818,9.887241],[2.899323,2.800843,4.153408],[4.84918,2.371097,5.489686],[-3.035087,-0.9401883,3.331022],[5.219366,2.129823,5.725201],[-3.275616,2.877387,4.473144],[1.346292,9.667265,9.811652],[9.499433,-8.969755,13.10327],[-1.215968,2.830204,3.238616],[-7.64644,3.633471,8.524679],[5.66618,0.6697079,5.79259],[5.546146,-9.375782,10.93915],[7.142503,-1.368971,7.340942],[8.766945,-9.774076,13.16783],[-9.867897,-7.895108,12.67707],[-2.833221,2.83088,4.128078],[-6.134891,9.111846,11.03008],[-0.06419513,9.07541,9.130564],[-9.295587,-2.840542,9.771214],[0.495623,2.029551,2.316186],[-6.160142,8.467386,10.51874],[4.948281,3.558643,6.176522],[2.179674,-5.345382,5.858676],[-5.873667,-2.838712,6.599867],[8.115447,-1.417761,8.298826],[-8.273927,-5.048388,9.743925],[2.562944,3.539474,4.482919],[-8.172512,1.985257,8.469427],[6.368373,0.436836,6.461192],[-0.8734065,-0.2534992,1.351703],[7.623536,-7.200106,10.53375],[-1.199864,6.344592,6.534028],[-2.314285,8.35978,8.731658],[-1.854185,-3.868324,4.404763],[-5.225898,1.703422,5.58674],[8.162394,7.273822,10.97876],[-4.53942,-4.512261,6.478181],[-6.27558,-6.751106,9.27148],[-6.430559,-1.377687,6.652075],[8.490353,6.328948,10.63681],[-9.783928,0.1828912,9.8366],[-0.4069764,-3.467886,3.63206],[0.9609663,6.996747,7.132876],[5.680158,4.494083,7.311702],[7.81978,7.879951,11.14642],[3.696755,0.2708317,3.839186],[6.753907,-6.788193,9.627815],[1.644007,-1.456293,2.413203],[-0.4227698,9.492949,9.554832],[1.54039,-2.564134,3.153979],[0.6216,-3.085399,3.302434],[0.4085452,-5.1043,5.217354],[-8.58848,3.80134,9.44522],[-3.201164,5.323492,6.291821],[3.609313,6.981125,7.922327],[-9.920141,-8.02449,12.7985],[-1.795236,5.646761,6.009058],[-7.231731,-2.202168,7.62545],[-7.576827,-7.799958,10.92006],[-4.150321,-1.925343,4.683173],[3.054161,1.029894,3.374697],[2.167949,4.965298,5.509463],[8.492682,-7.338418,11.26845],[4.036742,8.873575,9.799776],[1.238819,-2.10774,2.641447],[5.53414,1.267627,5.764858],[5.209553,2.438741,5.838399],[-5.804586,-4.038043,7.141358],[3.302181,7.341842,8.112154],[0.4445978,-6.78814,6.875792],[-8.180478,-4.373916,9.330132],[9.765917,-2.499104,10.13009],[3.551634,-2.781557,4.620732],[-0.9520357,-5.446328,5.618617],[2.539263,0.6036174,2.795033],[-2.575618,-4.461639,5.24786],[3.921132,-3.903223,5.622315],[0.08766359,3.615979,3.752731],[6.45107,-0.2302451,6.532176],[-8.230626,-1.789837,8.482141],[1.486085,3.349394,3.798274],[1.458294,0.6839184,1.895881],[-8.351471,-5.954769,10.30564],[-1.242331,6.992116,7.171686],[-8.512174,2.379429,8.894874],[-5.539298,-2.208616,6.046636],[9.883962,1.407377,10.03361],[-2.948647,9.470695,9.969381],[-2.004489,8.043962,8.350048],[4.317366,4.297235,6.172996],[0.7603397,5.014472,5.169434],[-1.03705,2.728283,3.085288],[-9.839784,2.771929,10.27156],[-5.109217,5.539526,7.602003],[1.258254,2.065607,2.617238],[-7.590696,-7.505877,10.72179],[-4.392595,0.8461297,4.583757],[4.731454,-4.357172,6.509348],[4.531736,-3.610459,5.879799],[9.422297,-2.922859,9.915784],[6.549312,-6.505973,9.285536],[-8.341627,2.976743,8.913122],[3.28631,-6.167225,7.059355],[1.212044,7.902057,8.05677],[4.718584,-3.008519,5.684736],[2.576645,5.792448,6.418065],[3.622828,-2.033818,4.273324],[9.028528,0.02504543,9.083775],[9.457362,-6.049624,11.27119],[-9.7353,-8.349508,12.86431],[4.581546,1.633982,4.96593],[-3.218571,-4.515484,5.634607],[-8.826765,-5.317999,10.3534],[7.906235,7.640223,11.04],[-8.593467,-4.036452,9.54676],[-6.583915,-8.07842,10.46942],[0.01355066,8.782938,8.839693],[7.725895,-0.0001989026,7.790344],[2.526086,5.634169,6.254995],[-9.814339,-5.97842,11.53528],[-3.017777,-7.206866,7.876922],[-4.998536,-7.955591,9.44864],[-4.89253,1.426206,5.193353],[7.91209,-1.596901,8.133343],[-7.134211,5.966797,9.354124],[-7.623719,-9.315269,12.07871],[7.411042,3.109991,8.099111],[-9.956186,0.9988199,10.05601],[3.685178,-6.114183,7.20859],[-7.470538,2.709942,8.00954],[-8.710026,1.586316,8.909598],[-3.242844,2.730262,4.355498],[8.37985,-9.908996,13.01576],[-7.993836,-0.0338022,8.056212],[2.334594,8.231142,8.614059],[-1.287272,-3.404085,3.774237],[-6.653809,-0.7597569,6.771292],[-2.318077,0.1299019,2.527915],[-3.566018,-3.487973,5.087479],[5.030667,-8.7639,10.15448],[-0.4099193,2.558668,2.777556],[-1.607957,0.4748652,1.952184],[1.13855,8.227063,8.365457],[2.278853,8.162086,8.533043],[-8.028608,-3.276806,8.729033],[-4.942138,-8.991689,10.30899],[4.583262,0.4444767,4.712097],[-6.995543,-2.365533,7.452072],[-8.738464,0.4322799,8.806113],[3.36752,-0.4778393,3.545211],[-4.933935,-5.145176,7.198371],[1.661434,6.284175,6.576566],[2.247256,2.871178,3.780717],[8.341209,-1.446343,8.524534],[8.640924,-8.922264,12.46083],[9.148058,9.873708,13.4973],[7.048779,-1.982616,7.390268],[-0.4077353,-0.6749324,1.273492],[-3.441042,7.156018,8.003085],[7.109136,-5.202234,8.865836],[-8.51216,5.200642,10.02515],[0.110187,-8.898687,8.955378],[8.19046,-2.312141,8.569109],[-0.7714875,8.794306,8.884538],[5.291195,2.369772,5.883244],[2.290035,-3.823735,4.567846],[8.365217,-3.792887,9.239202],[3.914841,-0.9055041,4.140763],[0.7836291,-4.072394,4.265967],[-5.321082,5.440382,7.675393],[7.224546,9.662083,12.10578],[-6.269326,2.559697,6.84518],[5.010461,-3.619702,6.261546],[-8.167981,-3.235171,8.842072],[-6.41914,-6.093175,8.906859],[-1.836858,-3.862442,4.392323],[-6.693031,-6.203371,9.180331],[6.487657,8.567205,10.7929],[1.227096,1.886894,2.462952],[-5.739096,2.43817,6.315211],[0.05224429,-9.298409,9.352172],[6.49795,3.007214,7.22957],[-2.603885,6.166858,6.768334],[-6.383596,-3.87799,7.535855],[-8.229796,-1.752917,8.473621],[-1.121638,6.225115,6.403915],[-4.676883,-5.206608,7.069795],[8.607748,4.204974,9.631986],[1.449065,-3.403363,3.831798],[5.60028,2.50017,6.214015],[9.470739,-7.789732,12.30345],[6.275753,-1.104111,6.450127],[-3.176027,-3.635947,4.93024],[7.188738,9.777526,12.17694],[-8.702328,-3.398479,9.395752],[-0.4473366,-0.2466115,1.12291],[9.868558,5.265122,11.22987],[-7.603662,5.264144,9.301983],[1.539964,7.210263,7.440388],[4.820233,9.385746,10.59844],[-1.083886,8.240823,8.371737],[3.130792,8.966063,9.549459],[7.783566,3.700001,8.676053],[-3.753003,-9.901081,10.63562],[5.966269,1.42789,6.215724],[4.992762,-5.363596,7.395663],[-8.597426,-9.813262,13.08495],[-8.804231,-5.534162,10.44708],[-9.456139,-8.532581,12.7759],[-4.848384,-9.449927,10.66808],[3.437904,-1.69978,3.963387],[-7.806082,3.604158,8.655915],[-7.816789,1.291956,7.985696],[2.09898,5.968791,6.405636],[9.902542,5.849402,11.54452],[-0.7812575,3.470397,3.695135],[6.131123,-1.289393,6.344541],[-2.280894,0.4838836,2.537049],[1.392283,0.3589211,1.751364],[5.262171,7.133982,8.920995],[-4.395662,-7.302435,8.581806],[3.679153,5.143423,6.402419],[-6.922775,-3.049582,7.630516],[-7.011956,6.381824,9.533897],[2.749953,-1.380092,3.235258],[-3.169292,-2.804842,4.348741],[9.347296,-1.025422,9.456396],[5.317281,-8.352743,9.951974],[7.877711,-5.106984,9.441378],[-6.817086,1.664488,7.088242],[1.221074,-8.525931,8.670786],[5.847062,-1.533443,6.126956],[-7.332918,-1.228738,7.502099],[-9.506182,-7.040963,11.87193],[0.8802206,-6.052578,6.197458],[-7.718095,5.147407,9.330852],[6.239284,2.275912,6.716282],[2.088439,-3.296391,4.028371],[-4.885619,1.261187,5.143915],[2.29446,-5.161702,5.736525],[-4.315888,-0.427763,4.450828],[5.178983,7.132309,8.870833],[6.416711,-8.093504,10.37685],[-4.668981,7.766009,9.116485],[-6.567901,-6.958319,9.620578],[-5.68936,6.735595,8.873391],[-7.544181,-5.604098,9.450956],[2.374298,-9.859146,10.19019],[3.701605,7.812054,8.702302],[0.4328016,8.734986,8.802688],[-2.707819,9.823438,10.23876],[-7.015729,1.687685,7.28483],[7.128049,-4.114085,8.290645],[-4.506529,2.17575,5.103204],[-2.176206,6.114533,6.56684],[-8.853323,0.4589798,8.921434],[8.761862,-5.620487,10.45754],[-7.432215,5.14282,9.093207],[-7.975033,6.203743,10.1532],[-2.271354,7.525535,7.924186],[-7.913134,-6.978258,10.59782],[8.853261,-2.216847,9.181211],[-4.262769,1.346763,4.580936],[7.881732,3.646243,8.741669],[3.767196,-5.796587,6.985141],[9.663768,-3.231267,10.23863],[-0.838113,5.440815,5.595078],[0.01758915,9.283054,9.336777],[7.18789,4.472132,8.52442],[-4.496712,-5.753802,7.370662],[0.4737675,3.855541,4.011191],[-6.395357,-3.903155,7.558784],[-2.071427,0.2525704,2.314002],[-0.583937,-1.922532,2.244351],[6.40242,-7.378621,9.820134],[-9.901577,4.32034,10.84927],[-9.15307,8.737449,12.69337],[-9.562348,7.469228,12.17489],[-8.502807,5.833885,10.36011],[-2.103853,1.51075,2.77643],[0.6278409,6.974176,7.073423],[-3.76777,8.062303,8.955267],[1.628041,-0.4235212,1.957009],[1.834822,-1.520167,2.584082],[0.8803871,8.174051,8.28192],[-2.864612,3.478539,4.615868],[8.502811,-2.216656,8.843719],[-4.292207,5.367061,6.944665],[3.665708,2.218157,4.399731],[-3.282311,-3.815122,5.131152],[-7.334166,8.857468,11.54317],[-2.108232,7.612644,7.962223],[-7.679197,-2.714946,8.206156],[6.280956,2.978266,7.022854],[8.575798,4.261388,9.628278],[-1.176049,7.865938,8.015988],[-6.910703,9.303584,11.63248],[6.801386,5.664369,8.907521],[-0.1826481,-4.544296,4.656607],[-0.6958897,-8.236432,8.326048],[-6.372269,5.869689,8.721184],[1.164694,3.561532,3.878275],[-6.570951,-2.050509,6.955717],[3.26231,5.26378,6.272961],[9.028268,-1.615678,9.226052],[-4.494044,-4.23007,6.252193],[8.239147,-9.330477,12.48765],[-7.753506,-0.6346967,7.84345],[9.594437,-1.063586,9.704866],[8.27592,-2.296895,8.646767],[-1.527382,-7.628835,7.844235],[3.135595,6.379183,7.178156],[1.281318,1.168977,2.00207],[3.354242,5.211828,6.278064],[6.797619,1.093817,6.957303],[-3.8277,-6.451684,7.568058],[-4.386517,9.210752,10.25083],[-5.361558,-2.414589,5.964608],[4.077808,1.411443,4.429525],[-0.3566394,8.445326,8.511799],[0.4379539,-3.192174,3.373689],[2.614228,1.382634,3.121837],[-4.522262,2.579446,5.301357],[4.852912,9.44384,10.66475],[-8.778483,1.609621,8.980682],[0.7209939,2.097156,2.432672],[-1.352192,-9.621225,9.767108],[2.43314,0.1913206,2.63757],[5.420979,-7.458278,9.274316],[-4.339304,-0.6493874,4.500141],[-6.930082,4.164631,8.146791],[-7.389629,-6.509341,9.898391],[-6.570888,-5.396851,8.561691],[2.457463,6.661322,7.170239],[6.415688,0.9102389,6.556644],[6.140507,-4.608752,7.742508],[7.654753,7.662013,10.87666],[8.578066,-8.720964,12.27349],[-1.410535,9.227992,9.38858],[-6.033105,-5.634171,8.315181],[8.878057,1.694841,9.093535],[2.216245,-6.124629,6.5896],[-8.759223,2.822046,9.256778],[-8.151454,-1.765532,8.400197],[-4.614935,-4.416363,6.465438],[0.580506,-2.147643,2.439131],[-7.554883,9.255985,11.98956],[8.039597,6.923586,10.65698],[1.294822,6.451134,6.655351],[0.5737064,1.917551,2.237441],[-0.2189846,7.62392,7.692341],[7.815151,4.466852,9.057006],[-2.230784,0.767038,2.562176],[1.735363,-1.973717,2.811947],[-0.8918424,2.172843,2.552769],[6.036754,-8.545058,10.51001],[5.601614,2.373123,6.165208],[-7.324876,1.060416,7.468486],[-1.284686,7.063678,7.24886],[-6.125092,6.264773,8.818398],[-3.673465,5.753618,6.899164],[-1.094023,-9.221831,9.340185],[7.110207,3.60023,8.032228],[-4.654252,-7.49536,8.879329],[-5.158544,9.769765,11.09319],[-6.350033,-2.80598,7.014018],[-1.342552,-6.484852,6.697443],[4.64857,9.549127,10.66747],[2.224705,5.961815,6.441471],[-8.019834,0.1616334,8.083555],[-5.288442,-7.935313,9.588368],[0.9381549,4.47694,4.682214],[2.146832,6.086014,6.530579],[-6.820617,8.808073,11.18494],[1.073914,5.406758,5.60235],[7.458989,3.638862,8.359296],[-3.144139,-1.922775,3.818727],[5.676257,2.776208,6.397439],[-1.884646,-6.292674,6.644519],[2.678828,-8.324955,8.802329],[-6.694368,4.006064,7.865311],[-8.935085,0.57959,9.009532],[-2.037726,2.866172,3.656127],[-1.378841,-3.47156,3.866902],[7.195556,-1.368278,7.392443],[9.189773,-1.452697,9.35747],[-2.250911,3.380278,4.182449],[9.742982,0.1139229,9.794829],[-2.002393,-1.46029,2.672456],[5.452972,-3.386751,6.496537],[-7.601783,3.563102,8.45475],[0.8048021,8.852252,8.944836],[-6.482598,-3.680233,7.521183],[9.397771,-7.175802,11.86635],[-9.992327,9.915764,14.11272],[2.761565,2.253851,3.702173],[1.579102,-2.120915,2.826985],[-6.466414,7.706761,10.10983],[6.286707,-0.967658,6.43887],[-9.807929,3.922649,10.6105],[0.1588592,-9.758213,9.810605],[0.9871626,7.163414,7.299931],[-5.67415,8.944867,10.63986],[-8.878024,2.234915,9.20946],[1.425817,9.028227,9.194663],[5.05016,-6.314926,8.14754],[-0.01506646,9.522015,9.574392],[7.66149,-9.332746,12.11605],[9.517308,-2.339504,9.85152],[3.339201,-7.324954,8.112041],[0.3698737,-9.019782,9.082582],[8.68538,-6.043775,10.62841],[-7.102218,-1.712586,7.373903],[-3.639417,8.757776,9.536456],[-3.060217,3.024961,4.417614],[6.136308,-3.579839,7.174226],[-6.271366,-5.309749,8.27789],[-0.02308518,7.861613,7.924992],[3.804981,2.848308,4.857029],[-4.973302,4.514073,6.790478],[-6.193156,9.112843,11.06341],[-3.326665,7.388991,8.164796],[9.344018,-9.423883,13.30865],[-1.189527,4.666159,4.918131],[6.26242,-5.895755,8.658974],[-7.232396,5.33037,9.039933],[-9.910393,-5.096123,11.18867],[8.586522,-2.835258,9.09764],[7.568013,9.056513,11.84463],[9.585442,6.698571,11.73676],[4.386148,-7.281508,8.559127],[0.1045197,-0.3170305,1.054245],[6.643882,-1.7954,6.954468],[6.332513,-8.384063,10.5543],[-7.816001,-1.40168,8.00341],[-3.013778,-5.9563,6.749843],[4.089204,5.651335,7.046926],[4.808768,-5.773274,7.579904],[-0.1417471,-6.633094,6.709548],[9.291702,9.527193,13.34553],[5.056002,-4.849937,7.077079],[-4.577136,8.676988,9.861049],[-1.037334,3.789022,4.053733],[-1.030892,5.166046,5.361975],[-8.145175,-1.940997,8.432755],[6.990674,-8.781493,11.26872],[-6.316634,5.520186,8.448214],[8.731665,-3.539948,9.474873],[-6.39118,0.858331,6.525635],[0.9387794,-8.334926,8.447028],[9.591467,-6.420211,11.58513],[5.207167,6.152712,8.12222],[-9.841431,-0.8555321,9.929033],[1.319028,8.397737,8.559311],[-9.530882,8.586939,12.86753],[8.269636,-0.8072497,8.368902],[3.162071,-7.132278,7.865627],[9.175628,7.813777,12.09327],[4.769326,7.845448,9.235666],[7.822712,-3.488331,8.623415],[0.4605006,6.694979,6.784895],[8.513551,-1.144994,8.648211],[1.738642,-7.085684,7.364088],[-3.055666,5.637854,6.490184],[2.717718,-7.909828,8.423264],[5.122303,-9.203934,10.58066],[2.243475,-7.662331,8.046397],[-8.278835,-7.205778,11.021],[-4.389565,-9.346236,10.37403],[-5.575451,0.1534223,5.666497],[7.921429,-9.041698,12.06239],[2.471972,-3.363847,4.292564],[9.804796,8.126839,12.77417],[-5.766611,-2.140164,6.231702],[2.815241,4.238675,5.185744],[1.195815,-2.386564,2.850555],[2.328872,4.043753,4.772377],[3.174418,-3.627871,4.923248],[1.194175,-5.215586,5.443196],[7.99209,-9.440923,12.40986],[7.33506,-1.372923,7.529145],[-3.875699,9.138737,9.976851],[-1.093602,-0.02147418,1.482035],[0.5499411,-8.830357,8.903798],[7.357243,3.140731,8.061836],[-9.647991,5.660707,11.23064],[2.167501,-0.9435672,2.566784],[1.973622,-9.091444,9.356791],[-0.5766166,6.523121,6.62447],[-3.645449,3.816537,5.371709],[3.028321,-3.199838,4.517708],[-4.259498,-2.37993,4.980702],[1.573673,-3.574972,4.031981],[9.928385,1.432314,10.08089],[5.668653,-4.023506,7.022978],[-3.015524,-6.736376,7.447963],[8.741213,-0.08338013,8.798622],[-4.664971,-3.820917,6.112394],[-8.833728,8.225812,12.11193],[-9.139448,0.8278351,9.231187],[8.898272,-0.04355759,8.954392],[9.644762,9.09332,13.29323],[-3.540021,9.289449,9.991277],[5.201117,-7.782029,9.413373],[-8.826631,-8.85052,12.53958],[5.789567,-2.985605,6.590365],[2.089862,6.618273,7.012065],[3.846438,-1.821747,4.371939],[8.172541,6.143625,10.273],[8.927029,4.608509,10.09605],[-4.411568,-4.609733,6.458449],[-5.735736,-8.023814,9.91364],[8.031073,6.148355,10.16368],[5.306867,-1.68213,5.656183],[2.992514,1.44649,3.470947],[-9.735484,9.176992,13.41629],[-8.850456,3.223975,9.472307],[5.359424,-4.959246,7.370044],[-6.850361,-7.549759,10.24335],[1.390629,-6.773109,6.986334],[6.527538,-1.252648,6.721449],[7.724956,-0.6646653,7.817719],[-9.097095,5.674687,10.76844],[-7.881074,-0.4224228,7.955487],[-8.019501,1.667505,8.251846],[0.6658719,5.052963,5.193826],[-3.719379,7.9034,8.791901],[9.648835,-7.592917,12.31878],[-8.421021,3.082421,9.023022],[9.322601,-7.927401,12.27822],[2.746874,4.97211,5.767772],[-6.321476,-6.113978,8.85109],[-8.943552,-7.790077,11.90262],[3.735515,-7.85208,8.75267],[-2.271551,4.570876,5.201235],[-5.290196,-9.966021,11.3273],[-8.550871,3.290435,9.216527],[-4.704498,8.763965,9.996968],[4.555178,7.043856,8.447813],[-1.393412,-6.027716,6.266973],[6.197787,-2.924438,6.92567],[-7.470422,2.070308,7.816226],[-9.002561,-5.003296,10.3479],[4.021904,-8.849232,9.771623],[-4.371559,-0.5781825,4.521595],[-1.394048,7.005796,7.212805],[5.151854,5.657102,7.716502],[0.783473,-5.586363,5.728987],[-3.893329,6.327384,7.496252],[-4.106097,-0.03957271,4.226298],[-0.6583478,-4.874676,5.019551],[-4.296177,-3.784023,5.81171],[2.435812,1.556926,3.058954],[-2.280282,-3.362389,4.183939],[5.392343,-6.210691,8.285532],[3.902527,2.041545,4.516373],[-6.583261,4.775305,8.194076],[-6.536058,-6.967498,9.605524],[-8.436217,6.876896,10.92984],[-3.190295,6.07799,6.936854],[-1.69324,9.154835,9.363656],[-6.37684,-8.837935,10.94409],[3.026399,1.806178,3.663519],[-4.392001,-7.701091,8.921685],[9.019797,-6.558826,11.1971],[-7.64614,-7.250592,10.58464],[-8.776855,-9.227051,12.77387],[8.761072,2.539946,9.176476],[4.766675,9.24598,10.45033],[7.529822,-4.864341,9.019979],[3.063493,-8.498067,9.088572],[2.558691,-5.823748,6.439172],[-5.859529,-0.5744239,5.971939],[-0.108496,8.145699,8.207568],[4.801093,-2.700113,5.598312],[-0.3307176,-0.9465841,1.41612],[9.975103,-5.641694,11.50354],[-6.956714,-1.798866,7.254777],[3.258018,-4.504288,5.6483],[-8.936925,-8.976037,12.70582],[-8.590981,-6.434268,10.77983],[7.084657,-6.528081,9.685464],[-5.187905,2.153354,5.705374],[-0.977889,-7.181745,7.316674],[-6.881467,4.318819,8.185767],[8.238012,-7.040447,10.88268],[2.254118,9.661015,9.97077],[-3.607501,-0.7381164,3.81561],[8.605948,4.560483,9.79083],[-5.889764,4.122805,7.258571],[9.474532,-1.491558,9.64321],[-8.784048,-2.765752,9.263308],[-0.931716,-4.695389,4.890274],[-5.544023,5.933953,8.182175],[-7.793084,1.441166,7.988061],[-5.456663,-8.520458,10.16727],[-4.4685,7.284647,8.604278],[-0.4957681,0.5219213,1.232148],[8.088213,-6.67859,10.53673],[7.223837,-9.917446,12.31014],[-4.832296,-5.756984,7.582476],[-4.832038,9.176151,10.41875],[6.330236,-2.287717,6.804818],[-2.143137,-1.943817,3.061284],[7.444915,0.04331499,7.5119],[-4.550204,2.669474,5.3694],[6.337031,3.287456,7.208698],[-6.49394,1.94681,6.852834],[-3.683126,-4.98513,6.278291],[3.347065,0.6817288,3.559157],[8.116897,7.779674,11.28749],[-6.062768,9.108914,10.9877],[5.748271,6.935249,9.063128],[-5.349812,1.882278,5.758772],[-2.015656,-8.373559,8.670604],[2.711586,0.1076271,2.892106],[-6.189147,-3.397052,7.130603],[6.975004,-5.159878,8.733557],[-4.958733,8.991615,10.31689],[7.526658,6.973432,10.30919],[8.567953,-0.9580713,8.679154],[9.975401,1.598837,10.15209],[6.041345,-9.076909,10.94934],[0.720969,1.088291,1.644437],[-9.08722,-1.188121,9.218958],[-5.056805,-6.935944,8.641678],[-9.950916,1.643418,10.13517],[-1.629544,2.149483,2.87675],[-1.049801,-1.286897,1.938604],[-3.32688,2.700846,4.400307],[-5.421143,6.374725,8.427687],[4.811043,2.814184,5.662665],[-8.4691,-5.725873,10.27187],[-6.151237,1.877115,6.508554],[5.530796,-4.031139,6.916631],[7.170907,0.4671847,7.255354],[7.624505,-3.574853,8.480132],[6.535028,-6.724841,9.430275],[-8.015165,0.3090455,8.083217],[-9.763792,-8.852089,13.21708],[-4.966646,9.179585,10.48486],[8.264064,-9.493368,12.62612],[0.792177,4.757202,4.925293],[-7.544542,0.2121989,7.613484],[7.418346,9.358835,11.98414],[0.6238882,-3.55512,3.745412],[4.611703,-0.7680504,4.780974],[2.755252,8.299909,8.802267],[2.670033,4.995256,5.751667],[4.142551,-3.068704,5.251445],[-0.6547215,4.372291,4.532724],[5.292856,2.901232,6.118127],[-2.095939,3.986862,4.613895],[4.869792,7.657408,9.129664],[-6.415383,7.651375,10.03497],[-4.856083,6.020186,7.798986],[-9.141344,4.07512,10.05837],[0.01039044,-4.947556,5.047616],[-1.786279,8.403962,8.649704],[8.07268,8.399835,11.69296],[6.13905,3.425942,7.101058],[5.234118,-4.162678,6.761944],[4.071911,1.868322,4.590325],[3.754162,-2.850402,4.81856],[-1.249953,-4.171374,4.467969],[-9.665374,4.591088,10.74698],[-1.538786,-9.451509,9.628027],[3.459909,-6.121176,7.102096],[-8.718925,8.848378,12.46248],[-1.721624,-1.390285,2.428349],[-9.694456,-7.900251,12.54577],[-1.322114,-4.182253,4.498803],[1.961344,-7.368572,7.690431],[7.04318,-4.892656,8.633914],[-3.476261,1.772403,4.028127],[5.383033,4.757623,7.253415],[6.691832,3.184972,7.478279],[7.325391,-6.526647,9.861971],[1.287592,3.93273,4.25726],[-3.995226,-8.459115,9.408424],[-1.07775,-1.074826,1.821207],[-5.844141,2.550693,6.454457],[-8.540077,-2.51236,8.957951],[8.675397,-7.796515,11.70676],[7.784962,8.286751,11.41385],[-0.7455349,-0.6947383,1.427755],[-5.850648,6.104971,8.514738],[-3.595155,-2.165032,4.314221],[-3.448289,4.770738,5.970816],[-1.969774,0.1866018,2.216941],[5.405728,-2.552995,6.061327],[6.135682,2.741385,6.794246],[-5.114353,6.120713,8.03864],[-3.072687,-3.195878,4.544782],[0.9492246,0.9237468,1.659619],[-9.078024,1.475537,9.251364],[-3.824661,1.910584,4.390714],[-9.181222,5.48505,10.74154],[-7.568254,2.270387,7.964491],[1.718228,4.593277,5.005047],[-6.650422,7.433678,10.02435],[-5.867073,6.514933,8.824222],[6.994081,7.994789,10.66929],[-8.978614,-7.681035,11.85807],[4.02502,-5.457002,6.854171],[-1.260802,0.7790422,1.787884],[-7.328275,-2.503732,7.808475],[0.7039286,-6.548667,6.661873],[0.2769133,5.23121,5.333127],[8.422505,-8.56652,12.05503],[-7.559659,-9.154775,11.91463],[2.533987,2.124338,3.454547],[1.38788,7.290708,7.4887],[7.434774,-8.968016,11.69193],[5.004758,1.641106,5.361048],[-6.868747,1.726777,7.152723],[4.383638,-3.968881,5.997358],[7.916744,-9.508169,12.4129],[-7.566618,-5.755548,9.559291],[-4.09993,5.438019,6.88342],[9.56525,5.999496,11.33525],[-7.727228,9.550424,12.32561],[8.463212,-4.984123,9.87256],[-5.784758,-4.220578,7.230263],[-1.402976,5.879892,6.12711],[1.783133,-2.683408,3.373461],[3.910944,-9.39973,10.22988],[3.278206,5.724522,6.67209],[9.024782,-9.307799,13.00315],[-5.218582,-6.007776,8.02041],[-9.922637,0.9718275,10.02014],[1.338774,1.172606,2.041401],[9.999704,-0.6919439,10.07337],[-7.05016,1.017414,7.193045],[4.560174,-6.900802,8.331641],[-1.84303,-2.004161,2.90059],[8.999015,5.252398,10.46757],[-3.200114,6.322969,7.156862],[-6.804102,-1.431819,7.024664],[-6.464664,-5.204144,8.359127],[7.274015,-0.9632553,7.405347],[-2.695135,-8.431828,8.908395],[-1.406049,-5.427601,5.695246],[9.552304,-4.591383,10.64553],[4.168467,-2.778877,5.108647],[7.764115,-6.202511,9.987623],[6.151637,-4.665113,7.78498],[-7.912083,0.8083904,8.015894],[3.668775,7.283665,8.216549],[-8.057689,-2.170293,8.404553],[-0.6293378,0.4412126,1.261243],[-7.865199,-3.892821,8.832633],[5.529031,-7.950016,9.735139],[7.529834,-6.324677,9.884328],[1.222427,-8.840517,8.980483],[-2.99599,-1.362308,3.439744],[2.805335,-4.665936,5.535419],[5.860879,-6.081991,8.505323],[9.51487,0.1619051,9.568645],[3.147053,3.428732,4.760268],[7.807281,7.948546,11.18629],[6.140275,-7.819922,9.992705],[3.439135,4.902809,6.071671],[-2.432884,-2.605927,3.702672],[6.162295,-6.15423,8.766323],[4.279854,-4.92226,6.598924],[1.969283,-2.876886,3.62692],[0.4123502,5.364861,5.472821],[-8.826879,-4.462129,9.941046],[4.154095,-5.021147,6.593059],[3.627757,8.607985,9.394575],[0.7112548,3.729844,3.926528],[-1.993246,-2.53299,3.374769],[-5.544368,6.877139,8.890166],[3.851819,-0.129206,3.981608],[-0.8989518,-4.910514,5.091293],[-8.890851,-7.50208,11.67598],[4.23878,-6.034218,7.44171],[3.805643,-2.373648,4.595336],[2.823805,-1.868224,3.530458],[-2.574398,8.30294,8.75022],[-9.212195,9.744407,13.44686],[-9.554458,8.467089,12.80544],[-5.332816,0.0009258324,5.425765],[-0.9352102,9.504898,9.603005],[8.596414,-2.744575,9.079153],[-3.264433,-5.331437,6.330936],[-3.293745,-7.557235,8.304249],[-8.672374,4.857015,9.990028],[-1.314236,0.8300928,1.848316],[-7.240139,-4.592278,8.631838],[-9.172147,4.696818,10.35318],[8.774025,9.635779,13.07026],[2.92818,-2.39205,3.911027],[0.8265426,-2.424136,2.749474],[7.053493,-9.82773,12.13821],[2.448555,-1.017013,2.833679],[-5.855382,0.8942896,6.0071],[-5.791171,0.3031885,5.884691],[7.696047,-2.334458,8.104248],[-5.895906,1.339982,6.128397],[5.642587,2.923238,6.433048],[-1.476298,-9.705963,9.868393],[1.199568,-9.519218,9.646475],[-6.060663,-6.947519,9.2736],[7.246321,3.176785,7.975032],[2.540748,-9.233099,9.62837],[-9.383009,5.4703,10.90711],[3.183561,-9.36678,9.943421],[8.713938,0.4855293,8.784557],[-3.936769,-2.939817,5.014048],[2.288825,-4.950212,5.544666],[6.049231,5.801808,8.441218],[-3.310434,3.704368,5.067674],[-8.692095,9.634013,13.0141],[8.133362,9.612084,12.63106],[0.2367634,6.693104,6.771535],[-3.289825,0.8122867,3.533095],[2.159907,7.402987,7.776208],[3.632519,-6.466856,7.484345],[-7.370052,-3.641711,8.281288],[8.692006,1.778056,8.928183],[3.560568,-7.012141,7.927659],[-0.08190703,-8.172306,8.233668],[-3.67495,0.5479887,3.847798],[-0.9451037,4.791742,4.98538],[5.55002,4.012274,6.92106],[1.400135,-6.18532,6.420168],[-1.488609,-8.145972,8.341032],[-2.17092,-3.973686,4.637141],[-1.560124,5.96696,6.248088],[-4.682688,8.050237,9.366637],[-5.294898,4.340503,6.919242],[-2.170773,-6.510557,6.935388],[4.201369,1.870418,4.706375],[-5.858255,0.4665318,5.961276],[-4.277745,-8.112552,9.225649],[-5.506805,4.040856,6.903145],[-2.374068,1.03564,2.776463],[3.265788,-5.774648,6.709093],[-8.118911,-5.426979,9.816762],[-0.7242826,7.706254,7.804546],[7.914749,3.395502,8.670218],[8.893354,2.763813,9.366451],[7.144872,-3.62748,8.075134],[5.121372,-7.210237,8.900334],[-8.18878,-7.486156,11.13996],[-7.578259,-9.826997,12.44989],[-6.090753,-3.248941,6.975162],[-3.956844,5.246691,6.647133],[8.056703,8.81746,11.98574],[6.687154,-8.964364,11.22844],[5.290472,1.380838,5.558399],[3.720553,-7.417311,8.358171],[-1.965638,9.277436,9.535961],[-4.075217,-9.418243,10.31071],[-7.758079,0.8201176,7.865137],[-5.283161,-8.855816,10.36037],[6.283272,6.495422,9.092305],[-8.565834,-6.997786,11.10597],[9.785643,8.870648,13.24565],[3.527998,-9.785699,10.4502],[5.018051,-9.675211,10.94489],[9.153428,-7.888123,12.12467],[-4.045822,5.620289,6.99688],[-8.835875,0.1141629,8.893014],[-6.659587,-1.741175,6.955702],[-5.539852,-0.6598487,5.667924],[-1.904704,0.4626071,2.200433],[3.239472,-3.512659,4.8819],[-7.717862,-3.63658,8.590117],[6.68314,0.725951,6.796423],[1.892064,8.29648,8.56805],[-9.954054,-2.807493,10.39063],[2.279047,1.640579,2.980865],[-4.346849,-9.43797,10.43889],[3.346283,8.527424,9.21491],[1.97928,-8.949726,9.220366],[-5.320354,-6.095983,8.15274],[-3.523422,-9.8565,10.51499],[-3.336872,-5.108668,6.1833],[-4.722311,-5.293189,7.163662],[0.0736405,-1.856202,2.109717],[2.243383,5.463639,5.990335],[-0.711205,0.7465994,1.436392],[4.693914,-0.9250808,4.887597],[-6.137969,-0.1914783,6.221842],[-8.098494,-5.737525,9.975209],[-3.131614,3.321217,4.67306],[9.953545,7.441852,12.46813],[6.453738,4.913034,8.172432],[8.816946,9.869316,13.27185],[-5.342802,5.281577,7.578957],[8.221967,-8.2325,11.67796],[4.191467,-2.613358,5.039647],[-3.5353,-9.21453,9.919975],[1.803554,8.237839,8.492043],[-9.466272,-8.853731,12.99996],[-1.522986,3.313057,3.780983],[-8.295181,5.374096,9.934331],[-0.5882798,5.95962,6.071502],[-3.389609,-9.942723,10.55212],[-2.939474,-8.926507,9.451087],[-1.414424,-3.136452,3.583005],[-7.412049,4.524833,8.741428],[-2.014258,-0.3502819,2.275947],[6.236513,6.371331,8.971508],[-1.581704,6.24323,6.517646],[0.5131085,0.612433,1.279982],[-0.7568983,-6.707337,6.823582],[-1.393305,8.373328,8.547159],[3.101409,-8.435834,9.043343],[-3.08002,1.424913,3.537923],[8.464237,-6.593049,10.77551],[0.9276742,1.240577,1.843803],[-4.879653,7.17759,8.736636],[6.035867,-0.02667818,6.118203],[7.412883,-2.26095,7.814265],[-4.928728,7.198143,8.78098],[-8.565647,7.523532,11.44438],[6.303196,7.135712,9.573331],[-4.296556,6.112131,7.537807],[8.503992,-0.9757835,8.618007],[7.671831,-0.779765,7.775926],[8.785624,3.955292,9.686666],[2.497396,-5.616624,6.227636],[0.7682014,-1.883587,2.266723],[-3.402829,9.304497,9.957556],[2.87362,-4.376971,5.330626],[3.828948,-7.480005,8.462347],[5.714533,5.295773,7.855004],[-6.797708,9.181025,11.46735],[8.640022,-1.152795,8.773764],[0.03164914,-8.938588,8.994407],[8.499027,8.124404,11.79998],[3.160216,1.088999,3.488966],[3.928357,-8.372392,9.302094],[7.342614,5.456699,9.202692],[9.45303,2.720441,9.887395],[-3.230067,-2.176985,4.021517],[-1.342089,-8.329427,8.495914],[1.230018,4.268458,4.553316],[-8.531075,-7.460285,11.37696],[-8.198013,9.449995,12.55029],[5.159562,-0.2912545,5.26364],[-8.352279,-9.535994,12.71596],[-8.431263,6.84596,10.90657],[9.131682,-0.1572132,9.187619],[-8.313376,1.304107,8.474251],[0.4557128,-2.061645,2.336248],[1.700167,-3.839171,4.316226],[3.884583,-5.498198,6.805892],[5.773306,3.497218,6.823606],[7.41402,-4.240247,8.599267],[-2.942019,2.405438,3.929582],[-3.302036,-4.220436,5.451195],[-4.294317,-0.3947437,4.426848],[-1.241212,7.682295,7.845907],[6.907396,-4.056017,8.072386],[0.9935642,0.5196707,1.502407],[-3.637551,-3.440467,5.105741],[6.775263,-3.066444,7.503817],[1.265105,2.632757,3.087377],[-8.00841,1.299985,8.174631],[-0.2172656,7.461051,7.530901],[-9.778321,-4.549777,10.83125],[-0.6144368,6.911464,7.010411],[7.844446,-1.856752,8.122984],[-0.1265412,1.148571,1.528145],[1.0589,-9.837396,9.944628],[6.392903,-3.821232,7.514721],[3.944023,9.161726,10.0246],[-7.657671,-2.919956,8.256274],[-9.008518,-0.7284731,9.093079],[-3.868994,8.203053,9.124648],[1.618934,-6.773863,7.036062],[-0.9805614,2.878367,3.201015],[7.956067,4.648369,9.268567],[5.056073,-8.352242,9.81447],[5.104334,-6.208344,8.099244],[3.617593,-0.6865909,3.815545],[9.172071,7.334178,11.78631],[-9.912702,0.6206676,9.982328],[0.7140785,8.723207,8.809328],[-4.645825,-7.40603,8.799601],[7.955473,-3.644958,8.807683],[-8.705577,-0.5884283,8.782557],[-6.85583,1.86845,7.175898],[7.758985,1.574982,7.980126],[9.500363,5.356575,10.95216],[-7.086548,9.562607,11.94415],[6.524494,9.888138,11.88883],[-7.78012,8.737201,11.74176],[7.2475,-5.788234,9.328982],[3.910177,-7.044117,8.118441],[3.375908,9.453892,10.08825],[-8.426269,5.493369,10.10837],[-4.907058,-3.54528,6.135815],[-6.25844,7.458004,9.787231],[2.795188,1.345079,3.259189],[-3.092826,-2.883555,4.345165],[-1.881412,-2.403254,3.211751],[3.913477,-9.071575,9.930196],[-8.622274,0.2618967,8.68402],[-0.7237727,-4.273581,4.448296],[-2.745336,6.618085,7.234357],[1.608967,-6.53871,6.807607],[0.8165967,1.655607,2.099491],[2.062727,-6.638676,7.023309],[4.187735,-1.109694,4.446184],[-5.060585,3.411557,6.184517],[-4.555074,1.884863,5.03005],[-5.793154,1.701741,6.120176],[-5.690922,-6.000469,8.330199],[-1.156797,-4.580519,4.82901],[-0.7730715,-7.014665,7.127634],[-0.5024052,-8.783413,8.854421],[-9.239203,8.027947,12.28051],[1.125546,-0.3600941,1.54807],[-0.3860157,4.944257,5.059119],[1.757559,1.945758,2.806241],[1.284689,-4.376941,4.669908],[-8.001415,-2.577013,8.465439],[1.573816,1.489577,2.386575],[3.392912,-3.719109,5.132604],[-3.465774,0.5096141,3.642979],[1.421127,-2.903837,3.384061],[8.839343,6.99048,11.31374],[-4.118101,-1.985243,4.679737],[-8.086476,0.6165342,8.171365],[-1.621389,6.415367,6.692222],[6.521732,-8.601357,10.8405],[-2.816975,-8.856125,9.346994],[-7.756667,4.020779,8.793893],[-0.1591686,-8.545244,8.60503],[-4.190453,0.8239931,4.386213],[5.078227,-6.090568,7.992709],[-6.396691,-2.76756,7.041097],[9.088747,5.826882,10.84241],[-3.593327,0.7480958,3.804162],[-1.961426,4.122888,4.673906],[3.413041,-9.669045,10.30239],[-4.711914,-3.124594,5.741535],[-5.072612,5.135477,7.287284],[-3.018935,-5.932097,6.730806],[-2.230787,-3.689917,4.426274],[5.274719,7.409838,9.15032],[-7.073236,1.381904,7.27601],[2.044195,-8.662006,8.955953],[-7.896027,-2.239056,8.268047],[9.178702,1.314698,9.326146],[-3.869946,-6.842029,7.924005],[-5.868618,-8.027988,9.994462],[1.466221,0.5784342,1.866652],[-3.897302,-0.04512352,4.023804],[1.378318,-3.488064,3.881539],[1.096284,-8.56199,8.689621],[9.181385,-0.08216676,9.236048],[9.646998,9.204898,13.37142],[-3.861032,-1.322104,4.201848],[7.402659,-0.3434136,7.477787],[-5.103405,-2.59179,5.810518],[5.346535,4.985034,7.378075],[7.395342,-6.546007,9.926796],[-2.752658,9.401275,9.846883],[0.9455354,6.609884,6.751637],[8.459929,-9.250417,12.5754],[-5.667725,1.091076,5.857778],[-7.681514,5.589438,9.552354],[3.97705,-4.89619,6.386674],[0.6886155,7.557859,7.654765],[-8.172928,-7.013574,10.81605],[-4.647655,-8.440685,9.687407],[1.972289,0.5841413,2.28717],[-4.61665,7.949775,9.247291],[8.160316,-3.868616,9.086085],[5.690356,-4.399115,7.261706],[3.956909,8.911331,9.801477],[-0.1949912,4.317881,4.436453],[5.862401,7.821668,9.825794],[3.105625,-1.675085,3.667536],[-4.032102,-9.683037,10.53656],[-1.313643,9.163334,9.310873],[-3.825372,-7.171228,8.189016],[9.458878,3.039787,9.985523],[7.832903,-5.233397,9.473269],[6.392731,-1.636837,6.674297],[7.161544,-3.373244,7.979128],[4.956673,6.113239,7.933493],[8.604821,7.649539,11.55675],[-2.944657,-9.340287,9.844388],[3.42394,6.831837,7.706968],[8.656447,-7.75092,11.66237],[9.708483,3.337308,10.31466],[-4.535065,-3.024656,5.542144],[9.336043,0.3869878,9.397418],[2.933005,3.685217,4.814909],[2.71901,-1.106613,3.101227],[2.132254,4.725662,5.279999],[-6.66116,-5.696234,8.821459],[7.975079,5.153825,9.547974],[7.036418,-8.186328,10.84099],[-1.742679,8.4203,8.656695],[-0.4990797,-2.765613,2.9829],[7.775609,2.484919,8.224045],[-4.103766,5.366756,6.829566],[2.275336,-3.08777,3.963771],[5.349056,-6.786456,8.698758],[1.755314,7.147473,7.427482],[4.533585,-6.465954,7.960022],[4.773693,2.207315,5.353539],[-1.824898,-2.552813,3.293494],[-9.11576,-5.125728,10.50572],[4.807636,-8.137838,9.504619],[5.00536,-1.630387,5.358338],[-1.539052,-3.082522,3.587565],[5.812063,3.763015,6.995738],[-4.002409,-1.400848,4.356793],[3.371857,-4.558335,5.757416],[5.732923,0.06749072,5.819876],[4.851419,-6.608312,8.258696],[3.692251,-0.7638254,3.900788],[9.116572,1.798292,9.345895],[4.769442,8.844292,10.09797],[4.547272,-9.749542,10.80422],[9.863743,8.947307,13.35469],[4.352216,-8.853247,9.915733],[-9.889881,-0.5347189,9.954681],[9.309204,-9.932529,13.64978],[5.673664,-1.703194,6.007607],[-0.5956883,-4.122004,4.283195],[3.836388,6.140979,7.309548],[8.18159,6.457518,10.47081],[8.332293,8.397147,11.87178],[-4.209202,0.5073656,4.356008],[-7.775099,-6.186798,9.986423],[-0.9909362,8.003054,8.125936],[-8.425021,-1.822727,8.677749],[-9.978869,8.070058,12.87259],[-7.481596,0.1370275,7.549375],[-2.495724,6.290946,6.841392],[-1.090874,8.683802,8.808996],[8.785177,-9.453193,12.94381],[-2.183738,2.365319,3.370971],[9.343595,-4.953338,10.62254],[-8.337256,9.915389,12.99326],[6.92354,-7.58918,10.32139],[-2.601277,6.513917,7.085038],[-2.105528,-9.924904,10.19495],[6.467261,0.3042453,6.551185],[8.813169,-4.258473,9.839032],[8.997934,3.259588,9.622252],[3.9247,4.507986,6.060133],[-5.553974,1.132167,5.75573],[9.076445,-2.747355,9.535712],[-3.265733,-8.678535,9.326413],[8.652298,5.265437,10.17777],[-2.305889,9.93751,10.25043],[5.63966,9.888097,11.42717],[-3.467799,-7.871281,8.659255],[-2.520455,1.654067,3.17626],[2.402801,-0.9343385,2.76522],[-5.986443,2.140195,6.435677],[5.106097,-2.257247,5.67163],[4.336025,-1.898591,4.83795],[-9.067716,-5.574539,10.69107],[1.102418,-7.011524,7.16776],[-0.1045234,8.041994,8.104604],[-2.607809,-7.580831,8.078965],[1.627296,0.8914442,2.107787],[9.53703,2.528933,9.917179],[-7.954943,-0.9258319,8.070829],[0.9461079,5.620181,5.786324],[7.831127,-3.684137,8.712027],[-8.347431,-3.432151,9.080709],[-9.555387,-5.241915,10.94455],[4.215145,4.907137,6.545796],[5.462389,8.148331,9.860679],[9.101638,-7.973267,12.14137],[-2.201456,-3.678858,4.402318],[-3.395812,-8.970531,9.643753],[3.906039,5.375623,6.719707],[3.611562,2.168442,4.32961],[1.232855,-7.173386,7.346931],[3.978066,7.842439,8.850359],[2.780921,-2.756804,4.041471],[4.105873,-7.507646,8.615273],[-0.5040063,-3.702314,3.867964],[4.220964,3.276284,5.436044],[-1.579259,-5.149135,5.477924],[-0.5051661,-0.6325328,1.286581],[3.260526,7.675269,8.398856],[6.053866,9.304554,11.14558],[4.734195,2.156251,5.29736],[3.267867,3.622482,4.980093],[-0.6172408,0.9578243,1.516052],[-4.697692,7.119747,8.588312],[-2.688927,-8.328742,8.808988],[-3.196027,-4.617124,5.70372],[2.427577,-5.541341,6.13185],[-2.150837,9.529118,9.819887],[-1.320893,-8.228212,8.393345],[2.067204,5.898046,6.329319],[5.339082,5.582166,7.788863],[4.512345,-9.911207,10.93587],[7.45048,-1.079528,7.594408],[1.417976,-3.187948,3.629555],[-6.291715,2.334022,6.784787],[-6.080217,-5.88175,8.518452],[-9.147601,-4.219443,10.12335],[3.750812,-1.620354,4.206439],[0.3241295,-3.177555,3.346926],[-6.179671,4.504398,7.712194],[-6.08424,7.664049,9.836444],[8.316434,-9.808362,12.89833],[1.923817,-9.146397,9.399875],[-4.468796,1.134504,4.717758],[2.356532,8.779983,9.145563],[0.8285567,-6.005584,6.144391],[-6.992067,7.235173,10.11122],[5.15147,1.684553,5.511384],[9.296889,-4.551106,10.39927],[-1.847455,-6.277057,6.619255],[8.059102,0.06481592,8.121165],[-3.889967,0.2191482,4.022421],[5.410938,-3.943845,6.769946],[2.6709,-8.405908,8.876541],[9.638663,8.062624,12.60594],[7.054151,2.757605,7.639727],[-0.09691031,-6.591172,6.667304],[-9.656896,2.050215,9.922651],[-4.325148,-7.947258,9.103066],[-0.9131834,-6.018161,6.168644],[-6.497045,-0.3163068,6.581158],[5.372506,2.394362,5.966305],[5.930819,8.248613,10.20854],[5.689508,-4.188765,7.135562],[-2.567661,7.567753,8.053804],[5.127347,-4.922411,7.177731],[8.42156,-2.324631,8.793553],[2.029898,2.409332,3.305354],[-8.415689,-4.036847,9.387223],[-5.967998,-7.034938,9.279405],[8.781093,-2.130947,9.091124],[4.463771,6.290485,7.777883],[0.08966027,9.374501,9.428113],[1.150754,-0.3150168,1.55675],[7.970949,-7.072564,10.70314],[5.182117,5.436411,7.576867],[-6.135428,7.959688,10.09951],[2.610349,-9.97078,10.35521],[7.076854,1.20478,7.24799],[-6.481087,-1.162195,6.659968],[1.895936,7.711941,8.004287],[-6.667254,4.935257,8.35518],[-5.296018,2.445245,5.918364],[-0.5569239,-8.435905,8.513205],[0.1368878,8.119997,8.182487],[-3.482261,-3.570922,5.087005],[3.966849,-3.264455,5.23379],[-2.449784,2.049441,3.346886],[4.240481,-7.701833,8.848724],[0.6043187,7.837348,7.923965],[-8.278341,8.984646,12.25785],[-6.161838,4.471349,7.67862],[-5.487699,5.878633,8.103898],[2.370713,8.315944,8.704896],[-3.844634,8.736897,9.597634],[6.007084,5.301047,8.073794],[7.791266,-8.603666,11.65019],[6.70668,-6.936698,9.700378],[1.358993,-4.007942,4.348616],[5.277054,8.271118,9.861982],[-4.599102,-5.383015,7.150426],[7.969011,-4.946503,9.432551],[-7.164047,-3.292398,7.947544],[-1.916175,-6.108802,6.479907],[9.444728,-3.534904,10.13402],[-6.236247,-0.7823165,6.36418],[-1.491895,-7.745686,7.951188],[2.211836,-3.187731,4.006725],[1.084832,9.602946,9.715629],[-7.796807,6.173684,9.995227],[-7.631226,-7.459112,10.71793],[9.66765,3.641966,10.37918],[-7.335021,3.068607,8.013669],[-1.090342,1.170644,1.886599],[-7.582307,6.215689,9.855261],[3.180058,-0.6471342,3.395813],[2.933642,4.250311,5.260361],[-2.732233,-1.409439,3.232896],[2.301126,3.853234,4.598107],[-3.981205,6.310432,7.528051],[0.8219969,0.3833369,1.350047],[9.932079,-1.64114,10.1163],[8.292449,-6.574603,10.62968],[8.287403,4.382566,9.42804],[-9.848344,0.9773747,9.947118],[9.762671,-3.059087,10.27948],[9.74319,3.064892,10.26271],[7.743554,-9.675999,12.43332],[-3.020437,3.808021,4.962264],[-9.522521,3.212657,10.09948],[-4.564081,-7.022442,8.434781],[-1.069704,-1.786735,2.310127],[5.018824,-2.62077,5.749524],[0.1771266,1.236532,1.60012],[-8.939567,0.03542064,8.995394],[8.677719,-1.215152,8.819263],[7.328279,-7.001496,10.18453],[-4.281662,-7.557171,8.743196],[-7.810601,-1.636136,8.042539],[-6.456265,5.897349,8.801255],[-7.107834,-9.0732,11.56911],[4.959996,-1.091613,5.176212],[-5.300412,6.672563,8.580062],[-9.575088,8.472,12.82408],[5.361254,8.754495,10.31427],[-4.373386,-2.996463,5.394932],[-3.531559,9.871614,10.53189],[5.410559,-7.411525,9.230647],[-4.617996,2.156578,5.193912],[-8.544888,-1.239986,8.692105],[-1.830522,8.194637,8.455938],[-9.675075,-1.651282,9.865789],[3.906553,1.44596,4.283918],[-9.606283,-2.653954,10.01619],[3.768788,9.509815,10.27815],[-4.01157,6.817894,7.973479],[-8.408295,-8.283887,11.84577],[-9.262033,5.830087,10.98977],[2.604995,-7.574852,8.072447],[-5.482109,-5.821367,8.05865],[5.420883,-8.258027,9.928796],[7.201836,5.135153,8.901474],[3.350782,-6.56102,7.434697],[-0.1416261,7.293512,7.363109],[-7.454583,-5.13459,9.106855],[-0.3547255,8.396651,8.463427],[-5.425512,5.441562,7.748986],[9.317666,9.670343,13.46605],[7.863175,-2.718521,8.37973],[2.88214,2.472476,3.926814],[3.309585,-2.386222,4.200882],[8.763013,8.741624,12.41799],[-7.972623,-1.486765,8.171486],[-2.923015,6.109588,6.846246],[-0.03019246,-3.511232,3.650981],[1.793492,-7.185418,7.473074],[-5.713451,2.561077,6.340555],[-5.431937,-2.668948,6.134266],[9.065042,-3.073666,9.624055],[6.120029,9.066875,10.98467],[-9.213751,-3.279983,9.831149],[-4.052538,5.73931,7.096672],[-1.886525,2.995689,3.67874],[0.9600481,-6.303245,6.453882],[-3.295546,-5.545703,6.52805],[8.648563,0.4917242,8.720059],[9.667921,-4.706967,10.79927],[-6.913043,-0.8779425,7.039954],[3.22879,8.885417,9.506614],[-3.586608,-9.519876,10.22212],[1.552099,3.108614,3.61559],[-4.85323,-1.697058,5.237733],[-0.3494792,-1.275768,1.658228],[-9.065651,7.923227,12.08154],[6.161022,-2.578762,6.753385],[1.517751,-1.669029,2.467636],[5.015759,8.398524,9.833262],[9.134568,1.402121,9.295498],[-5.784668,-1.840547,6.152236],[-6.487428,4.764674,8.111032],[9.721638,-1.810348,9.939196],[7.995433,2.74288,8.511777],[6.96268,-0.8029315,7.079803],[-7.067807,-4.134576,8.249158],[0.08312746,0.3273996,1.05551],[6.319187,6.652804,9.229948],[2.744331,-9.593898,10.02867],[6.340886,-5.302591,8.326122],[-6.533237,-9.609505,11.66301],[-7.847186,3.975194,8.853276],[9.276086,-9.42109,13.25906],[-2.567999,-2.066199,3.444386],[1.507167,5.727664,6.006471],[-4.377569,-2.875891,5.332341],[7.277967,-1.319549,7.463914],[-2.462607,-7.705773,8.15128],[-7.871677,-5.854197,9.860778],[3.859045,-6.647035,7.750825],[1.396679,7.714634,7.903562],[-2.547669,-2.898422,3.986411],[-5.288947,5.996678,8.058108],[-6.080444,7.658011,9.829391],[0.05636749,7.471814,7.538646],[-8.501192,9.417291,12.72618],[-0.7431832,0.2273937,1.266503],[-3.544834,1.713762,4.062367],[-6.611368,-9.969915,12.00456],[3.658901,-6.477952,7.506758],[-1.937926,3.586168,4.197161],[7.396028,-4.337619,8.632275],[-6.790079,-8.161118,10.66344],[4.299419,-9.880243,10.82147],[4.858232,-8.468619,9.814271],[-0.6557209,4.207897,4.374513],[6.667487,8.118443,10.55294],[0.5836998,1.159444,1.638601],[-6.488178,-6.974902,9.578398],[2.494824,-2.234003,3.494985],[3.418528,-5.970246,6.951991],[-5.434537,5.846418,8.044551],[-2.984448,1.861628,3.656855],[8.471598,0.356605,8.537865],[1.022493,-9.687049,9.792058],[3.30419,-9.289217,9.909955],[4.736517,8.961623,10.18554],[4.489333,5.908398,7.487541],[7.227281,-6.272296,9.621605],[0.1291775,-6.084775,6.167753],[-3.429725,-1.406945,3.839597],[5.127396,3.732782,6.42058],[-1.607248,-3.821559,4.264688],[1.977532,5.87404,6.278134],[-3.763825,7.213028,8.197205],[8.385225,9.06651,12.39006],[-3.180408,8.625345,9.247247],[-6.858607,-2.488442,7.364295],[0.0793945,2.436035,2.634496],[9.708916,-2.579351,10.09535],[-0.3738451,3.868665,4.013269],[1.51049,-7.787713,7.995627],[-2.040331,1.895448,2.958999],[-9.041863,6.453693,11.15372],[3.129939,-4.861234,5.867547],[0.9626269,4.778871,4.976369],[-2.750983,-2.521036,3.863099],[-2.193401,8.803934,9.127994],[5.204125,7.932351,9.53966],[1.154904,8.249766,8.39002],[4.943279,-4.522977,6.774462],[2.000416,7.748052,8.064364],[0.7200015,-9.387012,9.467544],[-3.897451,-4.355172,5.929388],[3.905285,5.09567,6.497469],[-8.482476,0.3613054,8.548857],[8.510964,6.877924,10.98828],[-0.743115,0.9819769,1.586348],[8.531843,-5.302313,10.09489],[6.217713,4.946216,8.007809],[1.633975,4.232572,4.645917],[-4.717134,-6.524675,8.113122],[-6.248425,-6.603447,9.145947],[4.59359,6.371948,7.918509],[-1.281706,-3.075392,3.478622],[5.410211,4.201266,6.9225],[-0.1705726,7.193582,7.264758],[-1.385369,9.24435,9.400918],[-4.40051,-1.268039,4.687474],[-9.958926,-4.981251,11.18003],[8.818188,-1.761084,9.047753],[-7.385609,5.140891,9.054059],[-3.32,-1.631031,3.831796],[-7.729207,-4.905153,9.208755],[-2.461737,2.965813,3.981984],[9.031911,-3.402581,9.703245],[-8.864042,-3.734098,9.670301],[-9.782412,5.585088,11.30879],[8.7902,8.92574,12.56728],[-9.818149,5.329083,11.21585],[8.925853,3.06943,9.491693],[-8.117188,8.771082,11.99252],[-9.397514,-4.039348,10.27763],[6.94065,-9.412771,11.73767],[-5.480756,-6.987105,8.936349],[3.096016,8.265113,8.882421],[8.752796,4.648306,9.960833],[-3.05083,-1.48056,3.535481],[-8.737484,-6.986384,11.23179],[-2.19229,2.750065,3.656363],[-1.080054,-6.289668,6.459601],[-2.671841,-2.975795,4.122388],[4.840734,6.740364,8.358541],[-2.715479,2.035814,3.53813],[9.194938,-0.0387273,9.249236],[-9.348006,-6.445932,11.39892],[2.019125,7.586775,7.914292],[2.274336,-1.356365,2.830606],[-1.594979,8.413983,8.622011],[-7.005465,-4.985167,8.656121],[0.9863861,-5.400437,5.580114],[-0.2337728,0.7345783,1.262638],[-8.779316,2.177715,9.100486],[1.630014,-6.452364,6.729781],[-3.666707,-4.291912,5.732822],[7.657179,-9.014653,11.86998],[9.951585,-8.479381,13.11236],[-1.216702,6.454689,6.644048],[-3.417977,1.30143,3.791607],[1.323982,3.030725,3.455173],[-3.007133,-1.294123,3.423099],[2.046302,6.532696,6.918343],[-2.579396,3.102377,4.156684],[9.671131,-5.627823,11.23402],[-0.4459549,-3.59103,3.754247],[-5.913123,-0.6204882,6.029099],[5.300447,-6.734713,8.628505],[6.883778,-6.678185,9.64285],[-6.016191,-1.448033,6.268282],[-1.129391,4.797773,5.02933],[-1.562144,-2.601199,3.194766],[-8.516914,5.644536,10.26638],[8.637904,-4.944825,10.00323],[1.84358,-5.70316,6.07658],[-0.3087837,4.514617,4.63434],[0.9342974,-9.003939,9.107349],[-2.630512,0.168294,2.819205],[-0.7690254,8.927218,9.015909],[-1.017821,9.267973,9.377168],[1.557907,9.513385,9.69183],[-4.86746,5.224019,7.209891],[4.000099,0.09461228,4.124287],[-8.753507,-1.648803,8.963394],[-6.55331,3.447573,7.472057],[-0.429565,1.317621,1.708991],[-8.199811,3.487483,8.966573],[-1.661065,-0.3228753,1.96555],[-8.716265,3.264412,9.361071],[2.647088,1.676738,3.289152],[9.638319,7.559939,12.29024],[-0.8764858,-7.066384,7.190411],[9.429744,0.7529694,9.512466],[9.847383,0.9298601,9.941609],[6.25801,-3.090954,7.051006],[2.462018,0.5678892,2.717358],[-3.838681,9.902414,10.66739],[-5.681934,8.464965,10.24402],[7.963779,-2.552094,8.42229],[-7.913896,-2.666478,8.410699],[6.1714,1.13241,6.353623],[2.744889,-2.826753,4.065089],[-7.316671,-9.798336,12.26952],[-8.941626,9.025511,12.74412],[8.116389,-6.092865,10.19798],[1.540646,-1.394418,2.306077],[8.462674,-1.353776,8.628416],[-5.331528,-4.361621,6.960526],[-5.127351,-7.991109,9.547123],[-8.982678,-7.297693,11.61658],[5.724152,0.6507374,5.847167],[-6.261768,-4.235706,7.625677],[-4.329756,9.172103,10.19187],[-6.218832,0.7122765,6.338865],[2.125359,4.520913,5.094684],[9.453441,-5.808335,11.14021],[-7.354702,7.084458,10.26066],[-2.505737,7.667425,8.12823],[2.32421,9.505895,9.836869],[3.741202,-4.286661,5.776855],[3.16271,2.165846,3.961518],[-9.434541,-8.305854,12.60943],[6.499024,4.276935,7.844074],[2.17498,-2.170666,3.231459],[2.405977,9.225812,9.586675],[-9.409488,-6.607799,11.54129],[-4.317581,9.747323,10.70756],[-2.293085,-5.737193,6.258883],[3.757021,2.519901,4.633046],[2.855282,-6.063993,6.776772],[8.343461,-6.098679,10.38303],[7.955764,-9.258489,12.24801],[-4.024394,-0.3841014,4.164526],[-7.083542,6.918124,9.951735],[4.833675,8.394437,9.73812],[5.05583,-2.940102,5.933432],[-1.082645,6.879809,7.035901],[-5.895067,3.754334,7.06023],[1.185799,-4.713655,4.962325],[-4.070023,-6.788394,7.977931],[8.315164,5.564547,10.05515],[-3.799644,-5.545042,6.795938],[-0.5556321,2.717601,2.948573],[0.697813,-9.862013,9.937114],[9.569399,4.85715,10.778],[6.351316,-7.730144,10.05457],[-0.2588238,7.303119,7.375807],[9.396074,7.172751,11.86316],[-8.984544,7.629799,11.82945],[-0.9248022,-6.441074,6.583517],[-1.947302,-5.151423,5.597245],[9.733056,-0.04013788,9.784374],[-5.630575,0.2609905,5.724639],[-8.08064,-0.9857324,8.201732],[5.64075,-1.339162,5.883147],[9.004729,5.28183,10.48727],[7.214299,-4.113894,8.364821],[9.563101,8.317739,12.71368],[3.910252,-0.5536227,4.073888],[-9.105172,-2.409622,9.47156],[4.256504,-6.586244,7.905468],[3.356564,-7.598704,8.367008],[-5.277678,3.32234,6.315998],[5.948592,-9.703931,11.42594],[1.619027,6.061372,6.353069],[-9.965819,-3.541204,10.62345],[-3.924944,8.475065,9.393185],[-4.330619,9.054907,10.0869],[-6.263229,-3.795506,7.391475],[-9.186654,2.210999,9.501743],[8.500083,-6.186438,10.56047],[-7.469251,1.400294,7.66489],[-9.455531,-4.467549,10.50552],[-4.28089,-7.482431,8.678294],[0.9238409,-3.029005,3.320896],[-7.879797,6.760355,10.43042],[1.082923,-9.571944,9.684773],[-7.324514,-1.56565,7.556439],[0.6990712,-5.407317,5.543264],[-1.723192,2.579258,3.259135],[-1.119988,-2.531333,2.943131],[7.441556,-1.44042,7.645362],[0.2856756,3.33345,3.491919],[-6.832502,-7.796116,10.41453],[8.321278,-8.687914,12.0716],[5.257184,8.648355,10.17015],[7.909396,-8.182192,11.42396],[6.925184,-3.29527,7.734144],[-6.423004,2.243521,6.876654],[7.252946,-5.434628,9.118137],[3.481032,-1.259807,3.83467],[7.387582,5.195416,9.086733],[-0.7358399,-5.006735,5.158377],[8.238169,3.628763,9.057336],[-7.849524,-7.095903,10.62859],[0.4084393,-2.939633,3.131815],[2.475915,-2.223273,3.474637],[2.484237,-4.228388,5.005067],[2.672413,0.7210938,2.943088],[5.559503,-8.497533,10.20373],[7.05869,7.164809,10.1074],[5.343874,1.415699,5.617935],[4.948767,1.932361,5.405951],[-8.889015,-9.877321,13.32577],[0.3097768,-5.704572,5.799836],[2.855456,0.6082838,3.086039],[-7.777306,5.970829,9.855825],[6.339019,2.823967,7.011273],[4.38928,2.375654,5.090138],[-0.8851867,-8.24357,8.351048],[-4.935272,-7.361081,8.918656],[-0.3454079,3.016747,3.196884],[-8.562421,-2.557631,8.992025],[5.967784,-3.541412,7.011138],[-5.973212,4.921917,7.804135],[-3.101011,-1.373932,3.536093],[-0.8852762,2.235829,2.604351],[-6.651644,0.2874828,6.732534],[-5.439477,-2.336606,6.003969],[-9.93841,9.883267,14.05172],[7.775953,-5.284043,9.454447],[3.118805,-5.941008,6.78399],[-6.627821,-3.398061,7.514974],[6.410026,-3.302818,7.279906],[9.532999,5.566576,11.08444],[4.654861,0.1477374,4.763356],[-1.353864,-7.756852,7.937361],[8.615736,8.607346,12.21955],[-8.804266,5.553329,10.45727],[5.190955,0.4593525,5.306318],[-0.8879007,7.361665,7.482144],[-6.94941,2.677927,7.514359],[4.618527,5.386382,7.165465],[4.805251,7.843675,9.252766],[0.4806498,2.432758,2.673824],[3.177994,-9.376887,9.951163],[-4.120412,5.589307,7.015565],[5.537644,8.3939,10.1056],[7.030504,-0.5521873,7.122703],[-3.002629,-0.6751144,3.235979],[3.615765,6.796436,7.763073],[1.177503,3.677153,3.98848],[-0.9363456,0.1661804,1.379985],[-2.195665,-8.031215,8.385783],[-9.185344,4.001422,10.06886],[-2.214412,5.284571,5.816383],[-9.868307,-3.693821,10.58432],[1.210772,-0.5270593,1.65643],[6.906159,-3.512191,7.812203],[-5.051567,0.9651259,5.239255],[-6.390922,-4.823406,8.069022],[8.361105,-4.991146,9.78875],[9.867305,-8.722146,13.20756],[-8.018476,-7.034538,10.71358],[-8.334758,-2.553138,8.774207],[-3.279745,8.879598,9.518613],[8.280025,9.32918,12.51369],[8.996838,2.416065,9.369123],[9.808517,-6.212262,11.65329],[-8.899082,7.206081,11.4944],[4.579411,-2.002954,5.097335],[-2.257466,-8.840346,9.178664],[5.661117,-4.171127,7.102573],[-9.706571,3.015545,10.21328],[9.247133,-2.036185,9.521319],[2.182299,9.573478,9.869848],[-9.361732,5.658043,10.98433],[9.290874,0.9936484,9.397216],[-9.304005,-6.980574,11.67446],[-8.858253,2.869965,9.365113],[6.853761,-5.708704,8.975709],[-9.56776,-8.77217,13.01895],[-4.722635,6.142035,7.812035],[-3.116193,-8.241952,8.867944],[-1.397599,-2.523278,3.052902],[2.346619,1.314062,2.869386],[-8.473837,-4.406491,9.603285],[1.900305,0.414188,2.186941],[5.424891,2.682921,6.134126],[9.547637,-2.248614,9.859697],[5.486035,6.287983,8.404482],[0.4940209,2.062742,2.344987],[-7.433337,-0.9041066,7.554595],[0.3947026,0.8626363,1.37838],[-3.759305,8.543415,9.387349],[8.295027,-0.5868743,8.375672],[8.455276,8.332938,11.91342],[4.394027,-7.110199,8.417981],[6.650087,5.329773,8.5808],[7.677896,7.868862,11.03943],[6.472447,9.353907,11.41876],[9.21455,-6.979062,11.60238],[-1.798314,8.037925,8.297117],[-1.914138,6.981771,7.308149],[9.024523,9.476378,13.12417],[-8.071795,8.306405,11.62541],[-1.077911,1.995729,2.478876],[-7.155457,-8.362831,11.05158],[-1.594652,-7.542224,7.773549],[6.378462,3.086607,7.15625],[8.818546,7.570209,11.66511],[5.489308,9.402327,10.93326],[-5.106024,-5.434169,7.523408],[8.750553,4.735037,9.999639],[8.284713,1.36708,8.456084],[-9.045237,-7.759229,11.95918],[9.025335,0.9987243,9.135324],[3.472428,5.418213,6.512664],[-2.210507,-9.815309,10.11072],[-8.794622,-8.910691,12.55969],[-6.100082,-7.784616,9.940385],[6.977111,0.06723807,7.04873],[9.907228,-8.108665,12.84148],[-9.012413,-8.793814,12.6315],[-2.468261,-2.407471,3.590018],[-2.934753,6.199063,6.931173],[2.613584,-5.414074,6.094507],[4.666449,1.776078,5.09217],[9.263367,7.876267,12.20023],[6.877477,-9.865987,12.06803],[4.711473,-2.280931,5.329223],[7.243971,-5.148477,8.943262],[-7.087565,0.719841,7.193869],[8.529179,-2.523646,8.950736],[-7.297935,4.273604,8.516075],[0.7457405,5.071621,5.222784],[-8.039132,7.635794,11.13252],[-2.51892,-2.633455,3.778895],[5.692453,3.141192,6.578078],[1.057502,-4.921965,5.132646],[-5.773711,-0.6523147,5.895867],[7.179937,-0.68894,7.281904],[-0.1919281,-2.592232,2.78505],[0.1538497,-0.2567036,1.043823],[0.03407677,-5.107383,5.204471],[-5.008235,-2.108597,5.52527],[-1.714579,-7.740038,7.990492],[2.250948,-5.439478,5.971154],[6.544402,-9.318598,11.4309],[-7.483016,-3.676786,8.397279],[9.136044,-0.8862554,9.233241],[7.460351,-3.955342,8.503033],[-9.107127,7.357878,11.75067],[-5.327539,-9.462308,10.90495],[-7.006706,-8.42065,11.00006],[-0.7100271,-3.661831,3.861754],[0.3741361,-5.086576,5.197425],[0.958651,-5.518213,5.689437],[5.427458,-2.096516,5.903615],[-4.576718,0.4855451,4.709788],[-2.346684,-0.09556466,2.552657],[8.365328,1.285214,8.522351],[5.09964,7.815167,9.385262],[2.974661,-0.2812568,3.150827],[-5.771164,-8.110549,10.00437],[0.7960012,1.23444,1.776924],[7.917594,-2.997341,8.524807],[8.270347,-1.308018,8.432648],[-5.954382,2.308191,6.463932],[1.307313,-1.271762,2.080011],[4.988663,3.867118,6.390724],[3.855012,-1.397443,4.220659],[4.124998,0.6837839,4.299206],[-5.439649,4.295053,7.002661],[-9.270196,-1.217458,9.403124],[9.105769,3.48705,9.801763],[9.65147,-1.255855,9.784072],[-3.063665,4.793519,5.776146],[-9.83349,-6.446624,11.8007],[0.4257514,-3.380692,3.551104],[-0.3831196,-5.636263,5.737094],[-7.355731,-5.95379,9.516007],[2.531679,7.629434,8.100473],[-0.4568526,5.181347,5.296704],[7.337766,1.873867,7.638991],[-7.102428,1.117481,7.259011],[-7.883779,6.161119,10.05551],[-8.987316,-3.924709,9.857748],[-0.515686,3.337099,3.52167],[8.071041,-4.655997,9.371233],[4.276924,4.736269,6.459437],[2.275834,8.343608,8.706045],[1.366447,7.73754,7.92065],[-4.559934,4.329689,6.36704],[-8.447812,-2.664625,8.914356],[6.323609,1.35744,6.544515],[-4.161082,4.542749,6.241087],[7.365647,6.737235,10.0321],[-1.817082,-2.847834,3.523059],[9.727248,2.873117,10.19187],[-5.960362,-9.897388,11.59673],[4.036988,-2.478403,4.841462],[-9.250761,-4.856565,10.49585],[-1.132484,2.118851,2.602316],[-8.344337,2.478607,8.761932],[5.077918,-5.682603,7.686171],[-8.413991,-7.507917,11.32096],[-0.7363116,-1.909726,2.277983],[-3.009999,-4.644807,5.62444],[-5.230474,-8.534769,10.05983],[6.851517,3.140395,7.602984],[-7.53136,-2.349385,7.952421],[-0.4390952,-5.659335,5.763756],[8.499512,3.481005,9.238998],[-3.196815,-7.186887,7.929122],[-3.130087,8.023958,8.670718],[3.970209,-3.335944,5.281201],[-8.186887,1.03436,8.312341],[-4.695768,-3.577466,5.987361],[-6.971942,-8.455154,11.00444],[1.778772,-4.321733,4.779268],[-2.702448,-0.6580647,2.955719],[-4.201855,8.575105,9.601459],[-2.191721,3.359065,4.133638],[4.384438,6.808509,8.1596],[-6.455182,-9.214353,11.29485],[-8.050876,2.949149,8.632154],[2.574652,-5.308305,5.98389],[0.8364726,2.764832,3.056793],[-8.078938,4.322987,9.217237],[3.041436,-8.48842,9.072134],[3.700261,7.359215,8.297588],[3.840542,-9.143258,9.967394],[-9.288502,-6.977803,11.66045],[5.543,3.40953,6.584053],[7.562519,9.385516,12.09461],[1.995383,-7.548348,7.871411],[0.1607402,8.672624,8.731565],[-7.317855,8.158545,11.00513],[-9.871163,-6.564285,11.89663],[-0.5437045,8.12111,8.20049],[-6.555251,-9.907906,11.92216],[-5.179623,6.540651,8.402893],[1.506051,-4.499685,4.849263],[-6.442364,-9.336809,11.38772],[-3.121041,-1.875009,3.775785],[6.176118,-0.7271423,6.298664],[4.145356,-3.882703,5.767093],[5.262579,-8.86446,10.35729],[1.039092,-3.858583,4.119269],[-2.30193,-3.201142,4.067701],[-9.220319,-4.424801,10.27585],[5.756808,1.359505,5.999091],[-8.91712,-5.040965,10.29205],[2.84891,-4.409703,5.344321],[-3.010178,4.648291,5.627413],[4.263401,-1.404338,4.598777],[-3.723753,4.29454,5.771431],[1.067071,8.691072,8.81325],[-3.89783,-0.8290271,4.108572],[-2.487268,-2.863611,3.922598],[-6.660068,5.747478,8.853812],[5.073587,-0.1774408,5.174242],[-9.616675,6.026206,11.39279],[8.820393,5.628447,10.51089],[2.569104,8.063481,8.521739],[-7.297857,2.993166,7.95096],[-3.803301,4.360104,5.871593],[8.732265,-3.495281,9.458828],[7.745537,1.674449,7.98731],[-9.834803,-0.8376865,9.920941],[-1.061328,-5.705637,5.889033],[-9.651464,-1.91869,9.891012],[-3.171285,9.903316,10.44666],[0.0563279,6.719602,6.793838],[9.395549,-2.880148,9.877833],[-9.020248,-8.094055,12.16053],[9.686136,5.317018,11.09468],[5.865831,-6.701387,8.961951],[-0.058093,5.389501,5.481796],[2.435114,9.618417,9.972147],[-9.568503,-2.100797,9.847315],[-8.739391,5.114414,10.17517],[-1.405107,-0.5090637,1.798185],[-4.989345,-8.920851,10.27011],[-1.671831,2.586332,3.237921],[-9.600349,-9.539162,13.57064],[4.902878,8.550792,9.907283],[7.549505,8.80075,11.63822],[6.146396,0.8615125,6.286524],[-7.67172,-1.63717,7.907946],[1.296679,5.000852,5.262119],[-8.552695,7.066598,11.13936],[-2.205918,-2.941158,3.81005],[8.983545,-6.291642,11.01312],[-7.524231,-2.989878,8.158029],[-8.577415,9.709098,12.99379],[-7.745353,-3.269469,8.466399],[-4.024871,-5.874642,7.191037],[-3.883561,4.448117,5.988972],[0.2720359,-3.040547,3.212309],[-0.2212987,9.267292,9.323715],[-3.274174,-2.113061,4.023088],[-6.614052,-8.68004,10.9585],[-5.657049,5.461384,7.926469],[-3.313859,-4.48102,5.662261],[1.382255,4.425842,4.74328],[-1.179213,5.114813,5.343394],[5.542331,-4.534362,7.230344],[-7.995508,-8.392656,11.63464],[8.871778,-2.000868,9.149422],[-7.521232,-1.386672,7.713092],[1.171833,6.72613,6.900291],[-8.170115,4.82053,9.538778],[3.359535,8.777493,9.4515],[-2.30741,-8.876643,9.225992],[3.150689,-5.484211,6.403391],[-4.27342,8.376751,9.456853],[8.452131,3.715903,9.286897],[-4.573811,6.525752,8.031511],[4.730932,-4.457972,6.576871],[3.958866,9.707193,10.53101],[0.1436155,2.801457,2.978051],[4.475457,-1.030222,4.700115],[8.353098,5.265168,9.924527],[4.056514,4.313087,6.004833],[-9.540756,-3.528484,10.22136],[4.457321,8.237097,9.418995],[-5.725309,9.117539,10.81243],[-8.418668,0.1974885,8.480151],[3.153896,0.3782662,3.330187],[-7.570054,-4.379375,8.802536],[-1.707923,-6.775105,7.058261],[-0.1217515,4.211831,4.330628],[6.060525,-7.78697,9.918007],[0.6945292,5.137134,5.279443],[4.200402,8.38931,9.435248],[-9.819388,-7.282794,12.26619],[7.560737,-2.054918,7.898571],[-5.710106,3.833735,6.950024],[-2.280266,3.984998,4.698917],[1.779105,9.735011,9.94664],[9.709916,-9.917535,13.91546],[9.813367,-4.770545,10.9572],[-4.358712,-8.01297,9.176386],[2.477756,8.376047,8.791896],[6.396301,-9.453664,11.45794],[-6.741276,-7.237201,9.940919],[8.94155,-8.540244,12.40512],[5.702064,-0.6252913,5.822759],[7.848375,-1.457853,8.045018],[5.084004,8.064917,9.585926],[-7.603934,9.789922,12.43633],[-0.4134541,1.000094,1.473476],[-6.675208,-6.295908,9.230215],[-5.096969,-9.818861,11.10807],[7.001001,-8.881346,11.35308],[5.881155,9.524705,11.23868],[-2.980885,0.6924743,3.219503],[1.042603,-8.80103,8.918809],[0.3174655,5.015146,5.123716],[8.959971,1.257269,9.102846],[-5.469697,6.930364,8.885242],[9.154168,7.294065,11.74743],[5.526065,-4.839615,7.413452],[8.117255,-8.715508,11.95198],[9.461662,8.866683,13.00543],[-2.414317,6.040425,6.581464],[0.05368252,9.735518,9.786889],[-9.349236,5.23558,10.76195],[-3.267764,1.324418,3.665019],[-3.894006,4.664481,6.157976],[1.267164,-6.683411,6.875586],[-0.2298904,-0.6323171,1.205269],[7.076639,8.347242,10.98887],[-1.885097,7.76726,8.055056],[-6.502326,9.749555,11.76155],[3.678898,-2.589121,4.608453],[-4.798928,-5.144135,7.105761],[3.03733,5.459948,6.327433],[9.670519,5.460569,11.15064],[3.877644,-4.350502,5.912951],[-6.437234,8.539247,10.74042],[1.069013,1.821746,2.336996],[-9.760849,0.1468956,9.81304],[3.611541,-9.904011,10.58927],[1.989929,9.131886,9.399529],[6.723576,9.037912,11.30886],[7.610737,-0.08751834,7.676651],[-2.140759,-3.664471,4.360183],[-2.159015,9.968084,10.24812],[-5.676472,9.424576,11.0474],[6.678045,-8.747154,11.05029],[-8.498635,2.365088,8.878088],[-4.073749,0.2831024,4.204233],[5.991421,8.64347,10.56441],[4.383621,6.577181,7.967147],[1.802249,6.228503,6.560667],[8.890569,3.798294,9.719529],[3.995889,-5.661216,7.001178],[9.778047,-5.933062,11.48092],[-5.48572,3.96854,6.844153],[-6.530749,-1.483261,6.771318],[6.606818,7.501803,10.04625],[5.417221,0.3155325,5.517776],[8.675779,8.12309,11.92702],[7.085883,4.390355,8.395532],[-7.256907,7.879449,10.75864],[-8.343482,-2.234293,8.695157],[6.697455,-1.52475,6.941237],[0.006046193,-5.435939,5.527157],[-0.6282765,-0.243369,1.205803],[-4.192049,9.323976,10.2718],[-4.070932,3.466961,5.439881],[-3.068352,7.229818,7.917389],[1.888747,-9.863667,10.09254],[-2.247804,-1.520722,2.892269],[-5.839566,7.98635,9.943958],[6.335563,9.79397,11.70731],[4.901665,-0.9386669,5.089932],[8.54783,-1.597589,8.753153],[2.812118,3.535459,4.626821],[-0.09436949,9.8448,9.895908],[-2.320908,1.299834,2.841862],[-7.114805,5.868664,9.276943],[4.727992,-6.032749,7.729681],[-0.8061469,-8.341537,8.439853],[-7.19957,-3.346104,8.001888],[9.983344,-2.989781,10.46929],[2.772044,-5.211648,5.987112],[2.986814,2.100066,3.785674],[-8.301259,-5.242431,9.868839],[9.180977,8.360816,12.45767],[2.217172,-9.364913,9.675611],[-7.042462,-9.893325,12.185],[4.494464,-3.126139,5.565335],[-9.062602,5.468362,10.63173],[9.814893,-7.994506,12.6982],[-3.322364,1.395748,3.739815],[-3.028989,2.933724,4.333764],[1.817012,-1.796409,2.743833],[8.227825,-4.385048,9.376873],[-2.050683,-2.450463,3.348144],[8.524322,-5.665154,10.28387],[-8.127131,-2.524321,8.56869],[0.9374768,-5.218096,5.395126],[-1.725521,-5.143566,5.516674],[-8.373802,5.612187,10.13002],[5.335165,-6.947776,8.816778],[0.3416752,-6.444231,6.530303],[5.866871,-0.3149542,5.959813],[2.089319,0.3375451,2.340767],[4.749691,-8.752966,10.00869],[1.776799,-8.748232,8.982681],[-1.621284,-4.996727,5.347507],[-1.683876,-7.576716,7.825731],[2.516871,4.978923,5.667832],[-1.767169,-4.030402,4.512984],[1.014109,-2.238785,2.653408],[-0.6826468,9.000643,9.081717],[-0.1900559,7.803774,7.86988],[6.085577,8.89895,10.82708],[4.495818,0.8640206,4.686033],[5.196871,-8.063067,9.644715],[-9.730252,-7.007998,12.03286],[1.798787,-7.877105,8.141524],[9.463124,-9.206002,13.24014],[1.467239,4.475827,4.815165],[-9.050193,5.327938,10.54955],[6.665246,2.402804,7.155346],[-8.013946,-3.918139,8.976365],[-6.486549,-7.495129,9.962543],[7.926046,4.590718,9.213951],[1.415715,-3.383583,3.801695],[1.627648,-3.870811,4.316528],[-7.144736,3.842362,8.1738],[7.02283,-5.275131,8.840088],[-6.722955,-5.197782,8.556581],[-2.561205,2.655481,3.82248],[0.2445347,-0.6587573,1.222194],[-7.60045,7.513473,10.73402],[2.932808,9.970504,10.4409],[-0.4193568,2.539734,2.761541],[2.763345,-8.258009,8.765317],[-7.818409,-1.436448,8.011923],[2.585891,-7.482489,7.979629],[7.946949,-2.699323,8.452238],[-9.423691,9.352599,13.31454],[-3.776657,7.28552,8.266918],[6.987598,-0.461293,7.073848],[5.485225,-6.808428,8.800135],[-7.51793,6.473206,9.971042],[3.347569,-3.142716,4.699243],[0.3795152,7.255478,7.333893],[8.353857,1.519195,8.549555],[7.884916,0.934549,8.00283],[0.7779847,-2.294533,2.621096],[9.332724,7.178868,11.81676],[5.264858,2.272345,5.820848],[8.682915,-6.8624,11.1124],[-7.559869,2.969754,8.183584],[-1.540396,-5.715488,6.003301],[8.905748,3.171482,9.506349],[-8.752276,-3.099894,9.33872],[5.259305,5.168972,7.441678],[-1.317375,6.53665,6.742645],[-7.07714,1.005969,7.217886],[0.04448572,-4.641075,4.747795],[5.548146,-3.282238,6.52342],[4.574085,6.998615,8.420384],[9.535464,2.508167,9.910398],[-2.542828,3.08561,4.121525],[-5.376565,-3.420582,6.450414],[0.5573341,7.275316,7.364839],[-3.597444,1.242685,3.935209],[-3.034036,-2.06164,3.80207],[-7.089075,-1.116634,7.245817],[-9.931342,8.164815,12.89557],[-1.103872,-8.039215,8.176033],[5.801492,-1.722915,6.133983],[3.585934,-9.767879,10.45325],[-0.2613812,3.058623,3.228544],[9.765979,5.928327,11.46819],[5.067154,-8.258693,9.740743],[9.833155,3.788495,10.58507],[2.705681,5.870229,6.540665],[-3.033623,-4.640895,5.633895],[0.8303958,-2.665977,2.965972],[1.43105,-9.579439,9.737226],[7.571125,3.08919,8.238023],[9.394603,7.092769,11.8138],[-1.108578,-1.271363,1.960946],[-7.270456,9.294601,11.84268],[8.079309,5.896771,10.05222],[0.1225288,-6.99917,7.071308],[4.202919,4.619568,6.324946],[-7.107224,9.033723,11.5378],[-6.447438,9.033813,11.14357],[-3.442085,7.897871,8.673195],[-2.892744,-6.92124,7.567796],[9.609443,-1.787272,9.82526],[2.531256,-8.320899,8.754691],[9.44817,9.898406,13.72029],[-0.4416661,9.36996,9.433516],[-3.212806,7.925436,8.610148],[-9.24003,-9.400044,13.21889],[6.185919,1.330051,6.405828],[0.6978033,-9.095151,9.17653],[1.534233,6.784938,7.027749],[7.258026,-8.493993,11.21726],[2.171502,6.701437,7.115102],[-0.9774226,5.637635,5.808466],[-6.110554,0.1837512,6.194565],[3.456857,-7.246337,8.090691],[-7.192708,-4.507564,8.547116],[8.175223,4.494832,9.382846],[-6.481455,9.852799,11.83583],[5.347162,-9.051738,10.56059],[1.128879,8.839611,8.967335],[6.334618,5.997565,8.780556],[-9.757228,6.237847,11.62387],[-7.047923,6.658036,9.746931],[0.3554588,-7.729477,7.801998],[-1.569239,-4.758078,5.108994],[-0.2862947,-6.050424,6.139185],[9.085858,1.666216,9.291346],[-1.328022,2.174167,2.736905],[-4.711282,7.376837,8.809875],[5.814189,8.94837,10.71812],[-6.467264,0.0189399,6.544147],[-2.98152,6.630057,7.33806],[-9.888857,6.256556,11.74453],[6.452325,6.837381,9.454221],[-4.870245,-8.893534,10.18893],[-3.822709,0.1115255,3.952916],[-4.019165,-8.500865,9.45613],[6.358894,-1.31986,6.570964],[-3.931115,-2.903332,4.988286],[9.09477,7.889347,12.08125],[5.410104,-6.150962,8.252488],[6.449808,6.77903,9.410381],[-9.103885,-4.740682,10.31285],[-1.787893,0.545781,2.120009],[8.988338,-6.252124,10.99451],[6.263295,-3.97896,7.487388],[-9.200784,-8.153522,12.33428],[-1.672017,-1.71004,2.592273],[3.543728,7.130653,8.025225],[-8.614675,4.021636,9.559611],[9.000957,2.153265,9.308801],[-8.109016,-1.139103,8.249467],[-3.146516,0.1120077,3.303499],[2.56628,1.581393,3.17594],[-3.455392,4.47798,5.74387],[-6.648096,8.058856,10.49487],[-9.35427,2.422297,9.714417],[-0.8878292,7.141387,7.265511],[8.781408,0.05091817,8.83831],[-0.6522209,-3.614348,3.806429],[4.352376,2.473576,5.105072],[-7.616106,7.600242,10.80596],[-0.2082577,7.686616,7.754189],[-7.147081,2.936642,7.791318],[-8.581828,1.149056,8.715968],[-5.158416,-6.69232,8.508608],[7.23813,3.773354,8.223669],[-1.14139,-8.561301,8.694748],[-3.914879,4.446796,6.00835],[2.866461,-9.521122,9.993416],[2.607465,9.941417,10.32621],[8.893923,-9.2972,12.90503],[0.3165582,1.885522,2.157638],[2.995021,-5.631243,6.456086],[-6.495008,0.2526748,6.576395],[5.975975,8.62524,10.54073],[-3.961747,-7.484834,8.527495],[8.43512,-1.002244,8.553113],[4.967566,-6.585263,8.309176],[1.45342,-5.02632,5.326943],[1.161562,0.4391521,1.594391],[-0.1787593,-6.628766,6.706153],[-2.047884,-3.052995,3.809804],[3.49123,4.908309,6.105751],[3.438717,-1.447785,3.862752],[-1.829091,-9.598876,9.822627],[-8.947632,6.020453,10.83079],[6.126363,1.458726,6.376536],[-4.202673,2.493553,4.988012],[-4.593033,-2.922624,5.535132],[-8.605408,-9.469215,12.83429],[9.337125,-6.337077,11.32874],[-4.753941,-4.067977,6.336276],[-5.301907,5.276918,7.546926],[-6.424314,-9.432957,11.45655],[-6.760554,-4.480027,8.171641],[7.629502,-3.776235,8.571421],[-6.235689,-1.586765,6.511654],[-1.136132,7.316988,7.471888],[2.742484,9.874319,10.29677],[9.720802,-9.777529,13.82368],[4.059032,5.529339,6.931763],[-5.001196,2.907128,5.87055],[-6.735356,-6.129835,9.161872],[6.354763,0.4406797,6.44804],[9.768438,0.1920953,9.821369],[-4.126097,6.183279,7.500507],[-2.401377,1.689492,3.101772],[4.724698,3.630991,6.042091],[6.708174,-4.079123,7.91447],[-5.457417,-2.731731,6.184315],[2.710279,-9.62024,10.04463],[8.724179,-6.580532,10.97336],[-8.453034,-0.6089199,8.533731],[3.95682,1.271329,4.274658],[5.32605,1.258763,5.563388],[-0.3532307,1.662414,1.9719],[-3.154928,-2.157859,3.95094],[6.873587,-7.360641,10.12053],[1.568869,5.943872,6.22824],[5.849385,9.656276,11.33397],[3.220246,-1.576554,3.722299],[1.001647,-6.673859,6.822293],[-2.500759,-8.906477,9.30479],[-8.969245,1.72965,9.189073],[8.322105,3.98756,9.282137],[-7.114598,-0.009739595,7.184539],[-8.638371,3.267454,9.289657],[-4.973714,4.143791,6.550483],[-3.057006,-7.038534,7.73862],[5.167496,6.676783,8.501908],[7.149391,-2.799904,7.742949],[3.479694,-4.998685,6.172124],[-2.57984,7.905544,8.375751],[-9.345629,-4.396735,10.37651],[-2.919431,-5.876326,6.63734],[-9.787966,-4.307847,10.74066],[-6.643061,-2.349177,7.116804],[-5.91692,-6.800507,9.069555],[-0.719407,9.700235,9.778145],[0.08967771,1.262665,1.613185],[9.308549,-9.968415,13.67547],[-8.082753,-4.721406,9.413957],[-8.245977,-1.416735,8.426344],[4.769582,-3.125557,5.789475],[3.891599,-4.107701,5.746108],[0.7644521,9.936826,10.01623],[-2.136371,-8.575363,8.893871],[4.598984,4.788814,6.714417],[-9.146389,9.111232,12.94878],[1.217703,2.65502,3.087383],[-3.284445,-7.923078,8.634973],[-8.471087,0.7460657,8.562472],[7.285286,-1.61646,7.529166],[8.389912,-6.863554,10.88572],[9.025734,0.5374122,9.09685],[-1.180281,8.293049,8.436096],[-1.376158,-4.843579,5.133621],[-9.705334,4.568925,10.77351],[-7.028446,4.149483,8.222971],[8.344283,2.479819,8.762223],[2.105952,4.307898,4.898267],[1.336579,1.37788,2.164485],[-0.9777235,5.332594,5.512939],[-1.434183,-9.695643,9.852024],[-8.548654,-7.203107,11.22338],[5.407984,6.711101,8.676702],[-5.230673,-2.134977,5.737427],[-1.551727,6.830709,7.075765],[2.156929,1.652314,2.895252],[-5.377623,6.035294,8.145159],[6.155353,-8.575088,10.60285],[1.799313,0.8816051,2.239365],[-0.08841647,-9.264,9.318235],[5.686598,-5.369621,7.88481],[8.336019,-9.838454,12.93385],[5.923166,0.5627087,6.033286],[-3.086261,-8.329372,8.938873],[1.952016,1.502451,2.658519],[7.033173,7.717665,10.48942],[7.365803,9.641955,12.17466],[7.177139,-2.878791,7.797356],[9.015049,5.749077,10.73885],[1.114204,6.516705,6.686471],[2.88548,-7.152742,7.777385],[7.909256,-2.233466,8.279173],[-8.236672,-5.266643,9.827527],[0.05108798,9.262949,9.316912],[8.053933,5.695292,9.914746],[3.461235,-1.086674,3.763112],[-2.362585,1.039229,2.767997],[-6.676546,-7.085462,9.786728],[9.728827,1.869912,9.957242],[6.592977,6.372203,9.223465],[-3.103691,5.101085,6.054252],[-3.525172,6.865457,7.782116],[4.359023,9.0679,10.11078],[-1.630761,-7.875064,8.104074],[5.995005,-5.937586,8.496765],[6.186158,-6.858808,9.290415],[-7.556148,9.927172,12.51575],[8.841129,-8.880692,12.57109],[-9.45821,8.99941,13.09378],[7.053361,8.899969,11.39997],[-6.520841,8.818426,11.01299],[6.645033,-9.071957,11.28968],[9.19484,-8.112863,12.30299],[2.432688,0.8865439,2.775595],[-7.227872,2.353487,7.666879],[8.720655,-8.925817,12.51879],[7.566294,-3.470777,8.384217],[2.886936,-6.156673,6.873065],[6.367738,-7.604116,9.968484],[-4.769611,0.3506742,4.885915],[4.303854,-1.584648,4.694068],[9.965402,4.13878,10.83692],[-8.512636,7.053324,11.1002],[-2.203094,-4.816769,5.390258],[-3.588014,-2.250613,4.351908],[8.801136,2.202357,9.127452],[-2.558908,-6.282612,6.857056],[-9.308473,9.584642,13.39825],[-5.991918,-2.5914,6.604425],[0.5623642,6.222707,6.327586],[-0.6571087,-8.876643,8.956929],[-9.622786,-0.683592,9.698727],[7.648746,0.4170209,7.725104],[9.36125,-4.955989,10.63931],[2.750953,-8.432299,8.925884],[-2.785273,-3.030056,4.235444],[4.209073,3.331489,5.460321],[8.442007,6.878038,10.93503],[-1.007618,-7.805453,7.933498],[-9.179793,-1.923817,9.432374],[9.988093,0.6577474,10.05956],[-0.2984064,-9.596757,9.653331],[6.122307,-2.987323,6.885255],[0.2787552,-4.312842,4.436024],[7.276969,7.039121,10.17367],[-9.319702,-6.869609,11.62103],[-1.397509,0.5471867,1.803454],[-1.59139,1.929206,2.693392],[3.586832,8.910706,9.657434],[-6.80502,-2.046342,7.176058],[4.229972,2.642745,5.086921],[3.822251,-9.495728,10.28487],[7.708927,9.746727,12.46701],[-2.193439,8.969528,9.287821],[-1.089389,-5.339833,5.540811],[-3.71719,7.513505,8.442172],[4.330385,2.46795,5.083602],[8.706561,-1.413793,8.877106],[2.725078,-6.346901,6.979198],[2.886431,-9.144958,9.641666],[-5.701901,-1.516788,5.984339],[9.125038,-0.1299488,9.180589],[-6.35394,2.689199,6.971681],[-2.022666,5.041821,5.523689],[-7.977414,9.369103,12.34582],[7.975768,6.842147,10.55594],[3.228816,-1.648467,3.760678],[6.227774,7.339038,9.67712],[-8.176476,-8.50629,11.8411],[-5.391586,-7.096285,8.96808],[-5.930111,3.837183,7.133736],[9.538124,-1.851256,9.767444],[-6.45814,0.7773647,6.581176],[-5.491472,-5.770138,8.028123],[-2.364752,-7.059816,7.512194],[6.733721,7.232028,9.932031],[9.976801,4.176486,10.86184],[4.379723,3.549239,5.725301],[-7.962668,3.026568,8.576958],[6.291411,-8.900658,10.94548],[-0.6730929,7.737399,7.830734],[8.074352,5.507142,9.824652],[1.724149,-2.291896,3.037347],[8.594369,-3.557268,9.35507],[-6.462662,-6.583423,9.27941],[1.541058,-0.0315937,1.837351],[-0.3924099,-6.33608,6.4265],[-2.9456,-1.561783,3.480765],[3.558359,-2.937077,4.721053],[-5.768657,1.546993,6.055624],[2.444483,1.632054,3.104689],[3.778471,2.283903,4.526926],[2.535997,6.272407,6.839179],[-9.940829,3.06462,10.45045],[8.284396,-9.410103,12.57701],[3.343011,5.572033,6.574441],[1.397865,-1.863913,2.535389],[-3.337565,5.737277,6.712353],[8.09921,-4.957433,9.548473],[-7.24733,1.242639,7.420778],[9.024527,5.593688,10.66449],[7.877152,4.904559,9.332964],[2.986919,1.031407,3.314436],[-9.677384,5.677188,11.2642],[9.230555,-0.4082595,9.293536],[5.113084,5.734038,7.747439],[-7.471643,-7.778997,10.83228],[7.739139,-4.310154,8.91469],[8.200529,-0.03130321,8.261334],[-5.986668,3.295812,6.906704],[9.900335,-3.551423,10.56547],[4.935203,5.698133,7.604271],[6.679972,-6.38053,9.291566],[7.635226,-0.9659994,7.760788],[7.315536,-1.581691,7.551081],[2.87922,4.31873,5.285956],[8.166503,2.875467,8.715508],[-9.76383,9.195138,13.44927],[1.727226,3.843928,4.331177],[-5.012064,4.65922,6.91586],[5.136998,-4.923748,7.185544],[0.2537997,-2.018875,2.267216],[3.602673,1.834312,4.164607],[-4.750686,-7.978419,9.339389],[8.19693,-8.433284,11.80296],[4.231031,-5.941548,7.362311],[-9.704782,-9.70581,13.76174],[8.215248,2.755201,8.722467],[1.861456,4.764777,5.212304],[-7.958154,-8.126219,11.41786],[-3.395,9.26473,9.917725],[8.95293,5.626122,10.62112],[3.61915,1.516023,4.049268],[-1.974789,-9.494608,9.749224],[4.02235,-5.576089,6.947811],[-8.254399,8.762672,12.07972],[9.188443,-6.129087,11.09023],[-2.467096,-7.478632,7.938293],[-6.018831,-6.90888,9.217319],[-0.9246269,-8.153808,8.266771],[-8.681619,-5.875215,10.53037],[5.007256,4.028661,6.504054],[8.868887,-1.117399,8.994761],[-6.903319,9.552498,11.82819],[4.445877,7.22493,8.541981],[7.497507,1.72818,7.758815],[1.589667,8.083004,8.298312],[-3.373686,6.959839,7.798789],[2.497607,-9.614494,9.983814],[0.3086573,-0.3875857,1.116016],[-6.66591,3.472776,7.582514],[7.535618,3.081184,8.202392],[-6.416052,1.917197,6.770626],[9.819991,9.529582,13.72025],[9.186749,-5.752191,10.88504],[0.08668057,4.448658,4.560491],[-4.397523,-5.128588,6.829394],[3.137731,3.527668,4.82595],[7.030819,7.043518,10.00218],[-2.232719,-0.0212764,2.446526],[-9.795912,-2.826422,10.24444],[1.392345,0.03213672,1.714543],[2.7033,2.463945,3.791946],[-5.16085,-3.483878,6.306488],[4.433296,0.5097912,4.573183],[-1.137217,3.851399,4.138422],[5.03279,-2.616686,5.759863],[-8.867289,6.931862,11.29954],[-2.764022,4.76533,5.598946],[2.802314,-6.402376,7.059985],[-2.121371,7.877507,8.219205],[-1.590977,9.286454,9.474673],[-8.816046,-0.6699641,8.897838],[-1.511165,7.416255,7.634427],[2.674253,5.386782,6.096642],[-2.450363,0.8564859,2.781699],[8.622622,9.274506,12.70299],[-8.027534,6.630534,10.4597],[9.735451,-8.171624,12.74968],[-0.7665257,1.042218,1.63517],[-7.717442,-3.626319,8.585401],[4.574441,-2.108612,5.135344],[8.7827,-0.3181025,8.845168],[4.919716,-6.695331,8.368457],[9.190126,1.021974,9.300691],[5.341414,8.733173,10.28586],[4.44009,3.365635,5.660556],[5.662347,-3.5732,6.76978],[1.651683,8.187069,8.411668],[-9.590072,-6.125689,11.42338],[4.541033,-6.895466,8.316756],[-2.195767,5.466592,5.975368],[-9.308893,7.49321,11.99182],[3.285446,-5.056606,6.112563],[-1.483943,-2.245931,2.871636],[-2.51213,-1.271853,2.988044],[1.16326,8.189335,8.331769],[-5.721913,0.4359732,5.824977],[-2.771897,3.870532,4.86461],[9.995577,8.279236,13.01758],[-1.880692,-4.626969,5.093706],[8.364713,6.612001,10.7092],[-5.335392,-8.135152,9.779933],[6.793189,7.190413,9.942307],[9.12758,1.118447,9.250061],[-2.168602,2.384898,3.374993],[1.555179,-2.796114,3.352139],[6.265238,6.2058,8.874973],[5.782384,-2.465317,6.365042],[4.70989,2.887762,5.614466],[-8.6951,6.518927,10.91335],[1.430928,-2.350515,2.927879],[1.996061,0.4968781,2.28717],[-8.098332,8.560811,11.82668],[4.086466,-0.5491541,4.242732],[9.133516,-7.440553,11.82298],[-6.092925,5.946148,8.572071],[-9.833797,2.244525,10.13615],[-0.1375821,-8.464081,8.52406],[9.551702,-5.100273,10.87418],[9.633879,7.942067,12.5255],[-2.577386,5.476626,6.134847],[-3.698343,-7.947155,8.822416],[6.758987,5.999448,9.092706],[6.157216,7.992589,10.13868],[-2.724678,-2.892081,4.097317],[-0.208334,-4.237328,4.35871],[-4.247097,1.053309,4.488574],[-4.107233,-7.94145,8.996443],[3.821431,-1.193496,4.126472],[8.021711,-4.04451,9.039132],[6.607568,-4.007018,7.792056],[-8.146136,0.5319188,8.224504],[3.762984,-4.321376,5.81673],[2.536395,3.067427,4.10395],[-1.669106,6.22882,6.525651],[-0.06719616,-2.348359,2.553293],[4.93654,7.260211,8.836294],[-7.985876,-1.190074,8.135754],[8.974156,-5.601364,10.62595],[3.222809,-8.075059,8.751747],[-4.947002,5.72608,7.632878],[9.239606,9.281075,13.13426],[1.102046,-5.040099,5.255198],[3.736293,-6.597278,7.647481],[0.9232691,5.910995,6.065665],[4.290192,0.8101659,4.479075],[5.860132,1.041424,6.035371],[-9.987302,-7.913805,12.7818],[7.493257,0.1628384,7.561442],[-8.959002,-9.892291,13.38361],[4.858397,-9.849054,11.0276],[-4.517028,-2.604101,5.308944],[0.2754277,8.509908,8.572887],[1.933881,-7.178392,7.50128],[-4.499557,4.019481,6.115737],[8.273193,-1.567188,8.479493],[7.50668,-8.972349,11.74109],[1.729105,-7.248565,7.518743],[-5.081962,-3.962711,6.521458],[-9.400064,1.025065,9.508521],[-6.964202,-6.613914,9.656292],[2.643056,4.718027,5.499593],[1.401744,-9.659814,9.812078],[0.3279737,4.702696,4.819016],[5.917629,-0.03103761,6.001608],[-3.005771,-6.84682,7.544111],[9.658606,7.529014,12.28718],[0.9521365,-4.4302,4.640391],[2.281769,9.334179,9.660919],[8.873232,-9.168705,12.79841],[-5.429903,-8.828409,10.41272],[-7.624604,-3.280442,8.360376],[-4.926692,1.922196,5.382112],[-3.413877,9.663634,10.29759],[1.722756,-3.406496,3.94615],[9.752666,6.841311,11.95483],[2.368335,-4.466183,5.153232],[3.056464,2.812296,4.272117],[-8.36513,1.766829,8.607966],[-3.109254,5.397729,6.308957],[7.545311,4.411655,8.79741],[-5.879078,6.736372,8.996792],[2.104112,7.277143,7.640949],[5.607522,7.028739,9.04696],[7.212952,-4.548128,8.585577],[-3.598277,-8.186029,8.997705],[-8.967277,2.428614,9.343994],[-4.014578,-9.444023,10.3105],[1.005391,8.127069,8.249853],[-4.016331,-6.88186,8.030624],[3.200202,8.972046,9.578043],[-5.649019,8.783291,10.49084],[1.166136,-9.565235,9.687807],[9.387866,-1.83944,9.618502],[2.103669,9.518847,9.799687],[2.327504,-3.712446,4.494389],[3.508416,-4.948055,6.147539],[-7.675206,5.53871,9.517673],[-5.920835,2.441941,6.482234],[2.458244,-3.195405,4.153742],[4.273759,-0.1761405,4.392726],[0.4734253,3.540241,3.709102],[7.449901,7.0696,10.31893],[-3.378277,-4.854597,5.998322],[-2.602253,0.01921602,2.787847],[5.212995,0.7491653,5.36065],[6.749562,1.988581,7.107112],[-9.035769,-8.103188,12.17813],[-6.138606,0.5411665,6.243023],[-5.960494,7.580759,9.695123],[-6.355009,1.782804,6.675667],[4.925193,4.638552,6.839129],[1.001271,6.501187,6.653419],[-3.332393,-9.905554,10.4988],[-8.99926,-6.484553,11.13715],[2.969697,-0.4663791,3.16806],[-4.348559,-8.389857,9.502614],[-0.8416639,-9.91759,10.00335],[4.142294,7.84552,8.928089],[-0.07629377,-3.213423,3.36629],[4.850558,8.469481,9.81122],[-1.200345,-5.942368,6.144312],[4.976585,-7.066247,8.700475],[7.769999,-7.284159,10.69728],[-6.881454,5.5069,8.870195],[0.831877,-2.660263,2.961253],[2.342125,0.6807148,2.636081],[5.510836,-3.155994,6.428811],[-9.445552,-0.5687637,9.515353],[-4.661745,6.335663,7.929218],[-5.721874,-2.923104,6.502644],[-5.456796,0.6304021,5.583371],[9.70949,5.253099,11.08464],[2.766297,0.3725084,2.964989],[-2.697643,7.53075,8.061604],[9.037578,-1.044074,9.152481],[-1.105024,3.466515,3.773302],[4.981349,-0.3324722,5.091598],[-1.258563,7.078665,7.258889],[-1.326885,-4.352457,4.65881],[7.080709,4.938965,8.690789],[-4.172138,-3.929942,5.818176],[5.865673,1.400075,6.1128],[5.779924,-7.636144,9.62903],[-7.033457,-9.869432,12.16039],[7.106021,0.2346319,7.179874],[-1.802717,-2.165224,2.989646],[-4.059247,-2.699739,4.976553],[4.066726,-1.749202,4.538498],[1.512066,-0.00311173,1.81283],[7.164791,-4.746936,8.652609],[-1.814752,-0.8564341,2.242053],[4.325398,-8.249094,9.36785],[-3.541405,-3.075261,4.795704],[0.3424194,7.053362,7.132123],[-9.223133,2.958226,9.737416],[8.395632,1.758363,8.635883],[4.52785,4.677428,6.586331],[-0.2945758,0.7672083,1.294366],[0.5397591,0.836157,1.41085],[-0.0756669,-5.620868,5.70963],[-4.499197,-9.694325,10.73419],[-2.728481,-6.284762,6.924077],[9.595323,4.621696,10.69721],[-2.703941,-8.741934,9.205036],[-5.876688,-1.279203,6.09687],[7.911356,-5.286425,9.567437],[-4.932025,1.046085,5.139958],[7.606833,-3.059488,8.259805],[-3.174063,0.2862914,3.340156],[-1.745291,-2.417857,3.145166],[2.43953,-8.260708,8.671252],[6.550898,-4.358201,7.931468],[-0.1887187,-6.087053,6.171534],[-5.813763,-8.49759,10.34451],[5.140632,4.287502,6.768218],[-0.823933,5.903151,6.043679],[4.662083,-7.369172,8.777227],[-4.819428,-0.2963587,4.930995],[9.863588,-6.758821,11.99883],[0.219812,9.871028,9.923986],[-6.208371,-1.612859,6.491932],[-4.005644,5.337673,6.748032],[3.466263,6.335796,7.290905],[6.925158,7.291933,10.10594],[-9.521317,4.004838,10.37758],[-4.058208,2.781599,5.020592],[4.024561,0.05919294,4.14736],[0.9981067,6.914491,7.057365],[-2.990633,6.932555,7.616049],[-3.742201,-9.852781,10.58685],[5.384836,-9.57519,11.0309],[1.336231,-3.552461,3.924983],[-1.152076,7.141569,7.302691],[-8.123281,-1.050738,8.251772],[1.340369,-9.391962,9.539682],[1.302532,-8.588085,8.743671],[-6.195855,3.711069,7.291135],[-0.4900052,-8.054438,8.131056],[-0.384175,4.239369,4.372624],[4.671275,7.024289,8.494789],[8.028484,8.067039,11.42513],[-9.545441,2.180571,9.842273],[3.586485,-3.863996,5.365942],[-2.693748,5.407343,6.123368],[-9.756125,-5.685719,11.3362],[-8.323656,-6.730685,10.75106],[6.715686,3.911841,7.836002],[-0.5331475,-0.1605588,1.144563],[5.426686,-0.4125124,5.533452],[-3.091128,-3.183656,4.548707],[1.226619,8.797573,8.938786],[-1.643916,-7.290673,7.540317],[-8.73745,3.160568,9.345171],[4.609836,9.810455,10.88557],[-5.145938,-7.250679,8.947235],[3.539937,-9.170094,9.880374],[-1.181191,0.0220716,1.547804],[-7.842568,-2.282121,8.228848],[-9.340626,-6.973419,11.6994],[-7.59757,-8.456885,11.41236],[8.824499,0.8088707,8.917738],[-9.292551,-5.981956,11.09663],[2.679304,-6.723249,7.306213],[-9.773339,-3.998205,10.60678],[-1.231197,4.324911,4.606593],[2.647211,-9.205343,9.630476],[-3.647274,-4.101408,5.578903],[-7.99578,-2.528956,8.445597],[3.346361,2.129011,4.090333],[-7.206783,0.3089704,7.282388],[5.662699,-5.16792,7.731336],[-6.622648,-2.856936,7.28159],[8.087882,-2.338845,8.478445],[7.950572,-8.582885,11.74213],[-5.678049,0.04600915,5.765618],[-6.245456,-4.856081,7.974161],[8.224795,-0.877148,8.331666],[-2.677323,-7.949429,8.447573],[-8.083122,4.858543,9.483792],[-4.492181,3.573102,5.826384],[9.660138,-2.292943,9.97877],[-8.395582,4.852795,9.748611],[3.511057,-1.045152,3.79735],[-5.656837,9.748077,11.31481],[-4.425598,2.004305,4.960156],[-7.981494,-3.708503,8.85761],[8.10352,7.087535,10.81204],[2.535029,-8.444074,8.872923],[1.205052,-3.708745,4.025784],[-4.685283,7.180394,8.631913],[-8.130489,3.515763,8.914339],[5.918495,-6.840697,9.100754],[0.5885875,1.728997,2.082274],[7.349698,-7.711635,10.69988],[3.457119,1.87545,4.058199],[0.7777522,-0.7449881,1.469662],[-3.812182,-3.85495,5.513019],[-6.263912,-0.8475741,6.399607],[9.803756,-3.740854,10.54076],[7.468621,-1.014493,7.603256],[-1.970203,-5.905822,6.305587],[4.240612,6.131842,7.522119],[6.221618,8.014585,10.1952],[-1.50288,-8.030087,8.230489],[-0.4661069,-4.893843,5.016667],[3.192487,-0.9891909,3.48862],[-7.773095,-2.618037,8.262877],[-1.550666,-3.357833,3.831398],[1.418966,-2.742363,3.245616],[-9.234354,2.845675,9.714482],[5.805252,4.190757,7.229342],[-6.437168,0.8579642,6.570634],[-7.715916,-5.376927,9.457627],[-9.270813,-7.393644,11.90017],[1.991181,-4.085368,4.653497],[-7.930744,-2.391779,8.343699],[3.486256,-8.814675,9.531657],[-4.281052,-9.027506,10.04108],[7.867065,5.990905,9.938896],[0.1555087,-1.632546,1.920779],[-1.381177,2.201705,2.784808],[4.546048,0.04054535,4.654911],[7.846531,5.429543,9.594164],[-2.628924,-5.448638,6.131795],[1.400033,3.004203,3.461983],[-8.32434,-5.573797,10.06786],[-7.490795,3.721234,8.423752],[-6.838972,7.89924,10.49617],[-9.674415,-6.652012,11.78319],[7.311444,-2.915282,7.934487],[9.077174,-7.231893,11.64884],[7.693282,2.890139,8.278858],[-6.52078,-1.5308,6.772291],[-5.212349,7.728436,9.375357],[3.317288,-2.135614,4.070042],[-1.628779,9.846558,10.03033],[-0.05581499,-9.821731,9.872664],[3.464185,-6.748479,7.65131],[1.035163,-5.880105,6.053692],[-1.79543,-6.76516,7.070428],[1.410998,4.414864,4.741512],[-4.471932,9.654901,10.68715],[3.524183,8.008802,8.806859]],"ignoreExtent":false,"flags":4096},"9":{"id":9,"type":"text","material":{"lit":false},"vertices":[[0.000690937,-13.38884,-1.203414]],"colors":[[0,0,0,1]],"texts":[["x"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[0.000690937,-13.38884,-1.203414]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"10":{"id":10,"type":"text","material":{"lit":false},"vertices":[[-13.38799,-0.001586914,-1.203414]],"colors":[[0,0,0,1]],"texts":[["y"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-13.38799,-0.001586914,-1.203414]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"11":{"id":11,"type":"text","material":{"lit":false},"vertices":[[-13.38799,-13.38884,7.568857]],"colors":[[0,0,0,1]],"texts":[["z"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-13.38799,-13.38884,7.568857]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"5":{"id":5,"type":"light","vertices":[[0,0,1]],"colors":[[1,1,1,1],[1,1,1,1],[1,1,1,1]],"viewpoint":true,"finite":false},"4":{"id":4,"type":"background","material":{"fog":true},"colors":[[0.2980392,0.2980392,0.2980392,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"6":{"id":6,"type":"background","material":{"lit":false,"back":"lines"},"colors":[[1,1,1,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"8":{"id":8,"type":"bboxdeco","material":{"front":"lines","back":"lines"},"vertices":[[-5,"NA","NA"],[0,"NA","NA"],[5,"NA","NA"],["NA",-5,"NA"],["NA",0,"NA"],["NA",5,"NA"],["NA","NA",2],["NA","NA",4],["NA","NA",6],["NA","NA",8],["NA","NA",10],["NA","NA",12],["NA","NA",14]],"colors":[[0,0,0,1]],"draw_front":true,"newIds":[19,20,21,22,23,24,25]},"1":{"id":1,"type":"subscene","par3d":{"antialias":8,"FOV":30,"ignoreExtent":false,"listeners":1,"mouseMode":{"left":"trackball","right":"zoom","middle":"fov","wheel":"pull"},"observer":[0,0,76.26823],"modelMatrix":[[0.8998286,0,0,-0.0006217249],[0,0.3077923,1.29054,-9.767421],[0,-0.8456525,0.469718,-79.82481],[0,0,0,1]],"projMatrix":[[2.665751,0,0,0],[0,3.732051,0,0],[0,0,-3.863704,-274.9382],[0,0,-1,0]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,0.3420201,0.9396926,0],[0,-0.9396926,0.3420201,0],[0,0,0,1]],"scale":[0.8998286,0.8999246,1.373363],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":1,"bbox":[-9.998322,9.999704,-9.999534,9.99636,1.017497,14.12022],"windowRect":[100,100,772,580],"family":"sans","font":1,"cex":1,"useFreeType":false,"fontname":"TT Arial","maxClipPlanes":8,"glVersion":4.5},"embeddings":{"viewport":"replace","projection":"replace","model":"replace"},"objects":[6,8,7,9,10,11,5,19,20,21,22,23,24,25],"subscenes":[],"flags":6736},"19":{"id":19,"type":"lines","material":{"lit":false},"vertices":[[-5,-10.29947,0.8209564],[5,-10.29947,0.8209564],[-5,-10.29947,0.8209564],[-5,-10.81437,0.4835614],[0,-10.29947,0.8209564],[0,-10.81437,0.4835614],[5,-10.29947,0.8209564],[5,-10.81437,0.4835614]],"colors":[[0,0,0,1]],"centers":[[0,-10.29947,0.8209564],[-5,-10.55692,0.6522589],[0,-10.55692,0.6522589],[5,-10.55692,0.6522589]],"ignoreExtent":true,"origId":8,"flags":64},"20":{"id":20,"type":"text","material":{"lit":false},"vertices":[[-5,-11.84415,-0.1912287],[0,-11.84415,-0.1912287],[5,-11.84415,-0.1912287]],"colors":[[0,0,0,1]],"texts":[["-5"],["0"],["5"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-5,-11.84415,-0.1912287],[0,-11.84415,-0.1912287],[5,-11.84415,-0.1912287]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"21":{"id":21,"type":"lines","material":{"lit":false},"vertices":[[-10.29829,-5,0.8209564],[-10.29829,5,0.8209564],[-10.29829,-5,0.8209564],[-10.81324,-5,0.4835614],[-10.29829,0,0.8209564],[-10.81324,0,0.4835614],[-10.29829,5,0.8209564],[-10.81324,5,0.4835614]],"colors":[[0,0,0,1]],"centers":[[-10.29829,0,0.8209564],[-10.55577,-5,0.6522589],[-10.55577,0,0.6522589],[-10.55577,5,0.6522589]],"ignoreExtent":true,"origId":8,"flags":64},"22":{"id":22,"type":"text","material":{"lit":false},"vertices":[[-11.84314,-5,-0.1912287],[-11.84314,0,-0.1912287],[-11.84314,5,-0.1912287]],"colors":[[0,0,0,1]],"texts":[["-5"],["0"],["5"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-11.84314,-5,-0.1912287],[-11.84314,0,-0.1912287],[-11.84314,5,-0.1912287]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"23":{"id":23,"type":"lines","material":{"lit":false},"vertices":[[-10.29829,-10.29947,2],[-10.29829,-10.29947,14],[-10.29829,-10.29947,2],[-10.81324,-10.81437,2],[-10.29829,-10.29947,4],[-10.81324,-10.81437,4],[-10.29829,-10.29947,6],[-10.81324,-10.81437,6],[-10.29829,-10.29947,8],[-10.81324,-10.81437,8],[-10.29829,-10.29947,10],[-10.81324,-10.81437,10],[-10.29829,-10.29947,12],[-10.81324,-10.81437,12],[-10.29829,-10.29947,14],[-10.81324,-10.81437,14]],"colors":[[0,0,0,1]],"centers":[[-10.29829,-10.29947,8],[-10.55577,-10.55692,2],[-10.55577,-10.55692,4],[-10.55577,-10.55692,6],[-10.55577,-10.55692,8],[-10.55577,-10.55692,10],[-10.55577,-10.55692,12],[-10.55577,-10.55692,14]],"ignoreExtent":true,"origId":8,"flags":64},"24":{"id":24,"type":"text","material":{"lit":false},"vertices":[[-11.84314,-11.84415,2],[-11.84314,-11.84415,4],[-11.84314,-11.84415,6],[-11.84314,-11.84415,8],[-11.84314,-11.84415,10],[-11.84314,-11.84415,12],[-11.84314,-11.84415,14]],"colors":[[0,0,0,1]],"texts":[["2"],["4"],["6"],["8"],["10"],["12"],["14"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-11.84314,-11.84415,2],[-11.84314,-11.84415,4],[-11.84314,-11.84415,6],[-11.84314,-11.84415,8],[-11.84314,-11.84415,10],[-11.84314,-11.84415,12],[-11.84314,-11.84415,14]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":8,"flags":2064},"25":{"id":25,"type":"lines","material":{"lit":false},"vertices":[[-10.29829,-10.29947,0.8209564],[-10.29829,10.2963,0.8209564],[-10.29829,-10.29947,14.31676],[-10.29829,10.2963,14.31676],[-10.29829,-10.29947,0.8209564],[-10.29829,-10.29947,14.31676],[-10.29829,10.2963,0.8209564],[-10.29829,10.2963,14.31676],[-10.29829,-10.29947,0.8209564],[10.29967,-10.29947,0.8209564],[-10.29829,-10.29947,14.31676],[10.29967,-10.29947,14.31676],[-10.29829,10.2963,0.8209564],[10.29967,10.2963,0.8209564],[-10.29829,10.2963,14.31676],[10.29967,10.2963,14.31676],[10.29967,-10.29947,0.8209564],[10.29967,10.2963,0.8209564],[10.29967,-10.29947,14.31676],[10.29967,10.2963,14.31676],[10.29967,-10.29947,0.8209564],[10.29967,-10.29947,14.31676],[10.29967,10.2963,0.8209564],[10.29967,10.2963,14.31676]],"colors":[[0,0,0,1]],"centers":[[-10.29829,-0.001586914,0.8209564],[-10.29829,-0.001586914,14.31676],[-10.29829,-10.29947,7.568857],[-10.29829,10.2963,7.568857],[0.000690937,-10.29947,0.8209564],[0.000690937,-10.29947,14.31676],[0.000690937,10.2963,0.8209564],[0.000690937,10.2963,14.31676],[10.29967,-0.001586914,0.8209564],[10.29967,-0.001586914,14.31676],[10.29967,-10.29947,7.568857],[10.29967,10.2963,7.568857]],"ignoreExtent":true,"origId":8,"flags":64}},"snapshot":"data:image/png;base64,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oCEM3AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDzAzek6/988vgMA94P/5wD3xyo66g7wJfBfHeArEF1H3QG+B/63A3wFCDzAt8H/doBvAXUH+Cr4Dw/wFWDBA3wb/G8H+AoQeIBvg//tAPeHLHqAL4T/6gA3h3XwAN8J/88BAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMANQeABAABuCAIPAABwQxB4AACAG4LAAwAA3BAEHgAA4IYg8AAAADcEgQcAALghCDwAAMAN+T8S5XoAuyK99gAAAABJRU5ErkJggg==","width":673,"height":481,"sphereVerts":{"vb":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1],[0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0]],"it":[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270],[17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288],[18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271]],"primitivetype":"triangle","material":null,"normals":null,"texcoords":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1]]},"context":{"shiny":false,"rmarkdown":"html_document"},"crosstalk":{"key":[],"group":[],"id":[],"options":[]}});
unnamed_chunk_3rgl.prefix = "unnamed_chunk_3";
</script>
<p id="unnamed_chunk_3debug">
You must enable Javascript to view this page properly.</p>
<script>unnamed_chunk_3rgl.start();</script>

```r
pd.crd <- apply(xyz,1,my.pdisk.coords)
plot(t(pd.crd),pch=20,cex=0.1,asp=TRUE)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-4-1.png" width="672" />

双曲面の外側にいくほど、円板上で点が密になる。

円板上では、短い距離にみえても、ミンコフスキー空間上では遠い距離に相当するからである。

### ミンコフスキー内積と点間の双曲幾何的距離


ベクトルの長さのようなもの(二次形式)が以下のように与えられているとき、
$$
Q(x_1,x_2,...) = -(x_1^2+...+x_{n}^2) + x_{n+1}^2
$$

双曲面は$Q(x_1,...) = 1$上の点の集まりのことである。

この双曲面上の２点に相当するベクトル$u,v$があったとき、その内積は

$$
<u,v> = \frac{Q(u+v) - (Q(u)+Q(v))}{2} = -(u_1v_2 + ... + u_n v_n) + u_{n+1}v_{n+1}
$$
で与えられる。

この２点の双曲幾何的距離は
$$
d(u,v) = arcosh(<u,v>)
$$

### 距離を視覚化

単位円板上に点を取り、そのうちのある１点からの双曲幾何的距離を視覚化してみる。


```r
n <- 10000
x. <- runif(n,min=-1,max=1)
y. <- runif(n,min=-1,max=1)
r2 <- x.^2 + y.^2
s <- which(r2 < 1)
x. <- x.[s]
y. <- y.[s]
plot(x.,y.,pch=20,cex=0.1,asp=TRUE)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-5-1.png" width="672" />


```r
hyp.cds <- t(apply(cbind(x.,y.),1,my.pdisk.coords.inv))
plot3d(hyp.cds)
```

<div id="unnamed_chunk_6div" class="rglWebGL"></div>
<script type="text/javascript">
var unnamed_chunk_6div = document.getElementById("unnamed_chunk_6div"),
unnamed_chunk_6rgl = new rglwidgetClass();
unnamed_chunk_6div.width = 673;
unnamed_chunk_6div.height = 481;
unnamed_chunk_6rgl.initialize(unnamed_chunk_6div,
unnamed_chunk_6rgl.prefix = "unnamed_chunk_6";
</script>
<p id="unnamed_chunk_6debug">
You must enable Javascript to view this page properly.</p>
<script>unnamed_chunk_6rgl.start();</script>


```r
my.hyp.dist <- function(u,v){
  acosh(- sum(u*v) + 2*u[length(u)]*v[length(v)])
}
d <- apply(hyp.cds,1,my.hyp.dist,hyp.cds[1,])
col <- ceiling(d/max(d)*10)
plot(x.,y.,pch=20,cex=1,asp=TRUE,col=col)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-7-1.png" width="672" />

### 直線・測地戦

双曲幾何的直線は、双曲面と３次元平面との交わり。

$$
<u,u> = -1\\
<v,v> = 1\\
<u,v> = 0
$$

という２つのベクトルが、パラメタ$w$を用いて
$$
u \cosh{w} + v \sinh{w}
$$
で表される点がそうなるのだと言う。

$u = (u_1,u_2,u_3),v = (v_1,v_2,v_3)$のうち、$u_1,u_2,v_1$を指定し残り３変数を計算する関数を以下のように作る。

上記の制約を満足する$v_2,v_3$の取り方は４通りあるが、そのうち、双曲面上の曲線となるのは、以下に示すように２通り。

点$u$は、双曲面上の点、点$v$はベクトル$uv$が双曲面の接面となるような点


```r
my.hyp.line.vec <- function(u1,u2,v1){
  u3 <- sqrt(u1^2+u2^2+1)
  tmp <- sqrt(u1^2*v1^2*u2^2 + (u1^2+v1^2+1)*(u1^2+u2^2+1))
  v2a <- (u1*v1*u2 + tmp )/(u1^2+v1^2+1)
  v2b <- (u1*v1*u2 - tmp )/(u1^2+v1^2+1)
  u <- c(u1,u2,u3)
  v3a <- sqrt(v1^2+v2a^2-1)
  v3b <- sqrt(v1^2+v2b^2-1)
  va <- c(v1,v2a,v3a)
  va. <- c(v1,v2a,-v3a)
  vb <- c(v1,v2b,v3b)
  vb. <- c(v1,v2b,-v3b)
  return(list(u=u,va=va,va.=va.,vb=vb,vb.=vb.))
}
```


```r
u1 <- 0.5
u2 <- 0.7
v1 <- -0.3
uv <- my.hyp.line.vec(u1,u2,v1)
uv
```

```
## $u
## [1] 0.500000 0.700000 1.319091
## 
## $va
## [1] -0.3000000  1.0638534  0.4709397
## 
## $va.
## [1] -0.3000000  1.0638534 -0.4709397
## 
## $vb
## [1] -0.30000 -1.22057  0.76144
## 
## $vb.
## [1] -0.30000 -1.22057 -0.76144
```


```r
w <- seq(from=-2,to=2,length=1000)
coshw <- cosh(w)
sinhw <- sinh(w)
cshw <- rbind(coshw,sinhw)
matplot(t(cshw),type="l")
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-10-1.png" width="672" />

```r
u <- uv$u
va <- uv$va 
va.  <- uv$va.
vb <- uv$vb 
vb. <- uv$vb. 

lina <- cbind(u,va) %*% cshw
lina. <- cbind(u,va.) %*% cshw

linb<- cbind(u,vb) %*% cshw
linb. <- cbind(u,vb.) %*% cshw
```

<div id="unnamed_chunk_10div" class="rglWebGL"></div>
<script type="text/javascript">
var unnamed_chunk_10div = document.getElementById("unnamed_chunk_10div"),
unnamed_chunk_10rgl = new rglwidgetClass();
unnamed_chunk_10div.width = 673;
unnamed_chunk_10div.height = 481;
unnamed_chunk_10rgl.initialize(unnamed_chunk_10div,
{"material":{"color":"#000000","alpha":1,"lit":true,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":false,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less","isTransparent":false},"rootSubscene":1,"objects":{"26":{"id":26,"reuse":"unnamed_chunk_6div"},"28":{"id":28,"reuse":"unnamed_chunk_6div"},"29":{"id":29,"reuse":"unnamed_chunk_6div"},"30":{"id":30,"reuse":"unnamed_chunk_6div"},"5":{"id":5,"reuse":"unnamed_chunk_3div"},"6":{"id":6,"reuse":"unnamed_chunk_3div"},"27":{"id":27,"type":"bboxdeco","material":{"front":"lines","back":"lines"},"vertices":[[0,"NA","NA"],[5000,"NA","NA"],[10000,"NA","NA"],["NA",0,"NA"],["NA",5000,"NA"],["NA",10000,"NA"],["NA","NA",5000],["NA","NA",10000],["NA","NA",15000]],"colors":[[0,0,0,1]],"draw_front":true,"newIds":[52,53,54,55,56,57,58]},"1":{"id":1,"type":"subscene","par3d":{"antialias":8,"FOV":30,"ignoreExtent":false,"listeners":1,"mouseMode":{"left":"trackball","right":"zoom","middle":"fov","wheel":"pull"},"observer":[0,0,69131.88],"modelMatrix":[[1.013613,0,0,-3601.799],[0,0.3467469,0.9154097,-9325.506],[0,-0.9526793,0.3331819,-67358.25],[0,0,0,1]],"projMatrix":[[2.665751,0,0,0],[0,3.732051,0,0],[0,0,-3.863703,-249212.4],[0,0,-1,0]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,0.3420201,0.9396926,0],[0,-0.9396926,0.3420201,0],[0,0,0,1]],"scale":[1.013613,1.01382,0.9741587],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":1,"bbox":[-4492.564,11599.42,-3254.372,12834.32,1.000455,16744.72],"windowRect":[100,100,772,580],"family":"sans","font":1,"cex":1,"useFreeType":false,"fontname":"TT Arial","maxClipPlanes":8,"glVersion":4.5},"embeddings":{"viewport":"replace","projection":"replace","model":"replace"},"objects":[6,27,26,28,29,30,5,52,53,54,55,56,57,58],"subscenes":[],"flags":6736},"52":{"id":52,"type":"lines","material":{"lit":false},"vertices":[[0,-3495.703,-250.1553],[10000,-3495.703,-250.1553],[0,-3495.703,-250.1553],[0,-3909.986,-681.3061],[5000,-3495.703,-250.1553],[5000,-3909.986,-681.3061],[10000,-3495.703,-250.1553],[10000,-3909.986,-681.3061]],"colors":[[0,0,0,1]],"centers":[[5000,-3495.703,-250.1553],[0,-3702.844,-465.7307],[5000,-3702.844,-465.7307],[10000,-3702.844,-465.7307]],"ignoreExtent":true,"origId":27,"flags":64},"53":{"id":53,"type":"text","material":{"lit":false},"vertices":[[0,-4738.554,-1543.608],[5000,-4738.554,-1543.608],[10000,-4738.554,-1543.608]],"colors":[[0,0,0,1]],"texts":[["0"],["5000"],["10000"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[0,-4738.554,-1543.608],[5000,-4738.554,-1543.608],[10000,-4738.554,-1543.608]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":27,"flags":2064},"54":{"id":54,"type":"lines","material":{"lit":false},"vertices":[[-4733.944,0,-250.1553],[-4733.944,10000,-250.1553],[-4733.944,0,-250.1553],[-5148.313,0,-681.3061],[-4733.944,5000,-250.1553],[-5148.313,5000,-681.3061],[-4733.944,10000,-250.1553],[-5148.313,10000,-681.3061]],"colors":[[0,0,0,1]],"centers":[[-4733.944,5000,-250.1553],[-4941.129,0,-465.7307],[-4941.129,5000,-465.7307],[-4941.129,10000,-465.7307]],"ignoreExtent":true,"origId":27,"flags":64},"55":{"id":55,"type":"text","material":{"lit":false},"vertices":[[-5977.05,0,-1543.608],[-5977.05,5000,-1543.608],[-5977.05,10000,-1543.608]],"colors":[[0,0,0,1]],"texts":[["0"],["5000"],["10000"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-5977.05,0,-1543.608],[-5977.05,5000,-1543.608],[-5977.05,10000,-1543.608]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":27,"flags":2064},"56":{"id":56,"type":"lines","material":{"lit":false},"vertices":[[-4733.944,-3495.703,5000],[-4733.944,-3495.703,15000],[-4733.944,-3495.703,5000],[-5148.313,-3909.986,5000],[-4733.944,-3495.703,10000],[-5148.313,-3909.986,10000],[-4733.944,-3495.703,15000],[-5148.313,-3909.986,15000]],"colors":[[0,0,0,1]],"centers":[[-4733.944,-3495.703,10000],[-4941.129,-3702.844,5000],[-4941.129,-3702.844,10000],[-4941.129,-3702.844,15000]],"ignoreExtent":true,"origId":27,"flags":64},"57":{"id":57,"type":"text","material":{"lit":false},"vertices":[[-5977.05,-4738.554,5000],[-5977.05,-4738.554,10000],[-5977.05,-4738.554,15000]],"colors":[[0,0,0,1]],"texts":[["5000"],["10000"],["15000"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-5977.05,-4738.554,5000],[-5977.05,-4738.554,10000],[-5977.05,-4738.554,15000]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"origId":27,"flags":2064},"58":{"id":58,"type":"lines","material":{"lit":false},"vertices":[[-4733.944,-3495.703,-250.1553],[-4733.944,13075.65,-250.1553],[-4733.944,-3495.703,16995.88],[-4733.944,13075.65,16995.88],[-4733.944,-3495.703,-250.1553],[-4733.944,-3495.703,16995.88],[-4733.944,13075.65,-250.1553],[-4733.944,13075.65,16995.88],[-4733.944,-3495.703,-250.1553],[11840.8,-3495.703,-250.1553],[-4733.944,-3495.703,16995.88],[11840.8,-3495.703,16995.88],[-4733.944,13075.65,-250.1553],[11840.8,13075.65,-250.1553],[-4733.944,13075.65,16995.88],[11840.8,13075.65,16995.88],[11840.8,-3495.703,-250.1553],[11840.8,13075.65,-250.1553],[11840.8,-3495.703,16995.88],[11840.8,13075.65,16995.88],[11840.8,-3495.703,-250.1553],[11840.8,-3495.703,16995.88],[11840.8,13075.65,-250.1553],[11840.8,13075.65,16995.88]],"colors":[[0,0,0,1]],"centers":[[-4733.944,4789.973,-250.1553],[-4733.944,4789.973,16995.88],[-4733.944,-3495.703,8372.859],[-4733.944,13075.65,8372.859],[3553.427,-3495.703,-250.1553],[3553.427,-3495.703,16995.88],[3553.427,13075.65,-250.1553],[3553.427,13075.65,16995.88],[11840.8,4789.973,-250.1553],[11840.8,4789.973,16995.88],[11840.8,-3495.703,8372.859],[11840.8,13075.65,8372.859]],"ignoreExtent":true,"origId":27,"flags":64}},"snapshot":"data:image/png;base64,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","width":673,"height":481,"sphereVerts":{"reuse":"unnamed_chunk_3div"},"context":{"shiny":false,"rmarkdown":"html_document"},"crosstalk":{"key":[],"group":[],"id":[],"options":[]}});
unnamed_chunk_10rgl.prefix = "unnamed_chunk_10";
</script>
<p id="unnamed_chunk_10debug">
You must enable Javascript to view this page properly.</p>
<script>unnamed_chunk_10rgl.start();</script>




```r
n <- 10000
x <- runif(n,min=-2,max=2)
y <- runif(n,min=-2,max=2)
z <- sqrt(x^2 + y^2 + 1)

xyz <- cbind(x,y,z*sample(c(-1,1),n,replace=TRUE))

plot3d(xyz)
spheres3d(t(lina ),radius=0.05,col=2)
spheres3d(t(lina. ),radius=0.05,col=3) # 双曲面に乗らない
spheres3d(t(linb ),radius=0.05,col=4) # 双曲面に乗らない
spheres3d(t(linb. ),radius=0.05,col=5) #
spheres3d(u,radius=0.1,col=6)
spheres3d(va,radius=0.1,col=2)
spheres3d(va.,radius=0.1,col=3)
spheres3d(vb,radius=0.1,col=4)
spheres3d(vb.,radius=0.1,col=5)
segments3d(rbind(rep(0,3),u))

segments3d(rbind(rep(0,3),va))
segments3d(rbind(rep(0,3),va. ))

segments3d(rbind(rep(0,3),vb ))
segments3d(rbind(rep(0,3),vb. ))
```

<div id="unnamed_chunk_11div" class="rglWebGL"></div>
<script type="text/javascript">
var unnamed_chunk_11div = document.getElementById("unnamed_chunk_11div"),
unnamed_chunk_11rgl = new rglwidgetClass();
unnamed_chunk_11div.width = 673;
unnamed_chunk_11div.height = 481;
unnamed_chunk_11rgl.initialize(unnamed_chunk_11div,
unnamed_chunk_11rgl.prefix = "unnamed_chunk_11";
</script>
<p id="unnamed_chunk_11debug">
You must enable Javascript to view this page properly.</p>
<script>unnamed_chunk_11rgl.start();</script>

この直線を単位円板に射影する。


```r
n <- 10000
x. <- runif(n,min=-1,max=1)
y. <- runif(n,min=-1,max=1)
r2 <- x.^2 + y.^2
s <- which(r2 < 1)
x. <- x.[s]
y. <- y.[s]
plot(x.,y.,pch=20,cex=0.01,asp=TRUE)

pd.crd.a  <- apply(t(lina ),1,my.pdisk.coords)
points(t(pd.crd.a),pch=20,cex=1,col=2)
pd.crd.b.  <- apply(t(linb. ),1,my.pdisk.coords)
points(t(pd.crd.b.),pch=20,cex=1,col=5)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-12-1.png" width="672" />

## ポアンカレ・ディスク

この節は[Poincare Disk model (Wikipedia記事)](https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model)
に従てって書いている。

### 直線

単位円板の辺縁の点と直交するすべてのユークリッド幾何の意味での円周は、双曲幾何での直線に相当する。


```r
t <- seq(from=0,to=1,length=100)*2*pi
plot(cos(t),sin(t),type="l",asp=1)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-13-1.png" width="672" />

単位円板のある点$(\cos{\theta},\sin{\theta})$を通り、単位円周と直交する円の中心は、
その点を通る、単位円の接線上の点である。

接線の方向ベクトルは$(-\sin{\theta},\cos{\theta})$であるから、

そのような円の座標は
$$
ctr_x = \cos{\theta} - r \sin{\theta}\\
ctr_y = \sin{\theta} + r \cos{\theta}
$$
であり、その半径は、$r$である。


```r
# 円周上の点座標ptと直交円の半径rを与え
# 中心を返す関数
my.orth.circle1 <- function(pt,r){
  ctr_x <- pt[1] - r * pt[2]
  ctr_y <- pt[2] + r * pt[1]
  return(c(ctr_x,ctr_y))
}
```

```r
plot(cos(t),sin(t),type="l",asp=1)

n.line <- 10

for(i in 1:n.line){
  theta <- runif(1) * 2 * pi
  r <- 10^runif(1)
  ctr <- my.orth.circle1(c(cos(theta),sin(theta)),r)
  
  pnts <- cbind(r*cos(t)+ctr[1],r*sin(t)+ctr[2])
  points(pnts,type="l",col=i+1)
}
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-15-1.png" width="672" />

#### Inversion 反転を使う方法

平面上の任意の２点を通る円で、単位円と直交するものを求める。

ある点がその円周上にあるとき、その点を単位円に関して反転した点も通ることを利用する。


```r
my.inversion <- function(x){
  r <- sqrt(sum(x^2))
  R <- 1/r
  return(x/r*R)
}
# P, Qは任意の２点の座標
my.orth.circle2 <- function(P,Q){
  P. <- my.inversion(P)
  Q. <- my.inversion(Q)
  M <- (P + P.)/2
  N <- (Q + Q.)/2
  
  k. <- (M[1]^2+M[2]^2-M[1]*N[1]-M[2]*N[2])/(M[1]*N[2]-M[2]*N[1])
  ctr_x <- N[1] + k.*N[2]
  ctr_y <- N[2] - k.*N[1]
  r <- sqrt((ctr_x-P[1])^2+(ctr_y-P[2])^2)
  
  return(list(ctr=c(ctr_x,ctr_y),r=r,P=P,Q=Q,P.=P.,Q.=Q.))
}
```



```r
plot(cos(t),sin(t),type="l",asp=1,xlim=c(-5,5),ylim=c(-5,5))

n.line <- 4

for(i in 1:n.line){
  P <- rnorm(2) * 1
  Q <- rnorm(2) *1
  #P <- c(2,0)
  #Q <- c(0,2)
  tmp <- my.orth.circle2(P,Q)
  r <- tmp$r
  ctr <- tmp$ctr
  pnts <- cbind(r*cos(t)+ctr[1],r*sin(t)+ctr[2])
  points(pnts,type="l",col=i+1)
  points(P[1],P[2],col=i+1,pch=20,cex=2)
  points(Q[1],Q[2],col=i+1,pch=20,cex=2)
  points(tmp$P.[1],tmp$P.[2],col=i+1,pch=20,cex=2)
  points(tmp$Q.[1],tmp$Q.[2],col=i+1,pch=20,cex=2)
}
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-17-1.png" width="672" />

### 距離

上記で描いた、ポアンカレ・ディスク上の曲線が測地線。
測地戦をたどって「距離」を測る。

円板内の２点を通る測地線を延ばして円周との交点を取る。

順番に、A,P,Q,B(A,Bは円周上の点、P,Qは円板内の点)とする。

双曲幾何的なPQ距離$d_{hyp}(P,Q)$は、ユークリッド的距離$|XY|$を使って、以下のように定まる。

$$
d_{hyp}(P,Q) = \log{\frac{|AQ||PB|}{|AP||QB|} }
$$

Pが円周上に近づくと、AP間のユークリッド距離が0に近づくので、上の式の分母が0に近づき、分数自体は無限大に向かう。

逆に、PとQとが近づくと、分数が1に近づくので、双曲幾何的距離も0に近づく。

別法もある。
ディスクの半径をrとし(今はr=1と固定して考えている）、ポアンカレ・ディスクの中心点をOとしたときに以下で得られる。
$$
d_{hyp}(P,Q) = arcosh(1 + \frac{2 |PQ|^2r^2}{(r^2-|OP|^2)(r^2-|OQ|^2)})\\
=arcosh(1+\frac{|PQ|^2}{(1-|OP|^2)(1-|OQ|^2)})
$$

２つの方法で計算して結果を比較してみる。

```r
my.hyp.dist1 <- function(P,Q){
  out <- my.orth.circle2(P,Q)
  r <- out$r
  d.ctr <- sqrt(sum(out$ctr^2)) # 原点から測地線円の中心までの距離
  theta <- acos(1/d.ctr)
  phi <- Arg(out$ctr[1] + 1i * out$ctr[2]) # 測地線円方向の角度
  AB1 <- c(cos(phi-theta),sin(phi-theta))
  AB2 <- c(cos(phi+theta),sin(phi+theta))

  AQ <- sqrt(sum((AB1-Q)^2))
  PB <- sqrt(sum((P-AB2)^2))
  AP <- sqrt(sum((AB1-P)^2))
  QB <- sqrt(sum((Q-AB2)^2))
  # A,BをAB1,AB2のどちらにするか、２通りで距離を算出する
  dhyp1 <- log(AQ*PB/(AP*QB))
  dhyp2 <- log(QB*AP/(PB*AQ))
  
  return(list(dhyp1=dhyp1,dhyp2=dhyp2,A=AB1,B=AB2,circle=out))
}
my.hyp.dist2 <- function(P,Q,r=1){
  #r <- 1
  PQ <- sqrt(sum((P-Q)^2))
  OP <- sqrt(sum((P)^2))
  OQ <- sqrt(sum((Q)^2))
  
  dhyp <- acosh(1+2*PQ^2*r^2/((r^2-OP^2)*(r^2-OQ^2)))
  return(dhyp)
}
```


```r
r4 <- runif(4)
#r4[3] <- r4[1]+rnorm(1)*0.0001
#r4[4] <- r4[2] + rnorm(1)*0.0001
P <- r4[1] * c(cos(r4[2]*2*pi),sin(r4[2]*2*pi))
Q <- r4[3] * c(cos(r4[4]*2*pi),sin(r4[4]*2*pi))

out1 <- my.hyp.dist1(P,Q)
out2 <- my.hyp.dist2(P,Q)
out1$dhyp1
```

```
## [1] -1.003046
```

```r
out1$dhyp2
```

```
## [1] 1.003046
```

```r
out2
```

```
## [1] 1.003046
```

念のため、測地線円と単位円との交点を描図しておく。


```r
plot(cos(t),sin(t),type="l",asp=1)
points(P[1],P[2],pch=20,col=2)
points(Q[1],Q[2],pch=20,col=3)
points(out1$A[1],out1$A[2],pch=20,col=4)
points(out1$B[1],out1$B[2],pch=20,col=5)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-20-1.png" width="672" />

ちなみに、双曲線模型では、円盤上の２点間距離を、いったん、双曲面上の座標に戻し、
そこでのミンコフスキー内積を計算し、それに基づいて計算していた。

その方法も再実装する。


```r
my.hyp.dist3 <- function(P,Q){
  P.hyp.cds <- my.pdisk.coords.inv(P)
  Q.hyp.cds <- my.pdisk.coords.inv(Q)
  return(my.hyp.dist(P.hyp.cds,Q.hyp.cds))
}
```


一致することを確認する。

```r
n <- 10000
x. <- runif(n,min=-1,max=1)
y. <- runif(n,min=-1,max=1)
r2 <- x.^2 + y.^2
s <- which(r2 < 1)
x. <- x.[s]
y. <- y.[s]
plot(x.,y.,pch=20,cex=0.01,asp=TRUE)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-22-1.png" width="672" />

```r
hypD1 <- apply(cbind(x.,y.),1,my.hyp.dist2,c(x.[1],y.[1]))
hypD2 <- apply(cbind(x.,y.),1,my.hyp.dist3,c(x.[1],y.[1]))

plot(hypD1,hypD2)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-22-2.png" width="672" />


### 円～等距離点集合

ポアンカレ・ディスク上の円は、ある点からの測地線的距離が一定の点の集合。


```r
n <- 10000
x. <- runif(n,min=-1,max=1)
y. <- runif(n,min=-1,max=1)
r2 <- x.^2 + y.^2
s <- which(r2 < 1)
x. <- x.[s]
y. <- y.[s]

hypD1 <- apply(cbind(x.,y.),1,my.hyp.dist2,c(x.[1],y.[1]))

col <- ceiling(hypD1/max(hypD1)*10)
plot(x.,y.,pch=20,cex=1,asp=TRUE,col=col)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-23-1.png" width="672" />


### 円ではないが、円のようなもの

双曲幾何では、いわゆる円の他に、Hypercycle・Hypocycleと言う２つの異なるものが、円と同様に閉じた点集合となる

* 円は、ある点から双曲幾何的距離が等距離の点の集合。
* Hypercycleは双曲幾何的直線から等距離にある点の集合。ただし、点と直線との距離は、直線上の点との距離のうち最小のもののこと。ポアンカレ・ディスクの円周と直交しない円が作る曲線がこれに相当する。
* Hypocycyleは、そのcycleに直行する測地線が漸近的に１点に集まるようなサイクルである。ポアンカレ・ディスクの内部の円で、ディスクに内接する。

### Hypercycle

ある直線から等距離の点の集合は、ユークリッド幾何では、平行な直線になる。
ただし、ユークリッド幾何での平行な直線は、「交点を持たない」直線のことでもある。

双曲幾何では、交わらないことと等距離であることが、同じことにならない。

ある距離を持って等距離な線はある線となり、それは（片側に）１本引けるだけだが、交点を持たない直線は２本以上（たいていは無数に）引ける。

これも視覚化してみよう。

２点を通る測地線を引き、その線との最短距離でディスクを色分けしてみる。

まずは測地線。特に、円板上の点。

```r
P <- rnorm(2) * 1
Q <- rnorm(2) *1
tmp <- my.orth.circle2(P,Q)
r <- tmp$r
ctr <- tmp$ctr
t <- seq(from=0,to=1,length=1000)*2*pi
pnts <- cbind(r*cos(t)+ctr[1],r*sin(t)+ctr[2])
# このうち、円板上の点のみを選ぶ
s <- which(apply(pnts^2,1,sum) < 1)
pnts <- pnts[s,]
plot(cos(t),sin(t),type="l",asp=1)
points(pnts,type="l",col=2)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-24-1.png" width="672" />
円板上の乱点と、測地線上の点との双曲幾何的距離を計算し、その中の最小値を取り出す。



```r
n <- 10000
x. <- runif(n,min=-1,max=1)
y. <- runif(n,min=-1,max=1)
r2 <- x.^2 + y.^2
s <- which(r2 < 1)
x. <- x.[s]
y. <- y.[s]

hypD.geod <- matrix(0,length(s),length(pnts[,1]))

for(i in 1:length(pnts[,1])){
  tmp <- apply(cbind(x.,y.),1,my.hyp.dist2,pnts[i,])
  hypD.geod[,i] <- tmp
}
hypD.geod.min <- apply(hypD.geod,1,min)
```

描図すると、測地線からの等距離点集合は円弧をなすが、それは単位円周と直交していないことがわかる。


```r
col <- ceiling(hypD.geod.min/max(hypD.geod.min)*10)
plot(x.,y.,pch=20,cex=1,asp=TRUE,col=col)
points(pnts,type="l",col=2,lwd=2)
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-26-1.png" width="672" />

### Hypocycle

可視化してみる。

円周上の点を取り、その点で内接する円を描く。これがHypocycle。

Hypocycle上の任意の点と、内接点とを通る測地線(円)が、Hypocycleと直交していることを示すことで、視覚化したものとする。

Hypocycleがしめしているのは、双曲幾何での直線に直交する直線は、互いに交差することがないが、無限遠では１点に集中する、ということである。


```r
# 内接点を定める角度
theta <- runif(1)*2*pi
# 内接円の半径
r <- runif(1)*0.5 + 0.25
# 内接円の中心座標
ctr <- (1-r) * c(cos(theta),sin(theta))

t <- seq(from=0,to=1,length=1000)*2*pi
pnts <- cbind(r*cos(t)+ctr[1],r*sin(t)+ctr[2])

plot(cos(t),sin(t),type="l",asp=1)
points(pnts,type="l",col=2)

n.line <- 10

t.samples <- sample(1:length(t),n.line)


for(i in 1:n.line){
  P <- c(cos(theta),sin(theta)) # 内接点
  Q <- pnts[t.samples[i],]
  tmp <- my.orth.circle2(P,Q)
  r <- tmp$r
  ctr <- tmp$ctr
  tmp.pnts <- cbind(r*cos(t)+ctr[1],r*sin(t)+ctr[2])
  points(tmp.pnts,type="l",col=i+1)
  
}
```

<img src="PoincareDisk_files/figure-html/unnamed-chunk-27-1.png" width="672" />


### ミンコフスキー空間としてとらえる

#### ミンコフスキー空間を双曲面の積み重ねに分解する

この節は、[Minkowski space(Wikipedia記事)のGeometryの節](https://en.wikipedia.org/wiki/Minkowski_space#Geometry)節に
に従って書いてある。

ミンコフスキー空間では、n個の正の固有値とm個の負の固有値によってミンコフスキー内積が定まる。

そのうちn=3,m=1のものは、物理学・相対性理論で重要なものであり、資料が多い。

この$M^{n=3,m=1}$では、内積が
$$
v_x u_x + v_y u_y + v_z u_z - v_t u_t
$$
と計算される。

この空間の各点を距離を保ったままユークリッド空間に埋め込むことはできない。

球面幾何の場合には、次元を一つ上げることで埋め込める（それが球面）のに対して、
双曲幾何では、次元を一つあげてもうまくいかない、ということである。

双曲幾何空間全体を埋め込むことはできないが、双曲幾何空間を

$$
x^2 + y^2 + z^2 - t^2 = R^2; R\in [0,\infty)
$$
という（超）双曲面に分けると、個々の（超）双曲面は、半径Rの円板に埋め込める。

埋め込みのモデルとして、「双曲面モデル」「ポアンカレ・ディスク（ポアンカレ・ボール）モデル」「ポアンカレの半平面モデル」など複数のものが知られている。

この（超）双曲面は、一定負のガウス曲率$-1/R^2$を持った「定負曲率空間」である。

ミンコフスキー空間自体はリーマン多様体ではなく、擬リーマン多様体であるが、（超）双曲面を埋め込んだ各モデルはリーマン多様体になっていて、そのリーマン計量が定まっている。

今、時間軸の座標を$t$、空間軸の座標を$\mathbf(x)$とすると

$$
(t,\mathbf{x}) \to \mathbf{u}= \frac{R \mathbf{x}}{(R+t)}\\
\mathbf{u} \to (t,\mathbf{x}) = (R\frac{R^2+|\mathbf{u}|^2}{R^2-|\mathbf{u}|^2},\frac{2R^2 \mathbf{u}}{R^2-|\mathbf{u}|^2})
$$



これは、双曲面モデル、ポアンカレ・ディスクモデルのところでは、$R=1$に限定して扱ってきた式の、一般化したものに相当する。

ディスクへの写像座標は(0,-R)と双曲面上の点とを結んだ線によって得られる。



```r
# ベクトルxの最後の座標が時間軸とする
my.pdisk.coords.R <- function(x){
  t <- x[length(x)]
  x. <- x[-length(x)]
  R <- sqrt(-sum(x^2) + 2 * t^2)
  u <- R*x./(R+t)
  return(list(u=u,R=R))
}
my.pdisk.coords.inv.R <- function(u,R){
  t <- R * (R^2+sum(u^2))/(R^2-sum(u^2))
  x. <- 2*R^2*u/(R^2-sum(u^2))
  return(c(x.,t))
}
```

念のため、R=1の例について検算する

```r
x <- 0.4
y <- 0.3
z <- sqrt(x^2+y^2+1)
u1 <- my.pdisk.coords(c(x,y,z))
u2 <- my.pdisk.coords.R(c(x,y,z))
print(u1)
```

```
## [1] 0.1888544 0.1416408
```

```r
print(u2)
```

```
## $u
## [1] 0.1888544 0.1416408
## 
## $R
## [1] 1
```

```r
my.pdisk.coords.inv(u1)
```

```
## [1] 0.400000 0.300000 1.118034
```

```r
my.pdisk.coords.inv.R(u2$u,u2$R)
```

```
## [1] 0.400000 0.300000 1.118034
```
