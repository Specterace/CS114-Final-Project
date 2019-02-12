/**
 * The functions in this file create models in an
 * IFS format that can be drawn using gl.drawElements
 * with primitive type gl.TRIANGLES.  Objects have
 * vertex coordinates, normal vectors, and texture
 * coordinates for each vertex, plus a list of indicies
 * for the element array buffer.  The return value
 * of each function is an object, model, with properties:
 * 
 *    model.vertexPositions -- the vertex coordinates;
 *    model.vertexNormals -- the normal vectors;
 *    model.vertexTextureCoords -- the texture coordinates;
 *    model.indices -- the face indices.
 *
 * The first three properties are of type Float32Array, while
 * model.indicesis of type Uint16Array.
 */

 
 /**
  * Create a model of a cube, centered at the origin.
  * @side the length of a side of the cube.  If not given, the value will be 1.
  */
function cube(side) {
   var s = (side || 1)/2;
   var coords = [];
   var normals = [];
   var texCoords = [];
   var indices = [];
   function face(xyz, nrm) {
      var start = coords.length/3;
      var i;
      for (i = 0; i < 12; i++) {
         coords.push(xyz[i]);
      }
      for (i = 0; i < 4; i++) {
         normals.push(nrm[0],nrm[1],nrm[2]);
      }
      texCoords.push(0,0,1,0,1,1,0,1);
      indices.push(start,start+1,start+2,start,start+2,start+3);
   }
   face( [-s,-s,s, s,-s,s, s,s,s, -s,s,s], [0,0,1] );
   face( [-s,-s,-s, -s,s,-s, s,s,-s, s,-s,-s], [0,0,-1] );
   face( [-s,s,-s, -s,s,s, s,s,s, s,s,-s], [0,1,0] );
   face( [-s,-s,-s, s,-s,-s, s,-s,s, -s,-s,s], [0,-1,0] );
   face( [s,-s,-s, s,s,-s, s,s,s, s,-s,s], [1,0,0] );
   face( [-s,-s,-s, -s,-s,s, -s,s,s, -s,s,-s], [-1,0,0] );
   return {
      vertexPositions: new Float32Array(coords),
      vertexNormals: new Float32Array(normals),
      vertexTextureCoords: new Float32Array(texCoords),
      indices: new Uint16Array(indices)
   }
}

/////////////////////////////////////////////////////////////////////START OF WATER MESH CODE
/*
//////////NOTE: This water mesh/heightmap implementation was specifically coded and created from scratch for the 114 Final Project
*/
// variables that were originally created with the intent of implementing changeable vertex positions and therefore moving water surface
var WMVertexPosition, WMVertexNormals;
var WMVertexVelocity, WMVertexForce;
var WMDimensions;
var waterDrawn = 0;
var timeStep = 0.01;
var timetotal = 0.0;
var per = 2.0*Math.PI;

////////////////////////////////////////////////////////////////Getters/Setters
/////////These functions were originally intended to be used for implementation of moving water surface
/////////Though the final water mesh implementation still uses many of these, they were not used to their initial potential
function getPosition(i, j, mRes) {
    var id = i*mRes + j;
	var pos = vec3.create();
	pos[0] = WMVertexPosition[3*id];
	pos[1] = WMVertexPosition[3*id + 1];
	pos[2] = WMVertexPosition[3*id + 2];
	return pos;
}

function setPosition(i, j, x, mRes) {
    var id = i*mRes + j;
    WMVertexPosition[3*id] = x[0]; WMVertexPosition[3*id + 1] = x[1]; WMVertexPosition[3*id + 2] = x[2];
}

function getNormal(i, j, mRes) {
    var id = i*mRes + j;
	var norm = vec3.create();
	norm[0] = WMVertexNormals[3*id];
	norm[1] = WMVertexNormals[3*id + 1];
	norm[2] = WMVertexNormals[3*id + 2];
	return norm;
}

function getVelocity(i, j, mRes) {
    var id = i*mRes + j;
	var velo = vec3.create();
	velo[0] = WMVertexVelocity[id][0];
	velo[1] = WMVertexVelocity[id][1];
	velo[2] = WMVertexVelocity[id][2];
	return velo;
}

function setVelocity(i, j, v, mRes) {
    var id = i*mRes + j;
	var veloSet = vec3.create();
	veloSet[0] = v[0];
	veloSet[1] = v[1];
	veloSet[2] = v[2];
	WMVertexVelocity[id] = veloSet;
}

function crossProd(a, b) {
	var d = a[0], e = a[1], k = a[2];
	var g = b[0], f = b[1], h = b[2];
	var c = vec3.create();
	c[0] = e * h - k * f;
	c[1] = k * g - d * h;
	c[2] = d * f - e * g;
	return c;
}



/////////////////////////////////////////////////////////////Computing normals

function computeNormals(lCount, wCount, mRes) {
	var dxa = [1, 0, 0, 1, 1, 0], dya = [1, 1, 0, 0, 1, 0];		//offset to find vertices for the "a" vectors (of a x b)
	var dxb = [0, 0, 1, 1, 0, 1], dyb = [1, 0, 1, 1, 0, 0];		//offset to find vertices for the "b" vectors (of a x b)
	var odx = [0, 1, 0, 0, 1, 1], ody = [0, 1, 1, 0, 0, 1];		//add this to the origin points to calculate proper vertices
	
    var e1, e2, e3, c1, c2, c3, n1, p0;
	var i1, i2, j1, j2;
    var i, j, t, k = 0;
	var iO, jO;
	c1 = vec3.create();
	c2 = vec3.create();
	c3 = vec3.create();
	e1 = vec3.create();
	e2 = vec3.create();
	e3 = vec3.create();
	n1 = vec3.create();
    for ( i = 0; i < (lCount - 1); ++i )
        for ( j = 0; j < (wCount - 1); ++j ) {
			var norms = [];
			
			p0 = vec3.create();
			
			for (t = 0; t < 6; ++t) {
				p0 = getPosition((i + ody[t]),(j + odx[t]),mRes);
				i1 = (i + dya[t]); 
				j1 = (j + dxa[t]);
				i2 = (i + dyb[t]);
				j2 = (j + dxb[t]);
				e1 = getPosition(i1, j1, mRes);
                e2 = getPosition(i2, j2, mRes);
				
				c1[0] = e1[0] - p0[0];
				c1[1] = e1[1] - p0[1];
				c1[2] = e1[2] - p0[2];
				c2[0] = e2[0] - p0[0];
				c2[1] = e2[1] - p0[1];
				c2[2] = e2[2] - p0[2];
				
				c3 = crossProd(c1, c2);
				n1 = vec3.normalize(n1, c3);
				norms.push(n1);
			}			
			
            for ( t = 0; t < 6; ++t ) {
				WMVertexNormals[((18*k) + (3*t))] = norms[t][0];
				WMVertexNormals[((18*k) + (3*t)) + 1] = norms[t][1];
				WMVertexNormals[((18*k) + (3*t)) + 2] = norms[t][2];
			}
            ++k;
        }
}

/////////////////////////////////////////////////////////////Creating the actual mesh

function waterMesh(xSide, ySide, meshResolution) {
   //meshResolution is how many vertices per 1.0 length on a side do you want to render
   //xSide and ySide must be integers
   var coords = [];
   var normals = [];
   var texCoords = [];
   var indices = [];
   var i, j, k;
   var lVertCount, wVertCount;
   var vRes;
   
   lVertCount = (xSide * meshResolution) + 1;		//number of vertices on length side (how many rows)
   wVertCount = (ySide * meshResolution) + 1;		//number of vertices on width side	(how many columns)
   vRes = (1.0 / meshResolution);			//length between vertices horizontally and vertically
  
   if (waterDrawn == 0) {
		WMVertexPosition = new Array(lVertCount * wVertCount*3);
		WMVertexNormals = new Array(lVertCount * wVertCount*9);		//modify number of values for this normal array
		WMVertexVelocity = new Array(lVertCount * wVertCount);
		WMVertexForce = new Array(lVertCount * wVertCount*3);
		WMDimensions = new Array(2);
		WMDimensions[0] = lVertCount;
		WMDimensions[1] = wVertCount;
   }
   
   texCoords.push(0,0,1,0,1,1,0,1);
   
   //specifying what position the triangle vertices are in
   timetotal = timetotal + timeStep;
   
	for (i = 0; i < lVertCount; ++i) {
		for (j = 0; j < wVertCount; ++j) {
				setPosition(i, j, [(j*vRes), (i*vRes), 0.0], wVertCount);			//set position vertex
				setVelocity(i, j, vec3.create(), wVertCount);
		}
	}
	
	//compute normals
	computeNormals(lVertCount, wVertCount, wVertCount);

	k = 0;		//specifying what vertices make the triangles of the mesh  (consider the -1 for vertex count)
	for (i = 0; i < lVertCount - 1; ++i) {
		for (j = 0; j < wVertCount - 1; ++j) {
			indices[6*k] = i*wVertCount + j;
			indices[6*k + 1] = (i + 1)*wVertCount + j + 1;
			indices[6*k + 2] = (i + 1)*wVertCount + j;
			indices[6*k + 3] = i*wVertCount + j;
			indices[6*k + 4] = i*wVertCount + j + 1;
			indices[6*k + 5] = (i + 1)*wVertCount + j + 1;
			++k;
		}
	}
   return {
	   
	  //convert vertex positions and normals into a format that can be used by these functions
      vertexPositions: new Float32Array(WMVertexPosition),
      vertexNormals: new Float32Array(WMVertexNormals),
	  
	  //CREATE A VERTEXVELOCITY ARRAY
	  //vertexVelocity: new Float32Array(WMVertexVelocity),
	  
      vertexTextureCoords: new Float32Array(texCoords),
      indices: new Uint16Array(indices)
   }
}

///////////////////An attempt to write functions to implement periodic movement for the water mesh.
/*
function waterMove(timeStep) {
	var tempVec1 = vec3.create();
	var tempVec2 = vec3.create();
	var forceVec = vec3.create();
	var id;
	var length = WMDimensions[0];		//for i, or rows
	var width = WMDimensions[1];		//for j, or columns
	for (i = 0; i < length; ++i) {
		for (j = 0; j < width; ++j) {
			forceVec = calculateWaterForce(i,j);
			id = (i*width) + j;
			WMVertexForce[3*id] = forceVec[0];
			WMVertexForce[(3*id) + 1] = forceVec[1];
			WMVertexForce[(3*id) + 2] = forceVec[2];
		}
	}
	
	for (i = 0; i < length; ++i) {
		for (j = 0; j < width; ++j) {
			id = (i*width) + j;
			//calculating and setting velocity
			tempVec1[0] = WMVertexForce[3*id] * timeStep;
			tempVec1[1] = WMVertexForce[3*id + 1] * timeStep;
			tempVec1[2] = WMVertexForce[3*id + 2] * timeStep;
			tempVec2 = getVelocity(i, j, width);
			tempVec1[0] = tempVec1[0] + tempVec2[0];
			tempVec1[1] = tempVec1[1] + tempVec2[1];
			tempVec1[2] = tempVec1[2] + tempVec2[2];
			setVelocity(i,j,tempVec1,width);	
			
			//calculating and setting position
			tempVec1 = getVelocity(i,j,width);
			tempVec2 = getPosition(i,j,width);
			tempVec1[0] = tempVec1[0] * timeStep;
			tempVec1[1] = tempVec1[1] * timeStep;
			tempVec1[2] = tempVec1[2] * timeStep;
			tempVec2[0] = tempVec2[0] + tempVec1[0];
			tempVec2[1] = tempVec2[1] + tempVec1[1];
			tempVec2[2] = tempVec2[2] + tempVec1[2];
			setPosition(i,j,tempVec2,width);
		}
	}
	computeNormals(length, width, width);
}

function calculateWaterForce(i, j) {
	var period = 2*Math.PI;
	var force = vec3.create();
	force[0] = 0.0;
	force[1] = 0.0;
	//var xPer = (period * (j / 10.0)) + (period * (timeStep));
	//var yPer = (period * (i / 10.0)) + (period * (timeStep));
	var xPer = (period * (j / 10.0));
	var yPer = (period * (i / 10.0));
	var zPos = 10.0 * (Math.sin(yPer) + Math.cos(xPer));
	force[2] = zPos;
	return force;
}

*/


/////////////////////////////////////////////////////////////////////////////////////////END WATER MESH CODE


/**
 * Creates a model of an annulus or disk lying in the xy plane,
 * centered at the origin.
 * @param innerRadius the radius of the hole in the radius; a value of
 *    zero will give a disk rather than a ring.  If not present,
 *    the default value is 0.25.
 * @param outerRadius the radius of the ring, from the center to the
 *    outer edge.  Must be greater than innerRadius.  If not provided,
 *    the default value is 2*innerRadius or is 0.5 if innerRadius is 0.
 * @slices the number of radial subdivisions in the circular approximation
 *    of an annulus.  If not provided, the value will be 32.
 */
function ring(innerRadius, outerRadius, slices) {
   if (arguments.length == 0)
      innerRadius = 0.25;
   outerRadius = outerRadius || innerRadius * 2 || 0.5;
   slices = slices || 32;
   var vertexCount, vertices, normals, texCoords, indices, i;
   vertexCount = (innerRadius == 0)? slices + 1 : slices * 2;
   vertices = new Float32Array( 3*vertexCount );
   normals = new Float32Array( 3* vertexCount );
   texCoords = new Float32Array( 2*vertexCount );
   indices = new Uint16Array( innerRadius == 0 ?  3*slices : 3*2*slices );
   var d = 2*Math.PI/slices;
   var k = 0;
   var t = 0;
   var n = 0;
   if (innerRadius == 0) {
      for (i = 0; i < slices; i++) {
         c = Math.cos(d*i);
         s = Math.sin(d*i);
         vertices[k++] = c*outerRadius;
         vertices[k++] = s*outerRadius;
         vertices[k++] = 0;
         texCoords[t++] = 0.5 + 0.5*c;
         texCoords[t++] = 0.5 + 0.5*s;
         indices[n++] = slices;
         indices[n++] = i;
         indices[n++] = i == slices-1 ? 0 : i + 1;
      }
      vertices[k++] = vertices[k++] = vertices[k++] = 0;
      texCoords[t++] = texCoords[t++] = 0;
   }
   else {
      var r = innerRadius / outerRadius;
      for (i = 0; i < slices; i++) {
         c = Math.cos(d*i);
         s = Math.sin(d*i);
         vertices[k++] = c*innerRadius;
         vertices[k++] = s*innerRadius;
         vertices[k++] = 0;
         texCoords[t++] = 0.5 + 0.5*c*r;
         texCoords[t++] = 0.5 + 0.5*s*r;
         vertices[k++] = c*outerRadius;
         vertices[k++] = s*outerRadius;
         vertices[k++] = 0;
         texCoords[t++] = 0.5 + 0.5*c;
         texCoords[t++] = 0.5 + 0.5*s;
      }
      for (i = 0; i < slices - 1; i++) {
         indices[n++] = 2*i;
         indices[n++] = 2*i+1;
         indices[n++] = 2*i+3;
         indices[n++] = 2*i;
         indices[n++] = 2*i+3;
         indices[n++] = 2*i+2;
      }
      indices[n++] = 2*i;
      indices[n++] = 2*i+1;
      indices[n++] = 1;
      indices[n++] = 2*i;
      indices[n++] = 1;
      indices[n++] = 0;
   }
   for (i = 0; i < vertexCount; i++) {
      normals[3*i] = normals[3*i+1] = 0;
      normals[3*i+2] = 1;
   }
   return {
       vertexPositions: vertices,
       vertexNormals: normals,
       vertexTextureCoords: texCoords,
       indices: indices
   };
}

/**
 * Create a model of a sphere.  The z-axis is the axis of the sphere,
 * with the north pole on the positive z-axis and the center at (0,0,0).
 * @param radius the radius of the sphere, default 0.5 if not specified.
 * @param slices the number of lines of longitude, default 32
 * @param stacks the number of lines of latitude plus 1, default 16.  (This 
 *    is the number of vertical slices, bounded by lines of latitude, the
 *    north pole and the south pole.)
 */
function uvSphere(radius, slices, stacks) {
   radius = radius || 0.5;
   slices = slices || 32;
   stacks = stacks || 16;
   var vertexCount = (slices+1)*(stacks+1);
   var vertices = new Float32Array( 3*vertexCount );
   var normals = new Float32Array( 3* vertexCount );
   var texCoords = new Float32Array( 2*vertexCount );
   var indices = new Uint16Array( 2*slices*stacks*3 );
   var du = 2*Math.PI/slices;
   var dv = Math.PI/stacks;
   var i,j,u,v,x,y,z;
   var indexV = 0;
   var indexT = 0;
   for (i = 0; i <= stacks; i++) {
      v = -Math.PI/2 + i*dv;
      for (j = 0; j <= slices; j++) {
         u = j*du;
         x = Math.cos(u)*Math.cos(v);
         y = Math.sin(u)*Math.cos(v);
         z = Math.sin(v);
         vertices[indexV] = radius*x;
         normals[indexV++] = x;
         vertices[indexV] = radius*y;
         normals[indexV++] = y;
         vertices[indexV] = radius*z;
         normals[indexV++] = z;
         texCoords[indexT++] = j/slices;
         texCoords[indexT++] = i/stacks;
      } 
   }
   var k = 0;
   for (j = 0; j < stacks; j++) {
      var row1 = j*(slices+1);
      var row2 = (j+1)*(slices+1);
      for (i = 0; i < slices; i++) {
          indices[k++] = row1 + i;
          indices[k++] = row2 + i + 1;
          indices[k++] = row2 + i;
          indices[k++] = row1 + i;
          indices[k++] = row1 + i + 1;
          indices[k++] = row2 + i + 1;
      }
   }
   return {
       vertexPositions: vertices,
       vertexNormals: normals,
       vertexTextureCoords: texCoords,
       indices: indices
   };
}


/**
 * Create a model of a torus (surface of a doughnut).  The z-axis goes through the doughnut hole,
 * and the center of the torus is at (0,0,0).
 * @param outerRadius the distance from the center to the outside of the tube, 0.5 if not specified.
 * @param innerRadius the distance from the center to the inside of the tube, outerRadius/3 if not
 *    specified.  (This is the radius of the doughnut hole.)
 * @param slices the number of lines of longitude, default 32.  These are slices parallel to the
 * z-axis and go around the tube the short way (through the hole).
 * @param stacks the number of lines of latitude plus 1, default 16.  These lines are perpendicular
 * to the z-axis and go around the tube the long way (arouind the hole).
 */
function uvTorus(outerRadius, innerRadius, slices, stacks) {
   outerRadius = outerRadius || 0.5;
   innerRadius = innerRadius || outerRadius/3;
   slices = slices || 32;
   stacks = stacks || 16;
   var vertexCount = (slices+1)*(stacks+1);
   var vertices = new Float32Array( 3*vertexCount );
   var normals = new Float32Array( 3* vertexCount );
   var texCoords = new Float32Array( 2*vertexCount );
   var indices = new Uint16Array( 2*slices*stacks*3 );
   var du = 2*Math.PI/slices;
   var dv = 2*Math.PI/stacks;
   var centerRadius = (innerRadius+outerRadius)/2;
   var tubeRadius = outerRadius - centerRadius;
   var i,j,u,v,cx,cy,sin,cos,x,y,z;
   var indexV = 0;
   var indexT = 0;
   for (j = 0; j <= stacks; j++) {
      v = -Math.PI + j*dv;
      cos = Math.cos(v);
      sin = Math.sin(v);
      for (i = 0; i <= slices; i++) {
         u = i*du;
         cx = Math.cos(u);
         cy = Math.sin(u);
         x = cx*(centerRadius + tubeRadius*cos);
         y = cy*(centerRadius + tubeRadius*cos);
         z = sin*tubeRadius;
         vertices[indexV] = x;
         normals[indexV++] = cx*cos;
         vertices[indexV] = y
         normals[indexV++] = cy*cos;
         vertices[indexV] = z
         normals[indexV++] = sin;
         texCoords[indexT++] = i/slices;
         texCoords[indexT++] = j/stacks;
      } 
   }
   var k = 0;
   for (j = 0; j < stacks; j++) {
      var row1 = j*(slices+1);
      var row2 = (j+1)*(slices+1);
      for (i = 0; i < slices; i++) {
          indices[k++] = row1 + i;
          indices[k++] = row2 + i + 1;
          indices[k++] = row2 + i;
          indices[k++] = row1 + i;
          indices[k++] = row1 + i + 1;
          indices[k++] = row2 + i + 1;
      }
   }
   return {
       vertexPositions: vertices,
       vertexNormals: normals,
       vertexTextureCoords: texCoords,
       indices: indices
   };
}

/**
 * Defines a model of a cylinder.  The axis of the cylinder is the z-axis,
 * and the center is at (0,0,0).
 * @param radius the radius of the cylinder
 * @param height the height of the cylinder.  The cylinder extends from -height/2
 * to height/2 along the z-axis.
 * @param slices the number of slices, like the slices of an orange.
 * @param noTop if missing or false, the cylinder has a top; if set to true,
 *   the cylinder has a top. The top is a disk at the positive end of the cylinder.
 * @param noBottom if missing or false, the cylinder has a bottom; if set to true,
 *   the cylinder has a bottom. The bottom is a disk at the negtive end of the cylinder.
 */
function uvCylinder(radius, height, slices, noTop, noBottom) {
   radius = radius || 0.5;
   height = height || 2*radius;
   slices = slices || 32;
   var vertexCount = 2*(slices+1);
   if (!noTop)
      vertexCount += slices + 2;
   if (!noBottom)
      vertexCount += slices + 2;
   var triangleCount = 2*slices;
   if (!noTop)
      triangleCount += slices;
   if (!noBottom)
      triangleCount += slices; 
   var vertices = new Float32Array(vertexCount*3);
   var normals = new Float32Array(vertexCount*3);
   var texCoords = new Float32Array(vertexCount*2);
   var indices = new Uint16Array(triangleCount*3);
   var du = 2*Math.PI / slices;
   var kv = 0;
   var kt = 0;
   var k = 0;
   var i,u;
   for (i = 0; i <= slices; i++) {
      u = i*du;
      var c = Math.cos(u);
      var s = Math.sin(u);
      vertices[kv] = c*radius;
      normals[kv++] = c;
      vertices[kv] = s*radius;
      normals[kv++] = s;
      vertices[kv] = -height/2;
      normals[kv++] = 0;
      texCoords[kt++] = i/slices;
      texCoords[kt++] = 0;
      vertices[kv] = c*radius;
      normals[kv++] = c;
      vertices[kv] = s*radius;
      normals[kv++] = s;
      vertices[kv] = height/2;
      normals[kv++] = 0;
      texCoords[kt++] = i/slices;
      texCoords[kt++] = 1;
   }
   for (i = 0; i < slices; i++) {
          indices[k++] = 2*i;
          indices[k++] = 2*i+3;
          indices[k++] = 2*i+1;
          indices[k++] = 2*i;
          indices[k++] = 2*i+2;
          indices[k++] = 2*i+3;
   }
   var startIndex = kv/3;
   if (!noBottom) {
      vertices[kv] = 0;
      normals[kv++] = 0;
      vertices[kv] = 0;
      normals[kv++] = 0;
      vertices[kv] = -height/2;
      normals[kv++] = -1;
      texCoords[kt++] = 0.5;
      texCoords[kt++] = 0.5; 
      for (i = 0; i <= slices; i++) {
         u = 2*Math.PI - i*du;
         var c = Math.cos(u);
         var s = Math.sin(u);
         vertices[kv] = c*radius;
         normals[kv++] = 0;
         vertices[kv] = s*radius;
         normals[kv++] = 0;
         vertices[kv] = -height/2;
         normals[kv++] = -1;
         texCoords[kt++] = 0.5 - 0.5*c;
         texCoords[kt++] = 0.5 + 0.5*s;
      }
      for (i = 0; i < slices; i++) {
         indices[k++] = startIndex;
         indices[k++] = startIndex + i + 1;
         indices[k++] = startIndex + i + 2;
      }
   }
   var startIndex = kv/3;
   if (!noTop) {
      vertices[kv] = 0;
      normals[kv++] = 0;
      vertices[kv] = 0;
      normals[kv++] = 0;
      vertices[kv] = height/2;
      normals[kv++] = 1;
      texCoords[kt++] = 0.5;
      texCoords[kt++] = 0.5; 
      for (i = 0; i <= slices; i++) {
         u = i*du;
         var c = Math.cos(u);
         var s = Math.sin(u);
         vertices[kv] = c*radius;
         normals[kv++] = 0;
         vertices[kv] = s*radius;
         normals[kv++] = 0;
         vertices[kv] = height/2;
         normals[kv++] = 1;
         texCoords[kt++] = 0.5 + 0.5*c;
         texCoords[kt++] = 0.5 + 0.5*s;
      }
      for (i = 0; i < slices; i++) {
         indices[k++] = startIndex;
         indices[k++] = startIndex + i + 1;
         indices[k++] = startIndex + i + 2;
      }
   }
   return {
       vertexPositions: vertices,
       vertexNormals: normals,
       vertexTextureCoords: texCoords,
       indices: indices
   };
}


/**
 * Defines a model of a cone.  The axis of the cone is the z-axis,
 * and the center is at (0,0,0).
 * @param radius the radius of the cone
 * @param height the height of the cone.  The cone extends from -height/2
 * to height/2 along the z-axis, with the tip at (0,0,height/2).
 * @param slices the number of slices, like the slices of an orange.
 * @param noBottom if missing or false, the cone has a bottom; if set to true,
 *   the cone has a bottom. The bottom is a disk at the wide end of the cone.
 */
function uvCone(radius, height, slices, noBottom) {
   radius = radius || 0.5;
   height = height || 2*radius;
   slices = slices || 32;
   var fractions = [ 0, 0.5, 0.75, 0.875, 0.9375 ];
   var vertexCount = fractions.length*(slices+1) + slices;
   if (!noBottom)
      vertexCount += slices + 2;
   var triangleCount = (fractions.length-1)*slices*2 + slices;
   if (!noBottom)
      triangleCount += slices;
   var vertices = new Float32Array(vertexCount*3);
   var normals = new Float32Array(vertexCount*3);
   var texCoords = new Float32Array(vertexCount*2);
   var indices = new Uint16Array(triangleCount*3);
   var normallength = Math.sqrt(height*height+radius*radius);
   var n1 = height/normallength;
   var n2 = radius/normallength; 
   var du = 2*Math.PI / slices;
   var kv = 0;
   var kt = 0;
   var k = 0;
   var i,j,u;
   for (j = 0; j < fractions.length; j++) {
      var uoffset = (j % 2 == 0? 0 : 0.5);
      for (i = 0; i <= slices; i++) {
         var h1 = -height/2 + fractions[j]*height;
         u = (i+uoffset)*du;
         var c = Math.cos(u);
         var s = Math.sin(u);
         vertices[kv] = c*radius*(1-fractions[j]);
         normals[kv++] = c*n1;
         vertices[kv] = s*radius*(1-fractions[j]);
         normals[kv++] = s*n1;
         vertices[kv] = h1;
         normals[kv++] = n2;
         texCoords[kt++] = (i+uoffset)/slices;
         texCoords[kt++] = fractions[j];
      }
   }
   var k = 0;
   for (j = 0; j < fractions.length-1; j++) {
      var row1 = j*(slices+1);
      var row2 = (j+1)*(slices+1);
      for (i = 0; i < slices; i++) {
          indices[k++] = row1 + i;
          indices[k++] = row2 + i + 1;
          indices[k++] = row2 + i;
          indices[k++] = row1 + i;
          indices[k++] = row1 + i + 1;
          indices[k++] = row2 + i + 1;
      }
   }
   var start = kv/3 - (slices+1);
   for (i = 0; i < slices; i++) { // slices points at top, with different normals, texcoords
      u = (i+0.5)*du;
      var c = Math.cos(u);
      var s = Math.sin(u);
      vertices[kv] = 0;
      normals[kv++] = c*n1;
      vertices[kv] = 0;
      normals[kv++] = s*n1;
      vertices[kv] = height/2;
      normals[kv++] = n2;
      texCoords[kt++] = (i+0.5)/slices;
      texCoords[kt++] = 1;
   }
   for (i = 0; i < slices; i++) {
      indices[k++] = start+i;
      indices[k++] = start+i+1;
      indices[k++] = start+(slices+1)+i;
   }
   if (!noBottom) {
      var startIndex = kv/3;
      vertices[kv] = 0;
      normals[kv++] = 0;
      vertices[kv] = 0;
      normals[kv++] = 0;
      vertices[kv] = -height/2;
      normals[kv++] = -1;
      texCoords[kt++] = 0.5;
      texCoords[kt++] = 0.5; 
      for (i = 0; i <= slices; i++) {
         u = 2*Math.PI - i*du;
         var c = Math.cos(u);
         var s = Math.sin(u);
         vertices[kv] = c*radius;
         normals[kv++] = 0;
         vertices[kv] = s*radius;
         normals[kv++] = 0;
         vertices[kv] = -height/2;
         normals[kv++] = -1;
         texCoords[kt++] = 0.5 - 0.5*c;
         texCoords[kt++] = 0.5 + 0.5*s;
      }
      for (i = 0; i < slices; i++) {
         indices[k++] = startIndex;
         indices[k++] = startIndex + i + 1;
         indices[k++] = startIndex + i + 2;
      }
   } 
   return {
       vertexPositions: vertices,
       vertexNormals: normals,
       vertexTextureCoords: texCoords,
       indices: indices
   };   
}

