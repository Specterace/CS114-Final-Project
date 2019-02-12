"use strict";

var gl;   // The webgl context.

var a_coords_loc;       // Location of the a_coords attribute variable in the shader program.
var a_normal_loc;       // Location of a_normal attribute.

var a_coords_buffer;    // Buffer for a_coords.
var a_normal_buffer;    // Buffer for a_normal.
var index_buffer;       // Buffer for indices.

var u_diffuseColor;     // Locations of uniform variables in the shader program
var u_specularColor;
var u_specularExponent;

//used to record color values of various objects in 114 final project
var u_rSphereColor;
var u_ySphereColor;
var u_bSphereColor;
var u_oSphereColor;
var u_wallColor;

//old 112 codebase/environment necesary variables
var u_modelview;
var u_projection;
var u_normalMatrix;

//used to record normals of walls and floor of pool in 114 final project
var u_pFloorNMatrix;
var u_nWallNMatrix;
var u_sWallNMatrix;
var u_eWallNMatrix;
var u_wWallNMatrix; 

//used to record reference positions of various objects in 114 final project
var u_sunPosition;
var u_redSpherePos;
var u_yellowSpherePos;
var u_blueSpherePos;
var u_orangeSpherePos;


//pool wall positions (to calculate if the projection hits the walls)
var u_poolNWAPos;	//A is top left
var u_poolNWBPos;	//B is top right
var u_poolNWCPos;	//C is bottom right
var u_poolNWDPos;	//D is bottom left
var u_poolEWBPos;
var u_poolEWCPos;
var u_poolSWBPos;
var u_poolSWCPos;

//various positions to record plane geometry
var u_poolFloorPos;
var u_nWallPos;
var u_sWallPos;
var u_eWallPos;
var u_wWallPos;

// flags for various conditions in 114 final project
var u_isDay;
var u_isLightSource;
var u_isWater;
var u_isUnderwater;

var projection = mat4.create();    // projection matrix
var modelview;                     // modelview matrix; value comes from rotator
var normalMatrix = mat3.create();  // matrix, derived from modelview matrix, for transforming normal vectors

//mat3 objects used to save normals of wall and floor of pool in 114 final project
var floorMatrix = mat3.create();
var nWallMatrix = mat3.create();
var sWallMatrix = mat3.create();
var eWallMatrix = mat3.create();
var wWallMatrix = mat3.create();


var rotator;                  // A TrackballRotator to implement rotation by mouse.

var mvMatrixStack = [];

//Various variables used for the specific 114 final project scene
var timeInt = 0;
var rotateDeg = 0;
var daylight = 1;

var heightDiff = 1.0;
var xDiff = 0.0;
var yDiff = 0.0;

/*
//was to be used for water movement/animation
var timeStep = 0.01;
*/

var objects = [                     // Objects for display, selected by popup menu
    uvTorus(2,1,64,32),
    uvCylinder(1.0,2.0),
    uvSphere(1),
	cube(1),
	ring(1.3,2,32),
	uvCone(1,3,32),
	waterMesh(2.0,2.0,10),
];

var currentModelNumber;             // contains data for the current object

//translate function
function translate(oMtx, iMtx, tVec)
{
	return mat4.translate(oMtx, iMtx, tVec);
}

//rotate function
function rotate(oMtx, iMtx, aRad, rAxis)
{
	return mat4.rotate(oMtx, iMtx, aRad, rAxis);
}

//scale function
function scale(oMtx, iMtx, sVec)
{
	return mat4.scale(oMtx, iMtx, sVec);
}

function degToRad(degrees) {
	return degrees * Math.PI / 180;
}

/**
* Routine for pushing a current model view matrix to a stack for hieroarchial modeling
* @return None
*/
function mvPushMatrix() {
    var copy = mat4.clone(modelview);
    mvMatrixStack.push(copy);
}

/**
* Routine for popping a stored model view matrix from stack for hieroarchial modeling
* @return None
*/
function mvPopMatrix() {
    if (mvMatrixStack.length == 0) {
    	throw "Invalid popMatrix!";
    }
    modelview = mvMatrixStack.pop();
}

/**
* Routine for updating all buffers, effectively "drawing" the shape by "drawing" the given vectors/faces/normals/etc. and sending them to the shaders
* @return None
*/
function update_uniform(modelview,projection,currentModelNumber){

    /* Get the matrix for transforming normal vectors from the modelview matrix,
       and send matrices to the shader program*/
    mat3.normalFromMat4(normalMatrix, modelview);
    
    gl.uniformMatrix3fv(u_normalMatrix, false, normalMatrix);
    gl.uniformMatrix4fv(u_modelview, false, modelview );
    gl.uniformMatrix4fv(u_projection, false, projection );   
    gl.drawElements(gl.TRIANGLES, objects[currentModelNumber].indices.length, gl.UNSIGNED_SHORT, 0);
}

//function to draw the water tank at a given position (component-by-component)
//114 final project new content
function drawWaterTank(){
	
	//floor
	mvPushMatrix();	//drawing base of tank
	//wall stuff
	//place the point of reference for the floor at the floor's middle
	mat3.normalFromMat4(floorMatrix, modelview);
    gl.uniformMatrix3fv(u_pFloorNMatrix, false, floorMatrix);		//saving normalmatrix for pool floor
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,0.1,1],modelview);
	gl.uniform4fv(u_poolFloorPos,transformedPos);					//saving position of pool floor
	scale(modelview, modelview,[6.0, 6.0, 0.2]);
	installModel(objects[3]);
	currentModelNumber = 3;
	gl.uniform4f(u_diffuseColor, 0, 0.5, 0.4, 1);
	update_uniform(modelview,projection,currentModelNumber);
	gl.uniform1i(u_isLightSource, 0);
	mvPopMatrix();
	
	//north wall
	mvPushMatrix();
	translate(modelview, modelview, [0,2.9,1.1]);
	rotate(modelview, modelview, ((Math.PI * 90.0) / 180.0),[1,0,0]);
	//wall stuff
	mat3.normalFromMat4(nWallMatrix, modelview);
    gl.uniformMatrix3fv(u_nWallNMatrix, false, nWallMatrix);		//saving normalmatrix for north wall
	//saving position of north wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,0.1,1],modelview);
	gl.uniform4fv(u_nWallPos,transformedPos);
	//Getting pool wall boundaries relative to north wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [-2.8,1.0,0.1,1],modelview);		//getting position of topleft of north wall
	gl.uniform4fv(u_poolNWAPos,transformedPos);					//saving position of topleft of north wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [2.8,1.0,0.1,1],modelview);		//getting position of topright of north wall
	gl.uniform4fv(u_poolNWBPos,transformedPos);					//saving position of topright of north wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [2.8,-0.9,0.1,1],modelview);		//getting position of bottomright of north wall
	gl.uniform4fv(u_poolNWCPos,transformedPos);					//saving position of bottomright of north wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [-2.8,-0.9,0.1,1],modelview);		//getting position of bottomleft of north wall
	gl.uniform4fv(u_poolNWDPos,transformedPos);					//saving position of bottomleft of north wall
	scale(modelview, modelview,[6.0, 2.0, 0.2]);
	installModel(objects[3]);
	currentModelNumber = 3;
	gl.uniform4f(u_diffuseColor, 0, 0.5, 0.4, 1);
	update_uniform(modelview,projection,currentModelNumber);
	gl.uniform1i(u_isLightSource, 0);
	mvPopMatrix();
	
	//east wall
	mvPushMatrix();
	translate(modelview, modelview, [2.9,0,1.1]);
	rotate(modelview, modelview, ((Math.PI * 90.0) / 180.0),[1,0,0]);
	rotate(modelview, modelview, ((Math.PI * 90.0) / 180.0),[0,1,0]);
	//wall stuff
	mat3.normalFromMat4(eWallMatrix, modelview);
    gl.uniformMatrix3fv(u_eWallNMatrix, false, eWallMatrix);		//saving normalmatrix for east wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,-0.1,1],modelview);
	gl.uniform4fv(u_eWallPos,transformedPos);					//saving position of east wall
	//Getting pool wall boundaries relative to east wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [-2.8,1.0,-0.1,1],modelview);		//getting position of topright of east wall
	gl.uniform4fv(u_poolEWBPos,transformedPos);					//saving position of topright of north wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [-2.8,-0.9,-0.1,1],modelview);		//getting position of bottomright of east wall
	gl.uniform4fv(u_poolEWCPos,transformedPos);					//saving position of bottomright of east wall
	scale(modelview, modelview,[6.0, 2.0, 0.2]);
	installModel(objects[3]);
	currentModelNumber = 3;
	gl.uniform4f(u_diffuseColor, 0, 0.5, 0.4, 1);
	update_uniform(modelview,projection,currentModelNumber);
	gl.uniform1i(u_isLightSource, 0);
	mvPopMatrix();
	
	//south wall
	mvPushMatrix();
	translate(modelview, modelview, [0,-2.9,1.1]);
	rotate(modelview, modelview, ((Math.PI * 90.0) / 180.0),[1,0,0]);
	//wall stuff
	mat3.normalFromMat4(sWallMatrix, modelview);
    gl.uniformMatrix3fv(u_sWallNMatrix, false, sWallMatrix);		//saving normalmatrix for south wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,-0.1,1],modelview);		
	gl.uniform4fv(u_sWallPos,transformedPos);					//saving position of south wall
	//Getting pool wall boundaries relative to south wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [-2.8,1.0,-0.1,1],modelview);		//getting position of topright of south wall
	gl.uniform4fv(u_poolSWBPos,transformedPos);					//saving position of topright of south wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [-2.8,-0.9,-0.1,1],modelview);		//getting position of bottomright of south wall
	gl.uniform4fv(u_poolSWCPos,transformedPos);					//saving position of bottomright of south wall
	scale(modelview, modelview,[6.0, 2.0, 0.2]);
	installModel(objects[3]);
	currentModelNumber = 3;
	gl.uniform4f(u_diffuseColor, 0, 0.5, 0.4, 0);
	update_uniform(modelview,projection,currentModelNumber);
	gl.uniform1i(u_isLightSource, 0);
	mvPopMatrix();
	
	//west wall
	mvPushMatrix();
	translate(modelview, modelview, [-2.9,0,1.1]);
	rotate(modelview, modelview, ((Math.PI * 90.0) / 180.0),[1,0,0]);
	rotate(modelview, modelview, ((Math.PI * 90.0) / 180.0),[0,1,0]);
	//wall stuff
	mat3.normalFromMat4(wWallMatrix, modelview);
    gl.uniformMatrix3fv(u_wWallNMatrix, false, wWallMatrix);		//saving normalmatrix for west wall
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,0.1,1],modelview);
	gl.uniform4fv(u_wWallPos,transformedPos);					//saving position of west wall
	scale(modelview, modelview,[6.0, 2.0, 0.2]);
	installModel(objects[3]);
	currentModelNumber = 3;
	gl.uniform4f(u_diffuseColor, 0, 0.5, 0.4, 1);
	update_uniform(modelview,projection,currentModelNumber);
	gl.uniform1i(u_isLightSource, 0);
	mvPopMatrix();
	
	
}

//function to draw the water "mesh"/"heightmap" which will serve as the water surface.
//114 final project new content
function drawWaterSurface(waterDrawn) {
	mvPushMatrix();
	scale(modelview, modelview,[1.0, 1.0, 1.0]);
	installModel(objects[6]);
	currentModelNumber = 6;
	gl.uniform4f(u_diffuseColor, 0.1, 0.1, 0.1, 1.0);
	gl.uniform1i(u_isWater, 1);
	update_uniform(modelview,projection,currentModelNumber);
	gl.uniform1i(u_isLightSource, 0);
	gl.uniform1i(u_isWater, 0);
	mvPopMatrix();
}


function draw() { 
	gl.clearColor(0,0,0,1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    mat4.perspective(projection,Math.PI/5,1,10,20);
    modelview = rotator.getViewMatrix();
	
	
	if (document.getElementById("animate").checked) {
		timeInt += 1;
		rotateDeg = timeInt % 360;
    }else{
		timeInt = 0;
		rotateDeg = timeInt % 360;
    }
	
	//Doing rotations so the scene is tilted at the right angles
	rotate(modelview, modelview, ((Math.PI * 0.0) / 180.0),[1,0,0]);
	translate(modelview, modelview, [0,0,-1]);
	
	daylight = 1;
	gl.uniform1i(u_isDay, 1)
	
	
	//drawing water tank that will be filled with water
	//114 final project new content
	mvPushMatrix();
	drawWaterTank();
	mvPopMatrix();
	
	
	//drawing red ball
	//114 final project new content
	mvPushMatrix();
	rotate(modelview, modelview, degToRad(rotateDeg * (-2)),[0,0,1]);
	heightDiff = Math.sin(degToRad(rotateDeg * (-10))) + Math.cos(degToRad(rotateDeg * (-10)));
	translate(modelview, modelview, [1.5,1.5,0.78 + (heightDiff * 0.1)]);
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,0,1],modelview);
	gl.uniform4fv(u_redSpherePos,transformedPos);
	scale(modelview, modelview,[0.7,0.7,0.7]);
	installModel(objects[2]);
    currentModelNumber = 2;
	gl.uniform4f(u_diffuseColor, 1, 0, 0, 1);
	gl.uniform3f(u_rSphereColor, 1.0, 0.0, 0.0);
	update_uniform(modelview,projection,currentModelNumber);
	mvPopMatrix();
	
	//drawing yellow ball
	//114 final project new content
	mvPushMatrix();
	rotate(modelview, modelview, degToRad(rotateDeg * (-2)),[0,0,1]);
	heightDiff = Math.sin(degToRad(rotateDeg * (-10))) + Math.cos(degToRad(rotateDeg * (-10)));
	translate(modelview, modelview, [-1.5,-1.5,0.78 - (heightDiff * 0.1)]);
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,0,1],modelview);
	gl.uniform4fv(u_yellowSpherePos,transformedPos);
	scale(modelview, modelview,[0.7,0.7,0.7]);
	installModel(objects[2]);
    currentModelNumber = 2;
	gl.uniform4f(u_diffuseColor, 1, 1, 0, 1);
	gl.uniform3f(u_ySphereColor, 1.0, 1.0, 0.0);
	update_uniform(modelview,projection,currentModelNumber);
	mvPopMatrix();
	
	//drawing blue ball (completely underwater)
	//114 final project new content
	mvPushMatrix();
	heightDiff = Math.sin(degToRad(rotateDeg * (-10)));
	translate(modelview, modelview, [0.0,0.0,0.5 - (heightDiff * 0.3)]);
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,0,1],modelview);
	gl.uniform4fv(u_blueSpherePos,transformedPos);
	scale(modelview, modelview,[0.3,0.3,0.3]);
	installModel(objects[2]);
    currentModelNumber = 2;
	gl.uniform4f(u_diffuseColor, 0, 0, 1, 1);
	gl.uniform3f(u_bSphereColor, 0.0, 0.0, 1.0);
	update_uniform(modelview,projection,currentModelNumber);
	mvPopMatrix();
	
	
	//drawing orange ball (completely overwater, starts on top of the far northernwall, but will swing around both to the left and the right periodically)
	//114 final project new content
	mvPushMatrix();
	xDiff = Math.sin(degToRad(rotateDeg * (2.0)));
	yDiff = Math.cos(degToRad(rotateDeg * (2.0)));
	yDiff = Math.abs(yDiff);
	translate(modelview, modelview, [(0.0 + (xDiff * 3.0)),(3.0 * yDiff),3.35]);
	var transformedPos = new Float32Array(4);
	vec4.transformMat4(transformedPos, [0,0,0,1],modelview);
	gl.uniform4fv(u_orangeSpherePos,transformedPos);
	scale(modelview, modelview,[1.2,1.2,1.2]);
	installModel(objects[2]);
    currentModelNumber = 2;
	gl.uniform4f(u_diffuseColor, 0.8, 0.4, 0, 1);
	gl.uniform3f(u_oSphereColor, 0.8, 0.4, 0.0);
	update_uniform(modelview,projection,currentModelNumber);
	mvPopMatrix();
	
	
	//drawing water surface 
	//114 final project new content
	
	mvPushMatrix();
	translate(modelview, modelview, [-3.0,-3.0,1.2]);
	scale(modelview, modelview,[2.9,2.9,1.0]);
	drawWaterSurface();
	mvPopMatrix();
	
	mvPushMatrix();		//Drawing specular light (mostly based on sun implementation of the 112 final project)
	translate(modelview, modelview, [0,0.0,4.0]);
	if (daylight == 1) {
		var transformedPos = new Float32Array(4);
		vec4.transformMat4(transformedPos, [0,0,0,1],modelview);
		gl.uniform4fv(u_sunPosition,transformedPos);
		gl.uniform3f(u_specularColor, 1, 1, 1);
		gl.uniform4f(u_diffuseColor, 0.1, 0.1, 0.1, 0.1);
		gl.uniform1f(u_specularExponent, 10);
	}
	mvPopMatrix();
	
	/*
	//Attempts to make a function that draws an ever-moving water surface
	//waterMove(timeStep);
	
	if (document.getElementById("animate").checked) {
		waterMove(timeStep);
    }
	*/
	
}


/* 
 * Called and data for the model are copied into the appropriate buffers, and the 
 * scene is drawn.
 */
function installModel(modelData) {
     gl.bindBuffer(gl.ARRAY_BUFFER, a_coords_buffer);
     gl.bufferData(gl.ARRAY_BUFFER, modelData.vertexPositions, gl.STATIC_DRAW);
     gl.vertexAttribPointer(a_coords_loc, 3, gl.FLOAT, false, 0, 0);
     gl.enableVertexAttribArray(a_coords_loc);
     gl.bindBuffer(gl.ARRAY_BUFFER, a_normal_buffer);
     gl.bufferData(gl.ARRAY_BUFFER, modelData.vertexNormals, gl.STATIC_DRAW);
     gl.vertexAttribPointer(a_normal_loc, 3, gl.FLOAT, false, 0, 0);
     gl.enableVertexAttribArray(a_normal_loc);
     gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,index_buffer);
     gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, modelData.indices, gl.STATIC_DRAW);
}


/* Initialize the WebGL context.  Called from init() 
	NOTE: All uncommented lines are implementations that were part of the 112 project. The lines added for the 114 Final Project are commented on.
*/
function initGL() {
    var prog = createProgram(gl,"vshader-source","fshader-source");
    gl.useProgram(prog);
    a_coords_loc =  gl.getAttribLocation(prog, "a_coords");
    a_normal_loc =  gl.getAttribLocation(prog, "a_normal");
    u_modelview = gl.getUniformLocation(prog, "modelview");
    u_projection = gl.getUniformLocation(prog, "projection");
    u_normalMatrix =  gl.getUniformLocation(prog, "normalMatrix");
	
	u_pFloorNMatrix =  gl.getUniformLocation(prog, "pFloorNMatrix");	//for saving the normal of the floor
	u_nWallNMatrix =  gl.getUniformLocation(prog, "nWallNMatrix");		//for saving the normal of the north wall
	u_sWallNMatrix =  gl.getUniformLocation(prog, "sWallNMatrix");		//for saving the normal of the south wall
	u_eWallNMatrix =  gl.getUniformLocation(prog, "eWallNMatrix");		//for saving the normal of the east wall
	u_wWallNMatrix =  gl.getUniformLocation(prog, "wWallNMatrix");		//for saving the normal of the west wall
	
    u_diffuseColor =  gl.getUniformLocation(prog, "diffuseColor");
    u_specularColor =  gl.getUniformLocation(prog, "specularColor");
    u_specularExponent = gl.getUniformLocation(prog, "specularExponent");	
	u_rSphereColor =  gl.getUniformLocation(prog, "rSphereColor");		//for saving the color of the red sphere
    u_ySphereColor = gl.getUniformLocation(prog, "ySphereColor");		//for saving the color of the yellow sphere
	u_bSphereColor = gl.getUniformLocation(prog, "bSphereColor");		//for saving the color of the blue sphere
	u_oSphereColor = gl.getUniformLocation(prog, "oSphereColor");		//for saving the color of the orange sphere
	u_wallColor = gl.getUniformLocation(prog, "wallColor");			//for saving the color of the walls
	
	u_isDay = gl.getUniformLocation(prog,"isDay");	 //
	u_isLightSource = gl.getUniformLocation(prog,"isLightSource");
	u_isWater = gl.getUniformLocation(prog,"isWater");				//isWater flag
	u_isUnderwater = gl.getUniformLocation(prog,"isUnderwater");	//isUnderwater flag
	
	u_sunPosition = gl.getUniformLocation(prog, "sunPosition");		//using the "sunPosition" as the position of the overall scene light
	
	u_redSpherePos = gl.getUniformLocation(prog, "redSpherePos");	//for the position of the red sphere
	u_yellowSpherePos = gl.getUniformLocation(prog, "yellowSpherePos");	//for the position of the yellow sphere
	u_blueSpherePos = gl.getUniformLocation(prog, "blueSpherePos");	//for the position of the blue sphere
	u_orangeSpherePos = gl.getUniformLocation(prog, "orangeSpherePos");	//for the position of the orange sphere
	u_poolFloorPos = gl.getUniformLocation(prog, "poolFloorPos");	//for the position of the floor (for plane calculations)
	u_nWallPos = gl.getUniformLocation(prog, "nWallPos");	//for the position of the north wall (for plane calculations)
	u_sWallPos = gl.getUniformLocation(prog, "sWallPos");	//for the position of the south wall (for plane calculations)
	u_eWallPos = gl.getUniformLocation(prog, "eWallPos");	//for the position of the east wall (for plane calculations)
	u_wWallPos = gl.getUniformLocation(prog, "wWallPos");	//for the position of the west wall (for plane calculations)
	
	u_poolNWAPos = gl.getUniformLocation(prog, "poolNWAPos");;	//saving position of far top left corner of pool wall
	u_poolNWBPos = gl.getUniformLocation(prog, "poolNWBPos");;	//saving position of far top right corner of pool wall
	u_poolNWCPos = gl.getUniformLocation(prog, "poolNWCPos");;	//saving position of far bottom right corner of pool wall
	u_poolNWDPos = gl.getUniformLocation(prog, "poolNWDPos");;	//saving position of far bottom left corner of pool wall
	u_poolEWBPos = gl.getUniformLocation(prog, "poolEWBPos");;	//saving position of near top right corner of pool wall
	u_poolEWCPos = gl.getUniformLocation(prog, "poolEWCPos");;	//saving position of near bottom right corner of pool wall
	u_poolSWBPos = gl.getUniformLocation(prog, "poolSWBPos");;	//saving position of near top left corner of pool wall
	u_poolSWCPos = gl.getUniformLocation(prog, "poolSWCPos");;	//saving position of near bottom left corner of pool wall
	
	
    a_coords_buffer = gl.createBuffer();
    a_normal_buffer = gl.createBuffer();
    index_buffer = gl.createBuffer();
    gl.enable(gl.DEPTH_TEST);
	
    gl.uniform3f(u_specularColor, 1.0, 1.0, 1.0);   
    gl.uniform4f(u_diffuseColor, 0, 0, 0, 0);       
    gl.uniform1f(u_specularExponent, 10);
	gl.uniform3f(u_wallColor, 0.0, 0.5, 0.4);
}

/* Creates a program for use in the WebGL context gl, and returns the
 * identifier for that program.  If an error occurs while compiling or
 * linking the program, an exception of type String is thrown.  The error
 * string contains the compilation or linking error.  If no error occurs,
 * the program identifier is the return value of the function.
 *    The second and third parameters are the id attributes for <script>
 * elementst that contain the source code for the vertex and fragment
 * shaders.
 */
function createProgram(gl, vertexShaderID, fragmentShaderID) {
    function getTextContent( elementID ) {
            // This nested function retrieves the text content of an
            // element on the web page.  It is used here to get the shader
            // source code from the script elements that contain it.
        var element = document.getElementById(elementID);
        var node = element.firstChild;
        var str = "";
        while (node) {
            if (node.nodeType == 3) // this is a text node
                str += node.textContent;
            node = node.nextSibling;
        }
        return str;
    }
    try {
        var vertexShaderSource = getTextContent( vertexShaderID );
        var fragmentShaderSource = getTextContent( fragmentShaderID );
    }
    catch (e) {
        throw "Error: Could not get shader source code from script elements.";
    }
    var vsh = gl.createShader( gl.VERTEX_SHADER );
    gl.shaderSource(vsh,vertexShaderSource);
    gl.compileShader(vsh);
    if ( ! gl.getShaderParameter(vsh, gl.COMPILE_STATUS) ) {
        throw "Error in vertex shader:  " + gl.getShaderInfoLog(vsh);
     }
    var fsh = gl.createShader( gl.FRAGMENT_SHADER );
    gl.shaderSource(fsh, fragmentShaderSource);
    gl.compileShader(fsh);
    if ( ! gl.getShaderParameter(fsh, gl.COMPILE_STATUS) ) {
       throw "Error in fragment shader:  " + gl.getShaderInfoLog(fsh);
    }
    var prog = gl.createProgram();
    gl.attachShader(prog,vsh);
    gl.attachShader(prog, fsh);
    gl.linkProgram(prog);
    if ( ! gl.getProgramParameter( prog, gl.LINK_STATUS) ) {
       throw "Link error in program:  " + gl.getProgramInfoLog(prog);
    }
    return prog;
}


/**
 * initialization function that will be called when the page has loaded
 */
function init() {
    try {
        var canvas = document.getElementById("myGLCanvas");
        gl = canvas.getContext("webgl") || 
                         canvas.getContext("experimental-webgl");
        if ( ! gl ) {
            throw "Browser does not support WebGL";
        }
    }
    catch (e) {
        document.getElementById("canvas-holder").innerHTML =
            "<p>Sorry, could not get a WebGL graphics context.</p>";
        return;
    }
    try {
        initGL();  // initialize the WebGL graphics context
    }
    catch (e) {
        document.getElementById("canvas-holder").innerHTML =
            "<p>Sorry, could not initialize the WebGL graphics context:" + e + "</p>";
        return;
    }

    document.getElementById("animate").checked = false;
    rotator = new TrackballRotator(canvas, draw, 15);
    tick();
}

function tick() {
	requestAnimFrame(tick);
    draw();
	
	
}



