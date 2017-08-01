/***********************************************************************************************************
 *                                                                                                         *
 *         *         *         *         *     Rendering functions     *        *        *        *        *
 *                                                                                                         *
 ***********************************************************************************************************/


/**
 * @param colorArray: the array of vertices of the mesh to be constructed
 * 
 * Works the same way as the loadOverlayMesh function
 */
function constructOverlayMesh(colorArray) {
	clearOverlayMesh();
	overlayMesh = {};
	
	// Add to the mesh vertices black and white color points
	var increasedColorArray = [];
	
	for (var i = 0; i < colorArray.length; i++) {
		increasedColorArray.push(colorArray[i]);
	}
	increasedColorArray.push(whiteColor);
	increasedColorArray.push(blackColor);
	
	var geometry = new THREE.Geometry();
	
	/// "Dumb" faces which are each a triplet of positions.
	// Make them smart (faces should be indices into hull.vs).
	var smart_faces = [];
	// Iterate over the faces.
	for( var fi = 0; fi < colorArray.length; ++fi ) { // number of faces equal to the number reference colors
		var dumb_face = [blackColor, whiteColor, colorArray[fi]];
		smart_faces.push( [] );
		// Each vertex of the face is a 3D position. Find the index of this position in hull.vs;
		for( var fvi = 0; fvi < 3; ++fvi ) {
			var vert = dumb_face[fvi];
			
			// Find the index of vert in hull.vs.
			var hull_vertex_index = -1;
			for( var vi = 0; vi < increasedColorArray.length; ++vi ) {
				if( vert[0] == increasedColorArray[vi][0] && vert[1] == increasedColorArray[vi][1] && vert[2] == increasedColorArray[vi][2] ) {
					hull_vertex_index = vi;
					break;
				}
			}
			if( -1 === hull_vertex_index ) { console.error("Overlay mesh faces use a non-existent vertex.") }
			smart_faces[fi].push( hull_vertex_index );
		}
	}
	
	
	/// "Smart" faces which are indices into vertices.
	// smart_faces = hull.faces;
	// Copy vertices
	for( var i = 0; i < increasedColorArray.length; ++i ) {
		var vert = increasedColorArray[i];
		geometry.vertices.push( new THREE.Vector3( vert[0]/255., vert[1]/255., vert[2]/255. ) );
	}
	// Copy faces
	for( var i = 0; i < colorArray.length; ++i ) {
		var face = smart_faces[i]; // hull.faces[i];
		geometry.faces.push( new THREE.Face3( face[0], face[1], face[2] ) );
	}
	
	geometry.verticesNeedUpdate = true;
	geometry.elementsNeedUpdate = true;
	
	var edge_material = new THREE.MeshBasicMaterial({
		color: 0xFFFFFF,
		wireframe: true,
		wireframeLinewidth: 2
	});
	
	overlayMesh.edges = new THREE.Mesh( geometry, edge_material );
	scene.add( overlayMesh.edges );
	
	overlayMesh.vertices = new THREE.Object3D();
	overlayMesh.verticesRims = new THREE.Object3D();
	// We will keep a copy of the original position in case we later move them.
	overlayMesh.originalPositions = [];
	
	for (var i = 0; i < increasedColorArray.length; i++) {
		var color = increasedColorArray[i];
		var col = new THREE.Vector3(color[0]/255., color[1]/255., color[2]/255.);
		col.clampScalar(0., 1.);
		var circleMaterial = new THREE.MeshBasicMaterial({
			color: new THREE.Color( col.x, col.y, col.z )
		});
		var radius = 5./255.;
		var segments = 32;
		var circleGeometry = new THREE.SphereGeometry( radius, segments, segments );
		var circle = new THREE.Mesh( circleGeometry, circleMaterial );
		circle.position.set( color[0]/255., color[1]/255., color[2]/255. );
		// Make a copy of the position to have in case we later move it.
		overlayMesh.originalPositions.push( circle.position.clone() );
		
		overlayMesh.vertices.add( circle );
		
		
		// HACK: A white rim.
		circleGeometry = new THREE.CircleGeometry( radius*1.2, segments );
		circle = new THREE.Mesh( circleGeometry, new THREE.MeshBasicMaterial({ color: 0xFFFFFF }) );
		circle.position.set( color[0]/255., color[1]/255., color[2]/255. );
		
		overlayMesh.verticesRims.add( circle );
		
		scene.add( overlayMesh.vertices );
		scene.add( overlayMesh.verticesRims );
	}
	
	render();
}

/**
 * From original index.html
 */
function createParticles( params ) {
	if( params === undefined ) params = {};
	// By default, only create particles for unique colors.
	if( !( 'only_unique_pixels' in params ) ) params.only_unique_pixels = true;
	// By default, clear the overlay mesh.
	if( !( 'clear_overlay_mesh' in params ) ) params.clear_overlay_mesh = true;

	if( particleSystem !== undefined ) scene.remove( particleSystem );

	if (kmeans_particles !== undefined) {
		for (var i = 0; i < kmeans_particles.length; i++) {
			var cloud = kmeans_particles[i];
			scene.remove(cloud);
		}
	}

	// Also clear the overlay mesh if there is one.
	if( params.clear_overlay_mesh ) clearOverlayMesh();

	// Get the particles from the image's pixels.
	// From: http://stackoverflow.com/questions/1041399/how-to-use-javascript-or-jquery-to-read-a-pixel-of-an-image
	var canvas = document.createElement('canvas');
	canvas.width = img.naturalWidth;
	canvas.height = img.naturalHeight;
	var ctx = canvas.getContext("2d");
	ctx.drawImage( img, 0, 0, img.naturalWidth, img.naturalHeight );
	// getImageData returns an RGBA byte array.
	var pixels = ctx.getImageData( 0, 0, canvas.width, canvas.height ).data;
	
	var particles = pixels.length / 4;
	RGBPointsArray = [];
	
	for ( var i = 0; i < particles; i += 1 ) {
		// colors are also positions
		var r = pixels[ 4*i + 0 ]/255.;
		var g = pixels[ 4*i + 1 ]/255.;
		var b = pixels[ 4*i + 2 ]/255.;
		
		RGBPointsArray.push([r, g, b]);
	}
	
	if( params.only_unique_pixels ) {
		// Only take unique colors.
		pixels = uniquePixels( pixels );
	}

	particles = pixels.length / 4;

	var geometry = new THREE.BufferGeometry();

	var positions = new Float32Array( particles * 3 );
	var colors = new Float32Array( particles * 3 );

	var color = new THREE.Color();
	RGBUniquePointsArray = [];

	for ( var i = 0; i < particles; i += 1 ) {

		// colors are also positions
		var r = pixels[ 4*i + 0 ]/255.;
		var g = pixels[ 4*i + 1 ]/255.;
		var b = pixels[ 4*i + 2 ]/255.;
		
		RGBUniquePointsArray.push([r, g, b]);
		
		positions[ 3*i + 0 ] = r;
		positions[ 3*i + 1 ] = g;
		positions[ 3*i + 2 ] = b;

		colors[ 3*i + 0 ] = r;
		colors[ 3*i + 1 ] = g;
		colors[ 3*i + 2 ] = b;

	}

	document.getElementById('debug').innerHTML =
		'width: ' + img.naturalWidth + ', height: ' + img.naturalHeight
		+
		'<br>total pixels: ' + canvas.width*canvas.height
		+
		'<br>unique pixels: ' + particles
		;

	geometry.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
	geometry.addAttribute( 'color', new THREE.BufferAttribute( colors, 3 ) );

	geometry.computeBoundingSphere();

	//

	var material = new THREE.PointCloudMaterial( { size: 4/255., vertexColors: THREE.VertexColors } );

	particleSystem = new THREE.PointCloud( geometry, material );
	// return particleSystem;
	scene.add( particleSystem );

	render();
}

/**
 * @param colorArray: the RGB data array to be represented as a point cloud in an RGB space
 * @param globalColor: the color assigned to every point of the data set
 */
function createKmeansParticles(colorArray, globalColor)
{
	var particles = colorArray.length;

	var geometry = new THREE.BufferGeometry();

	var positions = new Float32Array( particles * 3 );
	var colors = new Float32Array( particles * 3 );

	var color = new THREE.Color();

	for ( var i = 0; i < particles; i += 1 ) {
		
		// colors are also positions
		var r = colorArray[i][0];
		var g = colorArray[i][1];
		var b = colorArray[i][2];
		
		positions[ 3*i + 0 ] = r;
		positions[ 3*i + 1 ] = g;
		positions[ 3*i + 2 ] = b;
		
		if (globalColor === null) {
			colors[ 3*i + 0 ] = r;
			colors[ 3*i + 1 ] = g;
			colors[ 3*i + 2 ] = b;
		} else { // We want the color equal to the position (see restorePointsColor)
			colors[ 3*i + 0 ] = globalColor[0]/255.;
			colors[ 3*i + 1 ] = globalColor[1]/255.;
			colors[ 3*i + 2 ] = globalColor[2]/255.;
		}

	}

	geometry.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
	geometry.addAttribute( 'color', new THREE.BufferAttribute( colors, 3 ) );

	geometry.computeBoundingSphere();

	//

	var material = new THREE.PointCloudMaterial( { size: 4/255., vertexColors: THREE.VertexColors } );

	var particleCloud = new THREE.PointCloud( geometry, material );
	kmeans_particles.push(particleCloud);
	// return particleSystem;
	scene.add( particleCloud );
}

/**
 * @param index: the index of the text field where the color to save as a reference color is specified
 */
function saveColor(index) {
	var components;
	var color;
	
	if (index == 0) {
		components = ['r_comp0', 'g_comp0', 'b_comp0'];
	}
	else if (index == 1) {
		components = ['r_comp1', 'g_comp1', 'b_comp1'];
	}
	else {
		components = ['r_comp2', 'g_comp2', 'b_comp2'];
	}
	
	color = [
		document.getElementById(components[0]).value,
		document.getElementById(components[1]).value,
		document.getElementById(components[2]).value
	];
	
	if (!(color[0] >= 0 && color[0] < 256 && color[1] >= 0 && color[1] < 256 && color[2] >= 0 && color[2] < 256) || color[0] == '' || color[1] == '' || color[2] == '') {
		console.warn('[WARNING] Some color components are incorrect');
	} else {
		tabRefColors.push(color);
	}
}

/*
 * Draw mesh defined by reference colors specified by user
 */
function drawReferenceColors() {
	tabRefColors = [];
	
	saveColor(0);
	saveColor(1);
	saveColor(2);
	
	constructOverlayMesh(tabRefColors);
}

/*
 * Draw kmeans-computed centroids from initial specified color points, and replace them in RGB space
 */
function drawCentroids() {
	if (tabRefColors.length < 1) {
		console.warn('[WARNING] Nonexistent reference colors. Please specify first the reference colors to build the corresponding mesh');
	}
	
	var kmeansResult = kmeans(RGBUniquePointsArray, tabRefColors.length, tabRefColors);
	tabRefColors = kmeansResult.centroids;
	tabClusters = kmeansResult.clusters;
		
	
	constructOverlayMesh(tabRefColors);
}

function drawTriangleVertices() {
	if (tabRefColors.length < 1) {
		console.warn('[WARNING] Nonexistent reference colors. Please specify first the reference colors to build the corresponding mesh');
	}
	
	var ditrimeansResult = ditrimeans(RGBUniquePointsArray, tabRefColors.length, tabRefColors);
	tabRefColors = ditrimeansResult.centroids;
	
// 	console.log(tabRefColors[0]);
// 	console.log(tabRefColors[1]);
	
	tabClusters = ditrimeansResult.clusters;
	
	constructOverlayMesh(tabRefColors);
	
	recolor_data = {};
	saveWeights();
}

function colorizeClusters() {
	if (tabClusters.length < 1) {
		console.warn('[WARNING] No Kmeans centroids drawn. Please compute Kmeans centroids to determine clusters to colorize.');
	}
	
	if( kmeans_particles !== undefined ) {
		for (var i = 0; i < kmeans_particles.length; i++) {
			scene.remove(kmeans_particles[i]);
		}
	}
	
	scene.remove(particleSystem);
	
	for (var i = 0; i < tabClusters.length; i++) {
// 				var projectedPoints = findTriangleThirdVertex(blackColor, whiteColor, tabRefColors[i], tabClusters[i]);
		createKmeansParticles(tabClusters[i], tabRefColors[i]);
	}
	
	render();
}

function restorePointsColor() {
	if( kmeans_particles !== undefined ) {
		for (var i = 0; i < kmeans_particles.length; i++) {
			scene.remove(kmeans_particles[i]);
		}
	}
	
	scene.remove(particleSystem);
	createKmeansParticles(RGBUniquePointsArray, null);
}

//

function animate() {
	requestAnimationFrame( animate );
	controls.update();
	// If we have dynamic movement, we need to call render every frame.
	if( !controls.staticMoving ) render();
}

function render() {
	// HACK: Update overlay mesh vertices' white rims.
	if( overlayMesh !== null && overlayMesh.verticesRims !== null ) {
		for( var i = 0; i < overlayMesh.verticesRims.children.length; ++i )
		{
			overlayMesh.verticesRims.children[i].lookAt( camera.position );
		}
	}
	
	// Always set the plane's orientation to look at the camera.
	// picking_data.plane.lookAt( camera.position );
	// UPDATE: That is not necessarily the look vector.
	picking_data.plane.quaternion.setFromUnitVectors(new THREE.Vector3( 0, 0, -1 ), getLookVector());
	
	renderer.render( scene, camera );
	stats.update();
}
