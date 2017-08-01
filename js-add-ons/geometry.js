/***********************************************************************************************************
 *                                                                                                         *
 *         *         *         *        *  Generic geometric functions  *       *        *        *        *
 *                                                                                                         *
 ***********************************************************************************************************/


function clone(obj){
	
	try{
		var copy = JSON.parse(JSON.stringify(obj));
	} catch(ex){
		console.error("[ERROR] Incorrect object type, or your browser does not support JSON format");
	}
	return copy;
}

/**
 * @return the distance between two points
 * @param p0: the first point as an array
 * @param p1: the second point as an array
 */
function distance(p0, p1) {
	var vec0 = new THREE.Vector3().fromArray(p0);
	var vec1 = new THREE.Vector3().fromArray(p1);
	
	return vec0.distanceTo(vec1);
}

/**
 * @return a THREE.Vector3 object corresponding to a cloned version the array given in parameter
 * @param obj: the original array vector to be cloned
 *
 * From http://www.finalclap.com/faq/371-javascript-clone-dupliquer-objet
 */
function cloneToVector(obj){
	
	try{
		var copy = JSON.parse(JSON.stringify(obj));
	} catch(ex){
		console.error("[ERROR] Incorrect object type, or your browser does not support JSON format");
	}
	return new THREE.Vector3(copy[0], copy[1], copy[2]);
}

/**
 * @return the distance between a and the segment a-b
 * @param v: the point from which the distance to the segment ab has to be computed
 * 
 * From https://stackoverflow.com/questions/4858264/find-the-distance-from-a-3d-point-to-a-line-segment
 */
function distanceToSegment(v, a, b)
{
	var ab = cloneToVector(b.toArray());
	ab.sub(a);
	var av = cloneToVector(v.toArray());
	av.sub(a);

	if (av.dot(ab) <= 0.0) {           // Point is lagging behind start of the segment, so perpendicular distance is not viable.
		return av.length();         // Use distance to start of segment instead.
	}
	
	var bv = cloneToVector(v.toArray());
	bv.sub(b);

	if (bv.dot(ab) >= 0.0) {           // Point is advanced past the end of the segment, so perpendicular distance is not viable.
		return bv.length();         // Use distance to end of the segment instead.
	}
	
	return (ab.cross(av)).length();       // Perpendicular distance of point to segment.
}

/**
 * @return the orthogonal vector to two other ones
 * @param v0: the first reference vector, as a THREE.Vector3
 * @param v1: the second reference vector, as a THREE.Vector3
 */
function computeNormal(v0, v1) {
	return cloneToVector(v0.toArray()).cross(v1).normalize();
}

/**
 * @return the barycentric position of a point on a triangle
 * @param p: the point on the triangle as an array
 * @param A.B.C: the triangle vertices as arrays
 */
function solveBarycentricEquations(A, B, C, p) {
	var vecA = new THREE.Vector3(A[0], A[1], A[2]);
	var vecB = new THREE.Vector3(B[0], B[1], B[2]);
	var vecC = new THREE.Vector3(C[0], C[1], C[2]);
	var vecP = new THREE.Vector3(p[0], p[1], p[2]);
	
	var vecDist, dist = Number.MAX_VALUE, t;
	
	
	var
		ba = cloneToVector(B),
		ca = cloneToVector(C),
		bp = cloneToVector(B),
		cp = cloneToVector(C),
		ap = cloneToVector(A)
	;
	
	ba.sub(vecA);
	ca.sub(vecA);
	bp.sub(vecP);
	cp.sub(vecP);
	ap.sub(vecP);
	
	var normal = computeNormal(ba, ca);
	
	var
		normalABC = cloneToVector(normal.toArray()),
		normalPBC = cloneToVector(normal.toArray()),
		normalPCA = cloneToVector(normal.toArray())
	;
	
	var areaABC = normalABC.dot(ba.cross(ca));
	var areaPBC = normalPBC.dot(bp.cross(cp));
	var areaPCA = normalPCA.dot(cp.cross(ap));
	
	var alpha = areaPBC / areaABC;
	var beta = areaPCA / areaABC;
	var gamma = 1.0 - alpha - beta;
	
	return {
		normal: normal,
		alphaPCB: alpha,
		betaPCA: beta,
		gammaPBA: gamma
	};
}

/**
 * @return the distance between a point and a triangle
 * @param p: the point whose distance to the plane has to be computed, as an array
 * @param A.B.C: the triangle vertices as arrays
 * 
 * From https://math.stackexchange.com/questions/588871/minimum-distance-between-point-and-face
 */
function pointTriangleDistance(A, B, C, p) {
	
	var vecA = new THREE.Vector3(A[0], A[1], A[2]);
	var vecB = new THREE.Vector3(B[0], B[1], B[2]);
	var vecC = new THREE.Vector3(C[0], C[1], C[2]);
	var vecP = new THREE.Vector3(p[0], p[1], p[2]);
	
	var result = solveBarycentricEquations(A, B, C, p);
	
	var alpha = result.alphaPCB;
	var beta = result.betaPCA;
	var gamma = result.gammaPBA;
	var normal = result.normal; // THREE.Vector3
	
	if (alpha < 0. || alpha > 1. || beta < 0. || beta > 1. || gamma < 0. || gamma > 1.) {
		if (alpha >= 0. && beta < 0. && gamma < 0.) {
			dist = vecP.distanceTo(cloneToVector(vecA.toArray()));
		} else if (beta >= 0. && alpha < 0. && gamma < 0.) {
			dist = vecP.distanceTo(cloneToVector(vecB.toArray()));
		} else if (gamma >= 0. && alpha < 0. && beta < 0.) {
			dist = vecP.distanceTo(cloneToVector(vecC.toArray()));
		} else if (gamma < 0. && alpha >= 0. && beta >= 0.) {
			dist = distanceToSegment(cloneToVector(vecP.toArray()), cloneToVector(vecA.toArray()), cloneToVector(vecB.toArray()));
		} else if (beta < 0. && alpha >= 0. && gamma >= 0.) {
			dist = distanceToSegment(cloneToVector(vecP.toArray()), cloneToVector(vecA.toArray()), cloneToVector(vecC.toArray()));
		} else if (alpha < 0. && beta >= 0. && gamma >= 0.) {
			dist = distanceToSegment(cloneToVector(vecP.toArray()), cloneToVector(vecB.toArray()), cloneToVector(vecC.toArray()));
		}
	} else {
		dist_aux = cloneToVector(normal.toArray());
		vecDist = cloneToVector(dist_aux.toArray());
		t = vecDist.dot(vecA) - normal.dot(vecP);
		var p0 = cloneToVector(p);
		p0.add(normal.multiplyScalar(t));
		dist = p0.distanceTo(vecP);
	}
	
	return dist;
}

/**
 * @return the orthogonal distance of a point to the plane defined by a triangle
 * @param p: the point whose distance to the plane has to be computed, as an array
 * @param A.B.C: the triangle vertices as arrays
 */
function pointPlaneDistance(A, B, C, p) {
	
	var vecA = new THREE.Vector3(A[0], A[1], A[2]);
	var vecB = new THREE.Vector3(B[0], B[1], B[2]);
	var vecC = new THREE.Vector3(C[0], C[1], C[2]);
	var vecP = new THREE.Vector3(p[0], p[1], p[2]);
	
	var ba = cloneToVector(B);
	var ca = cloneToVector(C);
	ba.sub(vecA);
	ca.sub(vecA);
	
	var normal = computeNormal(ba, ca); // THREE.Vector3
	
	dist_aux = cloneToVector(normal.toArray());
	vecDist = cloneToVector(dist_aux.toArray());
	t = vecDist.dot(vecA) - normal.dot(vecP);
	var p0 = cloneToVector(p);
	p0.add(normal.multiplyScalar(t));
	dist = p0.distanceTo(vecP);
	
	return dist;
}

/* Basic K-means algorithm by Jonathan Spicer
* URL: http://www.mymessedupmind.co.uk/index.php/javascript-k-mean-algorithm/
*/
function kmeans( arrayToProcess, Clusters, initialCentroids ) { // arrayToProcess is an array of 3D points

	var Groups = new Array();
	var Centroids = initialCentroids;
	var oldCentroids = new Array();
	var changed = false;

	var count = 0;

	do {
		for( j=0; j < Clusters; j++ ) {
			Groups[j] = [];
		}

		changed=false;

		/* Clusterization */
		for( i=0; i < arrayToProcess.length; i++ ) {

			Distance=-1;
			oldDistance=-1

			for( j=0; j < Clusters; j++ ) {
				var vecCentroid = new THREE.Vector3(Centroids[j][0], Centroids[j][1], Centroids[j][2]);
				var vecParticle = new THREE.Vector3(arrayToProcess[i][0], arrayToProcess[i][1], arrayToProcess[i][2]);
				
				distance = vecCentroid.distanceTo(vecParticle);
				var newGroup;

				if ( oldDistance==-1 ) {
					oldDistance = distance;
					newGroup = j;
				} else if ( distance <= oldDistance ) {
					newGroup=j;
					oldDistance = distance;
				}

			}

			Groups[newGroup].push( arrayToProcess[i] );  
		}

		oldCentroids=Centroids.slice(0);

		/* New centroids computation */
		for ( j=0; j < Clusters; j++ ) {
			total = new THREE.Vector3(0, 0, 0);

			for( i=0; i < Groups[j].length; i++ ) {
				var arrayVector = Groups[j][i];
				var v = new THREE.Vector3().fromArray(arrayVector);
				total.add(v);
			}
			var newCentroid = total.divideScalar(Groups[j].length);

			Centroids[j] = newCentroid.toArray();
		}

		/* Compare new centroids to the older ones */
		for( j=0; j < Clusters; j++ ) {
			vecCentroid.fromArray(Centroids[j]);
			var vecOldCentroid = new THREE.Vector3(oldCentroids[j][0], oldCentroids[j][1], oldCentroids[j][2]);
			
			if ( vecCentroid.distanceTo(vecOldCentroid) > 1e-4 ) {
				changed = true;
			}
		}

	} while(changed);

	return {
		centroids: Centroids.map(function(x) x.map(function(y) y * 255)),
		clusters: Groups
	};
}

/**
 * @return an object which contains the final cluster centers positions and the point clusters themselves
 * @param arrayToProcess: the data array on which the k-means algorithm has to be applied
 * @param Clusters: the number of clusters we would like to finally obtain
 * @param initialCentroids: the set of centroids initial positions
 * 
 * Adaptation of the Jonathan Spicer's K-means algorithm defined above
 * Heuristic has changed here : we take into account the 
 */
function ditrimeans( arrayToProcess, Clusters, initialCentroids ) { // arrayToProcess is an array of 3D points

	var Groups = new Array();
	var Centroids = initialCentroids;
	var oldCentroids = new Array();
	var changed = false;

	do {
		for( j=0; j < Clusters; j++ ) {
			Groups[j] = [];
		}
		
		changed = false;
		
		// Clusterization
		for( i=0; i < arrayToProcess.length; i++ ) {

			/* Assignation of each point of the array to one of the clusters, depending on its distance from the triangle formed with the centroid of the current cluster */
			
			Distance=-1;
			oldDistance=-1

			for( j=0; j < Clusters; j++ ) {
				var centroid = [parseInt(Centroids[j][0])/255., parseInt(Centroids[j][1])/255., parseInt(Centroids[j][2])/255.];
				var particle = arrayToProcess[i];
				
				distance = pointTriangleDistance(blackColor, whiteColor, centroid, particle);
				
				var newGroup;

				if ( oldDistance==-1 ) {
					oldDistance = distance;
					newGroup = j;
				} else if ( distance <= oldDistance ) {
					newGroup=j;
					oldDistance = distance;
				}

			}

			Groups[newGroup].push( arrayToProcess[i] );
		}
		
		oldCentroids=Centroids.slice(0);
		
		/* New centroids determination, computing the mean orientation vector */
		for ( j=0; j < Clusters; j++ ) {
			total = new THREE.Vector3(0, 0, 0);
			var refDist = distanceToSegment(cloneToVector(Centroids[j]), cloneToVector(blackColor), cloneToVector([1., 1., 1.]));

			/* For each point of the cluster, compute the vector which passes through it and orthogonal to the chromatic axis */
			
			for( i=0; i < Groups[j].length; i++ ) {
				var arrayPoint = Groups[j][i];
				var v = new THREE.Vector3().fromArray(arrayPoint);
				v.normalize();
				total.add(v);
			}
			
			/* Get the mean orientation unit vector */
			
			var newCentroid = total.divideScalar(Groups[j].length);
			
			/* Determine the new centroid position along the computed vector at the original refDist distance */
			
			newCentroid.multiplyScalar(refDist);
			Centroids[j] = newCentroid.toArray();
		}

		/* Compare new centroids to the older ones */
		for( j=0; j < Clusters; j++ ) {
			var vecCentroid = new THREE.Vector3().fromArray(Centroids[j]);
			var vecOldCentroid = new THREE.Vector3().fromArray(oldCentroids);
			
			if ( vecCentroid.distanceTo(vecOldCentroid) > 1e-4 ) {
				changed = true;
			}
		}

		for (var k = 0; k < Clusters; k++) {
			var white = [1., 1., 1.];
			var result = findTriangleThirdVertex(blackColor, white, [parseInt(Centroids[k][0])/255., parseInt(Centroids[k][1])/255., parseInt(Centroids[k][2])/255.], Groups[k]);
			var recomputedCentroid = result.thirdPoint;
			Centroids[k] = [recomputedCentroid.x * 255., recomputedCentroid.y * 255., recomputedCentroid.z * 255.];
			//Groups[k] = result.projectedPoints;
		}
	} while(changed);

	return { centroids: Centroids, clusters: Groups };
}

/**
 * @param p: the point which has to be projected on the plane defined by a normal and an orientation vector, as an array
 * @param v1: one of the two orientation vectors of the plane, as an array
 * @return the new projected point p
 */
function projectPointOnPlane(p, normal, v1) {
	var p0 = new THREE.Vector3().fromArray(p);
	var p0_copy = new THREE.Vector3().fromArray(p);
	var n = new THREE.Vector3().fromArray(normal);
	var n_copy = new THREE.Vector3().fromArray(normal);
	var v = new THREE.Vector3().fromArray([v1[0], v1[1], v1[2]]);
	
	var t = - n_copy.dot(p0);
	var s = p0_copy.add(new THREE.Vector3().fromArray(normal).multiplyScalar(t));
	
	return s.toArray();
}

/**
 * @param arrayToProcess: an array of points which have to be projected on the plane defined by a normal and an orientation vector
 * @param v1: one of the two orientation vectors of the plane, as an array
 * @return the projected points array
 */
function projectAllPointsOnPlane(arrayToProcess, normal, v1) {
	var projectedPoints = [];
	
	for (var i = arrayToProcess.length - 1; i >= 0; i--) {
		projectedPoints.push(projectPointOnPlane(arrayToProcess[i], normal, v1));
	}
	
	return projectedPoints;
}

/**
 * @param blackColor: a segment first vertex (here, should be the predefined black color point), segment which defines the first edge of the triangle
 * @param whiteColor: the segment second vertex (here, should be the predefined white color point)
 */
function findTriangleThirdVertex(blackColor, whiteColor, meanColor, arrayToProcess) {
	var blackPoint = clone(blackColor);
	var whitePoint = clone(whiteColor);
	var meanPoint = clone(meanColor);
	
	/* For each point in arrayToProcess, compute the angle formed with the chromatic axis */
	
	var blackToMean = new THREE.Vector3().fromArray([meanPoint[0] - blackPoint[0], meanPoint[1] - blackPoint[1], meanPoint[2] - blackPoint[2]]);
	var blackToWhite = new THREE.Vector3().fromArray([whitePoint[0] - blackPoint[0], whitePoint[1] - blackPoint[1], whitePoint[2] - blackPoint[2]]);
	var whiteToBlack = new THREE.Vector3().fromArray([-blackToWhite.x, -blackToWhite.y, -blackToWhite.z]);
	
	var normal = computeNormal(blackToWhite, blackToMean);
	
	var projectedPoints = projectAllPointsOnPlane(arrayToProcess, normal.toArray(), blackToMean.toArray());
	
	var maxAngleWBP = Number.MIN_VALUE;
	var maxAngleBWP = Number.MIN_VALUE;
	var maxBlackToPoint = new THREE.Vector3();
	var maxWhiteToPoint = new THREE.Vector3();
	
	for (var i = 0; i < projectedPoints.length; i++) {
		var point = projectedPoints[i];
		var blackToPoint = new THREE.Vector3().fromArray([point[0] - blackPoint[0], point[1] - blackPoint[1], point[2] - blackPoint[2]]);
		var whiteToPoint = new THREE.Vector3().fromArray([point[0] - whitePoint[0], point[1] - whitePoint[1], point[2] - whitePoint[2]]);
		
		var angleWBP = blackToWhite.angleTo(blackToPoint);
		var angleBWP = whiteToBlack.angleTo(whiteToPoint);
		
		if (maxAngleWBP < angleWBP) {
			maxAngleWBP = angleWBP;
			maxBlackToPoint = blackToPoint;
		}
		if (maxAngleBWP < angleBWP) {
			maxAngleBWP = angleBWP;
			maxWhiteToPoint = whiteToPoint;
		}
	}
	
	/* Determine the third point position, given the chromatic axis length and the two computed angles */
	
	var angleBPW = Math.PI - maxAngleWBP - maxAngleBWP;
	
	var distBlackToPoint = blackToWhite.length() / Math.sin(angleBPW);
	distBlackToPoint *= Math.sin(maxAngleBWP);
	
	var thirdPoint = new THREE.Vector3().fromArray(blackColor);
	thirdPoint.add(maxBlackToPoint.multiplyScalar(distBlackToPoint / maxBlackToPoint.length()));
	
	return {
		thirdPoint: thirdPoint,
		projectedPoints: projectedPoints
	};
}
