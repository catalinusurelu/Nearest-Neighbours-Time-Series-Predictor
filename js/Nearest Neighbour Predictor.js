// Catalin Constantin Usurelu
// 341 C1

/* Theory:
	According to the Takens’ theorem [1] a sufficient embedding dimension to reconstruct
	a phase space topologically equivalent to the original space is 2m+1, where m is the
	dimension of the attractor of the system, and the embedding can be done by delayed
	coordinates of a single observable. Takens’ theorem states the sufficient condition of the
	embedding dimension, however, for many systems the embedding can be reached in fewer
	dimensions.

	Rtol - Threashold for measuring the separation (distance) between 2 points.
           If it is too small, some true false nearest neighbours will be counted as false,
           while if it's too big, some false nearest neighoburs will be counted as true.

	Atol (loneliness tolerance) - There is a problem with the above formula: what if two
	                              points are nearest neighbours, but the distance between
	                              them is big? The Atol criteria is used.
	Ra - Some attractor's size. It is set as the standard deviation of the input samples.
*/

function minimumEmbeddingDimension(y, Rtol, Atol, f_limit, maxM) {

    // Handle default values
    if(typeof(maxM)==='undefined') maxM = 10;
    if(typeof(Rtol)==='undefined') Rtol = 10;
    if(typeof(Atol)==='undefined') Atol = 2;
    if(typeof(f_limit)==='undefined') f_limit = 0.1;

    var Ra = standardDeviation(y);
    var m = 2;
    var N = y.length - 1;

    var minF = 1;
    var bestM = 2;

    do {
        var Y = new Array(N);
        var Ynext = new Array(N); // Ynext - are embedding dimension mai mare cu 1 decat Y
        var falseNearestNeighbours = 0;
        var nearestNeighbours = 0;
        var f;

        // Y[i] = (y[i], y[i-1], ... , y[i-(m-1)]); i = m - 1 ... N
        for(var i = m - 1; i <= N - 1; i++) {
            Y[i] = y.slice(i - (m - 1), i + 1);
        }

        // Ynext[i] = (y[i], y[i-1], ... , y[i-(m-1)]); i = m - 1 ... N
        for(var i = m; i <= N - 1; i++) {
            Ynext[i] = y.slice(i - m, i + 1);
        }

        // Incepem de la m ca sa cuprindem si Ynext
        for(var i = m; i <= N - 1; i++)
        {
            var R = -1;
            var nearestNeighbour;
            var Rnext;

            // Find nearest neighbour.
            for(var j = m; j <= N - 1; j++)
                if(i != j) {
                	// Initiate.
                	if(R == -1) {
                		nearestNeighbour = j;
                        R = dist(Y[i], Y[j]);
                        continue;
                	}

                    if(R > dist(Y[i], Y[j])) {  
                        nearestNeighbour = j;
                        R = dist(Y[i], Y[j]);
                    }
                }

            // Distance between neigbhours in m + 1 dimensions.
            Rnext = dist(Ynext[i], Ynext[nearestNeighbour]);

	        var Crieriu_Rtol = Math.sqrt(((Math.pow(Rnext, 2) - Math.pow(R, 2)) / Math.pow(R, 2)));
	        var Criteriu_Atol = Rnext / Ra;

	        if(Crieriu_Rtol > Rtol)
	            falseNearestNeighbours++;

	        nearestNeighbours++;
	    }

	        
        // Every point has exactly one nearest neighbour. Includes duplicates (x,y), (y,x).
        //nearestNeighbours = (N - 1) - m + 1;
        f = falseNearestNeighbours / nearestNeighbours;
        m++;

        // Find the minnimum embedding dimnesion for which f is minimum.
        if(f < minF) {
        	minF = f;
        	bestM = m;
        }

        if(m > maxM) {
        	break;
        }
    }
    while(f >= f_limit);

    // m - 1 to take into account the artificial increased in embedding dimension caused by noise.
    return bestM - 1;
}

//  phi = maximum angle between trajectories.
//  k = number of nearest neighbours to be considered. 
//      k should be < (m - (N - 1) + 1) (that is cardinal of [m ... N - 1])).
function nearestNeighboursPredictor(trainingSamples, nrOfPredictions, k, phi, rtol, atol, f_lim) {

    if(typeof(k)==='undefined') k = 25;
    if(typeof(phi)==='undefined') phi = 30;
    
    var m = minimumEmbeddingDimension(trainingSamples, rtol, atol, f_lim);
    var y = trainingSamples.slice(); // clone array.
    var N = y.length - 1;
    var Y = [];
    var delta_y = [];

    for(var i = m - 1; i <= N - 1; i++) {

    	// Y[i] = (y[i], y[i-1], ... , y[i-(m-1)]); i = m - 1 ... N
        Y[i] = y.slice(i - (m - 1), i + 1);

        // Delta_y - are 'labels' associated with each Y[i]. They represent 
        //           'points' in the velocity vector field .
        delta_y[i] = y[i + 1] - y[i];
    }

    // Y[N] = (y[N], y[N-1], ... , y[N-(m-1)])
    Y[N] = y.slice(N - (m - 1), N + 1);

    // Actual algorithm start.

    for(predictions = 0; predictions < nrOfPredictions; predictions++) {

        // Find k nearest neighbours.
        var k_nearest = new MaxKNeighboursQueue(k);
        for(var i = m; i <= N - 1; i++) {


            if(angle( (substract(Y[i - 1], Y[i])),
                      (substract(Y[N - 1], Y[N])) ) > phi)
                continue;

            // Store only the last k nearest neighbours.
            k_nearest.push(i, dist(Y[N], Y[i]));
        }

        // If we don't have enough neighbours try to increase the angle and find other neighbours.
        var originalPhi = phi;
        if(k_nearest.length() < k) {
            while(k_nearest.length() < k) {
                phi++;

                var k2 = k - k_nearest.length();
                var k_nearest2 = new MaxKNeighboursQueue(k2);

                for(var i = m; i <= N - 1; i++) {
                    if(Math.abs(angle( (substract(Y[i - 1], Y[i])),
                              (substract(Y[N - 1], Y[N])) ) ) > phi ||  
                        Math.abs(angle( (substract(Y[i - 1], Y[i])),
                              (substract(Y[N - 1], Y[N]))) ) < originalPhi
                        ) 
                        continue;

                    // Store only the last (best) k nearest neighbours.
                    k_nearest2.push(i, dist(Y[N], Y[i]));
                }

                for(var j = 0; j < k_nearest2.length(); j++) {
                    k_nearest.push(k_nearest2.get(j).neighbour, k_nearest2.get(j).distance);
                }
            }
        }

        // Compute weights.
        var q = [];


        var maxDistanceSquared = Math.pow(k_nearest.getMax().distance, 2);

        for(var i = 0; i < k_nearest.length(); i++) {
            q[i] = Math.pow(1 - Math.pow(k_nearest.get(i).distance, 2)/
                                maxDistanceSquared, 
                            2);
        }

        // Compute delta_y[n + 1]

        var numerator = 0;
        var denominator = 0;

        for(var i = 0; i < k_nearest.length(); i++) {
            numerator += q[i] * delta_y[k_nearest.get(i).neighbour];
            denominator += q[i];
        }

        delta_y[N] = numerator / denominator;

        // Compute next predicted y.
        y[N + 1] = Y[N][Y[N].length - 1] + delta_y[N];

        // Increase number of samples (original + predicted).
        N++;

        Y[N] = y.slice(N - (m - 1), N + 1);
    }

    // Return the predicted part of y.
    return y.slice(trainingSamples.length, trainingSamples.length + 1 + nrOfPredictions);
}




/*
-----------------------------------------------------------
     Max k queue - queue holding k values with maximum size.

*/

MaxKNeighboursQueue = function(k) {
    this.queue = [];

    this.push = function(neighbour, distance) {
        if(this.queue.length < k) {
            this.queue.push({neighbour: neighbour, distance: distance});
        }
        else if(this.queue.length >= k) { 

            this.queue.push( {neighbour: neighbour, distance: distance} );

            this.queue.sort(function(x, y) { return x.distance - y.distance })

            this.queue = this.queue.slice(0, this.queue.length - 1); // Remove k + 1th element
        }
    }

    this.contains = function(neighbour) {
    	// Could have done a binary search ...
    	for(var i = 0; i < this.queue.length; i++) {
    		if(this.queue[i].neighbour == neighbour) {
    			return true;
    		}
    	}
    	return false;
    }

    this.getMax = function() {
        return this.queue[this.queue.length - 1];
    }

    this.getMaxKNeighbours = function() {
        return this.queue;
    }

    this.get = function(i) {
        return this.queue[i];
    }

    this.length = function() {
    	return this.queue.length;
    }
}



/*
-----------------------------------------------------------
     Math stuff
*/

function dist(v1, v2) {
    var sum = 0;
    for(var i = 0; i < v1.length; i++) {
        sum += Math.pow((v1[i] - v2[i]), 2);
    }

    return Math.sqrt(sum);
}

function standardDeviation(values) {
  var avg = average(values);
  
  var squareDiffs = values.map(function(value) {
    var diff = value - avg;
    var squareDiff = diff * diff;
    return squareDiff;
  });
  
  var avgSquareDiff = average(squareDiffs);
 
  var stdDev = Math.sqrt(avgSquareDiff);
  return stdDev;
}
 
function average(values) {
  var sum = values.reduce(function(sum, x) {
    return sum + x;
  }, 0);
 
  var avg = sum / values.length;
  return avg;
}


function dotProduct(v1, v2) {
    var sum = 0;
    for(var i = 0; i < v1.length; i++) {
        sum += v1[i] * v2[i];
    }

    return sum;
}

function length(v) {
    var sum = 0;
    for(var i = 0; i < v.length; i++) {
        sum += Math.pow(v[i], 2);
    }

    return Math.sqrt(sum);
}

function angle(v1, v2) {
    var radians = Math.acos(dotProduct(v1, v2) / (length(v1) * length(v2)));
    var degrees = radians * 180 / Math.PI;
    return degrees;
}

function substract(v1, v2) {
    var rez = [];
    for(var i = 0; i < v1.length; i++) {
        rez[i] = v1[i] - v2[i];
    }

    return rez;
}

function meanSquaredError(v1, v2) {
    var rez = 0;
    for(var i = 0; i < v1.length; i++) {
        rez += Math.pow(v1[i] - v2[i], 2);
    }
    rez /= v1.length;

    return rez;
}