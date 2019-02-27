window.onload = function() {
	var ctx;
    var problem;
    var lastTimestamp = 0;
    var interval = 1000;
    var start = false;
    var T_MID = 50;
    var T_MAX = 250;
    var T_MIN = -273;
    var gradient;
    var UNIT;
    
    function multiply(m, v) {
        var r = [];
        for(var i = 0; i < m.length; i++) {
            r[i] = 0;
            for(var j = 0; j < m[i].length; j++) {
                r[i] += m[i][j] * v[j];
            }
        }
        return r;
    }
    
    /**
     * Compute matrix system AX = B (A tridiagonal) using Thomas algorithm (O(n) complexity)
     * A represented by vectors a,b,c, and B represented by vector d
     */
    function thomasAlgorithm(a, b, c, d) {
        var cPrime = [];
        var dPrime = [];
        var n = a.length;
        var x = [];
        
        cPrime[0] = c[0] / b[0];
        dPrime[0] = d[0] / b[0];
        for(var i = 1; i < n; i++) {
            cPrime[i] = c[i] / (b[i] - a[i] * cPrime[i-1]);
            dPrime[i] = (d[i] - a[i] * dPrime[i-1]) / (b[i] - a[i] * cPrime[i-1]);
        }
                
        x[n-1] = dPrime[n-1];
        for(var i = n-2; i >= 0; i--) {
            x[i] = dPrime[i] - cPrime[i] * x[i+1];
        }
                
        return x;
    }
    
    var functions = {
        dirach: function(x) {
            if(x == 66) {
                return T_MAX;
            }
            return 0;
        },
        
        sinusoid: function(x, y, t) {
            var d2 = (x - 50) * (x - 50) + (y - 50) * (y - 50);
            var r2 = 10*10;
            if(d2 < r2) {
                return 5 + Math.cos(t * (2 * Math.PI) / 10) * (T_MID / 2) * (1 - d2/r2) * (1 - d2/r2);
            }
            return 0;
        }
    };
    
    var utils = {
        
        init: function(n, val) {
            var r = [];
            for(var i = 0; i < n; i++) {
                r[i] = val;
            }
            return r;
        },
        
        copy: function(v) {
            var r = [];
            for(var i = 0; i < v.length; i++) {
                r[i] = v[i];
            }
            return r;
        },
        
        multiply: function(v, scalar) {
            var r = [];
            for(var i = 0; i < v.length; i++) {
                r[i] = v[i] * scalar;
            }
            return r;
        },
        
        add: function(v1, v2) {
            var r = [];
            for(var i = 0; i < v1.length; i++) {
                r[i] = v1[i] + v2[i];
            }
            return r;
        },
        
        discretise: function(fct, min, max, n) {
            var r = [];
            var dx = (max - min) / n;
            for(var i = 0; i < n; i++) {
                r[i] = fct(min + i * dx);
            }
            return r;
        },
        
        discretise2d: function(fct, min, max, n, t) {
            var r = [];
            var dx = (max.x - min.x) / n.x;
            var dy = (max.y - min.y) / n.y;
            for(var i = 0; i < n.x; i++) {
                r[i] = [];
                for(var j = 0; j < n.y; j++) {
                    r[i][j] = fct(min.x + i * dx, min.y + j * dy, t);
                }
            }
            return r;
        },
        
        tridiagonalMult: function(a, b, c, v) {
            var r = [];
            r[0] = b[0] * v[0] + c[0] * v[1];
            for(var i = 1; i < v.length-1; i++) {
                r[i] = a[i] * v[i-1] + b[i] * v[i] + c[i] * v[i+1];
            }
            r[v.length-1] = a[v.length-1] * v[v.length-2] + b[v.length-1] * v[v.length-1];
            return r;
        },
        
        column: function(m, j) {
            var r = [];
            for(var i = 0; i < m.length; i++) {
                r[i] = m[i][j];
            }
            return r;
        },
        
        setColumn(m, j, v) {
            for(var i = 0; i < v.length; i++) {
                m[i] = m[i] || [];
                m[i][j] = v[i];
            }
        }
    };
    
    var matrix = {
        scalarMult: function(m, s) {
            var r = [];
            for(var i = 0; i < m.length; i++) {
                r[i] = [];
                for(var j = 0; j < m[i].length; j++) {
                    r[i][j] = m[i][j] * s;
                }
            }
            return r;
        },
        
        transpose: function(m) {
            var r = [];
            for(var i = 0; i < m.length; i++) {
                for(var j = 0; j < m[i].length; j++) {
                    r[j] = r[j] || [];
                    r[j][i] = m[i][j];
                }
            }
            return r;
        }
    };
    
    class ColorGradient {
        constructor(r1, g1, b1, r2, g2, b2) {
            this.ar = (r2 - r1);
            this.br = r1;
            this.ag = (g2 - g1);
            this.bg = g1;
            this.ab = (b2 - b1);
            this.bb = b1;
        }
        
        get(x) {
            var r = this.ar * x + this.br;
            var g = this.ag * x + this.bg;
            var b = this.ab * x + this.bb;
            return 'rgb(' + r + ',' + g + ',' + b + ')';
        }
        
        set(x, r, g, b) {
            var i = 0;
            while(i < this.steps.length && x > this.steps[i]) {
                i++;
            }
            steps.splice(i, 0, x);
        }
    }
    
    class Diffusion {
        constructor(dt, aFct, T0, sFct) {
            this.dt = dt;
            this.T0 = T0;
            this.T = T0;
            this.sFct = sFct;
            this.aFct = aFct;
            this.t = 0;
            this.step = 0;
        }
    }
    
    class Diffusion1D extends Diffusion {
        
        constructor(dt, L, aFct, T0, sFct) {
            super(dt, aFct, T0, sFct);
            
            this.N = T0.length;
            this.L = L;
            this.dx = L / this.N;
            
            this.lambda = this.alpha * this.dt / (this.dx * this.dx);
            
            this.a = utils.init(this.N - 2, -this.lambda);
            this.b = utils.init(this.N - 2, 1 + 2 * this.lambda);
            this.c = utils.init(this.N - 2, -this.lambda);
            
            //console.log(this.a, this.b, this.c);
        }
        
        nextStep() {
            // Discretise heat source function
            this.S = utils.discretise(this.sFct, 0, this.L, this.N);
            
            // Main term
            var d1 = utils.copy(this.T).slice(1, this.N - 1);
                
            // Edge conditions
            var d2 = utils.init(this.N - 2, 0);
            d2[0] = this.T0[0];
            d2[d2.length - 1] = this.T0[this.T0.length - 1];
            d2 = utils.multiply(d2, this.lambda);
            
            // Heat source term
            var d3 = utils.copy(this.S).slice(1, this.N - 1);
                
            var d = utils.add(utils.add(d1, d2), d3);
            
            this.T = thomasAlgorithm(this.a, this.b, this.c, d);
            this.T.splice(0, 0, this.T0[0]);
            this.T.push(this.T0[this.N - 1]);
                
            this.t += this.dt;
        }
    }
    
    class Diffusion2D extends Diffusion {
        constructor(dt, L, H, aFct, T0, sFct) {
            super(dt, aFct, T0, sFct);
            
            // Col number
            this.N = T0[0].length;
            // Row number
            this.M = T0.length;
            
            this.H = H;
            this.L = L;
            
            var dy = H / this.M;
            var dx = L / this.N;
            
            // Dispersal coefficient function
            this.alpha = utils.discretise2d(this.aFct, {x: 0, y: 0}, {x: this.H, y: this.L}, {x: this.M, y: this.N});
            
            var lambdaX = matrix.scalarMult(this.alpha, this.dt / (dx * dx));
            var lambdaY = matrix.scalarMult(this.alpha, this.dt / (dy * dy));
                                    
            // Row system (M * N)
            this.aRow = [];
            this.bRow = [];
            this.cRow = [];
            this.dRow = [];
            this.eRow = [];
            this.fRow = [];
            for(var i = 0; i < this.M; i++) {
                this.aRow[i] = utils.multiply(lambdaX[i], 1 / this.dt).slice(1, this.N - 1);
                this.bRow[i] = utils.add(utils.multiply(lambdaX[i], -2 / this.dt), utils.init(this.N, -2 / this.dt)).slice(1, this.N - 1);
                this.cRow[i] = utils.multiply(lambdaX[i], 1 / this.dt).slice(1, this.N - 1);
                this.dRow[i] = utils.multiply(lambdaY[i], -1 / this.dt);
                this.eRow[i] = utils.add(utils.multiply(lambdaY[i], 2 / this.dt), utils.init(this.N, -2 / this.dt));
                this.fRow[i] = utils.multiply(lambdaY[i], -1 / this.dt);
            }
            
            // Column system (N * M)
            this.aCol = [];
            this.bCol = [];
            this.cCol = [];
            this.dCol = [];
            this.eCol = [];
            this.fCol = [];
            for(var j = 0; j < this.N; j++) {
                this.aCol[j] = utils.multiply(utils.column(lambdaY, j), 1 / this.dt).slice(1, this.N - 1);
                this.bCol[j] = utils.add(utils.multiply(utils.column(lambdaY, j), -2 / this.dt), utils.init(this.N, -2 / this.dt)).slice(1, this.N - 1);
                this.cCol[j] = utils.multiply(utils.column(lambdaY, j), 1 / this.dt).slice(1, this.N - 1);
                this.dCol[j] = utils.multiply(utils.column(lambdaX, j), -1 / this.dt);
                this.eCol[j] = utils.add(utils.multiply(utils.column(lambdaX, j), 2 / this.dt), utils.init(this.N, -2 / this.dt));
                this.fCol[j] = utils.multiply(utils.column(lambdaX, j), -1 / this.dt);
            }
        }
        
        nextStep() {
            this.S = utils.discretise2d(this.sFct, {x: 0, y: 0}, {x: this.H, y: this.L}, {x: this.M, y: this.N}, this.t);
            
            this.TPrime = [];
            for(var i = 1; i < this.M - 1; i++) {
                var Si = utils.copy(this.S[i]);
                var Ti = utils.copy(this.T[i]);
                var d = utils.add(utils.tridiagonalMult(this.dRow[i], this.eRow[i], this.fRow[i], Ti), utils.multiply(Si, -1)).slice(1, this.N - 1);
                this.TPrime[i] = thomasAlgorithm(this.aRow[i], this.bRow[i], this.cRow[i], d);
                this.TPrime[i].splice(0, 0, this.T0[i][0]);
                this.TPrime[i].push(this.T0[i][this.N - 1]);
            }
            this.TPrime[0] = this.T0[0];
            this.TPrime[this.M - 1] = this.T0[this.M - 1];
            
            //console.log(this.TPrime);
            
            this.t += (this.dt / 2);
            this.S = utils.discretise2d(this.sFct, {x: 0, y: 0}, {x: this.H, y: this.L}, {x: this.M, y: this.N}, this.t);
            
            this.T = [];
            for(var j = 1; j < this.N - 1; j++) {
                var TjPrime = utils.column(this.TPrime, j);var Sj = utils.column(this.S, j);
                var d = utils.add(utils.tridiagonalMult(this.dCol[j], this.eCol[j], this.fCol[j], TjPrime), utils.multiply(Sj, -1)).slice(1, this.M - 1);
                var Tcol = thomasAlgorithm(this.aCol[j], this.bCol[j], this.cCol[j], d);
                Tcol.splice(0, 0, this.T0[0][j]);
                Tcol.push(this.T0[this.M - 1][j]);
                utils.setColumn(this.T, j, Tcol);
            }
            utils.setColumn(this.T, 0, utils.column(this.T0, 0));
            utils.setColumn(this.T, this.N - 1, utils.column(this.T0, this.N - 1));
            
            this.t += (this.dt / 2);
            this.step++;
        }
    }
	
	function init() {
		ctx = document.getElementById("myCanvas").getContext("2d");
        
        gradient = new ColorGradient(0, 188, 212, 255, 0, 0);
        
        var N = 100;
        var L = 100;
        var dx = L / N;
        UNIT = 760 / N;
                
        // time step complying with CFL conditions
        var dt = 0.9 * (0.5 * dx * dx);
        
        var aFct = function(x, y, t) {
            var d = x*x + y*y;
            return d < 100*100 && d > 50*50 ? 5 : 1;
        };
        
        //var T = utils.discretise(sFct, 0, L, N);
        var T = utils.discretise2d(functions.sinusoid, {x: 0, y: 0}, {x: 100, y: 100}, {x: 100, y: 100}, 0);
        
        problem = new Diffusion2D(dt, 100, 100, aFct, T, functions.sinusoid);
        //problem = new Diffusion1D(dt, L, alpha, T, sFct);
        
        document.getElementById("clear").onclick = function() {
            //problem = new Diffusion1D(dt, L, alpha, T, sFct);
            problem = new Diffusion2D(dt, 100, 100, aFct, T, functions.sinusoid);
        }
        
        document.getElementById("start").onclick = function() {
            interval = Number.parseInt(document.getElementById("interval").value);
            start = !start;
        };
        
        document.getElementById("step").onclick = function() {
            problem.nextStep();
        };
        		
		window.requestAnimationFrame(draw);
	}
    
	function draw(timestamp) {
		if(start && timestamp - lastTimestamp >= interval) {
            problem.nextStep();
            document.getElementById("number").innerHTML = problem.step;
            lastTimestamp = timestamp;
        }
        
        ctx.clearRect(0, 0, 1000, 800);
                
        for(var i = 0; i < problem.T.length; i++) {
            for(var j = 0; j < problem.T[i].length; j++) {
                
                var pMid = Math.min(T_MID, problem.T[i][j]) / T_MID;
                var pMax = Math.min(T_MAX - T_MID, problem.T[i][j] - T_MID) / (T_MAX - T_MID);
                ctx.globalCompositeOperation = "lighter";
                
                ctx.fillStyle = gradient.get(pMid);
                ctx.globalAlpha = Math.sqrt(pMid);
                ctx.beginPath();
                //ctx.arc(20 + j * (760 / problem.T[i].length), 20 + i * (760 / problem.T.length), 10, 0, 2 * Math.PI);
                ctx.fillRect(20 + j * UNIT, 20 + i * UNIT, UNIT, UNIT);
                ctx.fill();
            }
        }
		
		window.requestAnimationFrame(draw);
	}
    
    init();
};

