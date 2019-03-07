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
        
        sinusoid: function(cx, cy, r, p, aMin, aMax, x, y, t) {
            var d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
            var r2 = r*r;
            if(d2 < r2) {
                return aMin + Math.cos(t * (2 * Math.PI) / p) * aMax * (1 - d2/r2) * (1 - d2/r2);
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
        },
        
        minFilter(m, v) {
            for(var i = 0; i < m.length; i++) {
                for(var j = 0; j < m[i].length; j++) {
                    m[j][i] = Math.max(0, m[i][j]);
                }
            }
        }
    };
    
    class LinearColorGradient {
        constructor(c1, c2) {
            this.X = [];
            this.C = [];
            this.A = [];
            this.B = [];
            
            this.set(0, c1);
            this.set(1, c2);
        }
        
        get(x) {
            var i = 1;
            while(i < this.X.length && x > this.X[i]) {
                i++;
            }
            var r = this.A[i-1].r * x + this.B[i-1].r;
            var g = this.A[i-1].g * x + this.B[i-1].g;
            var b = this.A[i-1].b * x + this.B[i-1].b;
            return 'rgb(' + r + ',' + g + ',' + b + ')';
        }
        
        set(x, c) {
            var i = 0;
            while(i < this.X.length && x > this.X[i]) {
                i++;
            }
            this.X.splice(i, 0, x);
            this.C.splice(i, 0, c);
            // Compute linear coefficients for each segment (a, b)
            while(i > 0 && i < this.X.length) {
                this.A[i-1] = {
                    r: (this.C[i].r - this.C[i-1].r) / (this.X[i] - this.X[i-1]), 
                    g: (this.C[i].g - this.C[i-1].g) / (this.X[i] - this.X[i-1]), 
                    b: (this.C[i].b - this.C[i-1].b) / (this.X[i] - this.X[i-1])
                };
                this.B[i-1] = {
                    r: this.C[i-1].r - this.A[i-1].r * this.X[i-1], 
                    g: this.C[i-1].g - this.A[i-1].g * this.X[i-1], 
                    b: this.C[i-1].b - this.A[i-1].b * this.X[i-1]
                };
                i++;
            }
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
    
    class Diffusion2D extends Diffusion {
        constructor(dt, dx, dy, aFct, T0, sFct) {
            super(dt, aFct, T0, sFct);
            
            // Col number
            this.N = T0[0].length;
            // Row number
            this.M = T0.length;
            
            this.H = this.M * dy;
            this.L = this.N * dx;
            
            var dy = dy;
            var dx = dx;
            
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
            //this.S = utils.discretise2d(this.sFct, {x: 0, y: 0}, {x: this.H, y: this.L}, {x: this.M, y: this.N}, this.t);
            
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
            
            //matrix.minFilter(this.T, 0);
            
            this.t += (this.dt / 2);
            this.step++;
        }
    }
	
	function init() {
		ctx = document.getElementById("myCanvas").getContext("2d");
        
        gradient = new LinearColorGradient({r: 0, g: 188, b: 212}, {r:255,g:235,b:59});
        gradient.set(0.2, {r: 146, g: 138, b: 0});
        gradient.set(0.6, {r: 255, g: 0, b: 0});
        console.log(gradient);
        
        var N = 100;
        var M = 100;
        var H = 100;
        var L = 100;
        UNIT = 760 / N;

        var dy = H / M;
        var dx = L / N;
                
        // time step complying with CFL conditions
        var dt = 0.9 * (0.5 * dx * dx);
        
        var aFct = function(x, y, t) {
            var d = x*x + y*y;
            return d < 100*100 && d > 50*50 ? 10 : 0.1;
        };
        
        var sFct = function(x, y, t) {
            return functions.sinusoid(30, 30, 10, 10, 1, T_MID / 3, x, y, t) + functions.sinusoid(75, 75, 10, 10, 5, T_MID / 2, x, y, t);
        };
        
        var T = utils.discretise2d(sFct, {x: 0, y: 0}, {x: L, y: L}, {x: L, y: L}, 0);
        
        problem = new Diffusion2D(dt, dx, dy, aFct, T, sFct);
        
        document.getElementById("clear").onclick = function() {
            //problem = new Diffusion1D(dt, L, alpha, T, sFct);
            problem = new Diffusion2D(dt, dx, dy, aFct, T, sFct);
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
        
        ctx.clearRect(0, 0, 800, 800);
                
        for(var i = 0; i < problem.T.length; i++) {
            for(var j = 0; j < problem.T[i].length; j++) {
                
                var p = Math.min(T_MAX, problem.T[i][j]) / T_MAX;
                ctx.globalCompositeOperation = "lighter";
                
                ctx.fillStyle = gradient.get(p);
                ctx.globalAlpha = 0.8 * Math.sqrt(p);
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

