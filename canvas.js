window.onload = function() {
	var ctx;
    var problem;
    var lastTimestamp = 0;
    var interval = 1000;
    var start = false;
    var T_MID = 50;
    var T_MAX = 250;
    var T_MIN = -273;
    
    var sFct = function(x) {
        if(x == 66) {
            return T_MAX;
        }
        return 0;
    };
    
    function sFct2d(x, y, t) {
        if((x - 50) * (x - 50) + (y - 50) * (y - 50) < 10*10) {
            return 1 + Math.cos(t * (2 * Math.PI) / 10) * T_MAX;
        }
        return 0;
    }
    
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
    
    class Diffusion {
        constructor(dt, alpha, T0, sFct) {
            this.dt = dt;
            this.alpha = alpha;
            this.T0 = T0;
            this.T = T0;
            this.sFct = sFct;
            this.t = 0;
            this.step = 0;
        }
    }
    
    class Diffusion1D extends Diffusion {
        
        constructor(dt, L, alpha, T0, sFct) {
            super(dt, alpha, T0, sFct);
            
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
        constructor(dt, L, H, alpha, T0, sFct) {
            super(dt, alpha, T0, sFct);
            
            // Col number
            this.N = T0[0].length;
            // Row number
            this.M = T0.length;
            
            this.H = H;
            this.L = L;
            
            var dy = H / this.M;
            var dx = L / this.N;
            
            this.lambdaX = this.alpha * this.dt / (dx * dx);
            this.lambdaY = this.alpha * this.dt / (dy * dy);
            
            // Row system
            this.aRow = utils.init(this.N - 2, this.lambdaX / this.dt);
            this.bRow = utils.init(this.N - 2, -2 * (this.lambdaX + 1) / this.dt);
            this.cRow = utils.init(this.N - 2, this.lambdaX / this.dt);
            this.dRow = utils.init(this.N, -this.lambdaY / this.dt);
            this.eRow = utils.init(this.N, 2 * (this.lambdaY - 1) / this.dt);
            this.fRow = utils.init(this.N, -this.lambdaY / this.dt);
            
            // Column system
            this.aCol = utils.init(this.M - 2, this.lambdaY / this.dt);
            this.bCol = utils.init(this.M - 2, -2 * (this.lambdaY + 1) / this.dt);
            this.cCol = utils.init(this.M - 2, this.lambdaY / this.dt);
            this.dCol = utils.init(this.M, -this.lambdaX / this.dt);
            this.eCol = utils.init(this.M, 2 * (this.lambdaX - 1) / this.dt);
            this.fCol = utils.init(this.M, -this.lambdaX / this.dt);
        }
        
        nextStep() {
            this.S = utils.discretise2d(this.sFct, {x: 0, y: 0}, {x: this.H, y: this.L}, {x: this.M, y: this.N}, this.t);
            
            this.TPrime = [];
            for(var i = 1; i < this.M - 1; i++) {
                var Si = utils.copy(this.S[i]);
                var Ti = utils.copy(this.T[i]);
                var d = utils.add(utils.tridiagonalMult(this.dRow, this.eRow, this.fRow, Ti), utils.multiply(Si, -1)).slice(1, this.N - 1);
                this.TPrime[i] = thomasAlgorithm(this.aRow, this.bRow, this.cRow, d);
                this.TPrime[i].splice(0, 0, this.T0[i][0]);
                this.TPrime[i].push(this.T0[i][this.N - 1]);
            }
            this.TPrime[0] = this.T0[0];
            this.TPrime[this.M - 1] = this.T0[this.M - 1];
            
            //console.log(this.TPrime);
            
            this.T = [];
            for(var j = 1; j < this.N - 1; j++) {
                var TjPrime = utils.column(this.TPrime, j);
                var Sj = utils.column(this.S, j);
                var d = utils.add(utils.tridiagonalMult(this.dCol, this.eCol, this.fCol, TjPrime), utils.multiply(Sj, -1)).slice(1, this.M - 1);
                var Tcol = thomasAlgorithm(this.aCol, this.bCol, this.cCol, d);
                Tcol.splice(0, 0, this.T0[0][j]);
                Tcol.push(this.T0[this.M - 1][j]);
                utils.setColumn(this.T, j, Tcol);
            }
            utils.setColumn(this.T, 0, utils.column(this.T0, 0));
            utils.setColumn(this.T, this.N - 1, utils.column(this.T0, this.N - 1));
            
            this.t += this.dt;
            this.step++;
        }
    }
	
	function init() {
		ctx = document.getElementById("myCanvas").getContext("2d");
        
        var N = 100;
        var L = 100;
        var alpha = 1;
        var dx = L / N;
                
        // time step complying with CFL conditions
        var dt = 0.9 * (0.5 * dx * dx) / alpha;
        
        //var T = utils.discretise(sFct, 0, L, N);
        var T = utils.discretise2d(sFct2d, {x: 0, y: 0}, {x: 100, y: 100}, {x: 100, y: 100}, 0);
        //T[0] = T_MAX;
        //T[N-1] = T_MAX;
        
        problem = new Diffusion2D(dt, 100, 100, alpha, T, sFct2d);
        //problem = new Diffusion1D(dt, L, alpha, T, sFct);
        
        document.getElementById("clear").onclick = function() {
            //problem = new Diffusion1D(dt, L, alpha, T, sFct);
            problem = new Diffusion2D(dt, 100, 100, alpha, T, sFct2d);
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
                var pMax = Math.min(T_MAX, problem.T[i][j]) / T_MAX;
                ctx.globalCompositeOperation = "lighter";
                
                ctx.fillStyle = "#ff0000";
                ctx.globalAlpha = pMid;
                ctx.beginPath();
                ctx.arc(20 + j * (760 / problem.T[i].length), 20 + i * (760 / problem.T.length) , 10, 0, 2 * Math.PI);
                ctx.fill();
                
                if(problem.T[i][j] > T_MID) {
                    
                    ctx.fillStyle = "#f1d800";
                    ctx.globalAlpha = pMax;
                    ctx.beginPath();
                    ctx.arc(20 + j * (760 / problem.T[i].length), 20 + i * (760 / problem.T.length) , 10, 0, 2 * Math.PI);
                    ctx.fill();
                    //ctx.globalCompositeOperation = "source-over";
                }
            }
        }
		
		window.requestAnimationFrame(draw);
	}
    
    init();
};

