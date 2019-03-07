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
            return d < 100*100 && d > 50*50 ? 10 : 0;
        };
        
        var sFct = function(x, y, t) {
            return Functions.sinusoid(30, 30, 10, 10, 1, T_MID / 3, x, y, t) + Functions.sinusoid(75, 75, 10, 10, 5, T_MID / 2, x, y, t);
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

