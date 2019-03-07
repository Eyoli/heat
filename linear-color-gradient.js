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

