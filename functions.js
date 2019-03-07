var Functions = {
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