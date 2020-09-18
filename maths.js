/* Adapted from https://github.com/erich666/GraphicsGems/blob/master/gems/Roots3And4.c */

function isZero(x) {
    return Math.abs(x) < 1e-9;
}

function solveQuadratic(c0, c1, c2) {
    var s = [];
    var root = Math.sqrt(c1*c1-4*c0*c2);
    var denom = 2*c0;
    if (!isNaN(root)) {
        s.push((-c1+root)/denom);
        if (root != 0) {
            s.push((-c1-root)/denom);
        }
    }
    return s;
}

function solveCubic(c0, c1, c2, c3) {

    /* normal form: x^3 + Ax^2 + Bx + C = 0 */

    var A = c1/c0;
    var B = c2/c0;
    var C = c3/c0;

    /*  substitute x = y - A/3 to eliminate quadatric term:
	x^3 +px + q = 0 */

    var A2 = A * A;
    var p = 1.0 / 3 * (- 1.0 / 3 * A2 + B);
    var q = 1.0 / 2 * (2.0 / 27 * A * A2 - 1.0 / 3 * A * B + C);

    /* use Cardano's formula */

    var p3 = p * p * p;
    var D = q * q + p3;

    var s = [];

    if (isZero(D)) {
        if (isZero(q)) /* one triple solution */ {
            s = [0];
        } else /* one single and one double solution */ {
            let u = Math.cbrt(-q);
            s = [2 * u, -u];
        }
    } else if (D < 0) /* Casus irreducibilis: three real solutions */ {
        var phi = 1.0 / 3 * Math.acos(-q / Math.sqrt(-p3));
        var t = 2 * Math.sqrt(-p);

        s = [t * Math.cos(phi),
            -t * Math.cos(phi + Math.PI / 3),
            -t * Math.cos(phi - Math.PI / 3)];

    } else /* one real solution */ {
        var sqrt_D = Math.sqrt(D);
        var u = Math.cbrt(sqrt_D - q);
        var v = -Math.cbrt(sqrt_D + q);

        s = [u + v];

    }

    /* resubstitute */

    var sub = 1.0 / 3 * A;

    for (var i = 0; i < s.length; ++i) {
        s[i] -= sub;
    }

    return s;
}

function solveQuartic(c0, c1, c2, c3, c4) {
    var s = [];
    var coeffs = [];
    var num;

    /* Normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 */

    var A = c1/c0;
    var B = c2/c0;
    var C = c3/c0;
    var D = c4/c0;

    /*  Substitute x = y - A/4 to eliminate cubic term:
	x^4 + px^2 + qx + r = 0 */

    var A2 = A * A;
    var p = - 3.0/8 * A2 + B;
    var q = 1.0/8 * A2 * A - 1.0/2 * A * B + C;
    var r = - 3.0/256*A2*A2 + 1.0/16*A2*B - 1.0/4*A*C + D;

    if (isZero(r)) {
        /* No absolute term: y(y^3 + py + q) = 0 */
        s = solveCubic(1, 0, p, q);
        s.push(0);
    } else {
        /* solve the resolvent cubic ... */

        var x3 = 1.0/2 * r * p - 1.0/8 * q * q;
        var x2 = -r;
        var x1 = -1.0/2 * p;
        var x0 = 1;

        s = solveCubic(x0, x1, x2, x3);

        /* ... and take the one real solution ... */

        var z = s[0];

        /* ... to build two quadratic equations */

        var u = z * z - r;
        var v = 2 * z - p;

        if (isZero(u)) {
            u = 0;
        } else if (u > 0) {
            u = Math.sqrt(u);
        } else {
            return [];
        }

        if (isZero(v)) {
            v = 0;
        } else if (v > 0) {
            v = Math.sqrt(v);
        } else {
            return [];
        }

        x2 = z - u;
        x1 = q < 0 ? -v : v;
        x0 = 1;

        s = solveQuadratic(x0, x1, x2);

        x2 = z + u;
        x1 = q < 0 ? v : -v;
        x0 = 1;

        s = s.concat(solveQuadratic(x0, x1, x2));
    }

    /* Resubstitute */

    var sub = 1.0/4 * A;

    for (var i = 0; i < s.length; ++i) {
        s[i] -= sub;
    }

    return s;
}

class Mat4x4 {
    constructor(a) {
        this.a = a;
    }
    
    times(m) {
        if (m instanceof Vec3) {
            return this.times_v3(m);
        } else {
            return this.times_m4x4(m);
        }
    }
    
    times_m4x4(m) {
        var ra = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
        for (var i = 0; i < 4; i++) {
            for (var j = 0; j < 4; j++) {
                for (var k = 0; k < 4; k++) {
                    ra[i][j] += this.a[i][k] * m.a[k][j];
                }
            }
        }
        return new Mat4x4(ra);
    }
    
    times_v3(v) {
        var va = [v.x, v.y, v.z, 1];
        var ra = [0, 0, 0, 0];
        for (var i = 0; i < 4; i++) {
            for (var k = 0; k < 4; k++) {
                ra[i] += this.a[i][k] * va[k];
            }
        }
        return new Vec3(ra[0], ra[1], ra[2]);
    }
    
    get inverse() {
        var a = this.a;
        var a2323 = a[2][2] * a[3][3] - a[2][3] * a[3][2];
        var a1323 = a[2][1] * a[3][3] - a[2][3] * a[3][1];
        var a1223 = a[2][1] * a[3][2] - a[2][2] * a[3][1];
        var a0323 = a[2][0] * a[3][3] - a[2][3] * a[3][0];
        var a0223 = a[2][0] * a[3][2] - a[2][2] * a[3][0] ;
        var a0123 = a[2][0] * a[3][1] - a[2][1] * a[3][0] ;
        var a2313 = a[1][2] * a[3][3] - a[1][3] * a[3][2] ;
        var a1313 = a[1][1] * a[3][3] - a[1][3] * a[3][1] ;
        var a1213 = a[1][1] * a[3][2] - a[1][2] * a[3][1] ;
        var a2312 = a[1][2] * a[2][3] - a[1][3] * a[2][2] ;
        var a1312 = a[1][1] * a[2][3] - a[1][3] * a[2][1] ;
        var a1212 = a[1][1] * a[2][2] - a[1][2] * a[2][1] ;
        var a0313 = a[1][0] * a[3][3] - a[1][3] * a[3][0] ;
        var a0213 = a[1][0] * a[3][2] - a[1][2] * a[3][0] ;
        var a0312 = a[1][0] * a[2][3] - a[1][3] * a[2][0] ;
        var a0212 = a[1][0] * a[2][2] - a[1][2] * a[2][0] ;
        var a0113 = a[1][0] * a[3][1] - a[1][1] * a[3][0] ;
        var a0112 = a[1][0] * a[2][1] - a[1][1] * a[2][0] ;

        var det = a[0][0] * ( a[1][1] * a2323 - a[1][2] * a1323 + a[1][3] * a1223 ) 
            - a[0][1] * ( a[1][0] * a2323 - a[1][2] * a0323 + a[1][3] * a0223 ) 
            + a[0][2] * ( a[1][0] * a1323 - a[1][1] * a0323 + a[1][3] * a0123 ) 
            - a[0][3] * ( a[1][0] * a1223 - a[1][1] * a0223 + a[1][2] * a0123 ) ;
        det = 1 / det;

        return new Mat4x4([
            [
                det *   ( a[1][1] * a2323 - a[1][2] * a1323 + a[1][3] * a1223 ),
                det * - ( a[0][1] * a2323 - a[0][2] * a1323 + a[0][3] * a1223 ),
                det *   ( a[0][1] * a2313 - a[0][2] * a1313 + a[0][3] * a1213 ),
                det * - ( a[0][1] * a2312 - a[0][2] * a1312 + a[0][3] * a1212 )
            ], [
                det * - ( a[1][0] * a2323 - a[1][2] * a0323 + a[1][3] * a0223 ),
                det *   ( a[0][0] * a2323 - a[0][2] * a0323 + a[0][3] * a0223 ),
                det * - ( a[0][0] * a2313 - a[0][2] * a0313 + a[0][3] * a0213 ),
                det *   ( a[0][0] * a2312 - a[0][2] * a0312 + a[0][3] * a0212 )
            ], [
                det *   ( a[1][0] * a1323 - a[1][1] * a0323 + a[1][3] * a0123 ),
                det * - ( a[0][0] * a1323 - a[0][1] * a0323 + a[0][3] * a0123 ),
                det *   ( a[0][0] * a1313 - a[0][1] * a0313 + a[0][3] * a0113 ),
                det * - ( a[0][0] * a1312 - a[0][1] * a0312 + a[0][3] * a0112 )
            ], [
                det * - ( a[1][0] * a1223 - a[1][1] * a0223 + a[1][2] * a0123 ),
                det *   ( a[0][0] * a1223 - a[0][1] * a0223 + a[0][2] * a0123 ),
                det * - ( a[0][0] * a1213 - a[0][1] * a0213 + a[0][2] * a0113 ),
                det *   ( a[0][0] * a1212 - a[0][1] * a0212 + a[0][2] * a0112 )
            ]
        ]);
    }
    
    get transpose() {
        var ra = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
        for (var i = 0; i < 4; i++) {
            for (var j = 0; j < 4; j++) {
                ra[j][i] = this.a[i][j];
            }
        }
        return new Mat4x4(ra);
    }
}

class Vec3 {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    
    add(v) {
        return new Vec3(this.x+v.x, this.y+v.y, this.z+v.z);
    }
    
    subtract(v) {
        return new Vec3(this.x-v.x, this.y-v.y, this.z-v.z);
    }
    
    times(p) {
        return new Vec3(this.x*p, this.y*p, this.z*p);
    }
    
    get length() {
        return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    }
    
    get normalised() {
        return new Vec3(this.x/this.length, this.y/this.length, this.z/this.length);
    }
    
    get negated() {
        return new Vec3(-this.x, -this.y, -this.z);
    }
    
    dot(v) {
        return this.x*v.x+this.y*v.y+this.z*v.z;
    }
    
    distanceTo(v) {
        return this.subtract(v).length;
    }
}
