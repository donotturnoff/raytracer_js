
class Entity {    
    constructor(material) {
        this.material = material;
        this.translationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.xRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.yRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.zRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.rotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.scalingMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.invTranslationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.invXRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.invYRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.invZRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.invRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.invScalingMatrix = new Mat4x4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
    }
    
    computeIntersection(ray, x) {
        var d = this.invScalingMatrix.times(this.invRotationMatrix).times(ray.d);
        var o = this.invScalingMatrix.times(this.invRotationMatrix).times(this.invTranslationMatrix).times(ray.o);
        return this.translationMatrix.times(this.rotationMatrix).times(this.scalingMatrix).times(o.add(d.times(x)));
    }
        
    set translation(v) {
        this.translationMatrix = new Mat4x4([[1, 0, 0, v.x], [0, 1, 0, v.y], [0, 0, 1, v.z], [0, 0, 0, 1]]);
        this.invTranslationMatrix = this.translationMatrix.inverse;
    }
    
    set xRotation(theta) {
        var s = Math.sin(theta);
        var c = Math.cos(theta);
        this.xRotationMatrix = new Mat4x4([[1, 0, 0, 0], [0, c, -s, 0], [0, s, c, 0], [0, 0, 0, 1]]);
        this.invXRotationMatrix = this.xRotationMatrix.inverse;
        this.updateRotationMatrix();
    }
    
    set yRotation(theta) {
        var s = Math.sin(theta);
        var c = Math.cos(theta);
        this.yRotationMatrix = new Mat4x4([[c, 0, s, 0], [0, 1, 0, 0], [-s, 0, c, 0], [0, 0, 0, 1]]);
        this.invYRotationMatrix = this.yRotationMatrix.inverse;
        this.updateRotationMatrix();
    }
    
    set zRotation(theta) {
        var s = Math.sin(theta);
        var c = Math.cos(theta);
        this.zRotationMatrix = new Mat4x4([[c, -s, 0, 0], [s, c, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        this.invZRotationMatrix = this.zRotationMatrix.inverse;
        this.updateRotationMatrix();
    }
    
    updateRotationMatrix() {
        this.rotationMatrix = this.yRotationMatrix.times_m4x4(this.xRotationMatrix).times_m4x4(this.zRotationMatrix);
        this.invRotationMatrix = this.rotationMatrix.inverse;
    }
    
    set scaling(v) {
        this.scalingMatrix = new Mat4x4([[v.x, 0, 0, 0], [0, v.y, 0, 0], [0, 0, v.z, 0], [0, 0, 0, 1]]);
        this.invScalingMatrix = this.scalingMatrix.inverse;
    }
}

class Sphere extends Entity {
    constructor(r, material) {
        super(material);
        this.r = r; // Radius
    }

    computeIntersectionParam(ray) {
        var d = this.invScalingMatrix.times(this.invRotationMatrix).times(ray.d);
        var o = this.invScalingMatrix.times(this.invRotationMatrix).times(this.invTranslationMatrix).times(ray.o);
        
        var k0 = d.dot(d);
		var k1 = 2 * o.dot(d);
		var k2 = o.dot(o) - this.r*this.r;
		
		var discriminant = k1*k1 - 4*k0*k2;
		
		if (discriminant < 0) {
			return Infinity;
		} else {
			var s = Math.sqrt(discriminant);
			var t1 = (-k1 + s) / (2*k0);
			var t2 = (-k1 - s) / (2*k0);
			
			if (t1 < 0) {
				t1 = Infinity;
			}
			
			if (t2 < 0) {
				t2 = Infinity;
			}
			
			return Math.min(t1, t2);
		}
    }
    
    computeNormal(p) {
        return this.invScalingMatrix.times(this.invRotationMatrix).transpose.times(this.invScalingMatrix.times(this.invRotationMatrix).times(this.invTranslationMatrix).times(p)).normalised;
    }
}

class Torus extends Entity {
    constructor(R, r, material) {
        super(material);
        this.R = R; // Major radius
        this.r = r; // Minor radius
    }

    computeIntersectionParam(ray) {
        var d = this.invScalingMatrix.times(this.invRotationMatrix).times(ray.d);
        var o = this.invScalingMatrix.times(this.invRotationMatrix).times(this.invTranslationMatrix).times(ray.o);
        
        //http://hugi.scene.org/online/hugi24/coding%20graphics%20chris%20dragan%20raytracing%20shapes.htm
        
        var R2 = this.R*this.R;
        var r2 = this.r*this.r;
        
        var k0 = d.dot(d);
        var k1 = d.dot(o);
        var k2 = o.dot(o);

        var p = k2 - R2 - r2;

        var c4 = p*p - 4*R2*(r2 - o.y*o.y);
        var c3 = 4*k1*p + 8*R2*o.y*d.y;
        var c2 = 2*k0*p + 4*k1*k1 + 4*R2*d.y*d.y;
        var c1 = 4*k0*k1;
        var c0 = k0*k0;

        var s = solveQuartic(c0, c1, c2, c3, c4);
        s = s.filter(t => t > EPSILON);
        return Math.min(...s);
    }
    
    computeNormal(p) {
        //http://cosinekitty.com/raytrace/chapter13_torus.html
        
        var point = this.invScalingMatrix.times(this.invRotationMatrix).times(this.invTranslationMatrix).times(p);

        var x = point.x;
        var y = point.y;
        var z = point.z;
            
        var a = this.R/Math.sqrt(x*x+z*z);
        var N = new Vec3((1-a)*x, y, (1-a)*z);

        return this.invScalingMatrix.times(this.invRotationMatrix).transpose.times(N).normalised;
    }
}
