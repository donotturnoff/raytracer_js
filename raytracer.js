/*
 * To do:
 * Progressive rendering
 * Transparency and refraction
 */

const EPSILON = 0.0001;

var information = [];

class Ray {
    constructor(scene, o, d, bounces) {
        this.scene = scene;
        this.o = o;
        this.d = d;
        this.entity = null;
        this.intersectionParam = Infinity;
        this.intersection = null;
        this.bounces = bounces;
    }
    
    cast() {
        var minX = Infinity;
        for (var i = 0; i < this.scene.entities.length; i++) {
            var x = this.scene.entities[i].computeIntersectionParam(this);
            if (x < this.intersectionParam) {
                minX = x;
                this.entity = this.scene.entities[i];
                this.intersection = this.entity.computeIntersection(this, x);
                this.intersectionParam = this.intersection.distanceTo(this.o);
            }
        }
    }
    
    get colour() {
        if (this.entity == null) {
            return this.scene.backgroundColour;
        }
        var material = this.entity.material;
        var colour = {r: this.scene.ambientLight.r*material.ka.r, g: this.scene.ambientLight.g*material.ka.g, b: this.scene.ambientLight.b*material.ka.b};
        
        for (var i = 0; i < this.scene.lights.length; i++) {
            var l = this.scene.lights[i];
            
            var L_ = l.location.subtract(this.intersection);
            var L = L_.normalised;
            var N = this.entity.computeNormal(this.intersection);
            var R = N.times(2*(N.dot(L))).subtract(L).normalised;
            
            var shadowRayOrigin = this.intersection.add(L.times(EPSILON));
            
            var shadowRay = new Ray(this.scene, shadowRayOrigin, L, this.bounces);
            shadowRay.cast();
            
            if (shadowRay.entity == null || shadowRay.intersectionParam > L_.length) {
            
                var V = this.intersection.negated.normalised;
                
                var diff = Math.max(N.dot(L), 0);
                var spec = Math.pow(Math.max(R.dot(V), 0), material.se);
                
                colour.r += l.colour.r*material.kd.r*diff;
                colour.r += l.colour.r*material.ks.r*spec;
                colour.g += l.colour.g*material.kd.g*diff;
                colour.g += l.colour.g*material.ks.g*spec;
                colour.b += l.colour.b*material.kd.b*diff;
                colour.b += l.colour.b*material.ks.b*spec;
            }
            
        }
            
        if (material.refl > 0 && this.bounces < this.scene.maxBounces) {
            colour.r *= 1-material.refl;
            colour.g *= 1-material.refl;
            colour.b *= 1-material.refl;
            
            var reflectedRayDirection = this.d.subtract(N.times(2*(N.dot(this.d)))).normalised;
            var reflectedRayOrigin = this.intersection.add(reflectedRayDirection.times(EPSILON));
            var reflectedRay = new Ray(this.scene, reflectedRayOrigin, reflectedRayDirection, this.bounces+1);
            reflectedRay.cast();
            var reflectedColour = reflectedRay.colour;
            colour.r += reflectedColour.r*material.refl;
            colour.g += reflectedColour.g*material.refl;
            colour.b += reflectedColour.b*material.refl;
        }
        return colour;
    }
}

class RayTracer {
    constructor(ctx, scene) {
        this.ctx = ctx;
        this.scene = scene;
        
        var w = this.ctx.canvas.clientWidth;
        var h = this.ctx.canvas.clientHeight;
        
        this.info = []
        for (var y = 0; y < h; y++) {
            this.info.push([]);
            for (var x = 0; x < w; x++) {
                this.info[y].push({});
            }
        }
        
    }
    
    getInfo(x, y) {
        return this.info[y][x];
    }

    trace() {
        var w = this.ctx.canvas.clientWidth;
        var h = this.ctx.canvas.clientHeight;
        var blockSize = (this.scene.blockSize) ? this.scene.blockSize : 1;
        var samples = (this.scene.samples) ? this.scene.samples : 1;
        
        var o = new Vec3(0, 0, 0);
        
        for (var y = 0; y < h; y+=blockSize) {
            for (var x = 0; x < w; x+=blockSize) {
                var colours = [];
                var ss = Math.floor(Math.sqrt(samples));
                var intersections = [];
                for (var i = 0; i < ss; i++) {
                    for (var j = 0; j < ss; j++) {
                        var d = new Vec3((x+blockSize*i/ss-w/2)/w, (y+blockSize*j/ss-h/2)/w, 1);
                        var ray = new Ray(this.scene, o, d, 0);
                        ray.cast();
                        colours.push(ray.colour);
                        intersections.push(ray.intersection);
                    }
                }
                var colour = {r: 0, g: 0, b: 0};
                for (var i = 0; i < samples; i++) {
                    colour.r += colours[i].r;
                    colour.g += colours[i].g;
                    colour.b += colours[i].b;
                }
                colour.r /= samples;
                colour.g /= samples;
                colour.b /= samples;
                this.ctx.fillStyle = "rgb(" + colour.r*255 + ", " + colour.g*255 + ", " + colour.b*255 + ")";
                this.ctx.fillRect(x, h-y, blockSize, blockSize);
                for (var i = 0; i < blockSize; i++) {
                    for (var j = 0; j < blockSize; j++) {
                        this.info[h-(y+j)-1][x+i]={colour: colour, intersections: intersections};
                    }
                }
            }
        }
    }
}
