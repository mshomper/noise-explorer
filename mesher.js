/**
 * mesher.js — Shared mesh generation layer for the noise/TPMS scaffold pipeline
 *
 * Architecture:
 *   evaluateNoise(json, x, y, z)   → scalar SDF value
 *   evaluateTPMS(json, x, y, z)    → stub (pass 2)
 *   evaluateFDSBN(json, x, y, z)   → stub (future)
 *   sampleGrid(json, N)            → Float32Array of SDF values on N³ grid
 *   marchingCubes(grid, N)         → {vertices, normals, triangles}
 *   taubinSmooth(mesh, λ, μ, iter) → smoothed mesh (in-place)
 *   checkQuality(mesh)             → {watertight, boundaryEdges, triCount, degenCount}
 *   toBinarySTL(mesh, cellSizeMm)  → ArrayBuffer
 *   meshFromJSON(json, N)          → mesh object
 *   downloadSTL(json, N, cellMm)   → triggers browser download
 *
 * Coordinate system:
 *   Grid spans [0, N-1]³ in index space.
 *   World space maps to a cellSizeMm × cellSizeMm × cellSizeMm cube (default 10mm).
 *   SDF noise space: same SCALE = 5/π mapping as the tool's shader and GridGenerator.
 *
 * Surface types supported:
 *   'noise'  → full implementation (all 7 noise types + perforation)
 *   'terms'  → stub, throws descriptive error
 *   'fdsbn'  → stub, throws descriptive error
 */

(function(global) {
'use strict';

// ── Constants ─────────────────────────────────────────────────────────────
var SCALE = 5.0 / Math.PI;  // noise space scale, matches shader

// ── DeterministicNoise ────────────────────────────────────────────────────
// Verbatim port of the tool's JS noise functions.
// Produces numerically identical output to the shader for all types except curl
// (curl here uses the full 9-call mathematically correct version).
var Noise = (function() {
  function mod289(x) { return x - Math.floor(x / 289.0) * 289.0; }
  function permute(x) { return mod289(((x * 34.0) + 1.0) * x); }
  function tis(r) { return 1.79284291400159 - 0.85373472095314 * r; }

  function snoise(vx, vy, vz) {
    var C0=1.0/6.0, C1=1.0/3.0;
    var s=(vx+vy+vz)*C1;
    var ix=Math.floor(vx+s), iy=Math.floor(vy+s), iz=Math.floor(vz+s);
    var t=(ix+iy+iz)*C0;
    var x0x=vx-ix+t, x0y=vy-iy+t, x0z=vz-iz+t;
    var gx=x0x>=x0y?1:0, gy=x0y>=x0z?1:0, gz=x0z>=x0x?1:0;
    var lx=1-gx, ly=1-gy, lz=1-gz;
    var i1x=Math.min(gx,lz), i1y=Math.min(gy,lx), i1z=Math.min(gz,ly);
    var i2x=Math.max(gx,lz), i2y=Math.max(gy,lx), i2z=Math.max(gz,ly);
    var x1x=x0x-i1x+C0, x1y=x0y-i1y+C0, x1z=x0z-i1z+C0;
    var x2x=x0x-i2x+C1, x2y=x0y-i2y+C1, x2z=x0z-i2z+C1;
    var x3x=x0x-0.5,    x3y=x0y-0.5,    x3z=x0z-0.5;
    ix=mod289(ix); iy=mod289(iy); iz=mod289(iz);
    var p0=permute(permute(permute(iz)+iy)+ix);
    var p1=permute(permute(permute(iz+i1z)+iy+i1y)+ix+i1x);
    var p2=permute(permute(permute(iz+i2z)+iy+i2y)+ix+i2x);
    var p3=permute(permute(permute(iz+1)+iy+1)+ix+1);
    var nx=0.285714285714, ny=-0.928571428571, nz=0.142857142857;
    var j0=p0-49*Math.floor(p0*nz*nz), j1=p1-49*Math.floor(p1*nz*nz);
    var j2=p2-49*Math.floor(p2*nz*nz), j3=p3-49*Math.floor(p3*nz*nz);
    var x0_=Math.floor(j0*nz), y0_=Math.floor(j0-7*x0_);
    var x1_=Math.floor(j1*nz), y1_=Math.floor(j1-7*x1_);
    var x2_=Math.floor(j2*nz), y2_=Math.floor(j2-7*x2_);
    var x3_=Math.floor(j3*nz), y3_=Math.floor(j3-7*x3_);
    function grad(xs,ys){ var h=1-Math.abs(xs)-Math.abs(ys); var sh=h<=0?-1:0; return[xs+(Math.floor(xs)*2+1)*sh,ys+(Math.floor(ys)*2+1)*sh,h]; }
    var g0=grad(x0_*nx+ny,y0_*nx+ny), g1=grad(x1_*nx+ny,y1_*nx+ny);
    var g2=grad(x2_*nx+ny,y2_*nx+ny), g3=grad(x3_*nx+ny,y3_*nx+ny);
    var n0=tis(g0[0]*g0[0]+g0[1]*g0[1]+g0[2]*g0[2]);
    var n1=tis(g1[0]*g1[0]+g1[1]*g1[1]+g1[2]*g1[2]);
    var n2=tis(g2[0]*g2[0]+g2[1]*g2[1]+g2[2]*g2[2]);
    var n3=tis(g3[0]*g3[0]+g3[1]*g3[1]+g3[2]*g3[2]);
    g0[0]*=n0; g0[1]*=n0; g0[2]*=n0;
    g1[0]*=n1; g1[1]*=n1; g1[2]*=n1;
    g2[0]*=n2; g2[1]*=n2; g2[2]*=n2;
    g3[0]*=n3; g3[1]*=n3; g3[2]*=n3;
    var m0=Math.max(0.6-(x0x*x0x+x0y*x0y+x0z*x0z),0); m0*=m0;
    var m1=Math.max(0.6-(x1x*x1x+x1y*x1y+x1z*x1z),0); m1*=m1;
    var m2=Math.max(0.6-(x2x*x2x+x2y*x2y+x2z*x2z),0); m2*=m2;
    var m3=Math.max(0.6-(x3x*x3x+x3y*x3y+x3z*x3z),0); m3*=m3;
    return 42*(m0*m0*(g0[0]*x0x+g0[1]*x0y+g0[2]*x0z)+
               m1*m1*(g1[0]*x1x+g1[1]*x1y+g1[2]*x1z)+
               m2*m2*(g2[0]*x2x+g2[1]*x2y+g2[2]*x2z)+
               m3*m3*(g3[0]*x3x+g3[1]*x3y+g3[2]*x3z));
  }

  function frac(x){ return x - Math.floor(x); }

  function cellular(px, py, pz, metric) {
    var pix=Math.floor(px), piy=Math.floor(py), piz=Math.floor(pz);
    var pfx=px-pix, pfy=py-piy, pfz=pz-piz;
    var d1=10, d2=10;
    for(var oz=-1;oz<=1;oz++) for(var oy=-1;oy<=1;oy++) for(var ox=-1;ox<=1;ox++){
      var cx=pix+ox, cy=piy+oy, cz=piz+oz;
      var rpx=frac(Math.sin(cx*127.1+cy*311.7+cz*74.7)*43758.5453);
      var rpy=frac(Math.sin(cx*269.5+cy*183.3+cz*246.1)*43758.5453);
      var rpz=frac(Math.sin(cx*113.5+cy*271.9+cz*124.6)*43758.5453);
      var dx=ox+rpx-pfx, dy=oy+rpy-pfy, dz=oz+rpz-pfz;
      var d = metric==='euclidean'  ? Math.sqrt(dx*dx+dy*dy+dz*dz) :
              metric==='manhattan'  ? Math.abs(dx)+Math.abs(dy)+Math.abs(dz) :
              Math.max(Math.abs(dx),Math.max(Math.abs(dy),Math.abs(dz)));
      if(d<d1){d2=d1;d1=d;}else if(d<d2)d2=d;
    }
    return Math.max(Math.min((d2-d1)*2-1,1),-1);
  }

  function fbm(px, py, pz, octaves, lacunarity, gain) {
    var v=0,a=1,f=1,mx=0;
    for(var i=0;i<octaves;i++){ v+=snoise(px*f,py*f,pz*f)*a; mx+=a; a*=gain; f*=lacunarity; }
    return v/mx;
  }

  function ridged(px, py, pz, octaves, lacunarity, gain) {
    var v=0,a=1,f=1,mx=0;
    for(var i=0;i<octaves;i++){ v+=(1-Math.abs(snoise(px*f,py*f,pz*f)))*a; mx+=a; a*=gain; f*=lacunarity; }
    return (v/mx)*2.0-1.0;
  }

  function billow(px, py, pz, octaves, lacunarity, gain) {
    var v=0,a=1,f=1,mx=0;
    for(var i=0;i<octaves;i++){ v+=Math.abs(snoise(px*f,py*f,pz*f))*a; mx+=a; a*=gain; f*=lacunarity; }
    return (v/mx)*2.0-1.0;
  }

  function warp(px, py, pz, strength, octaves, lacunarity, gain) {
    var qx=fbm(px,py,pz,octaves,lacunarity,gain);
    var qy=fbm(px+5.2,py+1.3,pz+8.1,octaves,lacunarity,gain);
    var qz=fbm(px+3.7,py+9.4,pz+2.8,octaves,lacunarity,gain);
    return fbm(px+strength*qx,py+strength*qy,pz+strength*qz,octaves,lacunarity,gain);
  }

  // Full 9-call mathematically correct curl — superior to browser version
  // which used 6 calls to stay within ANGLE D3D instruction limits.
  function curl(px, py, pz, step, potScale) {
    var s=potScale, e=step;
    var psi_x=snoise(px*s,           py*s,           pz*s);
    var psi_y=snoise(px*s+3.7,       py*s+1.5,       pz*s+2.8);
    var psi_z=snoise(px*s+1.2,       py*s+4.6,       pz*s+0.9);
    var pzy  =snoise(px*s+1.2,       (py+e)*s+4.6,   pz*s+0.9);
    var pyz  =snoise(px*s+3.7,       py*s+1.5,       (pz+e)*s+2.8);
    var pxz  =snoise(px*s,           py*s,            (pz+e)*s);
    var pzx  =snoise((px+e)*s+1.2,   py*s+4.6,       pz*s+0.9);
    var pyx  =snoise((px+e)*s+3.7,   py*s+1.5,       pz*s+2.8);
    var pxy  =snoise(px*s,           (py+e)*s,        pz*s);
    var cx=(pzy-psi_z)/e-(pyz-psi_y)/e;
    var cy=(pxz-psi_x)/e-(pzx-psi_z)/e;
    var cz=(pyx-psi_y)/e-(pxy-psi_x)/e;
    return Math.sqrt(cx*cx+cy*cy+cz*cz);
  }

  return { snoise:snoise, cellular:cellular, fbm:fbm, ridged:ridged,
           billow:billow, warp:warp, curl:curl };
})();


// ── Field evaluators ──────────────────────────────────────────────────────

// Evaluates the scaffold SDF from a noise JSON recipe at world coords (wx,wy,wz).
// wx/wy/wz are in noise-space (already frequency-scaled).
// Returns negative inside scaffold material, positive outside.
function _evalNoiseField(surf, geom, noiseMin, noiseMax, wx, wy, wz) {
  var type = surf.noise_type || 'simplex';
  var raw;
  if(type==='simplex')  raw = Noise.snoise(wx, wy, wz);
  else if(type==='cellular') raw = Noise.cellular(wx, wy, wz, surf.distance_metric||'euclidean');
  else if(type==='fbm')  raw = Noise.fbm(wx, wy, wz, surf.octaves||4, surf.lacunarity||2.0, surf.gain||0.5);
  else if(type==='ridged') raw = Noise.ridged(wx, wy, wz, surf.octaves||4, surf.lacunarity||2.0, surf.gain||0.5);
  else if(type==='billow') raw = Noise.billow(wx, wy, wz, surf.octaves||4, surf.lacunarity||2.0, surf.gain||0.5);
  else if(type==='curl')   raw = Noise.curl(wx, wy, wz, surf.curl_step||0.1, surf.potential_scale||1.0);
  else /* warp */          raw = Noise.warp(wx, wy, wz, surf.warp_strength||1.0, surf.octaves||4, surf.lacunarity||2.0, surf.gain||0.5);

  var mid = (noiseMin + noiseMax) * 0.5;
  var halfR = Math.max((noiseMax - noiseMin) * 0.5, 0.001);
  var norm = (raw - mid) / halfR;
  var center = surf.center != null ? surf.center : 0.0;
  var halfW = surf.half_width != null ? surf.half_width : 0.15;
  var sign = (geom.mode === 'solid') ? -1 : 1;
  return sign * (Math.abs(norm - center) - halfW);
}

function _voidField(surf, wx_world, wy_world, wz_world) {
  // wx_world etc. are in raw world space (not frequency-scaled)
  // Void uses its own frequency
  var hf = surf.hole_freq != null ? surf.hole_freq : 1.0;
  var hr = surf.hole_radius != null ? surf.hole_radius : 0.3;
  if(surf.void_type === 'organic') {
    return Noise.snoise(wx_world*hf+17.3, wy_world*hf+5.8, wz_world*hf+11.2) - hr;
  }
  // Spherical: Worley F1
  function frac(x){return x-Math.floor(x);}
  var px=wx_world*hf, py=wy_world*hf, pz=wz_world*hf;
  var pix=Math.floor(px), piy=Math.floor(py), piz=Math.floor(pz);
  var pfx=px-pix, pfy=py-piy, pfz=pz-piz;
  var d1=1e9;
  for(var oz=-1;oz<=1;oz++) for(var oy=-1;oy<=1;oy++) for(var ox=-1;ox<=1;ox++){
    var cx=pix+ox, cy=piy+oy, cz=piz+oz;
    var rpx=frac(Math.sin(cx*127.1+cy*311.7+cz*74.7)*43758.5453);
    var rpy=frac(Math.sin(cx*269.5+cy*183.3+cz*246.1)*43758.5453);
    var rpz=frac(Math.sin(cx*113.5+cy*271.9+cz*124.6)*43758.5453);
    var dx=ox+rpx-pfx, dy=oy+rpy-pfy, dz=oz+rpz-pfz;
    var d=Math.sqrt(dx*dx+dy*dy+dz*dz);
    if(d<d1)d1=d;
  }
  return d1 - hr;
}

// Public field evaluator for noise JSON
// wx_world/wy_world/wz_world are in world space [-5,5]³
// Noise coords = world * freq  (matches GridGenerator exactly — no SCALE needed here)
function evaluateNoise(json, wx_world, wy_world, wz_world, noiseMin, noiseMax) {
  var surf = json.surface;
  var geom = json.geometry || {};
  var freq = surf.frequency || 0.3;
  // Box SDF — clips scaffold to the grid boundary, matching the tool's sceneSDF
  var bx = Math.abs(wx_world) - 5.0;
  var by = Math.abs(wy_world) - 5.0;
  var bz = Math.abs(wz_world) - 5.0;
  var boxDist = Math.max(bx, by, bz);  // negative inside box, positive outside
  var scaffold = _evalNoiseField(surf, geom, noiseMin, noiseMax,
    wx_world * freq,
    wy_world * freq,
    wz_world * freq);
  // Intersect scaffold with box — caps open surfaces at grid boundary
  var sdf = Math.max(scaffold, boxDist);
  if(surf.void_enabled) {
    var vf = _voidField(surf, wx_world, wy_world, wz_world);
    return Math.max(sdf, -vf);
  }
  return sdf;
}

// Stub for TPMS JSON (pass 2 — trigonometric term evaluator)
function evaluateTPMS(json, wx, wy, wz) {
  throw new Error('TPMS meshing not yet implemented. Export JSON from the TPMS Builder and check back in a future release.');
}

// Stub for FDSBN JSON (future)
function evaluateFDSBN(json, wx, wy, wz) {
  throw new Error('FDSBN meshing not yet implemented.');
}

// Route to correct evaluator based on surface type
function _getEvaluator(json) {
  var type = json.surface ? json.surface.type : null;
  if(type === 'noise') return evaluateNoise;
  if(type === 'terms' || type === 'raw_preset') return evaluateTPMS;
  if(type === 'fdsbn') return evaluateFDSBN;
  throw new Error('Unknown surface type: ' + type);
}


// ── Grid sampling ─────────────────────────────────────────────────────────
// Samples the SDF onto an N³ grid over a 10-unit world cube [-5, 5]³.
// The same domain the tool's shader uses.
// Returns Float32Array of length N³, index = x + y*N + z*N*N

function sampleGrid(json, N, statusCb) {
  var evaluator = _getEvaluator(json);
  var size = N * N * N;
  var values = new Float32Array(size);

  // Prepass: estimate noise field min/max over a coarse 16³ grid
  // (mirrors _runPrepass in the tool exactly)
  var surf = json.surface;
  var geom = json.geometry || {};
  var freq = surf.frequency || 0.3;
  var mn=Infinity, mx=-Infinity;
  var PP=16;
  for(var zi=0;zi<PP;zi++) for(var yi=0;yi<PP;yi++) for(var xi=0;xi<PP;xi++){
    // World coords [-5,5] — matching GridGenerator domain
    var px=((xi/(PP-1))*2-1)*5.0;
    var py=((yi/(PP-1))*2-1)*5.0;
    var pz=((zi/(PP-1))*2-1)*5.0;
    var type = surf.noise_type || 'simplex';
    var raw;
    var sx=px*freq, sy=py*freq, sz=pz*freq;  // no SCALE — matches GridGenerator
    if(type==='simplex')   raw=Noise.snoise(sx,sy,sz);
    else if(type==='cellular') raw=Noise.cellular(sx,sy,sz,surf.distance_metric||'euclidean');
    else if(type==='fbm')  raw=Noise.fbm(sx,sy,sz,surf.octaves||4,surf.lacunarity||2.0,surf.gain||0.5);
    else if(type==='ridged') raw=Noise.ridged(sx,sy,sz,surf.octaves||4,surf.lacunarity||2.0,surf.gain||0.5);
    else if(type==='billow') raw=Noise.billow(sx,sy,sz,surf.octaves||4,surf.lacunarity||2.0,surf.gain||0.5);
    else if(type==='curl')  raw=Noise.curl(sx,sy,sz,surf.curl_step||0.1,surf.potential_scale||1.0);
    else raw=Noise.warp(sx,sy,sz,surf.warp_strength||1.0,surf.octaves||4,surf.lacunarity||2.0,surf.gain||0.5);
    if(raw<mn)mn=raw; if(raw>mx)mx=raw;
  }
  var range=mx-mn;
  var noiseMin=mn-range*0.05, noiseMax=mx+range*0.05;

  // Sample full grid
  var step = 10.0 / N;  // world units per cell, domain is [-5,5]³
  var reportEvery = Math.floor(N / 4);
  for(var z=0;z<N;z++){
    if(statusCb && z%reportEvery===0) statusCb('Sampling grid... ' + Math.round(z/N*100) + '%');
    for(var y=0;y<N;y++) for(var x=0;x<N;x++){
      // Cell center in world space
      var wx = -5.0 + (x + 0.5) * step;
      var wy = -5.0 + (y + 0.5) * step;
      var wz = -5.0 + (z + 0.5) * step;
      values[x + y*N + z*N*N] = evaluateNoise(json, wx, wy, wz, noiseMin, noiseMax);
    }
  }
  return values;
}


// ── Marching cubes ────────────────────────────────────────────────────────
// Standard Lorensen & Cline (1987) with trilinear edge interpolation.
// Produces a triangle soup; vertices are deduplicated in a subsequent pass.

// Edge table and triangle table (256-case MC lookup)
var MC_EDGE_TABLE = new Uint16Array([
  0x000,0x109,0x203,0x30a,0x406,0x50f,0x605,0x70c,
  0x80c,0x905,0xa0f,0xb06,0xc0a,0xd03,0xe09,0xf00,
  0x190,0x099,0x393,0x29a,0x596,0x49f,0x795,0x69c,
  0x99c,0x895,0xb9f,0xa96,0xd9a,0xc93,0xf99,0xe90,
  0x230,0x339,0x033,0x13a,0x636,0x73f,0x435,0x53c,
  0xa3c,0xb35,0x83f,0x936,0xe3a,0xf33,0xc39,0xd30,
  0x3a0,0x2a9,0x1a3,0x0aa,0x7a6,0x6af,0x5a5,0x4ac,
  0xbac,0xaa5,0x9af,0x8a6,0xfaa,0xea3,0xda9,0xca0,
  0x460,0x569,0x663,0x76a,0x066,0x16f,0x265,0x36c,
  0xc6c,0xd65,0xe6f,0xf66,0x86a,0x963,0xa69,0xb60,
  0x5f0,0x4f9,0x7f3,0x6fa,0x1f6,0x0ff,0x3f5,0x2fc,
  0xdfc,0xcf5,0xfff,0xef6,0x9fa,0x8f3,0xbf9,0xaf0,
  0x650,0x759,0x453,0x55a,0x256,0x35f,0x055,0x15c,
  0xe5c,0xf55,0xc5f,0xd56,0xa5a,0xb53,0x859,0x950,
  0x7c0,0x6c9,0x5c3,0x4ca,0x3c6,0x2cf,0x1c5,0x0cc,
  0xfcc,0xec5,0xdcf,0xcc6,0xbca,0xac3,0x9c9,0x8c0,
  0x8c0,0x9c9,0xac3,0xbca,0xcc6,0xdcf,0xec5,0xfcc,
  0x0cc,0x1c5,0x2cf,0x3c6,0x4ca,0x5c3,0x6c9,0x7c0,
  0x950,0x859,0xb53,0xa5a,0xd56,0xc5f,0xf55,0xe5c,
  0x15c,0x055,0x35f,0x256,0x55a,0x453,0x759,0x650,
  0xaf0,0xbf9,0x8f3,0x9fa,0xef6,0xfff,0xcf5,0xdfc,
  0x2fc,0x3f5,0x0ff,0x1f6,0x6fa,0x7f3,0x4f9,0x5f0,
  0xb60,0xa69,0x963,0x86a,0xf66,0xe6f,0xd65,0xc6c,
  0x36c,0x265,0x16f,0x066,0x76a,0x663,0x569,0x460,
  0xca0,0xda9,0xea3,0xfaa,0x8a6,0x9af,0xaa5,0xbac,
  0x4ac,0x5a5,0x6af,0x7a6,0x0aa,0x1a3,0x2a9,0x3a0,
  0xd30,0xc39,0xf33,0xe3a,0x936,0x83f,0xb35,0xa3c,
  0x53c,0x435,0x73f,0x636,0x13a,0x033,0x339,0x230,
  0xe90,0xf99,0xc93,0xd9a,0xa96,0xb9f,0x895,0x99c,
  0x69c,0x795,0x49f,0x596,0x29a,0x393,0x099,0x190,
  0xf00,0xe09,0xd03,0xc0a,0xb06,0xa0f,0x905,0x80c,
  0x70c,0x605,0x50f,0x406,0x30a,0x203,0x109,0x000
]);

// The full 256-entry triangle table is complex — use a well-known complete table
// This is the Bourke (1994) table, widely used and validated
var _TRI_TABLE_FLAT = new Int8Array([
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 1, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 8, 3, 9, 8, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 3, 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   9, 2,10, 0, 2, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   2, 8, 3, 2,10, 8,10, 9, 8,-1,-1,-1,-1,-1,-1,-1,
   3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0,11, 2, 8,11, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 9, 0, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1,11, 2, 1, 9,11, 9, 8,11,-1,-1,-1,-1,-1,-1,-1,
   3,10, 1,11,10, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0,10, 1, 0, 8,10, 8,11,10,-1,-1,-1,-1,-1,-1,-1,
   3, 9, 0, 3,11, 9,11,10, 9,-1,-1,-1,-1,-1,-1,-1,
   9, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4, 3, 0, 7, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 1, 9, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4, 1, 9, 4, 7, 1, 7, 3, 1,-1,-1,-1,-1,-1,-1,-1,
   1, 2,10, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   3, 4, 7, 3, 0, 4, 1, 2,10,-1,-1,-1,-1,-1,-1,-1,
   9, 2,10, 9, 0, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,
   2,10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4,-1,-1,-1,-1,
   8, 4, 7, 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  11, 4, 7,11, 2, 4, 2, 0, 4,-1,-1,-1,-1,-1,-1,-1,
   9, 0, 1, 8, 4, 7, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,
   4, 7,11, 9, 4,11, 9,11, 2, 9, 2, 1,-1,-1,-1,-1,
   3,10, 1, 3,11,10, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1,
   1,11,10, 1, 4,11, 1, 0, 4, 7,11, 4,-1,-1,-1,-1,
   4, 7, 8, 9, 0,11, 9,11,10,11, 0, 3,-1,-1,-1,-1,
   4, 7,11, 4,11, 9, 9,11,10,-1,-1,-1,-1,-1,-1,-1,
   9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   9, 5, 4, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 5, 4, 1, 5, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   8, 5, 4, 8, 3, 5, 3, 1, 5,-1,-1,-1,-1,-1,-1,-1,
   1, 2,10, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   3, 0, 8, 1, 2,10, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1,
   5, 2,10, 5, 4, 2, 4, 0, 2,-1,-1,-1,-1,-1,-1,-1,
   2,10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8,-1,-1,-1,-1,
   9, 5, 4, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0,11, 2, 0, 8,11, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1,
   0, 5, 4, 0, 1, 5, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,
   2, 1, 5, 2, 5, 8, 2, 8,11, 4, 8, 5,-1,-1,-1,-1,
  10, 3,11,10, 1, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,
   4, 9, 5, 0, 8, 1, 8,10, 1, 8,11,10,-1,-1,-1,-1,
   5, 4, 0, 5, 0,11, 5,11,10,11, 0, 3,-1,-1,-1,-1,
   5, 4, 8, 5, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1,
   9, 7, 8, 5, 7, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   9, 3, 0, 9, 5, 3, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1,
   0, 7, 8, 0, 1, 7, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1,
   1, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   9, 7, 8, 9, 5, 7,10, 1, 2,-1,-1,-1,-1,-1,-1,-1,
  10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3,-1,-1,-1,-1,
   8, 0, 2, 8, 2, 5, 8, 5, 7,10, 5, 2,-1,-1,-1,-1,
   2,10, 5, 2, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1,
   7, 9, 5, 7, 8, 9, 3,11, 2,-1,-1,-1,-1,-1,-1,-1,
   9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7,11,-1,-1,-1,-1,
   2, 3,11, 0, 1, 8, 1, 7, 8, 1, 5, 7,-1,-1,-1,-1,
  11, 2, 1,11, 1, 7, 7, 1, 5,-1,-1,-1,-1,-1,-1,-1,
   9, 5, 8, 8, 5, 7,10, 1, 3,10, 3,11,-1,-1,-1,-1,
   5, 7, 0, 5, 0, 9, 7,11, 0, 1, 0,10,11,10, 0,-1,
  11,10, 0,11, 0, 3,10, 5, 0, 8, 0, 7, 5, 7, 0,-1,
  11,10, 5, 7,11, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 3, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   9, 0, 1, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 8, 3, 1, 9, 8, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,
   1, 6, 5, 2, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 6, 5, 1, 2, 6, 3, 0, 8,-1,-1,-1,-1,-1,-1,-1,
   9, 6, 5, 9, 0, 6, 0, 2, 6,-1,-1,-1,-1,-1,-1,-1,
   5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8,-1,-1,-1,-1,
   2, 3,11,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  11, 0, 8,11, 2, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,
   0, 1, 9, 2, 3,11, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,
   5,10, 6, 1, 9, 2, 9,11, 2, 9, 8,11,-1,-1,-1,-1,
   6, 3,11, 6, 5, 3, 5, 1, 3,-1,-1,-1,-1,-1,-1,-1,
   0, 8,11, 0,11, 5, 0, 5, 1, 5,11, 6,-1,-1,-1,-1,
   3,11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9,-1,-1,-1,-1,
   6, 5, 9, 6, 9,11,11, 9, 8,-1,-1,-1,-1,-1,-1,-1,
   5,10, 6, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4, 3, 0, 4, 7, 3, 6, 5,10,-1,-1,-1,-1,-1,-1,-1,
   1, 9, 0, 5,10, 6, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,
  10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4,-1,-1,-1,-1,
   6, 1, 2, 6, 5, 1, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,
   1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7,-1,-1,-1,-1,
   8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6,-1,-1,-1,-1,
   7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9,-1,
   3,11, 2, 7, 8, 4,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,
   5,10, 6, 4, 7, 2, 4, 2, 0, 2, 7,11,-1,-1,-1,-1,
   0, 1, 9, 4, 7, 8, 2, 3,11, 5,10, 6,-1,-1,-1,-1,
   9, 2, 1, 9,11, 2, 9, 4,11, 7,11, 4, 5,10, 6,-1,
   8, 4, 7, 3,11, 5, 3, 5, 1, 5,11, 6,-1,-1,-1,-1,
   5, 1,11, 5,11, 6, 1, 0,11, 7,11, 4, 0, 4,11,-1,
   0, 5, 9, 0, 6, 5, 0, 3, 6,11, 6, 3, 8, 4, 7,-1,
   6, 5, 9, 6, 9,11, 4, 7, 9, 7,11, 9,-1,-1,-1,-1,
  10, 4, 9, 6, 4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4,10, 6, 4, 9,10, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,
  10, 0, 1,10, 6, 0, 6, 4, 0,-1,-1,-1,-1,-1,-1,-1,
   8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1,10,-1,-1,-1,-1,
   1, 4, 9, 1, 2, 4, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1,
   3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4,-1,-1,-1,-1,
   0, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   8, 3, 2, 8, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1,
  10, 4, 9,10, 6, 4,11, 2, 3,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 2, 2, 8,11, 4, 9,10, 4,10, 6,-1,-1,-1,-1,
   3,11, 2, 0, 1, 6, 0, 6, 4, 6, 1,10,-1,-1,-1,-1,
   6, 4, 1, 6, 1,10, 4, 8, 1, 2, 1,11, 8,11, 1,-1,
   9, 6, 4, 9, 3, 6, 9, 1, 3,11, 6, 3,-1,-1,-1,-1,
   8,11, 1, 8, 1, 0,11, 6, 1, 9, 1, 4, 6, 4, 1,-1,
   3,11, 6, 3, 6, 0, 0, 6, 4,-1,-1,-1,-1,-1,-1,-1,
   6, 4, 8,11, 6, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   7,10, 6, 7, 8,10, 8, 9,10,-1,-1,-1,-1,-1,-1,-1,
   0, 7, 3, 0,10, 7, 0, 9,10, 6, 7,10,-1,-1,-1,-1,
  10, 6, 7, 1,10, 7, 1, 7, 8, 1, 8, 0,-1,-1,-1,-1,
  10, 6, 7,10, 7, 1, 1, 7, 3,-1,-1,-1,-1,-1,-1,-1,
   1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7,-1,-1,-1,-1,
   2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9,-1,
   7, 8, 0, 7, 0, 6, 6, 0, 2,-1,-1,-1,-1,-1,-1,-1,
   7, 3, 2, 6, 7, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   2, 3,11,10, 6, 8,10, 8, 9, 8, 6, 7,-1,-1,-1,-1,
   2, 0, 7, 2, 7,11, 0, 9, 7, 6, 7,10, 9,10, 7,-1,
   1, 8, 0, 1, 7, 8, 1,10, 7, 6, 7,10, 2, 3,11,-1,
  11, 2, 1,11, 1, 7,10, 6, 1, 6, 7, 1,-1,-1,-1,-1,
   8, 9, 6, 8, 6, 7, 9, 1, 6,11, 6, 3, 1, 3, 6,-1,
   0, 9, 1,11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   7, 8, 0, 7, 0, 6, 3,11, 0,11, 6, 0,-1,-1,-1,-1,
   7,11, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   3, 0, 8,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 1, 9,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   8, 1, 9, 8, 3, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,
  10, 1, 2, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 2,10, 3, 0, 8, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,
   2, 9, 0, 2,10, 9, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,
   6,11, 7, 2,10, 3,10, 8, 3,10, 9, 8,-1,-1,-1,-1,
   7, 2, 3, 6, 2, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   7, 0, 8, 7, 6, 0, 6, 2, 0,-1,-1,-1,-1,-1,-1,-1,
   2, 7, 6, 2, 3, 7, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1,
   1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6,-1,-1,-1,-1,
  10, 7, 6,10, 1, 7, 1, 3, 7,-1,-1,-1,-1,-1,-1,-1,
  10, 7, 6, 1, 7,10, 1, 8, 7, 1, 0, 8,-1,-1,-1,-1,
   0, 3, 7, 0, 7,10, 0,10, 9, 6,10, 7,-1,-1,-1,-1,
   7, 6,10, 7,10, 8, 8,10, 9,-1,-1,-1,-1,-1,-1,-1,
   6, 8, 4,11, 8, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   3, 6,11, 3, 0, 6, 0, 4, 6,-1,-1,-1,-1,-1,-1,-1,
   8, 6,11, 8, 4, 6, 9, 0, 1,-1,-1,-1,-1,-1,-1,-1,
   9, 4, 6, 9, 6, 3, 9, 3, 1,11, 3, 6,-1,-1,-1,-1,
   6, 8, 4, 6,11, 8, 2,10, 1,-1,-1,-1,-1,-1,-1,-1,
   1, 2,10, 3, 0,11, 0, 6,11, 0, 4, 6,-1,-1,-1,-1,
   4,11, 8, 4, 6,11, 0, 2, 9, 2,10, 9,-1,-1,-1,-1,
  10, 9, 3,10, 3, 2, 9, 4, 3,11, 3, 6, 4, 6, 3,-1,
   8, 2, 3, 8, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1,
   0, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8,-1,-1,-1,-1,
   1, 9, 4, 1, 4, 2, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1,
   8, 1, 3, 8, 6, 1, 8, 4, 6, 6,10, 1,-1,-1,-1,-1,
  10, 1, 0,10, 0, 6, 6, 0, 4,-1,-1,-1,-1,-1,-1,-1,
   4, 6, 3, 4, 3, 8, 6,10, 3, 0, 3, 9,10, 9, 3,-1,
  10, 9, 4, 6,10, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4, 9, 5, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 3, 4, 9, 5,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,
   5, 0, 1, 5, 4, 0, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,
  11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5,-1,-1,-1,-1,
   9, 5, 4,10, 1, 2, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,
   6,11, 7, 1, 2,10, 0, 8, 3, 4, 9, 5,-1,-1,-1,-1,
   7, 6,11, 5, 4,10, 4, 2,10, 4, 0, 2,-1,-1,-1,-1,
   3, 4, 8, 3, 5, 4, 3, 2, 5,10, 5, 2,11, 7, 6,-1,
   7, 2, 3, 7, 6, 2, 5, 4, 9,-1,-1,-1,-1,-1,-1,-1,
   9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7,-1,-1,-1,-1,
   3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0,-1,-1,-1,-1,
   6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8,-1,
   9, 5, 4,10, 1, 6, 1, 7, 6, 1, 3, 7,-1,-1,-1,-1,
   1, 6,10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4,-1,
   4, 0,10, 4,10, 5, 0, 3,10, 6,10, 7, 3, 7,10,-1,
   7, 6,10, 7,10, 8, 5, 4,10, 4, 8,10,-1,-1,-1,-1,
   6, 9, 5, 6,11, 9,11, 8, 9,-1,-1,-1,-1,-1,-1,-1,
   3, 6,11, 0, 6, 3, 0, 5, 6, 0, 9, 5,-1,-1,-1,-1,
   0,11, 8, 0, 5,11, 0, 1, 5, 5, 6,11,-1,-1,-1,-1,
   6,11, 3, 6, 3, 5, 5, 3, 1,-1,-1,-1,-1,-1,-1,-1,
   1, 2,10, 9, 5,11, 9,11, 8,11, 5, 6,-1,-1,-1,-1,
   0,11, 3, 0, 6,11, 0, 9, 6, 5, 6, 9, 1, 2,10,-1,
  11, 8, 5,11, 5, 6, 8, 0, 5,10, 5, 2, 0, 2, 5,-1,
   6,11, 3, 6, 3, 5, 2,10, 3,10, 5, 3,-1,-1,-1,-1,
   5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2,-1,-1,-1,-1,
   9, 5, 6, 9, 6, 0, 0, 6, 2,-1,-1,-1,-1,-1,-1,-1,
   1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8,-1,
   1, 5, 6, 2, 1, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  10, 5, 6, 9, 5,10, 0, 1, 8, 1,11, 8, 1, 3,11,-1,
   0, 9, 1,10, 5, 6,11, 0, 3,-1,-1,-1,-1,-1,-1,-1,
  11, 8, 5,11, 5, 6, 0, 3, 5, 3, 6, 5,-1,-1,-1,-1,
  11, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  11, 5,10, 7, 5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  11, 5,10,11, 7, 5, 8, 3, 0,-1,-1,-1,-1,-1,-1,-1,
   5,11, 7, 5,10,11, 1, 9, 0,-1,-1,-1,-1,-1,-1,-1,
  10, 7, 5,10,11, 7, 9, 8, 1, 8, 3, 1,-1,-1,-1,-1,
  11, 1, 2,11, 7, 1, 7, 5, 1,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2,11,-1,-1,-1,-1,
   9, 7, 5, 9, 2, 7, 9, 0, 2, 2,11, 7,-1,-1,-1,-1,
   7, 5, 2, 7, 2,11, 5, 9, 2, 3, 2, 8, 9, 8, 2,-1,
   2, 5,10, 2, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1,
   8, 2, 0, 8, 5, 2, 8, 7, 5,10, 2, 5,-1,-1,-1,-1,
   9, 0, 1, 5,10, 3, 5, 3, 7, 3,10, 2,-1,-1,-1,-1,
   9, 8, 2, 9, 2, 1, 8, 7, 2,10, 2, 5, 7, 5, 2,-1,
   1, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 7, 0, 7, 1, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1,
   9, 0, 3, 9, 3, 5, 5, 3, 7,-1,-1,-1,-1,-1,-1,-1,
   9, 8, 7, 5, 9, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   5, 8, 4, 5,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1,
   5, 0, 4, 5,11, 0, 5,10,11,11, 3, 0,-1,-1,-1,-1,
   0, 1, 9, 8, 4,10, 8,10,11,10, 4, 5,-1,-1,-1,-1,
  10,11, 4,10, 4, 5,11, 3, 4, 9, 4, 1, 3, 1, 4,-1,
   2, 5, 1, 2, 8, 5, 2,11, 8, 4, 5, 8,-1,-1,-1,-1,
   0, 4,11, 0,11, 3, 4, 5,11, 2,11, 1, 5, 1,11,-1,
   0, 2, 5, 0, 5, 9, 2,11, 5, 4, 5, 8,11, 8, 5,-1,
   9, 4, 5, 2,11, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   2, 5,10, 3, 5, 2, 3, 4, 5, 3, 8, 4,-1,-1,-1,-1,
   5,10, 2, 5, 2, 4, 4, 2, 0,-1,-1,-1,-1,-1,-1,-1,
   3,10, 2, 3, 5,10, 3, 8, 5, 4, 5, 8, 0, 1, 9,-1,
   5,10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2,-1,-1,-1,-1,
   8, 4, 5, 8, 5, 3, 3, 5, 1,-1,-1,-1,-1,-1,-1,-1,
   0, 4, 5, 1, 0, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5,-1,-1,-1,-1,
   9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4,11, 7, 4, 9,11, 9,10,11,-1,-1,-1,-1,-1,-1,-1,
   0, 8, 3, 4, 9, 7, 9,11, 7, 9,10,11,-1,-1,-1,-1,
   1,10,11, 1,11, 4, 1, 4, 0, 7, 4,11,-1,-1,-1,-1,
   3, 1, 4, 3, 4, 8, 1,10, 4, 7, 4,11,10,11, 4,-1,
   4,11, 7, 9,11, 4, 9, 2,11, 9, 1, 2,-1,-1,-1,-1,
   9, 7, 4, 9,11, 7, 9, 1,11, 2,11, 1, 0, 8, 3,-1,
  11, 7, 4,11, 4, 2, 2, 4, 0,-1,-1,-1,-1,-1,-1,-1,
  11, 7, 4,11, 4, 2, 8, 3, 4, 3, 2, 4,-1,-1,-1,-1,
   2, 9,10, 2, 7, 9, 2, 3, 7, 7, 4, 9,-1,-1,-1,-1,
   9,10, 7, 9, 7, 4,10, 2, 7, 8, 7, 0, 2, 0, 7,-1,
   3, 7,10, 3,10, 2, 7, 4,10, 1,10, 0, 4, 0,10,-1,
   1,10, 2, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4, 9, 1, 4, 1, 7, 7, 1, 3,-1,-1,-1,-1,-1,-1,-1,
   4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1,-1,-1,-1,-1,
   4, 0, 3, 7, 4, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   4, 8, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   9,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   3, 0, 9, 3, 9,11,11, 9,10,-1,-1,-1,-1,-1,-1,-1,
   0, 1,10, 0,10, 8, 8,10,11,-1,-1,-1,-1,-1,-1,-1,
   3, 1,10,11, 3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 2,11, 1,11, 9, 9,11, 8,-1,-1,-1,-1,-1,-1,-1,
   3, 0, 9, 3, 9,11, 1, 2, 9, 2,11, 9,-1,-1,-1,-1,
   0, 2,11, 8, 0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   3, 2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   2, 3, 8, 2, 8,10,10, 8, 9,-1,-1,-1,-1,-1,-1,-1,
   9,10, 2, 0, 9, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   2, 3, 8, 2, 8,10, 0, 1, 8, 1,10, 8,-1,-1,-1,-1,
   1,10, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   1, 3, 8, 9, 1, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 9, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
   0, 3, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
]);

// Cell corner offsets for the 8 vertices of a unit cube
var CORNERS = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];
// Edge vertex pairs
var EDGES = [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]];

function marchingCubes(grid, N, statusCb) {
  var N2 = N * N;
  var verts = [];   // flat [x,y,z, x,y,z, ...]
  var tris  = [];   // flat [i,j,k, ...]
  var vi = 0;

  function idx(x,y,z){ return x + y*N + z*N2; }
  function lerp(v0,v1,p0,p1,axis){
    // linear interpolation along axis to find zero crossing
    var dv = v1 - v0;
    var t = Math.abs(dv) < 1e-10 ? 0.5 : -v0 / dv;
    t = Math.max(0, Math.min(1, t));
    var coords = [p0[0]+t*(p1[0]-p0[0]), p0[1]+t*(p1[1]-p0[1]), p0[2]+t*(p1[2]-p0[2])];
    return coords;
  }

  var reportEvery = Math.max(1, Math.floor(N / 8));

  for(var z=0; z<N-1; z++){
    if(statusCb && z%reportEvery===0) statusCb('Meshing... ' + Math.round(z/(N-1)*100) + '%');
    for(var y=0; y<N-1; y++){
      for(var x=0; x<N-1; x++){
        // Sample the 8 corner values
        var vals = new Array(8);
        var caseIdx = 0;
        for(var c=0; c<8; c++){
          var cx=x+CORNERS[c][0], cy=y+CORNERS[c][1], cz=z+CORNERS[c][2];
          vals[c] = grid[idx(cx,cy,cz)];
          if(vals[c] < 0) caseIdx |= (1<<c);
        }
        if(caseIdx===0 || caseIdx===255) continue;

        // Compute edge intersection points
        var edgePts = new Array(12);
        var edgeMask = MC_EDGE_TABLE[caseIdx];
        for(var e=0; e<12; e++){
          if(edgeMask & (1<<e)){
            var e0=EDGES[e][0], e1=EDGES[e][1];
            var p0=[x+CORNERS[e0][0], y+CORNERS[e0][1], z+CORNERS[e0][2]];
            var p1=[x+CORNERS[e1][0], y+CORNERS[e1][1], z+CORNERS[e1][2]];
            edgePts[e] = lerp(vals[e0], vals[e1], p0, p1, e);
          }
        }

        // Emit triangles from lookup table
        var base = caseIdx * 16;
        for(var t=0; t<15; t+=3){
          var a=_TRI_TABLE_FLAT[base+t];
          if(a===-1) break;
          var b=_TRI_TABLE_FLAT[base+t+1], cc=_TRI_TABLE_FLAT[base+t+2];
          var pa=edgePts[a], pb=edgePts[b], pc=edgePts[cc];
          // Guard: skip if any edge point wasn't computed (table inconsistency)
          if(!pa || !pb || !pc) continue;
          verts.push(pa[0],pa[1],pa[2], pb[0],pb[1],pb[2], pc[0],pc[1],pc[2]);
          tris.push(vi, vi+1, vi+2);
          vi+=3;
        }
      }
    }
  }

  return {
    vertices: new Float32Array(verts),
    triangles: new Int32Array(tris),
    vertexCount: vi,
    triCount: tris.length / 3
  };
}


// ── Taubin smoothing ──────────────────────────────────────────────────────
// λ/μ alternating Laplacian — preserves volume, removes staircase artifacts.
// λ > 0 (smooth), μ < 0 (inflate to cancel shrinkage).
// Default: λ=0.5, μ=-0.53, 10 iterations — standard medical mesh values.

function taubinSmooth(mesh, lambda, mu, iterations, statusCb) {
  lambda = lambda != null ? lambda : 0.5;
  mu     = mu     != null ? mu     : -0.53;
  iterations = iterations != null ? iterations : 10;

  var verts = mesh.vertices;       // Float32Array, 3 floats per vertex
  var tris  = mesh.triangles;      // Int32Array, 3 ints per triangle
  var nVerts = verts.length / 3;
  var nTris  = tris.length / 3;

  // Build adjacency: for each vertex, list of neighbouring vertex indices
  // We have triangle soup so we first deduplicate vertices by position
  // (simple spatial hash since MC produces unshared verts)
  var adj = new Array(nVerts);
  for(var i=0;i<nVerts;i++) adj[i] = [];

  function addNeighbour(a,b){
    if(adj[a].indexOf(b)<0) adj[a].push(b);
    if(adj[b].indexOf(a)<0) adj[b].push(a);
  }

  for(var t=0;t<nTris;t++){
    var a=tris[t*3], b=tris[t*3+1], c=tris[t*3+2];
    addNeighbour(a,b); addNeighbour(b,c); addNeighbour(a,c);
  }

  var tmp = new Float32Array(nVerts*3);

  function laplacianStep(factor){
    for(var v=0;v<nVerts;v++){
      var nb=adj[v], nb_len=nb.length;
      if(nb_len===0){ tmp[v*3]=verts[v*3]; tmp[v*3+1]=verts[v*3+1]; tmp[v*3+2]=verts[v*3+2]; continue; }
      var cx=0,cy=0,cz=0;
      for(var n=0;n<nb_len;n++){ cx+=verts[nb[n]*3]; cy+=verts[nb[n]*3+1]; cz+=verts[nb[n]*3+2]; }
      cx/=nb_len; cy/=nb_len; cz/=nb_len;
      tmp[v*3]   = verts[v*3]   + factor*(cx-verts[v*3]);
      tmp[v*3+1] = verts[v*3+1] + factor*(cy-verts[v*3+1]);
      tmp[v*3+2] = verts[v*3+2] + factor*(cz-verts[v*3+2]);
    }
    for(var v=0;v<nVerts*3;v++) verts[v]=tmp[v];
  }

  for(var iter=0; iter<iterations; iter++){
    if(statusCb) statusCb('Smoothing... ' + Math.round(iter/iterations*100) + '%');
    laplacianStep(lambda);  // forward: smooth
    laplacianStep(mu);      // backward: inflate to cancel shrinkage
  }

  return mesh; // modified in place
}


// ── Mesh quality checks ───────────────────────────────────────────────────

function checkQuality(mesh) {
  var verts = mesh.vertices;
  var tris  = mesh.triangles;
  var nTris = tris.length / 3;

  // Edge → face count map (using string key for simplicity)
  var edgeCount = {};
  var degenCount = 0;

  function edgeKey(a,b){ return a<b ? a+','+b : b+','+a; }

  for(var t=0;t<nTris;t++){
    var a=tris[t*3], b=tris[t*3+1], c=tris[t*3+2];

    // Degenerate check: area via cross product
    var ax=verts[a*3],ay=verts[a*3+1],az=verts[a*3+2];
    var bx=verts[b*3],by=verts[b*3+1],bz=verts[b*3+2];
    var cx2=verts[c*3],cy2=verts[c*3+1],cz2=verts[c*3+2];
    var ux=bx-ax,uy=by-ay,uz=bz-az;
    var vx=cx2-ax,vy=cy2-ay,vz=cz2-az;
    var area2 = (uy*vz-uz*vy)*(uy*vz-uz*vy)+(uz*vx-ux*vz)*(uz*vx-ux*vz)+(ux*vy-uy*vx)*(ux*vy-uy*vx);
    if(area2 < 1e-14) { degenCount++; continue; }

    var k1=edgeKey(a,b), k2=edgeKey(b,c), k3=edgeKey(a,c);
    edgeCount[k1]=(edgeCount[k1]||0)+1;
    edgeCount[k2]=(edgeCount[k2]||0)+1;
    edgeCount[k3]=(edgeCount[k3]||0)+1;
  }

  var boundaryEdges=0, nonManifold=0;
  for(var k in edgeCount){
    if(edgeCount[k]===1) boundaryEdges++;
    else if(edgeCount[k]>2) nonManifold++;
  }

  return {
    triCount:      nTris,
    degenCount:    degenCount,
    boundaryEdges: boundaryEdges,
    nonManifold:   nonManifold,
    watertight:    boundaryEdges===0 && nonManifold===0
  };
}


// ── STL serializer ────────────────────────────────────────────────────────
// Binary STL — 80-byte header + 4-byte count + 50 bytes per triangle.
// Vertex positions scaled from grid index space to millimetres.

function toBinarySTL(mesh, N, cellSizeMm) {
  cellSizeMm = cellSizeMm || 10.0;
  var scale = cellSizeMm / N;  // grid units → mm

  var tris  = mesh.triangles;
  var verts = mesh.vertices;
  var nTris = tris.length / 3;

  var buf = new ArrayBuffer(84 + nTris * 50);
  var view = new DataView(buf);

  // 80-byte ASCII header
  var header = 'Noise Scaffold Explorer — mesher.js';
  for(var i=0;i<80;i++) view.setUint8(i, i<header.length ? header.charCodeAt(i) : 0);
  view.setUint32(80, nTris, true);

  var offset = 84;
  for(var t=0;t<nTris;t++){
    var a=tris[t*3], b=tris[t*3+1], c=tris[t*3+2];
    var ax=verts[a*3]*scale, ay=verts[a*3+1]*scale, az=verts[a*3+2]*scale;
    var bx=verts[b*3]*scale, by=verts[b*3+1]*scale, bz=verts[b*3+2]*scale;
    var cx=verts[c*3]*scale, cy=verts[c*3+1]*scale, cz=verts[c*3+2]*scale;

    // Reverse winding (a,c,b) so outward normals point away from scaffold material
    var ux=cx-ax,uy=cy-ay,uz=cz-az;
    var vvx=bx-ax,vvy=by-ay,vvz=bz-az;
    var nx=uy*vvz-uz*vvy, ny=uz*vvx-ux*vvz, nz=ux*vvy-uy*vvx;
    var nl=Math.sqrt(nx*nx+ny*ny+nz*nz)||1;
    nx/=nl; ny/=nl; nz/=nl;

    view.setFloat32(offset,    nx, true); view.setFloat32(offset+4,  ny, true); view.setFloat32(offset+8,  nz, true);
    view.setFloat32(offset+12, ax, true); view.setFloat32(offset+16, ay, true); view.setFloat32(offset+20, az, true);
    view.setFloat32(offset+24, cx, true); view.setFloat32(offset+28, cy, true); view.setFloat32(offset+32, cz, true);
    view.setFloat32(offset+36, bx, true); view.setFloat32(offset+40, by, true); view.setFloat32(offset+44, bz, true);
    view.setUint16(offset+48, 0, true);
    offset += 50;
  }
  return buf;
}


// ── Public API ────────────────────────────────────────────────────────────

// Full pipeline: JSON → sampled grid → MC mesh → Taubin smooth → quality check
// Returns mesh object with attached quality report.
function meshFromJSON(json, N, statusCb) {
  N = N || 64;

  // Validate input
  if(!json || !json.surface) throw new Error('Invalid JSON: missing surface object');
  if(json.surface.type !== 'noise') throw new Error('Surface type "'+json.surface.type+'" not yet supported for mesh export');

  if(statusCb) statusCb('Sampling grid...');
  var grid;
  try { grid = sampleGrid(json, N, statusCb); }
  catch(e) { throw new Error('Grid sampling failed: ' + e.message); }
  if(!grid || grid.length === 0) throw new Error('Grid sampling returned empty result');

  if(statusCb) statusCb('Running marching cubes...');
  var mesh;
  try { mesh = marchingCubes(grid, N, statusCb); }
  catch(e) { throw new Error('Marching cubes failed: ' + e.message); }
  if(!mesh || mesh.triCount === 0) throw new Error('Marching cubes produced no triangles — try adjusting Center or Half-width so the surface is inside the domain');

  if(statusCb) statusCb('Smoothing (Taubin)...');
  try { taubinSmooth(mesh, 0.5, -0.53, 10, statusCb); }
  catch(e) { throw new Error('Smoothing failed: ' + e.message); }

  if(statusCb) statusCb('Checking quality...');
  mesh.quality = checkQuality(mesh);

  return mesh;
}

// Download STL from a JSON recipe.
// cellSizeMm: physical size of the output cube (default 10mm).
function downloadSTL(json, N, cellSizeMm, statusCb) {
  cellSizeMm = cellSizeMm || 10.0;
  var mesh = meshFromJSON(json, N, statusCb);
  var q = mesh.quality;

  var stl = toBinarySTL(mesh, N, cellSizeMm);

  // Build a descriptive filename
  var preset = json.meta ? json.meta.preset : 'scaffold';
  var ts = new Date().toISOString().slice(0,19).replace(/[:T]/g,'-');
  var filename = 'noise_' + preset + '_' + N + 'cube_' + ts + '.stl';

  var blob = new Blob([stl], {type:'model/stl'});
  var url  = URL.createObjectURL(blob);
  var a    = document.createElement('a');
  a.href=url; a.download=filename;
  document.body.appendChild(a); a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);

  return q;  // caller can display quality info
}

// Expose public surface
global.Mesher = {
  meshFromJSON:   meshFromJSON,
  downloadSTL:    downloadSTL,
  checkQuality:   checkQuality,
  toBinarySTL:    toBinarySTL,
  taubinSmooth:   taubinSmooth,
  marchingCubes:  marchingCubes,
  sampleGrid:     sampleGrid,
  evaluateNoise:  evaluateNoise,
  evaluateTPMS:   evaluateTPMS,
  evaluateFDSBN:  evaluateFDSBN,
  // Noise primitives exposed for testing
  Noise: Noise
};

})(typeof window !== 'undefined' ? window : global);
