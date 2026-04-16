# Noise Scaffold Explorer

**Part of the [F13LD](https://f13ld.app) tool suite by [Not a Robot Engineering](https://notarobot-eng.com)**

🔗 **[Launch → F13LD.noise](https://mshomper.github.io/f13ld.noise)**

A browser-based implicit field explorer for generating stochastic scaffold geometry from mathematical noise fields. Seven noise types share a common design pipeline: interactive 2D section view with X/Y/Z slicing, WebGL raymarched 3D preview, MIL-HS mechanical homogenization, STL mesh export, and JSON export for downstream validation.

All geometry is fully implicit — no meshes, no STL intermediates at the design stage. Structure is defined as a scalar field evaluated on demand at any point in space.

---

## Overview

Periodic lattices (TPMS, BCC, FCC) offer predictable but rigid geometry. Real biological structures — trabecular bone, open-cell cartilage, vascular networks — are stochastic: they achieve mechanical efficiency through statistical rather than crystallographic order. This tool lets you explore that design space.

The seven noise types occupy distinct positions in the stochastic geometry landscape:

| Noise type | Mechanism | Natural topology | Beam character | Directional control |
|---|---|---|---|---|
| Simplex | Gradient-hash lattice | Sheet / half-solid | Low | None |
| Cellular (F2−F1) | Voronoi distance field | Cell wall network | None | Via metric |
| FBM | Octave-layered simplex | Sheet / half-solid | Low–moderate | None |
| Ridged FBM | `1 − |simplex|` per octave | Sharp ridges / spicules | Moderate | None |
| Billow | `|simplex|` per octave | Rounded lobes | Low | None |
| Domain Warped | FBM-displaced FBM | Flowing organic | Low | None |
| Curl | Divergence-free ∇×Ψ | Tube networks | High | Via potential scale |

---

## Usage

Open `index.html` in Chrome or Edge alongside `mesher.js` in the same directory. No server, no build step, no dependencies.

**Recommended workflow:**

1. Select a noise type and set axis scaling for any desired anisotropy
2. Explore parameters in **Section 2D** — fast, interactive, full resolution; toggle X/Y/Z slice axis to read the 3D field structure
3. Switch to **Sheet 3D** to read spatial character and surface quality
4. Adjust topology (Sheet / Solid) and scaffold geometry parameters
5. Run **homogenization** to get estimated mechanical properties
6. **Export STL** for physical validation prints, or **Export JSON** for downstream pipeline handoff

---

## Noise Types

### Simplex

```
φ(x) = snoise(x · S)
```

Standard 3D simplex noise (Gustavson 2012). Smooth, continuous gradient field with no preferred direction. The foundational noise type — use it as a baseline before exploring the fractal variants. At mid volume fraction produces interconnected open-pore sheets; at low volume fraction the minority phase forms loosely connected blobs.

**No additional parameters.**

---

### Cellular (F2−F1)

```
φ(x) = F₂(x) − F₁(x)    [Worley distance field]
```

The difference between the second and first nearest Voronoi cell distances. Produces a field that is zero along cell boundaries and positive everywhere else — the natural analog of foam cell walls or honeycomb networks.

**Additional parameter:**
- **Distance metric** — Euclidean (smooth rounded cells), Manhattan (square-edged cells), Chebyshev (cubic cells). Controls the shape of the Voronoi geometry without changing the topology.

---

### FBM (Fractal Brownian Motion)

```
φ(x) = (1/Z) · Σᵢ aᵢ · snoise(x · f · lacunarityⁱ)
```

Octave-layered simplex: each octave adds detail at a higher spatial frequency, weighted by `gainⁱ`. Produces multi-scale texture with consistent statistical character across scales.

**Additional parameters:**
- **Octaves** — number of stacked noise layers (1–8). Live uniform — no recompile on change.
- **Lacunarity** — frequency ratio between successive octaves. 2.0 = classic doubling.
- **Gain** — amplitude falloff per octave. 0.5 = each octave contributes half the previous amplitude.

---

### Ridged FBM

```
φ(x) = (1/Z) · Σᵢ aᵢ · (1 − |snoise(x · f · lacunarityⁱ)|)
```

Replaces each octave's raw noise with `1 − |noise|`, inverting the minima into sharp ridges. Produces prominent connected ridge networks — a close match to trabecular spicules and coral branching geometry.

Parameters identical to FBM.

---

### Billow

```
φ(x) = (1/Z) · Σᵢ aᵢ · |snoise(x · f · lacunarityⁱ)|
```

Replaces each octave with `|noise|`, folding the negative half of the distribution upward. Produces rounded, cloud-like lobes — the natural analog of open-cell foam and adipose tissue.

Parameters identical to FBM.

---

### Domain Warped

```
q(x) = fbm(x)
φ(x) = fbm(x + strength · q(x))
```

Applies FBM displacement to the input coordinates before evaluating the final FBM. The q-vector is evaluated at 2 octaves for performance; the final field uses full octave count. Produces flowing, non-repeating geometry with no straight features.

**Additional parameter:**
- **Warp strength** — magnitude of the displacement field. 0 = no warp. 1–2 = organic flow. 3+ = extreme folding.

---

### Curl

```
∇×Ψ(x)    where Ψ is a simplex vector potential
cx = (∂Ψz/∂y − ∂Ψy/∂z)
```

Takes the curl of a simplex-derived vector potential via forward-difference finite differences (9 snoise calls per evaluation). Produces a divergence-free field whose magnitude forms connected tube networks — the closest noise analog to vascular trees and Haversian canals.

**Additional parameters:**
- **Curl step** — finite-difference step size. Smaller = more accurate; larger = softer tube geometry.
- **Potential scale** — spatial frequency of the underlying vector potential.

> **Performance note:** Curl requires 9 snoise calls per ray step vs 1 for simplex. Shader compilation is longer and the 3D view is slower. Section 2D is recommended for parameter exploration.

---

## Axis Scaling

Independent per-axis frequency scaling applied before the noise evaluation:

```
φ(x) = noiseField(x · diag(scaleX, scaleY, scaleZ) · freq)
```

Stretches or compresses the noise field independently along each world axis, introducing geometric anisotropy without changing the noise type or global frequency.

| Scale < 1.0 | Scale = 1.0 | Scale > 1.0 |
|---|---|---|
| Features compressed along axis | Isotropic (default) | Features elongated along axis |
| Denser pores, thicker walls | — | Coarser pores, thinner walls |

Range: 0.25–4.0. The Lipschitz constant used by the raymarcher automatically accounts for the maximum axis scale. Homogenization reads the anisotropic field and reports directional moduli that reflect the scaling — providing a direct design handle on stiffness anisotropy.

---

## Scaffold Geometry

**Topology** controls how the scalar field is converted to solid/void geometry:

| Mode | SDF expression | Character |
|---|---|---|
| Sheet | `\|φ − center\| < half-width` | Two surfaces bounding a hollow channel |
| Solid | `\|φ − center\| > half-width` (inverted) | Filled slab with open voids |

- **Center** — shifts the isosurface threshold. Primary volume fraction control.
- **Half-width** — slab thickness. Low values produce thin-wall open-pore scaffolds.
- **Smoothing** — separable box-filter blur before isosurface extraction. Most useful with Cellular noise.

---

## Perforation (Secondary Porosity)

Carves a second population of discrete voids through the scaffold via Boolean subtraction:

```
SDF_final = max(SDF_scaffold, −SDF_void)
```

- **Spherical** — Voronoi F1 nearest-cell field. Spherical voids at randomized positions, analogous to the lacunar-canalicular network in cortical bone.
- **Organic** — offset simplex field. Irregular cave-like voids with smooth boundaries.

**Parameters:** Hole frequency and radius are live uniforms — no shader recompile on change.

> Perforation is disabled in Section 2D for performance. Fully visible in Sheet 3D.

---

## Homogenization (MIL-HS)

Estimates effective elastic properties from the voxelized scalar field using the Mean Intercept Length method coupled to the Hashin-Shtrikman upper bound.

**Outputs:**
- **Volume fraction ρ** — solid phase fraction
- **Ex, Ey, Ez** — directional Young's moduli (GPa)
- **Gxy** — shear modulus (GPa)
- **νeff** — effective Poisson's ratio
- **Zener A** — anisotropy ratio. A = 1 is perfectly isotropic
- **Stiffness ellipsoid** — rotating WebGL visualization of the directional stiffness surface

**Grid resolution** — 32³, 48³, or 64³.

**E_s (GPa)** — solid material Young's modulus (Ti6Al4V ≈ 114 GPa, PEEK ≈ 3.6 GPa, 316L ≈ 193 GPa).

> MIL-HS assumes statistical homogeneity. Results improve at higher spatial frequency and higher grid resolution. Noise fields are not periodic — results are statistical approximations.

---

## Mesh Export

`mesher.js` must be present in the same directory as `index.html`.

| Resolution | Approx. time (simplex) | Approx. time (curl) | Typical triangle count |
|---|---|---|---|
| 64³ | ~1s | ~4s | ~50K |
| 128³ | ~8s | ~30s | ~200K |
| 256³ | ~60s | ~4min | ~800K |

Pipeline: scalar field sampling → marching cubes → Taubin λ/μ smoothing (λ=0.5, μ=−0.53, 10 iterations) → quality check → binary STL download.

**Quality report** — triangle count, watertight status, boundary edges, non-manifold edges.

> 256³ is not recommended for Curl in the browser.

---

## Export JSON

The exported JSON carries the complete parameter set and, if homogenization has been run, the mechanical results. Compatible with `mesher.js` for offline mesh generation and forward-compatible with Lattica.mesh (in development).

```json
{
  "meta": {
    "version": "1.0.0",
    "tool": "noise-scaffold-explorer",
    "timestamp": "...",
    "preset": "simplex"
  },
  "surface": {
    "type": "noise",
    "noise_type": "simplex",
    "frequency": 0.3,
    "scale_x": 1.0,
    "scale_y": 1.0,
    "scale_z": 1.0,
    "center": 0.0,
    "half_width": 0.15,
    "smoothing": 0.0,
    "topology": "sheet"
  },
  "geometry": {
    "cell_scale": 3.333,
    "mode": "shell",
    "wall_thickness": 0.15
  },
  "homogenization": {
    "volume_fraction": 18.4,
    "Ex_GPa": 12.3,
    "method": "MIL-HS",
    "note": "statistical_homogeneity_assumed"
  }
}
```

---

## Part of the F13LD Suite

| Tool | Field type | Topology | Directional control |
|---|---|---|---|
| [F13LD.tpms](https://mshomper.github.io/f13ld.tpms) | Periodic trigonometric | Sheet / strut | Crystallographic symmetry |
| **F13LD.noise** | Isotropic stochastic | Sheet / solid | Per-axis frequency scaling |
| **F13LD.grain**| Anisotropic stochastic | Sheet / half-solid | VMF + orthotropic |

---

## Technical Notes

**Fully implicit principle** — all geometry is defined as a scalar field evaluated pointwise. No mesh generation or STL intermediates exist during design. The STL export invokes marching cubes intentionally as a fabrication step.

**Sheet 3D** — WebGL1 sphere tracer. The Lipschitz constant is estimated from the noise gradient bound and current axis scaling. Internal render resolution is capped at 420px; CSS upscaling handles display. Shader compilation is asynchronous (KHR_parallel_shader_compile); the shader indicator shows amber while compiling and cyan when ready. Switching between noise types that require different GLSL programs triggers a recompile.

**Section 2D** — CPU rasterization of the scalar field at the chosen axis/position. Accurate for all noise types at all parameter values. Toggle X/Y/Z to read field anisotropy introduced by axis scaling. Recommended for parameter exploration, particularly for Curl and Domain Warp.

**Browser compatibility** — Chrome and Edge recommended. Firefox supported. WebGL1 required for Sheet 3D. Section 2D functions without WebGL.

---

## License

MIT — see LICENSE file.

*Not a Robot Engineering · notarobot-eng.com*
