f13ld.noise

Browser-based design tool for noise-derived metamaterial scaffolds

🔗 **[Launch the tool](https://mshomper.github.io/f13ld.noise)**

No installation. No server. Open the link and start designing.

---

What it does
field.noise generates three-dimensional porous scaffold geometry from mathematical noise fields — the same underlying mathematics used to simulate natural textures and terrains in computer graphics. Unlike crystallographic lattice structures (TPMS), noise-based scaffolds produce organic, statistically homogeneous geometries that closely resemble biological tissues such as trabecular bone, open-cell foam, and vascular networks.
The tool runs entirely in your browser using WebGL. You can explore designs visually in real time across a 2D cross-section view and an interactive 3D raymarched render, then export a JSON parameter file compatible with the companion TPMS Builder sweep tool for downstream mechanical analysis.

---

Noise types
Type	Character	Biological analog
Simplex	Smooth, rolling hills and valleys	Generic open-pore network
Cellular (F2−F1)	Voronoi cell boundaries	Honeycomb wall networks
FBM (Fractal)	Multi-scale layered texture	Fine-grained porous matrix
Ridged FBM	Sharp connected ridges	Trabecular spicules, coral branching
Billow	Rounded puffy lobes	Open-cell foam, adipose tissue
Domain Warped	Flowing, non-repeating organic forms	Biological soft tissue
Curl	Divergence-free tube networks	Vascular trees, Haversian canals
---
Controls
View
Section 2D — cross-section through the noise field. Blue = solid, red = void, white boundary = scaffold wall. Updates instantly.
Sheet 3D — interactive raymarched render. Drag to rotate, scroll to zoom.
Scaffold geometry
Topology — Sheet: thin-wall open-pore structure. Solid: bulk matrix with voids carved through it.
Center — slides the scaffold surface through the noise field. Primary volume fraction control.
Half-width — wall thickness, independent of center position.
Frequency — spatial scale of the noise. Higher = denser, finer pores.
Smoothing — spatial blur applied before surface extraction. Most useful with Cellular noise.
Fractal controls (FBM, Ridged, Billow, Domain Warp)
Octaves — number of stacked noise layers (1 = smooth, 8 = richly detailed).
Lacunarity — frequency ratio between octaves (default 2.0).
Gain — amplitude falloff between octaves (lower = subtler fine detail).
Warp strength — Domain Warp only. Controls how aggressively space is twisted.
Perforation (secondary porosity)
Punches a second population of discrete voids through the entire volume — analogous to the lacunar-canalicular network within a trabecular scaffold.
Spherical — Voronoi-based spherical holes at random positions.
Organic — Noise-field-derived irregular cave-like voids.
Hole frequency / radius — density and size of secondary pores, both live uniforms (no recompile on slider change).
---
Homogenization
Press Run homogenization to estimate mechanical properties via the MIL-HS method:
The noise field is voxelized at the selected grid resolution (32³, 48³, or 64³).
Mean Intercept Length (MIL) intercept counts give directional fabric tensors.
Hashin-Shtrikman upper bounds give isotropic stiffness from volume fraction.
Directional moduli Ex, Ey, Ez are scaled by normalized MIL weights.
Results include volume fraction, Ex/Ey/Ez (GPa), Gxy (GPa), Zener anisotropy index A, and effective Poisson's ratio. A live rotating directional stiffness ellipsoid visualizes the anisotropy.
> **Note:** noise fields are statistically homogeneous but not periodic. MIL-HS results are approximate and improve as frequency increases (more pore cycles in the domain).
---
JSON export
Export JSON produces a parameter file describing the scaffold. If homogenization has been run, mechanical results are included. The schema is structured to be compatible with the TPMS Builder sweep tool:
```json
{
  "meta": { "version": "1.0.0", "tool": "noise-scaffold-explorer", "preset": "simplex" },
  "surface": {
    "type": "noise",
    "noise_type": "simplex",
    "frequency": 0.3,
    "center": 0.0,
    "half_width": 0.15,
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
Browser compatibility
Browser	Section 2D	Sheet 3D
Chrome / Edge (desktop)	✅	✅ (Low quality recommended)
Safari / iOS	✅	✅
Firefox	✅	✅
The 3D view uses WebGL 1.0 raymarching with KHR_parallel_shader_compile for non-blocking shader compilation. Switching noise types shows a brief "Compiling..." status while the GPU driver compiles the new shader — the browser stays responsive throughout. On integrated GPUs, use Low quality (128 march iterations) for best performance.
---
Companion tools
This tool is part of an open-source computational materials design pipeline:
TPMS Builder — crystallographic TPMS surface design with PI-TPMS mode, FFT-CG homogenization sweep, and stiffness ellipsoid visualization.
JSON exports from the Noise Explorer use the same schema as the TPMS Builder, enabling direct comparison of stochastic and periodic scaffold architectures.
---
License
MIT — see LICENSE.
---
Developed by Matt Shomper / Not a Robot Engineering
