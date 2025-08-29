# 🚀 Laval Nozzle Designer (Python GUI)

A **fully interactive graphical application** for designing and analyzing **converging-diverging nozzles (Laval nozzles)**.  
The program provides:

- **Cross-sectional geometry preview** (real-time update as parameters change).  
- **Isentropic flow calculation** (exit Mach number, choking condition).  
- **3D model generation** with CAD/mesh export (STL, OBJ, PLY, OFF, GLTF, STEP).  
- **User-friendly GUI** with tooltips, annotations, and warnings.  

---

## 📖 Background

The **Laval nozzle** accelerates a compressible fluid to supersonic speeds by passing through:
1. A **converging section** (accelerates to sonic at the throat).
2. A **throat** (minimum cross-sectional area, Mach = 1).
3. A **diverging section** (further accelerates to supersonic).

This program allows you to **design, visualize, and export** such nozzles interactively.

---

## ✨ Features

- **GUI with dark theme** (Tkinter + Matplotlib).  
- **Parameter input fields** with live tooltips.  
- **Cross-section visualization** with:
  - Segment coloring (`Flange`, `Converging`, `Throat`, `Diverging`).  
  - Annotated diameters, lengths, and expansion angle.  
- **Flow analysis**:
  - Computes **exit Mach number** using area–Mach relation.  
  - Checks choking conditions (p_b ≤ p_crit).  
  - Reports both **subsonic and supersonic branches**.  
- **Warnings**:
  - Alerts if **divergence half-angle** is outside 12°–18°.  
- **3D model export**:
  - STL, OBJ, PLY, OFF, GLTF (mesh-based).
  - STEP (solid CAD format, via pythonOCC).  

---

## 🛠️ Implementation Details

### 1. GUI (Tkinter)
- **Tkinter Frames** split layout:  
  - Left → Matplotlib canvas (nozzle preview).  
  - Right → Input panel with all parameters.  
- **ToolTip class** implemented manually:  
  - `Toplevel` popup with black background + white text.  
  - Appears on hover, disappears on leave.  
- **Dark mode styling**:
  - All `Label`, `Entry`, `Button`, `Text` use black background, white foreground.

---

### 2. Geometry Generation

#### (a) **Converging Section – Vitoshinski Curve**
```python
generate_vitoshinski_curve(r_start, r_end, length, step=0.01)
```
- Based on **Vitoshinski empirical formula** for smooth contraction.  
- Ensures smooth slope continuity.  
- Adaptive point selection: adds a point if:
  - Radial/axial difference > threshold (`min_distance`, `min_radius_diff`, `min_x_diff`).  

#### (b) **Throat**
- Straight cylindrical section of user-defined length.  
- Maintains constant radius = D_throat / 2.

#### (c) **Diverging Section – Parabolic Expansion**
```python
generate_parabolic_expansion(r_start, r_end, x_start, length, curvature=1.0)
```
- Creates a **parabolic curve** between throat and exit.  
- Linear interpolation included → blend controlled by `curvature`.  

#### (d) **Outer Wall + Flange**
- Outer wall: constant radius (`outer_dia/2`).  
- Flange: cylinder section before converging.  

---

### 3. Flow Calculation

#### (a) **Critical Pressure**
```math
p_crit = p₀ * (2 / (γ + 1))^(γ / (γ - 1))
```

#### (b) **Area–Mach Relation**
```math
A/A* = (1/M) * [(2/(γ+1)) * (1 + (γ-1)/2 * M²)]^((γ+1)/(2(γ-1)))
```
- Solved numerically using `scipy.optimize.fsolve`.  
- Two solutions:
  - **M_sub (< 1)** → subsonic branch.  
  - **M_sup (> 1)** → supersonic branch.  

#### (c) **Branch Selection**
- If `p_b ≤ p_crit` → supersonic (choked).  
- Else → subsonic (not choked).  

---

### 4. Visualization (Matplotlib)

- **Filled color regions**:
  - Flange → orange (#FFB74D).  
  - Converging → teal (#4DB6AC).  
  - Throat → purple (#7986CB).  
  - Diverging → red (#E57373).  
- **Annotations**:
  - Diameters → with arrows and Ø symbol.  
  - Lengths → double arrows with labels.  
  - Divergence angle → placed at midpoint with arc arrow.  
- **Dark theme styling**:
  - Axes, ticks, labels all white on black.  
- **Auto-scaling**:
  - Axis limits dynamically adjusted with margins.  
- **Symmetry**:
  - Plots both +r and -r to show full cross-section.  

---

### 5. File Export

#### (a) **Mesh-based (Trimesh)**
```python
trimesh.creation.revolve(profile, angle, sections=rotres)
```
- Revolves 2D profile around X-axis.  
- Supports STL, OBJ, PLY, OFF, GLTF.  
- Optionally auto-opens STL file (Windows).  

#### (b) **STEP Export (pythonOCC)**
Steps:
1. Convert profile points → `gp_Pnt`.  
2. Remove duplicate points (tolerance 1e-7).  
3. Build edges → `BRepBuilderAPI_MakeEdge`.  
4. Build wire → `BRepBuilderAPI_MakeWire`.  
5. Build face → `BRepBuilderAPI_MakeFace`.  
6. Revolve face → `BRepPrimAPI_MakeRevol`.  
7. Write STEP file via `STEPControl_Writer`.  

---

## 📂 File Structure

```
nozzle_gui.py    # Main program with GUI, geometry, solver, export
README.md        # Project documentation
```

---

## ▶️ Usage

```bash
python nozzle_gui.py
```

- Modify inputs on the right panel.  
- Cross-section preview updates instantly.  
- Check calculated Mach numbers and warnings.  
- Export 3D model → choose file format.  

---

## 📌 Parameters (Full List)

| Parameter | Role |
|-----------|------|
| Inlet diameter | Controls nozzle entry size |
| Throat diameter | Minimum area → sonic flow location |
| Exit diameter | Controls expansion & Mach number |
| Converging length | Length of contraction curve |
| Throat straight section | Stabilizes throat flow |
| Diverging length | Expansion distance |
| Outer wall diameter | Controls housing profile |
| γ (Gamma) | Specific heat ratio |
| Total pressure p₀ | Chamber pressure |
| Total temperature T₀ | Chamber temperature |
| Back pressure p_b | Ambient pressure |
| Flange length | Length before nozzle contraction |
| Flange outer diameter | Outer flange size |
| Rotation angle | Revolve sweep (e.g., 180° = half nozzle) |
| Rotation resolution | Number of revolution segments |

---

## 📦 Dependencies

- `numpy`
- `matplotlib`
- `scipy`
- `trimesh`
- `pythonocc-core`
- `tkinter` (default in Python)

Install with:
```bash
pip install numpy matplotlib scipy trimesh pythonocc-core
```

---

## 📜 License

MIT License.  
Free for personal, academic, and commercial use.  

---

## 👨‍💻 Author

Developed by **[Your Name]**.  
If you find this project useful, please ⭐ the repo!
