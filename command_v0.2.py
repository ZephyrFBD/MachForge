import numpy as np
import matplotlib.pyplot as plt
import trimesh

def get_user_input(prompt, default, cast_type=float):
    try:
        value = input(f"{prompt} [Default: {default}]: ").strip()
        return cast_type(value) if value else default
    except:
        print("Invalid input, using default value")
        return default

def generate_vitoshinski_curve(r_start, r_end, length, step=0.01, min_distance=0.5, min_radius_diff=0.05, min_x_diff=0.05):
    """Generate Vitoshinski contraction curve."""
    x = []
    r = []
    prev_x = None
    prev_r = None
    for xi in np.arange(0, length + step, step):
        p1 = 1 - (r_end / r_start)**2
        p2 = (1 - (xi / length)**2)**2
        p3 = (1 + (1/3) * (xi / length)**2)**3
        current_r = r_end / np.sqrt(1 - p1 * p2 / p3)

        add_point = False
        if prev_r is None:
            add_point = True
        else:
            dx = xi - prev_x
            dr = current_r - prev_r
            distance = np.sqrt(dx**2 + dr**2)
            if (distance >= min_distance or
                abs(dr) >= min_radius_diff or
                abs(dx) >= min_x_diff):
                add_point = True
        if xi >= length:
            add_point = True

        if add_point:
            x.append(xi)
            r.append(current_r)
            prev_x = xi
            prev_r = current_r
    return np.array(x), np.array(r)

def generate_parabolic_expansion(r_start, r_end, x_start, length, curvature=1.0, num_points=100):
    """Generate parabolic expansion section."""
    x = np.linspace(x_start, x_start + length, num_points)
    linear = r_start + (r_end - r_start)*(x - x_start)/length
    a = (r_end - r_start) / (length**2)  # Parabolic coefficient
    parabolic = a * (x - x_start)**2 + r_start
    return x, curvature*parabolic + (1-curvature)*linear

def calculate_exit_mach(d_throat, d_exit, gamma):
    """
    Calculate the exit Mach number for an isentropic nozzle flow,
    assuming choked flow at the throat and exit area ratio Ae/At.
    """
    # Area ratio A_e / A_t
    A_t = np.pi * (d_throat / 2.0)**2
    A_e = np.pi * (d_exit   / 2.0)**2
    area_ratio = A_e / A_t
    
    # Isentropic area relation
    # M_e = sqrt( (2/(gamma-1)) * [ (area_ratio)^((gamma-1)/gamma) - 1 ] )
    exponent = (gamma - 1.0) / gamma
    term = (area_ratio**exponent) - 1.0
    M_e = np.sqrt((2.0 / (gamma - 1.0)) * term)
    
    return M_e, area_ratio

def generate_laval_nozzle():
    print("\n=== Laval Nozzle Configuration ===")
    
    # Basic parameters
    D_inlet = get_user_input("Inlet diameter (mm)", 20.0)
    D_throat = get_user_input("Throat diameter (mm)", 10.0)
    D_exit = get_user_input("Exit diameter (mm)", 25.0)
    L_converging = get_user_input("Converging length (mm)", 20.0)
    L_throat = get_user_input("Throat straight section (mm)", 5.0)
    L_exit = get_user_input("Diverging length (mm)", 30.0)
    outer_dia = get_user_input("Outer wall diameter (mm)", 35.0)
    
    # Prompt for gas properties
    gamma = get_user_input("Ratio of specific heats gamma", 1.4)
    
    # Expansion angle validation
    delta_r = (D_exit - D_throat)/2
    expansion_angle = np.degrees(np.arctan(delta_r / L_exit))
    print(f"Expansion semi-angle {expansion_angle:.1f}° (recommended range 12-18°)")
    if not (12 <= expansion_angle <= 18):
        print(f"\nWarning: Expansion semi-angle {expansion_angle:.1f}° is outside the recommended range (12-18°).")
        print("Consider adjusting exit diameter or diverging length.")
    
    # Flange parameters
    L_flange = get_user_input("Flange length (mm)", 5.0)
    D_flange = get_user_input("Flange outer diameter (mm)", 40.0)
    
    # Generation parameters
    angle = get_user_input("Rotation angle (°)", 360.0)
    rotres = int(get_user_input("Rotation res", 35))
    stl_filename = input("STL filename [default: nozzle.stl]: ").strip() or "nozzle.stl"
    
    # ========== Generate inner profile ==========
    # Flange section
    x_flange = np.linspace(0, L_flange, 10)
    r_flange = np.full_like(x_flange, D_inlet/2)
    
    # Converging section
    x_conv, r_conv = generate_vitoshinski_curve(
        r_start=D_inlet/2,
        r_end=D_throat/2,
        length=L_converging
    )
    x_conv += L_flange  # coordinate shift
    
    # Throat section
    x_throat = np.linspace(x_conv[-1], x_conv[-1] + L_throat, 5)
    r_throat = np.full_like(x_throat, D_throat/2)
    
    # Parabolic expansion
    x_para, r_para = generate_parabolic_expansion(
        r_start=D_throat/2,
        r_end=D_exit/2,
        x_start=x_throat[-1],
        length=L_exit,
        curvature=1.0,
        num_points=100
    )
    
    # Combine inner profile
    x_inner = np.concatenate([x_flange, x_conv, x_throat, x_para])
    r_inner = np.concatenate([r_flange, r_conv, r_throat, r_para])
    
    # ========== Generate outer profile ==========
    x_outer = x_inner.copy()
    r_outer = np.concatenate([
        np.full_like(x_flange, D_flange/2),
        np.full(len(x_conv)+len(x_throat)+len(x_para), outer_dia/2)
    ])
    
    # ========== Generate 3D model ==========
    profile = np.vstack([
        np.column_stack((r_inner, x_inner)),
        np.column_stack((r_outer[::-1], x_outer[::-1]))
    ])
    nozzle = trimesh.creation.revolve(
        profile,
        angle=np.radians(angle),
        cap_ends=True,
        sections=rotres
    )
    nozzle.export(stl_filename)
    print(f"\nSTL file saved: {stl_filename}")
    
    # ========== Compute exit Mach number (isentropic assumption) ==========
    mach_exit, area_ratio = calculate_exit_mach(d_throat=D_throat, d_exit=D_exit, gamma=gamma)
    
    # ========== Print the calculation steps in English ==========
    print("\n=== Exit Mach Number Calculation ===")
    print("We assume an ideal isentropic flow with choked conditions at the throat.\n")
    print(f"1) Throat diameter, D_throat = {D_throat:.2f} mm")
    print(f"   Exit diameter,   D_exit   = {D_exit:.2f} mm")
    print(f"   Ratio of specific heats, gamma = {gamma:.3f}")
    print()
    print("2) Compute throat area (A_t) and exit area (A_e):")
    print("   A_t = π * (D_throat/2)^2")
    print("   A_e = π * (D_exit/2)^2")
    A_t = np.pi * (D_throat/2)**2
    A_e = np.pi * (D_exit/2)**2
    print(f"   A_t = {A_t:.4f} mm^2")
    print(f"   A_e = {A_e:.4f} mm^2")
    print()
    print(f"   => Area ratio (A_e / A_t) = {area_ratio:.4f}")
    print()
    print("3) Use the isentropic area–Mach relation for an exit Mach number M_e:")
    print("   M_e = sqrt[ (2 / (gamma - 1)) * { (A_e/A_t)^((gamma-1)/gamma) - 1 } ]")
    exponent = (gamma - 1) / gamma
    term = (area_ratio**exponent) - 1
    print(f"   => exponent = (gamma - 1)/gamma = ({gamma:.3f} - 1) / {gamma:.3f} = {exponent:.4f}")
    print(f"   => term = (A_e/A_t)^exponent - 1 = {area_ratio:.4f}^{exponent:.4f} - 1 = {term:.4f}")
    print(f"   => M_e = {mach_exit:.3f}")
    
    # ========== Plot cross-section for visualization ==========
    plt.figure(figsize=(12, 6))
    ax = plt.gca()
    
    # Color scheme
    colors = {
        'Flange': '#FFB74D',
        'Converging': '#4DB6AC',
        'Throat': '#7986CB',
        'Diverging': '#E57373'
    }
    
    def plot_segment(x, r_in, r_out, name, color):
        ax.fill_between(x, r_in, r_out, color=color, alpha=0.6)
        ax.fill_between(x, -r_in, -r_out, color=color, alpha=0.6)
        ax.text(np.mean(x), np.mean(r_out)+1, name,
                ha='center', va='bottom', fontsize=10,
                bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='gray'))
    
    # Boolean masks for each segment
    flange_mask = x_inner <= L_flange
    conv_mask = (x_inner > L_flange) & (x_inner <= L_flange+L_converging)
    throat_mask = (x_inner > L_flange+L_converging) & (x_inner <= L_flange+L_converging+L_throat)
    exit_mask = x_inner > L_flange+L_converging+L_throat
    
    plot_segment(x_inner[flange_mask], r_inner[flange_mask], r_outer[flange_mask], 
                 'Flange', colors['Flange'])
    plot_segment(x_inner[conv_mask], r_inner[conv_mask], r_outer[conv_mask], 
                 'Converging', colors['Converging'])
    plot_segment(x_inner[throat_mask], r_inner[throat_mask], r_outer[throat_mask], 
                 'Throat', colors['Throat'])
    plot_segment(x_inner[exit_mask], r_inner[exit_mask], r_outer[exit_mask], 
                 'Diverging', colors['Diverging'])
    
    # Plot contours
    ax.plot(x_inner, r_inner, 'k-', lw=1)
    ax.plot(x_inner, r_outer, 'k-', lw=1)
    ax.plot(x_inner, -r_inner, 'k-', lw=1)
    ax.plot(x_inner, -r_outer, 'k-', lw=1)
    ax.axhline(0, color='gray', lw=0.8)
    
    # Annotate diameters
    def annotate_diameter(x, r, name):
        ax.annotate(f'{name}\nØ{2*r:.1f}mm',
                    xy=(x, r), xytext=(5, 5),
                    textcoords='offset points',
                    arrowprops=dict(arrowstyle='->', lw=1),
                    fontsize=9)
    
    annotate_diameter(x_flange[0], D_flange/2, 'Flange OD')
    annotate_diameter(x_flange[-1], D_inlet/2, 'Inlet')
    annotate_diameter(x_throat[0], D_throat/2, 'Throat')
    annotate_diameter(x_para[-1], D_exit/2, 'Exit')
    
    # Annotate lengths
    def annotate_length(x1, x2, y, name):
        ax.annotate('', xy=(x1, y), xytext=(x2, y),
                    arrowprops=dict(arrowstyle='<->', lw=1))
        ax.text((x1+x2)/2, y+1, f'{name}\n{abs(x2-x1):.1f}mm',
                ha='center', va='bottom', fontsize=9)
    
    annotate_length(0, L_flange, -D_flange/2-2, 'Flange')
    annotate_length(L_flange, L_flange+L_converging, outer_dia/2+3, 'Converging')
    annotate_length(x_throat[0], x_throat[-1], outer_dia/2+3, 'Throat')
    annotate_length(x_para[0], x_para[-1], outer_dia/2+3, 'Diverging')
    
    # Annotate expansion angle
    def annotate_angle(x_start, length, r_start, r_end):
        """Annotate expansion angle."""
        mid_x = x_start + length/2
        mid_r = (r_start + r_end)/2
        ax.annotate(f'{expansion_angle:.1f}°', 
                    xy=(mid_x, mid_r),
                    xytext=(mid_x+5, mid_r+5),
                    arrowprops=dict(arrowstyle='->', connectionstyle="arc3"),
                    fontsize=10,
                    bbox=dict(boxstyle="round", fc="white"))
    
    annotate_angle(x_para[0], L_exit, D_throat/2, D_exit/2)
    
    ax.set_xlabel('Axial Position (mm)')
    ax.set_ylabel('Radial Radius (mm)')
    ax.set_title('Laval Nozzle Cross-Section')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.axis('equal')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    generate_laval_nozzle()