import tkinter as tk
import tkinter.messagebox
from tkinter import ttk, filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import trimesh
import sys
import os
from scipy.optimize import fsolve  # 用于求解面积–马赫关系
import math

# 新增：pythonOCC 相关导入，用于 STEP 输出
from OCC.Core.gp import gp_Pnt, gp_Ax1, gp_Dir
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire,
                                     BRepBuilderAPI_MakeFace)
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeRevol
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.Core.IFSelect import IFSelect_RetDone

# ---------------- ToolTip Class ----------------
class ToolTip:
    """
    It creates a tooltip for a given widget as the mouse goes on it.
    """
    def __init__(self, widget, text='widget info'):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self.id = None
        self.widget.bind("<Enter>", self.show_tip)
        self.widget.bind("<Leave>", self.hide_tip)
    
    def show_tip(self, event=None):
        "Display text in tooltip window"
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 25
        y = y + self.widget.winfo_rooty() + 20
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                         background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)
    
    def hide_tip(self, event=None):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

# ---------------- Helper Functions ----------------

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
            if (distance >= min_distance or abs(dr) >= min_radius_diff or abs(dx) >= min_x_diff):
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
    linear = r_start + (r_end - r_start) * (x - x_start) / length
    a = (r_end - r_start) / (length**2)  # Parabolic coefficient
    parabolic = a * (x - x_start)**2 + r_start
    return x, curvature * parabolic + (1 - curvature) * linear

def calculate_exit_mach(d_throat, d_exit, gamma):
    """
    Calculate the exit Mach number for an isentropic nozzle flow,
    assuming choked flow at the throat and exit area ratio Ae/At.
    (此函数在新版本中不再使用)
    """
    A_t = np.pi * (d_throat / 2.0)**2
    A_e = np.pi * (d_exit / 2.0)**2
    area_ratio = A_e / A_t
    exponent = (gamma - 1.0) / gamma
    term = (area_ratio**exponent) - 1.0
    M_e = np.sqrt((2.0 / (gamma - 1.0)) * term)
    return M_e, area_ratio

# ---------------- GUI Application ----------------

class NozzleGUI:
    def __init__(self, master):
        self.master = master
        master.title("Laval Nozzle Designer")
        
        # 左侧预览区
        self.left_frame = tk.Frame(master)
        self.left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # 右侧参数及计算区
        self.right_frame = tk.Frame(master)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=10, pady=10)
        
        # 嵌入 Matplotlib 图形
        self.figure, self.ax = plt.subplots(figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.left_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # 参数列表（已去除“Output Format”项）
        self.param_vars = {}
        param_list = [
            ("Inlet diameter (mm)", "20.0", float, "内流道入口直径，控制喷管初始入口尺寸"),
            ("Throat diameter (mm)", "10.0", float, "喷管喉部直径，决定临界流动时的最小截面"),
            ("Exit diameter (mm)", "25.0", float, "喷管出口直径，影响膨胀段和出口马赫数"),
            ("Converging length (mm)", "20.0", float, "收敛段长度，控制入口至喉部的渐缩"),
            ("Throat straight section (mm)", "5.0", float, "喉部直段长度，保证喉部稳定达到 M=1"),
            ("Diverging length (mm)", "30.0", float, "发散段长度，决定气体膨胀和加速的空间"),
            ("Outer wall diameter (mm)", "35.0", float, "外壳直径，用于构建外轮廓"),
            ("Ratio of specific heats gamma", "1.4", float, "比热比，反映气体热力学性质（例如空气常用1.4）"),
            ("Total pressure p0 (Pa)", "500000", float, "入口总压（燃烧室全压），决定气体能量状态"),
            ("Total temperature T0 (K)", "3000", float, "入口总温，决定气体总能量"),
            ("Back pressure p_b (Pa)", "101325", float, "背压（环境压强），用于判断是否达到壅塞和超音速流动"),
            ("Flange length (mm)", "5.0", float, "法兰长度，用于定义连接部分尺寸"),
            ("Flange outer diameter (mm)", "40.0", float, "法兰外径，控制法兰外轮廓尺寸"),
            ("Rotation angle (°)", "360.0", float, "旋转角度，用于生成三维模型（360° 表示全转）"),
            ("Rotation resolution", "35", int, "旋转时的分段数，决定三维模型的平滑程度")
        ]
        
        for i, (label_text, default, cast_type, tooltip_text) in enumerate(param_list):
            label = tk.Label(self.right_frame, text=label_text)
            label.grid(row=i, column=0, sticky="w", padx=5, pady=2)
            var = tk.StringVar(value=default)
            entry = tk.Entry(self.right_frame, textvariable=var)
            entry.grid(row=i, column=1, padx=5, pady=2)
            self.param_vars[label_text] = (var, cast_type)
            ToolTip(entry, tooltip_text)
            var.trace_add("write", self.update_preview)
        
        # 文本框显示 Mach 数计算结果
        self.calc_text = tk.Text(self.right_frame, height=15, width=40)
        self.calc_text.grid(row=len(param_list), column=0, columnspan=2, padx=5, pady=5)
        self.calc_text.config(state="disabled")
        
        # 文本框显示警告信息
        self.warning_text = tk.Text(self.right_frame, height=5, width=40)
        self.warning_text.grid(row=len(param_list)+1, column=0, columnspan=2, padx=5, pady=5)
        self.warning_text.config(state="disabled")
        
        # “Generate File” 按钮
        self.stl_button = tk.Button(self.right_frame, text="Generate File", command=self.generate_file)
        self.stl_button.grid(row=len(param_list)+2, column=0, columnspan=2, padx=5, pady=10)
        
        self.generate_preview()
    
    def get_param(self, label_text):
        var, cast_type = self.param_vars[label_text]
        try:
            return cast_type(var.get())
        except Exception:
            return None
    
    def update_preview(self, *args):
        self.generate_preview()
    
    def generate_preview(self):
        try:
            D_inlet      = self.get_param("Inlet diameter (mm)")
            D_throat     = self.get_param("Throat diameter (mm)")
            D_exit       = self.get_param("Exit diameter (mm)")
            L_converging = self.get_param("Converging length (mm)")
            L_throat     = self.get_param("Throat straight section (mm)")
            L_exit       = self.get_param("Diverging length (mm)")
            outer_dia    = self.get_param("Outer wall diameter (mm)")
            gamma        = self.get_param("Ratio of specific heats gamma")
            p0           = self.get_param("Total pressure p0 (Pa)")
            T0           = self.get_param("Total temperature T0 (K)")
            p_b          = self.get_param("Back pressure p_b (Pa)")
            L_flange     = self.get_param("Flange length (mm)")
            D_flange     = self.get_param("Flange outer diameter (mm)")
            angle        = self.get_param("Rotation angle (°)")
            rotres       = self.get_param("Rotation resolution")
        except Exception:
            return
        
        if None in [D_inlet, D_throat, D_exit, L_converging, L_throat, L_exit, outer_dia, gamma, p0, T0, p_b, L_flange, D_flange, angle, rotres]:
            return
        
        delta_r = (D_exit - D_throat) / 2
        expansion_angle = np.degrees(np.arctan(delta_r / L_exit))
        
        warning_str = ""
        if not (12 <= expansion_angle <= 18):
            warning_str = (
                f"WARNING: Expansion half-angle {expansion_angle:.1f}° is outside the recommended range (12°-18°).\n"
                "Consider adjusting the exit diameter or diverging length."
            )
        self.warning_text.config(state="normal")
        self.warning_text.delete(1.0, tk.END)
        self.warning_text.insert(tk.END, warning_str)
        self.warning_text.config(state="disabled")
        
        x_flange = np.linspace(0, L_flange, 10)
        r_flange = np.full_like(x_flange, D_inlet / 2)
        x_conv, r_conv = generate_vitoshinski_curve(
            r_start=D_inlet / 2,
            r_end=D_throat / 2,
            length=L_converging
        )
        x_conv += L_flange
        x_throat = np.linspace(x_conv[-1], x_conv[-1] + L_throat, 5)
        r_throat = np.full_like(x_throat, D_throat / 2)
        x_para, r_para = generate_parabolic_expansion(
            r_start=D_throat / 2,
            r_end=D_exit / 2,
            x_start=x_throat[-1],
            length=L_exit,
            curvature=1.0,
            num_points=100
        )
        
        x_inner = np.concatenate([x_flange, x_conv, x_throat, x_para])
        r_inner = np.concatenate([r_flange, r_conv, r_throat, r_para])
        
        x_outer = x_inner.copy()
        r_outer = np.concatenate([
            np.full_like(x_flange, D_flange / 2),
            np.full(len(x_conv) + len(x_throat) + len(x_para), outer_dia / 2)
        ])
        
        self.ax.clear()
        colors = {
            'Flange': '#FFB74D',
            'Converging': '#4DB6AC',
            'Throat': '#7986CB',
            'Diverging': '#E57373'
        }
        def plot_segment(x, r_in, r_out, name, color):
            self.ax.fill_between(x, r_in, r_out, color=color, alpha=0.6)
            self.ax.fill_between(x, -r_in, -r_out, color=color, alpha=0.6)
            self.ax.text(np.mean(x), np.mean(r_out) + 1, name,
                         ha='center', va='bottom', fontsize=10,
                         bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='gray'))
        
        flange_mask = x_inner <= L_flange
        conv_mask   = (x_inner > L_flange) & (x_inner <= L_flange + L_converging)
        throat_mask = (x_inner > L_flange + L_converging) & (x_inner <= L_flange + L_converging + L_throat)
        exit_mask   = x_inner > L_flange + L_converging + L_throat
        
        plot_segment(x_inner[flange_mask], r_inner[flange_mask], r_outer[flange_mask], 'Flange', colors['Flange'])
        plot_segment(x_inner[conv_mask], r_inner[conv_mask], r_outer[conv_mask], 'Converging', colors['Converging'])
        plot_segment(x_inner[throat_mask], r_inner[throat_mask], r_outer[throat_mask], 'Throat', colors['Throat'])
        plot_segment(x_inner[exit_mask], r_inner[exit_mask], r_outer[exit_mask], 'Diverging', colors['Diverging'])
        
        self.ax.plot(x_inner, r_inner, 'k-', lw=1)
        self.ax.plot(x_inner, r_outer, 'k-', lw=1)
        self.ax.plot(x_inner, -r_inner, 'k-', lw=1)
        self.ax.plot(x_inner, -r_outer, 'k-', lw=1)
        self.ax.axhline(0, color='gray', lw=0.8)
        
        def annotate_diameter(x, r, name):
            self.ax.annotate(f'{name}\nØ{2*r:.1f} mm',
                             xy=(x, r), xytext=(5, 5),
                             textcoords='offset points',
                             arrowprops=dict(arrowstyle='->', lw=1),
                             fontsize=9)
        annotate_diameter(x_flange[0], D_flange/2, 'Flange OD')
        annotate_diameter(x_flange[-1], D_inlet/2, 'Inlet')
        annotate_diameter(x_throat[0], D_throat/2, 'Throat')
        annotate_diameter(x_para[-1], D_exit/2, 'Exit')
        
        def annotate_length(x1, x2, y, name):
            self.ax.annotate('', xy=(x1, y), xytext=(x2, y),
                             arrowprops=dict(arrowstyle='<->', lw=1))
            self.ax.text((x1+x2)/2, y+1, f'{name}\n{abs(x2-x1):.1f} mm',
                         ha='center', va='bottom', fontsize=9)
        annotate_length(0, L_flange, -D_flange/2-2, 'Flange')
        annotate_length(L_flange, L_flange+L_converging, outer_dia/2+3, 'Converging')
        annotate_length(x_throat[0], x_throat[-1], outer_dia/2+3, 'Throat')
        annotate_length(x_para[0], x_para[-1], outer_dia/2+3, 'Diverging')
        
        def annotate_angle(x_start, length, r_start, r_end):
            mid_x = x_start + length/2
            mid_r = (r_start + r_end)/2
            self.ax.annotate(f'{expansion_angle:.1f}°',
                             xy=(mid_x, mid_r),
                             xytext=(mid_x+5, mid_r+5),
                             arrowprops=dict(arrowstyle='->', connectionstyle="arc3"),
                             fontsize=10,
                             bbox=dict(boxstyle="round", fc="white"))
        annotate_angle(x_para[0], L_exit, D_throat/2, D_exit/2)
        
        self.ax.set_xlabel('Axial Position (mm)')
        self.ax.set_ylabel('Radius (mm)')
        self.ax.set_title('Laval Nozzle Cross-section')
        self.ax.grid(True, linestyle=':', alpha=0.6)
        
        x_min, x_max = np.min(x_inner), np.max(x_inner)
        x_margin = 0.1 * (x_max - x_min)
        self.ax.set_xlim(x_min - x_margin, x_max + x_margin)
        max_radius = max(np.max(r_outer), np.max(np.abs(r_outer)))
        y_margin = 0.1 * max_radius
        self.ax.set_ylim(-max_radius - y_margin, max_radius + y_margin)
        
        self.ax.axis('equal')
        self.figure.tight_layout()
        self.canvas.draw()
        
        A_t = np.pi * (D_throat/2)**2
        A_e = np.pi * (D_exit/2)**2
        area_ratio = A_e / A_t
        
        p_crit = p0 * (2/(gamma+1))**(gamma/(gamma-1))
        
        def area_mach(M):
            return (1.0/M) * ((2.0/(gamma+1))*(1 + (gamma-1)/2.0*M**2))**((gamma+1)/(2*(gamma-1)))
        
        mach_sub = fsolve(lambda M: area_mach(M) - area_ratio, 0.3)[0]
        mach_sup = fsolve(lambda M: area_mach(M) - area_ratio, 2.0)[0]
        
        if p_b <= p_crit:
            flow_condition = "Choked (supersonic branch)"
            mach_exit = mach_sup
        else:
            flow_condition = "Not choked (subsonic branch)"
            mach_exit = mach_sub
        
        calc_str = (
            "=== Exit Mach Number Calculation ===\n"
            "Assumptions: Ideal gas, isentropic flow, choked condition at throat if p_b <= p_crit.\n\n"
            f"1) Throat diameter, D_throat = {D_throat:.2f} mm\n"
            f"   Exit diameter,   D_exit   = {D_exit:.2f} mm\n"
            f"   Gamma = {gamma:.3f}\n"
            f"   Total pressure, p0 = {p0:.0f} Pa\n"
            f"   Back pressure,  p_b = {p_b:.0f} Pa\n\n"
            "2) Compute areas:\n"
            "   A_t = π * (D_throat/2)^2\n"
            "   A_e = π * (D_exit/2)^2\n"
            f"   A_t = {A_t:.4f} mm²\n"
            f"   A_e = {A_e:.4f} mm²\n"
            f"   => Area ratio (A_e/A_t) = {area_ratio:.4f}\n\n"
            "3) Critical pressure at throat:\n"
            "   p_crit = p0*(2/(gamma+1))^(gamma/(gamma-1))\n"
            f"   => p_crit = {p_crit:.0f} Pa\n\n"
            "4) Based on back pressure:\n"
            f"   p_b = {p_b:.0f} Pa => {flow_condition}\n\n"
            "5) Solve area-Mach relation:\n"
            "   For subsonic branch (M < 1):\n"
            f"      M_sub = {mach_sub:.4f}\n"
            "   For supersonic branch (M > 1):\n"
            f"      M_sup = {mach_sup:.4f}\n\n"
            f"6) Selected exit Mach number, M_e = {mach_exit:.4f}\n"
        )
        self.calc_text.config(state="normal")
        self.calc_text.delete(1.0, tk.END)
        self.calc_text.insert(tk.END, calc_str)
        self.calc_text.config(state="disabled")
    
    def generate_file(self):
        """生成 3D 喷管模型并导出文件，调用 Windows 的保存文件对话框，通过“Save as type”选择输出格式"""
        try:
            D_inlet      = self.get_param("Inlet diameter (mm)")
            D_throat     = self.get_param("Throat diameter (mm)")
            D_exit       = self.get_param("Exit diameter (mm)")
            L_converging = self.get_param("Converging length (mm)")
            L_throat     = self.get_param("Throat straight section (mm)")
            L_exit       = self.get_param("Diverging length (mm)")
            outer_dia    = self.get_param("Outer wall diameter (mm)")
            gamma        = self.get_param("Ratio of specific heats gamma")
            L_flange     = self.get_param("Flange length (mm)")
            D_flange     = self.get_param("Flange outer diameter (mm)")
            angle        = self.get_param("Rotation angle (°)")
            rotres       = self.get_param("Rotation resolution")
        except Exception:
            return
        
        x_flange = np.linspace(0, L_flange, 10)
        r_flange = np.full_like(x_flange, D_inlet/2)
        x_conv, r_conv = generate_vitoshinski_curve(
            r_start=D_inlet/2,
            r_end=D_throat/2,
            length=L_converging
        )
        x_conv += L_flange
        x_throat = np.linspace(x_conv[-1], x_conv[-1] + L_throat, 5)
        r_throat = np.full_like(x_throat, D_throat/2)
        x_para, r_para = generate_parabolic_expansion(
            r_start=D_throat/2,
            r_end=D_exit/2,
            x_start=x_throat[-1],
            length=L_exit,
            curvature=1.0,
            num_points=100
        )
        x_inner = np.concatenate([x_flange, x_conv, x_throat, x_para])
        r_inner = np.concatenate([r_flange, r_conv, r_throat, r_para])
        x_outer = x_inner.copy()
        r_outer = np.concatenate([
            np.full_like(x_flange, D_flange/2),
            np.full(len(x_conv) + len(x_throat) + len(x_para), outer_dia/2)
        ])
        
        # 调用保存文件对话框，通过“Save as type”选择输出格式
        filetypes = [
            ("STL files", "*.stl"),
            ("STEP files", "*.step"),
            ("OBJ files", "*.obj"),
            ("PLY files", "*.ply"),
            ("OFF files", "*.off"),
            ("GLTF files", "*.gltf"),
            ("All files", "*.*")
        ]
        filename = filedialog.asksaveasfilename(title="Save file", defaultextension=".stl", filetypes=filetypes)
        if not filename:
            return
        
        ext = os.path.splitext(filename)[1].lower()
        
        if ext in [".stl", ".obj", ".ply", ".off", ".gltf"]:
            profile = np.vstack([
                np.column_stack((r_inner, x_inner)),
                np.column_stack((r_outer[::-1], x_outer[::-1]))
            ])
            try:
                nozzle = trimesh.creation.revolve(
                    profile,
                    angle=np.radians(angle),
                    cap_ends=True,
                    sections=rotres
                )
                nozzle.export(filename)
                tk.messagebox.showinfo("File Export", f"{ext.upper()} file saved: {filename}")
                if ext == ".stl":
                    self.master.after(0, lambda: os.startfile(filename))
            except Exception as e:
                tk.messagebox.showerror("File Export Error", str(e))
        elif ext == ".step":
            if D_flange/2 <= D_inlet/2:
                tk.messagebox.showerror("File Export Error", "Flange outer diameter must be larger than inlet diameter.")
                return
            if outer_dia/2 <= max(r_inner):
                tk.messagebox.showerror("File Export Error", "Outer wall diameter must be larger than all inner radii.")
                return
            try:
                pts = []
                for xi, ri in zip(x_inner, r_inner):
                    pts.append(gp_Pnt(xi, ri, 0))
                for xi, ri in zip(x_outer[::-1], r_outer[::-1]):
                    pts.append(gp_Pnt(xi, ri, 0))
                
                tolerance = 1e-7
                filtered_pts = []
                if pts:
                    filtered_pts.append(pts[0])
                    for pt in pts[1:]:
                        last_pt = filtered_pts[-1]
                        dx = pt.X() - last_pt.X()
                        dy = pt.Y() - last_pt.Y()
                        dz = pt.Z() - last_pt.Z()
                        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                        if dist > tolerance:
                            filtered_pts.append(pt)
                    first_pt = filtered_pts[0]
                    last_pt = filtered_pts[-1]
                    dx = first_pt.X() - last_pt.X()
                    dy = first_pt.Y() - last_pt.Y()
                    dz = first_pt.Z() - last_pt.Z()
                    if math.sqrt(dx*dx + dy*dy + dz*dz) < tolerance:
                        filtered_pts.pop()
                
                if len(filtered_pts) < 3:
                    tk.messagebox.showerror("File Export Error", "Not enough distinct points to form a valid profile.")
                    return
                
                edges = []
                for i in range(len(filtered_pts)):
                    p1 = filtered_pts[i]
                    p2 = filtered_pts[(i+1) % len(filtered_pts)]
                    maker_edge = BRepBuilderAPI_MakeEdge(p1, p2)
                    if not maker_edge.IsDone():
                        tk.messagebox.showerror("File Export Error", f"Edge creation failed between point {i} and {i+1}.")
                        return
                    edge = maker_edge.Edge()
                    edges.append(edge)
                
                wire_maker = BRepBuilderAPI_MakeWire()
                for edge in edges:
                    wire_maker.Add(edge)
                if not wire_maker.IsDone():
                    tk.messagebox.showerror("File Export Error", "Wire creation failed.")
                    return
                wire = wire_maker.Wire()
                face_maker = BRepBuilderAPI_MakeFace(wire)
                if not face_maker.IsDone():
                    tk.messagebox.showerror("File Export Error", "Face creation failed.")
                    return
                face = face_maker.Face()
                axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))
                revol_maker = BRepPrimAPI_MakeRevol(face, axis, np.radians(angle))
                solid = revol_maker.Shape()
                step_writer = STEPControl_Writer()
                step_writer.Transfer(solid, STEPControl_AsIs)
                status = step_writer.Write(filename)
                if status == IFSelect_RetDone:
                    tk.messagebox.showinfo("File Export", f"STEP file saved: {filename}")
                else:
                    tk.messagebox.showerror("File Export Error", "STEP export failed.")
            except Exception as e:
                tk.messagebox.showerror("File Export Error", str(e))
        else:
            tk.messagebox.showerror("File Export Error", f"Unknown file extension: {ext}")
        
def on_closing():
    root.destroy()
    sys.exit(0)

if __name__ == "__main__":
    root = tk.Tk()
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    root.geometry(f"{int(screen_width*0.8)}x{int(screen_height*0.8)}")
    root.protocol("WM_DELETE_WINDOW", on_closing)
    app = NozzleGUI(root)
    root.mainloop()
