import numpy as np
import matplotlib.pyplot as plt
import trimesh
from scipy.interpolate import make_interp_spline

def get_user_input(prompt, default, cast_type=float):
    try:
        value = input(f"{prompt} [默认: {default}]: ").strip()
        return cast_type(value) if value else default
    except:
        print("输入无效，使用默认值")
        return default

def generate_laval_nozzle_with_flange():
    print("=== 拉瓦尔喷嘴参数设置（含加厚段） ===")
    D_inlet = get_user_input("入口直径 D_inlet (mm)", 20)
    D_throat = get_user_input("喉口直径 D_throat (mm)", 10)
    D_exit = get_user_input("出口直径 D_exit (mm)", 25)
    Outer_diameter = get_user_input("外轮廓最大直径 Outer_diameter (mm)", 30)

    L_flange = get_user_input("加厚段长度 L_flange (mm)", 5)
    D_flange = get_user_input("加厚段外径 D_flange (mm)", 36)

    L_straight = get_user_input("入口直段长度 L_straight (mm)", 10)
    L_curve1 = get_user_input("收缩段长度 L_curve1 (mm)", 10)
    L_curve2 = get_user_input("进一步收缩段 L_curve2 (mm)", 10)
    L_throat = get_user_input("喉口直段长度 L_throat (mm)", 5)
    L_exit = get_user_input("扩张段长度 L_exit (mm)", 30)

    angle = get_user_input("旋转角度 angle (°)", 360)
    num_points = get_user_input("插值点数 num_points", 300, int)
    stl_filename = input("输出 STL 文件名（例如 nozzle.stl）[默认: laval_nozzle.stl]: ").strip()
    if not stl_filename:
        stl_filename = 'laval_nozzle.stl'

    # 构建内外壁轮廓（带加厚段）
    x_flange = np.linspace(0, L_flange, 5)
    r_flange_inner = np.full_like(x_flange, D_inlet / 2)
    r_flange_outer = np.full_like(x_flange, D_flange / 2)

    x_straight = np.linspace(L_flange, L_flange + L_straight, 20)
    r_straight = np.full_like(x_straight, D_inlet / 2)

    x_key = [
        L_flange + L_straight,
        L_flange + L_straight + L_curve1,
        L_flange + L_straight + L_curve1 + L_curve2,
        L_flange + L_straight + L_curve1 + L_curve2 + L_throat,
        L_flange + L_straight + L_curve1 + L_curve2 + L_throat + L_exit
    ]
    r_key = [
        D_inlet / 2,
        (D_inlet + D_throat) / 4,
        D_throat / 2,
        D_throat / 2,
        D_exit / 2
    ]

    x_spline = np.linspace(x_key[0], x_key[-1], num_points)
    spline = make_interp_spline(x_key, r_key, k=3)
    r_spline = spline(x_spline)

    x_inner = np.concatenate([x_flange, x_straight, x_spline[1:]])
    r_inner = np.concatenate([r_flange_inner, r_straight, r_spline[1:]])
    r_outer = np.concatenate([
        r_flange_outer,
        np.full_like(np.concatenate([x_straight, x_spline[1:]]), Outer_diameter / 2)
    ])
    x_outer = x_inner.copy()

    # 构造轮廓剖面用于旋转建模
    profile = np.vstack([
        np.column_stack((r_inner, x_inner)),
        np.column_stack((r_outer[::-1], x_outer[::-1]))
    ])
    nozzle_3d = trimesh.creation.revolve(profile, angle=np.radians(angle), cap_ends=False)

    # 保存 STL 文件
    nozzle_3d.export(stl_filename)
    print(f"\nSTL 文件已保存为: {stl_filename}")

    # ------------------------------
    # 绘制剖面图：彩色段标注 + 半径标注
    # ------------------------------
    fig, ax = plt.subplots(figsize=(12, 6))

    # 工程配色
    colors = {
        'Flange': '#D95319',
        'Straight': '#0072BD',
        'Curve1': '#EDB120',
        'Curve2': '#7E2F8E',
        'Throat': '#77AC30',
        'Exit': '#A2142F',
    }

    segments = [
        ('Flange', 0, L_flange),
        ('Straight', L_flange, L_flange + L_straight),
        ('Curve1', L_flange + L_straight, L_flange + L_straight + L_curve1),
        ('Curve2', L_flange + L_straight + L_curve1, L_flange + L_straight + L_curve1 + L_curve2),
        ('Throat', L_flange + L_straight + L_curve1 + L_curve2,
         L_flange + L_straight + L_curve1 + L_curve2 + L_throat),
        ('Exit', L_flange + L_straight + L_curve1 + L_curve2 + L_throat,
         L_flange + L_straight + L_curve1 + L_curve2 + L_throat + L_exit)
    ]

    for name, x_start, x_end in segments:
        mask = (x_inner >= x_start) & (x_inner <= x_end)
        xi = x_inner[mask]
        ri = r_inner[mask]
        ro = r_outer[mask]
        ax.fill_between(xi, ri, ro, facecolor=colors[name], alpha=0.6)
        ax.fill_between(xi, -ri, -ro, facecolor=colors[name], alpha=0.6)

        # 段名 + 长度标注放在 x 轴
        x_center = (x_start + x_end) / 2
        ax.annotate(f'{name}\nL = {x_end - x_start:.0f} mm',
                    xy=(x_center, 0),
                    xytext=(x_center, 0),
                    textcoords='offset points',
                    ha='center',
                    va='bottom',
                    fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='black', lw=0.5))

    # 上下边界线
    ax.plot(x_inner, r_inner, 'k-')
    ax.plot(x_inner, r_outer, 'k-')
    ax.plot(x_inner, -r_inner, 'k-')
    ax.plot(x_inner, -r_outer, 'k-')
    ax.axhline(0, color='gray', linewidth=0.8)

    # 半径标注
    def mark_radius(x, r, label):
        ax.annotate(f'{label} = {2*r:.0f} mm',
                    xy=(x, r),
                    xytext=(x + 1, r + 2),
                    arrowprops=dict(arrowstyle='->', lw=1),
                    fontsize=10)

    mark_radius(2, D_flange / 2, 'D_flange')
    mark_radius(8, D_inlet / 2, 'D_inlet')
    mark_radius(x_key[2], D_throat / 2, 'D_throat')
    mark_radius(x_key[-1], D_exit / 2, 'D_exit')

    ax.set_xlabel('X (Length along nozzle, mm)')
    ax.set_ylabel('Y (Radius from axis, mm)')
    ax.set_title('Cross Section of Laval Nozzle')
    ax.axis('equal')
    ax.grid(True)
    plt.tight_layout()
    plt.show()

# 运行程序
generate_laval_nozzle_with_flange()