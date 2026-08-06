import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib as mpl
from scipy.interpolate import splprep, splev
from scipy.interpolate import CubicSpline

import matplotlib.pyplot as plt

def interpolate_missing(traj):
    """
    Given a list/array of 2D points (each as [a, b]) where missing points
    are marked with [-1, -1], fill the missing points by linear interpolation.
    If missing values occur at the beginning or the end, they are replaced
    with the nearest valid point.
    """
    traj = np.array(traj, dtype=float)
    n = len(traj)
    valid = np.all(traj != -1, axis=1)
    
    # If all points are missing, return zero array.
    if not np.any(valid):
        return np.zeros_like(traj)
    
    # For each coordinate column, interpolate separately.
    for col in range(traj.shape[1]):
        vals = traj[:, col]
        indices = np.arange(n)

        # Get indices for valid points
        valid_idx = indices[valid]
        valid_vals = vals[valid]

        # If the first or last points are missing, fill them with first/last valid value.
        if valid_idx[0] > 0:
            vals[:valid_idx[0]] = valid_vals[0]
        if valid_idx[-1] < n - 1:
            vals[valid_idx[-1]+1:] = valid_vals[-1]
        # Interpolate for the missing points.
        missing = ~valid
        vals[missing] = np.interp(indices[missing], valid_idx, valid_vals)
        traj[:, col] = vals
    return traj

def plot_3d_traj(left_view, top_view, save_path):
    """
    绘制3d轨迹。输入：
      left_view: 左视图轨迹，假设每个点为 [y, z]
      top_view:  上视图轨迹，假设每个点为 [x, y]
    若某个坐标点为 [-1, -1]，则采用前后线性插值进行处理。
    该函数将绘制出三维轨迹 (x, y, z)，轨迹颜色随时间变化（彩虹色），且不添加 marker。
    """
    # 处理缺失值：用前后插值处理
    left_view_filled = interpolate_missing(left_view)
    top_view_filled = interpolate_missing(top_view)
    
    # 提取坐标：
    # 假设上视图提供 x 和 y 坐标，左视图提供 y 和 z 坐标。
    # 这里取 top_view 的 x 坐标，top_view 的 y 坐标，
    # 左视图的第二个值作为 z 坐标。
    x = top_view_filled[:, 0]
    y = top_view_filled[:, 1]
    z = left_view_filled[:, 1]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # # 构造线段以便多彩着色
    # points = np.array([x, y, z]).T.reshape(-1, 1, 3)
    # segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # 对轨迹进行平滑处理
    t = np.linspace(0, 1, len(x))
    cs_x = CubicSpline(t, x)
    cs_y = CubicSpline(t, y)
    cs_z = CubicSpline(t, z)
    t_new = np.linspace(0, 1, 10000)
    x_new = cs_x(t_new)
    y_new = cs_y(t_new)
    z_new = cs_z(t_new)

    points = np.array([x_new, y_new, z_new]).T.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    # 创建彩虹 colormap，颜色随时间变化
    norm = mpl.colors.Normalize(vmin=0, vmax=len(x_new) - 1)
    lc = Line3DCollection(segments, cmap='rainbow', norm=norm)
    lc.set_array(np.arange(len(x_new) - 1))
    lc.set_linewidth(2)
    ax.add_collection(lc)

    # 设置字体为times new roman
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 10

    # 自动缩放坐标轴以适应数据
    ax.set_xlim([np.min(x), np.max(x)])
    ax.set_ylim([np.min(y), np.max(y)])
    ax.set_zlim([np.min(z), np.max(z)])

    ax.set_xlabel("X", fontdict={'family': 'Times New Roman'})
    ax.set_ylabel("Y", fontdict={'family': 'Times New Roman'})
    ax.set_zlabel("Z", fontdict={'family': 'Times New Roman'})
    ax.set_title("3D Trajectory", fontdict={'family': 'Times New Roman'})

    # 添加颜色条
    cbar = fig.colorbar(lc, ax=ax, pad=0.2, shrink=0.5)
    cbar.set_label('Time (s)')
    # ticks = np.linspace(0, len(x_new) - 1, num=int(len(x) / 60))
    # cbar.set_ticks(ticks)
    # ratio = int(len(x_new) / len(x) * 60)
    # tick_labels = [f"{int(tick / ratio)}" for tick in ticks]
    tick_seconds = np.arange(0, len(x), 60)
    tick_locations = tick_seconds * (len(x_new) - 1) / (len(x) - 1)
    tick_labels = [str(int(sec / 60)) for sec in tick_seconds]
    cbar.set_ticks(tick_locations)
    cbar.set_ticklabels(tick_labels)
    # cbar.set_ticks([0, len(x) - 1])
    # cbar.set_ticklabels([0, len(x) - 1])
    # 去除坐标轴背景
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # 去除坐标轴刻度线上的数字，保留刻度线
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    # 设置初始角度
    ax.view_init(elev=-160, azim=-130, roll=2)

    plt.show() 
    # 保存为pdf
    fig.savefig(save_path, bbox_inches='tight', dpi=300)
    fig.savefig(save_path.replace('.pdf', '.svg'), bbox_inches='tight', dpi=300)

# 示例调用：
if __name__ == '__main__':
    # 示例数据：用 -1 表示缺失的点
    # 假设top_view中的每个点格式为 [x, y]
    src_data = '18-31'

    top_view = np.load(f'all/{src_data}/Top_2025-03-19_{src_data}.npy')
    # 假设left_view中的每个点格式为 [y, z]
    left_view = np.load(f'all/{src_data}/Left_2025-03-19_{src_data}.npy')

    # save_path = 'all/16-58/3D_2025-03-19_16-58.pdf'
    # save_path = 'all/17-00/3D_2025-03-19_17-00.pdf'
    save_path = f'all/{src_data}/3D_2025-03-19_{src_data}.pdf'

    min_length = min(len(top_view), len(left_view))
    top_view = top_view[:min_length]
    left_view = left_view[:min_length]

    # 将每条轨迹的初始坐标点设置为0
    # top_view = top_view - top_view[0]
    # left_view = left_view - left_view[0]
    
    plot_3d_traj(left_view, top_view, save_path)
