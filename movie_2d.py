import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import matplotlib
from matplotlib.lines import Line2D
import tqdm

# Use Agg backend for headless rendering
matplotlib.use("Agg")

# Simulation parameters
dG = -0.02
lmbda = 100.0

# Load simulation
simf = "build/output/pepa.pvd"
pvrd = pv.PVDReader(simf)
n_frames = pvrd.number_time_points

# Set PyVista theme
pv.set_plot_theme('document')
pv.global_theme.allow_empty_mesh = True

# Define colors for eta fields
etas = ["solution_0", "solution_1", "solution_2", "solution_3"]
eta_colors = {"solution_0": "red", "solution_1": "green", "solution_2": "blue", "solution_3": "gold"}

# Initialize plotter centered at (0,0)
plotter = pv.Plotter(off_screen=True, window_size=(600, 600))
plotter.view_xy()
plotter.camera.SetParallelProjection(True)
plotter.camera.SetPosition(25, 25, 1)      # Centered at (0, 0)
plotter.camera.SetFocalPoint(25, 25, 0)    # Centered at (0, 0)
plotter.camera.SetParallelScale(25)

# Initialize video writer
video_path = f"simulation-2d_dG_{dG}_lam_{lmbda}_elliptic2.mp4"
writer = imageio.get_writer(video_path, fps=10, codec="libx264")

# Reference domain area (square 50x50)
V = 50**2

# Step size for frames
frame_step = 2
frame_indices = range(0, n_frames, frame_step)

with tqdm.tqdm(total=len(frame_indices), desc="Saving video") as pbar:
    for i in frame_indices:
        # Read current timestep
        pvrd.set_active_time_point(i)
        mesh = pvrd.read()[0]
        time_val = pvrd.active_time_value

        # Clear and set up plotter
        plotter.clear()
        plotter.add_mesh(mesh, color="gray", opacity=0.2)

        # Collect statistics
        stats = {}
        for eta in etas:
            mesh.set_active_scalars(eta)
            avg = mesh.integrate_data()[eta].sum() / V
            mx = mesh[eta].max()
            mn = mesh[eta].min()
            stats[eta] = (avg, mn, mx)

            # Clip and color
            clipped = mesh.clip_scalar(scalars=eta, value=0.5, invert=False)
            plotter.add_mesh(clipped, color=eta_colors[eta])

        # Take PyVista screenshot
        img = plotter.screenshot(return_img=True)

        # Create matplotlib figure with 2 columns (image + legend panel)
        fig = plt.figure(figsize=(9, 6), dpi=100)
        gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])

        ax_img = fig.add_subplot(gs[0])
        ax_info = fig.add_subplot(gs[1])

        # Show PyVista image on left
        ax_img.axis("off")
        ax_img.imshow(img)

        # Right panel: stats + legend
        ax_info.axis("off")

        # Time text
        ax_info.text(0.5, 1.0, f"Time $t^*$ = {time_val:.3f}",
                     fontsize=12, ha="center", va="top",
                     transform=ax_info.transAxes)

        # Legend entries with colored squares
        handles = []
        for eta in etas:
            avg, mn, mx = stats[eta]
            greek = f"$\\eta_{{{eta[-1]}}}$"
            label = f"{greek}: avg={avg:.2f}, min={mn:.2f}, max={mx:.2f}"
            
            square = Line2D([0], [0], marker="s", color="w",
                            markerfacecolor=eta_colors[eta],
                            markersize=12, label=label)
            handles.append(square)

        ax_info.legend(handles=handles,
                       loc="upper left", bbox_to_anchor=(0, 0.92),
                       fontsize=10, frameon=False)

        # Render figure to image array
        fig.tight_layout()
        fig.canvas.draw()
        rgba_img = np.asarray(fig.canvas.buffer_rgba())
        fig_img = rgba_img[:, :, :3]  # Convert RGBA -> RGB
        plt.close(fig)

        # Write to video
        writer.append_data(fig_img)
        pbar.update(1)

writer.close()