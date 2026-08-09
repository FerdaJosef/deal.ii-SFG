import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

pvd_path = "output/shrink.pvd"  # adjust to your actual filename

reader = pv.get_reader(pvd_path)
times = np.array(reader.time_values)

radii = np.full(times.shape, np.nan)

for i, t in enumerate(times):
    reader.set_active_time_point(i)
    mesh = reader.read()

    # If MultiBlock (common for parallel .pvtu output), merge partitions
    if isinstance(mesh, pv.MultiBlock):
        mesh = mesh.combine()

    # "solution" is however you named add_data_vector(...) in output_results()
    # It's likely stored as an n-component array — pick the component for
    # the grain you're tracking (0-indexed: eta1 might be component 0 or 1
    # depending on how your FE system orders the fields)
    field = mesh.point_data["solution_0"]

    if field.ndim > 1:
        eta = field[:, 0]   # <-- adjust index to the correct grain component
    else:
        eta = field

    mesh.point_data["eta_grain"] = eta

    thresholded = mesh.threshold(0.5, scalars="eta_grain")

    if thresholded.n_points == 0:
        continue  # grain has fully vanished

    sized = thresholded.compute_cell_sizes(length=False, area=True, volume=False)
    area = sized["Area"].sum()

    radii[i] = np.sqrt(area / np.pi)

# --- Inspect before fitting ---
mask = (times >= 2) & (times <= 5) & ~np.isnan(radii)
slope, intercept = np.polyfit(times[mask], radii[mask]**2, 1)
M_gamma_measured = -slope / 2
print(f"Measured M*gamma = {M_gamma_measured:.5f}")
print(f"Paper's M_kappa * gamma_omega_omega = {1.0 * 0.022:.5f}")
print(f"Ratio = {M_gamma_measured / 0.022:.3f}")