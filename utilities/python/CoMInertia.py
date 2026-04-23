import trimesh

mesh = trimesh.load_mesh(".stl")

scale = 1 

mesh.apply_scale(scale)

rho = 1000

if not mesh.is_watertight:
    print("Mesh is not watertight. Results may be unreliable.")

mass_properties = mesh.mass_properties

volume = mass_properties["volume"]

mass = rho * volume

inertia_unit_density = mass_properties["inertia"]

inertia_with_density = rho * inertia_unit_density

print("=== Mesh Properties (with Scale and Density) ===")
print("Applied Scale Factor:", scale)
print("Volume [m^3]:", volume)
print("Mass [kg] (rho = {}):".format(rho), mass)
print("Center of Mass [m]:", mass_properties["center_mass"])
print("Inertia Tensor [kg·m²]:\n", inertia_with_density)
