using Plots

include("../src/ode-solver.jl")
include("../src/stability-region.jl")

n_points = 400

center = (-1.0, 0.0)
radius = 3.0

# radius = 1.5
region = StabilityRegion(center, radius, n_points, euler)
p1 = plot_stability_region(region, title="Explicit Euler")

# radius = 2.0
region = StabilityRegion(center, radius, n_points, midp)
p2 = plot_stability_region(region, title="Midpoint Method")

# radius = 2.0
region = StabilityRegion(center, radius, n_points, heun)
p3 = plot_stability_region(region, title="Heun Method")

# radius = 2.0
region = StabilityRegion(center, radius, n_points, rk2)
p4 = plot_stability_region(region, title="Runge-Kutta 2")

# radius = 3.0
region = StabilityRegion(center, radius, n_points, rk3)
p5 = plot_stability_region(region, title="Runge-Kutta 3")

# radius = 3.0
region = StabilityRegion(center, radius, n_points, rk4)
p6 = plot_stability_region(region, title="Runge-Kutta 4")

p = plot(p1, p2, p3, p4, p5, p6, layout=6, size=(900, 600))

display(p)
readline()
# savefig("../images/stability.svg")

