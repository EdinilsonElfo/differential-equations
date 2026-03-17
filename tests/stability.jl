using Plots

include("../src/ode-solver.jl")
include("../src/stability-region.jl")

n_points = 400

center = (-1.0, 0.0)

radius = 1.5
region = StabilityRegion(center, radius, euler)
p1 = plot_stability_region(region, title="Explicit Euler")

radius = 2.0
region = StabilityRegion(center, radius, midp)
p2 = plot_stability_region(region, title="Midpoint Method")

radius = 2.0
region = StabilityRegion(center, radius, heun)
p3 = plot_stability_region(region, title="Heun Method")

radius = 2.0
region = StabilityRegion(center, radius, rk2)
p4 = plot_stability_region(region, title="Runge-Kutta 2")

radius = 3.0
region = StabilityRegion(center, radius, rk4)
p5 = plot_stability_region(region, title="Runge-Kutta 4")

p = plot(p1, p2, p3, p4, p5, layout=6, size=(900, 600))

display(p)
readline()

