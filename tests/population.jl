using Plots

include("../src/ode-solver.jl")

struct PopulationGrowth <: OdeProblem
    growRatio :: Float64
    initial_value :: Float64
    time_step :: Float64
end

function (p::PopulationGrowth)(time::Float64, value::Float64)
    return p.growRatio * value
end

ratio = 0.5
u0 = 20.0
t_final = 1.0
nt = 10
dt = t_final / nt

# populationGrowth = OdeProblem(1.0, 2.0)
populationGrowth = PopulationGrowth(ratio, u0, dt)
display(populationGrowth)

t, u1 = solve(populationGrowth, dt, nt, euler)
t, u2 = solve(populationGrowth, dt, nt, midp)
t, u3 = solve(populationGrowth, dt, nt, rk2)
t, u4 = solve(populationGrowth, dt, nt, rk4)

p = plot(t -> u0*exp(ratio*t), label="Analytical Solution", xlims=(0, t_final))
scatter!(t, u1, ms=2.0, label="Explicit Euler")
scatter!(t, u2, ms=2.0, label="Midpoint")
scatter!(t, u3, ms=2.0, label="Runge-Kutta 2")
scatter!(t, u4, ms=2.0, label="Runge-Kutta 4")

display(p)
readline()
