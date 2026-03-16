using Plots

include("../src/ode-solver.jl")

struct VanDerPol <: OdeProblem
    coeff :: Float64
    initial_value :: Vector{Float64}
    time_step :: Float64
end

function (p::VanDerPol)(time::Float64, value::Vector{Float64})
    du1 = value[2]
    du2 = p.coeff*(1.0 - value[1]^2)*value[2] - value[1]
    return [du1, du2]
end

coeff = 5.0
u0 = [2.0, 0.0]
t_final = 40.0
nt = 2000
dt = t_final / nt

model = VanDerPol(coeff, u0, dt)

t, u = solve(model, dt, nt, euler)
u = reduce(hcat, u)'

p1 = plot(t, u, title="Van der Pol", label=["Position" "Velocity"])
p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
p_euler = plot(p1, p2, layout=2)
# plot!(t -> u0*exp(ratio*t))
# display(p_euler)
# readline()

t, u = solve(model, dt, nt, rk4)
u = reduce(hcat, u)'

p1 = plot(t, u, title="Van der Pol", label=["Position" "Velocity"])
p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
p_midp = plot(p1, p2, layout=2)
p = plot(p_euler, p_midp, layout=(2,1))
# plot!(t -> u0*exp(ratio*t))
display(p)
readline()
