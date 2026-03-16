using Plots

include("../src/ode-solver.jl")

struct Pendulum <: OdeProblem
    coeff :: Float64
    initial_value :: Vector{Float64}
    time_step :: Float64
end

function (p::Pendulum)(time::Float64, value::Vector{Float64})
    dx = value[2]
    dv = -p.coeff*sin(value[1])
    return [dx, dv]
end

coeff = 0.2
u0 = [3.0, 0.0]
t_final = 100.0
nt = 300
dt = t_final / nt

pendulum = Pendulum(coeff, u0, dt)

t, u = solve(pendulum, dt, nt, euler)
u = reduce(hcat, u)'

p1 = plot(t, u, title="Pendulum", label=["position" "velocity"])
p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
p_euler = plot(p1, p2, layout=2)
# plot!(t -> u0*exp(ratio*t))
# display(p_euler)
# readline()

t, u = solve(pendulum, dt, nt, rk4)
u = reduce(hcat, u)'

p1 = plot(t, u, title="Pedulum", label=["position" "velocity"])
p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
p_midp = plot(p1, p2, layout=2)
p = plot(p_euler, p_midp, layout=(2,1))
# plot!(t -> u0*exp(ratio*t))
display(p)
readline()
