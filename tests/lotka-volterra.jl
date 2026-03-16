using Plots

include("../src/ode-solver.jl")

struct LotkaVolterra <: OdeProblem
    coeff :: Vector{Float64}
    initial_value :: Vector{Float64}
    time_step :: Float64
end

function (p::LotkaVolterra)(time::Float64, value::Vector{Float64})
    du1 =   p.coeff[1]*value[1] - p.coeff[2]*value[1]*value[2]
    du2 = - p.coeff[3]*value[2] + p.coeff[4]*value[1]*value[2]
    return [du1, du2]
end

coeff = [1.5, 1.0, 1.0, 3.0]
u0 = [1.0, 2.0]
t_final = 20.0
nt = 400
dt = t_final / nt

model = LotkaVolterra(coeff, u0, dt)

t, u = solve(model, dt, nt, midp)
u = reduce(hcat, u)'

p1 = plot(t, u, title="Lotka-Volterra", label=["prey" "predator"])
p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
p_euler = plot(p1, p2, layout=2)
# plot!(t -> u0*exp(ratio*t))
# display(p_euler)
# readline()

t, u = solve(model, dt, nt, rk4)
u = reduce(hcat, u)'

p1 = plot(t, u, title="Lotka-Volterra", label=["prey" "predator"])
p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
p_midp = plot(p1, p2, layout=2)
p = plot(p_euler, p_midp, layout=(2,1))
# plot!(t -> u0*exp(ratio*t))
display(p)
readline()
