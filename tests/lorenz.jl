using Plots

include("../src/ode-solver.jl")

struct Lorenz <: OdeProblem
    coeff :: Vector{Float64}
    initial_value :: Vector{Float64}
    time_step :: Float64
end

function (p::Lorenz)(time::Float64, value::Vector{Float64})
    du1 = p.coeff[1]*(value[2] - value[1])
    du2 = value[1]*(p.coeff[2] - value[3]) - value[2]
    du3 = value[1]*value[2] - p.coeff[3]*value[3]
    return [du1, du2, du3]
end

coeff = [10.0, 28.0, 8/3]
u0 = [1.0, 0.0, 0.0]
t_final = 200.0
nt = 100000
dt = t_final / nt

model = Lorenz(coeff, u0, dt)

t, u = solve(model, dt, nt, rk4)
u = reduce(hcat, u)'

p = plot(t, u, title="Lorenz System", label=["x" "y" "z"])
display(p)
readline()

# p = plot(u[:,1], u[:,2], u[:,3], title="Phase Plane", label="Solution")
p = plot(u[:,1], u[:,2], u[:,3], xlabel="x", ylabel="y", zlabel="z",
legend=false, linewidth=0.5, size=(600, 600))
display(p)
readline()
savefig("../images/lorenz.svg")

