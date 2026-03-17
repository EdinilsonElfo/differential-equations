using Plots

include("../src/ode-solver.jl")

struct DoublePendulum <: OdeProblem
    coeff :: Float64
    initial_value :: Vector{Float64}
    time_step :: Float64
end

function (p::DoublePendulum)(time::Float64, value::Vector{Float64})
    q1 = value[1]
    q2 = value[2]
    p1 = value[3]
    p2 = value[4]

    a1 = p1*p2*sin(q1 - q2)
    a1 /= (1.0 + sin(q1 - q2)^2)
    a2 = (p1^2 + 2*p2^2 - 2*p1*p2*cos(q1 - q2)) * sin(q1 - q2) * cos(q1 - q2)
    a2 /= (1.0 + sin(q1 - q2)^2)^2

    dq1 =   (p1 - p2*cos(q1 - q2)) / (2.0 - cos(q1 - q2)^2)
    dq2 = (2*p2 - p1*cos(q1 - q2)) / (2.0 - cos(q1 - q2)^2)
    dp1 = -a1 + a2 - 2*sin(q1)
    dp2 =  a1 - a2 -   sin(q2)

    return [dq1, dq2, dp1, dp2]
end

function polar_to_cartesian(q::Vector{Float64})
    x1 = sin(q[1])
    y1 = -cos(q[1])
    x2 = x1 + sin(q[2])
    y2 = y1 - cos(q[2])

    return [x1, y1, x2, y2]
end

coeff = 1.0
u0 = [1.5, 2.0, 0.0, 0.0]
t_final = 250.0
nt = 10000
dt = t_final / nt

pendulum = DoublePendulum(coeff, u0, dt)

t, u = solve(pendulum, dt, nt, rk4)
r = polar_to_cartesian.(u)
u = reduce(hcat, u)'
r = reduce(hcat, r)'

p1 = plot(r[:,3], r[:,4], title="Double Pedulum", label="Position",
    xlims=(-2.1, 2.1), ylims=(-2.1, 2.1), size=(600,600))
# plot!(r[:,1], r[:,2], title="Double Pedulum", label="Position")
# p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
# p = plot(p1, p2, layout=2, size=(1000, 600))
display(p1)
readline()

