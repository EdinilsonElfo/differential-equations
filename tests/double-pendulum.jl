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
u0 = [2.0, 3.0, 0.0, 0.0]
t_final = 200.0
nt = 8000
dt = t_final / nt
step_animation = 10

pendulum = DoublePendulum(coeff, u0, dt)

t, u = solve(pendulum, dt, nt, rk4)
r = polar_to_cartesian.(u)
u = reduce(hcat, u)'
r = reduce(hcat, r)'

# for i in 1:10:nt
#     line_x = [0.0, r[i,1], r[i,3]]
#     line_y = [0.0, r[i,2], r[i,4]]
#     p = plot(line_x, line_y, title="Double Pedulum", legend=false, linewidth=5,
#         xlims=(-2.1, 2.1), ylims=(-2.1, 2.1), size=(600,600))
#     plot!(r[1:i, 3], r[1:i, 4], linewidth=1)
#     display(p)
# end
# readline()

function plot_animation!(r)
    anim = @animate for i in 1:step_animation:nt
        line_x = [0.0, r[i,1], r[i,3]]
        line_y = [0.0, r[i,2], r[i,4]]
        plot(r[1:i, 3], r[1:i, 4], linewidth=1,
            xlims=(-2.1, 2.1), ylims=(-2.1, 2.1), size=(600,600))
        plot!(line_x, line_y, title="Double Pedulum", legend=false, linewidth=5)
    end
    return anim
end

println("Generating Animation...")
anim = plot_animation!(r)
println("Animation generated.")
# plot!(r[:,1], r[:,2], title="Double Pedulum", label="Position")
# p2 = plot(u[:,1], u[:,2], title="Phase Plane", label="Solution")
# p = plot(p1, p2, layout=2, size=(1000, 600))
p = gif(anim, "../images/double-pendulum.gif", fps=20)
display(p)
# readline()

