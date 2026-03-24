abstract type OdeProblem end
abstract type OdeSolver end

const OdeFunction = Union{Float64, Vector{Float64}}

struct ExplicitEuler <: OdeSolver end
struct MidPoint <: OdeSolver end
struct Heun <: OdeSolver end

struct RungeKutta <: OdeSolver
    a :: Matrix{Float64}
    b :: Vector{Float64}
    c :: Vector{Float64}
    order :: Int64
end

struct EmbeddedRungeKutta <: OdeSolver
    a :: Matrix{Float64}
    b1 :: Vector{Float64}
    b2 :: Vector{Float64}
    c :: Vector{Float64}
    order :: Int64
end

include("ode-methods.jl")

function (euler::ExplicitEuler)(problem::OdeProblem, time::Float64, value::OdeFunction)
    return value + problem(time, value)*problem.time_step
end

function (mp::MidPoint)(problem::OdeProblem, time::Float64, value::OdeFunction)
    dt2 = 0.5 * problem.time_step
    prob_tu = problem(time, value)
    return value + problem(time + dt2, value + dt2*prob_tu)*problem.time_step
end

function (heun::Heun)(problem::OdeProblem, time::Float64, value::OdeFunction)
    val_next = value + problem(time, value)*problem.time_step
    return value + 0.5*problem.time_step * (problem(time, value) + problem(time+problem.time_step, val_next))
end

function (rk::RungeKutta)(problem::OdeProblem, time::Float64, value::OdeFunction)
    k = Vector{OdeFunction}(undef, rk.order)
    sum_k = zero(value)
    for i in 1:rk.order
        sum_ak = zero(value)
        for j in 1:(i-1)
            sum_ak += rk.a[i,j]*k[j]
        end
        k[i] = problem(time + rk.c[i]*problem.time_step, value + sum_ak*problem.time_step)
        sum_k += rk.b[i]*k[i]
    end
    return value + sum_k * problem.time_step
end

function (rk::EmbeddedRungeKutta)(problem::OdeProblem, time::Float64, value::OdeFunction)
    k = Vector{OdeFunction}(undef, rk.order)
    sum_k1 = zero(value)
    sum_k2 = zero(value)
    for i in 1:rk.order
        sum_ak = zero(value)
        for j in 1:(i-1)
            sum_ak += rk.a[i,j]*k[j]
        end
        k[i] = problem(time + rk.c[i]*problem.time_step, value + sum_ak*problem.time_step)
        sum_k1 += rk.b1[i]*k[i]
        sum_k2 += rk.b2[i]*k[i]
    end
    val1 = value + sum_k1 * problem.time_step
    val2 = value + sum_k2 * problem.time_step
    error = val1 - val2
    return val1
end

function solve(problem::OdeProblem, step::Float64, N::Int, solver::OdeSolver)
    time = Float64[]
    solution = OdeFunction[]
    u = problem.initial_value
    t = 0.0
    push!(time, t)
    push!(solution, u)
    for i in 1:N
        t_new = t + step
        u_new = solver(problem, t, u)
        t = t_new
        u = u_new
        push!(time, t)
        push!(solution, u)
    end
    
    return time, solution
end

