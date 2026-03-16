abstract type OdeProblem end

struct Pendulum <: OdeProblem
    initial_value :: Float64
    paramaters :: Vector{Float64}
end

function (problem::Pendulum)(t::Float64)
    return - sin(t)
end

