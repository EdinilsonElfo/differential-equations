struct StabilityRegion
    x :: Vector{Float64}
    y :: Vector{Float64}
    stability :: Matrix{Float64}
    center :: Tuple{Float64, Float64}
    radius :: Float64
    x_lims :: Tuple{Float64, Float64}
    y_lims :: Tuple{Float64, Float64}
end

function StabilityRegion(center, radius, method)
    x_start = center[1] - radius
    x_stop  = center[1] + radius
    y_start = center[2] - radius
    y_stop  = center[2] + radius

    xrange = range(x_start, x_stop, length=n_points)
    yrange = range(y_start, y_stop, length=n_points)

    x = collect(xrange)
    y = collect(yrange)

    stability = test_stability(x, y, method)

    return StabilityRegion(x, y, stability, center, radius, (x_start, x_stop), (y_start, y_stop))
end


struct StabilityProblem <: OdeProblem
    coeff :: Vector{Float64}
    initial_value :: Vector{Float64}
    time_step :: Float64
end

function (p::StabilityProblem)(time::Float64, value::Vector{Float64})
    du1 = p.coeff[1]*value[1] - p.coeff[2]*value[2]
    du2 = p.coeff[2]*value[1] + p.coeff[1]*value[2]
    return [du1, du2]
end

function test_stability(xrange, yrange, solver)
    u0 = [1.0, 0.0]
    t_final = 1.0
    nt = 10
    dt = t_final / nt
    tol = 0.95
    xsize = length(xrange)
    ysize = length(yrange)
    stability = zeros(xsize, ysize)

    for (i, x) in enumerate(xrange)
        for (j, y) in enumerate(yrange)
            coeff = [x/dt, y/dt]
            model = StabilityProblem(coeff, u0, dt)
            t, u = solve(model, dt, nt, solver)
            norm = u[end][1]^2 + u[end][2]^2
            if norm < tol
                stability[i,j] = 1.0
            else
                stability[i,j] = 0.0
            end
        end
    end
    return stability
end

function plot_stability_region(region::StabilityRegion; title="StabilityRegion"::String)
    p = heatmap(region.x, region.y, region.stability', title=title, 
        cbar=false, size=(600,600), color=cgrad([:white, :gray]),
        xlims=region.x_lims, ylims=region.y_lims)
    plot!([region.x_lims...], [0.0, 0.0], color=:black, legend=false)
    plot!([0.0, 0.0], [region.y_lims...], color=:black, legend=false)

    return p
end

