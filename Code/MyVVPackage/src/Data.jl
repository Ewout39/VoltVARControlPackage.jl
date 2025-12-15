function _calc_load_current_max(load::Dict{String,<:Any}, bus::Dict{String,<:Any})::Vector{Float64}
    if all([haskey(load, prop) for prop in ["pmax", "pmin", "qmax", "qmin"]]) && haskey(bus, "vmin")
        pabsmax = max.(abs.(load["pmin"]), abs.(load["pmax"]))
        qabsmax = max.(abs.(load["qmin"]), abs.(load["qmax"]))
        smax = sqrt.(pabsmax.^2 + qabsmax.^2)
        vmin = bus["vmin"][[findfirst(isequal(c), bus["terminals"]) for c in load["connections"]]]

        return smax./vmin
    else
        return fill(Inf, length(load["connections"]))
    end
end