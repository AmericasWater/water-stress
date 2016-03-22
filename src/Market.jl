# The market component
#
# Determines the available resource for consumption, as a balance between local
# production, imports, and exports.

using Mimi

@defcomp Market begin
    regions = Index()
    crops = Index()

    # External
    # Local production from Agriculture
    produced = Parameter(index=[regions, crops, time])

    # Imports and exports from Transportation
    regionimports = Parameter(index=[regions, crops, time])
    regionexports = Parameter(index=[regions, crops, time])

    # The balance of available resource
    available = Variable(index=[regions, crops, time])
end

"""
Compute the available local resource for consumption, `available`.
"""
function timestep(c::Market, tt::Int)
    v = c.Variables
    p = c.Parameters
    d = c.Dimensions

    for rr in d.regions
        for cc in d.crops
            v.available[rr, cc, tt] = p.produced[rr, cc, tt] + p.regionimports[rr, cc, tt] - p.regionexports[rr, cc, tt]
        end
    end
end

"""
Add a market component to the model.
"""
function initmarket(m::Model)
    market = addcomponent(m, Market)

    market[:produced] = repeat([0.], outer=[m.indices_counts[:regions], m.indices_counts[:crops], m.indices_counts[:time]])

    market
end

