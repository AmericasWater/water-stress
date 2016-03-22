# The consumption component
#
# Each region has an `demand`, and compares this to the amount of `available`
# resource.  The result is a `surplus`, which may be negative if `available` <
# `demand`.

using Mimi

@defcomp Consumption begin
    regions = Index()
    crops = Index()

    # External
    # Resource availability from Infrastructure
    wateravailable = Parameter(index=[regions, time])
    cropavailable = Parameter(index=[regions, crops, time])

    # Resource surplus over (or below) demand
    waterdemand = Variable(index=[regions, time])
    cropdemand = Variable(index=[regions, crops, time])

    # Internal
    watersurplus = Variable(index=[regions, time])
    cropsurplus = Variable(index=[regions, crops, time])
end

"""
Compute the `surplus` as `available` - `demand`.
"""
function timestep(c::Consumption, tt::Int)
    v = c.Variables
    p = c.Parameters
    d = c.Dimensions

    for rr in d.regions
        v.watersurplus[rr, tt] = p.wateravailable[rr, tt] - p.waterdemand[rr, tt]

        for cc in d.crops
            v.cropsurplus[rr, cc, tt] = p.cropavailable[rr, cc, tt] - p.cropdemand[rr, cc, tt]
        end
    end
end
