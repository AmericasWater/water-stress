# The waterdemand component
#
# Combines all of the sources of water demand, and determines the
# conjunctive use division between surface water and groundwater.

using Mimi
using DataFrames

@defcomp WaterDemand begin
    regions = Index()

    # External
    # Irrigation water (1000 m^3)
    totalirrigation = Parameter(index=[regions, time])
    # Combined water use for domestic sinks (1000 m^3)
    domesticuse = Parameter(index=[regions, time])

    # Optimized
    # How much is taking from groundwater
    gwportion = Parameter(index=[regions, time])

    # Internal
    # The cost in USD / 1000m^3 of pumping
    cost_pumping = Parameter(index=[regions, time])

    # Total water demand (1000 m^3)
    totaldemand = Variable(index=[regions, time])
    # Portion from surface water, in 1000 m^3
    swdemand = Parameter(index=[regions, time])
    # Groundwater to pump, in 1000 m^3
    gwdemand = Parameter(index=[regions, time])
    # The cost to pump it (USD)
    pumpingcost = Variable(index=[regions, time])
end

"""
Compute the amount extracted and the cost for doing it.
"""
function timestep(c::WaterDemand, tt::Int)
    v = c.Variables
    p = c.Parameters
    d = c.Dimensions

    for rr in d.regions
        # Sum all demands
        v.totaldemand[rr, tt] = p.totalirrigation[rr, tt] + p.domesticuse[rr, tt]

        # Split into surface and groundwater
        v.swdemand[rr, tt] = v.totaldemand[rr, tt] * (1 - p.gwportion[rr, tt])
        v.gwdemand[rr, tt] = v.totaldemand[rr, tt] * p.gwportion[rr, tt]

        # Total cost is pumping * cost-per-unit
        v.pumpingcost[rr, tt] = v.gwdemand[rr, tt] * p.cost_pumping[rr, tt]
    end
end

"""
Add a waterdemand component to the model.
"""
function initconjunctiveuse(m::Model, years)
    waterdemand = addcomponent(m, WaterDemand);

    # From http://www.oecd.org/unitedstates/45016437.pdf
    # Varies between 6.78 to 140 USD / 1000 m^3
    waterdemand[:cost_pumping] = 100. * ones(m.indices_counts[:regions], m.indices_counts[:time])

    waterdemand
end
