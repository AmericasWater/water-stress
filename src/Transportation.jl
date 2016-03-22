# The transporation component
#
# Each region is linked to other regions, accoring to the `regionnet`.  This
# component determines the total imports and exports for each region, based on the
# transport on each edge.

using Mimi

@defcomp Transportation begin
    regions = Index()
    edges = Index()
    crops = Index()

    # Internal
    # Cost per unit for transportation on a given edge
    cost_edge = Parameter(index=[edges, time])

    # Set by optimiation
    # Amount of resource imported on each link
    imported = Parameter(index=[edges, crops, time])

    # The costs for each edge's transportation
    cost = Variable(index=[edges, crops, time])

    # The total imported to and exported from each region
    regionimports = Variable(index=[regions, crops, time])
    regionexports = Variable(index=[regions, crops, time])
end

"""
Compute the amount imported and exported by region.
"""
function timestep(c::Transportation, tt::Int)
    v = c.Variables
    p = c.Parameters
    d = c.Dimensions

    # Costs are easy: just resource imported * cost-per-unit
    for ee in d.edges
        for cc in d.crops
            v.cost[ee, cc, tt] = p.imported[ee, cc, tt] * p.cost_edge[ee, tt]
        end
    end
    for cc in d.crops
        for ii in d.regions
            v.regionexports[ii, cc, tt] = 0.0
        end

        # Sum over all edges for each region to translate to region-basis
        edge1 = 1
        for ii in d.regions
            # Get the number of edges this county imports from
            numneighbors = out_degree(regverts[names[ii]], regionnet)

            # Sum over all *out-edges* to get import
            v.regionimports[ii, cc, tt] = sum(p.imported[edge1:edge1 + numneighbors - 1, cc, tt])

            # Sum over the edges that have this as an out-edge
            sources = get(sourceiis, ii, Int64[])
            for source in sources
                v.regionexports[source, cc, tt] += p.imported[edge1, cc, tt]
                edge1 += 1 # length(sources) == numneighbors
            end
        end
    end
end

"""
The objective of the transportation component is to minimize transport costs.
"""
function soleobjective_transportation(model::Model)
    sum(model[:Transportation, :cost])
end

"""
Add a transportation component to the model.
"""
function inittransportation(m::Model)
    transit = addcomponent(m, Transportation)

    # 10 USD / 1000 m^3
    transit[:cost_edge] = repmat([10.], m.indices_counts[:edges], m.indices_counts[:time])
    transit[:imported] = repeat([0.], outer=[m.indices_counts[:edges], m.indices_counts[:crops], m.indices_counts[:time]])

    transit
end

