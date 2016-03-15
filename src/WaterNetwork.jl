# The Water Network component
#
# Determines how flows added and removed from the network propogate through.

using Mimi

@defcomp WaterNetwork begin
    gauges = Index()

    # External
    added = Parameter(index=[gauges, time]) # Water added at node
    removed = Parameter(index=[gauges, time]) # Water removed from node

    inflows = Variable(index=[gauges, time]) # Sum of upstream outflows
    outflows = Variable(index=[gauges, time]) # inflow + added - removed
end

"""
Compute the inflows and outflows at each node
"""
function timestep(c::WaterNetwork, tt::Int)
    v = c.Variables
    p = c.Parameters
    d = c.Dimensions

    for rr in d.gauges
        gg = vertex_index(downstreamorder[rr])
        gauge = downstreamorder[rr].label
        println("Process $gauge at $gg")
        allflow = 0.
        for upstream in out_neighbors(wateridverts[gauge], waternet)
            println(upstream)
            println(vertex_index(upstream, waternet))
            println(v.outflows[vertex_index(upstream, waternet), tt])
            allflow += v.outflows[vertex_index(upstream, waternet), tt]
            println(allflow)
        end

        v.inflows[gg, tt] = allflow
        v.outflows[gg, tt] = allflow + p.added[gg, tt] - p.removed[gg, tt]
    end
end
