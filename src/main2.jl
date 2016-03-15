include("world.jl")
include("WaterNetwork.jl")

m = newmodel(1);

waternetwork = addcomponent(m, WaterNetwork);

## These are given in downstreamorder order
# Basically, 1 unit added at ends and 1 unit removed from junctions
waternetwork[:added] = repmat(Float64[1, 0, 1, 0, 1], 1, 1);
waternetwork[:removed] = repmat(Float64[0, 1, 0, 1, 0], 1, 1);

run(m)

getdataframe(m, :WaterNetwork, :outflows) # Should be all 1's
