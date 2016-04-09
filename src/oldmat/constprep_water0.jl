using OptiMimi

include("world.jl")
include("weather.jl")

include("WaterNetwork.jl")
include("Agriculture.jl")
include("ConjunctiveUse.jl")

m = newmodel();

agriculture = initagriculture(m); # optimization-only
conjunctiveuse = initconjunctiveuse(m); # dep. Agriculture, DomesticDemand
waternetwork = initwaternetwork(m); # dep. ConjunctiveUse

# Variables set here by optimization: Agriculture.irrigatedareas
# Constraint: available > 0

conjunctiveuse[:totalirrigation] = agriculture[:totalirrigation];
conjunctiveuse[:domesticuse] = zeros(numcounties, numsteps)

waternetwork[:removed] = conjunctiveuse[:swdemand];

println("Creating LP Constraints")

# Make a network constraint for gauge gg, time tt
constraints = repmat(Function[m -> 0], numgauges * numsteps)
for rr in 1:numcounties
    println("Gauge $rr")
    for gg in 1:numgauges
        for tt in 1:numsteps
            constraints[(gg - 1) * numsteps + tt] = model -> -model[:WaterNetwork, :outflows][gg, tt]
        end
    end
end

println("Calculating constraints")

crop0constraints = savelpconstraints(m, [:ConjunctiveUse, :Agriculture], [:pumping, :irrigatedareas], [0., 0.], [Inf, Inf], soleobjective_conjunctiveuse, constraints)
serialize(open(joinpath(todata, "crop0constraints-t$numsteps$suffix.jld"), "w"), crop0constraints)
