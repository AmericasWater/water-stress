using OptiMimi

include("world.jl")
include("weather.jl")

include("WaterNetwork.jl")
include("Agriculture.jl")
include("ConjunctiveUse.jl")

# Variables set here by optimization: ConjunctiveUse.pumping, Agriculture.irrigatedareas
# Constraint: outflow > 0

pumpingparameters = numregions * numsteps
areasparameters = numregions * numcrops * numsteps
outflowvariables = numgauges * numsteps
A = spzeros(outflowvariables, pumpingparameters + areasparameters)

A[:, 1:pumpingparameters] = gradients_conjunctiveuse_pumping_swdemand()
A[:, pumpingparameters+1:pumpingparameters + areasparameters] = gradients_agriculture_irrigatedareas_irrigation()

b = soleconstvector_conjunctiveuse_outflow()

f = gradients_conjunctiveuse_pumping_cost() * soleobjvector_conjunctiveuse_cost()

water0constraints = savelpconstraints(m, [:ConjunctiveUse, :Agriculture], [:pumping, :irrigatedareas], [0., 0.], [Inf, Inf], soleobjective_conjunctiveuse, constraints)
serialize(open(joinpath(todata, "water0constraints-t$numsteps$suffix.jld"), "w"), water0constraints)
