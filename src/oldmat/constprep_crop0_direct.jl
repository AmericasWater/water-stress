using OptiMimi

include("world.jl")

include("Market.jl")
include("Transportation.jl")

# Variables set here by optimization: Market.produced, Transportation.imported
# Constraint: available > 0

importedparameters = numedges * numcrops * numsteps
producedparameters = numcounties * numcrops * numsteps
availablevariables = numcounties * numcrops * numsteps
A = spzeros(availablevariables, importedparameters + producedparameters)

A[:, 1:importedparameters] += gradients_market_regionimports_available() * gradients_transportation_imported_regionimports()
A[:, 1:importedparameters] += gradients_market_regionexports_available() * gradients_transportation_imported_regionexports()
A[:, importedparameters+1:importedparameters + producedparameters] = gradients_market_produced_available()

b = soleconstvector_market_available()

f = gradients_transportation_imported_cost() * soleobjvector_transportation_cost()

crop0constraints = MatrixConstraintSet([:Transportation, :Market], [:imported, :produced], [0., 0.], [Inf, Inf], f, A, b)
serialize(open(joinpath(todata, "crop0constraints-t$numsteps$suffix.jld"), "w"), crop0constraints)
