using OptiMimi

include("world.jl")

include("Market.jl")
include("Transportation.jl")

m = newmodel();

transportation = inittransportation(m); # optimization-only
market = initmarket(m); # dep. Transporation, Agriculture

# Variables set here by optimization: Market.produced, Transportation.imported
# Constraint: available > 0

market[:regionimports] = transportation[:regionimports];
market[:regionexports] = transportation[:regionexports];

println("Creating LP Constraints")

# Make a network constraint for county rr, crop cc, time tt
constraints = repmat(Function[m -> 0], numcounties * numcrops * numsteps)
for rr in 1:numcounties
    println("County $rr")
    for cc in 1:numcrops
        for tt in 1:numsteps
            constraints[(rr - 1) * numcrops * numsteps + (cc - 1) * numsteps + tt] = model -> -model[:Market, :available][rr, cc, tt]
        end
    end
end

println("Calculating constraints")

crop0constraints = savelpconstraints(m, [:Market, :Transportation], [:produced, :imported], [0., 0.], [Inf, Inf], soleobjective_transportation, constraints)
serialize(open(joinpath(todata, "crop0constraints-t$numsteps$suffix.jld"), "w"), crop0constraints)
