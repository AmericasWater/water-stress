using OptiMimi

include("world.jl")
include("weather.jl")

include("Agriculture.jl")
include("ConjunctiveUse.jl")
include("Consumption.jl")
include("DomesticDemand.jl")
include("Market.jl")
include("Transportation.jl")
include("WaterNetwork.jl")

println("Creating model...")

# First solve entire problem in a single timestep
m = newmodel(1);

# Add all of the components
domesticdemand = initdomesticdemand(m, DomesticDemand); # exogenous
agriculture = initagriculture(m); # optimization-only
conjunctiveuse = initconjunctiveuse(m, years); # dep. Agriculture, DomesticDemand
waternetwork = addcomponent(m, WaterNetwork); # dep. ConjunctiveUse
transportation = inittransportation(m); # optimization-only
market = initmarket(m); # dep. Transporation, Agriculture
consumption = initconsumption(m, years); # dep. DomesticDemand, Market

# Connect up the components
agriculture[:irrigation] = farmerchoices[:irrigation];
agriculture[:areas] = farmerschoices[:areas];

waternetwork[:removed] = conjunctiveuse[:swdemand];

conjunctiveuse[:totalirrigation] = farmerchoices[:totalirrigation];
conjunctiveuse[:domesticuse] = domesticdemand[:waterdemand];

market[:produced] = extraction[:produced];
market[:regionimports] = transportation[:regionimports];
market[:regionexports] = transportation[:regionexports];

consumption[:waterdemand] = domesticdemand[:waterdemand];
consumption[:cropdemand] = domesticdemand[:cropdemand];
consumption[:cropavailable] = market[:available];

# Run it and time it!
@time run(m)

println("Testing:")
println(m[:Consumption, :surplus][1, 1])

println("Create linear optimization problem...")
# Make a network constraint for county rr, time tt
function makeconstraint(rr, tt)
    # The constraint function
    function constraint(model)
        -model[:Consumption, :surplus][rr, tt]
    end
end

# Set up the constraints
constraints = Function[]
for tt in 1:m.indices_counts[:time]
    constraints = [constraints; map(rr -> makeconstraint(rr, tt), 1:m.indices_counts[:regions])]
end

# Combine component-specific objectives
function objective(model::Model)
    soleobjective_extraction(model) + soleobjective_transportation(model)
end

# Create the OptiMimi optimization problem
optprob = problem(m, [:Extraction, :Transportation], [:pumping, :imported], [0., 0.], [Inf, Inf], objective, constraints=constraints, algorithm=:GUROBI_LINPROG);

println("Solving...")
@time sol = solution(optprob)
println(sol)

setparameters(m, [:Extraction, :Transportation], [:pumping, :imported], sol)
@time run(m)

df = DataFrame(fips=m.indices_values[:regions], demand=vec(m[:Consumption, :demand]),
               allotment=vec(m.components[:Extraction].Parameters.free_allotment),
               pumping=vec(m.components[:Extraction].Parameters.pumping),
               imports=vec(m[:Transportation, :regionimports]),
               exports=vec(m[:Transportation, :regionexports]))
writetable("results/counties$suffix.csv", df)


