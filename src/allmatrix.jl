using Mimi
include("linproghouse.jl")

redohouse = true
redogwwo = true

include("world.jl")
include("weather.jl")

include("Agriculture.jl")
include("ConjunctiveUse.jl")
include("DomesticDemand.jl")
include("Market.jl")
include("Transportation.jl")
include("WaterNetwork.jl")

println("Creating model...")

# First solve entire problem in a single timestep
m = newmodel();

# Add all of the components
domesticdemand = initdomesticdemand(m, m.indices_values[:time]); # exogenous
agriculture = initagriculture(m); # optimization-only
conjunctiveuse = initconjunctiveuse(m); # dep. Agriculture, DomesticDemand
waternetwork = initwaternetwork(m); # dep. ConjunctiveUse
transportation = inittransportation(m); # optimization-only
market = initmarket(m); # dep. Transporation, Agriculture

# Only include variables needed in constraints and parameters needed in optimization

paramcomps = [:ConjunctiveUse, :ConjunctiveUse, :Agriculture, :Agriculture, :Transportation, :Market]
parameters = [:pumping, :withdrawals, :rainfedareas, :irrigatedareas, :imported, :internationalsales]
constcomps = [:Agriculture, :WaterNetwork, :ConjunctiveUse, :Market, :Market]
constraints = [:allagarea, :outflows, :swbalance, :available, :domesticbalance]

if redohouse
    house = LinearProgrammingHouse(m, paramcomps, parameters, constcomps, constraints);

    # Optimize revenue_domestic + revenue_international - pumping_cost - transit_cost
    setobjective!(house, -varsum(grad_conjunctiveuse_pumping_cost(m)))
    setobjective!(house, -varsum(grad_transportation_imported_cost(m)))
    setobjective!(house, deriv_market_internationalsales_totalrevenue(m))
    setobjective!(house, deriv_market_produced_totalrevenue(m) * room_relabel(grad_agriculture_rainfedareas_production(m), :production, :Market, :produced))
    setobjective!(house, deriv_market_produced_totalrevenue(m) * room_relabel(grad_agriculture_irrigatedareas_production(m), :production, :Market, :produced))

    # Constrain agriculture < county area
    setconstraint!(house, grad_agriculture_irrigatedareas_allagarea(m)) # +
    setconstraint!(house, grad_agriculture_rainfedareas_allagarea(m)) # +
    setconstraintoffset!(house, constraintoffset_agriculture_allagarea(m))

    # Constrain outflows + runoff > 0, or -outflows < runoff
    if redogwwo
        gwwo = grad_waternetwork_withdrawals_outflows(m);
        serialize(open(joinpath(todata, "partialhouse$suffix.jld"), "w"), gwwo);
        cwro = constraintoffset_waternetwork_runoff(m);
        serialize(open(joinpath(todata, "partialhouse2$suffix.jld"), "w"), cwro);
    else
        gwwo = deserialize(open(joinpath(todata, "partialhouse$suffix.jld"), "r"));
        cwro = deserialize(open(joinpath(todata, "partialhouse2$suffix.jld"), "r"));
    end

    setconstraint!(house, -room_relabel_parameter(gwwo, :withdrawals, :ConjunctiveUse, :withdrawals)) # +
    setconstraintoffset!(house, cwro) # +

    # Constrain available market > 0
    setconstraint!(house, -grad_market_produced_available(m) * room_relabel(grad_agriculture_irrigatedareas_production(m), :production, :Market, :produced)) # -
    setconstraint!(house, -grad_market_produced_available(m) * room_relabel(grad_agriculture_rainfedareas_production(m), :production, :Market, :produced)) # -
    setconstraint!(house, -grad_market_internationalsales_available(m)) # +
    setconstraint!(house, -(grad_market_regionimports_available(m) * grad_transportation_imported_regionimports(m) +
                   grad_market_regionexports_available(m) * grad_transportation_imported_regionexports(m))) # +-

    # Constrain swdemand < swsupply, or irrigation + domestic < pumping + withdrawals, or irrigation - pumping - withdrawals < -domestic
    setconstraint!(house, grad_conjunctiveuse_totalirrigation_swbalance(m) * grad_agriculture_irrigatedareas_totalirrigation(m)) # +
    setconstraint!(house, -grad_conjunctiveuse_pumping_swbalance(m)) # -
    setconstraint!(house, -grad_conjunctiveuse_withdrawals_swbalance(m)) # - THIS IS SUPPLY
    setconstraintoffset!(house, -hall_relabel(constraintoffset_domesticdemand_waterdemand(m), :waterdemand, :ConjunctiveUse, :swbalance)) # -

    # Constrain domesticsales < domesticdemand
    # Reproduce 'available'
    setconstraint!(house, room_relabel(grad_market_produced_available(m) * room_relabel(grad_agriculture_irrigatedareas_production(m), :production, :Market, :produced), :available, :Market, :domesticbalance)) # +
    setconstraint!(house, room_relabel(grad_market_produced_available(m) * room_relabel(grad_agriculture_rainfedareas_production(m), :production, :Market, :produced), :available, :Market, :domesticbalance)) # +
    setconstraint!(house, room_relabel(grad_market_internationalsales_available(m), :available, :Market, :domesticbalance)) # -
    setconstraint!(house, room_relabel(grad_market_regionimports_available(m) * grad_transportation_imported_regionimports(m) +
                                       grad_market_regionexports_available(m) * grad_transportation_imported_regionexports(m), :available, :Market, :domesticbalance)) # +-
    setconstraintoffset!(house, -hall_relabel(constraintoffset_domesticdemand_cropinterest(m), :cropinterest, :Market, :domesticbalance)) # -

    # Clean up

    house.b[isnan(house.b)] = 0
    house.b[!isfinite(house.b)] = 0
    house.f[isnan(house.f)] = 0
    house.f[!isfinite(house.f)] = 0

    ri, ci, vv = findnz(house.A)
    for ii in find(isnan(vv))
        house.A[ri[ii], ci[ii]] = vv[ii]
    end
    for ii in find(!isfinite(vv))
        house.A[ri[ii], ci[ii]] = 1e9
    end

    serialize(open(joinpath(todata, "fullhouse$suffix.jld"), "w"), house)
else
    house = deserialize(open(joinpath(todata, "fullhouse$suffix.jld"), "r"));
end

using MathProgBase
@time sol = linprog(-house.f, house.A, '<', house.b, house.lowers, house.uppers) #1e8 * ones(length(house.uppers)))

# Look at parameter values
varlens = varlengths(m, house.paramcomps, house.parameters)
for ii in 1:length(house.parameters)
    println(house.parameters[ii])
    index1 = sum(varlens[1:ii-1]) + 1
    index2 = sum(varlens[1:ii])

    values = sol.sol[index1:index2][1:min(100, index2 - index1 + 1)]
    if (sum(values .!= 0) == 0)
        println("All zero.")
    else
        println(values)
        println("Sum: $(sum(values))")
    end
end

# Get constraint values
constvalues = house.A * sol.sol

varlens = varlengths(m, house.constcomps, house.constraints)
for ii in 1:length(house.constraints)
    println(house.constraints[ii])
    index1 = sum(varlens[1:ii-1]) + 1
    index2 = sum(varlens[1:ii])

    values = constvalues[index1:index2][1:min(100, index2 - index1 + 1)]
    if (sum(values .!= 0) == 0)
        println("All zero.")
    else
        println(values)
    end
end
