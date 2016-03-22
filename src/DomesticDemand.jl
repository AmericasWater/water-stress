# The domestic demand component

using Mimi
using DataFrames

populations = readtable("../../../data/county-pops.csv", eltypes=[Int64, UTF8String, UTF8String, Int64, Float64]);

function getpopulation(fips, year)
    pop = populations[(populations[:FIPS] .== parse(Int64, fips)) & (populations[:year] .== year), :population]
    if length(pop) != 1
        NA
    else
        pop[1]
    end
end

@defcomp DomesticDemand begin
    regions = Index()
    crops = Index()

    # Internal
    # Resource demands
    population = Parameter(index=[regions, time])

    waterdemandperperson = Parameter()
    cropdemandperperson = Parameter(index=[crops])

    # Resource surplus over (or below) demand
    waterdemand = Variable(index=[regions, time])
    cropdemand = Variable(index=[regions, crops, time])
end

"""
Compute the `surplus` as `available` - `demand`.
"""
function timestep(c::DomesticDemand, tt::Int)
    v = c.Variables
    p = c.Parameters
    d = c.Dimensions

    for rr in d.regions
        v.waterdemand[rr, tt] = p.population[rr, tt] * p.waterdemandperperson

        for cc in d.crops
            v.cropdemand[rr, cc, tt] = p.population[rr, tt] * p.cropdemandperperson[cc]
        end
    end
end

"""
Add a domesticdemand component to the model.
"""
function initdomesticdemand(m::Model, years)
    domesticdemand = addcomponent(m, DomesticDemand)

    # Blue water from http://waterfootprint.org/media/downloads/Hoekstra_and_Chapagain_2006.pdf
    domesticdemand[:waterdemandperperson] = 575 * 365.25 * .001 # m^3 / yr
    domesticdemand[:cropdemandperperson] = repmat(100., length(crops), 1)

    allpops = Matrix{Float64}(m.indices_counts[:regions], length(years))
    totalpop = 0
    for tt in 1:length(years)
        year = years[tt]
        for ii in 1:m.indices_counts[:regions]
            fips = m.indices_values[:regions][ii]
            pop = getpopulation(fips, year)
            if isna(pop) && mod(year, 10) != 0
                # Estimate from decade
                pop0 = getpopulation(fips, div(year, 10) * 10)
                pop1 = getpopulation(fips, (div(year, 10) + 1) * 10)
                pop = pop0 * (1 - mod(year, 10) / 10) + pop1 * mod(year, 10) / 10
            end
            if isna(pop)
                pop = 0.
            end
            allpops[ii, tt] = pop
            totalpop += pop
        end
    end

    domesticdemand[:population] = allpops

    domesticdemand
end

