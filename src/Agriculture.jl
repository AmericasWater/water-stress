using Mimi

@defcomp Agriculture begin
    regions = Index()
    crops = Index()

    # Yield base: combination of GDDs, KDDs, and intercept
    yield0 = Parameter(index=[regions, crops, time])

    # Coefficient on the effects of water deficits
    deficit_coeff = Paramter(index=[regions, crops])

    # Water requirements per unit area, in mm
    water_demand = Parameter(index=[crops])

    # Precipitation water per unit area, in mm
    precipitation = Parameter(index=[regions, time])

    # Irrigation water per unit area, in mm-- OPTIMIZED
    irrigation = Parameter(index=[regions, crops, time])

    # Land area appropriated to each crop-- OPTIMIZED
    areas = Parameter(index=[regions, crops, time])

    # Deficit after irrigation
    water_deficit = Variable(index[regions, crops, time])

    # Yield per hectare
    yield = Variable(index=[regions, crops, time])

    # Total production
    production = Variable(index=[regions, crops, time])
end

"""Simulates crop yields as a Cobb-Douglas model of water and energy."""
function timestep(state::agriculture, tt::Int)
    v = s.Variables
    p = s.Parameters
    d = s.Dimensions

    for rr in d.regions
        for cc in d.crops
            # Calculate deficit by crop
            v.water_deficit[rr, cc, tt] = p.water_demand[rr, cc, tt] - p.precipitation[rr, tt] - p.irrigation[rr, cc, tt]

            # Calculate total yield
            v.yield[rr, cc, tt] = p.yield0[rr, cc, tt] + p.deficit_coeff[rr, cc] * v.water_deficit[rr, cc, tt]

            # Calculate total production
            v.production[rr, cc, tt] = p.yield[rr, cc, tt] * p.areas[rr, cc, tt]
        end
    end
end
