using DataFrames
using Mimi

include("world.jl")

water_requirements = Dict("alfalfa" => 1.63961100235402, "otherhay" => 1.63961100235402,
                          "Barley" => 1.18060761343329, "Barley.Winter" => 1.18060761343329,
                          "Maize" => 1.47596435526564,
                          "Sorghum" => 1.1364914374721,
                          "Soybeans" => 1.37599595071683,
                          "Wheat" => 0.684836198198068, "Wheat.Winter" => 0.684836198198068) # in m

type StatisticalAgricultureModel
    intercept::Float64
    interceptse::Float64
    gdds::Float64
    gddsse::Float64
    kdds::Float64
    kddsse::Float64
    wreq::Float64
    wreqse::Float64
end

function StatisticalAgricultureModel(df::DataFrame, filter::Symbol, fvalue::UTF8String)
    interceptrow = (df[filter] .== fvalue) & (df[:coef] .== "intercept")
    gddsrow = (df[filter] .== fvalue) & (df[:coef] .== "gdds")
    kddsrow = (df[filter] .== fvalue) & (df[:coef] .== "kdds")
    wreqrow = (df[filter] .== fvalue) & (df[:coef] .== "wreq")

    if sum(interceptrow) == 1
        intercept = df[interceptrow, :mean]
        interceptse = df[interceptrow, :serr]
    else
        intercept = 0
        interceptse = 0
    end

    gdds = df[gddsrow, :mean]
    gddsse = df[gddsrow, :serr]
    kdds = df[kddsrow, :mean]
    kddsse = df[kddsrow, :serr]
    wreq = df[wreqrow, :mean]
    wreqse = df[wreqrow, :serr]

    StatisticalAgricultureModel(intercept, interceptse, gdds, gddsse, kdds, kddsse, wreq, wreqse)
end

function gaussianpool(mean1, sdev1, mean2, sdev2)
    (mean1 / sdev1^2 + mean2 / sdev2^2) / (1 / sdev1^2 + 1 / sdev2^2)
end

# Prepare all the agricultural models
agmodels = Dict{UTF8String, Dict{Int64, StatisticalAgricultureModel}}() # {crop: {fips: model}}
nationals = readtable("../data/nationals.csv")
for crop in crops
    agmodels[crop] = Dict{Int64, StatisticalAgricultureModel}()

    # Create the national model
    national = StatisticalAgricultureModel(nationals, :crop, crop)
    counties = readtable("../data/unpooled-$crop.csv")
    for fips in unique(counties[:fips])
        county = StatisticalAgricultureModel(counties, :fips, fips)
        # Construct a pooled combination
        gdds, gddsse = gaussianpool(national.gdds, national.gddsse, county.gdds, county.gddsse)
        kdds, kddsse = gaussianpool(national.kdds, national.kddsse, county.kdds, county.kddsse)
        wreq, wreqse = gaussianpool(national.wreq, national.wreqse, county.wreq, county.wreqse)
        agmodel = StatisticalAgricultureModel(gdds, gddsse, kdds, kddsse, wreq, wreqse)
        agmodels[crop][fips] = agmodel
    end
end

@defcomp Agriculture begin
    regions = Index()
    crops = Index()

    # Optimized
    # Land area appropriated to each crop, irrigated to full demand
    irrigatedareas = Parameter(index=[regions, crops, time])
    rainfedareas = Parameter(index=[regions, crops, time])

    # Internal
    # Yield base: combination of GDDs, KDDs, and intercept
    logirrigatedyield = Parameter(index=[regions, crops, time])

    # Coefficient on the effects of water deficits
    deficit_coeff = Parameter(index=[regions, crops])

    # Water requirements per unit area, in mm
    water_demand = Parameter(index=[crops])

    # Precipitation water per unit area, in mm
    precipitation = Parameter(index=[regions, time])

    # Computed
    # Land area appropriated to each crop
    totalareas = Variable(index=[regions, crops, time])

    # Deficit for any unirrigated areas, in mm
    water_deficit = Variable(index[regions, crops, time])

    # Total irrigation water (1000 m^3)
    totalirrigation = Variable(index=[regions, time])

    # Yield per hectare for rainfed (irrigated has irrigatedyield)
    lograinfedyield = Variable(index=[regions, crops, time])

    # Total production
    production = Variable(index=[regions, crops, time])
end

"""Simulates crop yields as a Cobb-Douglas model of water and energy."""
function timestep(state::agriculture, tt::Int)
    v = s.Variables
    p = s.Parameters
    d = s.Dimensions

    for rr in d.regions
        totalirrigation = 0.
        for cc in d.crops
            v.totalareas[rr, cc, tt] = p.irrigatedareas[rr, cc, tt] + p.rainfedareas[rr, cc, tt]

            # Calculate deficit by crop, for unirrigated areas
            v.water_deficit[rr, cc, tt] = p.water_demand[rr, cc, tt] - p.precipitation[rr, tt]

            # Calculate irrigation water, summed across all crops: 1 mm * Ha^2 = 10 m^3
            totalirrigation += v.water_deficit[rr, cc, tt] * p.irrigatedareas[rr, cc, tt] / 100

            # Calculate rainfed yield
            v.lograinfedyield[rr, cc, tt] = p.logirrigatedyield[rr, cc, tt] + p.deficit_coeff[rr, cc] * v.water_deficit[rr, cc, tt]

            # Calculate total production
            v.production[rr, cc, tt] = exp(p.logirrigatedyield[rr, cc, tt]) * p.irrigatedareas[rr, cc, tt] + exp(p.lograinfedyield[rr, cc, tt]) * p.rainfedareas[rr, cc, tt]
        end

        v.totalirrigation[rr, cc] = totalirrigation
    end
end

function initagriculture(m::Model)
    # precip_month loaded by weather.jl

    # Currently summing over all months
    precip_year = zeros(nrow(precip_month), ncol(precip_month) / 12)
    for year in 1:ncol(precip_month) / 12
        allcounties = zeros(nrow(precip_month))
        for month in 1:12
            allcounties += precip_month[:, (year - 1) * 12 + month]
        end

        precip_year[:, year] = allcounties
    end

    # Match up values by FIPS
    logirrigatedyield = zeros(numregions, numcrops, numsteps)
    deficit_coeff = zeros(numregions, numcrops)
    for rr in 1:numcounties
        for cc in 1:numcrops
            thismodel = agmodels[crops[cc]][names[rr]]
            logirrigatedyield[rr, cc, :] = repmat([thismodel.intercept], numsteps)
            deficit_coeff[rr, cc] = thismodel.wreq
        end
    end

    water_demand = zeros(numcrops)
    for cc in 1:numcrops
        water_demand[cc] = water_requirements[crops[cc]] * 1000
    end

    agriculture = addcomponent(m, Agriculture)

    agriculture[:logirrigatedyield] = repmat(logirrigatedyield, numsteps)
    agriculture[:deficit_coeff] = deficit_coeff
    agriculture[:water_demand] = water_demand
    agriculture[:precipitation] = precip_year
end
