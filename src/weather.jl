using NetCDF

statefips = ncread("../data/VIC_WB.nc", "state_fips")
countyfips = ncread("../data/VIC_WB.nc", "county_fips")
fips = map(fipsnum -> (fipsnum < 10000 ? "0$fipsnum" : "$fipsnum"), statefips * 1000 + countyfips)

function reorderfips(weather, fromfips, tofips)
    result = zeros(length(tofips), size(weather, 2))
    for rr in 1:length(tofips)
        result[rr, :] = weather[fromfips .== tofips[rr], :]
    end

    result
end

# Load data from the water budget
runoff = reorderfips(ncread("../data/VIC_WB.nc", "runoff"), fips, names)
precip_month = reorderfips(ncread("../data/VIC_WB.nc", "precip"), fips, names) # mm / month


