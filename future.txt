room_relabel(grad_agriculture_irrigatedareas_production(m), :production, :Market, :produced) * deriv_market_produced_totalrevenue(m)
PAR x VAR * VAR x 1 -> PAR x 1

gradient(Component, par1, var1)

function gradient(s::Agriculture, p::TotalAreaParameter, v::WaterDeficitVariable)
end


[grad/deriv]_component_parameter_variable -> VAR x PAR
[grad/deriv]_component_variable_parameter -> VAR x PAR

deriv_component_variable2_variable3 * grad_comp_var1_var2 * grad_component_parameter_variable1
 1 x VAR * VAR x PAR -> 1 x PAR


## REDO
function grad_irrigatedareas_totalirrigation(m::Model, s::Agriculture)
    p = s.Parameters

    gen(rr, tt) = min(0., p.water_demand[cc] - p.precipitation[rr, tt]) / 100
    roomdiagonal(m, :Agriculture, :totalirrigation, :irrigatedareas, gen)
end
