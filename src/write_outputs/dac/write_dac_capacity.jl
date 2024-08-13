"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

@doc raw"""
	write_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model))

Function for writing the diferent capacities for the different generation technologies (starting capacities or, existing capacities, retired capacities, and new-built capacities).
"""
function write_dac_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    ##############
	# write DAC capacity
	# scale factor
	## ModelScalingFactor = 1e3
	scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1
	# Capacity decisions
	dfDac = inputs["dfDac"]
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	#MultiStage = setup["MultiStage"]
	NEW_CAP = dfDac[dfDac[!,:NEW_CAP] .==1,:R_ID]
	DAC_COMMIT = dfDac[dfDac[!,:DAC_COMMIT] .==1,:R_ID]
	capdDac = zeros(size(dfDac[!,:Resource]))
	for i in  NEW_CAP
		capdDac[i] = value(EP[:vCAP_DAC][i])
	end

	dfCapDac = DataFrame(
		Resource = dfDac[!,:Resource],
		Zone = dfDac[!,:Zone],
		StartCap = dfDac[!,:Existing_Cap_CO2],
		NewCap = capdDac[:],
		EndCap = capdDac[:]
	)

	dfCapDac.StartCap = dfCapDac.StartCap * scale_factor
	dfCapDac.NewCap = dfCapDac.NewCap * scale_factor
	dfCapDac.EndCap = dfCapDac.EndCap * scale_factor

	#dfCapDac = vcat(dfCapDac, total)
	CSV.write(joinpath(path, "capacity_dac.csv"), dfCapDac)
	##############

	##############
	# write DAC energy consumption

	dfDac = inputs["dfDac"]
	n = length(dfDac[!,:DAC_ID])

    dfEnergy = DataFrame(Resource = repeat(dfDac[!,:Resource],3), 
                        Zone = repeat(dfDac[!,:Zone],3),
                        MWelec = repeat(["Total_Consumption", "Heat", "Electricity"], inner=n)
    )

    # total energy consumption includes electricity for heat demand (MWelec) and operational demand (eDAC_power, MWelec)
    dac_total_consumption = value.(EP[:eDAC_demand_tot]).*scale_factor^2
    dac_heat = value.(EP[:eDAC_heat_demand_MWelec]).*scale_factor^2
    dac_power = value.(EP[:eDAC_power]).*scale_factor^2

    energyDAC = DataFrame(vcat(dac_total_consumption, dac_heat, dac_power),:auto)

    dfEnergy = hcat(dfEnergy, energyDAC)

    CSV.write(joinpath(path, "dac_energy.csv"), dftranspose(dfEnergy, false), writeheader=false)

	#######################
	# write DAC cost 

	## Find electricity cost associated with DAC demand (heat, operation)
	## reference files: write_price, write_costs
	## right now, using full year T and no scaling factor so these don't modify cPowerBalance
	prices = dual.(EP[:cPowerBalance])./inputs["omega"]*scale_factor
	price_all_zones = []
	for z in 1:Z
		price_zone = prices[:,z] .* dac_total_consumption[z,:]
		price_zone_total = sum(price_zone[t] for t in 1:T)
		push!(price_all_zones, price_zone_total)
	end
	dfPrices = DataFrame(Costs = ["cElectricity_z$z" for z in 1:Z], Total = price_all_zones)
	
	hp_cap = Array(value.(EP[:vHP_CAP])).*scale_factor^2 ## heat pump capacity per zone, MW electric 
	dfHP = DataFrame(Costs = ["HP_cap_z$z" for z in 1:Z], Total = hp_cap)

	# Find other costs

	dfCost = DataFrame(Costs = ["cTotal", "cFix", "cFix_Energy", "cVar", "cCO2_seq", "cCO2_tax"])
	cVar = (!isempty(inputs["dfDac"]) ? value(EP[:eCTotalVariableDAC]*scale_factor^2) : 0.0)
	cFix = (!isempty(inputs["dfDac"]) ? value(EP[:eTotalCFixedDAC]*scale_factor^2) : 0.0)
	cFix_E = (!isempty(inputs["dfDac"]) ? value(EP[:eTotalCFixedDAC_Energy]*scale_factor) : 0.0)
    cCO2_seq =  (!isempty(inputs["dfDac"]) ? value(EP[:eCTotalCO2TS])*scale_factor^2 : 0.0)
	cCO2_tax =  ((setup["CO2Tax"]  > 0)  ? value.(EP[:eTotalCCO2TaxDAC])*scale_factor^2 : 0)

	## Adding cost of electricity for DAC demand to total variable costs 
	cVar += sum(price_all_zones)
	
	cDacTotal = 0 
	cDacTotal += (cVar + cFix + cFix_E+ cCO2_seq + cCO2_tax)

	dfCost[!,Symbol("Total")] = [cDacTotal, cFix, cFix_E, cVar, cCO2_seq, cCO2_tax]

	# dfCost = vcat(dfCost, dfPrices)
	dfCost = vcat(dfCost, dfPrices, dfHP)

	CSV.write(joinpath(path, "Dac_HP_costs.csv"), dfCost)
	#######################

	
end