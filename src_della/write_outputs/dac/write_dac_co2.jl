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
	write_power(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)

Function for writing the different values of co2 flow
"""
function write_dac_co2(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)


	dfDac = inputs["dfDac"]
	DAC_ID = dfDac[!,:DAC_ID]
	n = length(DAC_ID)

	dfCO2 = DataFrame(Resource = repeat(dfDac[!,:Resource],2), 
					Zone = repeat(dfDac[!,:Zone],2), 
					# CO2 = repeat(["CO2 heat", "CO2 gross", "CO2 net"], inner=n)
					CO2 = repeat(["CO2 gross", "CO2 net"], inner=n)
	)


	CO2_gross = -value.(EP[:vCO2_DAC])
	CO2_net = value.(EP[:eCO2_DAC_net])

	if setup["ParameterScale"] == 1
		CO2_gross *= ModelScalingFactor
		CO2_net *= ModelScalingFactor
	end

	# nuclear_penalty_from_dac = value.(EP[:ePowerReductionnuclearheat])
	# CSV.write(joinpath(path, "dac_nuke_penalty.csv"), DataFrame(nuclear_penalty_from_dac, :auto), writeheader=false)

	# CO2 = DataFrame(vcat(CO2_heat,CO2_gross, CO2_net),:auto)
	CO2 = DataFrame(vcat(CO2_gross, CO2_net),:auto)

	#dfPower.AnnualSum .= power * inputs["omega"]

	dfCO2= hcat(dfCO2, CO2)

	CSV.write(joinpath(path, "dac_co2.csv"), dftranspose(dfCO2, false), writeheader=false)

	#total_emissions = value.(EP[:eEmissionsTotalZoneYear])

	return dfCO2
end
