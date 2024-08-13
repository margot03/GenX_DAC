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
    write_dac_energy_consumption(path::AbstractString, inputs::Dict, setup::Dict, EP::Model))

Function for writing the energy consumption of DAC in each hour.
"""
function write_dac_energy_consumption(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    # println("hello world")
    # scale factor
    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1
    # Capacity decisions
    dfDac = inputs["dfDac"]
	n = length(dfDac[!,:DAC_ID])

    dfEnergy = DataFrame(Resource = repeat(dfDac[!,:Resource],3), 
                        Zone = repeat(dfDac[!,:Zone],3),
                        MWelec = repeat(["Total_Consumption", "Heat", "Electricity"], inner=n)
    )

    # total energy consumption includes electricity for heat demand (MWelec) and operational demand (eDAC_power, MWelec)
    dac_total_consumption = value.(EP[:eDAC_demand_tot])
    dac_heat = value.(EP[:eDAC_heat_demand_MWelec])
    dac_power = value.(EP[:eDAC_power])

    energyDAC = DataFrame(vcat(dac_total_consumption, dac_heat, dac_power),:auto)

    dfEnergy = hcat(dfEnergy, energyDAC)

    CSV.write(joinpath(path, "dac_energy.csv"), dftranspose(dfEnergy, false), writeheader=false)

end