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
	load_dac_data!(setup::Dict, path::AbstractString, inputs::Dict)

Read input data related to dac
"""

function load_dac_data!(setup::Dict, path::AbstractString, inputs::Dict)
    filename ="dac.csv"

    dac_in = load_dataframe(joinpath(path, filename))

	# Add Resource IDs after reading to prevent user errors - THIS CODE HAS BEEN CHANGED
	dac_in[!,:DAC_ID] = dac_in.R_ID

	# Store DataFrame of dac input data for use in model
	inputs["dfDac"] = dac_in

    # Number of resources
	#inputs["D"] = length(collect(skipmissing(dac_in[!,:R_ID])))

	# Set indices for internal use
	#D = inputs["D"]   # Number of DAC resources 

    # Zones resources are located in
	#zones = collect(skipmissing(dac_in[!,:Zone][1:inputs["D"]]))
	# Resource identifiers by zone (just zones in resource order + resource and zone concatenated)
	#inputs["D_ZONES"] = zones
	#inputs["dfDac"][!,:DAC_ZONES] = dac_in["RESOURCES"] .* "_z" .* string.(zones)

    # the capacity of DAC is expressed as t CO2/h 
    # the unit of fixed cost is $/t CO2
    # the unit of non-fuel VOM is $/t CO2
    # the unit of heat consumption is MMBTU/t CO2
    # the unit of electricity consumption is MWh/t CO2

    if setup["ParameterScale"] == 1  # Parameter scaling turned on - adjust values of subset of parameter values
        # keep the unit of CO2 to be t through out of the calculations...
        # convert the cost from $/t CO2/h to million $/t CO2/h
        # I dont think we have exiting capacity for DAC..so I will just do the conversion for Fixed Cost, VOM, heat, electricity consumption
        inputs["dfDac"][!, :Fix_Cost_per_tCO2_yr] = dac_in[!, :Fix_Cost_per_tCO2_yr]/ModelScalingFactor #0809 

        inputs["dfDac"][!, :Var_OM_Cost_per_tCO2] = dac_in[!, :Var_OM_Cost_per_tCO2]/ModelScalingFactor #0809 

        inputs["dfDac"][!, :CO2_Transport_Storage_Per_t] = dac_in[!, :CO2_Transport_Storage_Per_t]/ModelScalingFactor
        
        inputs["dfDac"][!, :Electricity_MWh_per_CO2_metric_ton] = dac_in[!, :Electricity_MWh_per_CO2_metric_ton]/ModelScalingFactor
        
        inputs["dfDac"][!, :Deployment] = dac_in[!, :Deployment]/ModelScalingFactor

    end

    println(filename * " Successfully Read!")


    filename ="dac_params.csv"

    df_dac_params = load_dataframe(joinpath(path, filename))

    dac_params = Dict{AbstractString, Float64}()

    for i = 1:length(df_dac_params.key)
        dac_params[df_dac_params.key[i]] = df_dac_params.value[i]
    end

	# Store DataFrame of dac input data for use in model
	inputs["Dac_params"] = dac_params

    # heat extraction per percent steam diversion is in MW thermal

    println(filename * " Successfully Read!")


    ## Heat Pump Params loading
    filename ="heat_pump_params.csv"

    hp_in = load_dataframe(joinpath(path, filename))

    hp_params = Dict{AbstractString, Float64}()

    for i = 1:length(hp_in.key)
        hp_params[hp_in.key[i]] = hp_in.value[i]
    end

	# Store DataFrame of dac input data for use in model
	inputs["HP_params"] = hp_params

    if setup["ParameterScale"] == 1
        inputs["HP_params"].heatpump_cost_perMW = inputs["HP_params"].heatpump_cost_perMW / ModelScalingFactor
    end
    
    println(filename * " Successfully Read!")

   # return inputs
end




