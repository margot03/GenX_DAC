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

@doc raw""" CO2 emissions and CO2 capture"""
function co2!(EP::Model, inputs::Dict, setup::Dict)

    println("C02 Module")

    dfGen = inputs["dfGen"]
    G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones
    #dfDac = inputs["dfDac"]
    #G_DAC = length(collect(skipmissing(dfDac[!,:R_ID])))
    
    scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    dfGen.BECCS = "BECCS" in names(dfGen) ? dfGen.BECCS : zeros(Int, nrow(dfGen))

    ### Expressions ###
    # CO2 emissions from power plants in "Generator_data.csv"
    if setup["CO2Capture"] == 0
        @expression(EP, eEmissionsByPlant[y=1:G, t=1:T], 
            ((EP[:vFuel][y, t] + EP[:eStartFuel][y, t]) * 
                inputs["fuel_CO2"][dfGen[y,:Fuel]]))
    else # setup["CO2Capture"] == 1
        @expression(EP, eEmissionsByPlant[y=1:G, t=1:T],
            ((1-dfGen.BECCS[y]) - dfGen[!, :CO2_Capture_Rate][y]) * 
            ((EP[:vFuel][y, t] + EP[:eStartFuel][y, t]) * 
                inputs["fuel_CO2"][dfGen[y,:Fuel]]))
        # CO2  captured from power plants in "Generator_data.csv"
        @expression(EP, eEmissionsCaptureByPlant[y=1:G, t=1:T],
            (dfGen[!, :CO2_Capture_Rate][y]) * 
            ((EP[:vFuel][y, t] + EP[:eStartFuel][y, t]) * 
                inputs["fuel_CO2"][dfGen[y,:Fuel]]))
        
        @expression(EP, eEmissionsCaptureByPlantYear[y=1:G], 
            sum(inputs["omega"][t] * eEmissionsCaptureByPlant[y, t] 
                for t in 1:T))
        @expression(EP, eEmissionsCaptureByZone[z=1:Z, t=1:T], 
            sum(eEmissionsCaptureByPlant[y, t] 
                for y in dfGen[(dfGen[!, :Zone].==z), :R_ID]))
        @expression(EP, eEmissionsCaptureByZoneYear[z=1:Z], 
            sum(eEmissionsCaptureByPlantYear[y] 
                for y in dfGen[(dfGen[!, :Zone].==z), :R_ID]))
    
    
        # add CO2 sequestration cost to objective function
        # when scale factor is on tCO2/MWh = > kt CO2/GWh
        @expression(EP, ePlantCCO2Sequestration[y=1:G], 
            sum(inputs["omega"][t] * eEmissionsCaptureByPlant[y, t] * 
                dfGen[y, :CO2_Capture_Cost_per_Metric_Ton] /scale_factor for t in 1:T))
    
        @expression(EP, eZonalCCO2Sequestration[z=1:Z], 
            sum(ePlantCCO2Sequestration[y] 
                for y in dfGen[(dfGen[!, :Zone].==z), :R_ID]))
    
        @expression(EP, eTotaleCCO2Sequestration, 
            sum(eZonalCCO2Sequestration[z] for z in 1:Z))
    
        add_to_expression!(EP[:eObj], EP[:eTotaleCCO2Sequestration])
    end

    @expression(EP, eEmissionsByPlantYear[y = 1:G], 
        sum(inputs["omega"][t] * eEmissionsByPlant[y, t] for t in 1:T))
    #@expression(EP, eEmissionsByZone[z = 1:Z, t = 1:T], 
    #    sum(eEmissionsByPlant[y, t] for y in dfGen[(dfGen[!, :Zone].==z), :R_ID]))
    

    #
    @expression(EP, eEmissionsByZoneAll[z=1:Z, t=1:T], 0)
    
    # emissions from generator_data.csv
    @expression(EP, eEmissionsByZone[z = 1:Z, t = 1:T], 
        sum(eEmissionsByPlant[y, t] for y in dfGen[(dfGen[!, :Zone].==z), :R_ID]))

    @expression(EP, eEmissionsPowerPlants, 
        sum(sum(eEmissionsByZone[z, t] for z in 1:Z) for t = 1:T))

    eEmissionsByZoneAll += eEmissionsByZone

    # CO2 emissions from FLECCS 
	if setup["FLECCS"] >= 1
        gen_ccs = inputs["dfGen_ccs"]
	    @expression(EP, eEmissionsByZoneFleccs[z=1:Z, t=1:T], 
        sum(EP[:eEmissionsByPlantFLECCS][y,t] for y in unique(gen_ccs[(gen_ccs[!,:Zone].==z),:R_ID])))
        eEmissionsByZoneAll += eEmissionsByZoneFleccs
    end


    ## when modeling scale factor is on, electricity from MW => GW, EmissionsByPlant = GW * (kt/GW) = kt CO2
    
    if setup["DAC"] == 1
        dfDac = inputs["dfDac"]
	    @expression(EP, eEmissionsByZoneDAC[z=1:Z, t=1:T], 
        sum(EP[:eCO2_DAC_net][y,t] for y in unique(dfDac[(dfDac[!,:Zone].==z),:DAC_ID])))
        eEmissionsByZoneAll += eEmissionsByZoneDAC
    end

    @expression(EP, eEmissionsByZoneYear[z = 1:Z], 
        sum(eEmissionsByZoneAll[z, t] for t in 1:T))

    @expression(EP, eEmissionsTotalZoneYear, 
        sum(eEmissionsByZoneYear[z] for z in 1:Z))

    #@constraint(EP, cNetRemovals, eEmissionsTotalZoneYear <= 84087512 - sum(dfDac[y, :Deployment] for y in 1:G_DAC))

    # CO2 zontal constraint
    if setup["NetZero"] == 1
        @constraint(EP, cNetZeroGrid, eEmissionsTotalZoneYear == 0)
    end

    if setup["CarbonBudget"] != 0
        @constraint(EP, CarbonBudget, eEmissionsTotalZoneYear <= setup["CarbonBudget"])
    end

    return EP

end
