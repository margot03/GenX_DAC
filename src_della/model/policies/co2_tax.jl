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
    co2_tax!(EP::Model, inputs::Dict, setup::Dict)
"""
function co2_tax!(EP::Model, inputs::Dict, setup::Dict)

    println("C02 Tax Module")

    dfGen = inputs["dfGen"]
    SEG = inputs["SEG"]  # Number of lines
    G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones

    ### Expressions ###
    #CO2 Tax
    # Sum CO2 Tax to plant level
    @expression(EP, ePlantCCO2Tax[y = 1:G], sum(inputs["omega"][t] * EP[:eEmissionsByPlant][y, t] for t in 1:T) * inputs["dfCO2Tax"][dfGen[y, :Zone], "CO2Tax"])
    # Sum CO2 Tax to zonal level
    @expression(EP, eZonalCCO2Tax[z = 1:Z], EP[:vZERO] + sum(EP[:ePlantCCO2Tax][y] for y in dfGen[(dfGen[!, :Zone].==z), :R_ID]))
    # Sum CO2 Tax to system level
    @expression(EP, eTotalCCO2Tax, sum(EP[:eZonalCCO2Tax][z] for z in 1:Z))

    if setup["FLECCS"] > 0
        G_F = inputs["G_F"] # Number of FLECCS generator
        dfGen_ccs = inputs["dfGen_ccs"]
        @expression(EP, ePlantCCO2TaxFLECCS[y = 1:G_F], sum(EP[:eEmissionsByPlantFLECCS][y, t] for t in 1:T) * inputs["dfCO2Tax"][dfGen_ccs[y, :Zone], "CO2Tax"])
        # Sum CO2 Tax to zonal level
        @expression(EP, eZonalCCO2TaxFLECCS[z = 1:Z],  sum(EP[:ePlantCCO2TaxFLECCS][y] for y in unique(dfGen_ccs[(dfGen_ccs[!, :Zone].==z), :R_ID])))
        # Sum CO2 Tax to system level
        @expression(EP, eTotalCCO2TaxFLECCS, sum(EP[:eZonalCCO2TaxFLECCS][z] for z in 1:Z))
        add_to_expression!(EP[:eObj],  EP[:eTotalCCO2TaxFLECCS])
    end    



    if setup["DAC"] > 0
        dfDac = inputs["dfDac"]
        DAC_ID  = dfDac[!,:DAC_ID] # Number of DAC
        @expression(EP, ePlantCCO2TaxDAC[y = DAC_ID], sum(EP[:eCO2_DAC_net][y, t] for t in 1:T) * inputs["dfCO2Tax"][dfDac[y, :Zone], "CO2Tax"])
        # Sum CO2 Tax to zonal level
        @expression(EP, eZonalCCO2TaxDAC[z = 1:Z],  sum(EP[:ePlantCCO2TaxDAC][y] for y in unique(dfDac[(dfDac[!, :Zone].==z), :R_ID])))
        # Sum CO2 Tax to system level
        @expression(EP, eTotalCCO2TaxDAC, sum(EP[:eZonalCCO2TaxDAC][z] for z in 1:Z))
        add_to_expression!(EP[:eObj],  EP[:eTotalCCO2TaxDAC])
    end    



    add_to_expression!(EP[:eObj], EP[:eTotalCCO2Tax])

end
