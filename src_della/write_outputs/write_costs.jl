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
	write_costs(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)

Function for writing the costs pertaining to the objective function (fixed, variable O&M etc.).
"""
function write_costs(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	## Cost results
	dfGen = inputs["dfGen"]
	dfDac =  (setup["DAC"] > 0 ? inputs["dfDac"] : 0.0)
	SEG = inputs["SEG"]  # Number of lines
	Z = inputs["Z"]     # Number of zones
	T = inputs["T"]     # Number of time steps (hours)

	dfCost = DataFrame(Costs = ["cTotal", "cFix", "cVar", "cNSE", "cStart", "cUnmetRsv", "cNetworkExp", "cUnmetPolicyPenalty", "cCO2_seq", "cCO2_tax"])
	cVar = value(EP[:eTotalCVarOut])+ value.(EP[:eTotalCFuelOut])+ (!isempty(inputs["STOR_ALL"]) ? value(EP[:eTotalCVarIn]) : 0.0) + (!isempty(inputs["FLEX"]) ? value(EP[:eTotalCVarFlexIn]) : 0.0)+ (dfDac != 0 ? value(EP[:eCTotalVariableDAC]) : 0.0)
	cFix = value(EP[:eTotalCFix]) + (!isempty(inputs["STOR_ALL"]) ? value(EP[:eTotalCFixEnergy]) : 0.0) + (!isempty(inputs["STOR_ASYMMETRIC"]) ? value(EP[:eTotalCFixCharge]) : 0.0)+(dfDac != 0 ? value(EP[:eTotalCFixedDAC]) : 0.0) + (dfDac != 0 ? value(EP[:eTotalCFixedDAC_Energy]) : 0.0)
    
	cCO2_seq = (setup["CO2Capture"] >0) ? (value(EP[:eTotaleCCO2Sequestration]) + (dfDac != 0 ? value(EP[:eCTotalCO2TS]) : 0.0)) : 0
	cCO2_tax = ((setup["CO2Tax"]  > 0)  ?  (value.(EP[:eTotalCCO2Tax]) ) : 0.0) + ((setup["DAC"] > 0 & setup["CO2Tax"]  > 0)  ?  ( value.(EP[:eTotalCCO2TaxDAC])) : 0.0)


	dfCost[!,Symbol("Total")] = [objective_value(EP), cFix, cVar,  value(EP[:eTotalCNSE]), 0.0, 0.0, 0.0, 0.0, cCO2_seq, cCO2_tax]

	if setup["ParameterScale"] == 1
		dfCost.Total *= ModelScalingFactor^2
	end

	if setup["UCommit"]>=1
		dfCost[5,2] = value(EP[:eTotalCStart])
	end

	if setup["Reserves"]==1
		dfCost[6,2] = value(EP[:eTotalCRsvPen])
	end

	if setup["NetworkExpansion"] == 1 && Z > 1
		dfCost[7,2] = value(EP[:eTotalCNetworkExp])
	end

	if haskey(inputs, "dfCapRes_slack")
		dfCost[8,2] += value(EP[:eCTotalCapResSlack])
	end

	if haskey(inputs, "dfESR_slack")
		dfCost[8,2] += value(EP[:eCTotalESRSlack])
	end
	
	if haskey(inputs, "dfCO2Cap_slack")
		dfCost[8,2] += value(EP[:eCTotalCO2CapSlack])
	end
	
	if haskey(inputs, "MinCapPriceCap")
		dfCost[8,2] += value(EP[:eTotalCMinCapSlack])
	end	



	if setup["ParameterScale"] == 1
		dfCost[5,2] *= ModelScalingFactor^2
		dfCost[6,2] *= ModelScalingFactor^2
		dfCost[7,2] *= ModelScalingFactor^2
		dfCost[8,2] *= ModelScalingFactor^2
	end

	for z in 1:Z
		tempCTotal = 0.0
		tempCFix = 0.0
		tempCVar = 0.0
		tempCStart = 0.0
		tempCNSE = 0.0


		Y_ZONE = dfGen[dfGen[!,:Zone].==z,:R_ID]

		if dfDac != 0
			Y_ZONE_DAC = dfDac[dfDac[!,:Zone].==z,:R_ID]
		end

		STOR_ALL_ZONE = intersect(inputs["STOR_ALL"], Y_ZONE)
		STOR_ASYMMETRIC_ZONE = intersect(inputs["STOR_ASYMMETRIC"], Y_ZONE)
		FLEX_ZONE = intersect(inputs["FLEX"], Y_ZONE)
		COMMIT_ZONE = intersect(inputs["COMMIT"], Y_ZONE)

		eCFix = sum(value.(EP[:eCFix][Y_ZONE]))
		eCFix_DAC = ((setup["DAC"]  > 0)  ? sum(value.(EP[:eCFixed_DAC][Y_ZONE_DAC])) : 0)
		tempCFix += eCFix + eCFix_DAC
		tempCTotal += tempCFix

		tempCVar = sum(value.(EP[:eCVar_out][Y_ZONE,:]))
		CVar_DAC = ((setup["DAC"]  > 0)  ? sum(value.(EP[:eCTotalVariableDACT][Y_ZONE_DAC])) : 0 )
		CVar_fuel = sum(value.(EP[:ePlantCFuelOut][Y_ZONE,:]))
		tempCVar += CVar_DAC + CVar_fuel
		tempCTotal += tempCVar 

		tempCO2_seq =  (setup["CO2Capture"]  > 0)  ? sum(value.(EP[:ePlantCCO2Sequestration][Y_ZONE])) : 0
		CCO2seq_DAC = ((setup["DAC"]  > 0)  ? sum(value.(EP[:eCCO2_TS_ByPlant][Y_ZONE_DAC,:])) : 0 )
		tempCO2_seq += CCO2seq_DAC
		tempCTotal += tempCO2_seq

		tempCO2_tax = ((setup["CO2Tax"]  > 0)  ? sum(value.(EP[:ePlantCCO2Tax][Y_ZONE])) : 0 )
		CCO2_DAC = ((setup["CO2Tax"]  & setup["DAC"]  > 0)  ? sum(value.(EP[:ePlantCCO2TaxDAC][Y_ZONE_DAC])) : 0)
		tempCO2_tax += CCO2_DAC
		tempCTotal += tempCO2_tax


		if !isempty(STOR_ALL_ZONE)
			eCVar_in = sum(value.(EP[:eCVar_in][STOR_ALL_ZONE,:]))
			tempCVar += eCVar_in
			eCFixEnergy = sum(value.(EP[:eCFixEnergy][STOR_ALL_ZONE]))
			tempCFix += eCFixEnergy

			tempCTotal += eCVar_in + eCFixEnergy
		end
		if !isempty(STOR_ASYMMETRIC_ZONE)
			eCFixCharge = sum(value.(EP[:eCFixCharge][STOR_ASYMMETRIC_ZONE]))
			tempCFix += eCFixCharge
			tempCTotal += eCFixCharge
		end
		if !isempty(FLEX_ZONE)
			eCVarFlex_in = sum(value.(EP[:eCVarFlex_in][FLEX_ZONE,:]))
			tempCVar += eCVarFlex_in
			tempCTotal += eCVarFlex_in
		end

		if setup["UCommit"] >= 1
			eCStart = sum(value.(EP[:eCStart][COMMIT_ZONE,:]))
			tempCStart += eCStart
			tempCTotal += eCStart
		end

		tempCNSE = sum(value.(EP[:eCNSE][:,:,z]))
		tempCTotal += tempCNSE

		if setup["ParameterScale"] == 1
			tempCTotal *= ModelScalingFactor^2
			tempCFix *= ModelScalingFactor^2
			tempCVar *= ModelScalingFactor^2
			tempCNSE *= ModelScalingFactor^2
			tempCStart *= ModelScalingFactor^2
		end
		dfCost[!,Symbol("Zone$z")] = [tempCTotal, tempCFix, tempCVar, tempCNSE, tempCStart, "-", "-", "-",tempCO2_seq,tempCO2_tax]
	end
	CSV.write(joinpath(path, "costs.csv"), dfCost)
end
