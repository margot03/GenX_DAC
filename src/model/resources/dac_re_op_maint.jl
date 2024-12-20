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
	DAC RE
"""

function dac_re!(EP::Model, inputs::Dict, setup::Dict)

	println("DAC module")

    #######################
    # case details
    additional_flag = setup["Additional"] # Uses custom file from dac_additional_resources_dac.jl containing R_IDs of clean, additional gens. False: Uses R_IDs for clean generators in generator input data.
    matching = setup["Matching"] # Uses custom file from dac_additional_resources_dac.jl containing R_IDs of clean, additional gens. False: Uses R_IDs for clean generators in generator input data.
    policy_case = lowercase(setup["PolicyCase"]) # "hourly", "annual", "emissions", "none" (basecase, DAC no policy)
    overproc_penalty = 5
    #######################

    # Parameter Scaling factor
	## ModelScalingFactor = 1e3
	scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1
    omega = inputs["omega"] ## has to deal with subperiods and weights for time sampling

    T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

    hours_per_subperiod = inputs["hours_per_subperiod"]

    ## Heat Pump params data
    hp_params = inputs["HP_params"]
    heat_pump_coeff_perf = hp_params["coeff_performance"] # MW thermal / MW    ex. 2 = 2 MW thermal output / 1 MW electricity input 
    heatpump_cost_perMW = hp_params["heatpump_cost_perMW"] # $/MW electric, annualized (investment fixed cost)

    # load dac input data and coupling parameters
    dfDac = inputs["dfDac"]
    Dac_params = inputs["Dac_params"]

    maintenance_duration = dfDac[!,:Maintenance_duration][1] #number of hours for maintenance downtime

    maint_freq_years= dfDac[!,:Maint_freq_years][1] # once per year is default

    maint_begin_cadence = dfDac[!,:Maint_begin_cadence][1] # hour increments in which maintenance can begin, you don't want it to be able to begin in every hour of the year, default is every 100 hours

    maintenance_begin_hours = 1:maint_begin_cadence:T # creates the range 

    println(maintenance_duration)
    println(maint_freq_years)
    println(maint_begin_cadence)
    println(maintenance_begin_hours)


    # dfDac[:,:CF] = ones(length(dfDac[:,:CF])) .* 0.98
    # println("cf:", dfDac[:,:CF])
    # println(dfDac[:,:CF])
    # println("z", dfDac[dfDac.Deployment .== 1e7, :Zone])

    DAC_ID = dfDac[!,:DAC_ID]  # collect ids
    G_DAC = length(collect(skipmissing(dfDac[!,:R_ID])))  # number of DAC types
    
    re_dac = dfDac[dfDac.DAC_heat_resource.=="RE", :R_ID] #subset of electricity (RE) dac facilities

    dfGen = inputs["dfGen"]   #generator data

    duplicated_gens = []
    println("duplicated generators R_IDs")
    # G_Dup = []
    # for col in names(dfGen)
    #     if occursin("MaxCapTag_", col)
    #         gen_pair = dfGen[dfGen[!,Symbol(col)].==1, :]
    #         org_res = gen_pair[1,:]
    #         # push!(G_Dup, org_res[1,:R_ID])
    #         push!(G_Dup, org_res.R_ID)

    #         dup_res = gen_pair[2,:]
    #         # push!(duplicated_gens, dup_res[1,:R_ID]) 
    #         push!(duplicated_gens, dup_res.R_ID)
    #         println(col, ": ", org_res.R_ID, ", ", dup_res.R_ID)
    #     end
    # end
    
    duplicated_gens = [dfGen[dfGen.Resource.==a, :R_ID][1] for a in dfGen.Resource if occursin("RPSH_1", a)]
    
    # maxcap = 1:inputs["NumberOfMaxCapReqs"]], sum(EP[:eTotalCap][y] for y in dfGen[dfGen[!, Symbol("MaxCapTag_$maxcap")] .== 1

    # println("original gens:", G_Dup)
    println("duplicated_gens:", duplicated_gens)

    if additional_flag == "add"
        println("additional gens")
        ## ADDITIONAL RESOURCES, already filtered for new build, clean, and additional logic
        ad_df = load_dataframe("/Users/margotadam/Documents/GitHub/DAC/GenX/dac_additional_resources_29_26.csv")
        eligible_gens = ad_df.R_ID

    elseif additional_flag == "incr"
        # incremental: filtering for new build
        # println("incremental gens")
        new_build_eligible_gens = inputs["NEW_CAP"] # R_IDs of all gens eligible for new build
        # vre_gens = intersect(inputs["VRE"], new_build_eligible_gens) # all VRE gens (R_IDs) eligible for new build
        vre_gens = intersect(inputs["Wind"], new_build_eligible_gens)
        # hydro_gens = intersect(inputs["HYDRO_RES"], new_build_eligible_gens) # HYDRO_RES
        # geothermal_gens = intersect(inputs["geothermal"], new_build_eligible_gens) # all geothermal gens eligible for new build
        nuclear_gens = intersect(inputs["nuclear"], new_build_eligible_gens)
        eligible_gens = [vre_gens; nuclear_gens] #; hydro_gens R_IDs of all clean gens eligible for new build
        println("incremental gens: ", eligible_gens)
        eligible_gens = intersect(duplicated_gens, eligible_gens)

    elseif additional_flag == "exist"
        # existing only
        println("existing gens")
        existing_gens = dfGen[dfGen.Existing_Cap_MW.!=0,:R_ID]
        vre_gens = intersect(inputs["VRE"], existing_gens) # all VRE gens (R_IDs) eligible for new build
        hydro_gens = intersect(inputs["HYDRO_RES"], existing_gens) # HYDRO_RES

        small_hydro_gens = [12,23,31]

        # geothermal_gens = intersect(inputs["geothermal"], new_build_eligible_gens) # all geothermal gens eligible for new build
        nuclear_gens = intersect(inputs["nuclear"], existing_gens)
        eligible_gens = [vre_gens; hydro_gens; small_hydro_gens; nuclear_gens] # R_IDs of all clean gens eligible for new build
        eligible_gens = intersect(duplicated_gens, eligible_gens)

    elseif additional_flag == "all"
        # existing + new build, non-incremental 
        println("all gens")
        vre_gens = inputs["VRE"] # all VRE gens (R_IDs) eligible for new build
        # geothermal_gens = inputs["geothermal"]
        nuclear_gens = inputs["nuclear"]
        eligible_gens = [vre_gens; nuclear_gens] # R_IDs of all clean gens
        eligible_gens = intersect(duplicated_gens, eligible_gens)
    else 
        # noPolicy
        println("noPolicy")
        eligible_gens = dfGen.R_ID
        eligible_gens = intersect(duplicated_gens, eligible_gens)
    end

    # inputs["MinCapReq"] .* dfGen.R_ID --> this should produce 0's for non-eligible mincapreq and R_ID for eligible mincapreq
    # dfGen[dfGen[!, Symbol("MinCapTag_$mincap")]
    # we want to exclude all mincapreq generators, since we care about additionality being investment decisions so if these are required to gen a certain amount then we can ignore them
    req_RIDs = []
    if setup["MinCapReq"] == 1
        num_reqs = inputs["NumberOfMinCapReqs"]
        
        # dfGen[dfGen[!, Symbol("MinCapTag_$mincap")] .== 1, :R_ID]
        for mincap in 1:num_reqs
            append!(req_RIDs, [y for y in dfGen[dfGen[!, Symbol("MinCapTag_$mincap")] .== 1, :R_ID]])
        end
    end

    println(eligible_gens)
    # println("min cap req R_IDs:")
    # println(req_RIDs)
    eligible_gens = setdiff(eligible_gens, req_RIDs)
    println("eligible gens:")
    println(eligible_gens)
    inputs["dac_eligible_gens"] = eligible_gens

    MW_to_GJ_conversion = Dac_params["MW_to_GJ_conversion"]   #3.6
    
    #decision variables
    # vCO2_DAC: the amount of hourly capture by a DAC facility, metric ton CO2/h. 
    # vCAP_DAC: the ANNUAL removal capacity of a DAC facility, metric ton CO2
    # vHP_heat: the hourly heat output by heat pumps for each DAC facility, MW thermal 
    # vHP_CAP: maximum heat pump capacity for each DAC facility, MW electric 
    
    ## MW thermal --> difference is that MW is electricity consumed, MW thermal is heat output (takes into account coeff of performance!!!)
    ## assuming all DAC_IDs are electric DACs, in future can replace DAC_ID with DAC_RE_ID for HP variables
    @variables(EP, begin
        vCO2_DAC[y in DAC_ID,t = 1:T] >= 0
        vCAP_DAC[y in DAC_ID] >= 0
        vHP_heat[y in DAC_ID,t = 1:T] >= 0 
        vHP_CAP[y in DAC_ID] >= 0  
    end)

    # heat consumption from heat resources, GJ / hour 
    @expression(EP, eDAC_heat_consumption[y in DAC_ID, t = 1:T], vCO2_DAC[y,t] * dfDac[y,:Heat_GJ_per_CO2_metric_ton]) 
    ## GJ / hour = tCO2 removed / hour * GJ / tCO2

    # the electricity consumption for DAC, MWh/t CO2
    @expression(EP, eDAC_power[y in DAC_ID, t = 1:T], vCO2_DAC[y,t] * dfDac[y,:Electricity_MWh_per_CO2_metric_ton])
    ## tCO2 captured / hour * MWh / tCO2 = MW

    ## the heat demand (in MW electric) for DAC facility y in hour t, converted from GJ to MW thermal to MW electric, MWelec / hr  
    @expression(EP, eDAC_heat_demand_MWelec[y in DAC_ID, t = 1:T], EP[:eDAC_heat_consumption][y,t]*(1/MW_to_GJ_conversion)*(1/heat_pump_coeff_perf) / scale_factor)

    ## DAC total demand, annual demand for policy constraints (hourly, annual, emissions)
    @expression(EP, eDAC_demand_tot[z = 1:Z, t = 1:T], sum((eDAC_heat_demand_MWelec[y,t] + EP[:eDAC_power][y, t]) for y in dfDac[dfDac[!,:Zone].==z,:R_ID]))
    @expression(EP, eDAC_Annual_demand[z = 1:Z], sum(eDAC_demand_tot[z,t] * omega[t] for t in 1:T))
    
    ## Clean generation (RE) total generation, annual generation for policy constraints
    @expression(EP, eRE_output[z = 1:Z, t=1:T], sum(EP[:vP][y, t] for y in intersect(eligible_gens, dfGen[dfGen[!,:Zone].==z,:R_ID])))         
    @expression(EP, eRE_Annual_gen[z = 1:Z], sum(eRE_output[z,t] * omega[t] for t in 1:T)) # should be down here!!!

    
    for z in 1:Z
        # all DAC demand in this zone that is of type RE
        dac_z_re = intersect(re_dac, dfDac[dfDac[!,:Zone].==z,:R_ID])

        if policy_case == "hourly"
            println("hourly matching: ", z)
            
            @constraints(EP, begin
                [t in 1:T], eDAC_demand_tot[z,t] * (matching / 100) <= eRE_output[z,t]
            end)
            ## MWelec for the heatpumps + DAC power MWelec <= MWelec from RE supply

        elseif policy_case == "annual"
            println("annual matching: ", z)
            
            ## all DAC demand must be equal to current generation of all RE gens
            @constraint(EP, eDAC_Annual_demand[z] * (matching / 100) <= eRE_Annual_gen[z])
            ## MWelec for the heatpumps + DAC power MWelec <= MWelec from RE supply
            ## must be <= for each zone, otherwise could some some generation that's undeliverable meeting the dac demand

        else
            println("noPolicy")
        end

        ## Heat Pump Capacity Constraints

        ## HP heat output, heat pump output MW thermal == DAC heat demand, MW thermal, for t in 1:t
        ## convert MWelec to MW thermal with COP factor
        @constraints(EP, begin
            [y in dac_z_re, t = 1:T], vHP_heat[y,t] == eDAC_heat_demand_MWelec[y,t]*(heat_pump_coeff_perf)
        end)

        ## HP heat output (MW thermal)(COP for MW electric) <= HP capacity (MW electric)
        @constraints(EP, begin
            [y in dac_z_re, t = 1:T], vHP_heat[y,t] <= vHP_CAP[y]*(heat_pump_coeff_perf)
        end)

        ## DAC heat demand for hour t <= max HP capacity (implied by constraints above)

    end

    # # annual excess procurement
    # @expression(EP, eAnnualExcessProc, sum(eRE_Annual_gen[z] - eDAC_Annual_demand[z] for z in 1:Z))
    # # add the cost of overprocurement (annual excess procurement * cost $/MWh) to the objective function 
    # EP[:eObj] += eAnnualExcessProc * overproc_penalty

    # excess procurement: the difference between the eligible generaiton and the DAC demand in each hour
    @expression(EP, eExcessProcurement[t=1:T,z=1:Z], eRE_output[z,t] - eDAC_demand_tot[z,t]) # could we filter this somehow so we only find the excess procurement if DAC demand < RE output in that hour? 
    # EP[:ePowerBalance] = EP[:ePowerBalance] - eExcessProcurement #EP[:ePowerBalance] has dimensions T x Z
    # @expression(EP, ePositiveExcess, sum(eZonalExcess[t] for t in 1:T if value.(EP[:eZonalExcess][t]) >= 0)) 
    # --> can't find out if each hour is positive or negative because the optimization hasn't been solved yet! each hour is still unknown
    
    # @expression(EP, eExcessProcurement[t=1:T, z=1:Z], eRE_output[z,t] - eDAC_demand_tot[z,t]) # could we filter this somehow so we only find the excess procurement if DAC demand < RE output in that hour? 
    # @expression(EP, eZonalExcess[t=1:T], sum(eExcessProcurement[t,z] for z in 1:Z))
    # @constraint(EP, cPositiveExcess[t=1:T], eZonalExcess >= 0)
    
    # for t in 1:T
    #     excess_proc_t = sum(eExcessProcurement[t,z] for z in 1:Z)
    #     if excess_proc_t >= 0
    #         positive_costs += excess_proc_t
    #     end
    # end

    # @variable(EP, vEXCESS[t=1:T] >= 0)
    # @constraint(EP, cExcessProc[t=1:T], vEXCESS[t] >= sum(eRE_output[z,t] - eDAC_demand_tot[z,t] for z in 1:Z))
    # @expression(EP, eTotalExcess, sum(vEXCESS[t] for t in 1:T))
    
    # EP[:eObj] += eTotalExcess * overproc_penalty
    # the power used for DAC and HP must also go into a power balance equation
    ## heat pump electricity demand, MW: heat pump heat output (MW thermal) / heat_pump_coeff_perf (MW thermal / MW)
    # things in power balance epxression should be scaled by omega 
	@expression(EP, ePowerBalanceDAC[t=1:T, z=1:Z], sum((eDAC_power[y,t] + vHP_heat[y,t]*(1/heat_pump_coeff_perf))*omega[t] for y in intersect(dfDac[dfDac[!,:Zone].==z,:][!,:DAC_ID])))
    EP[:ePowerBalance] = EP[:ePowerBalance] - ePowerBalanceDAC
    
    # duplicated gens + original gens < max capacity
    # constraint: sum(dfGen[R_ID == duplicated gen ID (should be list of two)]) <= dfGen[R_ID == duplicated gen ID][1] 
    # @constraint(EP, cDuplicatedGens[t=1:T], [EP[:vP][y,t] for y in G_Dup] .+ [EP[:vP][y,t] for y in duplicated_gens] .<= [dfGen[dfGen.R_ID.==r, :Max_Cap_MW] for r in G_Dup])
    # @constraint(EP, cDuplicatedGens[t=1:T], [EP[:vP][y,t] for y in G_Dup] .<= [dfGen[dfGen.R_ID.==r, :Max_Cap_MW] for r in G_Dup])
    # @expression(EP, eOriginalGenerator_Gen[y in G_Dup, t=1:T], EP[:vP][y,t])
    # @expression(EP, eDuplicatedGenerator_Gen[y in duplicated_gens, t=1:T], EP[:vP][y,t])
    # @expression(EP, eMaxCap_OriginalGen[y in G_Dup], dfGen[dfGen.R_ID.==y, :Max_Cap_MW])
    # # @constraint(EP, cDuplicatedGens[y=1:length(G_Dup), t=1:T], eOriginalGenerator_Gen[y,t] + eDuplicatedGenerator_Gen[y,t] <= eMaxCap_OriginalGen[y])
    # @constraint(EP, cDuplicatedGens[t=1:T], eOriginalGenerator_Gen[:,t] .+ eDuplicatedGenerator_Gen[:,t] .<= eMaxCap_OriginalGen)
    # println("duplicated cap constraint")
    # println([dfGen[dfGen.R_ID.==r, :Max_Cap_MW] for r in G_Dup])

    #---------------------------------- add up cost ---------------------------------------
    # Fixed Cost
    # Combine CAPEX and FOM into annualized Fixed Cost
    # Fixed cost for a DAC y 
	@expression(EP, eCFixed_DAC[y in DAC_ID], dfDac[y,:Fix_Cost_per_tCO2_yr] * vCAP_DAC[y])   #Annualized CAPEX cost 
	# total fixed costs for all the DAC
	@expression(EP, eTotalCFixedDAC, sum(eCFixed_DAC[y] for y in DAC_ID))
	EP[:eObj] += eTotalCFixedDAC

    # Fixed DAC cost for energy (this covers capex of heat pump, and other energy related fixed costs for the DAC plant itself)
    ## should eCFixed_DAC_Energy be scaled by factor of 1/heat_pump_coeff_perf? Both heatpump_cost_perMW and vHP_CAP are in terms of MW-electric, not MW-thermal 
    @expression(EP, eCFixed_DAC_Energy[y in DAC_ID], heatpump_cost_perMW*vHP_CAP[y] + dfDac[y, :Energy_Fix_Cost_per_yr]) # heatpump_cost_perMW*vHP_CAP[y]*(1/heat_pump_coeff_perf)
    # total fixed costs for all the DAC energy demand 
	@expression(EP, eTotalCFixedDAC_Energy, sum(eCFixed_DAC_Energy[y] for y in DAC_ID))
	EP[:eObj] += eTotalCFixedDAC_Energy

    # Variable cost
    # the total variable cost (heat cost + non-fuel vom cost) for DAC y at time t, $ = $/t CO2 * CO2 captured in hour t
    @expression(EP, eCDAC_Variable[y in DAC_ID, t = 1:T], (dfDac[y,:Var_OM_Cost_per_tCO2])*vCO2_DAC[y,t])  
    # Cost associated with variable costs for DAC for the whole year
    @expression(EP, eCTotalVariableDACT[y in DAC_ID], sum(omega[t] * eCDAC_Variable[y,t] for t in 1:T ))
    # Total variable cost for all DAC facilities across the year
    @expression(EP, eCTotalVariableDAC, sum(eCTotalVariableDACT[y] for y in DAC_ID ))

    EP[:eObj] += eCTotalVariableDAC

    DAC_COMMIT = dfDac[dfDac[!,:DAC_COMMIT] .== 1, :R_ID]
    if setup["UCommit"] > 0 && !isempty(DAC_COMMIT)
        @variables(EP, begin
            vCOMMIT_DAC[y in DAC_COMMIT, t=1:T] >= 0 # #number of dac units committed in time t
            vSTART_DAC[y in DAC_COMMIT, t=1:T] >= 0 # #number of dac units sent for startup in time t
            vSHUT_DAC[y in DAC_COMMIT, t=1:T] >= 0 # #number of dac units sent for shutdown in time t
            vMDOWN_DAC[y in DAC_COMMIT,t=1:T] >= 0 # number of dac units down for maintenance in time t
            vMSHUT_DAC[y in DAC_COMMIT,t in maintenance_begin_hours] >= 0 # dac units for maintenance shutdown selected only from maintenance begin hours
        end)

         # set the unit commitment variables to integer is UC = 1
        for y in DAC_COMMIT
            if setup["UCommit"] == 1
                set_integer.(vCOMMIT_DAC[y,:])
                set_integer.(vSTART_DAC[y,:])
                set_integer.(vSHUT_DAC[y,:])
                set_integer.(vMDOWN_DAC[y,:])
                set_integer.(vMSHUT_DAC[y,:])
            end
        end

        ## units committed are less than max possible number of units
        @constraints(EP, begin
            [y in DAC_COMMIT, t=1:T], vCOMMIT_DAC[y,t] <= (EP[:vCAP_DAC][y]) / dfDac[y,:Cap_Size]   
            [y in DAC_COMMIT, t=1:T], vSTART_DAC[y,t] <= (EP[:vCAP_DAC][y]) / dfDac[y,:Cap_Size]    
            [y in DAC_COMMIT, t=1:T], vSHUT_DAC[y,t] <= (EP[:vCAP_DAC][y]) / dfDac[y,:Cap_Size]     
        end)

        #### this ensures number of units down for maintenance in possible mainteance hours are less than max possible units  
        @constraint(EP, [y in DAC_COMMIT, t in maintenance_begin_hours], vMSHUT_DAC[y, t] .<= (EP[:vCAP_DAC][y]) / dfDac[y,:Cap_Size])

        ### Units down for maintenance cannot be committed
        @constraint(EP, [y in DAC_ID, t in 1:T], vMDOWN_DAC[y,t] + vCOMMIT_DAC[y, t] <= (EP[:vCAP_DAC][y]) / dfDac[y,:Cap_Size])

        # Commitment state constraint linking startup and shutdown decisions (Constraint #1)
        p = hours_per_subperiod
        @constraints(EP, begin
            [y in DAC_ID, t in 1:T], vCOMMIT_DAC[y,t] == vCOMMIT_DAC[y, hoursbefore(p, t, 1)] + vSTART_DAC[y,t] - vSHUT_DAC[y,t]
        end)

        ### specific maintenance constraints

        #this ensures that number of units down in time t for maintenance is the sum of units down in the previous hours of maintenance duration # int, int, vector, steprange
        @constraint(EP, [y in DAC_ID, t in 1:T], vMDOWN_DAC[y,t] .==sum(vMSHUT_DAC[y,t] for t in controlling_maintenance_start_hours(p, t, maintenance_duration, maintenance_begin_hours)) )  

        #### this ensures that across all possible maintenance start hours, the sum of units shut down are the required number needing maintenance per year
        @constraint(EP, [y in DAC_ID], sum(vMSHUT_DAC[y,t] for t in maintenance_begin_hours)>=(EP[:vCAP_DAC][y] / dfDac[y,:Cap_Size]) / maint_freq_years)


        #max and min output constraints
        @constraints(EP, begin
        # Minimum negative CO2 per DAC "y" at hour "t" > Min stable CO2 capture
            [y in DAC_COMMIT, t=1:T], EP[:vCO2_DAC][y,t] >= dfDac[y,:Min_DAC]*EP[:vCOMMIT_DAC][y,t]*(dfDac[y,:Cap_Size]) /(dfDac[y, :CF]*sum(omega[t] for t in 1:T))
        # Maximum negative CO2 per DAC "y"  "y" at hour "t" < Max capacity per hour
            [y in DAC_COMMIT, t=1:T], EP[:vCO2_DAC][y,t] <= (dfDac[y,:Cap_Size]) *EP[:vCOMMIT_DAC][y,t] /(dfDac[y, :CF]*sum(omega[t] for t in 1:T))
         end)

         ## For Start Hours
        # Links last time step with first time step, ensuring position in hour 1 is within eligible ramp of final hour position
        # rampup constraints
        @constraint(EP,cRampUpDAC[y in DAC_COMMIT, t in 1:T],
                EP[:vCO2_DAC][y,t] - EP[:vCO2_DAC][y, hoursbefore(p, t, 1)] <= (dfDac[y,:Ramp_Up_Percentage]*dfDac[y,:Cap_Size]*(EP[:vCOMMIT_DAC][y,t]-EP[:vSTART_DAC][y,t]) 
            + min(max(dfDac[y,:Min_DAC],dfDac[y,:Ramp_Up_Percentage]))*dfDac[y,:Cap_Size]*EP[:vSTART_DAC][y,t] - dfDac[y,:Min_DAC]*dfDac[y,:Cap_Size]*EP[:vSHUT_DAC][y,t]) / (dfDac[y, :CF]*sum(omega[t] for t in 1:T)) )

        # rampdown constraints
        @constraint(EP,cRampDnDAC[y in DAC_COMMIT, t in 1:T],
            EP[:vCO2_DAC][y, hoursbefore(p,t,1)] - EP[:vCO2_DAC][y,t] <= (dfDac[y,:Ramp_Dn_Percentage]*dfDac[y,:Cap_Size]*(EP[:vCOMMIT_DAC][y,t]-EP[:vSTART_DAC][y,t])
        - dfDac[y,:Min_DAC]*dfDac[y,:Cap_Size]*EP[:vSTART_DAC][y,t] + min(max(dfDac[y,:Min_DAC],dfDac[y,:Ramp_Dn_Percentage]))*dfDac[y,:Cap_Size]*EP[:vSHUT_DAC][y,t]) / (dfDac[y, :CF]*sum(omega[t] for t in 1:T)) )



        ### Minimum up and down times once fully switched off or fully switched on
        Up_Time = zeros(Int, nrow(dfDac))
        Up_Time[DAC_COMMIT] .= Int.(floor.(dfDac[DAC_COMMIT,:Up_Time]))
        @constraint(EP, [y in DAC_COMMIT, t in 1:T],
            EP[:vCOMMIT_DAC][y,t] >= sum(EP[:vSTART_DAC][y, hoursbefore(p, t, 0:(Up_Time[y] - 1))])  #minimum number of plants online are those that were online uptime hours before
        )

        Down_Time = zeros(Int, nrow(dfDac))
        Down_Time[DAC_COMMIT] .= Int.(floor.(dfDac[DAC_COMMIT,:Down_Time]))
        @constraint(EP, [y in DAC_COMMIT, t in 1:T],
        EP[:vCAP_DAC][y]/dfDac[y,:Cap_Size] - EP[:vCOMMIT_DAC][y,t] - vMDOWN_DAC[y,t] >= sum(EP[:vSHUT_DAC][y, hoursbefore(p, t, 0:(Down_Time[y] - 1))])  #minimum number of plants offline are those that were offline downtime hours before
        )
    end

    ### Ani's postOpMaint CF constraints
    # annual removals for DAC
    @expression(EP, eDAC_removals_annual[y in DAC_ID], sum(vCO2_DAC[y, t] for t in 1:T))

    # annual DAC removals must be greater than the minimum DAC annual deployment required for each DAC type
    @constraint(EP, cDACRemoval[y in DAC_ID], eDAC_removals_annual[y] .>=  dfDac[y, :Deployment])
    # can't remove more CO2 than the capacity, should force new build capacity
    @constraint(EP, cDACRemovalCapacity[y in DAC_ID], sum(vCO2_DAC[y,t] for t in 1:T) .<= vCAP_DAC[y]) # don't need the CF multiplier because the CF should be 100% except for maintenance downtime

    # preOpMaint constraint capacity constraints
    ## max annual capacity constraint for DAC
    ##  Constraint on maximum annual DAC removal (if applicable) [set input to -1 if no constraint on maximum capacity]
	##  DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MW is >= Max_Cap_MW and lead to infeasabilty
    # @expression(EP, DAC_removals_hourly[y = 1:G_DAC], sum(vCO2_DAC[y, t] * omega[t] for t in 1:T))
    ##  @constraint(EP, cDAC_CF, sum(DAC_removals_hourly[y] for y in 1:G_DAC) ==  sum(vCAP_DAC[y] * dfDac[y, :CF]  for y in 1:G_DAC))
    # @constraint(EP, cDAC_CF[y = 1:G_DAC], DAC_removals_hourly[y] .==  (vCAP_DAC[y] .* dfDac[y, :CF]))
    ##  this cDAC_removal enforces one zone having capacity and one zone having none (if Deployment = 1Mt for zone 1, = 0Mt for zone 2)
    ## @constraint(EP, cDAC_removal[y = 1:G_DAC], DAC_removals_hourly[y] ==  dfDac[y, :Deployment])
    # @constraint(EP, cDAC_removal, sum(DAC_removals_hourly[y] for y in 1:G_DAC) >=  sum(dfDac[y, :Deployment] for y in 1:G_DAC))

    
    # get the CO2 balance 
    # the net negative CO2 for each DAC y at each hour t, CO2 emissions from heat consumption minus CO2 captured by DAC = net negative emissions
    ## remove eDAC_heat_CO2[y,t] since we are using heat pump and RE electricity for heat, not natural gas / fuel 
    ## should we add grid emissions in net CO2 expression? 
    @expression(EP, eCO2_DAC_net[y in DAC_ID, t = 1:T], 0.0 - vCO2_DAC[y,t] )  
    # the net negative CO2 from all the DAC facilities
    @expression(EP, eCO2_DAC_net_ByZoneT[z = 1:Z, t = 1:T], 
        sum(eCO2_DAC_net[y, t] for y in dfDac[(dfDac[!, :Zone].==z), :DAC_ID]))  
    # the net negative CO2 from all DAC facilities during the whole year
    @expression(EP, eCO2_DAC_net_ByZone[z = 1:Z], 
        sum(eCO2_DAC_net_ByZoneT[z, t] * omega[t] for t in 1:T))
    # sum of net CO2 across all the zone
    @expression(EP, eCO2_ToT_DAC_net, sum(eCO2_DAC_net_ByZone[z] for z in 1:Z))


    # separately account for the amount of CO2 that is captured.
    # actually the eCO2_net should be the total sequestration carbon. since eCO2_DAC_net should be a negative value, put minus sign in front of it..
    # costs associated with co2 transport & storage ($/(t CO2/h)) = captured co2 (t CO2/h) * Co2 transport and storage cost ($/t CO2)
    @expression(EP, eCCO2_TS_ByPlant[y in DAC_ID, t = 1:T], vCO2_DAC[y, t]* dfDac[y, :CO2_Transport_Storage_Per_t])
    # the sequestrated CO2 from all DAC facilities 
    @expression(EP, eCCO2_TS_ByZoneT[z = 1:Z, t = 1:T], 
        sum(eCCO2_TS_ByPlant[y, t] for y in dfDac[(dfDac[!, :Zone].==z), :DAC_ID]))    
    # the sequestrated CO2 from all DAC facilities during the whole year ($/t CO2)
    @expression(EP, eCCO2_TS_ByZone[z = 1:Z], 
        sum(eCCO2_TS_ByZoneT[z, t] * omega[t] for t in 1:T))
    # sum of CO2 sequestration costs.
    @expression(EP, eCTotalCO2TS, sum(eCCO2_TS_ByZone[z] for z in 1:Z))


    EP[:eObj] += eCTotalCO2TS

    @expression(EP, eDACcostTOTAL, eCTotalVariableDAC + eTotalCFixedDAC + eTotalCFixedDAC_Energy + eCTotalCO2TS)

    return EP
end