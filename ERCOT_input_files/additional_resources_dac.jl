using JuMP
using DataFrames
using CSV

function compare_generation()
    
    path = "/Users/margotadam/Documents/GitHub/DAC/GenX/ERCOT_input_files/Results_"
    # basecase_results_num = "65"
    # policy_results_num = "63"
    # path = "/Users/margotadam/Downloads/Results_"
    basecase_results_num = "79"
    policy_results_num = "80"
    basecase_path = path * basecase_results_num * "/capacity.csv"
    policy_path = path * policy_results_num * "/capacity.csv"
    
    generator_path = "/Users/margotadam/Documents/GitHub/DAC/GenX/ERCOT_input_files/Generators_data.csv"
    # generator_path = "/Users/margotadam/Downloads/DAC research project/DAC data files/conus-dac-v2/base/Inputs/Inputs_p1/Generators_data.csv"

    # read in data
    basecase_df = DataFrame(CSV.File(basecase_path))
    policy_df = DataFrame(CSV.File(policy_path))
    gen_in = DataFrame(CSV.File(generator_path))

    # find Resources that are clean
    # clean generator R_IDs
    vre_gens = gen_in[gen_in.VRE.>=1,:R_ID]
    # geotherm_gens = gen_in[(gen_in.Geothermal.==1),:R_ID]
    nuclear_gens = gen_in[(gen_in.Nuclear.==1),:R_ID]
    # nuclear_gens = gen_in[(gen_in.technology.=="Nuclear"),:R_ID] # needed for CONUS data
    clean_RIDs = [vre_gens; nuclear_gens]

    # filter for only new-build eligible generators in the clean subset
    new_build_gens_input = intersect(gen_in[gen_in.New_Build.==1,:R_ID], gen_in[gen_in.Max_Cap_MW.!=0,:R_ID])
    
    # combine new build gens with clean gens for all qualifying R_IDs
    new_clean_RIDs = intersect(new_build_gens_input, clean_RIDs) 

    # find Resources with associated clean_RIDs
    # clean_resources = gen_in[gen_in.R_ID=clean_RIDs, :Resource]

    clean_resources = []
    for r in new_clean_RIDs
        push!(clean_resources, gen_in[gen_in.R_ID.==r, :Resource][1])
    end

    clean_gens_df = DataFrame(Resource = clean_resources, R_ID = new_clean_RIDs)

    # find new build resources, can do this straight to output data because it's an output field
    # BUG: this doesn't filter for generators eligible for new build, it only selects generators that had added capacity!
    # new_gen_bc = basecase_df[basecase_df.NewCap.>0, :]
    # new_gen_p = policy_df[policy_df.NewCap.>0, :]

    # go through generators output data, find gens to add to DAC resources (logic)
    dac_resources = DataFrame(Resource = [], R_ID = [], built_in_DAC = [], not_built = [], retired_in_bc = [])
    for r in basecase_df.Resource
        if r in clean_resources
            # logic for DAC resource or not 
            # built_in_policy = policy_df[policy_df.Resource.==r, :EndCap] - basecase_df[basecase_df.Resource.==r, :EndCap] >= 5.0
            built_in_policy = policy_df[policy_df.Resource.==r, :EndCap] > basecase_df[basecase_df.Resource.==r, :EndCap]
            not_build = (policy_df[policy_df.Resource.==r, :NewCap][1] < 5.0) && (basecase_df[basecase_df.Resource.==r, :NewCap][1] < 5.0)
            retired_in_bc = basecase_df[basecase_df.Resource.==r, :RetCap] > policy_df[policy_df.Resource.==r, :RetCap]

            if built_in_policy || not_build || retired_in_bc
                # add to DAC resources
                push!(dac_resources.Resource, r)
                push!(dac_resources.R_ID, clean_gens_df[clean_gens_df.Resource.==r, :R_ID][1])
                push!(dac_resources.built_in_DAC,built_in_policy)
                push!(dac_resources.not_built,not_build)
                push!(dac_resources.retired_in_bc,retired_in_bc)
                
            end
        end

    end

    # write output files + results

    # results_df = DataFrame(Resource = basecase_df.Resource, 
    #                         BaseCase_NewCap = basecase_df.NewCap,
    #                         Policy_NewCap = policy_df.NewCap) 
    # # results_df.NewCap_Diff = diff_cap
    # results_df.DAC_Resource = added_col

    CSV.write("dac_additional_resources_" * basecase_results_num * "_" * policy_results_num * ".csv", dac_resources)

    return dac_resources, basecase_df, policy_df
end

compare_generation()