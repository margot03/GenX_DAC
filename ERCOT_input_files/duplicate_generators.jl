using JuMP
using DataFrames
using CSV

function duplicate_gens()
    println("duplicate_gens")

    # load generators dataframe 
    generator_path = "/Users/margotadam/Documents/GitHub/DAC/GenX_DAC/ERCOT_input_files/Generators_data.csv"
    gen_in = DataFrame(CSV.File(generator_path))
    
    # find all clean generators 
    vre_gens = gen_in[gen_in.VRE.>=1,:R_ID]
    # geotherm_gens = gen_in[(gen_in.Geothermal.==1),:R_ID]
    nuclear_gens = gen_in[(gen_in.Nuclear.==1),:R_ID]
    # nuclear_gens = gen_in[(gen_in.technology.=="Nuclear"),:R_ID] # needed for CONUS data
    clean_RIDs = sort([vre_gens; nuclear_gens])

    # println(names(gen_in))
    clean_resources = DataFrame([c => [] for c in names(gen_in)])
    for r in clean_RIDs
        append!(clean_resources, gen_in[gen_in.R_ID.==r, :])
    end
    println(clean_resources.R_ID)
    # increase R_ID numbers so there aren't any repeats
    println("new R_IDs")
    max_R_ID = maximum(gen_in.R_ID)
    for i in 1:length(clean_resources.R_ID)
        # println(max_R_ID + i)
        clean_resources.R_ID[i] = max_R_ID + i
        clean_resources.Resource[i] = clean_resources.Resource[i] + "_DAC"
    end

    # create DAC column (0 = non-duplicate, R_ID from clean_RIDs (non-zero) = duplicate)
    non_duplicate_gens = zeros(length(gen_in.R_ID))
    
    res_df = vcat(gen_in, clean_resources)

    res_df[!, :DAC_R_ID] = vcat(non_duplicate_gens, clean_RIDs)

    dfGen_DAC = res_df[res_df.DAC_R_ID .!= 0,:]
    dfGen_GRID = res_df[res_df.DAC_R_ID .== 0,:]
    initial_RIDs_for_duplicated_gens = intersect(dfGen_DAC.DAC_R_ID, dfGen_GRID.R_ID)
    
    println(dfGen_DAC.DAC_R_ID)
    # println(initial_RIDs_for_duplicated_gens)
    
    # edit duplicated gens for only new_builds
    # new_build_gens = intersect(clean_resources[clean_resources.New_Build.==1,:R_ID], clean_resources[clean_resources.Max_Cap_MW.!=0,:R_ID])
    # dfGen_DAC = DataFrame([c => [] for c in names(clean_resources)])
    # for r in new_build_gens
    #     append!(dfGen_DAC, clean_resources[clean_resources.R_ID.==r, :])
    # end
    # println(dfGen_DAC.R_ID)

    # write out new CSV file of generators 
    CSV.write("/Users/margotadam/Documents/GitHub/DAC/GenX_DAC/ERCOT_input_files/Duplicated_Generators_data.csv", res_df)

end

duplicate_gens()