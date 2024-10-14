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
	load_dac_params!(setup::Dict, path::AbstractString, inputs::Dict)

Read input parameters related to dac
"""

function load_dac_params!(setup::Dict, path::AbstractString, inputs::Dict)
    filename ="dac_params.csv"

    dac_params = load_dataframe(joinpath(path, filename))

	# Store DataFrame of dac input data for use in model
	inputs["dfDac_params"] = dac_params

    # heat extraction per percent steam diversion is in MW thermal

    println(filename * " Successfully Read!")

   # return inputs
end




