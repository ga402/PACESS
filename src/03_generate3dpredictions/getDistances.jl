using CSV: makeunique
using DataFrames: make_unique
using Base: _mask1_uint128, final_shred!
using DataFrames
using Distances
using CSV
using ArgParse
using DataFrames
using Queryverse
using Pipe:@pipe
using Statistics
include("getDistancesUtils.jl")
using .getDistancesUtils # loads local module


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--file"
            help = "input csv file"
            arg_type= String
            default="file.csv"
            required=true
        "--output"
            arg_type=String
            help = "output csv file"
        "--pixelsXYZ"
            help = "xyz pixel dimensions, eg. 0.87, 0.87, 5"
            nargs = '+'
            arg_type= Float64
            default= [0.87,0.87, 5]
        "--cells"
            help = "cells to select"
            nargs = '+'
            arg_type= Int64
            default= [0,2, 4]
        "--intensity_match"
            help = "cells to select"
            nargs = '+'
            arg_type= String
            default= ["c2","c4", "c0"]
    end

    return parse_args(s)
end


function getArgs()
    parsed_args = parse_commandline()
    return parsed_args["file"], parsed_args["output"], parsed_args["pixelsXYZ"], parsed_args["cells"] , parsed_args["intensity_match"] 
end


"""
Quick function to remove white space in lists
"""
ff(a)= filter(x -> !isspace(x), a)



"""
Main function
"""
function main()
    #=
    Load the file names..
    =#
    input, output ,xyz, cells, intensity_match = getArgs()

    # example parameters...    
    #input =  "myinputdata.csv"
    #output = "myoutput.csv"
    #xyz = [0.41, 0.41, 5] - the relative distances for each dimension
    #cells = [0, 2, 4] - cells you want to select from the dataset
    #intensity_match = ["c2", "c4", "c0"] - the names of your cells you want 

    # create a dictionary of cells to intensities in dataframe
    intensity_match=[ff(x) for x in intensity_match]
    dict_cells = Dict(zip(cells, intensity_match))


    # Check the parameters in the model 
    println("input = $input")
    println("output = $output")
    println("cells = $cells")
    println("intensity_match = $intensity_match")
    println("dict_cells = $dict_cells")

    # convert the xyz pixels to individual values
    x, y, z= Tuple(x for x in xyz)
    println("pixel dimensions x = $x, y = $y, z = $z")

    # read the dataset
    df = CSV.File("$input") |> DataFrame;
    #df = df |> @select(-:Column1) |> DataFrame; # remove any unwanted column..

    #=
    Run as a loop through each cell type.
    =#
    for i in cells
        ddf = df |> @filter(_.cell == i)  |> DataFrame;
        
        # main function.
        df_out = identityEuclideanGroups(ddf, x,y, z, dict_cells, i) 

        groupdf = df_out |>
        @mutate(wmin = _.x - (_.w/2),
                wmax = _.x + (_.w/2),
                hmin = _.y - (_.h/2),
                hmax = _.y + (_.h/2),
                area = _.h * _.w) |>
        @groupby(_.FinalCell) |>
        @map({cell_number=key(_), 
                zmin=minimum(_.z),
                zmean=mean(_.z),
                zmax=maximum(_.z),
                hmax=quantile(_.hmax, 0.50),
                hmin=quantile(_.hmin, 0.50),
                wmax=quantile(_.wmax, 0.50),
                wmin=quantile(_.wmin, 0.50),
                x = mean(_.x),
                y =mean(_.y)
                }) |>
        DataFrame;

        CSV.write("$(output[1:end-4])$i.csv", groupdf)
    end
end

# run.
main()





