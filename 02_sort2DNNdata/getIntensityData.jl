using CSV: makeunique
using DataFrames: make_unique
using Base: _mask1_uint128, final_shred!
using DataFrames
using Distances
using CSV
using ArgParse
using DataFrames
using NPZ
using Queryverse
using Pipe:@pipe
using Statistics


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--data"
            help = "input csv file"
            arg_type= String
            default="file.csv"
            required=true
        "--image"
            help = "image"
            arg_type=String
            default="image"
            required=true
        "--output"
            arg_type=String
            help = "output csv file"
            default="out.csv"
        "--channel_order"
            help = "order of channels for cells"
            nargs = '+'
            arg_type= Int
            default= [1,2, 3, 4]
        "--color_order"
            help = "order of channels for cells"
            nargs = '+'
            arg_type= String
            default= ["dapi","gfp", "bone", "dTomato"]
        "--select_colors"
            help = "selections"
            nargs = '+'
            arg_type= String
            default= ["dTomato","gfp", "all"]  
    end

    return parse_args(s)
end



"""
Quick function to remove white space in lists
"""
ff(a)= filter(x -> !isspace(x), a)



function getArgs()
    parsed_args = parse_commandline()
    return parsed_args["data"], parsed_args["image"], parsed_args["output"], parsed_args["channel_order"], parsed_args["color_order"], parsed_args["select_colors"] 
end
 




"""
Identify the colour from the options
"""
function select_channels(col, dict_colors)
    if col=="all"
        return 0
    end
    try
        v = dict_colors[col]
        return v
    catch e
        println("color not present within options")
    end
end

#select_channels("dapi", dict_colors)

function rgb2gray3D(rgb)
    gray = mean(rgb, dims=4)
    return gray[:, :, :, 1]
end



function image_select(img, col, dict_colors)
    c = select_channels(col, dict_colors)
    if c==0
        return rgb2gray3D(img)
    end
    return img[:, :, :, c]
end


function readNPZImage(input)
    img = npzread(input)
    img = Int.(img)
    return img
end


"""
gets Z level beginning
"""
function getZlevel(path_name)
    x = match(r"Z\d+", path_name).match
    x = split(x, "Z")[2]
    return parse(Int64, x)
end


function getXlevel(path_name)
    x = match(r"X\d+", path_name).match
    x = split(x, "X")[2]
    return parse(Int64, x)
end


function getYlevel(path_name)
    x = match(r"Y\d+", path_name).match
    x = split(x, "Y")[2]
    return parse(Int64, x)
end


"""
Adds the index position for each image in the dataset value
"""
function addBegins(df)
    @assert "image_name" in names(df) "df must contain 'image_name'"
    df[!, :X_begin] = getXlevel.(df[!, :image_name])
    df[!, :Y_begin] = getYlevel.(df[!, :image_name])
    df[!, :Z_begin] = getZlevel.(df[!, :image_name])
    return (df)
end



"""
Calculates the MFI for each location (box) in the image.
"""
function calcMFI(img,df)

    @assert (img |> size |> length) == 3 "image should be only 3 dimensions"
    
    xmin = df.x           #+ #df.X_begin 
    xmax = df.x + df.w    #+ #df.X_begin + df.w 
    ymin = df.y           #+ #df.Y_begin 
    ymax = df.y + df.h    #+ #df.Y_begin + df.h 
    zlevel = df.z
    v = zeros(Float64,size(df)[1])
    for i ∈ 1:size(df)[1]
        z = zlevel[i] + 1 ## 1 index not 0
        yx = ymin[i] + 1
        yy = ymax[i] 
        xx = xmin[i] + 1
        xy = xmax[i]
        v[i]=mean(img[z, yx:yy, xx:xy])
    end
    return(v)
end




function main()
    println("starting...")
    data_path, image_path, output, channel_order, color_order, select_colors = getArgs()
    color_order=[ff(x) for x in color_order]
    select_colors=[ff(x) for x in select_colors]
    println("data_path = $data_path")
    println("image_path = $image_path")
    println("channel order = $channel_order")
    println("color_order = $color_order")
    println("selected_channels = $select_colors")

    # read the dataset
    df = CSV.File(data_path) |> DataFrame;
    df = df |> @select(-:Column1) |> DataFrame; # remove any unwanted column..
    df = addBegins(df);

    # read the image
    img = readNPZImage(image_path);
    println("image size = $(size(img))")

    # dictionary of colours and channels
    dict_colors = Dict(zip(color_order, channel_order))


    for col in select_colors
        img_sub = image_select(img, col, dict_colors)
        if col=="all"
            v = 0
        else
            v = dict_colors[col]
        end
        df[!, ("c"*"$v" * "Intensity")]  = calcMFI(img_sub,df)
    end

    CSV.write("$output",df)

    
end
    
# run.
main()




