module getDistancesUtilsNEW

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
using NearestNeighbors


export sortMaxArea, createMatrix, remove_diagonal, getMinDistancesFromMatrix, 
getDistancesBelowMaxDiam, determineFinalCellPos, identityEuclideanGroups, getZlevel, getXlevel, getYlevel, addBegins, sortMaxArea_Intensity, calcMFI


"""
gets Z level
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
    sortMaxArea(df_m)
calculates maxArea and sorts by it.
"""
function sortMaxArea(df_m)
    df_m[!, :maxArea]=  df_m.w .* df_m.h
    
    return sort!(df_m, :maxArea, rev= true)
end



"""
    sortMaxArea_Intensity(df_m)
calculates maxArea and sorts by area and intensity!
cell = 0,2, 4 etc..
"""
function sortMaxArea_Intensity(df, dict_cells, cell)
    #@assert "mean_intensity" in names(df) "df must contain 'mean_intensity'"
    df[!, :maxArea]=  df.w .* df.h
    return sort!(df, [(dict_cells[cell]*"Intensity"), "maxArea"], rev= true)
    #return sort!(df, [:mean_intensity, :maxArea], rev= true)
end



"""
Calculates the MFI for each location in the image...
"""
function calcMFI(img,df)

    @assert (img |> size |> length) == 3 "image should be only 3 dimensions"
    
    xmax = df.x + df.X_begin + df.w 
    xmin = df.x + df.X_begin 
    xmax = df.x + df.X_begin + df.w 
    ymin = df.y + df.Y_begin 
    ymax = df.y + df.Y_begin + df.h 
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



"""
create a Matrix from x, y, z cordinates
"""
function createMatrix(df_m, x_correction, y_correction, z_correction)
    mat = df_m  |> @select(:x, :y, :z) |> 
                  @mutate(x = _.x *x_correction, 
                  y = _.y*y_correction, 
                  z= _.z*z_correction) |> 
                  DataFrame |> Matrix |> transpose;
    return float.(mat)
end



"""
Remove diagonal from matrix
 - NO LONGER USED
"""
function remove_diagonal(x)
    mat = Array{Float64}(undef, size(x, 1), size(x, 1) - 1)
    for i = 1:size(mat, 1)
        for j ∈ 1:size(mat, 2)
            if i > j
                mat[i, j] = x[i, j]
            else
                mat[i, j] = x[i, j+1]
            end
        end    
    end
    return mat
end



"""
Extract the minimum distances from the distance matrix
 - NO LONGER USED
"""
function getMinDistancesFromMatrix(R_out)
    R_min = Float64[]
    for i in 1:size(R_out)[1]
        append!(R_min, minimum(R_out[i, :]))
    end
    return R_min
end



"""
Select only those distances below a certain value
"""
function getDistancesBelowMaxDiam(ddf, x_correction, y_correction, z_correction)

    mat = createMatrix(ddf, x_correction, y_correction, z_correction)

    # create a KDTree with 10 leaves
    balltree = BallTree(mat)
    
    for j in 1:size(ddf)[1]
        # get the max diameter
        r = maximum([((ddf[j, :w]*x_correction)) ((ddf[j, :h]*y_correction))])
        if r > 250
            continue
        end
        idxs = inrange(balltree, mat[:, j], r, false)
        # find position already taken...
        taken_pos  = Set(findall(x -> x > 0, ddf[!, :UpdateValue]))
        candidate_pos = Set(idxs)
        ## returns the elements in first set but not sec
        idx = Int.(setdiff(candidate_pos, taken_pos))
        if isempty(idx)
            continue
        end
        ddf[idx, :UpdateValue] .= ddf[j, :CellValueBase]
    end
    return ddf
end


"""
Identify the end position of the value
"""
function determineFinalCellPos(df_out)
    fv = Int64[]
    for i in 1:size(df_out)[1]
        t = df_out[i, :UpdateValue] ==0 ? df_out[i, :CellValueBase] : df_out[i, :UpdateValue]
        append!(fv, t)
    end
    df_out[!, :FinalCell]  = fv
    return df_out
end


"""
Main function for identifying groups based on the euclidean distances between them. Uses the above functions.
"""
function identityEuclideanGroups(df, x_correction, y_correction, z_correction, dict_cells, cell)   

    df =sortMaxArea_Intensity(df, dict_cells,cell)
    df[!, :CellValueBase] = 1:size(df)[1]
    
    try
        df[!, :UpdateValue] .=0
    catch
        insertcols!(df ,8, :UpdateValue => 0, makeunique=true)
    end

    #mat = createMatrix(df, x_correction, y_correction, z_correction)
    
    #R = @pipe pairwise(Euclidean(), mat, dims=1) |> remove_diagonal(_);

    #R = remove_diagonal(R)

    df = getDistancesBelowMaxDiam(df, x_correction, y_correction,z_correction)

    df = determineFinalCellPos(df)

    return df
end


end

