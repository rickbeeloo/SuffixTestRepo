using Random 
Random.seed!(1234)
using BenchmarkTools

const LIBSAIS = "/home/codegodz/tools/fais/libsais/libsais.so.2"
const MASK = Int32(1<<30) 

using Profile
using ProfileView
using StatProfilerHTML

struct RefLoc 
    ref_id::Int32 
    loc::Int64 
    reversed::Bool
end

struct Color
    origin::Vector{RefLoc}
    max_len::Vector{Int32}
end

function create_suffix_array(in_vector::Vector{Int32}, free_space::Int32)
    out_vector = zeros(Int32, length(in_vector) + free_space)
    d = Dict(Iterators.map(reverse,pairs(sort(unique(in_vector)))))
    in_vector = Int32.(get.(Ref(d), in_vector, 0))
    n = length(in_vector)
    k = maximum(in_vector) +1 
    @ccall LIBSAIS.libsais_int(in_vector::Ptr{Int32}, out_vector::Ptr{Int32}, n::Int32, k::Int32, free_space::Int32)::Int32
    return out_vector
end

flipnode(n::Int32) = n ⊻ MASK
isflipped(n::Int32) = ifelse(n & MASK != 0, true, false)
noflip(n::Int32) = n & (MASK - 1)
convert_node(n::Int32) = n < 0 ? flipnode(abs(n)) : noflip(n) 

function convert_nodes!(in_vector::Vector{Int32})
    for i in eachindex(in_vector)
        in_vector[i] = convert_node(in_vector[i])
    end
end

function concat_with_seperator(vectors::Vector{Vector{Int32}})
    # Probably a cleaner way to do this :) Like using map 
    # to flip the node ids in slices and copying the slices to the output
    total_size = sum(map(length, vectors))
    concat_arr = zeros(Int32, total_size + length(vectors)) 
    vect_id = -1 * length(vectors) # to have it decending for sorting
    for i in eachindex(vectors)
        convert_nodes!(vectors[i])
    end
    # Concat with seperator + store sign in most significant bit
    i = 1
    @inbounds for v in vectors
        for node_id in v 
            concat_arr[i] = node_id
            i +=1
        end 
        concat_arr[i] = vect_id 
        vect_id +=1
        i +=1
    end
    return concat_arr
end

function create_k_suffix_array(vectors::Vector{Vector{Int32}}, free_space::Int32)
    concat_array = concat_with_seperator(vectors)
    suffix_array = create_suffix_array(concat_array, free_space)
    suffix_array .+= 1
    return concat_array, suffix_array
end


function locate_insert_point(sa::Vector{Int32}, concat_arr::Vector{Int32}, ref::AbstractVector{Int32})
    low = 1 
    high = length(sa)
    suffix_arr_len = length(sa)
    @inbounds while low <= high 
        mid = low + ((high - low) >>> 0x01) 
        suffix_start = sa[mid]
        suffix = view(concat_arr, suffix_start:suffix_arr_len)
        if ref < suffix 
            high = mid - 1
        elseif ref > suffix
            low = mid + 1
        else 
            # then they should be equal, exact match
            return mid 
        end
    end
    return low
end

function matches_till(arr1::AbstractVector{Int32}, arr2::AbstractVector{Int32})
    smallest = min(length(arr1), length(arr2))
    @inbounds for i in 1:smallest
        if arr1[i] != arr2[i]
            return i - 1
        end 
    end
    return smallest
end

function get_suffix_match_location(sa::Vector{Int32}, concat_arr::Vector{Int32}, array_to_find::AbstractVector{Int32}, point_in_sa::Int)
    # Get the suffix start position 
    suffix_start = sa[point_in_sa]
    # Get a view of the whole suffix to iterate over
    suffix = view(concat_arr, suffix_start:length(concat_arr))
    # Check for how long the suffix is the same as the array_to_find
    match_size = matches_till(array_to_find, suffix)
    match_size == 0 && return nothing 
    return suffix_start:suffix_start+match_size-1
end


function scan_flanks(sa::Vector{Int32}, insert_point::Int, array_to_find::AbstractVector{Int32}, concat_arr::AbstractVector{Int32})
    # Now we know where to look in the suffix area, check if we can prefix match the suffix left/right from this point

    concat_matching_indexes = []

    # Scan left of insert point
    left_shift = insert_point
    while left_shift - 1 > 0
        left_match_location = get_suffix_match_location(sa, concat_arr, array_to_find, left_shift - 1)
        isnothing(left_match_location) && break
        push!(concat_matching_indexes, left_match_location)
        left_shift -= 1
        break
    end

    # Scan right of insert point (including the insert point index itself)
    right_shift = insert_point
    while right_shift <= length(sa)
        right_match_location = get_suffix_match_location(sa, concat_arr, array_to_find, right_shift)
        isnothing(right_match_location) && break 
        push!(concat_matching_indexes, right_match_location)
        right_shift += 1
        break
    end

    return concat_matching_indexes


end

function find_longest_match(ref::Vector{Int32}, position::Int64, sa::Vector{Int32}, query_concat::Vector{Int32})
    array_to_find = view(ref, position:length(ref))
    insert_point = locate_insert_point(sa, query_concat, array_to_find)
    # Scan the left and right flank of the suffix array for matches
    concat_matching_indexes = scan_flanks(sa, insert_point, array_to_find, query_concat)
    return concat_matching_indexes
end



function slide_over_ref(ref_id::Int32, ref::Vector{Int32}, sa::Vector{Int32}, query_concat::Vector{Int32}, query_colors::Color, size_map::Dict{Int32,Int64}, reversed::Bool)
    c = 0
    #println(ref_id, " reversed: ", reversed, " with length: ", length(ref))
    for i in 1:length(ref)
        # Find all matches with the Qs for this sub-ref region 
        query_match_locations = find_longest_match(ref, i, sa, query_concat)
        if length(query_match_locations) > 0
            # For each location obtain the nt sizes
            for m in query_match_locations
                q_match_region = view(query_concat, m)
                # Get the nt sizes from the size map 
                q_match_size = sum(get.(Ref(size_map), q_match_region, 0))
                # If the sizes are bigger (based on the colors) replace origin and max size
                for i in m
                    if query_colors.max_len[i] < q_match_size
                        query_colors.max_len[i] = q_match_size
                        query_colors.origin[i] = RefLoc(ref_id, i, reversed)
                    end
                end
            end
        end
    end
    return c
end

function find_longest_matches!(ref_id::Int32, ref::Vector{Int32}, sa::Vector{Int32}, query_concat::Vector{Int32}, query_colors::Color, size_map::Dict{Int32,Int64})
    # Convert the ref nodes to search forward matches
    c = 0
    convert_nodes!(ref)
    c += slide_over_ref(ref_id, ref, sa, query_concat, query_colors, size_map, false)

    # Then flip for reverse matches 
    ref = flipnode.(reverse!(ref))
    c += slide_over_ref(ref_id, ref, sa, query_concat, query_colors, size_map, true)
    
end

function test(queries::Vector{Vector{Int32}}, refs::Vector{Vector{Int32}})
    println("Concat and building suffix array")
    concat_arr, sa = @time create_k_suffix_array(queries, Int32(0))
    println("Q array size; ", length(concat_arr))
    println("\nAllocating colors")
    query_colors = @time Color(Vector{RefLoc}(undef,length(concat_arr)), zeros(Int32, length(concat_arr)))
    println("\nGetting unique nodes")
    unique_nodes = @time Set(reduce(vcat, queries))
    println("\nCreating size map")
    size_map = @time Dict(unique_nodes .=> ones(Int64, length(unique_nodes))) 
    println("Finding longest matches")
    @profview find_longest_matches!(Int32(2), refs[2], sa, concat_arr, query_colors, size_map)
    println("Origins: ",  query_colors.origin[1:20])
    println("Lens: ", query_colors.max_len[1:20])
end



function main() 
    #f = "/media/codegodz/TOSHIBA EXT/staph_cuttlefish_graph.gfa_reduced.cf_seq_head-50000"
    f = "sub_test.txt"
    queries = Vector{Vector{Int32}}()
    refs = Vector{Vector{Int32}}()
    count = 0
    qs_done = false
    q_include = 5
    max_count = 10
    

    last_genome_id = ""
    i = 0
    for line in eachline(f) 
        i +=1
        println(i)

        identifier, path = split(line, "\t")
        tag, genome_id, contig_id = split(identifier, "_")

        # Get the numbers from the path 
        nodes = Int32[]
        for node in split(path, " ")
            if node[end] == '+'
                mult = 1
            else 
                mult = -1
            end 
            node_id = parse(Int32, node[1:length(node)-1]) * mult
            push!(nodes, node_id)
        end

        if genome_id != last_genome_id
            count += 1
            

            if count == q_include !qs_done
                println("Querires done, switching to refs")
                qs_done = true
            end

            last_genome_id = genome_id

        end

        # Store them
        if !qs_done
            println(count, " Pushign q")
            push!(queries, nodes)
        else 
            println(count, " Pushing r")
            push!(refs, nodes)
        end

        count == max_count && break
        
    end
    # Call suffix 
    test(queries, refs)

end

main()


