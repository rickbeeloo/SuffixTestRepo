
const LIBSAIS = "/home/codegodz/tools/fais/libsais/libsais.so.2"
const MASK = Int32(1<<30) 

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

function create_suffix_array(in_vector::Vector{Int32}, free_space::Int32)
    out_vector = zeros(Int32, length(in_vector) + free_space)
    d = Dict(Iterators.map(reverse,pairs(sort(unique(in_vector)))))
    in_vector = Int32.(get.(Ref(d), in_vector, 0))
    n = length(in_vector)
    k = maximum(in_vector) +1 
    @ccall LIBSAIS.libsais_int(in_vector::Ptr{Int32}, out_vector::Ptr{Int32}, n::Int32, k::Int32, free_space::Int32)::Int32
    return out_vector
end

function create_k_suffix_array(vectors::Vector{Vector{Int32}}, free_space::Int32)
    concat_array = concat_with_seperator(vectors)
    suffix_array = create_suffix_array(concat_array, free_space)
    suffix_array .+= 1
    return concat_array, suffix_array
end

function build_lcp(sa::Vector{Int32}, V::Vector{Int32}, base::Integer=1)
    T = eltype(sa)
    pos = sa .+ T(1-base)
    n = length(pos)
    lcparr = similar(pos)
    rank = invperm(pos)
    h = 0
    for i in 1:n
        if rank[i] == 1
            continue
        end
        j = pos[rank[i]-1]
        maxh = n - max(i, j)
        while h <= maxh && V[i+h] == V[j+h]
            h += 1
        end
        lcparr[rank[i]] = h
        h = max(h-1, 0)
    end
    lcparr[1] = 0
    return lcparr
end

function locate_insert_point(sa::Vector{Int32}, concat_arr::Vector{Int32}, ref::AbstractVector{Int32})
    low = 1 
    high = length(sa)
    suffix_arr_len = length(sa)
    while low <= high 
        mid = low + ((high - low) >>> 0x01) 
        suffix_start = sa[mid]
        suffix = view(concat_arr, suffix_start:suffix_arr_len)
        if ref < suffix 
            high = mid - 1
        elseif ref > suffix
            low = mid + 1
        else 
            return mid  # then they should be equal, exact match
        end
    end
    return Int32(low)
end

function inverse_perm_sa(sa::Vector{Int32})
    inv_sa_perm = similar(sa)
    for i in 1:length(sa)
        inv_sa_perm[sa[i]] = i
    end
    return inv_sa_perm
end

function matches_till(ref::AbstractVector{Int32}, ref_start::Int32, ca::Vector{Int32}, q_start::Int32)
    (ref_start > length(ref) || q_start > length(ca)) && return 0
    smallest_n = min(length(ref)-ref_start+1, length(ca)-q_start+1)
    #println("Smalles tn: ", smallest_n)
    # println("Ref: ", view(ref,ref_start:length(ref)))
    # println("Q: ", view(ca,q_start:length(ref)))
    for i::Int32 in 1:smallest_n
        if ref[ref_start + i - 1] != ca[q_start+i-1]
            return Int32(i - 1)
        end 
    end 
   # println("full return")
    return smallest_n
end


function extend_from_point(ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, lcp::Vector{Int32}, point::Int32, forward::Bool)
    # When we have a seeding point we have to extend in both directions
    move_dir = forward ? 1 : -1
    lcp_dir  = forward ? 0 :  1
    matches = 0

    i = point += move_dir
    while i >= 1 && i < length(sa) && lcp[i + lcp_dir] > 0
        # We can skip the LCP part when extending
        start_check_from = lcp[i + lcp_dir]
        # Check the size of this match starting from +1 of the LCP value
        #println("Starting check form: ", start_check_from)
        match_size =  matches_till(ref, start_check_from + Int32(1), ca, sa[i] + start_check_from) 
        match_size += start_check_from
        println("(M): Match at: ", i, " of size: ", match_size)
        matches +=1
        i += move_dir
    end
    return matches
end

function check_this_point(ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32, point::Int32, skip::Int32)
    # Given a point in the suffix array, compare the suffix to the Ref 
    ca_suffix_index = sa[point]
    ca_start = ca_suffix_index + skip
    ref_start = ref_start + skip
    match_size = matches_till(ref, ref_start, ca, ca_start) 
    match_size += skip 
    match_size > 0 && println("(M): Match at: ", point, " of size: ", match_size)
    return match_size
end

function decide_on_seed(insert_point::Int32, ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32)
    # Check left for a match
    left_of_ip = insert_point > 1 ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point-Int32(1),  Int32(0)) : 0 
    left_of_ip > 0 && return insert_point-Int32(1), left_of_ip

    # Check right for a match 
    right_of_ip = check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point, Int32(0))
    right_of_ip > 0 && return insert_point, right_of_ip

    # Neither actually matches our Ref
    return 0, 0
end

function align(ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    ref_start = Int32(1)
    while ref_start <= length(ref)

        # Do binary search to locate the insert point
        println("Binary searching for: ", view(ref, ref_start:length(ref)))
        println("Insert lookup: ", view(ref, ref_start:length(ref)) )
        insert_point = locate_insert_point(sa, ca, view(ref, ref_start:length(ref)))
        max_match_index, max_match_size = decide_on_seed(insert_point, ca, sa, ref, ref_start)
        #println("Found max match at: ", max_match_index, " with size: ", max_match_size)
        println()
        # If we have a match keep using the linked list to extend 
        if max_match_size > 0 
            matches = true 
            while matches && ref_start <= length(ref)
                println("Working on: ", view(ref, ref_start:length(ref)))
                # Check the match size at this point
                max_match_size = check_this_point(ca, sa, ref, ref_start, max_match_index, Int32(max_match_size-1)) # skip k-1

                # If we don't have any match we don't have to check the flanks
                max_match_size == 0 && break 
                
                # Check up and down in suffix array for other matches
                up_matches = extend_from_point(ca, sa, ref, lcp, max_match_index, false)
                down_matches = extend_from_point(ca, sa, ref, lcp, max_match_index, true)

                # If we still have matches keep going
                # could also break here but this will invalidate with the while anyway
                matches = up_matches > 0 || down_matches > 0

                # Move to next location in suffix array for the second around
                max_match_index = inv_perm_sa[sa[max_match_index]+1]
                ref_start += Int32(1)
            end 
        else 
            ref_start += Int32(1)
        end
    end
end

function test()
    qs = [Int32[1,2], Int32[1,2,3, 4,1], Int32[1,2, 3, 100, 1, 2]]

    ref = Int32[100, 1]
    # #ref = Int32[1,2,3]
    ca, sa = create_k_suffix_array(qs, Int32(0))
    inv_sa_perm = inverse_perm_sa(sa)
    lcp = build_lcp(sa, ca)


    #println("R: ", ref)
    #println("Q: ", ca)
    println()
    
    println("Suffix array")
    for (i, si) in enumerate(sa)
        println(i, " -> ", view(ca, si:length(sa)))
    end 

    # println("\nLCP array")
    # for (i, l) in enumerate(lcp)
    #     println(i, " -> ", l)
    # end
  
    # println("\nInv permuted")
    # for i in eachindex(ca) 
    #     println("at: ", inv_sa_perm[i], " in the sa: ", view(ca, i:length(ca)))
    # # end
    # @time for (i, ref) in enumerate(refs)
    align(ca, sa, ref, inv_sa_perm, lcp)
    # end
   

end

test()

# function main() 
#     f = "/media/codegodz/TOSHIBA EXT/staph_cuttlefish_graph.gfa_reduced.cf_seq_head-50000"
#     #f = "/home/codegodz/SuffixTestRepo/sub_test.txt"
#     queries = Vector{Vector{Int32}}()
#     refs = Vector{Vector{Int32}}()
#     count = 0
#     qs_done = false
#     q_include = 100
#     max_count = 110
#     last_genome_id = ""
#     for line in eachline(f) 
#         identifier, path = split(line, "\t")
#         tag, genome_id, contig_id = split(identifier, "_")
#         # Get the numbers from the path 
#         nodes = Int32[]
#         for node in split(path, " ")
#             if node[end] == '+'
#                 mult = 1
#             else 
#                 mult = -1
#             end 
#             node_id = parse(Int32, node[1:length(node)-1]) * mult
#             push!(nodes, node_id)
#         end

#         if genome_id != last_genome_id
#             count += 1
#             if count == q_include
#                 qs_done = true
#             end
#             last_genome_id = genome_id
#         end

#         # Store them
#         if !qs_done
#             println("q: ", count)
#             push!(queries, nodes)
#         else 
#             println("r: ", count)
#             push!(refs, nodes)
#         end
#         println(count)
#         count == max_count && break
#     end
#     # Call suffix 
#     println("Start aln")
#     test(queries, refs)
# end

# main()

