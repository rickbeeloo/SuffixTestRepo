
const LIBSAIS = "/home/codegodz/tools/fais/libsais/libsais.so.2"
const MASK = Int32(1<<30) 


flipnode(n::Int32) = n ⊻ MASK
isflipped(n::Int32) = ifelse(n & MASK != 0, true, false)
noflip(n::Int32) = n & (MASK - 1)
convert_node(n::Int32) = n < 0 ? flipnode(abs(n)) : noflip(n) 

struct Color 
    origin::Vector{Int32}
    len::Vector{Int32}
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
    out_vector .+= 1
    return out_vector
end

function create_k_suffix_array(vectors::Vector{Vector{Int32}}, free_space::Int32)
    concat_array = concat_with_seperator(vectors)
    suffix_array = create_suffix_array(concat_array, free_space)
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
            return Int32(mid)  # then they should be equal, exact match
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

function update_color!(color::Color, ref_id::Int32, match_start::Int32, match_size::Int32)
    match_end = match_start+match_size-1
    for i in match_start:match_end
        if color.len[i] < match_size 
            color.len[i]  = match_size
            color.origin[i] = ref_id
        end
    end
end

function reverse_complement_ref!(ref::Vector{Int32})
    reverse!(ref)
    @inbounds for i in eachindex(ref)
        ref[i] = flipnode(ref[i])
    end
end

function convert_nodes!(in_vector::Vector{Int32})
    for i in eachindex(in_vector)
        in_vector[i] = convert_node(in_vector[i])
    end
end

# Check left and right from insert point to see if we have a match
function decide_on_seed(insert_point::Int32, ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32)
    # Check left for a match
    left_of_ip = insert_point > 1 ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point-Int32(1),  Int32(0)) : 0 
    left_of_ip > 0 && return insert_point-Int32(1), Int32(left_of_ip)

    # Check right for a match, no need to check if it's outside the bounds of the SA <= length(sa)
    right_of_ip = insert_point <= length(sa) ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point, Int32(0)) : 0
    right_of_ip > 0 && return insert_point, Int32(right_of_ip)

    # Neither actually matches our Ref, return 0,0 to move on to the next one
    return Int32(0), Int32(0)
end

function matches_till(ref::AbstractVector{Int32}, ref_start::Int32, ca::Vector{Int32}, q_start::Int32)
    (ref_start > length(ref) || q_start > length(ca)) && return 0
    smallest_n = min(length(ref)-ref_start+1, length(ca)-q_start+1)
    for i in 1:smallest_n
        if ref[ref_start + i - 1] != ca[q_start+i-1]
            return Int32(i - 1)
        end 
    end 
    return Int32(smallest_n)
end

function check_this_point(ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32, point::Int32, skip::Int32)
    # Given a point in the suffix array, compare the suffix to the Ref 
    ca_suffix_index = sa[point]
    ca_start = ca_suffix_index + skip
    ref_start = ref_start + skip
    match_size = matches_till(ref, ref_start, ca, ca_start) + skip
    return match_size
end

function extend_from_point!(ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, lcp::Vector{Int32}, point::Int32, forward::Bool, ref_start::Int32, match_size::Int32, ref_id::Int32, color::Color)
    move_dir = forward ? 1 : -1
    lcp_dir  = forward ? 0 :  1
    
    i = point += move_dir
    while i >= 1 && i < length(sa) && lcp[i + lcp_dir] > 0
        # We can skip the LCP part when extending, note though we also have to 
        # check the previous match size so min(lcp valu, prev match size)
        start_check_from = Int32(min(lcp[i + lcp_dir], match_size))
        # Check the size of this match starting from +1 of the LCP value)
        match_size = check_this_point(ca, sa, ref, ref_start, Int32(i), start_check_from )
        update_color!(color, ref_id, sa[i], Int32(match_size))
        i += move_dir        
    end
  
end

function align(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    max_match_index = Int32(0)
    max_match_size = Int32(0)
    ref_start = Int32(1)
    while ref_start <= length(ref)

        # Do binary search to locate the insert point
        insert_point = locate_insert_point(sa, ca, view(ref, ref_start:length(ref)))
        max_match_index, max_match_size = decide_on_seed(insert_point, ca, sa, ref, ref_start)

        # If we have a match keep using the linked list to extend 
        if max_match_size > 0 
            while ref_start <= length(ref)
                # Check the match size at this point 
                max_match_size = check_this_point(ca, sa, ref, ref_start, max_match_index, Int32(max_match_size-1)) # skip k-1
                
                # If we don't have any match we don't have to check the flanks
                max_match_size == 0 && break 
                update_color!(color, ref_id, sa[max_match_index], Int32(max_match_size))
                              
                # Check up and down in suffix array for other matches
                extend_from_point!(ca, sa, ref, lcp, max_match_index, false, ref_start, Int32(max_match_size), ref_id, color)
                extend_from_point!(ca, sa, ref, lcp, max_match_index, true, ref_start, Int32(max_match_size), ref_id, color)

                # # Move to next location in suffix array for the second around
                max_match_index = inv_perm_sa[sa[max_match_index]+1]
                ref_start += Int32(1)
            end 
        else 
            # No match at current point, move +1 to do a binary search again
            ref_start += Int32(1)
        end
    end
end

function align_forward_and_reverse(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    # First do the forward align 
    convert_nodes!(ref)
    align(ref_id, color, ca, sa, ref, inv_perm_sa,lcp)
    # Flip the nodes and reverse to do the reverse alignment 
    reverse_complement_ref!(ref)
    align(ref_id, color, ca, sa, ref, inv_perm_sa,lcp)
end

function test(qs, refs)
    println("Building data structures")
    println("Qs: ", length(qs))
    println("Refs: ", length(refs))
    ca, sa = @time create_k_suffix_array(qs, Int32(0))
    inv_sa_perm = @time inverse_perm_sa(sa)
    lcp = @time build_lcp(sa, ca)

     len = zeros(Int32, length(ca))
     ori = zeros(Int32, length(ca))
     color = Color(len, ori)

    println("Started aligning")
    for (ref_id, ref) in enumerate(refs)
        align_forward_and_reverse(Int32(ref_id), color, ca, sa, ref, inv_sa_perm, lcp)
    end 

    # Just to check if it did something :)
    println(count(x->x>0, color.len))
    println(count(x->x>0, color.origin))
end


function main() 
    f = "sub_test.txt"
    queries = Vector{Vector{Int32}}()
    refs = Vector{Vector{Int32}}()
    count = 0
    qs_done = false
    q_include = 5
    max_count = 11
    last_genome_id = ""
    for line in eachline(f) 
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
            if count == q_include
                qs_done = true
            end
            last_genome_id = genome_id
        end

        # Store them
        if !qs_done
            push!(queries, nodes)
        else 
            push!(refs, nodes)
        end
        count == max_count && break
    end
    # Call suffix 
    println("Start aln")
    test(queries, refs)
end

main()

