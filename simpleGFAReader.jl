
using StringViews
include("./NodeFlips.jl")

using Profile 
using ProfileView

mutable struct Entry 
    contig_id::Vector{UInt8}
    numbers::Vector{Int32}
    i::Int64
end

struct Reader
    f::IOStream 
    arr::Vector{UInt8}
    entry::Entry
end 

function eachEntry(gfa::String; buffer_size=10_000_000)
    # For the streaming
    h = open(gfa, "r")
    buffer = zeros(UInt8, buffer_size)
    # Storing the entry information
    id_buffer = zeros(UInt8,100)
    number_buffer = zeros(Int32, round(Int64, buffer_size/8)) # ~ bytes to int guess
    return Reader(h, buffer, Entry(id_buffer, number_buffer, 1))
end

function find_newline(bytes::Vector{UInt8})
    newline_loc = findfirst(isequal(UInt8('\n')), bytes)
    return newline_loc
end

function move_to_begin(arr::Vector{UInt8}, from::Int)
    # move bytes to the beginning of the buffer 
    #println("Slice to copy")
    #println(view(arr, from:length(arr)))
    j = 1 
    for i in from:length(arr)
        arr[j] = copy(arr[i])
        j +=1
    end
    # println("Till: ",j )
    for i in j:length(arr)
        arr[i] = 0x00 
    end
    
    return j
end

function add_number!(reader::Reader, number::Int32, number_i::Int64)
    # Store the number in the number array and increase its size by 20% 
    # if it does not fit 
    if number_i > length(reader.entry.numbers)
      #  println("Resizing number vector")
        increased_i = round(Int64, length(reader.entry.numbers) * 1.2) # 20% increase
        resize!(reader.entry.numbers, increased_i)
        reader.entry.numbers[number_i] = number
    else 
        reader.entry.numbers[number_i] = number 
    end
end

function find_line(reader::Reader, state::Int)
    # Move the unprocessed bytes to the fron of the array 
    if state > 0
        # Move the old bytes to the start 
        j = move_to_begin(reader.arr, state)
    else 
        j = 1
    end
    
    deposit_v = view(reader.arr, j:length(reader.arr))
    # Read a chunk of bytes 
    readbytes!(reader.f, deposit_v, length(deposit_v))
    # check if we found a newline already 
    newline_loc = find_newline(reader.arr)
    while isnothing(newline_loc)
        # We didn't find the newline yet so keep resizing the 
        # byte array till we do 
        i = length(reader.arr)
        increased_i = round(Int64, i * 1.2) # 20% increase
        resize!(reader.arr, increased_i)
        # Read more bytes 
        readbytes!(reader.f, view(reader.arr, i+1:increased_i), increased_i - i)
        # Check newline againn 
        newline_loc = find_newline(reader.arr)
    end
    return newline_loc
end

function parse_line(reader::Reader, newline_location::Int64)
    line = view(reader.arr, 1:newline_location-1)
    tab_loc = findfirst(isequal(0x09), line)
    # copy the id to the buffer
    for i in 1:tab_loc-1
        reader.entry.contig_id[i] = line[i]
    end 
    # Now we also have to parse the numbers from the byte vectors
    # and store that in the number array (dont forget last entry)
    prev_loc = tab_loc
    space_loc = findnext(isequal(UInt8(' ')), line, prev_loc + 1)
    number_i = 1
    while !isnothing(space_loc)
        number_bytes = view(line, prev_loc+1:space_loc-2) # Also exclude the +/1 sign 
        number = parse(Int32, StringView(number_bytes))
        number = line[space_loc-1] == 0x2d ? flipnode(number) : noflip(number)
        add_number!(reader, Int32(number), number_i)
        number_i +=1
        prev_loc = space_loc
        space_loc = findnext(isequal(UInt8(' ')), line, prev_loc + 1)
    end
    last_number =  parse(Int32, StringView(view(line, prev_loc+1:length(line)-1)))
    last_number = line[end] == 0x2d ? flipnode(last_number) : noflip(last_number)
    add_number!(reader, Int32(last_number), number_i)
   # println("Last number: ", last_number)
    reader.entry.i = number_i - 1
    return reader.entry
end


@inline function Base.iterate(reader::Reader)
    # Find the first newline char
    newline_loc = find_line(reader, 0)
    entry = parse_line(reader, newline_loc)
    # return StringView(view(reader.arr, 1:newline_loc-1)), newline_loc + 1
    return entry, newline_loc + 1
end

@inline function Base.iterate(reader::Reader, state::Int)
    # Find the first newline char
    if !eof(reader.f)
        newline_loc = find_line(reader, state)
        entry = parse_line(reader, newline_loc)
        entry, newline_loc + 1
       # return StringView(view(reader.arr, 1:newline_loc-1)), newline_loc + 1
    else 
        return nothing
    end
end

using BenchmarkTools

function test2()
    ProfileView.@profview for (i, entry) in enumerate(eachEntry("sub_test.txt"))
        path = view(entry.numbers, 1:entry.i+1)
        #println(StringView(entry.contig_id))
    end
end 

function test() 
   
    @btime begin 
        path_pre_alloc = zeros(Int32, 20_000_000)
        for line in eachline("sub_test.txt")
            identifier, path = split(line, "\t")
            i = 1
            for node in split(path, " ")
                number = parse(Int64, view(node, 1:length(node)-1))
                number = node[end] == "-" ? number * -1 : number 
                path_pre_alloc[i] = number
                i +=1
            end
        end
    end
end

function main() 
   # test()
    test2()
end

main()