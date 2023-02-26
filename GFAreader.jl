


using StringViews

struct BufferedReader
    io::IOStream
    buffer::Int64
    tot_alloc::Int64
    path_alloc::Vector{Int32}
    arr::Vector{UInt8}
end

struct GFAEntry 
    contig_id::String 
    path::Vector{Int32}
end


# Function to flip elements in an array to a specified offset(buffer size here)
function flip!(arr::Vector{UInt8}, buffer::Int64)
    @inbounds @simd for i in 1:buffer
        arr[i] = arr[i+buffer]
    end
end


function read_next_chunk!(reader::BufferedReader)
    # Move last read chunk to front of the array
    # (except in first iter)
    flip!(reader.arr, reader.buffer)
    
    # Store new chunk in second part of the array
    bytes_read::Int = readbytes!(reader.io, view(reader.arr, reader.buffer+1:reader.tot_alloc), reader.buffer)  

    # If we read less than the buffer size we have to reset the array
    # values after "bytes_read" as this is old data (previous read)
    if bytes_read < reader.buffer
        @inbounds for i in reader.buffer+bytes_read+1:reader.tot_alloc
            reader.arr[i] = 0x00
        end   
    end   
end

function find_newline(reader::BufferedReader, state::Int64)
    cur_stop = copy(state) + 1
    tab_location = findnext(isequal(UInt8('\n')), reader.arr, state+1)
    if !isnothing(tab_location)
        return cur_stop:tab_location, tab_location
    else 
        return 0:0, cur_stop
    end
end

function eachlineV(file_path::String; buffer_size::Int64=10_000)
    # Allocate buffer array
    tot_alloc = buffer_size * 2
    buffer_arr = zeros(UInt8, tot_alloc)     
    path_alloc = Vector{Int32}(undef, 10_000_000)
    
    # Open the file handle for streaming
    io =  open(file_path, "r")  
    
    # We will set up a buffered reader through which we 
    # stream the file bytes, >4x as fast as a regular reader
    reader = BufferedReader(io, buffer_size, buffer_size*2, buffer_arr, path_alloc)

    # Also populate the reader with the first chunk already 
    read_next_chunk!(reader)
    return reader
end


@inline function Base.iterate(reader::BufferedReader)
    # This is the first iter so only the last half of the array is filled now 
    # hence start reading from buffer + 1
    r, state = find_newline(reader, reader.buffer)
    if r.stop == 0
        guess_buffer_size(reader, r)
    end
    return StringView(view(reader.arr, r)), state
end

@inline function Base.iterate(reader::BufferedReader, state::Int64)
    r, state = find_newline(reader, state)
    if r.start == 0
        if !eof(reader.io)
            read_next_chunk!(reader)
            r, state = find_newline(reader, state - reader.buffer)
        else
            close(reader.io)
            return nothing  
        end 
    end
    # I twould be odd to not reach EOF but still not find 
    # a full line, throw warning
    r.stop == 0 && @warn ("Buffer probably too small")
    return StringView(view(reader.arr, r)), state
end

for line in eachlineV("sub_test.txt", buffer_size=100)
    println(line[1:100])
    break
end