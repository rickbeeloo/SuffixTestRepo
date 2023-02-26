
incude("./NodeFlips.jl")

struct Entry
    contig_id::String 
    path::Vector{Int32}
end

struct GFAreader 
    h::IOStream 
    entry::Entry 
end

@inline function Base.iterate(reader::GFAreader)
   # Get the first line 
   
end
