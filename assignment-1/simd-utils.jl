
@inline function range_chunks(r, n)
    rv = @inbounds r[1:n:end-(n-1)]
    rr = @inbounds r[length(rv)*n+1:end]
    rv, rr
end
