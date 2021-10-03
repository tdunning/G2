module G2

using DataStructures

"""
Computes the G² test statistic for a contingency table in a 2x2 matrix

This test is a χ² test that avoids the assumption of normal distribution
inherent in Pearson's χ² test. This is very helpful for sparse problems.

See http://tdunning.blogspot.com/2008/03/surprise-and-coincidence.html and https://www.researchgate.net/publication/2477641_Accurate_Methods_for_the_Statistics_of_Surprise_and_Coincidence for more information.
"""
function g2(m::Matrix)::Float64
    2 * sum(m) * (entropy(sum(m, dims=1)) + entropy(sum(m, dims=2)) - entropy(m))
end

"""
Compares the counts of the entries in two or more tables and returns a single
value indicating how incompatible the counts in each table are.
"""
function g2(d::Vector{AbstractDict}...)
    # copy to a matrix
    k = union(keys.(d...)...)
    m = zeros(Int32, length(k), length(d))
    i = 1
    for x in k
        j = 1
        for dx in d
            if k in dx
	        m[i, j] = dx[k]
            end
	    j += 1
	end
	i += 1
    end
    g2(m)
end

"""
Computes the signed square root of the G² statistic for a 2x2 matrix
of counts. The result has positive sign if the upper left count is
bigger than expected.

This is useful when looking primarily for things that are over-represented
against a background level of counts or when comparing frequencies between
two sources where we want to remember whether the rate was higher in one
or another collection and also want an indication of whether that was 
interesting.
"""
function g2root(m::Matrix)
    @assert size(m) == (2, 2)
    r = g2(m)
    s = sum(m, dims=1)
    if m[1,1] * (s[1] + s[2]) >= s[1] *  (m[1, 1] + m[1, 2])
        return sqrt(r)
    else
        return -sqrt(r)
    end
end

"""
Compares counts in two tables to find which elements have interestingly
different counts. The result has one entry for any key in either of 
the source tables and a value equal to the g2root score related to that
key relative to all other counts.

This is done by computing a 2x2 table for each key that contains the counts
for that key for each input table in one column and has the total counts
for all other keys in each table in the second column. The result is
roughly calibrated in terms of standard deviation (ISH!) and is positive
if the frequency of a particular key is greater in d1 than in d2 and negative
if the frequency is higher in d2.
"""
function compare(d1::AbstractDict, d2::AbstractDict)::Dict
    kx = union(keys(d1), keys(d2))
    n1 = sum(d1)
    n2 = sum(d2)
    Dict(
        let k11 = get(d1, k, 0)
            k21 = get(d2, k, 0);
            k => g2root([k11 n1 - k11; k21 n2 - k21])
	end
        for k in kx)
end
    

"""
Lazily compares counts in two tables to find which elements have interestingly
different counts. The result has one entry for any key in either of 
the source tables and a value equal to the g2root score related to that
key relative to all other counts.

This is done by computing a 2x2 table for each key that contains the counts
for that key for each input table in one column and has the total counts
for all other keys in each table in the second column. The result is
roughly calibrated in terms of standard deviation (ISH!) and is positive
if the frequency of a particular key is greater in d1 than in d2 and negative
if the frequency is higher in d2.
"""
function lazy(d1::AbstractDict{K,V}, d2::AbstractDict{K,V}) where {K, V <: Number}
    n1 = sum(d1)
    n2 = sum(d2)
    DefaultDict(passkey=true) do k
        k11 = get(d1, k, 0)
        k21 = get(d2, k, 0);
	if k11 == k21 == 0
	    throw(KeyError(k))
	else
	    g2root([k11 n1 - k11; k21 n2 - k21])
	end
    end
end
    




"""
Computes entropy of counts
"""
function entropy(v)
    s = sum(v)
    sum(-v ./ s .* log.(v/s + (v .== 0)))
end

end # module
