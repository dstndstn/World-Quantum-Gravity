module WorldQuantumGravity

using GraphPlot
using LightGraphs



# We can add comments like this
export Coord
struct Coord{D}
    c::NTuple{D,Int}
end

Base.getindex(c::Coord, i) = c.c[i]



export listGridSize
const listGridSize = (3,3,5)



export gridGraph
function gridGraph(m)
    sk = path_graph(listGridSize[length(listGridSize)])
    for i in length(listGridSize)-1:-1:1
       sk = cartesian_product(sk , path_graph(listGridSize[i]))
    end
    return sk
end



# grid object

export gridSize
function gridSize(j)
    prod(listGridSize[1:j])
end

export crd
function crd(m::Int)
    d = length(listGridSize)
    c = zeros(Int, d)
    m -= 1
    for i in 1:d
        c[i] = mod(m, listGridSize[i])
        m = fld(m, listGridSize[i])
    end
    @assert m == 0
    Coord(tuple(c...))
end

export lbl
function lbl(c::Coord{D}) where {D}
    m = 0
    for i in D:-1:1
        m = m * listGridSize[i] + c[i]
    end
    m + 1
end



export atom
function atom(c::Coord{D}) where {D}
    c1 = Coord{D}[]
    for offset in CartesianIndices(ntuple(i->2, D))
        push!(c1, Coord(ntuple(i -> c[i] + offset.I[i]-1, D)))
    end
    c1
end

function atom(m::Int)
    [lbl(c) for c in atom(crd(m))]
end



export Edgevv
"""refer to edge by (vertex, vertex)"""
struct Edgevv
    v1::Int
    v2::Int
end

export atomCorner

function atomCorner(m::Int, n::Int)
    # m labels atom, n labels vertex within atom
    d = length(listGridSize)
    crn = Edgevv[]
    for i in 1:2^d
        if sum(abs.([crd(atom(m)[n]).c...]-[crd(atom(m)[i]).c...]))==1
            push!(crn, Edgevv(atom(m)[n], atom(m)[i]))
        end
    end
    crn
end

function atomCorner(c::Coord{D}, n::Int) where {D}
    crn = Tuple{Coord{D},Coord{D}}[]
    for i in 1:2^D
        if sum(abs.([atom(c)[n].c...]-[atom(c)[i].c...]))==1
            push!(crn, (atom(c)[n], atom(c)[i]))
        end
    end
    crn
end



# amplitude

export vvmd
"""dspc = spatial dimension"""
function vvmd(s,c,dspc)
    if s == 0
        return 1.0
    elseif c == 0
        return dspc/s
    else
        return (sinc(s*sqrt(c/dspc)/pi))^(-dspc)
    end
end

export Edgevd
"""refer to edge by (vertex,direction)"""
struct Edgevd
    vert::Int
    dr::Int
end

export ampEdge
function ampEdge(svalue,cvalue,evd::Edgevd,L)
    dspc = length(listGridSize)-1
    m, dr = evd.vert, evd.dr
    if svalue[m,dr] == 0
        return 1
    else
        vd = vvmd(svalue[m,dr],cvalue[m,dr],dspc)
        return vd^(3/vd)*exp(-L*svalue[m,dr]/vd)
    end
end

export direction
function direction(e::Edgevv)
    d = length(listGridSize)
    mc = crd(e.v1)
    nc = crd(e.v2)
    dist = abs.(mc.c .- nc.c)
    for i in 1:d
        if dist == ntuple(j -> i==j, d)
            return i
        end
    end
    @assert false
end

export edgeVD
"""edge (vertices) to edge (vertex,direction)"""
function edgeVD(evv::Edgevv)
    return Edgevd(min(evv.v1, evv.v2),direction(evv))
end

export expCorner
"""(m,n) label corners. in 3d n=1,...,8"""
function expCorner(svalue, m::Int, n::Int, a)
    d = length(listGridSize)
    r = a
    for i in 1:d
        evdi = edgeVD(atomCorner(m,n)[i])
        r *= svalue[evdi.vert, evdi.dr]
    end
    r
end

export ampCorner
function ampCorner(svalue, cvalue, m::Int, n::Int,a,L)
    d = length(listGridSize)
    p = 1.0
    for i in 1:d
        evd = edgeVD(atomCorner(m,n)[i])
        ss = svalue[evd.vert, evd.dr]
        ae = ampEdge(svalue, cvalue, edgeVD(atomCorner(m,n)[i]),L)
        ec = expCorner(svalue, m, n, a)
        p *= ae^(ec / ss)
    end
    return p
end



# total amplitude
export ampVG
function ampVG(svalue, cvalue, a,L)
    d = length(listGridSize)
    lgs = listGridSize  
    las = lgs.- 1    
    p = 1.0
    for i in 1:prod(las), j in 1:2^d
        #albl=atomlabel[i]
        albl = lbl(Coord(CartesianIndices(las)[i].I.-1))
        p *= ampCorner(svalue, cvalue, albl, j, a, L)
    end
    return p
end

end
