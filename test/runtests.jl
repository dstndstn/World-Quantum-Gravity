using WorldQuantumGravity
using Test

@testset "Labels and coordinates" begin
    @test all(lbl(crd(m)) == m for m in 1:gridSize(3))
end

@testset "Atoms" begin
    @test atom(Coord((0,0,0))) == [
        Coord{3}((0, 0, 0)),
        Coord{3}((1, 0, 0)),
        Coord{3}((0, 1, 0)),
        Coord{3}((1, 1, 0)),
        Coord{3}((0, 0, 1)),
        Coord{3}((1, 0, 1)),
        Coord{3}((0, 1, 1)),
        Coord{3}((1, 1, 1)),
    ]

    @test atom(3) == [3, 4, 6, 7, 12, 13, 15, 16]
end

@testset "Edges" begin
    @test atomCorner(2,1) == [Edgevv(2, 3), Edgevv(2, 5), Edgevv(2, 11)]

    @test atomCorner(Coord((0,0,1)),1) == [
        (Coord{3}((0, 0, 1)), Coord{3}((1, 0, 1))),
        (Coord{3}((0, 0, 1)), Coord{3}((0, 1, 1))),
        (Coord{3}((0, 0, 1)), Coord{3}((0, 0, 2))),
    ]

    @test direction(Edgevv(1,10)) == 3

    @test edgeVD(Edgevv(1,4)) == Edgevd(1, 2)
end

@testset "Amplitudes" begin
    svalue=fill(1.0,(gridSize(length(listGridSize)),length(listGridSize)))
    cvalue=fill(0.0,(gridSize(length(listGridSize)),length(listGridSize)))

    @test ampEdge(svalue,cvalue,Edgevd(2,1),1) ≈ 1.7155277699214135

    @test expCorner(svalue,2,2,2) ≈ 2.0

    @test ampCorner(svalue,cvalue,2,2,1,0) ≈ 22.627416997969515

    @test ampVG(svalue, cvalue, 1,2) ≈ 4.209169856678294e6
end
