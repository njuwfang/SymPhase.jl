using SafeTestsets
using SymPhase

@safetestset "Misc" begin
    @safetestset "Aqua" begin include("Aqua.jl") end
end

@safetestset "SymPhase.jl" begin
    @safetestset "GHZ" begin include("test_GHZ.jl") end
end

