using SymPhase

c = parse_stim("./GHZ.stim")
sampler = Sampler(c)
r = sampler(64)
lens, nm = size(r)
tm = typemax(eltype(r))
for j2 in 1:nm-1
    for j1 in 1:lens
        @test r[j1,j2]+r[j1,j2+1] == tm
    end
end