L = 2
type=:d
G = BCS_G(L,type)
region = [[1,1],[2,1]]
ρ = RDM(G,region)
ρp = RDM_PH(G,region)

o = Op(4)
obs1 = ad(o,1)*ad(o,4)
obs2 = -ad(o,1)*a(o,4)

tr(ρ*obs1)
tr(ρp*obs2)
