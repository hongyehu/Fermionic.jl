function createRowSwapMatrix(n, i, j)
    # Create an n x n identity matrix
    E = Matrix{Float64}(I, n, n)
    
    # Swap rows i and j
    E[i, :], E[j, :] = E[j, :], E[i, :]
    
    return E
end
function shuffle_G(G,Is::Vector{Int64},N::Int64)
    """
    Is: contains the index of region that will not be traced out,
    N: total number of modes/regions
    Return:
    G: shuffled G matrix
    cs: shuffled c vector, new c_dag list
    """
    cs = [i for i in 1:N]
    region_size = size(Is)[1]
    swap_position = 1
    for i in 1:region_size
        if Is[i] != i
            U = createRowSwapMatrix(N,swap_position,Is[i])
            G = U*G*U
            cs = U*cs
            swap_position += 1
        else
            swap_position += 1
        end
    end
    return G,cs
end

function get_α_β(G,Is::Vector{Int64},N::Int64)
    # N is the number of spinless fermion
    # new_G,cs = shuffle_G(G,Is,N)
    new_G = G
    region_size = size(Is)[1]
    c = new_G[1:region_size,region_size+1:end]*inv(Matrix{Float64}(I, N-region_size, N-region_size)-real(new_G[region_size+1:end,region_size+1:end]))
    # print("debug: ","\n")
    # display(Matrix{Float64}(I, N-region_size, N-region_size)-new_G[region_size+1:end,region_size+1:end])
    # print("a11: ",  "\n")
    # display(round.(new_G[1:region_size,1:region_size],digits=4))
    # print("a12: ", "\n")
    # display(round.(new_G[1:region_size,region_size+1:end],digits=4))
    # α = (new_G[1:region_size,1:region_size]+c*new_G[region_size+1:end,region_size+1:end]*transpose(c)).*0.5
    α = (new_G[1:region_size,1:region_size]-c*new_G[region_size+1:end,region_size+1:end]*transpose(c)).*0.5
    β = c*transpose(c)
    return α, β
end
function logmat(A::Matrix{ComplexF64})
    """
    Compute the matrix logarithm of a complex hermitian matrix
    """
    V, U = eigen(A)
    newV = Vector{Float64}()
    for i in 1:size(V)[1]
        if real(V[i]) < 10^-10
            push!(newV,0.0)
        else
            push!(newV,log(real(V[i])))
        end
    end
    return U*diagm(newV)*inv(U)
end
function get_full_RDM(α::Matrix{ComplexF64},β::Matrix{ComplexF64},Nsub::Int)
    # Nsub is the number of spinless fermion
    logβ = logmat(β);
    o = Op(Nsub);
    term1 = zeros(ComplexF64,2^Nsub,2^Nsub);
    term2 = zeros(ComplexF64,2^Nsub,2^Nsub);
    term3 = zeros(ComplexF64,2^Nsub,2^Nsub);
    for i in 1:size(α)[1]
        for j in 1:size(α)[2]
            term1 += α[i,j]*ad(o,i)*ad(o,j)
            term2 += logβ[i,j]*ada(o,i,j)
            term3 -= α[i,j]*a(o,i)*a(o,j)
        end
    end
    unnorm_RDM = exp(Matrix(term1))*exp(Matrix(term2))*exp(Matrix(term3));
    unnorm_RDM = 0.5.*(unnorm_RDM+transpose(unnorm_RDM));
    return unnorm_RDM./real(tr(unnorm_RDM))
end

function BCS_G(L::Int64, type::Symbol;μ::Float64=0.5,Δ::Float64=5.0)
    G = zeros(ComplexF64,2*L^2,2*L^2)
    for i in 1:L
        for j in 1:L
            for ip in 1:L
                for jp in 1:L
                    x_coor = Cartesian2Index([i,j],[L,L],1)
                    y_coor = Cartesian2Index([ip,jp],[L,L],2)
                    Gtmp = 0.0+0.0im
                    # for kx in Int(-(L-1)/2):1:Int((L-1)/2)#1:L
                    #     for ky in Int(-(L-1)/2):1:Int((L-1)/2)
                    for kx in 1:L
                        for ky in 1:L
                            # print("kx: ", kx, " ky: ", ky, " acoef: ", a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ) ,"\n")
                            Gtmp += a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ)*exp(1im*(2π*kx/L*(i-ip)+2π*ky/L*(j-jp)))/(L^2)
                        end
                    end
                    G[x_coor,y_coor] = Gtmp/2.0
                    # See Minh's code
                    G[y_coor,x_coor] = -Gtmp/2.0
                end
            end
        end
    end
    return G 
end
function getBCSaij(L::Int64, type::Symbol;μ::Float64=0.5,Δ::Float64=5.0)
    # return the matrix a_ij in Eq. (3) of https://journals.aps.org/prb/pdf/10.1103/PhysRevB.38.931
    # transformed from Minh's code, huge thanks to Minh!
    N = 2 * L^2
    ## ak in momentum space
    # kxlist = 2 * π / L * collect(-(L-1)/2:1:(L-1)/2)
    kxlist = 2 * π / L * collect(1:L)
    kylist = kxlist
    ak = zeros(Complex{Float64}, L, L)
    for i = 1:L
        for j = 1:L
            kx = kxlist[i]
            ky = kylist[j]
            ak[i, j] = getak(kx, ky, type;μ=μ,Δ=Δ)
        end
    end
    ## convert to real space
    aij = zeros(Complex{Float64}, N, N)
    for i = 1:N
        for j = 1:N
            if i == j
                aij[i, j] = 0
            else
                xi, yi, si = index2coordinates(L, i)
                xj, yj, sj = index2coordinates(L, j)
                if si == 0 && sj == 1
                    rx = xi - xj
                    ry = yi - yj
    
                    a = 0
                    for kx in kxlist
                        for ky in kylist
                            a += getak(kx, ky, type;μ=μ,Δ=Δ) * exp(1im * (kx * rx + ky * ry)) / L^2
                        end
                    end
    
                    aij[i, j] = a / 2
                    aij[j, i] = -a / 2
                end
            end
        end
    end

    return aij
end
function index2coordinates(L, i)
    s = (i + 1) % 2
    x = ((i - s + 1) ÷ 2) % L
    x = x == 0 ? L : x
    y = ((i - s + 1) ÷ 2 - x) ÷ L + 1
    return x, y, s
end
function getak(kx::Float64, ky::Float64, type::Symbol ;μ::Float64=0.5,Δ::Float64=5.0)
    Delta = Δ
    t = 1.0
    mu = μ
    Ek = -2 * t * (cos(kx) + cos(ky)) - mu
    Deltak = if type == :d
                Delta * (cos(kx) - cos(ky))
            else
                Delta
            end
    if abs(Deltak)  < 10^-10
        ak = 0.0
    else
        ak = Deltak / (Ek + sqrt(Ek^2 + Deltak^2))
    end
    return ak
end

function Δfun(kx::Float64,ky::Float64, type::Symbol;Δ::Float64=5.0)
    if type == :s
        return Δ
    elseif type == :d
        return Δ*(cos(kx)-cos(ky))
    else
        error("type not implemented")
    end
end
function ξfun(kx::Float64,ky::Float64;μ::Float64=5.0)
    return -2*(cos(kx)+cos(ky))-μ
end
function a_coef(kx::Float64,ky::Float64, type::Symbol;μ::Float64=0.5,Δ::Float64=5.0)
    if abs(Δfun(kx,ky,type;Δ=Δ)) < 10^-10
        return 0.0
    else
        return Δfun(kx,ky,type;Δ=Δ)/(sqrt(Δfun(kx,ky,type;Δ=Δ)^2+ξfun(kx,ky;μ=μ)^2)+ξfun(kx,ky;μ=μ))
    end
end
function RDM(G::Matrix{ComplexF64},sites::Vector{Vector{Int64}})
    # L is the lattice size of the two dimensional lattice (spinfull-fermions)
    # sites is the coordinates of the sites that we want to keep[[x1,y1],[x2,y2],...]
    L = Int64(sqrt(size(G)[1]/2))
    Is = Vector{Int64}()
    for i in 1:size(sites)[1]
        push!(Is,Cartesian2Index([sites[i][1],sites[i][2]],[L,L],1))
        push!(Is,Cartesian2Index([sites[i][1],sites[i][2]],[L,L],2))
    end
    Gtmp,cs = shuffle_G(G,Is,2*L^2);
    α, β = get_α_β(Gtmp,Is,2*L^2);
    # print("α: ",round.(α,digits=4),"\n")
    # print("β: ",round.(β,digits=4),"\n")
    ρ = get_full_RDM(α,β,size(Is)[1]);
    return ρ
end

function RDM(G::Matrix{ComplexF64},Is::Vector{Int64})
    # L is the lattice size of the two dimensional lattice (spinfull-fermions)
    # sites is the coordinates of the sites that we want to keep[[x1,y1],[x2,y2],...]
    L = Int64(sqrt(size(G)[1]/2))
    Gtmp,cs = shuffle_G(G,Is,2*L^2);

    α, β = get_α_β(Gtmp,Is,2*L^2);
    ρ = get_full_RDM(α,β,size(Is)[1]);
    return ρ
end
####### Particle-Hole PH_transformation ############
# We need a function that takes the sites::Matrix{Int64} and returns the sites
# that needs to be applied to PH transformation together with the signs
function PH_transformation_indices(L::Int64,sites::Vector{Vector{Int64}})
    # L is the lattice size of the two dimensional lattice (spinfull-fermions)
    # sites is the coordinates of the sites that we want to keep[[x1,y1],[x2,y2],...]
    # Return:
    # Is: the indices of the sites that we want to keep (one dimensional indices)
    # Transform: 0 for no transformation, 1 for PH transformation (spin-up no transform, spin-down transform)
    # Signs: 1 for no transformation (or transform with 1.0 sign), -1 for PH transformation
    Is = Vector{Int64}()
    Transform = Vector{Int64}()
    Signs = Vector{Float64}()
    for i in 1:size(sites)[1]
        push!(Is,Cartesian2Index([sites[i][1],sites[i][2]],[L,L],1))
        push!(Transform,0) #Don't transform
        push!(Signs,1.0)
        push!(Is,Cartesian2Index([sites[i][1],sites[i][2]],[L,L],2))
        push!(Transform,1) #PH Transform
        if mod(sites[i][1]+sites[i][2],2) == 0
            push!(Signs,1.0)
        else
            push!(Signs,-1.0)
        end
    end
    return Is,Transform,Signs
end
function get_full_RDM_PH(α::Matrix{ComplexF64},β::Matrix{ComplexF64},Nsub::Int,Transform::Vector{Int64},Signs::Vector{Float64})
    # Nsub is the number of spinless fermion
    logβ = logmat(β);
    o = Op(Nsub);
    term1 = zeros(ComplexF64,2^Nsub,2^Nsub);
    term2 = zeros(ComplexF64,2^Nsub,2^Nsub);
    term3 = zeros(ComplexF64,2^Nsub,2^Nsub);
    for i in 1:size(α)[1]
        for j in 1:size(α)[2]
            ##
            if Transform[i]==0&&Transform[j]==0
                term1 += α[i,j]*ad(o,i)*ad(o,j)
                term2 += logβ[i,j]*ad(o,i)*a(o,j)
                term3 -= α[i,j]*a(o,i)*a(o,j)
            elseif Transform[i]==1&&Transform[j]==1
                term1 += α[i,j]*Signs[i]*Signs[j]*a(o,i)*a(o,j)
                term2 += logβ[i,j]*Signs[i]*Signs[j]*a(o,i)*ad(o,j)
                term3 -= α[i,j]*Signs[i]*Signs[j]*ad(o,i)*ad(o,j)
            elseif Transform[i]==0&&Transform[j]==1
                term1 += α[i,j]*Signs[j]*ad(o,i)*a(o,j)
                term2 += logβ[i,j]*Signs[j]*ad(o,i)*ad(o,j)
                term3 -= α[i,j]*Signs[j]*a(o,i)*ad(o,j)
            elseif Transform[i]==1&&Transform[j]==0
                term1 += α[i,j]*Signs[i]*a(o,i)*ad(o,j)
                term2 += logβ[i,j]*Signs[i]*a(o,i)*a(o,j)
                term3 -= α[i,j]*Signs[i]*ad(o,i)*a(o,j)
            end
            ##
        end
    end
    unnorm_RDM = exp(Matrix(term1))*exp(Matrix(term2))*exp(Matrix(term3));
    unnorm_RDM = 0.5.*(unnorm_RDM+transpose(unnorm_RDM));
    return unnorm_RDM./real(tr(unnorm_RDM))
end
"""
Particle-hole transformed reduced density matrix
"""
function RDM_PH(G::Matrix{ComplexF64},sites::Vector{Vector{Int64}})
    # L is the lattice size of the two dimensional lattice (spinfull-fermions)
    # sites is the coordinates of the sites that we want to keep[[x1,y1],[x2,y2],...]
    L = Int64(sqrt(size(G)[1]/2))
    Is, Transform, Signs = PH_transformation_indices(L, sites)
    Gtmp,cs = shuffle_G(G,Is,2*L^2);
    α, β = get_α_β(Gtmp,Is,2*L^2);
    ρ = get_full_RDM_PH(α,β,size(Is)[1],Transform,Signs);
    return ρ
end

# We need a function to generate the measurement observable in Attractive region
"""
left and right are the paddel index within flower: 1,2,3,4
"""
function two_paddle_obs(;left::Int64, right::Int64)
    o = Op(8);
    if left == 1
        Δdag_left = -ada(o,3,2)+ada(o,1,4)
    elseif left == 2
        Δdag_left = ada(o,1,4)-ada(o,3,2)
    elseif left == 3
        Δdag_left = ada(o,3,2)-ada(o,1,4)
    elseif left == 4
        Δdag_left = -ada(o,1,4)+ada(o,3,2)
    end
    O_left = Δdag_left+Δdag_left'
    if right == 1
        Δdag_right = -ada(o,7,6)+ada(o,5,8)
    elseif right == 2
        Δdag_right = ada(o,5,8)-ada(o,7,6)
    elseif right == 3
        Δdag_right = ada(o,7,6)-ada(o,5,8)
    elseif right == 4
        Δdag_right = -ada(o,5,8)+ada(o,7,6)
    end
    O_right = Δdag_right+Δdag_right'
    return O_left*O_right
end