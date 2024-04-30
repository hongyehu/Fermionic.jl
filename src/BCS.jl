
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
    new_G,cs = shuffle_G(G,Is,N)
    region_size = size(Is)[1]
    c = new_G[1:region_size,region_size+1:end]*inv(Matrix{Float64}(I, N-region_size, N-region_size)-new_G[region_size+1:end,region_size+1:end])
    α = (new_G[1:region_size,1:region_size]+c*new_G[region_size+1:end,region_size+1:end]*transpose(c)).*0.5
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
            # print("i: ",i, " j: ", j , " debug: ",α[i,j]*ad(o,i)*ad(o,j))
            term1 += α[i,j]*ad(o,i)*ad(o,j)
            term2 += logβ[i,j]*ada(o,i,j)
            # if β[i,j]!= 0.0
            #     term2 += log(β[i,j])*ada(o,i,j)
            # end
            term3 -= α[i,j]*a(o,i)*a(o,j)
        end
    end
    unnorm_RDM = exp(Matrix(term1))*exp(Matrix(term2))*exp(Matrix(term3));
    return unnorm_RDM, term1, term2, term3
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
                    for kx in 1:L
                        for ky in 1:L
                            # print("kx: ", kx, " ky: ", ky, " acoef: ", a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ) ,"\n")
                            Gtmp += (2.0/L^2)*a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ)*exp(1im*(2π*kx/L*(i-ip)+2π*ky/L*(j-jp)))
                        end
                    end
                    G[x_coor,y_coor] = Gtmp
                end
            end
        end
    end
    return G 
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

using PythonCall
plt = pyimport("matplotlib.pyplot")
L = 4
Gtmp = BCS_G(L,:s);
Gtmp,_ = shuffle_G(Gtmp,[1,2,3,4,19,20,27,28],2*L^2);
sα, sβ = get_α_β(Gtmp,[1,2,3,4,19,20,27,28],2*L^2);
Gtmp = BCS_G(L,:d);
Gtmp,_ = shuffle_G(Gtmp,[1,2,3,4,19,20,27,28],2*L^2);
dα, dβ = get_α_β(Gtmp,[1,2,3,4,19,20,27,28],2*L^2);


ρs, term1s, term2s, term3s = get_full_RDM(sα,sβ,8);
ρd, term1d, term2d, term3d = get_full_RDM(dα,dβ,8);



plt.imshow(real(ρd))
plt.show()

plt.imshow(real(get_full_RDM(α,β,8)))
plt.show()