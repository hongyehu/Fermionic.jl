
function index2coordinates(L, i)
    s = (i + 1) % 2
    x = ((i - s + 1) ÷ 2) % L
    x = x == 0 ? L : x
    y = ((i - s + 1) ÷ 2 - x) ÷ L + 1
    return x, y, s
end
function getak(kx, ky,type::Symbol)
    Delta = 5.0
    t = 1.0
    mu = 0.5
    Ek = -2 * t * (cos(kx) + cos(ky)) - mu
    Deltak = if type == :d
                Delta * (cos(kx) - cos(ky))
            else
                Delta
            end
    if Deltak == 0
        ak = 0
    else
        ak = Deltak / (Ek + sqrt(Ek^2 + Deltak^2))
    end
    return ak
end
function getBCSaij(L::Int64,type::Symbol)
    # return the matrix a_ij in Eq. (3) of https://journals.aps.org/prb/pdf/10.1103/PhysRevB.38.931
    N = 2 * L^2
    ## ak in momentum space
    kxlist = 2 * π / L * collect(-(L-1)/2:1:(L-1)/2)
    kylist = kxlist
    ak = zeros(Complex{Float64}, L, L)
    for i = 1:L
        for j = 1:L
            kx = kxlist[i]
            ky = kylist[j]
            ak[i, j] = getak(kx, ky, type)
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
                            a += getak(kx, ky, type) * exp(1im * (kx * rx + ky * ry)) / L^2
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



