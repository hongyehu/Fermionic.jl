using Fermionic
function Diag_real_skew(M, rand_perturbation::Int64=0)
    N = div(size(M, 1), 2)
  
    # Random perturbation before forcing skew symmetrisation
    if (rand_perturbation != 0)
      if (rand_perturbation == 1)
          random_M = rand(2N, 2N) * eps()
          random_M = (random_M - random_M') / 2.0
          M += random_M
        end
        if (rand_perturbation == 4)
          random_M = zeros(Complex{Float64}, 2N, 2N)
          random_M[1, N+2] = eps()
          random_M[2, N+1] = -random_M[1, N+2]
          random_M = (random_M - random_M') / 2.0
          M += random_M
        end
        if (rand_perturbation == 5)
          random_M = zeros(Complex{Float64}, 2N, 2N)
          r = eps()
          M[1, 2] += r
          M[2, 1] -= r
        end
    end
  
    Schur_object = LinearAlgebra.schur(M)
  
    Schur_ort_i = Schur_object.vectors
    Schur_blocks_i = Schur_object.Schur
  
    # Reorder so that all 2x2 blocks have positive value in top right, accounting for 1x1 blocks.
    swap_perm = Vector{Int64}(undef, 0)
    sizehint!(swap_perm, 2N)
    iiter = 1
    while iiter < 2N
        if abs(Schur_blocks_i[iiter, iiter]) >= abs(Schur_blocks_i[iiter+1, iiter]) # alternative: s_b_i[i+1, i] == 0
            # We have a 1x1 block - move it to the front. These always come in pairs, but not always sequentially.
            pushfirst!(swap_perm, iiter)
            iiter += 1
        elseif (Schur_blocks_i[iiter + 1, iiter] >= 0.0)
            # Flipped 2x2 block
            push!(swap_perm, iiter + 1, iiter)
            iiter += 2
        else
            # Unflipped 2x2 block
            push!(swap_perm, iiter, iiter + 1)
            iiter += 2
        end
    end
    if iiter == 2N  # catch the final 1x1 block if it exists
      pushfirst!(swap_perm, 2N)
    end
    
    M_temp = Schur_blocks_i[swap_perm, swap_perm]
    O_temp = Schur_ort_i[:, swap_perm]
  
    # Sort the blocks, λ_1>=λ_2>=...>=λ_N with λ_1 the coefficient in the upper left block
    psort = sortperm(diag(M_temp, 1)[begin:2:end], rev=true)
    full_psort = zeros(Int64, 2N)
    full_psort[begin:2:end] .= 2 .* psort .- 1
    full_psort[begin + 1:2:end] .= 2 .* psort
  
    M_f = M_temp[full_psort, full_psort]
    O_f = real(O_temp[:, full_psort])
  
    return M_f, O_f
  end

using Random

function random_antisymmetric_matrix(n::Int64)
    # Initialize a n x n matrix with zeros
    A = zeros(Float64, n, n)
    
    # Fill the upper triangular part above the diagonal with random values
    for i in 1:n
        for j in (i+1):n
            A[i, j] = randn()  # Generate a random number from a normal distribution
            A[j, i] = -A[i, j]  # Set the symmetric element to the negative
        end
    end
    
    return A
end
test = random_antisymmetric_matrix(200)
M, O = Diag_real_skew(test);
round.(M,digits=2)

Cartesian2Index([1,1],[2,2],2)

test = createRowSwapMatrix(5,1,2)
test*test
tmp = [1,2,3,4,5]
size(tmp)[1]