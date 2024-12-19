using Fermionic
using PythonCall
using LinearAlgebra
########### Compare with analytical result on <c*_{0,0,up}c*_{xi,yi,down}>
function analytical_two_fermion(xi::Int64, yi::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    for kx in 0:(L-1)
        for ky in 0:(L-1)
            res += conj(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))*exp(1im*(2π*(kx*xi+ky*yi)/L))/(abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
        end
    end
    return res/(L^2)
end

analytical_two_fermion(0,1,10,:d;μ=0.5,Δ=0.1)
analytical_two_fermion(0,-1,10,:d;μ=0.5,Δ=0.1)



########### Compare with analytical result on <c*_{0,0,up}c*_{xj,yj,down}c_{xk,yk,up}c_{xl,yl,down}>
function analytical_four_fermion(xj::Int64,yj::Int64,xk::Int64,yk::Int64,xl::Int64,yl::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    for kx in 0:(L-1)
        for ky in 0:(L-1)
            for px in 0:(L-1)
                for py in 0:(L-1)
                    coef1 = conj(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))/(abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
                    coef2 = a_coef(2π*px/L,2π*py/L,type;μ=μ,Δ=Δ)/(abs(a_coef(2π*px/L,2π*py/L,type;μ=μ,Δ=Δ))^2+1.0)
                    coef3 = (abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2)/(abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
                    coef4 = (abs(a_coef(2π*px/L,2π*py/L,type;μ=μ,Δ=Δ))^2)/(abs(a_coef(2π*px/L,2π*py/L,type;μ=μ,Δ=Δ))^2+1.0)
                    res += coef1*coef2*exp(1im*(2π*kx*xj/L+2π*ky*yj/L+2π*px*(xk-xl)/L+2π*py*(yk-yl)/L))
                    res += coef3*coef4*exp(1im*(2π*kx*xk/L+2π*ky*yk/L+2π*px*(xl-xj)/L+2π*py*(yl-yj)/L))
                end
            end
        end
    end
    return -res./(L^4)
end

ori = analytical_four_fermion(0,1,0,-5,0,-6,10,:d;μ=0.5,Δ=0.1)
singlet_sym = analytical_four_fermion(0,1,0,-5,0,-6,10,:d;μ=0.5,Δ=0.1)+analytical_four_fermion(0,1,0,-6,0,-5,10,:d;μ=0.5,Δ=0.1)+analytical_four_fermion(0,-1,0,-6,0,-7,10,:d;μ=0.5,Δ=0.1)+analytical_four_fermion(0,-1,0,-7,0,-6,10,:d;μ=0.5,Δ=0.1)
triplet_sym = analytical_four_fermion(0,1,0,-5,0,-6,10,:d;μ=0.5,Δ=0.1)-analytical_four_fermion(0,1,0,-6,0,-5,10,:d;μ=0.5,Δ=0.1)-analytical_four_fermion(0,-1,0,-6,0,-7,10,:d;μ=0.5,Δ=0.1)+analytical_four_fermion(0,-1,0,-7,0,-6,10,:d;μ=0.5,Δ=0.1)
triplet_sym
singlet_sym




########### Analytical result on <c*_{0,0,up}c_{0,0,up}>, and system has translational invariance
function analytical_n_up(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    for kx in 0:(L-1)
        for ky in 0:(L-1)
            res += abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2/(abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
        end
    end
    return res/(L^2)
end
########### Analytical result on <c*_{0,0,down}c_{0,0,down}>, and system has translational invariance
function analytical_n_down(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    for kx in 0:(L-1)
        for ky in 0:(L-1)
            res += abs(a_coef(-2π*kx/L,-2π*ky/L,type;μ=μ,Δ=Δ))^2/(abs(a_coef(-2π*kx/L,-2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
        end
    end
    return res/(L^2)
end
##### Analytical result on <c*(0,up)c*(x,y,down)+c*(x,y,up)c*(0,0,down)>
function SpinSymmetricPairWave(x::Int64,y::Int64,L::Int64,type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    for kx in 0:(L-1)
        for ky in 0:(L-1)
            res += conj(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))*2*cos(2π*(kx*x+ky*y)/L)/(abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
        end
    end
    return res/(L^2) 
end
##### Analytical result on <c*(0,up)c*(x,y,down)+c*(x,y,up)c*(0,0,down)>
# function SpinAntisymmetricPairWave(x::Int64,y::Int64,L::Int64,type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
#     res = 0.0+0.0im
#     for kx in 0:(L-1)
#         for ky in 0:(L-1)
#             res += conj(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))*2*cos(2π*(kx*x+ky*y)/L)/(abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
#         end
#     end
#     return res/(L^2) 
# end
## Plot the pairing function
# begin
#     L=31
#     μ=0.2
#     Δ=0.2
#     type = :s
#     Test_s = zeros(Complex{Float64},L,L)
#     for x in 1:L
#         for y in 1:L
#             Test_s[x,y] = SpinSymmetricPairWave(x,y,L,:s,μ=μ,Δ=Δ)
#         end
#     end
#     type=:d
#     Test_d = zeros(Complex{Float64},L,L)
#     for x in 1:L
#         for y in 1:L
#             Test_d[x,y] = SpinSymmetricPairWave(x,y,L,:d,μ=μ,Δ=Δ)
#         end
#     end
#     center = (L + 1) ÷ 2
#     # Calculate the shift needed to move [1, 1] to the center
#     shift_row = center 
#     shift_col = center 
#     # Use circshift to shift the matrix
#     shifted_Test_s = circshift(Test_s, (shift_row, shift_col))
#     shifted_Test_d = circshift(Test_d, (shift_row, shift_col))
#     begin
#         # vmin = minimum([minimum(abs.(s_real_space)), minimum(abs.(d_real_space))])
#         # vmax = maximum([maximum(abs.(s_real_space)), maximum(abs.(d_real_space))])
#         fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#         img1 = axs[0].imshow(real(shifted_Test_s),cmap="Oranges_r")
#         axs[0].set_title("s-wave")
#         img2 = axs[1].imshow(real(shifted_Test_d),cmap="Oranges_r")
#         axs[1].set_title("d-wave")
#         fig.colorbar(img1, ax=axs[0])
#         fig.colorbar(img2, ax=axs[1])
#         plt.show()
#         # fig.savefig("real_space_pairing_strong_limit.pdf", bbox_inches="tight")
#     end
# end

# analytical_two_fermion(0,0,23,:d;μ=0.2,Δ=0.2)*2
# SpinSymmetricPairWave(0,0,23,:d;μ=0.2,Δ=0.2)

# analytical_n_up(17,:s;μ=0.0,Δ=0.1)
# analytical_n_down(17,:s;μ=0.0,Δ=0.1)

# L=25
# μ=0.2
# Δ=0.2
# type=:s
# analytical_n_up(L,type;μ=μ,Δ=Δ)+analytical_n_up(L,type;μ=μ,Δ=Δ)
# so = 