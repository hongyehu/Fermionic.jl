using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")


"""
The following function calculate correlation of two pair of fermions that is
    either aligned or orthogonal to each other.
"""
function BCS_pairing(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0,rotation = true,connect=true)
    if type == :s
        G = BCS_G(L, :s;μ=μ,Δ=Δ);
    elseif type == :d
        G = BCS_G(L, :d;μ=μ,Δ=Δ);
    else
        error("type not implemented")
    end
    o = Op(8);
    Δcre = ad(o,1)*ad(o,4);
    Δana = a(o,5)*a(o,8);
    if rotation==false
        res = Vector{ComplexF64}(undef, L-1);
        obs = Δcre*Δana;
        for i in 1:L-1
            Is = [Cartesian2Index([1,1],[L,L],1),
                Cartesian2Index([1,1],[L,L],2),
                Cartesian2Index([2,1],[L,L],1),
                Cartesian2Index([2,1],[L,L],2),
                Cartesian2Index([1,1+i],[L,L],1),
                Cartesian2Index([1,1+i],[L,L],2),
                Cartesian2Index([2,1+i],[L,L],1),
                Cartesian2Index([2,1+i],[L,L],2)];
            ρs = RDM(G, Is);
            if connect==true
                res[i] = tr(ρs*obs)-tr(ρs*Δcre)*tr(ρs*Δana);
            else
                res[i] = tr(ρs*obs);
            end
        end
    else
        obs = Δcre*Δana;
        res = Vector{ComplexF64}(undef, L-2);
        for i in 1:L-2
            Is = [Cartesian2Index([1,1],[L,L],1),
                Cartesian2Index([1,1],[L,L],2),
                Cartesian2Index([2,1],[L,L],1),
                Cartesian2Index([2,1],[L,L],2),
                Cartesian2Index([1,1+i],[L,L],1),
                Cartesian2Index([1,1+i],[L,L],2),
                Cartesian2Index([1,1+i+1],[L,L],1),
                Cartesian2Index([1,1+i+1],[L,L],2)];
            ρs = RDM(G, Is);
            if connect==true
                res[i] = tr(ρs*obs)-tr(ρs*Δcre)*tr(ρs*Δana);
            else
                res[i] = tr(ρs*obs);
            end
        end
    end
    return res
end

begin
    μ=0.5
    Δ=0.2
    L = 15
    res_parallel_s = BCS_pairing(L, :s, rotation=false,connect=false,μ=μ,Δ=Δ);
    res_parallel_d = BCS_pairing(L, :d, rotation=false,connect=false,μ=μ,Δ=Δ);
    res_rot_s = BCS_pairing(L, :s, rotation=true,connect=false,μ=μ,Δ=Δ);
    res_rot_d = BCS_pairing(L, :d, rotation=true,connect=false,μ=μ,Δ=Δ);
end


begin
    plt.figure()
    plt.plot(1:L-1, real.(res_parallel_s), marker="o",label="parallel s")
    plt.plot(1:L-1, real.(res_parallel_d), marker="o",linestyle = "--", label="parallel d")
    plt.plot(1:L-2, real.(res_rot_s),  marker="o",label="rotated s")
    plt.plot(1:L-2, real.(res_rot_d),  marker="o",linestyle = "--",label="rotated d")
    plt.legend()
    plt.show()
end


G = G = BCS_G(L, :s;μ=0.5,Δ=5.0);
begin
    plt.imshow(real.(G))
    plt.show()
end


########### Compare with analytical result on <c*_{xi,yi,up}c*_{xi+l,yi,down}>
function analytical_two_fermion(l::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    for kx in 0:(L-1)
        for ky in 0:(L-1)
            res += conj(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))*exp(1im*(2π*kx*l/L))/(L^2)/(abs(a_coef(2π*kx/L,2π*ky/L,type;μ=μ,Δ=Δ))^2+1.0)
        end
    end
    return res
end
function numerical_two_fermion(l::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    G = BCS_G(L, type;μ=μ,Δ=Δ);
    Is = [Cartesian2Index([1,1],[L,L],1),
        Cartesian2Index([1+l,1],[L,L],1),
        Cartesian2Index([1,1],[L,L],2),
        Cartesian2Index([1+l,1],[L,L],2)];
    ρs = RDM(G, Is);
    o = Op(4);
    obs = ad(o,1)*ad(o,4);
    return tr(ρs*obs)
end
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
# function benchmark_four_fermion(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
#     G = BCS_G(L, type;μ=μ,Δ=Δ);
#     Is = [Cartesian2Index([1,1],[L,L],1),
#         Cartesian2Index([1,1],[L,L],2),
#         Cartesian2Index([1,1],[L,L],1),
#         Cartesian2Index([1,1],[L,L],2),
#         Cartesian2Index([1,1],[L,L],1),
#         Cartesian2Index([1,1],[L,L],2),
#         Cartesian2Index([1,1],[L,L],1),
#         Cartesian2Index([1,1],[L,L],2)];
#     ρs = RDM(G, Is);
#     o = Op(8);
#     obs = ad(o,1)*ad(o,4)*a(o,7)*a(o,8);
#     return tr(ρs*obs)
# end
function benchmark_four_fermion_Minhplot(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = zeros(Complex{Float64}, L, L)
    for i in 1:L
        for j in 1:L
            res[i,j] = analytical_four_fermion(i-1,j-1,0,1,0,0,L,type;μ=μ,Δ=Δ)
        end
    end
    return res
end
L = 13
mat_s=benchmark_four_fermion(L,:s,μ=0.5,Δ=5.0)
mat_d=benchmark_four_fermion(L,:d,μ=0.5,Δ=5.0)
center = (L + 1) ÷ 2
# Calculate the shift needed to move [1, 1] to the center
shift_row = center - 1
shift_col = center - 1

# Use circshift to shift the matrix
shifted_mat_s = circshift(mat_s, (shift_row, shift_col))
shifted_mat_d = circshift(mat_d, (shift_row, shift_col))
begin
	vmin = -0.045
	vmax = 0.03
	fig, axs = plt.subplots(1, 2, figsize=(10, 5))
	img1 = axs[0].imshow(real(shifted_mat_s),cmap="Oranges_r")
	axs[0].set_title("s-wave")
	img2 = axs[1].imshow(real(shifted_mat_d),cmap="Oranges_r")
	axs[1].set_title("d-wave")
    axs[1].set_xlabel("X-axis label")  # Set x-axis label
    axs[1].set_ylabel("Y-axis label") 
	img1.set_clim(vmin, vmax)
	img2.set_clim(vmin, vmax)

	# Show color bar
	fig.colorbar(img1, ax=axs[0])
	fig.colorbar(img2, ax=axs[1])
	plt.show()
end
shifted_mat_d


