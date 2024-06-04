using Fermionic
using PythonCall
using LinearAlgebra
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