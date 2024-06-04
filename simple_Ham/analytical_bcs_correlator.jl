using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")
include("./analytical_bcs.jl")

"""
The following function calculate true pairing correlation using analytical result.
"""
function s_pairing_correlation(l::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    res += analytical_four_fermion(1,0,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(1,0,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(1,0,mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(1,0,mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res += analytical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(mod(-1,L),0,mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(mod(-1,L),0,mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res += analytical_four_fermion(0,1,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,1,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,1,mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,1,mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res += analytical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,mod(-1,L),mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,mod(-1,L),mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    return -res
end
function single_site_s_pairing_correlation(l::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    res += analytical_four_fermion(0,0,mod(l,L),0,mod(l,L),0,L,type;μ=μ,Δ=Δ)
    return -res
end

"""
The following function calculate true pairing correlation using analytical result.
"""
function d_pairing_correlation(l::Int64, L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    res = 0.0+0.0im
    res += analytical_four_fermion(1,0,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(1,0,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(1,0,mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(1,0,mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res += analytical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(mod(-1,L),0,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(mod(-1,L),0,mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(mod(-1,L),0,mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res -= analytical_four_fermion(0,1,mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(0,1,mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,1,mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,1,mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    ##
    res -= analytical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l+1,L),0,L,type;μ=μ,Δ=Δ)
    res -= analytical_four_fermion(0,mod(-1,L),mod(l,L),0,mod(l-1,L),0,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,mod(-1,L),mod(l,L),0,l%L,1,L,type;μ=μ,Δ=Δ)
    res += analytical_four_fermion(0,mod(-1,L),mod(l,L),0,l%L,mod(-1,L),L,type;μ=μ,Δ=Δ)
    return -res
end

L = 23
s_single_correlation_4_swave = [single_site_s_pairing_correlation(l,L,:s;μ=3.0,Δ=0.1) for l in 1:Int((L-1)/2)]
s_single_correlation_4_dwave = [single_site_s_pairing_correlation(l,L,:d;μ=3.0,Δ=0.1) for l in 1:Int((L-1)/2)]
s_correlation_4_swave = [s_pairing_correlation(l,L,:s;μ=3.0,Δ=0.1) for l in 1:Int((L-1)/2)]
s_correlation_4_dwave = [s_pairing_correlation(l,L,:d;μ=3.0,Δ=0.1) for l in 1:Int((L-1)/2)]
d_correlation_4_swave = [d_pairing_correlation(l,L,:s;μ=0.2,Δ=0.1) for l in 1:Int((L-1)/2)]
d_correlation_4_dwave = [d_pairing_correlation(l,L,:d;μ=0.2,Δ=0.1) for l in 1:Int((L-1)/2)]
begin
    f = plt.figure(figsize=(4,2.5))
    plt.plot(1:Int((L-1)/2),abs.(real(s_single_correlation_4_swave)),marker = ".",label="s-wave")
    plt.plot(1:Int((L-1)/2),abs.(real(s_single_correlation_4_dwave)),marker = ".",label="d-wave")
    plt.xlabel("l")
    plt.ylabel(raw"$\langle \Delta^{\dagger}_{\text{onsite}}(0)\Delta_{\text{onsite}}(l)\rangle$")
    plt.yscale("log")
    plt.legend()
    plt.show()
    # f.savefig("onsite_s_correlation.pdf", bbox_inches="tight")
end

begin
    f = plt.figure(figsize=(4,2.5))
    plt.plot(1:Int((L-1)/2),abs.(real(s_correlation_4_swave)),marker = ".",label="s-wave")
    plt.plot(1:Int((L-1)/2),abs.(real(s_correlation_4_dwave)),marker = ".",label="d-wave")
    plt.xlabel("l")
    plt.ylabel(raw"$\langle \Delta^{\dagger}_{s}(0)\Delta_{s}(l)\rangle$")
    plt.yscale("log")
    plt.legend()
    plt.show()
    # f.savefig("s_correlation.pdf", bbox_inches="tight")
end

begin
    f = plt.figure(figsize=(4,2.5))
    plt.plot(1:Int((L-1)/2),abs.(real(d_correlation_4_swave)),marker = ".",label="s-wave")
    plt.plot(1:Int((L-1)/2),abs.(real(d_correlation_4_dwave)),marker = ".",label="d-wave")
    plt.xlabel("l")
    plt.ylabel(raw"$\langle \Delta^{\dagger}_{d}(0)\Delta_{d}(l)\rangle$")
    plt.yscale("log")
    plt.legend()
    plt.show()
    # f.savefig("d_correlation.pdf", bbox_inches="tight")
end

