using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")

function singlet_pairing_scaling(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    if type == :s
        G = BCS_G(L, :s;μ=μ,Δ=Δ);
    elseif type == :d
        G = BCS_G(L, :d;μ=μ,Δ=Δ);
    else
        error("type not implemented")
    end
    res = Vector{ComplexF64}(undef, L-1);
    o = Op(8);
    Δcre = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3);
    Δana = a(o,8)*a(o,5)-a(o,7)*a(o,6);
    # obs = Δcre*Δana+(Δcre*Δana)';
    obs = Δcre*Δana+Δcre'*Δana';
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
        res[i] = tr(ρs*obs);
    end
    return res
end



"""
The following function is used to calculated (Delta_i^dag Delta_j^dag+Delta_j Delta_i) for s-wave and d-wave pairing.
"""
function symmetrized_singlet_pairing_scaling(L::Int64, type::Symbol; μ::Float64=0.5,Δ::Float64=5.0)
    if type == :s
        G = BCS_G(L, :s;μ=μ,Δ=Δ);
    elseif type == :d
        G = BCS_G(L, :d;μ=μ,Δ=Δ);
    else
        error("type not implemented")
    end
    res = Vector{ComplexF64}(undef, L-1);
    o = Op(8);
    Δidag = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3);
    Δjdag = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7);
    obs = (Δidag+Δidag')*(Δjdag+Δjdag');
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
        res[i] = tr(ρs*obs);
    end
    return res
end


"""
The following function is used to calculated (Delta_i^dag Delta_j^dag+Delta_j Delta_i) for s-wave and d-wave pairing.
"""
function antisymmetrized_singlet_pairing_scaling(L::Int64, type::Symbol)
    if type == :s
        G = BCS_G(L, :s);
    elseif type == :d
        G = BCS_G(L, :d);
    else
        error("type not implemented")
    end
    res = Vector{ComplexF64}(undef, L-1);
    o = Op(8);
    Δidag = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3);
    Δjdag = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7);
    obs = -(Δidag-Δidag')*(Δjdag-Δjdag');
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
        res[i] = tr(ρs*obs);
    end
    return res
end

results = singlet_pairing_scaling(19, :s;Δ=5.0,μ=0.5);
resultd = singlet_pairing_scaling(19, :d;Δ=5.0,μ=0.5);
begin
    fig = plt.figure(figsize=(4, 3))
    plt.plot(collect(3:10), real(results[3:10]), label="swave","o")
    plt.plot(collect(3:10), real(resultd[3:10]), label="dwave","o")
    plt.plot(collect(3:10), real(results[3:10]), color = "C0")
    plt.plot(collect(3:10), real(resultd[3:10]), color = "C1")
    plt.legend(["swave", "dwave"])
    plt.xlabel("distance |i-j|",fontsize = 15)
    plt.ylabel(raw"$\langle \Delta_i^{\dagger}\Delta_j+\Delta_j^{\dagger}\Delta_i\rangle$", fontsize=15)
    plt.show()
    fig.savefig("pairing.pdf", bbox_inches="tight")
end


results_sym = symmetrized_singlet_pairing_scaling(19, :s;Δ=5.0,μ=0.5);
resultd_sym = symmetrized_singlet_pairing_scaling(19, :d;Δ=5.0,μ=0.5);

begin
    fig = plt.figure(figsize=(4, 3))
    plt.plot(collect(3:10), real(results_sym[3:10]), label="swave","o")
    plt.plot(collect(3:10), real(resultd_sym[3:10]), label="dwave","o")
    plt.plot(collect(3:10), real(results_sym[3:10]), color = "C0")
    plt.plot(collect(3:10), real(resultd_sym[3:10]), color = "C1")
    plt.legend(["swave", "dwave"])
    plt.xlabel("distance |i-j|",fontsize = 15)
    plt.ylabel(raw"$\langle (\Delta_i^{\dagger}+\Delta_i)(\Delta_j^{\dagger}+\Delta_j)\rangle$", fontsize=15)
    plt.show()
    fig.savefig("sym_pairing.pdf", bbox_inches="tight")
end

results_anti = antisymmetrized_singlet_pairing_scaling(15, :s);
resultd_anti = antisymmetrized_singlet_pairing_scaling(15, :d);

begin
    plt.plot(collect(1:10), real(results_anti[1:10]), label="swave","o")
    plt.plot(collect(1:10), real(resultd_anti[1:10]), label="dwave","o")
    plt.plot(collect(1:10), real(results_anti[1:10]), color = "C0")
    plt.plot(collect(1:10), real(resultd_anti[1:10]), color = "C1")
    plt.legend(["swave", "dwave"])
    plt.xlabel("distance |i-j|",fontsize = 15)
    plt.show()
end


#system_size_scaling
"""
The following function is used to calculated (Delta_i^dag Delta_j^dag+Delta_j Delta_i) for s-wave and d-wave pairing.
"""
function symmetrized_system_scaling(type::Symbol)
    Ls = collect(9:2:21);
    res = Vector{ComplexF64}(undef, length(Ls));
    o = Op(8);
    Δidag = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3);
    Δjdag = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7);
    obs = (Δidag+Δidag')*(Δjdag+Δjdag');
    i = 1;
    for L in Ls
        if type == :s
            G = BCS_G(L, :s);
        elseif type == :d
            G = BCS_G(L, :d);
        else
            error("type not implemented")
        end
        l = 5;
        Is = [Cartesian2Index([1,1],[L,L],1),
        Cartesian2Index([1,1],[L,L],2),
        Cartesian2Index([2,1],[L,L],1),
        Cartesian2Index([2,1],[L,L],2),
        Cartesian2Index([1,1+l],[L,L],1),
        Cartesian2Index([1,1+l],[L,L],2),
        Cartesian2Index([2,1+l],[L,L],1),
        Cartesian2Index([2,1+l],[L,L],2)];
        ρs = RDM(G, Is);
        res[i] = tr(ρs*obs);
        i += 1;
    end
    return Ls, res
end
Ls, res_s = symmetrized_system_scaling(:s);
Ls, res_d = symmetrized_system_scaling(:d);
begin
    plt.plot(Ls, real(res_s), label="swave")
    plt.plot(Ls, real(res_d), label="dwave")
    # plt.legend(["swave", "dwave"])
    plt.xlabel("system size",fontsize = 15)
    plt.show()
end

"""
Pairing of c_dag_(i,1,up)c_dag_(i,2,down)c_(j,1,up)c_(j,2,down) 
"""
function updown_updown(L::Int64,type::Symbol)
    if type == :s
        G = BCS_G(L, :s;Δ=10.0);
    elseif type == :d
        G = BCS_G(L, :d;Δ=10.0);
    else
        error("type not implemented")
    end
    res = Vector{ComplexF64}(undef, L-1);
    o = Op(8);
    Δidag = ad(o,1)*ad(o,4);
    Δj = a(o,5)*a(o,8);
    obs = Δidag*Δj;
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
        res[i] = tr(ρs*obs);
    end
    return res
end
"""
Pairing of c_dag_(i,1,up)c_dag_(i,2,down)c_(j,1,down)c_(j,2,up) 
"""
function updown_downup(L::Int64,type::Symbol)
    if type == :s
        G = BCS_G(L, :s;Δ=10.0);
    elseif type == :d
        G = BCS_G(L, :d;Δ=10.0);
    else
        error("type not implemented")
    end
    res = Vector{ComplexF64}(undef, L-1);
    o = Op(8);
    Δidag = ad(o,1)*ad(o,4);
    Δj = a(o,6)*a(o,7);
    obs = Δidag*Δj;
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
        res[i] = tr(ρs*obs);
    end
    return res
end

scaling1_s = updown_updown(15, :s);
scaling1_d = updown_updown(15, :d);
scaling2_s = updown_downup(15, :s);
scaling2_d = updown_downup(15, :d);

begin
    plt.plot(real(scaling1_s), label="swave udud")
    plt.plot(real(scaling2_s), label="swave uddu")
    plt.plot(real(scaling1_d), label="dwave udud")
    plt.plot(real(scaling2_d), label="dwave uddu")
    plt.legend()
    plt.xlabel("system size",fontsize = 15)
    plt.show()
end



scaling1_s
scaling2_s
scaling1_d
scaling2_d


o = Op(8);
Δidag = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3);
Δjdag = ad(o,5)*ad(o,8)-ad(o,6)*ad(o,7);
Δidag'*Δjdag≈Δjdag*Δidag'


o = Op(4);
Δidag = ad(o,1)*ad(o,4)-ad(o,2)*ad(o,3)
Matrix((Δidag+Δidag'))