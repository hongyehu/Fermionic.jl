using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")

function singlet_pairing_scaling(L::Int64, type::Symbol)
    if type == :s
        G = BCS_G(L, :s);
    elseif type == :d
        G = BCS_G(L, :d);
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
function singlet_pairing_scaling_extra(L::Int64, type::Symbol)
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
    Δjdag = a(o,5)*a(o,8)-a(o,6)*a(o,7);
    obs = Δidag*Δjdag+(Δidag'*Δjdag');
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

results = singlet_pairing_scaling(15, :s);
resultd = singlet_pairing_scaling(15, :d);
begin
    plt.plot(collect(1:10), real(results[1:10]), label="swave","o")
    plt.plot(collect(1:10), real(resultd[1:10]), label="dwave","o")
    plt.plot(collect(1:10), real(results[1:10]), color = "C0")
    plt.plot(collect(1:10), real(resultd[1:10]), color = "C1")
    plt.legend(["swave", "dwave"])
    plt.xlabel("distance |i-j|",fontsize = 15)
    # plt.ylabel(r"Re(\$\langle \Delta_i^{\dagger}\Delta_j\$)",fontsize = 15)
    plt.show()
end

results_extra = singlet_pairing_scaling_extra(15, :s);
resultd_extra = singlet_pairing_scaling_extra(15, :d);

begin
    plt.plot(collect(1:10), real(results_extra[1:10]), label="swave","o")
    plt.plot(collect(1:10), real(resultd_extra[1:10]), label="dwave","o")
    plt.plot(collect(1:10), real(results_extra[1:10]), color = "C0")
    plt.plot(collect(1:10), real(resultd_extra[1:10]), color = "C1")
    plt.legend(["swave", "dwave"])
    plt.xlabel("distance |i-j|",fontsize = 15)
    # plt.ylabel(r"Re(\$\langle \Delta_i^{\dagger}\Delta_j\$)",fontsize = 15)
    plt.show()
end