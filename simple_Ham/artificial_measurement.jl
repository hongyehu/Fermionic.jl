using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using PythonCall
plt = pyimport("matplotlib.pyplot")


o = Op(8)
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
obs = two_paddle_obs(left=1,right=1)
u, v = eigen(Matrix(obs))
tmp = v'*(Matrix(obs))*v

begin
    plt.imshow(real(tmp))
    plt.show()
end

ρmeasure!(ρ)