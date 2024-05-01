using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")

L = 11
G = BCS_G(L, :s)
ρs = RDM(G, [1 1; 1 2; 3 1]);
G
function get_correlator_matrix(L::Int64,type::Symbol)
    if type == :s
        G = getBCSaij(L, :s);
    elseif type == :d
        G = getBCSaij(L, :d);
    else
        error("type not implemented")
    end
    i_site = Cartesian2Index([Int((L+1)/2),Int((L+1)/2)],[L,L],1);
    k_site = Cartesian2Index([Int((L+1)/2),Int((L+1)/2+1)],[L,L],1);
    l_site = Cartesian2Index([Int((L+1)/2),Int((L+1)/2)],[L,L],2);#i_site<l_site<k_site
    res = Matrix{ComplexF64}(undef, L, L);
    for jx in 1:L
      for jy in 1:L
        if jy<(L+1)/2
			#j_site<i_site<l_site<k_site
        	j_site = Cartesian2Index([jx,jy],[L,L],2);
			Is = [j_site,i_site,l_site,k_site];
			ρs = RDM(G, Is);
			o = Op(4);
			obs = ad(o,2)*ad(o,1)*a(o,4)*a(o,3);
			res[jx,jy] = tr(ρs*obs);
        elseif jy>(L+1)/2
			#i_site<l_site<k_site<j_site
        	j_site = Cartesian2Index([jx,jy],[L,L],2);
			Is = [i_site,l_site,k_site,j_site];
			ρs = RDM(G, Is);
			o = Op(4);
			obs = ad(o,1)*ad(o,4)*a(o,3)*a(o,2);
			res[jx,jy] = tr(ρs*obs);
		elseif jy==(L+1)/2&&jx<(L+1)/2
			#j_site<i_site<l_site<k_site
        	j_site = Cartesian2Index([jx,jy],[L,L],2);
			Is = [j_site,i_site,l_site,k_site];
			ρs = RDM(G, Is);
			o = Op(4);
			obs = ad(o,2)*ad(o,1)*a(o,4)*a(o,3);
			res[jx,jy] = tr(ρs*obs);
		elseif jy==(L+1)/2&&jx>(L+1)/2
			#i_site<l_site<k_site<j_site
        	j_site = Cartesian2Index([jx,jy],[L,L],2);
			Is = [i_site,l_site,k_site,j_site];
			ρs = RDM(G, Is);
			o = Op(4);
			obs = ad(o,1)*ad(o,4)*a(o,3)*a(o,2);
			res[jx,jy] = tr(ρs*obs);
		elseif jy==(L+1)/2&&jx==(L+1)/2
			#i_site<l_site==j_site<k_site
			j_site = Cartesian2Index([jx,jy],[L,L],2);
			Is = [i_site,l_site,k_site];
			ρs = RDM(G, Is);
			o = Op(3);
			obs = ad(o,1)*ad(o,2)*a(o,3)*a(o,2);
			res[jx,jy] = tr(ρs*obs);
        end
      end
    end
    return res
end

tests = get_correlator_matrix(11,:s);
testd = get_correlator_matrix(11,:d);
begin
	vmin = -0.045
	vmax = 0.03
	fig, axs = plt.subplots(1, 2, figsize=(10, 5))
	img1 = axs[0].imshow(real(tests),cmap="Oranges_r")
	axs[0].set_title("s-wave")
	img2 = axs[1].imshow(real(testd),cmap="Oranges_r")
	axs[1].set_title("d-wave")
	img1.set_clim(vmin, vmax)
	img2.set_clim(vmin, vmax)

	# Show color bar
	fig.colorbar(img1, ax=axs[0])
	fig.colorbar(img2, ax=axs[1])
	plt.show()
end

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

ξfun(0.3, 0.1)
5.0*(cos(0.3)-cos(0.1))
getak(2*π/11*-5, 2*π/11*3, :d)
a_coef(2*π/11*0, 2*π/11*8, :d)
rho1 = getBCSaij(3,:s)
rho2 = BCS_G(3, :s)
begin
	fig, axs = plt.subplots(1, 2, figsize=(10, 5))
	axs[0].imshow(real(rho1))
	axs[0].set_title("Minh")
	axs[1].imshow(real(rho2))
	axs[1].set_title("HYH")
	plt.show()
end
real(rho1)[1:10,1:10]
real(rho2)[1:10,1:10]


tmp = 2 * π / L * collect(-(L-1)/2:1:(L-1)/2)
tmp
Cartesian2Index([1,2],[2,2],1)
collect(-(L-1)/2:1:(L-1)/2)