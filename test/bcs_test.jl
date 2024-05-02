using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")


function get_correlator_matrix(L::Int64,type::Symbol)
    if type == :s
        G = BCS_G(L, :s);
    elseif type == :d
        G = BCS_G(L, :d);
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
    axs[1].set_xlabel("X-axis label")  # Set x-axis label
    axs[1].set_ylabel("Y-axis label") 
	img1.set_clim(vmin, vmax)
	img2.set_clim(vmin, vmax)

	# Show color bar
	fig.colorbar(img1, ax=axs[0])
	fig.colorbar(img2, ax=axs[1])
	plt.show()
end
