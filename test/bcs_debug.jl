using Fermionic
using PythonCall
using LinearAlgebra
plt = pyimport("matplotlib.pyplot")

# Minh's code
using LinearAlgebra, SparseArrays
o=Op(3);
basis(o)
ad(o,1)
# G = BCS_G(3, :d)
index2coordinates(3,4)
G1 = getBCSaij(3,:s);
G2 = real(BCS_G(3,:s));
begin
	vmin = -0.045
	vmax = 0.03
	fig, axs = plt.subplots(1, 2, figsize=(10, 5))
	img1 = axs[0].imshow(real(G1),cmap="Oranges_r")
	axs[0].set_title("Minh")
	img2 = axs[1].imshow(real(G2),cmap="Oranges_r")
	axs[1].set_title("HYH")
	img1.set_clim(vmin, vmax)
	img2.set_clim(vmin, vmax)

	# Show color bar
	fig.colorbar(img1, ax=axs[0])
	fig.colorbar(img2, ax=axs[1])
	plt.show()
end


