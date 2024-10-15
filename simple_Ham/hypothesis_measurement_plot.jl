using LinearAlgebra
using SparseArrays
using Fermionic
using ProgressMeter
using Statistics
using ArgParse
using JSON3
using JLD2
using PythonCall
using Distributions
plt = pyimport("matplotlib.pyplot")
function gaussian_dist(x::Float64,mean::Float64,std::Float64)
    return 1/(std*sqrt(2*pi))*exp(-0.5*((x-mean)/std)^2)
end
# Here is how we use t-distribution to calculate confidence
filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/processed_measurement300_L40_type1_dis5_Delta0.1_mu0.5_samples10000000.json"
data = JSON3.read(filename)
signs = [1.0,-1.0,-1.0,1.0];
dwave = []
for bt in 1:10000
    tmp = 0.0
    for left in 1:4
        for right in 1:4
            tmp+= signs[left]*signs[right]*data["$(left)_$(right)"][bt]
        end
    end
    append!(dwave,tmp)
end
begin
    count, bins, _ = plt.hist(dwave,bins = 40,edgecolor="black",alpha=0.5)
    count = pyconvert(Vector,count)
    bins = pyconvert(Vector,bins)
    plt.show()
end


# Here we calculate confidence of all measurements
measurements = [50,100,150,200,300,400,500,600,1000]
confidence_01 = []
confidence_02 = []
threshold = 0.01
for m in measurements
    t_dist = TDist(m)
    filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/processed_measurement$(m)_L40_type1_dis5_Delta0.1_mu0.5_samples10000000.json"
    data = JSON3.read(filename)
    signs = [1.0,-1.0,-1.0,1.0];
    dwave = []
    for bt in 1:10000
        tmp = 0.0
        for left in 1:4
            for right in 1:4
                tmp+= signs[left]*signs[right]*data["$(left)_$(right)"][bt]
            end
        end
        append!(dwave,tmp)
    end
    shifted_scaled_cdf = x -> cdf(t_dist, (x - mean(dwave)) / std(dwave))
    confidence_ = 1-shifted_scaled_cdf(threshold)
    append!(confidence_01,confidence_)
end
for m in measurements
    t_dist = TDist(m)
    filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/processed_measurement$(m)_L40_type1_dis5_Delta0.2_mu0.5_samples10000000.json"
    data = JSON3.read(filename)
    signs = [1.0,-1.0,-1.0,1.0];
    dwave = []
    for bt in 1:10000
        tmp = 0.0
        for left in 1:4
            for right in 1:4
                tmp+= signs[left]*signs[right]*data["$(left)_$(right)"][bt]
            end
        end
        append!(dwave,tmp)
    end
    shifted_scaled_cdf = x -> cdf(t_dist, (x - mean(dwave)) / std(dwave))
    confidence_ = 1-shifted_scaled_cdf(threshold)
    append!(confidence_02,confidence_)
end

ScalarFormatter = pyimport("matplotlib.ticker").ScalarFormatter
begin
    f = plt.figure(figsize = (3.5,3.3))
    plt.plot(measurements,confidence_01,"o--",color="#1B78B2")
    plt.plot(measurements,confidence_02,"o--",color="#BEDCED")
    plt.legend([raw"$Δ_{\text{BCS}}=0.1$",raw"$Δ_{\text{BCS}}=0.2$"],fontsize = 12,loc=4)
    plt.ylabel("Success prob.",fontsize = 12)
    plt.xlabel("Number of measurements",fontsize = 12)
    # Set y-axis to display two decimal places
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=true))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    plt.show()
    f.savefig("./test_confidence_2.svg",bbox_inches="tight")
end

# Plot two histograms
filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/processed_measurement100_L40_type1_dis5_Delta0.1_mu0.5_samples10000000.json"
data = JSON3.read(filename)
signs = [1.0,-1.0,-1.0,1.0];
dwave_100 = []
for bt in 1:10000
    tmp = 0.0
    for left in 1:4
        for right in 1:4
            tmp+= signs[left]*signs[right]*data["$(left)_$(right)"][bt]
        end
    end
    append!(dwave_100,tmp)
end
filename = "/Users/hyhu/Git_Code/Fermionic.jl/simple_Ham/simple_direct_measurement_data/hypothesis_data/processed_measurement300_L40_type1_dis5_Delta0.1_mu0.5_samples10000000.json"
data = JSON3.read(filename)
signs = [1.0,-1.0,-1.0,1.0];
dwave_300 = []
for bt in 1:10000
    tmp = 0.0
    for left in 1:4
        for right in 1:4
            tmp+= signs[left]*signs[right]*data["$(left)_$(right)"][bt]
        end
    end
    append!(dwave_300,tmp)
end
begin
    f = plt.figure(figsize = (3.5,3.3))
    counts300, bins300, _ = plt.hist(dwave_300,bins = 35,edgecolor="black",alpha=0.4,label="300 Meas.")
    counts300 = pyconvert(Vector,counts300)
    bins300 = pyconvert(Vector,bins300)
    counts100, bins100, _ = plt.hist(dwave_100,bins = bins300,edgecolor="black",alpha=0.4,label="100 Meas.")
    counts100 = pyconvert(Vector,counts100)
    bins100 = pyconvert(Vector,bins100)
    # Add a vertical line at x=0.01
    plt.axvline(x=0.06, color="gray", linestyle="--", linewidth=1.5)
    # Add text near the vertical line
    plt.text(0.1, 750, raw"ϵ=0.05", fontsize=10.5, color="black")  # Adjust text position as needed
    plt.text(-0.4, 600, raw"No", fontsize=11.5, color="black")
    plt.text(-0.4, 530, raw"d-wave", fontsize=11.5, color="black")  # Adjust text position as needed
    plt.text(1.05, 600, raw"Exist", fontsize=11.5, color="black")
    plt.text(1.05, 530, raw"d-wave", fontsize=11.5, color="black")  # Adjust text position as needed
    plt.ylabel("Histogram",fontsize = 12)
    plt.xlabel("Estimated d-wave pairing",fontsize = 12)
    plt.legend(handlelength=1, handleheight=1, handletextpad=0.5, borderpad=0.3,frameon=false)
    plt.show()
    f.savefig("./hypothesis_test_histogram_2.svg",bbox_inches="tight")
end
mean(dwave_300)
0.05/mean(dwave_300)

## Here we show the t-distribution converge to Gaussian distribution when sample is large than 300
begin
    y = [gaussian_dist(x,mean(dwave),std(dwave)) for x in bins[1:end-1]]
    plt.plot(bins[1:end-1],count/(sum(count)*(bins[2]-bins[1])),"o")
    plt.plot(bins[1:end-1],y)
    plt.show()
end
mean(dwave)
std(dwave)
t_dist = TDist(300-1)
shifted_scaled_pdf = x -> pdf(t_dist, (x - mean(dwave)) / std(dwave)) / std(dwave)
y2 = [shifted_scaled_pdf(x) for x in bins[1:end-1]]
sum(y2)*(bins[2]-bins[1])
begin
    y = [gaussian_dist(x,mean(dwave),std(dwave)) for x in bins[1:end-1]]
    y2 = [shifted_scaled_pdf(x) for x in bins[1:end-1]]
    plt.plot(bins[1:end-1],count/(sum(count)*(bins[2]-bins[1])),"o")
    plt.plot(bins[1:end-1],y)
    plt.plot(bins[1:end-1],y2,"--")
    plt.show()
end

shifted_scaled_cdf = x -> cdf(t_dist, (x - mean(dwave)) / std(dwave))
shifted_scaled_cdf(0.0)
println("confidence: ",1-shifted_scaled_cdf(0.01))
# test it matches integral of gaussian dist and t-dist
y_gaussain_int = [gaussian_dist(x,mean(dwave),std(dwave)) for x in -3:0.01:0.0]
y_t_int = [shifted_scaled_pdf(x) for x in -3:0.01:0.0]
sum(y_gaussain_int)*0.01
sum(y_t_int)*0.01