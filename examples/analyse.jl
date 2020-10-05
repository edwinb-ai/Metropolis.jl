using Metropolis
using CSV
using Bootstrap
using Statistics

filedata = "npt_0.850_8.100.csv"

to_analyse = joinpath("results", filedata)
data = CSV.File(to_analyse; header = false)

n_boot = 200000

cil = 0.95
bs_mean = bootstrap(mean, data.Column2, BasicSampling(n_boot))
bci1 = confint(bs_mean, BasicConfInt(cil))

display(bs_mean)
display(bci1)
