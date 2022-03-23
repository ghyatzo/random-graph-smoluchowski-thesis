using JLD2
using CairoMakie, GraphMakie, NetworkLayout
CairoMakie.activate!(type = "pdf")

includet("../utils.jl")
includet("../RandomGraphs/random_graph_multicomp.jl")
## ------- Error Analysis Study

# We shall focus only on weak components.
# We shall pick a single time t to perform the calculations (at first, then we can do it for more time)
# We then generate che corresponding random graphs and build the solution.
# Then we compute the cumulative absolute error accross all entries and plot it againts
#	an increasing number of nodes.

# near to tc, but not too much.
t=0.39

## Code used to generate the error
# η = [0.5, 0.5]
# kmax = 100
# Nspan = [floor(Int, 2^(n/4)) for n=28:80]
# N_err = zeros(kmax-1, length(Nspan))
# for (ni, nn) = enumerate(Nspan)
# 	@info "running n: $nn"
# 	test_sol = avg_weak_cc_dist(t, kernel, nn, η; kmax=kmax-1) ./ nn
# 	N_err[:, ni] = test_sol
# end
# avg_err_N = (sum(abs.(N_err .- true_sol); dims=1) ./ 99)'
# cum_err_N = sum(abs.(N_err .- true_sol); dims=1)'

Nspan = load("../../data/cumulative_abs_error_N=1e6.jld2")["Nspan"]
cum_err_N = [load("../../data/cumulative_abs_error_N=1e6.jld2")["err"]...]
avg_err_N = [load("../../data/avg_abs_error_N=1e6.jld2")["err"]...]

ax_err_kw = (
	xscale=log10,
 	yminorticksvisible = true,
 	yminorgridvisible = true,
 	xminorticksvisible = true,
 	xminorgridvisible = true,
 	yminorticks = IntervalsBetween(5),
 	xminorticks = IntervalsBetween(10),
 	xlabel="N",aspect=1)
begin
F_Err = Figure(resolution=(1000, 500))
Ax_cum = Axis(F_Err[1, 1]; ylabel="err", ax_err_kw...)
Ax_avg = Axis(F_Err[1, 2]; yscale=log10, ax_err_kw..., ylabel="log err")

scatter!(Ax_cum, Nspan, cum_err_N; marker=:diamond, color=:darkorange)
scatter!(Ax_avg, Nspan, cum_err_N; marker=:xcross, color=:darkorange)
lines!(Ax_avg, Nspan, 0.214 .*Nspan.^(-0.5); color=:darkred, label=L"\frac{1}{5\sqrt{N}}")
axislegend(Ax_avg)

save("../../imgs/error_anal/cumulative_error_normallog_loglog.pdf", F_Err)
end

# regression
y = log10.(cum_err_N)
x = log10.(Nspan)
sx = sum(x)
sy = sum(y)
sxy = x'*y
sx2 = sum(x.^2)
n = length(Nspan)
m = (n*sxy - sx*sy)/(n*sx2 - sx^2)
b = 10^(sy/n - m*(sx/n))


## ------- Time Complexity (time/num of nodes)
tcd = load("../../data/time_complexity_data_Nmax=7e7-NEW.jld2")
tcd_tri = load("../../data/time_complexity_data_d=3_Nmax=7e7-NEWNEW.jld2")

Ntspan = tcd["Nspan"]
tcd_sec = tcd["time"]
tcd_gc = tcd["gc"]

tcd_tri_sec = tcd_tri["time"]
tcd_tri_gc = tcd_tri["gc"]

ax_time_kw = (
	xscale=log10,
 	xminorticksvisible = true,
 	xminorgridvisible = true,
 	xminorticks = IntervalsBetween(10),
 	yscale=log10,
 	yminorticksvisible = true,
 	yminorgridvisible = true,
 	yminorticks = IntervalsBetween(10),
 	xlabel="N", ylabel=L"runtime (s)", aspect=1
 )

F_timeN = Figure()
AxtimeN = Axis(F_timeN[1,1]; limits = (60, 1e8, nothing, nothing), ax_time_kw...)
AxgcN = Axis(F_timeN[1, 1]; xscale = log10, yscale=log10, yaxisposition=:right, aspect=1, ylabel="garbage collection (s)", limits = (60, 1e8, nothing, nothing))

lint = lines!(AxtimeN, Ntspan, tcd_sec; color=:darkorange, label="timed=2")
scat = scatter!(AxtimeN, Ntspan, tcd_sec; marker=:diamond, color=:darkorange, markersize=20, label="time d=2")
lint_tri = lines!(AxtimeN, Ntspan, tcd_tri_sec; color=:darkblue, label="time d=3")
scat_tri = scatter!(AxtimeN, Ntspan, tcd_tri_sec; marker=:diamond, color=:darkblue, markersize=20, label="time d=3")

lingc = lines!(AxgcN, Ntspan, nanmap(tcd_gc); color=:darkred, label="gc d=2")
scagc = scatter!(AxgcN, Ntspan, nanmap(tcd_gc); marker=:rect, color=:darkred, markersize=10, label="gc d=2")
lingc_tri = lines!(AxgcN, Ntspan, nanmap(tcd_tri_gc); color=:darkcyan, label="gc d=3")
scagc_tri = scatter!(AxgcN, Ntspan, nanmap(tcd_tri_gc); marker=:rect, color=:darkcyan, markersize=10, label="gc d=3")

axislegend(AxgcN,
    [[lint, scat], [lingc, scagc], [lint_tri, scat_tri], [lingc_tri, scagc_tri]],["time d=2", "gc d=2", "time d=3", "gc d=3"],
    tellheight = false,
    tellwidth = false,
    margin = (10, 100, 10, 10),
    halign=:left, valign=:top, framevisible=false)

save("../../imgs/time_comp/time_complexity_and_gc_times_new.pdf", F_timeN)

## ----- Time Complexity 2, Number of Iterations
iters2d_data = load("../../data/iterations_data_Nmax=7e7.jld2")
iters3d_data = load("../../data/iterations_data_d=3_Nmax=7e7.jld2")
Ntspan = iters2d_data["Nspan"]
iters2d = iters2d_data["iterations"]
iters3d = iters3d_data["iterations"]

ax_iter_kw = (
    xscale=log10,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
    yscale=log10,
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xlabel="N", ylabel="Iterations", aspect=1
 )

F_iterN = Figure()
AxiterN = Axis(F_iterN[1,1]; limits = (60, 1e8, nothing, nothing), ax_iter_kw...)

lint = lines!(AxiterN, Ntspan, iters2d; color=:darkorange, label="d=2")
scat = scatter!(AxiterN, Ntspan, iters2d; marker=:diamond, color=:darkorange, markersize=20, label="d=2")

lingc = lines!(AxiterN, Ntspan, iters3d; color=:darkblue, label="d=3")
scagc = scatter!(AxiterN, Ntspan, iters3d; marker=:diamond, color=:darkblue, markersize=20, label="d=3")
axislegend(AxiterN,
    [[lint, scat], [lingc, scagc]],
    ["d=2", "d=3"], framevisible = false,
    tellheight = false,
    tellwidth = false,
    margin = (10, 100, 10, 10),
    halign=:left, valign=:top)
save("../../imgs/time_comp/iterations_over_N.pdf", F_iterN)


## ----- Time complexity 3, O constant (time/(d+1)N)

Otcd_sec_2 = @. iters2d / ((2+0.39/0.526)*Ntspan)
Otcd_sec_3 = @. iters3d / ((3+0.39/0.526)*Ntspan)

ax_Oiter_kw = (
    xscale=log10,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(10),
    yminorticksvisible = true,
    yminorgridvisible = true,
    yminorticks = IntervalsBetween(10),
    xlabel="N", ylabel=L"Iter/(d+\frac{t}{t_c})N", aspect=1
 )

F_timeON = Figure()
AxtimeON = Axis(F_timeON[1,1];limits = (nothing, nothing, 0.5, 1.5), ax_Oiter_kw...)
lintO2 = lines!(AxtimeON, Ntspan, Otcd_sec_2; color=:darkorange, label="d=2")
scatO2 = scatter!(AxtimeON, Ntspan, Otcd_sec_2; marker=:rect, color=:darkorange, markersize=12, label="d=2")
lintO3 = lines!(AxtimeON, Ntspan, Otcd_sec_3; color=:darkblue, label="d=3")
scatO3 = scatter!(AxtimeON, Ntspan, Otcd_sec_3; marker=:xcross, color=:darkblue, markersize=14, label="d=3")

axislegend(AxtimeON,
    [[lintO2, scatO2], [lintO3, scatO3]],["d=2", "d=3"],
    tellheight = false,
    tellwidth = false,
    margin = (10, 100, 10, 10),
    halign=:left, valign=:top,
    framevisible = false)

save("../../imgs/time_comp/iterations_O_constant_new.pdf", F_timeON)