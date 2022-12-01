# Aluminium - Silica mixed slab from doi:10.1016/j.jcp.2022.111127
#
using DFTK
using MKL
setup_threading(n_blas=2)


full = true
if full
    file = joinpath(@__DIR__, "AlSiO2H_20.cif")
else
    file = joinpath(@__DIR__, "AlSiO2H_10.cif")
end

lattice = load_lattice(file)
Al = ElementPsp(:Al, psp=load_psp("hgh/pbe/Al-q3.hgh"))
H  = ElementPsp(:H,  psp=load_psp("hgh/pbe/h-q1.hgh"))
O  = ElementPsp(:O,  psp=load_psp("hgh/pbe/o-q6.hgh"))
Si = ElementPsp(:Si, psp=load_psp("hgh/pbe/Si-q4.hgh"))
atoms = map(atomic_symbol.(load_atoms(file))) do sym
    Dict(:Al => Al, :H => H, :O => O, :Si => Si)[sym]
end
positions = load_positions(file)

model  = model_PBE(lattice, atoms, positions;
                   temperature=1e-3, smearing=Smearing.Gaussian())
basis  = PlaneWaveBasis(model; Ecut=30, kgrid=[2, 2, 1])

DFTK.reset_timer!(DFTK.timer)
scfres = self_consistent_field(basis; is_converged=DFTK.ScfConvergenceDensity(1e-6))
println(DFTK.timer)
