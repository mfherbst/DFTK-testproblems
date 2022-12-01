# half-metallic ferromagnetic character
# https://sci-hub.se/10.1016/j.jmmm.2014.10.094
#
using DFTK
using LinearAlgebra
using MKL
setup_threading(n_blas=2)

lattice = 5.50855165328412 * [0 1 1;
                              1 0 1;
                              1 1 0]

Fe = ElementPsp(:Fe, psp=load_psp("hgh/pbe/fe-q16.hgh"))
Mn = ElementPsp(:Mn, psp=load_psp("hgh/pbe/mn-q15.hgh"))
Al = ElementPsp(:Al, psp=load_psp("hgh/pbe/al-q3.hgh"))
atoms     = [Mn, Al, Fe, Fe]
positions = [[0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]]

magnetic_moments = [
    5.0,  # Mn
    0.0,  # Al
    5.0,  # Fe
    5.0,  # Fe
]

model = model_PBE(lattice, atoms, positions;
                  temperature=1e-2, smearing=Smearing.Gaussian(), magnetic_moments)
full = true
if full
    basis = PlaneWaveBasis(model; kgrid=[13, 13, 13], Ecut=45)
    # basis = PlaneWaveBasis(model; kgrid=[13, 13, 13], Ecut=35)
else
    basis = PlaneWaveBasis(model; kgrid=[8, 8, 8], Ecut=25)
end

ρ0 = guess_density(basis, magnetic_moments)
DFTK.reset_timer!(DFTK.timer)
scfres = self_consistent_field(basis; ρ=ρ0, is_converged=DFTK.ScfConvergenceDensity(1e-6))
println(DFTK.timer)
