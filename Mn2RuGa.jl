# Tricky compensation between spin on nearest-neighbours Mn atoms
# https://arxiv.org/pdf/1506.03735.pdf
#
using DFTK
using LinearAlgebra
using MKL
setup_threading(n_blas=2)

lattice = 5.66917837387731 * [0 1 1;
                              1 0 1;
                              1 1 0]
Mn = ElementPsp(:Mn, psp=load_psp("hgh/pbe/mn-q15.hgh"))
Ru = ElementPsp(:Ru, psp=load_psp("hgh/pbe/ru-q16.hgh"))
Ga = ElementPsp(:Ga, psp=load_psp("hgh/pbe/ga-q3.hgh"))
atoms     = [Ru, Mn, Mn, Ga]
positions = [[0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [0.75, 0.75, 0.75], [0.25, 0.25, 0.25]]

magnetic_moments = [
    5.0,  # Ru
    5.0,  # Mn
    5.0,  # Mn
    0.0,  # Ga
]

model = model_PBE(lattice, atoms, positions;
                  temperature=1e-2, smearing=Smearing.Gaussian(), magnetic_moments)
full = true
if full
    basis = PlaneWaveBasis(model; kgrid=[13, 13, 13], Ecut=35)
else
    basis = PlaneWaveBasis(model; kgrid=[8, 8, 8], Ecut=25)
end

ρ0 = guess_density(basis, magnetic_moments)
DFTK.reset_timer!(DFTK.timer)
scfres = self_consistent_field(basis; ρ=ρ0, is_converged=DFTK.ScfConvergenceDensity(1e-6))
println(DFTK.timer)
