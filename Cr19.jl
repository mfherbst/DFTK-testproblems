# Chromium bulk with a defect
#
using DFTK
using MKL
setup_threading(n_blas=2)

Cr = ElementPsp(:Cr, psp=load_psp("hgh/pbe/cr-q14.hgh"))
atoms = fill(Cr, 19)
magnetic_moments = fill(5.0, 19)
lattice   = load_lattice("Cr19.cif")
positions = load_positions("Cr19.cif")

model = model_PBE(lattice, atoms, positions;
                  temperature=1e-2,
                  smearing=Smearing.Gaussian(),
                  magnetic_moments)


full = true
if full
    basis = PlaneWaveBasis(model; kgrid=[9, 9, 1], Ecut=35)
else
    basis = PlaneWaveBasis(model; kgrid=[3, 3, 1], Ecut=20)
end

ρ0 = guess_density(basis, magnetic_moments)
DFTK.reset_timer!(DFTK.timer)
scfres = self_consistent_field(basis; ρ=ρ0, is_converged=DFTK.ScfConvergenceDensity(1e-6))
println(DFTK.timer)
