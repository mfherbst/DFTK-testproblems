# Iron-Wolfram mixed slab from doi:10.1021/acs.jctc.1c00630
#
using DFTK
using MKL
setup_threading(n_blas=2)

function attach_pseudos(atoms::AbstractArray; pseudomap...)
    pseudomap = Dict(pseudomap)
    map(atoms) do element
        pspfile = get(pseudomap, element.symbol, nothing)
        ElementPsp(element.symbol; psp=load_psp(pspfile))
    end
end

function build_magnetic_moments(atoms::AbstractArray; magmoms...)
    magmoms = Dict(magmoms)
    map(atoms) do element
        magmoms[element.symbol]
    end
end

lattice   = load_lattice(joinpath(@__DIR__, "WFe.cif"))
positions = load_positions(joinpath(@__DIR__, "WFe.cif"))
atoms     = load_atoms(joinpath(@__DIR__, "WFe.cif"))
atoms     = attach_pseudos(atoms, Fe="hgh/pbe/fe-q16.hgh", W="hgh/pbe/w-q14.hgh")
magnetic_moments = build_magnetic_moments(atoms; Fe=3.0, W=2.0)

model = model_PBE(lattice, atoms, positions;
                  smearing=Smearing.Gaussian(), temperature=0.01, magnetic_moments)


full = true
if full
    basis = PlaneWaveBasis(model; Ecut=25, kgrid=(7, 7, 1), supersampling=1.5)
else
    basis = PlaneWaveBasis(model; Ecut=15, kgrid=(3, 3, 1))
end

scfargs = (; ρ=guess_density(basis, magnetic_moments),
             is_converged=DFTK.ScfConvergenceDensity(1e-6),
             mixing=KerkerMixing(),
          )
DFTK.reset_timer!(DFTK.timer)
scfres = self_consistent_field(basis; scfargs...);
println(DFTK.timer)
