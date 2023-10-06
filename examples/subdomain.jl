# First we load Ferrite, and some other packages we need.
using FerriteGmsh, Ferrite, FerriteCon2023
using SparseArrays, BlockArrays
using CairoMakie, Colors

grid = togrid("logo_refined.geo");

# Trial and test functions
dim = 2
order_u = 2 # quadratic interpolation
order_p = 1 # linear interpolation

qr = QuadratureRule{RefTriangle}(2) # 3 quadrature points in total
ip_u = Lagrange{RefTriangle, order_u}()^dim
cv_u = CellValues(qr, ip_u);
ip_p = Lagrange{RefTriangle, order_p}()
cv_p = CellValues(qr, ip_p)

cv_up = (cv_u, cv_p)

# ### Degrees of freedom
dh = DofHandler(grid)

# subdomains can only consist of a single element type
sdh_up = SubDofHandler(dh, getcellset(grid, "green"))
add!(sdh_up, :u, ip_u)
add!(sdh_up, :p, ip_p)

sdh_u = SubDofHandler(dh, setdiff!(Set(1:getncells(grid)), getcellset(grid, "green")))
add!(sdh_u, :u, ip_u)

close!(dh);

# ### Boundary conditions
ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), (x, t) -> 0.0, 2))
add!(ch, Dirichlet(:u, getfaceset(grid, "left"), (x, t) -> 0.0, 1))
add!(ch, Dirichlet(:u, getfaceset(grid, "top"), (x, t) -> 0.1, 2))
close!(ch);

########################
# material definitions #
########################
E = 200e3 # Young's modulus [MPa]
ν = 0.3 # Poisson's ratio [-]
compressible = Elasticity(E/2(1+ν), E/3(1-2ν));

ν = 0.5 # Poisson's ratio [-]
incompressible = Elasticity(E/2(1+ν), E/3(1-2ν));

materials = (incompressible, compressible)

#################################
# Pre-allocate element matrices #
#################################
ke_u = zeros(ndofs_per_cell(sdh_u), ndofs_per_cell(sdh_u))
re_u = zeros(ndofs_per_cell(sdh_u))

n = ndofs_per_cell(sdh_up)
blocks = [length(dof_range(sdh_up, :u)), length(dof_range(sdh_up, :p))]
ke_up = PseudoBlockArray(zeros(n, n), blocks, blocks)
re_up = PseudoBlockArray(zeros(n), blocks)

buffers = ((cv_up, ke_up, re_up), (cv_u, ke_u, re_u))

###################
# Global assembly #
###################
function assemble_global!(assembler, dh, buffers, materials)
    for (sdh, _buffers, material) in zip(dh.subdofhandlers, buffers, materials)
        assemble_subdofhandler!(assembler, sdh, _buffers, material)
    end
end

function assemble_subdofhandler!(assembler, sdh, buffers, material)
    cv, ke = buffers
    ## Loop over all cells
    for cell in CellIterator(sdh)
        ## Compute element contribution
        assemble_cell!(ke, cv, material, getcoordinates(cell))
        ## Assemble ke and fe into K and f
        assemble!(assembler, celldofs(cell), ke)
    end
end

##########################
# Solution of the system #
##########################
K = create_sparsity_pattern(dh)
f = zeros(ndofs(dh))
assembler = start_assemble(K, f)
assemble_global!(assembler, dh, buffers, materials);

apply!(K, f, ch)
u = K \ f

##################
# postprocessing #
##################
function compute_pressure_vonmises(cellvalues::CellValues, material, xe, ue)
    reinit!(cellvalues, xe)

    n_basefuncs = getnbasefunctions(cellvalues)

    p = zero(eltype(ue))
    σᵥₘ = 0.0
    area = 0.0
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        area += dΩ
        ## For each integration point, compute strain, stress and material stiffness
        ε = function_symmetric_gradient(cellvalues, q_point, ue)
        ε_3D = SymmetricTensor{2,3}((i,j)->i<3 && j<3 ? ε[i,j] : zero(eltype(ε)))
        σ, = material_routine(material, ε_3D)
        σᵥₘ += sqrt(3/2 * dev(σ) ⊡ dev(σ)) * dΩ
        p += -tr(σ)/3 * dΩ
    end
    return p/area, σᵥₘ/area
end

function compute_mixed_von_mises(cv_u, cv_p, material, xe, ue, pe)
    reinit!(cv_u, xe)
    reinit!(cv_p, xe)

    (; G, K) = material

    σᵥₘ = 0.0
    area = 0.0
    for q_point in 1:getnquadpoints(cv_u)
        dΩ = getdetJdV(cv_u, q_point)
        area += dΩ
        ## For each integration point, compute strain, stress and material stiffness
        p = function_value(cv_p, q_point, pe)
        ε = function_symmetric_gradient(cv_u, q_point, ue)
        ε_3D = SymmetricTensor{2,3}((i,j)->i<3 && j<3 ? ε[i,j] : zero(eltype(ε)))
        σ = 2G * dev(ε_3D) + p * one(ε_3D)
        σᵥₘ += sqrt(3/2 * dev(σ) ⊡ dev(σ)) * dΩ
    end
    return σᵥₘ/area
end

pressures = Vector{Float64}(undef, getncells(grid))
fill!(pressures, NaN)
von_mises_stresses = Vector{Float64}(undef, getncells(grid))
fill!(von_mises_stresses, NaN)

dofrange_u = dof_range(sdh_up, :u)
dofrange_p = dof_range(sdh_up, :p)
# ### Exporting to VTK
for cell in CellIterator(sdh_u)
    @views ue = u[cell.dofs]
    p, σᵥₘ = compute_pressure_vonmises(cv_u, compressible, cell.coords, ue)
    pressures[cellid(cell)] = p 
    von_mises_stresses[cellid(cell)] = σᵥₘ
end
for cell in CellIterator(sdh_up)
    @views ue = u[cell.dofs][dofrange_u]
    @views pe = u[cell.dofs][dofrange_p]
    σᵥₘ = compute_mixed_von_mises(cv_u, cv_p, compressible, cell.coords, ue, pe)
    von_mises_stresses[cellid(cell)] = σᵥₘ
end

isdir("paraview") || mkdir("paraview")
vtk_grid("paraview/incompressible_center_grain", dh) do vtk
    vtk_point_data(vtk, dh, u)
    vtk_cell_data(vtk, pressures, "pressure")
    vtk_cell_data(vtk, von_mises_stresses, "vonMises")
    vtk_cellset(vtk, grid)
end

##################
# Makie plotting #
##################
logo_colormaps = Dict(
        "1" => :jpurple,
        "2" => :jred,
        "3" => :jred,
        "4" => :jblue,
        "5" => :jpurple,
        "6" => :jgreen,
    )

u_nodal = Ferrite._evaluate_at_grid_nodes(dh, u, :u)
data = [[Point2f(grid.nodes[i].x + u_nodal[i]) for i in c.nodes] for c in grid.cells]
ux_min, ux_max = extrema(u->u[1], u_nodal)
uy_min, uy_max = extrema(u->u[2], u_nodal)

size_cm = (10.0, 9.51)
size_pt = 72 .* size_cm ./ 2.54
f = Figure(; resolution=size_pt, backgroundcolor=:transparent)
ax = Axis(f[1,1];
    backgroundcolor=:transparent,
    xtickalign=1,
    ytickalign=1,
    xticks=0.0:0.2:1.0,
    yticks=0.0:0.2:1.0,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    bottomspinecolor=:white,
    rightspinecolor=:white,
    topspinecolor=:white,
    leftspinecolor=:white,
    xtickcolor=:white,
    ytickcolor=:white,
   )

colsize!(f.layout, 1, Aspect(1,(1.0+ux_min)/(1.0+uy_max))) # correct aspect ratio
σᵥₘ_limits = extrema(von_mises_stresses)

# plot in julia logo colors
for (cellset_name, colormap) in pairs(logo_colormaps)
    cellset = collect(getcellset(grid, cellset_name))
    _data = data[cellset]
    poly!(_data; strokecolor=:black, strokewidth=1, color=von_mises_stresses[cellset], transparency=true, colormap, colorrange=σᵥₘ_limits)
end

cb_layout = GridLayout()
f.layout[1,2] = cb_layout
for (i, colormap) in pairs((:jgreen, :jred, :jblue, :jpurple))
    cb = Colorbar(cb_layout[1,i];
         colormap,
         ticklabelsvisible=false,
         tickalign=1,
         limits=σᵥₘ_limits,
         bottomspinecolor=:white,
         rightspinecolor=:white,
         topspinecolor=:white,
         leftspinecolor=:white,
         tickcolor=:white,
         labelcolor=:white,
    )
    i == 4 && (cb.label = "von Mises stress")
end
colgap!(cb_layout, 0.0)
resize_to_layout!(f)

isdir("graphics") || mkdir("graphics")
save("graphics/incompressible_vonMises_logo.svg", f, pt_per_unit=1)
