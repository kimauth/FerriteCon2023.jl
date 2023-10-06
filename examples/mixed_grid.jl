# First we load Ferrite, and some other packages we need.
using Ferrite, FerriteGmsh, SparseArrays
using FerriteCon2023
using Colors, CairoMakie

grid = togrid("logo_quads.geo");

# we will need cellsets of Triangles + Quadrilaterals
tris = findall(cell -> cell isa Triangle, grid.cells)
addcellset!(grid, "triangles", Set(tris))
quads = findall(cell -> cell isa Quadrilateral, grid.cells)
addcellset!(grid, "quadrilaterals", Set(quads))

# Trial and test functions
dim = 2
order = 1 # linear interpolation
ip_tri = Lagrange{RefTriangle, order}()
qr_tri = QuadratureRule{RefTriangle}(1) # 1 quadrature point
cv_tri = CellValues(qr_tri, ip_tri^dim);

ip_quad = Lagrange{RefQuadrilateral, order}()
qr_quad = QuadratureRule{RefQuadrilateral}(2) # 2x2 quadrature points
cv_quad = CellValues(qr_quad, ip_quad^dim)

cellvalues = (cv_tri, cv_quad)

# ### Degrees of freedom
dh = DofHandler(grid)

# subdomains can only consist of a single element type
sdh_tri = SubDofHandler(dh, getcellset(grid, "triangles"))
add!(sdh_tri, :u, ip_tri^dim)

sdh_quad = SubDofHandler(dh, getcellset(grid, "quadrilaterals"))
add!(sdh_quad, :u, ip_quad^dim)

close!(dh);

# ### Boundary conditions
ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfaceset(grid, "bottom"), (x, t) -> 0.0, 2))
add!(ch, Dirichlet(:u, getfaceset(grid, "left"), (x, t) -> 0.0, 1))
add!(ch, Dirichlet(:u, getfaceset(grid, "top"), (x, t) -> 0.1, 2))
close!(ch);

E = 200e3 # Young's modulus [MPa]
ν = 0.3 # Poisson's ratio [-]
material = Elasticity(E/2(1+ν), E/3(1-2ν));

# #### Global assembly
function assemble_global!(assembler, dh, cellvalues, material)
    ## Create an assembler
    for (sdh, cv) in zip(dh.subdofhandlers, cellvalues)
        assemble_subdofhandler!(assembler, sdh, cv, material)
    end
end

function assemble_subdofhandler!(assembler, sdh, cv, material)
    n_basefuncs = getnbasefunctions(cv)
    ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    ## Loop over all cells
    for cell in CellIterator(sdh)
        ## Compute element contribution
        assemble_cell!(ke, cv, material, cell.coords)
        ## Assemble ke and fe into K and f
        assemble!(assembler, cell.dofs, ke)
    end
end

# ### Solution of the system
K = create_sparsity_pattern(dh)
f = zeros(ndofs(dh))
assembler = start_assemble(K, f)
assemble_global!(assembler, dh, cellvalues, material);

apply!(K, f, ch)
u = K \ f;

# ### Exporting to VTK
isdir("paraview") || mkdir("paraview")
vtk_grid("paraview/linear_elasticity_quads", dh) do vtk
    vtk_point_data(vtk, dh, u)
    vtk_cellset(vtk, grid)
end

##################
# Makie plotting #
##################
logo_colors = Dict(
        "1" => Colors.JULIA_LOGO_COLORS.purple,
        "2" => Colors.JULIA_LOGO_COLORS.red,
        "3" => Colors.JULIA_LOGO_COLORS.red,
        "4" => Colors.JULIA_LOGO_COLORS.blue,
        "5" => Colors.JULIA_LOGO_COLORS.purple,
        "6" => Colors.JULIA_LOGO_COLORS.green,
    )

f = Figure(; backgroundcolor=:transparent)
ax = Axis(f[1,1]; backgroundcolor=:transparent)
colsize!(f.layout, 1, Aspect(1,1)) # correct aspect ratio
hidedecorations!(ax); hidespines!(ax)
resize_to_layout!(f)

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
    xticksmirrored=true,
    yticksmirrored=true,
   )

colsize!(f.layout, 1, Aspect(1,(1.0+ux_min)/(1.0+uy_max))) # correct aspect ratio

for (cellset_name, color) in pairs(logo_colors)
    cellset = collect(getcellset(grid, cellset_name))
    _data = data[cellset]
    poly!(_data; strokecolor=:black, strokewidth=1, color, transparency=true)
end
isdir("graphics") || mkdir("graphics")
save("graphics/mixed_grid.svg", f)
