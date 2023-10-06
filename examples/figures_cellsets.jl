using CairoMakie
using Colors

grid = togrid("logo_quads.geo");
data = [[Point2f(grid.nodes[i].x) for i in c.nodes] for c in grid.cells]

logo_colors = Dict(
        "1" => Colors.JULIA_LOGO_COLORS.purple,
        "2" => Colors.JULIA_LOGO_COLORS.red,
        "3" => Colors.JULIA_LOGO_COLORS.red,
        "4" => Colors.JULIA_LOGO_COLORS.blue,
        "5" => Colors.JULIA_LOGO_COLORS.purple,
        "6" => Colors.JULIA_LOGO_COLORS.green,
    )

isdir("graphics") || mkdir("graphics")

################
# base version #
################
f = Figure(; backgroundcolor=:transparent)
ax = Axis(f[1,1]; backgroundcolor=:transparent)
colsize!(f.layout, 1, Aspect(1,1)) # correct aspect ratio
hidedecorations!(ax); hidespines!(ax)
resize_to_layout!(f)

for (cellset_name, color) in pairs(logo_colors)
    cellset = getcellset(grid, cellset_name)
    _data = [[Point2f(grid.nodes[i].x) for i in grid.cells[cellidx].nodes] for cellidx in cellset]
    poly!(_data; strokecolor=:black, strokewidth=1, color, transparency=true)
end
save("graphics/base_logo.svg", f)

###################
# highlight quads #
###################
f = Figure(; backgroundcolor=:transparent)
ax = Axis(f[1,1]; backgroundcolor=:transparent)
colsize!(f.layout, 1, Aspect(1,1)) # correct aspect ratio
hidedecorations!(ax); hidespines!(ax)
resize_to_layout!(f)

for (cellset_name, color) in pairs(logo_colors)
    cellset = getcellset(grid, cellset_name)
    _data = [[Point2f(grid.nodes[i].x) for i in grid.cells[cellidx].nodes] for cellidx in cellset]
    poly!(_data; strokecolor=:grey, strokewidth=1, color, transparency=true)
end

highlight_cellset = getcellset(grid, "quadrilaterals")
other_cells = setdiff(1:getncells(grid), highlight_cellset)
poly!(data[collect(other_cells)]; color=(:black, 0.5), transparency=true)
poly!(data[collect(highlight_cellset)]; strokecolor=:black, strokewidth=2, color=:transparent, transparency=true)
save("graphics/highlight_quads.svg", f)

##########################
# highlight center grain #
##########################
f = Figure(; backgroundcolor=:transparent)
ax = Axis(f[1,1]; backgroundcolor=:transparent)
colsize!(f.layout, 1, Aspect(1,1)) # correct aspect ratio
hidedecorations!(ax); hidespines!(ax)
resize_to_layout!(f)

for (cellset_name, color) in pairs(logo_colors)
    cellset = getcellset(grid, cellset_name)
    _data = [[Point2f(grid.nodes[i].x) for i in grid.cells[cellidx].nodes] for cellidx in cellset]
    poly!(_data; strokecolor=:grey, strokewidth=1, color, transparency=true)
end

highlight_cellset = getcellset(grid, "green")
other_cells = setdiff(1:getncells(grid), highlight_cellset)
poly!(data[collect(other_cells)]; color=(:black, 0.5), transparency=true)
poly!(data[collect(highlight_cellset)]; strokecolor=:black, strokewidth=2, color=:transparent, transparency=true)
save("graphics/highlight_center.svg", f)

#######################
# highlight all other #
#######################
f = Figure(; backgroundcolor=:transparent)
ax = Axis(f[1,1]; backgroundcolor=:transparent)
colsize!(f.layout, 1, Aspect(1,1)) # correct aspect ratio
hidedecorations!(ax); hidespines!(ax)
resize_to_layout!(f)

for (cellset_name, color) in pairs(logo_colors)
    cellset = getcellset(grid, cellset_name)
    _data = [[Point2f(grid.nodes[i].x) for i in grid.cells[cellidx].nodes] for cellidx in cellset]
    poly!(_data; strokecolor=:grey, strokewidth=1, color, transparency=true)
end

other_cells = union(getcellset(grid, "quadrilaterals"), getcellset(grid, "6"))
highlight_cellset = setdiff(1:getncells(grid), other_cells)
poly!(data[collect(other_cells)]; color=(:black, 0.5), transparency=true)
poly!(data[collect(highlight_cellset)]; strokecolor=:black, strokewidth=2, color=:transparent, transparency=true)
save("graphics/highlight_all_other.svg", f)

##########################
# base version triangles #
##########################
grid = togrid("logo.geo");
data = [[Point2f(grid.nodes[i].x) for i in c.nodes] for c in grid.cells]

f = Figure(; backgroundcolor=:transparent)
ax = Axis(f[1,1]; backgroundcolor=:transparent)
colsize!(f.layout, 1, Aspect(1,1)) # correct aspect ratio
hidedecorations!(ax); hidespines!(ax)
resize_to_layout!(f)

for (cellset_name, color) in pairs(logo_colors)
    cellset = getcellset(grid, cellset_name)
    _data = [[Point2f(grid.nodes[i].x) for i in grid.cells[cellidx].nodes] for cellidx in cellset]
    poly!(_data; strokecolor=:black, strokewidth=1, color, transparency=true)
end
save("graphics/base_logo_tris.svg", f)
