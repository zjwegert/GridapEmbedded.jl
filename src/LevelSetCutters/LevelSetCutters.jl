module LevelSetCutters

using AbstractTrees
import GridapEmbedded.CSG
using GridapEmbedded.CSG: Node, Leaf, replace_data
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: similar_geometry
import GridapEmbedded.CSG: compatible_geometries

using GridapEmbedded.Interfaces
import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: cut_facets
import GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut
import GridapEmbedded.Interfaces: compute_bgfacet_to_inoutcut
using GridapEmbedded.Interfaces: Simplex
using GridapEmbedded.Interfaces: UNDEF
using GridapEmbedded.Interfaces: merge_sub_face_data
using GridapEmbedded.Interfaces: compute_inoutcut
using GridapEmbedded.Interfaces: SubCellTriangulation, SubFacetTriangulation, SubFacetBoundaryTriangulation

using LinearAlgebra
using MiniQhull

using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Arrays: testitem, return_type
using Gridap.Fields
using Gridap.Helpers
using Gridap.Geometry
using Gridap.CellData
using Gridap.Polynomials
using Gridap.Visualization
using Gridap.FESpaces
using Gridap.MultiField

using GridapDistributed: DistributedMultiFieldFESpace

export LevelSetCutter
export AnalyticalGeometry
export doughnut
export popcorn
export tube
export olympic_rings
export sphere
export disk
export cylinder
export plane
export square
export quadrilateral
export cube
export discretize
export DiscreteGeometry

include("AnalyticalGeometries.jl")

include("DiscreteGeometries.jl")

include("LookupTables.jl")

include("CutTriangulations.jl")

include("DifferentiableTriangulations.jl")

"""
    struct LevelSetCutter <: Cutter end

Cutter for `DiscreteGeometry` and `AnalyticalGeometry`.

## Usage

    cut(background::DiscreteModel,geom::AnalyticalGeometry)
    cut(background::DiscreteModel,geom::DiscreteGeometry)
    cut_facets(background::DiscreteModel,geom::AnalyticalGeometry)
    cut_facets(background::DiscreteModel,geom::DiscreteGeometry)

"""
struct LevelSetCutter <: Cutter end

function cut(cutter::LevelSetCutter,background::DiscreteModel,geom)
  data = _cut_ls(background,geom)
  EmbeddedDiscretization(background, data..., geom)
end

function cut(background::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  cut(cutter,background,geom)
end

function cut(background::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  cut(cutter,background,geom)
end

function cut(cutter::LevelSetCutter,background::Triangulation,geom)
  bgmodel = get_background_model(background)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(background) == D """
  Subdomain cutting requires a Triangulation of the same cell dimension ($D) as its background DiscreteModel.
  """
  data = _cut_ls(background,geom)
  EmbeddedDiscretization(bgmodel, data..., geom)
end

function cut(background::Triangulation,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  cut(cutter,background,geom)
end

function cut(background::Triangulation,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  cut(cutter,background,geom)
end

function _cut_ls(model::DiscreteModel,geom)
  grid = get_grid(model)
  _cut_ls(grid,geom)
end

function _cut_ls(trian::Triangulation,geom)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(trian) == D
  glue = get_glue(trian,Val(D))
  @check isa(glue,FaceToFaceGlue) && !isnothing(glue.mface_to_tface) """
  Subdomain cutting requires an injective body-fitted Triangulation of the background DiscreteModel.
  """
  tcell_to_bgcell = glue.tface_to_mface
  if isa(tcell_to_bgcell,IdentityVector) && length(tcell_to_bgcell) == num_cells(bgmodel)
    return _cut_ls(bgmodel,geom)
  end
  grid = _get_subdomain_grid(trian)
  ls_to_tcell_to_inoutcut, subcells_local, ls_to_cell_to_inout, subfacets_local, ls_to_facet_to_inout, oid_to_ls = _cut_ls(grid,geom)
  n_bgcells = num_cells(bgmodel)
  ls_to_bgcell_to_inoutcut = map(ls_to_tcell_to_inoutcut) do tcell_to_ioc
    bgcell_to_ioc = fill(Int8(UNDEF),n_bgcells)
    bgcell_to_ioc[tcell_to_bgcell] .= tcell_to_ioc
    bgcell_to_ioc
  end
  scell_to_bgcell = collect(Int32,lazy_map(Reindex(tcell_to_bgcell),subcells_local.cell_to_bgcell))
  subcells = SubCellData(
    subcells_local.cell_to_points,
    scell_to_bgcell,
    subcells_local.point_to_coords,
    subcells_local.point_to_rcoords)
  sfacet_to_bgcell = collect(Int32,lazy_map(Reindex(tcell_to_bgcell),subfacets_local.facet_to_bgcell))
  subfacets = SubFacetData(
    subfacets_local.facet_to_points,
    subfacets_local.facet_to_normal,
    sfacet_to_bgcell,
    subfacets_local.point_to_coords,
    subfacets_local.point_to_rcoords)
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, subfacets, ls_to_facet_to_inout, oid_to_ls
end

function _cut_ls(grid::Grid,geom)
  out = initial_sub_triangulation(grid,geom)
  subcells0, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = out
  out2 = cut_sub_triangulation_with_boundary_several_levelsets(subcells0,ls_to_point_to_value)
  subcells, ls_to_cell_to_inout, subfacets, ls_to_facet_to_inout  = out2
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, subfacets, ls_to_facet_to_inout, oid_to_ls
end

function cut_facets(cutter::LevelSetCutter,background::DiscreteModel,geom)
  data = _cut_ls_facets(background,geom)
  EmbeddedFacetDiscretization(background, data..., geom)
end

function cut_facets(background::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  cut_facets(cutter,background,geom)
end

function cut_facets(background::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  cut_facets(cutter,background,geom)
end

function cut_facets(cutter::LevelSetCutter,background::Triangulation,geom)
  bgmodel = get_background_model(background)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(background) == D """
  Subdomain facet cutting requires a Triangulation of the same cell dimension ($D) as its background DiscreteModel.
  """
  data = _cut_ls_facets(background,geom)
  EmbeddedFacetDiscretization(bgmodel, data..., geom)
end

function cut_facets(background::Triangulation,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  cut_facets(cutter,background,geom)
end

function cut_facets(background::Triangulation,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  cut_facets(cutter,background,geom)
end

function cut_facets(cut::EmbeddedDiscretization,geo::CSG.Geometry)
  if !isempty(cut.ls_to_bgcell_to_inoutcut)
    mask = map(!=(UNDEF),first(cut.ls_to_bgcell_to_inoutcut))
    if !all(mask)
      return cut_facets(Triangulation(cut.bgmodel,mask),geo)
    end
  end
  cut_facets(cut.bgmodel,geo)
end

function cut_facets(cut::EmbeddedFacetDiscretization,args...)
  cut
end

function _cut_ls_facets(model::DiscreteModel,geom)
  D = num_cell_dims(model)
  grid = Grid(ReferenceFE{D-1},model)
  _cut_ls_facets(grid,geom)
end

function _subdomain_facet_grid(trian::Triangulation)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(trian) == D
  glue = get_glue(trian,Val(D))
  @check isa(glue,FaceToFaceGlue) && !isnothing(glue.mface_to_tface) """
  Subdomain facet cutting requires an injective body-fitted Triangulation of the background DiscreteModel.
  """
  tcell_to_bgcell = glue.tface_to_mface
  topo = get_grid_topology(bgmodel)
  cell_to_facets = Table(get_faces(topo,D,D-1))
  n_bgfacets = num_facets(bgmodel)
  bgfacet_in_trian = fill(false,n_bgfacets)
  for bgcell in tcell_to_bgcell
    for p in cell_to_facets.ptrs[bgcell]:(cell_to_facets.ptrs[bgcell+1]-1)
      bgfacet_in_trian[cell_to_facets.data[p]] = true
    end
  end
  tfacet_to_bgfacet = findall(bgfacet_in_trian)
  bgfacet_grid = Grid(ReferenceFE{D-1},bgmodel)
  tfacet_grid = GridPortion(bgfacet_grid,tfacet_to_bgfacet)
  cgrid = _get_subdomain_grid(trian)
  if isa(cgrid,Geometry.GridPortion)
    @check tfacet_grid.node_to_parent_node == cgrid.node_to_parent_node
  end
  return tfacet_grid, tfacet_to_bgfacet, n_bgfacets
end

function _cut_ls_facets(trian::Triangulation,geom)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(trian) == D
  glue = get_glue(trian,Val(D))
  @check isa(glue,FaceToFaceGlue) && !isnothing(glue.mface_to_tface)
  tcell_to_bgcell = glue.tface_to_mface
  if isa(tcell_to_bgcell,IdentityVector) && length(tcell_to_bgcell) == num_cells(bgmodel)
    return _cut_ls_facets(bgmodel,geom)
  end
  tfacet_grid, tfacet_to_bgfacet, n_bgfacets = _subdomain_facet_grid(trian)
  ls_to_tfacet_to_inoutcut, subfacets_local, ls_to_subfacet_to_inout, oid_to_ls = _cut_ls_facets(tfacet_grid,geom)
  ls_to_bgfacet_to_inoutcut = map(ls_to_tfacet_to_inoutcut) do tfacet_to_ioc
    bgfacet_to_ioc = fill(Int8(UNDEF),n_bgfacets)
    bgfacet_to_ioc[tfacet_to_bgfacet] .= tfacet_to_ioc
    bgfacet_to_ioc
  end
  sfacet_to_bgfacet = collect(Int32,lazy_map(Reindex(tfacet_to_bgfacet),subfacets_local.cell_to_bgcell))
  subfacets = SubCellData(
    subfacets_local.cell_to_points,
    sfacet_to_bgfacet,
    subfacets_local.point_to_coords,
    subfacets_local.point_to_rcoords)
  ls_to_bgfacet_to_inoutcut, subfacets, ls_to_subfacet_to_inout, oid_to_ls
end

function _cut_ls_facets(grid::Grid,geom)
  out = initial_sub_triangulation(grid,geom)
  subcells0, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = out
  out2 = cut_sub_triangulation_several_levelsets(subcells0,ls_to_point_to_value)
  subcells, ls_to_cell_to_inout  = out2
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, oid_to_ls
end

function compute_bgcell_to_inoutcut(model::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  compute_bgcell_to_inoutcut(cutter,model,geom)
end

function compute_bgcell_to_inoutcut(model::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  compute_bgcell_to_inoutcut(cutter,model,geom)
end

function compute_bgcell_to_inoutcut(::LevelSetCutter,model::DiscreteModel,geom)
  grid = get_grid(model)
  _compute_bgcell_to_inoutcut(grid,geom)
end

function compute_bgcell_to_inoutcut(trian::Triangulation,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  compute_bgcell_to_inoutcut(cutter,trian,geom)
end

function compute_bgcell_to_inoutcut(trian::Triangulation,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  compute_bgcell_to_inoutcut(cutter,trian,geom)
end

function compute_bgcell_to_inoutcut(::LevelSetCutter,trian::Triangulation,geom)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(trian) == D
  glue = get_glue(trian,Val(D))
  @check isa(glue,FaceToFaceGlue) && !isnothing(glue.mface_to_tface)
  tcell_to_bgcell = glue.tface_to_mface
  if isa(tcell_to_bgcell,IdentityVector) && length(tcell_to_bgcell) == num_cells(bgmodel)
    return _compute_bgcell_to_inoutcut(get_grid(bgmodel),geom)
  end
  tcell_to_ioc = _compute_bgcell_to_inoutcut(_get_subdomain_grid(trian),geom)
  bgcell_to_ioc = fill(Int8(UNDEF),num_cells(bgmodel))
  bgcell_to_ioc[tcell_to_bgcell] .= tcell_to_ioc
  return bgcell_to_ioc
end

function compute_bgfacet_to_inoutcut(model::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,model,geom)
end

function compute_bgfacet_to_inoutcut(model::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,model,geom)
end

function compute_bgfacet_to_inoutcut(::LevelSetCutter,model::DiscreteModel,geom)
  D = num_cell_dims(model)
  grid = Grid(ReferenceFE{D-1},model)
  _compute_bgcell_to_inoutcut(grid,geom)
end

function compute_bgfacet_to_inoutcut(trian::Triangulation,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,trian,geom)
end

function compute_bgfacet_to_inoutcut(trian::Triangulation,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,trian,geom)
end

function compute_bgfacet_to_inoutcut(::LevelSetCutter,trian::Triangulation,geom)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(trian) == D
  glue = get_glue(trian,Val(D))
  @check isa(glue,FaceToFaceGlue) && !isnothing(glue.mface_to_tface)
  tcell_to_bgcell = glue.tface_to_mface
  if isa(tcell_to_bgcell,IdentityVector) && length(tcell_to_bgcell) == num_cells(bgmodel)
    return compute_bgfacet_to_inoutcut(LevelSetCutter(),bgmodel,geom)
  end
  tfacet_grid, tfacet_to_bgfacet, n_bgfacets = _subdomain_facet_grid(trian)
  tfacet_to_ioc = _compute_bgcell_to_inoutcut(tfacet_grid,geom)
  bgfacet_to_ioc = fill(Int8(UNDEF),n_bgfacets)
  bgfacet_to_ioc[tfacet_to_bgfacet] .= tfacet_to_ioc
  return bgfacet_to_ioc
end

function _compute_bgcell_to_inoutcut(grid,geom)
  ls_to_cell_to_inoutcut, tree, oid_to_ls = _compute_ls_to_bgcell_to_inoutcut(grid,geom)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = oid_to_ls[oid]
    cell_to_inoutcut = ls_to_cell_to_inoutcut[ls]
    cell_to_inoutcut, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end

function _compute_ls_to_bgcell_to_inoutcut(grid,geom::AnalyticalGeometry)
  _compute_ls_to_bgcell_to_inoutcut(grid,discretize(geom,grid))
end

function _compute_ls_to_bgcell_to_inoutcut(grid,geom::DiscreteGeometry)
  ugrid = UnstructuredGrid(grid)
  tree = get_tree(geom)
  ls_to_point_to_value, oid_to_ls = _find_unique_leaves(tree)
  ls_to_point_to_value = _adapt_point_values_to_grid(ls_to_point_to_value,grid)
  p = _check_and_get_polytope(ugrid)
  table = LookupTable(p)
  cell_to_points = get_cell_node_ids(ugrid)
  ls_to_cell_to_inoutcut = [
    _compute_in_out_or_cut(table,cell_to_points,point_to_value)
    for point_to_value in ls_to_point_to_value]
  ls_to_cell_to_inoutcut, tree, oid_to_ls
end

end #module
