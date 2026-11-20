module SubdomainEmbeddedTests

using Test
using FiniteDiff
using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.MultiField
using GridapEmbedded
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: DifferentiableTriangulation

# Helper to nudge level-set nodal values away from exact zero (avoiding degenerate cuts in FDM)
function nudge_zeros!(φh)
  x_φ = get_free_dof_values(φh)
  idx = findall(isapprox(0.0; atol=1e-10), x_φ)
  x_φ[idx] .+= 100 * eps(eltype(x_φ))
  return φh
end

############################################################################################
# 1. Subdomain cutting basics, CSG, GhostSkeleton, merge_triangulations, dimension checks
############################################################################################
@testset "Subdomain basics & CSG (2D)" begin
  model = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (10, 10)))
  sq_ids_quads = [34:37; 44:47; 54:57; 64:67]
  sq_ids = vcat([2*i-1:2*i for i in sq_ids_quads]...)
  remaining_ids = setdiff(1:num_cells(model), sq_ids)
  Ω_φ = Triangulation(model, sq_ids)
  Ω_r = Triangulation(model, remaining_ids)

  V_φ = TestFESpace(Ω_φ, ReferenceFE(lagrangian, Float64, 1))
  φh1 = nudge_zeros!(interpolate(x -> x[1] - x[2] - 0.1, V_φ))
  φh2 = nudge_zeros!(interpolate(x -> sqrt((x[1] - 0.5)^2 + (x[2] - 0.5)^2) - 0.17, V_φ))

  geo1 = DiscreteGeometry(φh1, Ω_φ)
  geo2 = DiscreteGeometry(φh2, Ω_φ)
  geo_csg = union(geo1, geo2)
  geo_diff = setdiff(geo1, geo2)
  geo_comp = !geo1
  geo_ana = disk(0.17, x0=Point(0.5, 0.5))

  # Dimension check: cutting a lower-dimensional triangulation must fail
  Γ_bdry = BoundaryTriangulation(model)
  @test_throws AssertionError cut(Γ_bdry, geo_ana)
  @test_throws AssertionError cut_facets(Γ_bdry, geo_ana)

  # Cut with DiscreteGeometry & AnalyticalGeometry
  for geo in (geo1, geo_csg, geo_diff, geo_comp, geo_ana)
    cutgeo = cut(Ω_φ, geo)
    cutfacet = cut_facets(Ω_φ, geo)

    Ω_φ_act = Triangulation(cutgeo, ACTIVE)
    Ω_φ_phys = Triangulation(cutgeo, PHYSICAL)
    Ω_φ_out = Triangulation(cutgeo, PHYSICAL_OUT)
    Γ = EmbeddedBoundary(cutgeo)
    Γg = GhostSkeleton(cutgeo)

    @test get_background_model(Ω_φ_act) === model
    @test get_background_model(Ω_φ_phys) === model
    @test get_background_model(Γ) === model
    @test get_background_model(Γg) === model

    # Volume conservation strictly within Ω_φ
    vol_in = sum(∫(1) * Measure(Ω_φ_phys, 2))
    vol_out = sum(∫(1) * Measure(Ω_φ_out, 2))
    vol_φ = sum(∫(1) * Measure(Ω_φ, 2))
    @test vol_in + vol_out ≈ vol_φ atol=1e-12

    # All active cells of cutgeo must be a subset of sq_ids
    act_cells = get_glue(Ω_φ_act, Val(2)).tface_to_mface
    @test issubset(act_cells, sq_ids)

    # All GhostSkeleton facets must be strictly interior to Ω_φ (both cells in sq_ids)
    if num_cells(Γg) > 0
      g_plus = get_glue(Γg.⁺, Val(2)).tface_to_mface
      g_minus = get_glue(Γg.⁻, Val(2)).tface_to_mface
      @test issubset(g_plus, sq_ids)
      @test issubset(g_minus, sq_ids)
    end

    # merge_triangulations builds an injective BodyFittedTriangulation
    Ωact = merge_triangulations(Ω_φ_act, Ω_r)
    @test get_background_model(Ωact) === model
    @test sort(collect(get_glue(Ωact, Val(2)).tface_to_mface)) == sort(union(act_cells, remaining_ids))
    @test !isnothing(get_glue(Ωact, Val(2)).mface_to_tface)

    # Cut facets via Interface(Ω_r, Ω_φ) and Boundary(Ω_r)
    Λ = Interface(Ω_r, Ω_φ)
    Γ_int_in = Boundary(Λ.⁺, cutfacet, PHYSICAL_IN)
    Γ_int_out = Boundary(Λ.⁺, cutfacet, PHYSICAL_OUT)
    len_in = num_cells(Γ_int_in) > 0 ? sum(∫(1) * Measure(Γ_int_in, 2)) : 0.0
    len_out = num_cells(Γ_int_out) > 0 ? sum(∫(1) * Measure(Γ_int_out, 2)) : 0.0
    len_tot = sum(∫(1) * Measure(Λ.⁺, 2))
    @test len_in + len_out ≈ len_tot atol=1e-12

    # Also test CompositeTriangulation input (Boundary(Ω_r))
    Γ_bdry_r_out = Boundary(Boundary(Ω_r), cutfacet, PHYSICAL_OUT)
    len_out_r = num_cells(Γ_bdry_r_out) > 0 ? sum(∫(1) * Measure(Γ_bdry_r_out, 2)) : 0.0
    @test len_out_r ≈ len_out atol=1e-12
  end

  # Also test cut_facets(cutgeo, geo) and DiscreteGeometry(φh_full, model) on subdomain Ω_φ
  V_full = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1))
  φh_full = nudge_zeros!(interpolate(x -> x[1] - x[2] - 0.1, V_full))
  geo_full = DiscreteGeometry(φh_full, model)
  cutgeo_full = cut(Ω_φ, geo_full)
  cutfacet_full = cut_facets(cutgeo_full, geo_full)
  @test sum(∫(1) * Measure(Triangulation(cutgeo_full, PHYSICAL), 2)) ≈
        sum(∫(1) * Measure(Triangulation(cut(Ω_φ, geo1), PHYSICAL), 2)) atol=1e-12
  @test sum(∫(1) * Measure(Boundary(Interface(Ω_r, Ω_φ).⁺, cutfacet_full, PHYSICAL_OUT), 2)) ≈
        sum(∫(1) * Measure(Boundary(Interface(Ω_r, Ω_φ).⁺, cut_facets(Ω_φ, geo1), PHYSICAL_OUT), 2)) atol=1e-12
end

############################################################################################
# 2. Subdomain CutFEM Poisson PDE + Full AD verification (2D)
############################################################################################
@testset "Subdomain CutFEM Poisson & AD (2D)" begin
  n = 10
  h = 1.0 / n
  model = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (n, n)))

  sq_ids_quads = [34:37; 44:47; 54:57; 64:67]
  sq_ids = vcat([2*i-1:2*i for i in sq_ids_quads]...)
  remaining_ids = setdiff(1:num_cells(model), sq_ids)
  Ω_φ = Triangulation(model, sq_ids)
  Ω_r = Triangulation(model, remaining_ids)

  V_φ = TestFESpace(Ω_φ, ReferenceFE(lagrangian, Float64, 1))
  φh = nudge_zeros!(interpolate(((x, y),) -> x - y - 0.13, V_φ))

  geo = DiscreteGeometry(φh, Ω_φ)
  cutgeo = cut(Ω_φ, geo)
  cutfacet = cut_facets(Ω_φ, geo)

  # Manufactured solution (linear -> exact for order 1!)
  u(x) = x[1] + 2 * x[2] - 0.5
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  Ω_φ_act = Triangulation(cutgeo, ACTIVE)
  Ω_φ_phys = Triangulation(cutgeo, PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)
  Λ_rφ = Interface(Ω_r, Ω_φ)
  Γ_∂Ωr = Boundary(Λ_rφ.⁺, cutfacet, PHYSICAL_OUT)
  Γg = GhostSkeleton(cutgeo)

  n_Γ = get_normal_vector(Γ)
  n_Γg = get_normal_vector(Γg)
  n_∂Ωr = get_normal_vector(Γ_∂Ωr) # already outwards from Ω_r into Ω_φ

  order = 1
  degree = 2 * order
  dΩ_φ_phys = Measure(Ω_φ_phys, degree)
  dΩ_r = Measure(Ω_r, degree)
  dΓg = Measure(Γg, degree)
  quad = Quadrature(duffy, degree)
  dΓ = Measure(Γ, quad)
  dΓ_∂Ωr = Measure(Γ_∂Ωr, degree)

  Ωact = merge_triangulations(Ω_φ_act, Ω_r)
  V = TestFESpace(Ωact, ReferenceFE(lagrangian, Float64, order), conformity=:H1, dirichlet_tags="boundary")
  U = TrialFESpace(V, ud)

  γd = 10.0
  γg = 0.1

  a(u, v) =
    ∫(∇(v) ⋅ ∇(u)) * dΩ_φ_phys + ∫(∇(v) ⋅ ∇(u)) * dΩ_r +
    ∫((γd / h) * v * u - v * (n_Γ ⋅ ∇(u)) - (n_Γ ⋅ ∇(v)) * u) * dΓ +
    ∫((γd / h) * v * u - v * (n_∂Ωr ⋅ ∇(u)) - (n_∂Ωr ⋅ ∇(v)) * u) * dΓ_∂Ωr +
    ∫((γg * h) * jump(n_Γg ⋅ ∇(v)) * jump(n_Γg ⋅ ∇(u))) * dΓg

  l(v) =
    ∫(v * f) * dΩ_φ_phys + ∫(v * f) * dΩ_r +
    ∫((γd / h) * v * ud - (n_Γ ⋅ ∇(v)) * ud) * dΓ +
    ∫((γd / h) * v * ud - (n_∂Ωr ⋅ ∇(v)) * ud) * dΓ_∂Ωr

  op = AffineFEOperator(a, l, U, V)
  uh = solve(op)

  e = u - uh
  el2 = sqrt(sum(∫(e * e) * dΩ_φ_phys + ∫(e * e) * dΩ_r))
  eh1 = sqrt(sum(∫(e * e + ∇(e) ⋅ ∇(e)) * dΩ_φ_phys + ∫(e * e + ∇(e) ⋅ ∇(e)) * dΩ_r))
  @test el2 < 1e-10
  @test eh1 < 1e-9

  # --- Automatic Differentiation verification against FiniteDiff ---
  # Use a non-linear smooth field wh on Ωact so every gradient/trace term varies spatially
  wh = interpolate(x -> sin(π * x[1]) * cos(π * x[2]) + x[1]^2 + 0.5 * x[2]^2, V)

  diff_Ω_φ_phys = DifferentiableTriangulation(Ω_φ_phys, V_φ)
  diff_Γ = DifferentiableTriangulation(Γ, V_φ)
  diff_Γ_∂Ωr = DifferentiableTriangulation(Γ_∂Ωr, V_φ)

  ddiff_Ω_φ_phys = Measure(diff_Ω_φ_phys, degree)
  ddiff_Γ = Measure(diff_Γ, quad)
  ddiff_Γ_∂Ωr = Measure(diff_Γ_∂Ωr, degree)

  # (a) Bulk term on diff_Ω_φ_phys
  J_bulk(φ) = ∫(∇(wh) ⋅ ∇(wh) + wh * wh) * ddiff_Ω_φ_phys
  g_bulk_ad = assemble_vector(gradient(J_bulk, φh), V_φ)
  g_bulk_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
    cg = cut(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
    dΩp = Measure(Triangulation(cg, PHYSICAL), degree)
    sum(∫(∇(wh) ⋅ ∇(wh) + wh * wh) * dΩp)
  end
  @test norm(g_bulk_ad - g_bulk_fd, Inf) < 1e-7

  # (b) Embedded boundary term on diff_Γ (including moving normal!)
  function J_gamma(φ)
    n = get_normal_vector(diff_Γ)
    ∫((γd / h) * wh * wh - 2 * wh * (n ⋅ ∇(wh))) * ddiff_Γ
  end
  g_gamma_ad = assemble_vector(gradient(J_gamma, φh), V_φ)
  g_gamma_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
    cg = cut(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
    Γ_i = EmbeddedBoundary(cg)
    n_i = get_normal_vector(Γ_i)
    dΓ_i = Measure(Γ_i, quad)
    sum(∫((γd / h) * wh * wh - 2 * wh * (n_i ⋅ ∇(wh))) * dΓ_i)
  end
  @test norm(g_gamma_ad - g_gamma_fd, Inf) < 1e-6

  # (c) Cut background facet term on diff_Γ_∂Ωr (SubFacetBoundaryTriangulation glued to Ω_r!)
  function J_facet(φ)
    n = get_normal_vector(diff_Γ_∂Ωr)
    ∫((γd / h) * wh * wh - 2 * wh * (n ⋅ ∇(wh))) * ddiff_Γ_∂Ωr
  end
  g_facet_ad = assemble_vector(gradient(J_facet, φh), V_φ)
  g_facet_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
    cf = cut_facets(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
    Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
    nf_i = get_normal_vector(Γf_i)
    dΓf_i = Measure(Γf_i, degree)
    sum(∫((γd / h) * wh * wh - 2 * wh * (nf_i ⋅ ∇(wh))) * dΓf_i)
  end
  @test norm(g_facet_ad - g_facet_fd, Inf) < 1e-6

  # (d) Hessian of SubFacetBoundaryTriangulation term
  H_facet_ad = assemble_matrix(hessian(J_facet, φh), V_φ, V_φ)
  H_facet_fd = FiniteDiff.finite_difference_hessian(get_free_dof_values(φh)) do φ_vec
    cf = cut_facets(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
    Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
    nf_i = get_normal_vector(Γf_i)
    dΓf_i = Measure(Γf_i, degree)
    sum(∫((γd / h) * wh * wh - 2 * wh * (nf_i ⋅ ∇(wh))) * dΓf_i)
  end
  @test norm(Matrix(H_facet_ad) - H_facet_fd, Inf) < 1e-4

  # (e) Full combined functional j(u, v, φ) (including dΩ_r and dΓg!)
  function j_full(u_f, v_f, φ)
    n_Γ_ad = get_normal_vector(diff_Γ)
    n_∂Ωr_ad = get_normal_vector(diff_Γ_∂Ωr)
    ∫(∇(v_f) ⋅ ∇(u_f)) * ddiff_Ω_φ_phys + ∫(∇(v_f) ⋅ ∇(u_f)) * dΩ_r +
    ∫((γd / h) * v_f * u_f - v_f * (n_Γ_ad ⋅ ∇(u_f)) - (n_Γ_ad ⋅ ∇(v_f)) * u_f) * ddiff_Γ +
    ∫((γd / h) * v_f * u_f - v_f * (n_∂Ωr_ad ⋅ ∇(u_f)) - (n_∂Ωr_ad ⋅ ∇(v_f)) * u_f) * ddiff_Γ_∂Ωr +
    ∫((γg * h) * jump(n_Γg ⋅ ∇(v_f)) * jump(n_Γg ⋅ ∇(u_f))) * dΓg
  end
  g_full_ad = assemble_vector(gradient(φ -> j_full(wh, wh, φ), φh), V_φ)
  g_full_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
    g_i = DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ)
    cg = cut(Ω_φ, g_i)
    cf = cut_facets(Ω_φ, g_i)
    Ωp_i = Triangulation(cg, PHYSICAL)
    Γ_i = EmbeddedBoundary(cg)
    Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
    n_Γ_i = get_normal_vector(Γ_i)
    n_f_i = get_normal_vector(Γf_i)
    sum(
      ∫(∇(wh) ⋅ ∇(wh)) * Measure(Ωp_i, degree) + ∫(∇(wh) ⋅ ∇(wh)) * dΩ_r +
      ∫((γd / h) * wh * wh - 2 * wh * (n_Γ_i ⋅ ∇(wh))) * Measure(Γ_i, quad) +
      ∫((γd / h) * wh * wh - 2 * wh * (n_f_i ⋅ ∇(wh))) * Measure(Γf_i, degree) +
      ∫((γg * h) * jump(n_Γg ⋅ ∇(wh)) * jump(n_Γg ⋅ ∇(wh))) * dΓg
    )
  end
  @test norm(g_full_ad - g_full_fd, Inf) < 1e-6

  # (f) MultiField AD with V_φ on Ω_φ and V on Ωact
  Y = MultiFieldFESpace([V_φ, V])
  φwh = interpolate([φh, wh], Y)
  g_mf = assemble_vector(gradient(((φ, w),) -> j_full(w, w, φ), φwh), Y)
  g_w_ad = assemble_vector(gradient(w -> j_full(w, w, φh), wh), V)
  n_φ = num_free_dofs(V_φ)
  @test norm(g_mf[1:n_φ] - g_full_ad, Inf) < 1e-12
  @test norm(g_mf[n_φ+1:end] - g_w_ad, Inf) < 1e-12
end

############################################################################################
# 3. 3D Subdomain CutFEM + SubFacetBoundaryTriangulation AD (TET & HEX)
############################################################################################
@testset "Subdomain CutFEM & AD (3D)" begin
  n = 4
  h = 1.0 / n
  model = simplexify(CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (n, n, n)))

  # Interior 2x2x2 cube of hexes -> 6 tets per hex
  # Select cells whose centroid lies in [0.25, 0.75]^3
  cell_coords = get_cell_coordinates(model)
  sq_ids = findall(cell_coords) do pts
    c = sum(pts) / length(pts)
    all(0.25 .<= Tuple(c) .<= 0.75)
  end
  remaining_ids = setdiff(1:num_cells(model), sq_ids)
  Ω_φ = Triangulation(model, sq_ids)
  Ω_r = Triangulation(model, remaining_ids)

  V_φ = TestFESpace(Ω_φ, ReferenceFE(lagrangian, Float64, 1))
  φh = nudge_zeros!(interpolate(x -> x[1] - 0.5 * x[2] + 0.3 * x[3] - 0.41, V_φ))

  geo = DiscreteGeometry(φh, Ω_φ)
  cutgeo = cut(Ω_φ, geo)
  cutfacet = cut_facets(Ω_φ, geo)

  Ω_φ_act = Triangulation(cutgeo, ACTIVE)
  Ω_φ_phys = Triangulation(cutgeo, PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)
  Λ_rφ = Interface(Ω_r, Ω_φ)
  Γ_∂Ωr = Boundary(Λ_rφ.⁺, cutfacet, PHYSICAL_OUT)
  Γg = GhostSkeleton(cutgeo)

  Ωact = merge_triangulations(Ω_φ_act, Ω_r)
  order = 1
  degree = 2 * order
  u(x) = x[1] + 2 * x[2] - 3 * x[3] + 0.25
  f(x) = -Δ(u)(x)
  V = TestFESpace(Ωact, ReferenceFE(lagrangian, Float64, order), conformity=:H1, dirichlet_tags="boundary")
  U = TrialFESpace(V, u)

  n_Γ = get_normal_vector(Γ)
  n_Γg = get_normal_vector(Γg)
  n_∂Ωr = get_normal_vector(Γ_∂Ωr)

  dΩ_φ_phys = Measure(Ω_φ_phys, degree)
  dΩ_r = Measure(Ω_r, degree)
  dΓ = Measure(Γ, Quadrature(duffy, degree))
  dΓ_∂Ωr = Measure(Γ_∂Ωr, degree)
  dΓg = Measure(Γg, degree)

  γd = 10.0
  γg = 0.1
  a(u, v) =
    ∫(∇(v) ⋅ ∇(u)) * dΩ_φ_phys + ∫(∇(v) ⋅ ∇(u)) * dΩ_r +
    ∫((γd / h) * v * u - v * (n_Γ ⋅ ∇(u)) - (n_Γ ⋅ ∇(v)) * u) * dΓ +
    ∫((γd / h) * v * u - v * (n_∂Ωr ⋅ ∇(u)) - (n_∂Ωr ⋅ ∇(v)) * u) * dΓ_∂Ωr +
    ∫((γg * h) * jump(n_Γg ⋅ ∇(v)) * jump(n_Γg ⋅ ∇(u))) * dΓg

  l(v) =
    ∫(v * f) * dΩ_φ_phys + ∫(v * f) * dΩ_r +
    ∫((γd / h) * v * u - (n_Γ ⋅ ∇(v)) * u) * dΓ +
    ∫((γd / h) * v * u - (n_∂Ωr ⋅ ∇(v)) * u) * dΓ_∂Ωr

  uh = solve(AffineFEOperator(a, l, U, V))
  e = u - uh
  el2 = sqrt(sum(∫(e * e) * dΩ_φ_phys + ∫(e * e) * dΩ_r))
  @test el2 < 1e-10

  # 3D AD verification on bulk, embedded boundary, and cut background facets
  wh = interpolate(x -> x[1]^2 + 2 * x[2]^2 - x[3]^2 + x[1] * x[2] * x[3], V)
  diff_Ω_φ_phys = DifferentiableTriangulation(Ω_φ_phys, V_φ)
  diff_Γ = DifferentiableTriangulation(Γ, V_φ)
  diff_Γ_∂Ωr = DifferentiableTriangulation(Γ_∂Ωr, V_φ)

  ddiff_Ω_φ_phys = Measure(diff_Ω_φ_phys, degree)
  ddiff_Γ = Measure(diff_Γ, Quadrature(duffy, degree))
  ddiff_Γ_∂Ωr = Measure(diff_Γ_∂Ωr, degree)

  function J_3d(φ)
    n1 = get_normal_vector(diff_Γ)
    n2 = get_normal_vector(diff_Γ_∂Ωr)
    ∫(wh * wh + ∇(wh) ⋅ ∇(wh)) * ddiff_Ω_φ_phys +
    ∫(wh * wh - wh * (n1 ⋅ ∇(wh))) * ddiff_Γ +
    ∫(wh * wh - wh * (n2 ⋅ ∇(wh))) * ddiff_Γ_∂Ωr
  end
  g_3d_ad = assemble_vector(gradient(J_3d, φh), V_φ)
  g_3d_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
    g_i = DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ)
    cg = cut(Ω_φ, g_i)
    cf = cut_facets(Ω_φ, g_i)
    Ωp_i = Triangulation(cg, PHYSICAL)
    Γ_i = EmbeddedBoundary(cg)
    Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
    n1 = get_normal_vector(Γ_i)
    n2 = get_normal_vector(Γf_i)
    sum(
      ∫(wh * wh + ∇(wh) ⋅ ∇(wh)) * Measure(Ωp_i, degree) +
      ∫(wh * wh - wh * (n1 ⋅ ∇(wh))) * Measure(Γ_i, Quadrature(duffy, degree)) +
      ∫(wh * wh - wh * (n2 ⋅ ∇(wh))) * Measure(Γf_i, degree)
    )
  end
  @test norm(g_3d_ad - g_3d_fd, Inf) < 1e-6

  # Also test SubFacetBoundaryTriangulation AD on full-domain Boundary(cutfacet, PHYSICAL_IN)
  model_full = CartesianDiscreteModel((0, 1, 0, 1), (6, 6))
  V_full = TestFESpace(model_full, ReferenceFE(lagrangian, Float64, 1))
  φh_full = nudge_zeros!(interpolate(x -> sqrt((x[1] - 0.5)^2 + (x[2] - 0.1)^2) - 0.31, V_full))
  cf_full = cut_facets(model_full, DiscreteGeometry(φh_full, model_full))
  Γ_full = Boundary(cf_full, PHYSICAL_IN)
  diff_Γ_full = DifferentiableTriangulation(Γ_full, V_full)
  dΓ_full_ad = Measure(diff_Γ_full, 2)
  wh_full = interpolate(x -> x[1]^2 + x[2]^2 + sin(x[1]), V_full)
  g_bdry_ad = assemble_vector(gradient(φ -> ∫(wh_full) * dΓ_full_ad, φh_full), V_full)
  g_bdry_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh_full)) do φ_vec
    cf_i = cut_facets(model_full, DiscreteGeometry(FEFunction(V_full, φ_vec), model_full))
    sum(∫(wh_full) * Measure(Boundary(cf_i, PHYSICAL_IN), 2))
  end
  @test norm(g_bdry_ad - g_bdry_fd, Inf) < 1e-6
end

end # module
