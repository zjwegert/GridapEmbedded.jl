using Gridap
using Gridap.ReferenceFEs, Gridap.Geometry
using GridapEmbedded
using GridapEmbedded.LevelSetCutters

n = 10
h = 1 / n
model = CartesianDiscreteModel((0,1,0,1),(n,n))
model = simplexify(model)

sq_ids_quads = [34:37; 44:47; 54:57; 64:67]
sq_ids = vcat([2*i-1:2*i for i in sq_ids_quads]...)
remaining_ids = setdiff(1:length(model.grid.cell_node_ids), sq_ids)
Ω_φ = Triangulation(model, sq_ids)
Ω_r = Triangulation(model, remaining_ids)

V_φ = TestFESpace(Ω_φ, ReferenceFE(lagrangian,Float64,1))
φh = interpolate(((x,y),) -> x - y - 0.13, V_φ)

geo = DiscreteGeometry(φh, Ω_φ)
cutgeo = cut(Ω_φ, geo)
cutfacet = cut_facets(Ω_φ, geo)

# Manufactured solution
u(x) = x[1] + x[2]
∇u(x) = ∇(u)(x)
f(x) = -Δ(u)(x)
ud(x) = u(x)

# Integration meshes
Ω_φ_act = Triangulation(cutgeo,ACTIVE)
Ω_φ_phys = Triangulation(cutgeo,PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)
Λ_rφ = Interface(Ω_r, Ω_φ)
Γ_∂Ωr_cap_Ω_φ_physᶜ = Boundary(Λ_rφ.⁺, cutfacet, PHYSICAL_OUT) # glued to Ω_r with outward normal pointing into Ω_φ
Γg = GhostSkeleton(cutgeo)

# Setup normal vectors
n_Γ = get_normal_vector(Γ)
n_Γg = get_normal_vector(Γg)
n_∂Ω_r_into_Ω_φ = get_normal_vector(Γ_∂Ωr_cap_Ω_φ_physᶜ)

# Setup Lebesgue measures
order = 1
degree = 2*order
dΩ_φ_phys = Measure(Ω_φ_phys,degree)
dΩ_r = Measure(Ω_r,degree)
dΓg = Measure(Γg,degree)
quad = Quadrature(duffy,degree)
dΓ = Measure(Γ,quad) # 4 integration points when order = 1
dΓ_∂Ωr_cap_Ω_φ_physᶜ = Measure(Γ_∂Ωr_cap_Ω_φ_physᶜ,degree)

# Setup FESpace
Ωact = merge_triangulations(Ω_φ_act, Ω_r)
V = TestFESpace(Ωact,ReferenceFE(lagrangian,Float64,order),conformity=:H1, dirichlet_tags="boundary")
U = TrialFESpace(V,ud)

# Weak form Nitsche + ghost penalty (CutFEM paper Sect. 6.1)
const γd = 10.0
const γg = 0.1

a(u,v) =
  ∫( ∇(v)⋅∇(u) ) * dΩ_φ_phys + ∫( ∇(v)⋅∇(u) ) * dΩ_r +
  ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ +
  ∫( (γd/h)*v*u  - v*(n_∂Ω_r_into_Ω_φ⋅∇(u)) - (n_∂Ω_r_into_Ω_φ⋅∇(v))*u ) * dΓ_∂Ωr_cap_Ω_φ_physᶜ +
  ∫( (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u)) ) * dΓg

l(v) =
  ∫( v*f ) * dΩ_φ_phys + ∫( v*f ) * dΩ_r +
  ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ +
  ∫( (γd/h)*v*ud - (n_∂Ω_r_into_Ω_φ⋅∇(v))*ud ) * dΓ_∂Ωr_cap_Ω_φ_physᶜ

# FE problem
op = AffineFEOperator(a,l,U,V)
uh = solve(op)

e = u - uh

# Postprocess
l2(u) = sqrt(sum( ∫( u*u )*dΩ_φ_phys + ∫( u*u )*dΩ_r ))
h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ_φ_phys + ∫( u*u + ∇(u)⋅∇(u) )*dΩ_r ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)
@show el2/ul2 eh1/uh1

# Differentiation
using GridapEmbedded.LevelSetCutters: DifferentiableTriangulation

diff_Ω_φ_phys = DifferentiableTriangulation(Ω_φ_phys, V_φ)
diff_Γ = DifferentiableTriangulation(Γ, V_φ)
diff_Γ_∂Ωr_cap_Ω_φ_physᶜ = DifferentiableTriangulation(Γ_∂Ωr_cap_Ω_φ_physᶜ, V_φ)

ddiff_Ω_φ_phys = Measure(diff_Ω_φ_phys,degree)
ddiff_Γ = Measure(diff_Γ,quad)
ddiff_Γ_∂Ωr_cap_Ω_φ_physᶜ = Measure(diff_Γ_∂Ωr_cap_Ω_φ_physᶜ,degree)

function j(u,v,φ)
  n_Γ_ad = get_normal_vector(diff_Γ)
  n_∂Ω_r_ad = get_normal_vector(diff_Γ_∂Ωr_cap_Ω_φ_physᶜ)
  ∫( ∇(v)⋅∇(u) ) * ddiff_Ω_φ_phys + ∫( ∇(v)⋅∇(u) ) * dΩ_r +
  ∫( (γd/h)*v*u  - v*(n_Γ_ad⋅∇(u)) - (n_Γ_ad⋅∇(v))*u ) * ddiff_Γ +
  ∫( (γd/h)*v*u  - v*(n_∂Ω_r_ad⋅∇(u)) - (n_∂Ω_r_ad⋅∇(v))*u ) * ddiff_Γ_∂Ωr_cap_Ω_φ_physᶜ +
  ∫( (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u)) ) * dΓg
end

grad = gradient(φ->j(uh,uh,φ), φh)
grad_vec = assemble_vector(grad, V_φ)
@show length(grad_vec) norm(grad_vec)

# Finite difference verification of AD (combined + per-term)
using FiniteDiff
using Test

# Use a non-linear field wh so gradients and traces vary non-trivially across cells
wh = interpolate(x -> sin(π*x[1])*cos(π*x[2]) + x[1]^2 + 0.5*x[2]^2, V)

# 1. Bulk term on diff_Ω_φ_phys
g_bulk_ad = assemble_vector(gradient(φ -> ∫( ∇(wh)⋅∇(wh) ) * ddiff_Ω_φ_phys, φh), V_φ)
g_bulk_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
  cg = cut(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
  sum(∫( ∇(wh)⋅∇(wh) ) * Measure(Triangulation(cg, PHYSICAL), degree))
end
@show norm(g_bulk_ad - g_bulk_fd, Inf)
@test norm(g_bulk_ad - g_bulk_fd, Inf) < 1e-6

# 2. Embedded boundary term on diff_Γ
function j_Γ(φ)
  n_Γ_ad = get_normal_vector(diff_Γ)
  ∫( (γd/h)*wh*wh - 2*wh*(n_Γ_ad⋅∇(wh)) ) * ddiff_Γ
end
g_Γ_ad = assemble_vector(gradient(j_Γ, φh), V_φ)
g_Γ_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
  cg = cut(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
  Γ_i = EmbeddedBoundary(cg)
  n_i = get_normal_vector(Γ_i)
  sum(∫( (γd/h)*wh*wh - 2*wh*(n_i⋅∇(wh)) ) * Measure(Γ_i, quad))
end
@show norm(g_Γ_ad - g_Γ_fd, Inf)
@test norm(g_Γ_ad - g_Γ_fd, Inf) < 1e-6

# 3. Cut background facet term on diff_Γ_∂Ωr_cap_Ω_φ_physᶜ (SubFacetBoundaryTriangulation)
function j_∂Ωr(φ)
  n_∂Ω_r_ad = get_normal_vector(diff_Γ_∂Ωr_cap_Ω_φ_physᶜ)
  ∫( (γd/h)*wh*wh - 2*wh*(n_∂Ω_r_ad⋅∇(wh)) ) * ddiff_Γ_∂Ωr_cap_Ω_φ_physᶜ
end
g_∂Ωr_ad = assemble_vector(gradient(j_∂Ωr, φh), V_φ)
g_∂Ωr_fd = FiniteDiff.finite_difference_gradient(get_free_dof_values(φh)) do φ_vec
  cf = cut_facets(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
  Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
  nf_i = get_normal_vector(Γf_i)
  sum(∫( (γd/h)*wh*wh - 2*wh*(nf_i⋅∇(wh)) ) * Measure(Γf_i, degree))
end
@show norm(g_∂Ωr_ad - g_∂Ωr_fd, Inf)
@test norm(g_∂Ωr_ad - g_∂Ωr_fd, Inf) < 1e-6

# 4. Combined functional j(wh, wh, φ)
g_full_ad = assemble_vector(gradient(φ -> j(wh, wh, φ), φh), V_φ)
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
    ∫( ∇(wh)⋅∇(wh) ) * Measure(Ωp_i, degree) + ∫( ∇(wh)⋅∇(wh) ) * dΩ_r +
    ∫( (γd/h)*wh*wh - 2*wh*(n_Γ_i⋅∇(wh)) ) * Measure(Γ_i, quad) +
    ∫( (γd/h)*wh*wh - 2*wh*(n_f_i⋅∇(wh)) ) * Measure(Γf_i, degree) +
    ∫( (γg*h)*jump(n_Γg⋅∇(wh))*jump(n_Γg⋅∇(wh)) ) * dΓg
  )
end
@show norm(g_full_ad - g_full_fd, Inf)
@test norm(g_full_ad - g_full_fd, Inf) < 1e-6