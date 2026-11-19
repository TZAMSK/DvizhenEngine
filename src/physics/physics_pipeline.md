# Physics Engine Pipeline — `sliyan_engine`

A code-oriented guide to the **structs, functions, execution order, and mathematics** of the current Rust rigid-body physics implementation.

> **Scope:** This describes the implementation shared in the conversation: `collider.rs`, `collider_builder.rs`, `world_collider.rs`, the GJK/EPA collision module, `rigid_body.rs`, and `physics_world.rs`. Some names retain their current spelling (`intertia`, `orentation`, `drag_coefficience`). The GJK/EPA description covers the algorithm and the previously shared implementation; it is not a proof that all edge cases are correct.

## 1. Architecture at a glance

```text
Scene::Shape
 ├── mesh.vertices                 (local geometry)
 ├── shape_type                    (Circle / Cube / other)
 └── transform.scale
          │
          ▼
make_collider_parts(shape)         collider_builder.rs
 ├── SphereCollider                collider.rs
 ├── BoxCollider                   collider.rs
 └── ConvexHullCollider            collider.rs
          │
          ▼
Vec<ColliderPart>
 ├── collider: Box<dyn Collider>
 ├── local_position: Vec3
 └── local_rotation: Quat
          │
          ▼
RigidBody                          rigid_body.rs
 ├── position / orentation
 ├── velocity / angular_velocity
 ├── force / inv_mass
 ├── inv_inertia_local
 └── parts: Vec<ColliderPart>
          │
          ▼
PhysicsWorld::add_body(body)       physics_world.rs
          │
          ▼
PhysicsWorld::step(dt)
 ├── integrate_velocities(dt)
 ├── find_body_contacts(...)
 │    ├── to_world(body, part) → WorldCollider
 │    ├── bounding_radius() broadphase
 │    └── collide(a, b)
 │         ├── GJK → intersecting simplex?
 │         └── EPA → Contact { point, normal, depth }
 ├── find_ground_contacts(...)
 │    └── WorldCollider::support(downward directions)
 ├── solve(contacts)
 │    ├── correct_penetration(...) once per contact
 │    └── resolve_contact(...) 8 passes
 └── integrate_positions(dt)
      ├── position update
      └── quaternion orientation update
```

### Coordinate systems

- **Mesh local space:** `shape.mesh.vertices` before scale is baked into the collider.
- **Collider local space:** the coordinates consumed by `Collider::support`; the factory bakes shape scale into these coordinates.
- **Rigid-body local space:** each `ColliderPart` has an offset and rotation relative to the body.
- **World space:** `RigidBody.position` and `orentation` determine where the collider is in the simulation. GJK/EPA and contact resolution operate on world-space geometry/contact data.

---

## 2. Collider geometry — `collider.rs`

### `trait Collider`

| Method | Output | Purpose |
|---|---|---|
| `support(direction: Vec3) -> Vec3` | Farthest local point | GJK/EPA and ground support sampling |
| `intertia(mass: f32) -> Mat3` | Local inertia tensor | Rigid-body mass properties |
| `bounding_radius() -> f32` | Origin-centered enclosing radius | Cheap broadphase rejection |

**Support mapping:**

$$
S_C(\mathbf d)=\underset{\mathbf p\in C}{\\mathrm{argmax}}\;(\mathbf p\cdot\mathbf d)
$$

It returns the point with the greatest projection along direction $\mathbf d$. The returned point is in the collider's **local space**.

### `SphereCollider { radius }`

**Support:**

$$
S(\mathbf d)=r\frac{\mathbf d}{\|\mathbf d\|},\qquad \|\mathbf d\|>0
$$

For nearly zero-length directions, the implementation returns the origin.

**Solid sphere inertia:**

$$
I=\frac25 mr^2\mathbf 1_3
$$

**Bounding radius:** $R=r$.

### `BoxCollider { half: Vec3 }`

For half-extents $(h_x,h_y,h_z)$, choose the corner whose coordinates have the signs of the search direction:

$$
S(\mathbf d)=\begin{pmatrix}h_x\\mathrm{sgn}(d_x)\\h_y\\mathrm{sgn}(d_y)\\h_z\\mathrm{sgn}(d_z)\end{pmatrix}
$$

In this code, a zero direction component chooses the **positive** half-extent.

**Solid cuboid inertia:**

$$
I=\begin{bmatrix}
\frac m3(h_y^2+h_z^2)&0&0\\
0&\frac m3(h_x^2+h_z^2)&0\\
0&0&\frac m3(h_x^2+h_y^2)
\end{bmatrix}
$$

**Bounding radius:** $R=\sqrt{h_x^2+h_y^2+h_z^2}$.

### `ConvexHullCollider { points, half, radius }`

`new(points)` computes:

1. Componentwise minimum/maximum vertex positions, `lo` and `hi`.
2. AABB half-extents $\mathbf h=(\mathbf{hi}-\mathbf{lo})/2$.
3. Origin-centered radius $R=\max_i\|\mathbf p_i\|$.

`support(direction)` scans all points, choosing the vertex with the largest dot product:

$$
S(\mathbf d)=\underset{i}{\\mathrm{argmax}}\;(\mathbf p_i\cdot\mathbf d)
$$

**Cost:** $O(V)$ for $V$ input vertices. It does not require an explicit hull-face mesh.

`intertia(mass)` **approximates** the convex hull as a solid AABB-sized box. It is not the exact inertia tensor of an arbitrary convex polyhedron.

---

## 3. Collider creation — `collider_builder.rs`

### `struct ColliderPart`

```rust
pub struct ColliderPart {
    pub collider: Box<dyn Collider>,
    pub local_position: glm::Vec3,
    pub local_rotation: glm::Quat,
}
```

This represents **one collision shape attached to a rigid body**. `Box<dyn Collider>` enables polymorphism: a rigid body can contain boxes, spheres, and convex hulls without knowing their concrete types. `ColliderPart::new` initializes its local position to zero and local rotation to the identity quaternion.

### `extract_vertex_points(shape)`

Copies `shape.mesh.vertices`, lexicographically sorts `(x,y,z)`, and removes adjacent duplicates. Duplicate triangle vertices do not help support mapping. Sorting costs $O(V\log V)$.

### `thicken_if_flat(points)`

Measures the Z extent:

$$
\Delta z=z_{\max}-z_{\min}
$$

If $\Delta z\le10^{-4}$, each vertex $(x,y,z)$ is duplicated as:

$$
\mathbf p_+=(x,y,z+h),\quad\mathbf p_-=(x,y,z-h),\quad h=0.05
$$

This gives an originally planar XY-oriented shape a nominal thickness of $0.1$ world/mesh units along local Z, reducing zero-volume degeneracy in 3D GJK/EPA. It does **not** thicken a planar mesh whose plane is tilted relative to local XY.

### `make_collider_parts(shape)`

The factory branches on `shape.shape_type`:

| Shape type | Result | Mathematics / assumption |
|---|---|---|
| `Circle` | `SphereCollider` | $r=\max(|s_x|,|s_y|,|s_z|)$; approximates a circle as a 3D sphere |
| `Cube` | `BoxCollider` | $\mathbf h=|\mathbf s|h_{mesh}$, with $h_{mesh}=1$ |
| Other types | `ConvexHullCollider` | Scale each vertex, optionally thicken, use all vertices as support candidates |

**Nonuniform scaling of mesh points:**

$$
\mathbf p'=\\mathrm{diag}(s_x,s_y,s_z)\mathbf p
$$

**Important:** This code bakes `scale` into collider geometry. The body transform should not apply the same scale again.

---

## 4. Local-to-world support — `world_collider.rs`

### `struct WorldCollider<'a>`

```rust
pub struct WorldCollider<'a> {
    pub collider: &'a dyn Collider,
    pub position: glm::Vec3,
    pub rotation: glm::Quat,
}
```

`WorldCollider` **borrows** collider geometry and pairs it with a world translation and orientation. It does not own or duplicate the collider.

### `WorldCollider::support(world_direction)`

The local collider understands only local-space directions. To query it in world space:

1. Rotate the direction into local space.
2. Call the collider's local support function.
3. Rotate the resulting point back into world space.
4. Add world translation.

$$
\boxed{S_{world}(\mathbf d)=\mathbf t+R\,S_{local}(R^T\mathbf d)}
$$

For unit quaternions, $R^{-1}=R^T$. Quaternion vector rotation is conceptually $q(0,\mathbf v)q^{-1}$.

### `to_world(body, part)` in `physics_world.rs`

Builds a `WorldCollider` from the body's transform and a collider part's local transform:

$$
\mathbf t_{world}=\mathbf x_{body}+R_{body}\mathbf t_{part}
$$

$$
q_{world}=q_{body}q_{part}
$$

This is why `ColliderPart.local_position` matters for compound rigid bodies: its offset rotates with the parent body.

---

## 5. GJK collision detection — `convex_collision.rs`

**Purpose:** Determine whether two convex sets intersect, using only support functions.

### Minkowski difference

$$
A-B=\{\mathbf a-\mathbf b\mid\mathbf a\in A,\mathbf b\in B\}
$$

The two shapes intersect if and only if the origin belongs to their Minkowski difference:

$$
A\cap B\ne\varnothing\iff\mathbf0\in A-B
$$

**Minkowski support:**

$$
S_{A-B}(\mathbf d)=S_A(\mathbf d)-S_B(-\mathbf d)
$$

The implementation's `minkowski(a,b,d)` also stores the individual witness points on A and B, not just their difference. Those witnesses are needed later by EPA to reconstruct a contact point.

### `Vertex`

A support vertex stores:

- `minkowski_different_point`: $\mathbf w=\mathbf a-\mathbf b$.
- `single_point_shape_a`: witness $\mathbf a$.
- `single_point_shape_b`: witness $\mathbf b$.

### `Simplex`

A simplex is the current set of support vertices being tested against the origin:

| Vertices | Simplex | Geometric test |
|---|---|---|
| 1 | Point | Search toward origin |
| 2 | Line segment | Determine which side/region contains the origin |
| 3 | Triangle | Check triangle edge regions and face direction |
| 4 | Tetrahedron | Determine whether the origin is inside all outward faces |

### GJK search loop

1. Choose an initial search direction (the code uses the displacement between collider positions).
2. Add the Minkowski support point in that direction.
3. Update the simplex and the direction toward the origin.
4. If a new support point fails to advance toward the origin, report no collision.
5. If a tetrahedron contains the origin, report intersection and hand the simplex to EPA.

The separating-direction test is based on:

$$
S_{A-B}(\mathbf d)\cdot\mathbf d\le0
$$

When this holds in the tested direction, the support point does not pass the origin-facing plane, so the origin cannot be enclosed by the convex Minkowski difference along that direction.

### `line_case`

For simplex vertices $A,B$:

$$
\mathbf{AB}=B-A,\qquad\mathbf{AO}=-A
$$

When the origin lies in the direction of segment AB, a new search direction perpendicular to AB and toward the origin can be formed using a vector triple product:

$$
\boxed{\mathbf d=(\mathbf{AB}\times\mathbf{AO})\times\mathbf{AB}}
$$

By the triple-product identity:

$$
(\mathbf a\times\mathbf b)\times\mathbf a=\|\mathbf a\|^2\mathbf b-(\mathbf a\cdot\mathbf b)\mathbf a
$$

If the origin is not in the segment's forward region, discard B and search toward $-A$.

### `triangle_case`

For triangle $A,B,C$:

$$
\mathbf{AB}=B-A,\quad\mathbf{AC}=C-A,\quad\mathbf n=\mathbf{AB}\times\mathbf{AC}
$$

The function checks whether the origin lies beyond the triangle's edges or above/below its plane. It either reduces the simplex to an edge or retains the triangle and chooses a search direction parallel to $\pm\mathbf n$.

### `tetrahedron_case`

For $A,B,C,D$, test whether the origin lies outside any tetrahedron face. If outside a face, reduce the simplex to that face and continue. If inside all four consistently oriented face half-spaces, the simplex encloses the origin and GJK succeeds.

**Implementation caveat:** Tetrahedron face normals must be oriented relative to their **opposite vertices**. The earlier implementation's normal-orientation logic deserves verification with unit tests; incorrect winding can produce false collision results.

---

## 6. EPA penetration estimation

**Purpose:** Given a GJK simplex enclosing the origin, estimate how deeply the convex shapes overlap and in what direction they should be separated.

### `Face`

Each EPA polytope face contains:

- Three indices into the support-vertex array.
- An outward unit normal $\mathbf n$.
- Its signed distance $d$ from the origin.

For triangle vertices $A,B,C$:

$$
\mathbf n=\frac{(B-A)\times(C-A)}{\|(B-A)\times(C-A)\|}
$$

$$
d=\mathbf n\cdot A
$$

The implementation flips the normal/winding if needed to make the face distance nonnegative.

### EPA expansion loop

1. Build a tetrahedral polytope from GJK's simplex.
2. Find the face closest to the origin (minimum $d$).
3. Query Minkowski support along that face's outward normal.
4. Compare the support distance $d_{new}$ with the face distance $d_{face}$.
5. If the difference is below tolerance (previous implementation: $10^{-3}$), treat the face as converged.
6. Otherwise, remove faces visible from the new support point, determine the horizon boundary, and connect the new point to the horizon edges.
7. Repeat until convergence or the iteration limit.

**Convergence test:**

$$
S_{A-B}(\mathbf n)\cdot\mathbf n-d_{face}<\varepsilon
$$

The converged face distance estimates penetration depth:

$$
\delta\approx d_{face}
$$

### Barycentric contact reconstruction

For the closest face with vertices $V_0,V_1,V_2$, compute barycentric weights $\lambda_i$ of the closest point to the origin:

$$
\lambda_0+\lambda_1+\lambda_2=1
$$

$$
\mathbf p_A=\lambda_0\mathbf a_0+\lambda_1\mathbf a_1+\lambda_2\mathbf a_2
$$

$$
\mathbf p_B=\lambda_0\mathbf b_0+\lambda_1\mathbf b_1+\lambda_2\mathbf b_2
$$

These are corresponding witness points on the original shapes. The previously shared code used the A witness for `Contact.point`; a midpoint $(\mathbf p_A+\mathbf p_B)/2$ is another common approximation for a single solver contact.

### `Contact`

```rust
pub struct Contact {
    pub point: glm::Vec3,
    pub normal: glm::Vec3,
    pub depth: f32,
}
```

- `point`: world-space location where the solver applies impulses.
- `normal`: collision direction; **the solver assumes it points from body A toward body B**.
- `depth`: overlap distance used for positional correction.

**Critical convention to verify:** EPA's face normal and the sign of `Contact.normal` must agree with the solver's A→B convention. With Minkowski difference **A−B**, the outward normal of the nearest EPA face generally gives the direction for separating A toward the opposite side from B in the simple sphere test; test the sign explicitly before trusting impulse resolution. An inverted normal can cause the solver to skip approaching contacts and push bodies further into each other.

---

## 7. Rigid-body state and mass properties — `rigid_body.rs`

### `struct RigidBody`

| Field | Meaning | Units / representation |
|---|---|---|
| `position` | Body reference position (assumed center of mass) | Vec3, meters |
| `orentation` | Body orientation | Unit quaternion |
| `velocity` | Linear velocity | m/s |
| `angular_velocity` | Angular velocity | rad/s |
| `force` | Accumulated external force | N |
| `volume` | Stored volume; not used by shown mass-property calculation | m³ if SI |
| `drag_coefficience` | Aerodynamic drag coefficient | Dimensionless |
| `cross_section` | Drag reference area | m² |
| `inv_mass` | $1/m$; zero for static bodies | kg⁻¹ |
| `restitution` | Normal collision elasticity | Dimensionless |
| `inv_inertia_local` | Inverse inertia tensor in body-local axes | (kg·m²)⁻¹ |
| `parts` | Attached collider geometry | `Vec<ColliderPart>` |

### `RigidBody::is_static()`

Checks `inv_mass == 0.0`. An immovable body has effectively infinite mass in the impulse equations:

$$
m^{-1}=0\implies\Delta\mathbf v=\mathbf Jm^{-1}=\mathbf0
$$

### `RigidBody::compute_mass_properties(mass)`

1. Set inverse mass $m^{-1}=1/m$ if positive; otherwise zero.
2. If static or without collider parts, leave inverse inertia zero.
3. Divide mass equally among $N$ parts: $m_i=m/N$.
4. For each part, calculate its local inertia $I_i$.
5. Rotate inertia into the body's local axes.
6. Apply the parallel-axis theorem for its offset.
7. Sum all part contributions and invert the result.

**Rotated inertia:**

$$
I_{rotated,i}=R_i I_iR_i^T
$$

**Parallel-axis theorem:**

$$
I_{offset,i}=m_i\left[(\mathbf p_i\cdot\mathbf p_i)\mathbf1_3-\mathbf p_i\mathbf p_i^T\right]
$$

**Total local inertia:**

$$
\boxed{I_{body}=\sum_i\left(R_iI_iR_i^T+m_i[(\mathbf p_i\cdot\mathbf p_i)\mathbf1_3-\mathbf p_i\mathbf p_i^T]\right)}
$$

Then:

$$
I_{body}^{-1}=\\mathrm{inverse}(I_{body})
$$

The current code uses `try_inverse().unwrap_or_else(glm::Mat3::zeros)` to avoid a failed matrix inversion.

**Assumptions:** Mass is distributed **equally per part**, not by collider volume; each part's local inertia is assumed to be about its own center of mass; and the rigid-body reference origin is assumed to be the aggregate center of mass. If these assumptions are violated, angular dynamics will be approximate.

### `RigidBody::inv_intertia_world()`

Transforms inverse inertia from body coordinates to world coordinates:

$$
\boxed{I_{world}^{-1}=R_{body}I_{local}^{-1}R_{body}^T}
$$

This is necessary because impulses and contact offsets are expressed in **world coordinates**, while the inertia tensor is stored in **body-local coordinates**.

---

## 8. Simulation world — `physics_world.rs`

### `BodyId`

```rust
pub enum BodyId { Body(usize), Ground }
```

`Body(usize)` indexes `PhysicsWorld.bodies`; `Ground` refers to a separately stored static rigid body. These indices remain valid only if the body vector is not reordered or removals are handled carefully.

### `Manifold`

```rust
struct Manifold {
    a: BodyId,
    b: BodyId,
    contact: Contact,
}
```

Each `Manifold` associates **one contact** with two bodies. Despite the name, this is not a multi-point contact manifold; the code can create several `Manifold` entries for a single body pair.

### `PhysicsWorld`

```rust
pub struct PhysicsWorld {
    bodies: Vec<RigidBody>,
    ground: RigidBody,
}
```

`new()` sets ground inverse mass and inverse inertia to zero, making it immovable. `add_body()` appends a body and returns its `BodyId`.

### `PhysicsWorld::step(dt)` — actual execution order

```rust
self.integrate_velocities(dt);
let mut contacts = Vec::new();
self.find_body_contacts(&mut contacts);
self.find_ground_contacts(&mut contacts);
self.solve(&contacts);
self.integrate_positions(dt);
```

**Important:** The current implementation detects collisions using the positions/orientations **before** this step's position integration. It updates velocities first, resolves detected contacts, then advances positions. This can delay detection of newly crossed contacts and can permit tunneling at large timesteps.

### `integrate_velocities(dt)`

For each non-static body:

**Speed:**

$$
v=\|\mathbf v\|=\sqrt{v_x^2+v_y^2+v_z^2}
$$

**Quadratic air drag:**

$$
\boxed{\mathbf F_{drag}=-\frac12\rho C_dA\|\mathbf v\|\mathbf v}
$$

The code uses $\rho=1.2$ kg/m³ as an approximate air density.

**Newton's second law:**

$$
\mathbf a=\mathbf g+(\mathbf F_{external}+\mathbf F_{drag})m^{-1}
$$

**Velocity update:**

$$
\boxed{\mathbf v_{new}=\mathbf v_{old}+\mathbf a\Delta t}
$$

Finally `force = Vec3::zeros()` so forces must be accumulated again next step.

**Notably:** The shared implementation does not integrate an accumulated **torque** field or angular acceleration from external torques; angular velocity changes mainly through collision impulses and damping.

### `find_body_contacts(out)`

For each unordered pair of bodies ($N(N-1)/2$ possible pairs):

1. Skip pairs where both bodies are static.
2. Iterate over every `ColliderPart` in each body.
3. Convert each part to `WorldCollider` using `to_world`.
4. Perform a cheap bounding-sphere overlap test.
5. If spheres may overlap, call `collide(&wa,&wb)` (GJK followed by EPA).
6. Store the returned `Contact` as a `Manifold`.

**Bounding-sphere rejection:**

$$
\|\mathbf c_B-\mathbf c_A\|>R_A+R_B\implies\text{no overlap}
$$

Body-pair traversal is $O(N^2)$ before accounting for the number of collider parts and narrowphase work. This is a simple brute-force broadphase, not a spatial tree/grid.

### `find_ground_contacts(out)`

For each non-static body and collider part:

1. Construct its world collider.
2. Query support points in five near-downward directions.
3. If a point lies below `GROUND_HEIGHT = 0`, generate a contact.

**Ground condition and penetration:**

$$
p_z<z_{ground},\qquad\delta=z_{ground}-p_z
$$

The contact normal is $(0,0,-1)$, consistent with the solver's convention because **A is the body and B is the ground**.

**Caveat:** Multiple support directions may return the **same point**, generating duplicate constraints and excessive correction. This ground check is also a plane-specific collision method, not GJK/EPA against a ground collider.

### `get_pair_mut(a,b)` / `get_two_mut(items,a,b)`

These functions are primarily about **Rust borrowing**, not physics. `get_two_mut` uses `split_at_mut` to obtain two non-overlapping mutable references to bodies; `get_pair_mut` additionally handles the separately stored ground.

### `solve(contacts)`

The current solver has two phases:

1. **Position correction:** call `correct_penetration` once per contact.
2. **Velocity resolution:** call `resolve_contact` for every contact in each of `SOLVER_ITERATIONS = 8` passes.

This is a **sequential impulse** approach: after each contact impulse, the bodies' velocities change, and subsequent constraints observe those changes.

---

## 9. Penetration correction — `correct_penetration`

Let $\delta$ be penetration depth, $\mathbf n$ the normal pointing **A→B**, and $m_A^{-1},m_B^{-1}$ inverse masses.

The code uses:

- `SLOP = 0.005`: allowed small overlap.
- `PERCENT = 0.4`: fraction of excess overlap corrected per step.

**Excess penetration:**

$$
\delta'=\max(\delta-\mathrm{SLOP},0)
$$

**Mass-weighted correction:**

$$
\mathbf C=\mathbf n\frac{\delta'}{m_A^{-1}+m_B^{-1}}\mathrm{PERCENT}
$$

**New positions:**

$$
\boxed{\mathbf x_A'=\mathbf x_A-m_A^{-1}\mathbf C}
$$

$$
\boxed{\mathbf x_B'=\mathbf x_B+m_B^{-1}\mathbf C}
$$

The lighter body (larger inverse mass) moves farther. If both inverse masses are zero, no correction is performed. This is **positional** correction: it does not directly change velocity.

---

## 10. Collision impulse — `resolve_contact`

**Purpose:** Change linear and angular velocities so contact points stop moving into each other, with optional bounce.

### Step 1 — contact offsets

For contact point $\mathbf p$:

$$
\mathbf r_A=\mathbf p-\mathbf x_A,\qquad\mathbf r_B=\mathbf p-\mathbf x_B
$$

These lever arms determine the angular effect of a contact impulse.

### Step 2 — velocity at each contact point

A point on a rotating rigid body has velocity:

$$
\boxed{\mathbf v_p=\mathbf v_{center}+\boldsymbol\omega\times\mathbf r}
$$

Thus:

$$
\mathbf v_A^p=\mathbf v_A+\boldsymbol\omega_A\times\mathbf r_A
$$

$$
\mathbf v_B^p=\mathbf v_B+\boldsymbol\omega_B\times\mathbf r_B
$$

### Step 3 — relative normal velocity

$$
\boxed{v_n=(\mathbf v_B^p-\mathbf v_A^p)\cdot\mathbf n}
$$

For an **A→B** normal:

- $v_n<0$: approaching → resolve collision.
- $v_n\ge0$: separating → return without an impulse.

### Step 4 — restitution

$e=\min(e_A,e_B)$. The code sets $e=0$ when $-v_n<1$ m/s to suppress tiny resting bounces.

### Step 5 — world inverse inertia

$$
I_A^{-1}=R_AI_{A,local}^{-1}R_A^T,\qquad I_B^{-1}=R_BI_{B,local}^{-1}R_B^T
$$

### Step 6 — effective inverse mass

$$
\boxed{K=m_A^{-1}+m_B^{-1}+\mathbf n\cdot\left[\left(I_A^{-1}(\mathbf r_A\times\mathbf n)\right)\times\mathbf r_A+\left(I_B^{-1}(\mathbf r_B\times\mathbf n)\right)\times\mathbf r_B\right]}
$$

The first two terms represent translation. The remaining terms represent how the impulse induces **rotation** about each body's center of mass.

### Step 7 — impulse magnitude and direction

$$
\boxed{j=-\frac{(1+e)v_n}{K}},\qquad\boxed{\mathbf J=j\mathbf n}
$$

### Step 8 — linear velocity updates

$$
\boxed{\mathbf v_A'=\mathbf v_A-m_A^{-1}\mathbf J}
$$

$$
\boxed{\mathbf v_B'=\mathbf v_B+m_B^{-1}\mathbf J}
$$

These are equal and opposite impulses (Newton's third law).

### Step 9 — angular velocity updates

$$
\boxed{\boldsymbol\omega_A'=\boldsymbol\omega_A-I_A^{-1}(\mathbf r_A\times\mathbf J)}
$$

$$
\boxed{\boldsymbol\omega_B'=\boldsymbol\omega_B+I_B^{-1}(\mathbf r_B\times\mathbf J)}
$$

**Example:** A cube struck off-center receives a nonzero $\mathbf r\times\mathbf J$, so it both translates and spins. A centered hit with $\mathbf r\parallel\mathbf J$ produces zero angular impulse.

**Current solver limitations:** The code solves only **normal** impulses, not tangential Coulomb friction. It does not cache accumulated impulses across steps (warm starting), nor does it explicitly enforce a nonnegative *accumulated* normal impulse across iterations.

---

## 11. Position and orientation integration — `integrate_positions(dt)`

For each non-static body:

### Damping

```rust
b.velocity *= 0.9995;
b.angular_velocity *= 0.998;
```

These factors are applied **per step**, so their effect depends on simulation frequency.

### Position

$$
\boxed{\mathbf x_{new}=\mathbf x_{old}+\mathbf v_{new}\Delta t}
$$

This is the position half of a semi-implicit Euler-style update: velocity is updated before position.

### Orientation

Angular velocity is represented as a pure quaternion:

$$
\Omega=(0,\omega_x,\omega_y,\omega_z)
$$

With angular velocity in world space:

$$
\frac{dq}{dt}=\frac12\Omega q
$$

Explicit quaternion integration:

$$
\boxed{q_{new}=\\mathrm{normalize}\left(q_{old}+\frac{\Delta t}{2}\Omega q_{old}\right)}
$$

Normalization prevents floating-point drift from breaking the unit-quaternion constraint.

---

## 12. End-to-end example — a falling cube

**Setup:** A cube with half-extents $(1,1,1)$, mass 2 kg, initial center at $z=3$, and no initial velocity.

1. **Collider construction:** `make_collider_parts` creates `BoxCollider { half: (1,1,1) }`.
2. **Mass properties:** `compute_mass_properties(2)` gives $m^{-1}=0.5$ kg⁻¹. Its diagonal inertia is $I_{xx}=I_{yy}=I_{zz}=\frac23(1^2+1^2)=\frac43$ kg·m².
3. **World insertion:** `PhysicsWorld::add_body` stores it and returns `BodyId::Body(index)`.
4. **Gravity integration:** over $\Delta t=1/60$ s, with zero drag/force, $\Delta v_z=-9.81/60\approx-0.1635$ m/s on the first step.
5. **Contact detection:** while the cube's bottom remains above $z=0$, no ground contact is generated.
6. **Ground penetration:** once a downward support point satisfies $p_z<0$, a ground `Contact` is added with downward normal and depth $-p_z$.
7. **Position correction:** the cube is moved upward by a fraction of its penetration (the ground does not move).
8. **Impulse:** the downward contact-point velocity yields $v_n<0$, so the solver applies an upward impulse to the cube, changing its linear velocity and possibly angular velocity if the contact is off-center.
9. **Position/orientation integration:** the cube advances using its corrected velocities.

**Note:** This example illustrates the intended flow; actual resting stability depends on contact uniqueness, contact-point selection, timestep, and solver convergence.

---

## 13. Important assumptions and verification checklist

| Area | Current behavior | Why it matters |
|---|---|---|
| Convex geometry | Arbitrary mesh vertices are treated as a convex hull | Concave meshes are not represented exactly |
| `Circle` collider | Uses a 3D sphere | A flat circle/disk has different collision geometry |
| Mass distribution | Equal mass per collider part | Compound-body inertia is approximate |
| Center of mass | Body origin treated as center of mass | Off-center geometry requires COM handling |
| EPA normal | Must match A→B solver convention | Wrong sign causes incorrect impulses/correction |
| GJK tetrahedron | Face orientation must be consistent | Wrong simplex classification causes false positives/negatives |
| Ground contacts | Five near-downward support queries | May produce duplicate contact constraints |
| Broadphase | All body pairs, bounding-sphere rejection | Scales quadratically with body count |
| Collision timing | Detection precedes position integration | Newly entered collisions can be detected late |
| Time stepping | No fixed-step accumulator shown | Large/variable `dt` can destabilize simulation |
| Contact impulses | Normal impulses only | No surface friction |
| Angular dynamics | No external torque integration shown | Angular acceleration is collision-driven |
| Shape scaling | Baked into collider points/extents | Avoid applying scale twice |

### Suggested focused tests

1. **Two overlapping spheres on X:** verify normal orientation, positive depth, and that correction separates them.
2. **Head-on moving spheres:** verify approaching relative velocity produces a separating impulse.
3. **One cube falling onto ground:** verify stable rest without excessive bouncing or sinking.
4. **Off-center cube impact:** verify induced angular velocity has the expected sign.
5. **Rotated, offset collider part:** verify `to_world` and `WorldCollider::support` agree with transformed geometry.
6. **Two non-overlapping convex shapes:** verify GJK returns no contact.
7. **Nearly coplanar geometry:** verify `thicken_if_flat` avoids degenerate GJK/EPA cases.

---

## 14. Mathematical references

- [Support functions](https://en.wikipedia.org/wiki/Support_function)
- [GJK collision detection](https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm)
- [Convex hulls](https://en.wikipedia.org/wiki/Convex_hull)
- [Rigid-body dynamics](https://en.wikipedia.org/wiki/Rigid_body_dynamics)
- [Moments of inertia](https://en.wikipedia.org/wiki/List_of_moments_of_inertia)
- [Inertia tensor](https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor)
- [Parallel-axis theorem](https://en.wikipedia.org/wiki/Parallel_axis_theorem)
- [Semi-implicit Euler integration](https://en.wikipedia.org/wiki/Semi-implicit_Euler_method)
- [Drag equation](https://en.wikipedia.org/wiki/Drag_equation)
- [Quaternion rotations](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation)
- [Erin Catto, *Sequential Impulses* (GDC 2006)](https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf)
- [Polyhedral mass properties (Geometric Tools)](https://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf)
