# Physics Pipeline — Visual Graphs (`sliyan_engine`)

This companion to `physics_pipeline.md` visualizes the **actual code architecture**, **function calls**, **data flow**, and **mathematics** in the Rust rigid-body engine. The diagrams use [Mermaid](https://mermaid.js.org/), which renders in compatible Markdown viewers (including GitHub).

> File names follow the proposed module layout (`collider_builder.rs`, `world_collider.rs`, `convex_collision.rs`); identifiers such as `orentation` and `intertia` retain their current spelling. A diagram arrow means *data flow* or *a function call*, as indicated in each section; it does not imply that every stage happens automatically.

## 1. Full pipeline — geometry to simulation

```mermaid
flowchart TD
    A["Scene Shape<br/>mesh.vertices · shape_type · transform.scale"] --> B["make_collider_parts(shape)"]
    B --> C{"ShapeType?"}
    C -->|Circle| D["SphereCollider<br/>radius = max(abs(scale))"]
    C -->|Cube| E["BoxCollider<br/>half = abs(scale) × mesh half extent"]
    C -->|Other| F["extract_vertex_points<br/>sort · dedup"]
    F --> G["Apply componentwise scale"]
    G --> H["thicken_if_flat"]
    H --> I["ConvexHullCollider"]
    D --> J["Vec of ColliderPart"]
    E --> J
    I --> J
    J --> K["RigidBody.parts"]
    K --> L["compute_mass_properties(mass)<br/>inv_mass · inv_inertia_local"]
    L --> M["PhysicsWorld.add_body(body)"]
    M --> N["PhysicsWorld.step(dt)"]
    N --> O["Velocity integration"]
    O --> P["Detect body and ground contacts"]
    P --> Q["Resolve contacts"]
    Q --> R["Position and quaternion integration"]
```

**Interpretation:** `make_collider_parts` returns colliders; the caller must attach them to a `RigidBody`, compute its mass properties, and add the body to `PhysicsWorld`.

## 2. Struct ownership and relationships

```mermaid
classDiagram
    class Shape {
        +shape_type
        +mesh.vertices
        +transform.scale
    }
    class Collider {
        <<trait>>
        +support(direction) Vec3
        +intertia(mass) Mat3
        +bounding_radius() f32
    }
    class SphereCollider
    class BoxCollider
    class ConvexHullCollider
    class ColliderPart {
        +collider Box~dyn Collider~
        +local_position Vec3
        +local_rotation Quat
    }
    class WorldCollider {
        +collider &dyn Collider
        +position Vec3
        +rotation Quat
        +support(world_direction) Vec3
    }
    class RigidBody {
        +position Vec3
        +orentation Quat
        +velocity Vec3
        +angular_velocity Vec3
        +force Vec3
        +inv_mass f32
        +inv_inertia_local Mat3
        +parts Vec~ColliderPart~
        +compute_mass_properties(mass)
        +inv_intertia_world() Mat3
    }
    class PhysicsWorld {
        -bodies Vec~RigidBody~
        -ground RigidBody
        +add_body(body) BodyId
        +step(dt)
    }
    class Contact {
        +point Vec3
        +normal Vec3
        +depth f32
    }
    Collider <|.. SphereCollider
    Collider <|.. BoxCollider
    Collider <|.. ConvexHullCollider
    ColliderPart --> Collider : owns trait object
    RigidBody "1" *-- "0..*" ColliderPart
    PhysicsWorld "1" *-- "0..*" RigidBody
    WorldCollider --> Collider : borrows
    Shape ..> ColliderPart : factory input
    PhysicsWorld ..> Contact : detects and resolves
```

`WorldCollider` is a **temporary world-space view**, not a second owner of the collider. `Contact` is the collision result used by the solver.

## 3. Exact order of one simulation step

```mermaid
sequenceDiagram
    participant Caller
    participant World as PhysicsWorld
    participant Bodies as RigidBody collection
    participant Collision as GJK / EPA
    participant Solver as Contact solver
    Caller->>World: step(dt)
    World->>Bodies: integrate_velocities(dt)
    Note over Bodies: v += [g + (F + F_drag)/m] dt
    World->>World: find_body_contacts(out)
    loop each body pair and collider-part pair
        World->>World: to_world(body, part)
        Note over World: bounding sphere rejection
        World->>Collision: collide(worldA, worldB)
        Collision-->>World: Option<Contact>
    end
    World->>World: find_ground_contacts(out)
    Note over World: downward support samples vs z = 0
    World->>Solver: solve(contacts)
    loop each contact (once)
        Solver->>Bodies: correct_penetration(A, B, contact)
    end
    loop 8 solver iterations
        loop each contact
            Solver->>Bodies: resolve_contact(A, B, contact)
        end
    end
    World->>Bodies: integrate_positions(dt)
    Note over Bodies: damping; x += v dt; normalize updated q
    World-->>Caller: return
```

**Important order detail:** contact detection happens **after velocity integration but before position integration** in this implementation. Position correction is applied before the eight impulse passes.

## 4. Body–body collision detection: broadphase → GJK → EPA

```mermaid
flowchart TD
    A["find_body_contacts"] --> B["For each unordered pair of bodies i, j"]
    B --> C{"Both static?"}
    C -->|Yes| X["Skip"]
    C -->|No| D["For each ColliderPart pair"]
    D --> E["to_world(A, partA)<br/>to_world(B, partB)"]
    E --> F{"Distance between centers ≤ rA + rB?"}
    F -->|No| X
    F -->|Yes| G["collide(worldA, worldB)"]
    G --> H["gjk: origin inside Minkowski difference?"]
    H -->|No| X
    H -->|Yes| I["epa: estimate contact geometry"]
    I --> J{"EPA returns contact?"}
    J -->|No| X
    J -->|Yes| K["Manifold { body A, body B, Contact }"]
    K --> L["Append to contacts"]
```

**Broadphase sphere test:** `||xA − xB|| ≤ RA + RB`. This is a rejection test; passing it does not prove intersection.

## 5. GJK — simplex decision loop

```mermaid
flowchart TD
    A["Initial direction d = centerB − centerA"] --> B["First Minkowski support vertex"]
    B --> C["d = −newest vertex"]
    C --> D["New support: S_A(d) − S_B(−d)"]
    D --> E{"support · d ≤ 0?"}
    E -->|Yes| N["No collision"]
    E -->|No| F["Push vertex into simplex"]
    F --> G{"Simplex size"}
    G -->|2| H["line_case<br/>keep segment or point; update d"]
    G -->|3| I["triangle_case<br/>keep face or edge; update d"]
    G -->|4| J["tetrahedron_case<br/>test origin against faces"]
    H --> K["Continue search"]
    I --> K
    J -->|Origin enclosed| Y["Collision detected"]
    J -->|Outside a face| K
    K --> L{"64 iterations exhausted?"}
    L -->|No| D
    L -->|Yes| N
```

**Minkowski support equation:** `S_(A−B)(d) = S_A(d) − S_B(−d)`.

**Line search direction:** `d = (AB × AO) × AB`.

**Triangle face normal:** `n = AB × AC`.

**GJK criterion:** intersection corresponds to `0 ∈ A−B` for intersecting convex sets. The code's numerical stopping conditions can produce false negatives in difficult or degenerate cases.

## 6. EPA — penetration and contact reconstruction

```mermaid
flowchart TD
    A["GJK tetrahedron"] --> B["Build four triangular Face records"]
    B --> C["Find face nearest Minkowski origin"]
    C --> D["Query Minkowski support along face normal"]
    D --> E["delta = n · support(n) − face.distance"]
    E --> F{"delta < 0.001?"}
    F -->|Yes| G["make_contact(face)"]
    F -->|No| H["Remove faces visible from new vertex"]
    H --> I["Collect horizon edges; cancel reversed edges"]
    I --> J{"Horizon empty?"}
    J -->|Yes| Z["EPA failure"]
    J -->|No| K["Add new vertex and connect horizon faces"]
    K --> L{"Iteration limit reached?"}
    L -->|No| C
    L -->|Yes| G
    G --> M["Project origin onto face: q = n × distance"]
    M --> N["Compute barycentric weights w0, w1, w2"]
    N --> O["Reconstruct points pA, pB from stored supports"]
    O --> P["Contact point = (pA + pB) / 2"]
    P --> Q["Contact normal = −face.normal<br/>depth = face.distance"]
```

**Barycentric interpolation:** `q = w0 A + w1 B + w2 C`, with `w0 + w1 + w2 = 1`.

**Contact reconstruction:** `pA = Σ wi · supportA_i`, `pB = Σ wi · supportB_i`, `point = (pA + pB)/2`.

The normal shown is the **current code's convention** (`-face.normal`), not a claim that its sign is correct for every collision configuration. The solver expects a normal from body A toward body B.

## 7. Ground contacts — separate from GJK/EPA

```mermaid
flowchart TD
    A["find_ground_contacts"] --> B["Each non-static body"]
    B --> C["Each ColliderPart"]
    C --> D["to_world(body, part)"]
    D --> E["Five downward and slightly tilted directions"]
    E --> F["p = WorldCollider.support(normalize(dir))"]
    F --> G{"p.z < GROUND_HEIGHT?"}
    G -->|No| H["No contact for sample"]
    G -->|Yes| I["Contact point = p"]
    I --> J["normal = (0, 0, −1)"]
    J --> K["depth = GROUND_HEIGHT − p.z"]
    K --> L["Append body–ground Manifold"]
```

The ground is represented by a **static `RigidBody` in the solver**, but ground contacts are generated analytically from support points and the plane `z = 0`, rather than by running GJK against a ground collider.

## 8. Contact solver — position and velocity

```mermaid
flowchart TD
    A["solve(contacts)"] --> B["Pass 1: correct_penetration for each contact"]
    B --> C["inverse mass sum = inv_mA + inv_mB"]
    C --> D["correction = n × max(depth − slop, 0) × percent / sum"]
    D --> E["xA −= correction × inv_mA<br/>xB += correction × inv_mB"]
    E --> F["Pass 2: 8 sequential impulse iterations"]
    F --> G["Compute contact offsets rA, rB"]
    G --> H["Point velocities v + ω × r"]
    H --> I["Normal relative speed vn = (vB − vA) · n"]
    I --> J{"vn < 0?"}
    J -->|No| K["Skip impulse"]
    J -->|Yes| L["World inverse inertias<br/>R I_local⁻¹ Rᵀ"]
    L --> M["Effective inverse mass K"]
    M --> N["Impulse scalar j = −(1 + e) vn / K"]
    N --> O["Impulse J = j n"]
    O --> P["Update both linear and angular velocities"]
    P --> Q["Next contact / iteration"]
    K --> Q
```

**Contact-point velocity:** `v_contact = v_center + ω × r`.

**Effective inverse mass:**

```text
K = inv_mA + inv_mB
  + n · [ (IA_inv (rA × n)) × rA
        + (IB_inv (rB × n)) × rB ]
```

**Impulse:** `j = −(1 + e) vn / K`, `J = j n`.

**Linear and angular updates:** `Δv = ±inv_m J`, `Δω = ±I_world⁻¹ (r × J)`.

There is no tangential/friction impulse in the current `resolve_contact` function.

## 9. Mass properties and inertia — calculation graph

```mermaid
flowchart TD
    A["RigidBody.compute_mass_properties(mass)"] --> B["inv_mass = 1 / mass, or 0 if static"]
    B --> C{"Static or no parts?"}
    C -->|Yes| X["Leave inverse inertia zero"]
    C -->|No| D["part_mass = total mass / number of parts"]
    D --> E["For each ColliderPart"]
    E --> F["I_local = collider.intertia(part_mass)"]
    F --> G["Rotate part inertia: R I_local Rᵀ"]
    G --> H["Parallel-axis offset: m (||p||² E − p pᵀ)"]
    H --> I["Accumulate both terms into total inertia"]
    I --> J["Inverse of total inertia; zero if singular"]
    J --> K["RigidBody.inv_inertia_local"]
    K --> L["During contact: inv_intertia_world()"]
    L --> M["I_world⁻¹ = R_body I_local⁻¹ R_bodyᵀ"]
```

**Assumptions:** every collider part receives equal mass; collider inertia is evaluated about its local origin; and the body's reference origin is treated as its center of mass. Convex hull inertia is approximated by its axis-aligned bounding box.

## 10. Per-frame dynamics equations

```mermaid
flowchart LR
    A["State: x, q, v, ω"] --> B["External force + aerodynamic drag"]
    B --> C["a = g + (F + F_drag) inv_mass"]
    C --> D["v ← v + a dt"]
    D --> E["Collision impulses modify v and ω"]
    E --> F["Damping: v × 0.9995; ω × 0.998"]
    F --> G["x ← x + v dt"]
    G --> H["q ← normalize(q + 0.5 dt Ωq)"]
    H --> I["New state: x, q, v, ω"]
```

```text
F_drag = −(1/2) ρ Cd A |v| v, with ρ = 1.2
Ω = quaternion(0, ωx, ωy, ωz)
q_next = normalize(q + (dt/2) Ωq)
```

The quaternion multiplication order `Ωq` assumes **world-space angular velocity**. Damping factors are per-step constants, so changing the simulation timestep changes their effective damping rate.

## 11. Code map — where each stage lives

| Stage | Structs / functions | Current source |
|---|---|---|
| Geometry | `Collider`, `SphereCollider`, `BoxCollider`, `ConvexHullCollider` | `collider.rs` |
| Collider creation | `ColliderPart`, `extract_vertex_points`, `thicken_if_flat`, `make_collider_parts` | `collider_builder.rs` |
| World-space support | `WorldCollider::support` | `world_collider.rs` |
| GJK | `minkowski`, `Simplex`, `line_case`, `triangle_case`, `tetrahedron_case`, `gjk` | `convex_collision.rs` |
| EPA | `Face`, `make_face`, `make_contact`, `find_closest_face`, `epa` | `convex_collision.rs` |
| Contact API | `Contact`, `collide` | `convex_collision.rs` |
| Mass properties | `RigidBody::compute_mass_properties`, `inv_intertia_world` | `rigid_body.rs` |
| Simulation | `PhysicsWorld::step`, `integrate_velocities`, `integrate_positions` | `physics_world.rs` |
| Contact generation | `find_body_contacts`, `find_ground_contacts`, `to_world` | `physics_world.rs` |
| Contact resolution | `solve`, `resolve_contact`, `correct_penetration` | `physics_world.rs` |

## 12. Implementation caveats worth keeping in mind

- **GJK tetrahedron face normals:** `tetrahedron_case` uses cross products without checking each normal against the opposite tetrahedron vertex. Incorrect orientation can cause false intersection results.
- **Contact normal convention:** EPA returns `-face.normal`, while the impulse solver expects A → B. Confirm the convention using simple test cases before relying on the resulting impulses.
- **Multiple ground samples:** nearby support directions may generate duplicate contacts, affecting the number of corrections and impulses.
- **No continuous collision detection:** contacts are sampled at discrete steps; fast bodies can tunnel through thin geometry.
- **No friction or persistent contact cache:** the solver handles normal impulses only and rebuilds contacts every step.

## References

- [GJK collision detection](https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm)
- [Minkowski addition and difference](https://en.wikipedia.org/wiki/Minkowski_addition)
- [Support functions](https://en.wikipedia.org/wiki/Support_function)
- [EPA walkthrough](https://winter.dev/articles/epa-algorithm/)
- [Sequential impulses — Erin Catto](https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf)
- [Parallel-axis theorem](https://en.wikipedia.org/wiki/Parallel_axis_theorem)
- [Rigid-body dynamics](https://en.wikipedia.org/wiki/Rigid_body_dynamics)
