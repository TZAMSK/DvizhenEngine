use core::panic;

use crate::physics::{
    collisions::{
        collider_factory::ColliderPart, convex_collision::collide, convex_collision::Contact,
        world_collider::WorldCollider,
    },
    rigid_body::RigidBody,
};

/// Identifies a rigid body in the physics simulation
///
/// Body(index) refers to a body stored in PhysicsWorld::bodies
/// Ground refers to the static ground plane
#[derive(Clone, Copy)]
pub enum BodyId {
    Body(usize),
    Ground,
}

/// Represents a collision constraint between two rigid bodies
///
/// A collision manifold associates
/// - Body A
/// - Body B
/// A contact point, normal and penetrative depth
///
/// The solver uses this info to calculate impulses
/// and correct interprenetration
struct Manifold {
    a: BodyId,
    b: BodyId,
    contact: Contact,
}

// Sequential impulse solver iterations
//
// Each iteration processes every collision constraint
//
// Increasing the iteration count generally improves
// convergence, particularly for stacked objects
//
// Reference:
// https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
const SOLVER_ITERATIONS: u32 = 8;

// Penetration tolerance
//
// Penetrations smaller than SLOP are ignored during
// positional correction to reduce numerical jitter
//
// Reference:
// https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
const SLOP: f32 = 0.005;

// Fraction of penetration corrected during a simulation step
//
// correction = PERCENT * max(depth - SLOP, 0)
//
// Correcting only a fraction can recuce abrupt positional
// changes and oscillation
const PERCENT: f32 = 0.4;

// Gravitational acceleration
// Gravity along the z axis
const GRAVITY: glm::Vec3 = glm::Vec3::new(0.0, 0.0, -9.81);

// Height of the horizontal ground plane
const GROUND_HEIGHT: f32 = 0.0;

// Directions used to sample support points near the
// downward direction
//
// The first direction finds a lowest point
//
// The remaining directions are slightly tilted to sample
// additional nearby extermal points
//
// Each direction is normalized before support mapping
// Reference:
// https://en.wikipedia.org/wiki/Support_function
const DIRS: [glm::Vec3; 5] = [
    glm::Vec3::new(0.0, 0.0, -1.0),
    glm::Vec3::new(0.05, 0.0, -1.0),
    glm::Vec3::new(-0.05, 0.0, -1.0),
    glm::Vec3::new(0.0, 0.05, -1.0),
    glm::Vec3::new(0.0, -0.05, -1.0),
];

/// Stores and simulates rigid bodies
/// 1. Velocity
/// 2. Collision Detection
/// 3. Collision resilution
/// 4. Position and orientation
///
/// Reference:
/// https://en.wikipedia.org/wiki/Rigid_body_dynamics
pub struct PhysicsWorld {
    bodies: Vec<RigidBody>,
    ground: RigidBody,
}

impl PhysicsWorld {
    pub fn new() -> Self {
        let mut ground = RigidBody::default();
        ground.inv_mass = 0.0;
        ground.inv_inertia_local = glm::Mat3::zeros();

        Self {
            bodies: Vec::new(),
            ground,
        }
    }

    pub fn add_body(&mut self, rigid_body: RigidBody) -> BodyId {
        let id = self.bodies.len();
        self.bodies.push(rigid_body);

        BodyId::Body(id)
    }

    /*
    fn remove_body(&mut self, rigid_body: RigidBody) {
        self.bodies.remove(index)
    }
    */

    pub fn step(&mut self, dt: f32) {
        self.integrate_velocities(dt);

        let mut contacts = Vec::new();

        self.find_body_contacts(&mut contacts);
        self.find_ground_contacts(&mut contacts);

        self.solve(&contacts);

        self.integrate_positions(dt);
    }

    fn integrate_velocities(&mut self, dt: f32) {
        for b in &mut self.bodies {
            if b.is_static() {
                continue;
            }

            // Speed:
            // |v| = sqrt(vx² + vy² + vz²)
            // Magnitude
            let speed = glm::length(&b.velocity);

            // Quadratic aerodynamic drag
            //
            // F_drag = -(1/2) * rho * Cd * A * |v| * v
            //
            // Where:
            // rho = fluid density (1.2 kg/m³ for approximate air)
            // |v| * v = produces a force whose magnitude
            // is proportional to speed squared
            //
            // Reference:
            // https://en.wikipedia.org/wiki/Drag_equation
            let drag = -0.5 * 1.2 * b.drag_coefficience * b.cross_section * speed * b.velocity;

            // Newton's second law:
            //
            // F = ma
            // a = F / m
            // using inverse = a = F / inv_mass
            // total acc = g + (F_extrnal + F_drag) * m_inv
            //
            // Semi implicit update:
            // v(t + dt) = v(t) + a * dt
            //
            // Reference:
            // https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
            b.velocity += (GRAVITY + (b.force + drag) * b.inv_mass) * dt;

            // Clear accumulated forces after integration
            //
            // Forces must be accumulated  again for the
            // next simulation step
            b.force = glm::Vec3::zeros();
        }
    }

    fn integrate_positions(&mut self, dt: f32) {
        for b in &mut self.bodies {
            if b.is_static() {
                continue;
            }

            // Numerical damping;
            //
            // v_new = 0.9995 * v
            // w_new = 0.998 * w
            //
            // Slightly reduces linear and angular velocies
            // to suppress persisten numerical oscillations
            //
            // These are per-step factors rather than
            // time-independant damping coefficients
            b.velocity *= 0.9995;
            b.angular_velocity *= 0.998;

            // Position intergration
            // dx / dt = v
            //
            // Discrete approx
            // x(t + dt) = x(t) + v(t + dt) * dt
            //
            // This uses the updated velocity
            // Reference:
            // https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
            b.position += b.velocity * dt;

            // Angulat velocity
            // w = (wx, wy, wz)
            // Represents the body's rotational velocity
            // in radians per second, in world coords
            let w = b.angular_velocity;

            // Pure quaterion representing angular velocity;
            //
            // Omega = (0, wx, wy, wz)
            //
            // The scalar component is zero because angular
            // velocity is represented by a 3d vector
            let spin = glm::Quat::new(0.0, w.x, w.y, w.z);

            // Quaterion orientation detrivative
            //
            // dq / dt = (1/2) * Omega * q
            //
            // Discrete intergation
            // q_new = q + (dt/2) * Omega * q
            //
            // Normalize:
            //
            // q_normalized = q_new / ||q__new||
            //
            // Normalization prevents numerical drift from
            // producing a non-unit orientation quaterion
            //
            // The multiplication order assumes angular
            // velocity is expressed in world coords
            //
            // Reference:
            // https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
            b.orentation = glm::quat_normalize(&(b.orentation + 0.5 * dt * (spin * b.orentation)));
        }
    }

    fn find_body_contacts(&self, out: &mut Vec<Manifold>) {
        // Test every unordered pair of rigid bodies
        //
        // For N bodies, the number of possible pairs is
        //
        // N(N - 1) / 2
        //
        // Therefore, this broadphase hasO (N²)
        // body-pair complexity
        for i in 0..self.bodies.len() {
            for j in (i + 1)..self.bodies.len() {
                let a = &self.bodies[i];
                let b = &self.bodies[j];

                if a.is_static() && b.is_static() {
                    continue;
                }

                // A rigid body can have multiple colliders
                //
                // Test every collider pair belonging to
                // the two rigid bodies
                for pa in &a.parts {
                    for pb in &b.parts {
                        let wa = to_world(a, pa);
                        let wb = to_world(b, pb);

                        // Bounding sphere overlaps test
                        // Reference:
                        // https://en.wikipedia.org/wiki/Bounding_sphere
                        let reach = wa.collider.bounding_radius() + wb.collider.bounding_radius();

                        if glm::distance(&wa.position, &wb.position) > reach {
                            continue;
                        }

                        // Narrowphase collision detection
                        //
                        // GJK determines whether the convex
                        // colliders intersect
                        //
                        // EPA calculates penetration depth,
                        // collision normal and contact point
                        //
                        // Reference:
                        // https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm
                        if let Some(contact) = collide(&wa, &wb) {
                            out.push(Manifold {
                                a: BodyId::Body(i),
                                b: BodyId::Body(j),
                                contact,
                            });
                        }
                    }
                }
            }
        }
    }

    fn find_ground_contacts(&self, out: &mut Vec<Manifold>) {
        for (i, b) in self.bodies.iter().enumerate() {
            if b.is_static() {
                continue;
            }

            for part in &b.parts {
                let w = to_world(&b, &part);
                for dir in DIRS {
                    // Support mapping
                    //
                    // S(d) = argmax(p · d), p in C
                    //
                    // A downward search direction selects
                    // an extermal point near the bottom
                    // of the collider
                    //
                    // d_normalized = d / ||d||
                    //
                    // Reference:
                    // https://en.wikipedia.org/wiki/Support_function
                    let p = w.support(&glm::normalize(&dir));

                    // Ground penetration
                    if p.z < GROUND_HEIGHT {
                        out.push(Manifold {
                            a: BodyId::Body(i),
                            b: BodyId::Ground,
                            contact: Contact {
                                point: p,
                                normal: glm::Vec3::new(0.0, 0.0, -1.0),
                                depth: GROUND_HEIGHT - p.z,
                            },
                        });
                    }
                }
            }
        }
    }

    fn get_pair_mut(&mut self, a: BodyId, b: BodyId) -> (&mut RigidBody, &mut RigidBody) {
        // Retrieve two distinct mutable rigid-body references
        //
        // Rust requires mutable references to refer to
        // non-overlapping memory
        //
        // The ground is stored separately from dynamic bodies
        match (a, b) {
            (BodyId::Body(i), BodyId::Body(j)) => get_two_mut(&mut self.bodies, i, j),

            (BodyId::Body(i), BodyId::Ground) => (&mut self.bodies[i], &mut self.ground),

            (BodyId::Ground, BodyId::Body(i)) => (&mut self.ground, &mut self.bodies[i]),

            (BodyId::Ground, BodyId::Ground) => {
                panic!("Ground cannot collide with itself")
            }
        }
    }

    fn solve(&mut self, contacts: &[Manifold]) {
        // Position correction
        //
        // Move penetrating bodies apart to reduce
        // overlap between colliders
        //
        // This is seperate from velocity impulses
        for m in contacts {
            let (a, b) = self.get_pair_mut(m.a, m.b);

            correct_penetration(a, b, &m.contact);
        }

        // Sequential impulse solver
        //
        // Process contacts repeatedly, updating velocities
        // after each impulse
        //
        // Later contacts see the velocity changes produced
        // by earlier
        //
        // This is an iterative approximation to solving
        // multiple coupled collision constraints
        //
        // Reference:
        // https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
        for _ in 0..SOLVER_ITERATIONS {
            for m in contacts {
                let (a, b) = self.get_pair_mut(m.a, m.b);
                resolve_contact(a, b, &m.contact);
            }
        }
    }
}

/// Returns two distinct mutable references from a slice
///
/// split_at_mut(k) divides the slice into
/// left = items[0..k]
/// right = items[k..]
///
/// The two slices cannot overlap
fn get_two_mut<T>(items: &mut [T], a: usize, b: usize) -> (&mut T, &mut T) {
    assert_ne!(a, b);

    if a < b {
        let (left, right) = items.split_at_mut(b);
        (&mut left[a], &mut right[0])
    } else {
        let (left, right) = items.split_at_mut(a);
        (&mut right[0], &mut left[b])
    }
}

// Converts a collider from rigid-body local space to world space
//
// World pos:
// p_World = p_body + R_body * p_local
//
// World orientation:
// q_world = q_body * q_local
//
// Reference:
// https://en.wikipedia.org/wiki/Rigid_transformation
fn to_world<'a>(rigid_body: &'a RigidBody, collider_part: &'a ColliderPart) -> WorldCollider<'a> {
    WorldCollider {
        collider: collider_part.collider.as_ref(),

        // Rotate the collider's local offset by the
        // rigid body's orientation, the translate it
        // by the rigid body's world position
        position: rigid_body.position
            + glm::quat_rotate_vec3(&rigid_body.orentation, &collider_part.local_position),

        // Quaterion composition
        // q_world = q_body * q_part
        //
        // Apply the collider's local orientation followed
        // by the rigid body's orientation
        rotation: rigid_body.orentation * collider_part.local_rotation,
    }
}

/// Resolves a collision using normal impulse
///
/// The impulse changes linear and angular velocities
/// so that the contact points stop approaching each other
///
/// Reference:
/// https://en.wikipedia.org/wiki/Collision_response
/// https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
fn resolve_contact(rigid_body_a: &mut RigidBody, rigid_body_b: &mut RigidBody, contact: &Contact) {
    // Collision normal
    //
    // ||n|| = 1
    //
    // Convention: points from A toward B
    let n = contact.normal;

    // Collision offsets
    //
    // r_A = p_contact - x_A
    // r_B = p_contact - x_B
    //
    // These vectors connect each body's center of mass
    // to the world-space contact point
    let (ra, rb) = (
        contact.point - rigid_body_a.position,
        contact.point - rigid_body_b.position,
    );

    // Velocity of a point
    //
    // v_point = v_center + w x r
    //
    // Where:
    // v_center = linear velocity
    // w        = angular velocity
    // r        = offset from center of pass
    //
    // The cross product accounts for the velocity
    // produced by rotation at the contact point
    //
    // Reference:
    // https://en.wikipedia.org/wiki/Rigid_body_dynamics
    let (va, vb) = (
        rigid_body_a.velocity + glm::cross(&rigid_body_a.angular_velocity, &ra),
        rigid_body_b.velocity + glm::cross(&rigid_body_b.angular_velocity, &rb),
    );

    // Relative normal velocity:
    //
    // v_rel = v_B - v_A
    // vn = v_rel · n
    //
    // vn < 0: contact points approaching
    // vn > 0: contact points separating
    //
    // Only approaching contacts require a
    // collision impulse
    let vn = glm::dot(&(vb - va), &n);

    if vn >= 0.0 {
        return;
    }

    // Coefficient of restitution
    //
    // e in [0, 1]
    //
    // e = 0: perfectly inelastic normal collision
    // e = 1: perfectly elastic normal collision
    //
    // Use the smaller restitution coefficient
    // of the two bodies
    //
    // Use the smaller restitution coefficient
    // of the bodies
    //
    // Reference:
    // https://en.wikipedia.org/wiki/Coefficient_of_restitution
    let mut e = rigid_body_a.restitution.min(rigid_body_b.restitution);

    // Disable restitution at low impact speeds
    //
    // If |vn| < 1 m/s, treat the impact as inelastic
    //
    // This reduces small persistent bounces
    // when objects settle on surfaces
    if -vn < 1.0 {
        e = 0.0
    }

    // World space inverse intertia tensors:
    //
    // I_world_inv = R * I_local_inv * R^T
    //
    // These tensors determine how angular velocity
    // responds to torque or angular impulse
    let (inv_ia, inv_ib) = (
        rigid_body_a.inv_intertia_world(),
        rigid_body_b.inv_intertia_world(),
    );

    // Effective inverse mass along the contact normal:
    //
    // K = mA_inv + mB_inv +
    //   + n · [
    //       (IA_inv * (rA x n)) x rA
    //       +
    //       (IB_inv * (rB x n)) x rB
    //   ]
    //
    // The first two terms describe linear motion
    //
    // The cross-product terms describe rotation
    // motion included by an impulse at the contact point
    //
    // A contact farther from the center of mass can
    // produce a larger angular response
    //
    // Reference:
    // https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
    let denom = rigid_body_a.inv_mass
        + rigid_body_b.inv_mass
        + glm::dot(
            &n,
            &(glm::cross(&(inv_ia * glm::cross(&ra, &n)), &ra)
                + glm::cross(&(inv_ib * glm::cross(&rb, &n)), &rb)),
        );

    // A nonpositive effective inverse mass cannot
    // produce a valid normal impulse
    if denom <= 0.0 {
        return;
    }

    // Scalar impulse magnitude
    //
    // j = -(1 + e) * vn / K
    //
    // Vector impulse
    //
    // J = j * n
    //
    // The impulse is directed along the collision normal
    //
    // The (1 + e) factor accounts for restitution
    //
    // Reference:
    // https://en.wikipedia.org/wiki/Impulse_(physics)
    let impulse_magnitude: f32 = -(1.0 + e) * vn / denom;
    let impulse: glm::Vec3 = n * impulse_magnitude;

    // Linear impulse-momentum equation:
    // Delta_v = J / m = J * m_inv
    //
    // Angular impulse-momentum equation:
    // Delta_w = I_inv * (r x J)
    //
    // Newton's third law gives equal and opposite
    // impulses on the two bodies
    //
    // Body A receives -J
    rigid_body_a.velocity -= impulse * rigid_body_a.inv_mass;
    rigid_body_a.angular_velocity -= inv_ia * glm::cross(&ra, &impulse);

    // Body B receives +J
    rigid_body_b.velocity += impulse * rigid_body_b.inv_mass;
    rigid_body_b.angular_velocity += inv_ib * glm::cross(&rb, &impulse);
}

/// Corrects positional interprenetration between 2 bodies
///
/// The correctoin is distributed according to inverse mass
///
/// Lighter bodies move more than heavier bodies
///
/// Reference:
/// https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
fn correct_penetration(
    rigid_body_a: &mut RigidBody,
    rigid_body_b: &mut RigidBody,
    contact: &Contact,
) {
    // Combined inverse mass
    //
    // M = mA_inv + mB_inv
    //
    // If M = 0, both bodies are static
    let inv_sum = rigid_body_a.inv_mass + rigid_body_b.inv_mass;
    if inv_sum == 0.0 {
        return;
    }

    // Corrected penetration
    //
    // d = max(depth - SLOP, 0)
    //
    // Correction vector
    //
    // C = n * (d / N) * PERCENT
    //
    // SLOP allows a small amount of penetration
    //
    // PERCENT controls the fraction corrected
    // during this simulation step
    let correction = contact.normal * (contact.depth - SLOP).max(0.0) / inv_sum * PERCENT;

    // Distribute correction by inverse mass
    //
    // xA_new = xA - C * mA_inv
    // xB_new = xB + C * mB_inv
    //
    // A static body has m_inv = 0 and does not move
    rigid_body_a.position -= correction * rigid_body_a.inv_mass;
    rigid_body_b.position += correction * rigid_body_b.inv_mass;
}
