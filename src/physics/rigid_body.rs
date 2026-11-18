use crate::physics::collisions::collider_factory::ColliderPart;

/// Represents the physical state of a rigid body
///
/// A rigid body is an an object whose internal distances
/// remain constant during motion
///
/// Its motion consists of:
/// - Translation of its center of mass
/// - Rotation around its center of mass
pub struct RigidBody {
    pub position: glm::Vec3,
    pub orentation: glm::Quat,
    // Linear velocity
    // v = dx / dt
    // m/s
    pub velocity: glm::Vec3,
    // Angular velocity
    // w = (wx, wy, wz)
    // rad/s
    pub angular_velocity: glm::Vec3,
    // Accumulated external forces
    // F_total = sum(F_i)
    // N
    // Cleared after each velocity intergration
    pub force: glm::Vec3,
    pub volume: f32,
    // Aerodynamics
    pub drag_coefficience: f32,
    // Reference cross-sectional area used for drag
    // m²
    // Larger area produces more drag
    // for the same velocity and drag_coefficience
    pub cross_section: f32,
    // Inverse mass
    // 1 / mass
    // Avoids repeated division
    // during force and impulse calculations
    pub inv_mass: f32,
    // e in [0, 1]
    // Determines the normal velocity retained
    // after a collision
    pub restitution: f32,
    // Inv inertia
    // I_inv = 1 / I = I⁻¹
    // Relats angular impulse to angular velocity
    // Delta_w = I_inv * Delta_L
    pub inv_inertia_local: glm::Mat3,
    // Multiple collider parts allow compound bodies
    // composed several geometric primitives
    pub parts: Vec<ColliderPart>,
}

impl Default for RigidBody {
    fn default() -> Self {
        Self {
            position: glm::Vec3::zeros(),
            orentation: glm::quat_identity(), // q = (1, 0, 0, 0)
            velocity: glm::Vec3::zeros(),
            angular_velocity: glm::Vec3::zeros(),
            force: glm::Vec3::zeros(),
            volume: 0.0,
            drag_coefficience: 0.0,
            cross_section: 0.0,
            inv_mass: 1.0, // Default mass 1 / 1 = 1
            restitution: 0.4,
            inv_inertia_local: glm::Mat3::zeros(),
            parts: vec![],
        }
    }
}

impl RigidBody {
    pub fn is_static(&self) -> bool {
        self.inv_mass == 0.0
    }

    // Calculates the inverse inertia tensor in world space
    //
    // Inertia tensor changes orentation when body rotates
    //
    // Tensor transformation
    // I_world = R * I_local * R^T
    //
    // For inverse tensor
    // I_world_inv = R * I_local_inv * R^T
    //
    // Where:
    // R = rotation matrix
    // R^T = transpose of R
    // Because R is orthogonal
    // R⁻¹ = R^T
    pub fn inv_intertia_world(&self) -> glm::Mat3 {
        // Convert the orientation quaterion to
        // a 3x3 rotation matrix
        let r = glm::quat_to_mat3(&self.orentation);

        // Rotate the local inverse inertia tensor
        // into world coords
        r * self.inv_inertia_local * glm::transpose(&r)
    }

    // Computes inverse mass and inverse inertia
    // from the body's colliders and total mass
    //
    // For compound rigid bodies, the total inertia
    // tensor is the sum of the individual collider
    // inertia tensors, transformed into the body's
    // local coord system
    //
    // References:
    // https://en.wikipedia.org/wiki/Parallel_axis_theorem
    // https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
    pub fn compute_mass_properties(&mut self, mass: f32) {
        self.inv_mass = if mass > 0.0 { 1.0 / mass } else { 0.0 };

        // Reset the inverse inertia tensor before
        // calculating new mass properties
        //
        // Also ensures static bodies have zero inverse inertia
        self.inv_inertia_local = glm::Mat3::zeros();

        if self.is_static() || self.parts.is_empty() {
            return;
        }

        // Equal mass distribution approx
        //
        // m_part = m_total / new
        //
        // Assumes every collider part contributes
        // an equal amount of mass, regardless of volume
        let part_mass = mass / self.parts.len() as f32;

        //      | 1 0 0 |
        // E  = | 0 1 0 |
        //      | 0 0 1 |
        //
        // Used in parallel axis theroem
        let identity = glm::Mat3::identity();

        // total inertia tensor
        let mut inertia = glm::Mat3::zeros();

        for part in &self.parts {
            // Local rotation matrix of a collider
            let r = glm::quat_to_mat3(&part.local_rotation);

            // Collider intertia tensor about its own
            // local origin
            //
            // Different collider geometries use
            // different intertia equations
            let local_inertia = part.collider.intertia(part_mass);

            // Rotate the collider inertia tensor into
            // the rigid body's local coord system
            //
            // Accounts for collider orientation
            // relative to the rigid body
            let rotated_intertia = r * local_inertia * glm::transpose(&r);

            // Collider offset from the rigid body's
            // local reference point
            //
            // p = (px, py, pz)
            let p = part.local_position;

            // Squared distance from the origin
            //
            // d² = p · p
            // d² = px² + py² + pz²
            let distance_squared = glm::dot(&p, &p);

            // Outer product
            //
            // P = p * p^T
            //
            //     | px²   pxpy  pxpz |
            // P = | pypx  py²   pypz |
            //     | pzpx  pzpy  pz²  |
            //
            // Unlike the dot product, the outer
            // product produces a 3x3 matrix
            //
            // Refence:
            // https://en.wikipedia.org/wiki/Outer_product
            let outer = p * glm::transpose(&p);

            // Parallel axis theorem in tensor form
            //
            // I_offset = m * (d² * E - p * p^T)
            //
            // Where:
            // m = collider mass
            // d = distance from ref origin
            // E = identity matrix
            // p = collider offset vector
            //
            // Accounts for additional rotational inertia
            // cause by placing mass away from the rigid
            // body's ref point
            //
            // Reference:
            // https://en.wikipedia.org/wiki/Parallel_axis_theorem
            let offset_inertia = part_mass * (distance_squared * identity - outer);

            // Total inertia
            inertia += rotated_intertia + offset_inertia;
        }

        // Invert total intertia tensor
        // I_inv = I_total⁻¹
        //
        // Inverse intertia tensor is used to
        // calculate changes in angular velocity
        //
        // Delta_w = I_inv * Delta_L
        //
        // A singular tensor cannot be converted
        //
        // This implementation uses direct inversion
        // and assumes the computed tensor is invertible
        self.inv_inertia_local = inertia.try_inverse().unwrap_or_else(glm::Mat3::zeros)
    }
}
