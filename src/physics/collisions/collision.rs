use crate::physics::collisions::collider::Collider;

pub struct WorldCollider<'a> {
    pub collider: &'a dyn Collider,
    pub position: glm::Vec3,
    pub rotation: glm::Quat,
}

impl<'a> WorldCollider<'a> {
    pub fn new(collider: &'a dyn Collider) -> Self {
        Self {
            collider,

            // t = (0, 0, 0)
            //
            // Initial world-space translation
            position: glm::Vec3::zeros(),

            // q_identity = (w: 1, x: 0, y: 0, z: 0)
            //
            // Initial world-space rotation
            //
            // Reference:
            // https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
            rotation: glm::quat_identity(),
        }
    }

    /// Returns the collider's support point in the world space
    ///
    /// Complete equation
    ///
    /// S_world(d) = t + R * S_local(R⁻¹ * d)
    ///
    /// Where:
    /// - d       = world-space search direction
    /// - t       = world-space position
    /// - R       = world-space rotation
    /// - S_local = local-space support function
    ///
    /// Reference:
    /// https://en.wikipedia.org/wiki/Support_function
    /// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
    pub fn support(&self, world_diretion: &glm::Vec3) -> glm::Vec3 {
        // q⁻¹ = conjugate(q) / ||q||²
        //
        // For a unit quaterion
        //
        // q⁻¹ = conjugate(q)
        //
        // Calculate the inverse rotation so that a vector
        // can be transformed from world space into a local space
        //
        // Reference:
        // https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
        let inverse_rotation = glm::quat_inverse(&self.rotation);

        // d_local = R⁻¹ * d_world
        //
        // Quaternion representation
        //
        // d_local = R⁻¹ * d_world * q
        //
        // Transform the search direction into the collider's
        // local coordinate system
        //
        // The vector is interepreted as a pure quaterion
        // during quaterion multiplication
        //
        // Reference:
        // https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
        let local_direction = glm::quat_rotate_vec3(&inverse_rotation, world_diretion);

        // p_local = S_local(d_local)
        //
        // Find the furthest point of the collider along
        // the local-space search direction
        //
        // S_local(d) = argmax(p · d), for p ∈ C
        //
        // Reference:
        // https://en.wikipedia.org/wiki/Support_function
        let local_support = self.collider.support(local_direction);

        // p_world = t + R * p_local
        //
        // Quaternion representation:
        //
        // p_world = t + R * p_local * q⁻¹
        //
        // Rotate the local support point into world space
        // then translate it by the collider's world position
        //
        // Reference:
        // https://en.wikipedia.org/wiki/Rigid_transformation
        self.position + glm::quat_rotate_vec3(&self.rotation, &local_support)
    }
}
