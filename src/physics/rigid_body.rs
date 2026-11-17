use crate::physics::collisions::collider_factory::ColliderPart;

pub struct RigidBody {
    pub position: glm::Vec3,
    pub orentation: glm::Quat,
    pub velocity: glm::Vec3,
    pub angular_velocity: glm::Vec3,
    pub force: glm::Vec3,
    pub volume: f32,
    pub drag_coefficience: f32,
    pub cross_section: f32,
    pub inv_mass: f32,
    pub restitution: f32,
    pub inv_inertia_local: glm::Mat3,
    pub parts: Vec<ColliderPart>,
}

impl Default for RigidBody {
    fn default() -> Self {
        Self {
            position: glm::Vec3::zeros(),
            orentation: glm::quat_identity(),
            velocity: glm::Vec3::zeros(),
            angular_velocity: glm::Vec3::zeros(),
            force: glm::Vec3::zeros(),
            volume: 0.0,
            drag_coefficience: 0.0,
            cross_section: 0.0,
            inv_mass: 1.0,
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

    pub fn inv_intertia_world(&self) -> glm::Mat3 {
        let r = glm::quat_to_mat3(&self.orentation);
        r * self.inv_inertia_local * glm::transpose(&r)
    }

    pub fn compute_mass_properties(&mut self, mass: f32) {
        self.inv_mass = if mass > 0.0 { 1.0 / mass } else { 0.0 };

        self.inv_inertia_local = glm::Mat3::zeros();

        if self.is_static() || self.parts.is_empty() {
            return;
        }

        let part_mass = mass / self.parts.len() as f32;
        let identity = glm::Mat3::identity();
        let mut inertia = glm::Mat3::zeros();

        for part in &self.parts {
            let r = glm::quat_to_mat3(&part.local_rotation);
            let local_inertia = part.collider.intertia(part_mass);
            let rotated_intertia = r * local_inertia * glm::transpose(&r);
            let p = part.local_position;

            let distance_squared = glm::dot(&p, &p);

            let outer = p * glm::transpose(&p);

            let offset_inertia = part_mass * (distance_squared * identity - outer);

            inertia += rotated_intertia + offset_inertia;
        }

        self.inv_inertia_local = inertia.try_inverse().unwrap_or_else(glm::Mat3::zeros)
    }
}
