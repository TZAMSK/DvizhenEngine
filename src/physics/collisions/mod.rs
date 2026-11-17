use crate::physics::collisions::collision::WorldCollider;

mod collider;
pub mod collider_factory;
pub mod collision;

#[derive(Default)]
pub struct Contact {
    // world-space contact position
    pub point: glm::Vec3,

    // Unit-length collision normal
    //
    // ||n|| = 1
    //
    // Its orientation depends on the collision convention
    // such as pointing from collider A toward collider B
    pub normal: glm::Vec3,

    // penetration depth between collider
    //
    // depth >= 0
    //
    // Used during collision resolution to separate
    // interesting bodies
    pub depth: f32,
}

#[derive(Default, Clone, Copy)]
struct Vertex {
    // p = p_A - p_B
    //
    // Support point in the minkowski difference
    minkowski_different_point: glm::Vec3,

    // Original world-space support point on collider A
    single_point_shape_a: glm::Vec3,

    // Original world-space support point on collider B
    single_point_shape_b: glm::Vec3,
}

/// Returns a support vertex of the minkowski difference A - B
///
/// minkowski difference:
///
/// A - B = { a - b | a ∈ A, b ∈ B }
///
/// Support mapping:
///
/// S_(A-B)(d) = S_A(d) - S_B(d)
///
/// Reference:
/// https://en.wikipedia.org/wiki/Minkowski_addition
/// https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm
fn minkowski(
    collider_a: &WorldCollider<'_>,
    collider_b: &WorldCollider<'_>,
    search_direction: glm::Vec3,
) -> Vertex {
    // p_A = S_A(d)
    //
    // Find the furthest point of the collider A along
    // the search_direction
    let point_a = collider_a.support(&search_direction);

    // p_B = S_B(-d)
    //
    // Find the furthest point of the collider B along
    // the opposite search_direction
    //
    // This is equivalent to finding the point of B
    // with the smallest projection along d
    let point_b = collider_b.support(&-search_direction);

    // p = p_A - p_B
    //
    // Calculate the support point of the minkowski difference
    //
    // Mathematical derivation
    //
    // max((a-b) · d)
    //     = max(a · d - b · d)
    //     = max(a · d) + max(b · (-d))
    //
    // Therefore
    //
    // S_(A-B)(d) = S_A(d) - S_B(-d)
    //
    // GJK uses this point to constrct a simplex and
    // determine whether the minkowski difference contains
    // the origin
    //
    // A ∩ B ≠ ∅  ⇔  0 ∈ (A - B)
    //
    // The original support points are preserved so EPA
    // can later reconstruct world-space contact positions
    //
    // Reference:
    // https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm
    Vertex {
        minkowski_different_point: point_a - point_b,
        single_point_shape_a: point_a,
        single_point_shape_b: point_b,
    }
}

fn is_acute(a: &glm::Vec3, b: &glm::Vec3) -> bool {
    glm::dot(&a, &b) > 0.0
}

#[derive(Default, Clone, Copy)]
struct Simplex {
    vertices: [Vertex; 4],
    size: usize,
}

impl Simplex {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set(&mut self, new_vertices: &[Vertex]) {
        assert!(new_vertices.len() <= 4);
        self.size = new_vertices.len();
        self.vertices[..self.size].copy_from_slice(new_vertices);
    }

    pub fn push(&mut self, new_vertex: Vertex) {
        for i in (1..=self.size.min(3)).rev() {
            self.vertices[i] = self.vertices[i - 1];
        }

        self.vertices[0] = new_vertex;
        self.size = (self.size + 1).min(4);
    }
}

fn line_case(simplex: &mut Simplex, search_direction: &mut glm::Vec3) -> bool {
    let newest_vertex = simplex.vertices[0];
    let newest = newest_vertex.minkowski_different_point;
    let older = simplex.vertices[1].minkowski_different_point;

    let edge = older - newest;
    let to_origin = -newest;

    if is_acute(&edge, &to_origin) {
        *search_direction = glm::cross(&glm::cross(&edge, &to_origin), &edge);

        if glm::dot(&search_direction, &search_direction) < 1e-12 {
            *search_direction = if edge.x.abs() < 0.9 * glm::length(&edge) {
                glm::cross(&edge, &glm::vec3(1.0, 0.0, 0.0))
            } else {
                glm::cross(&edge, &glm::vec3(0.0, 1.0, 0.0))
            };
        }
    } else {
        simplex.set(&[newest_vertex]);
        *search_direction = to_origin;
    }

    false
}

fn triangle_case(simplex: &mut Simplex, search_direction: &mut glm::Vec3) -> bool {
    let newest_vertex = simplex.vertices[0];
    let vertex_b = simplex.vertices[1];
    let vertex_c = simplex.vertices[2];

    let newest = newest_vertex.minkowski_different_point;
    let edge_to_b = vertex_b.minkowski_different_point - newest;
    let edge_to_c = vertex_c.minkowski_different_point - newest;
    let to_origin = -newest;
    let face_normal = glm::cross(&edge_to_b, &edge_to_c);

    if is_acute(&glm::cross(&face_normal, &edge_to_c), &to_origin) {
        if is_acute(&edge_to_c, &to_origin) {
            simplex.set(&[newest_vertex, vertex_c]);
            *search_direction = glm::cross(&glm::cross(&edge_to_c, &to_origin), &edge_to_c);
        } else {
            simplex.set(&[newest_vertex, vertex_b]);
            return line_case(simplex, search_direction);
        }
    } else if is_acute(&glm::cross(&edge_to_b, &face_normal), &to_origin) {
        simplex.set(&[newest_vertex, vertex_b]);
        return line_case(simplex, search_direction);
    } else if is_acute(&face_normal, &to_origin) {
        *search_direction = face_normal;
    } else {
        simplex.set(&[newest_vertex, vertex_c, vertex_b]);
        *search_direction = -face_normal;
    }

    false
}

fn tetrahedron_case(simplex: &mut Simplex, search_direction: &mut glm::Vec3) -> bool {
    let newest_vertex = simplex.vertices[0];
    let vertex_b = simplex.vertices[1];
    let vertex_c = simplex.vertices[2];
    let vertex_d = simplex.vertices[3];

    let newest = newest_vertex.minkowski_different_point;
    let edge_to_b = vertex_b.minkowski_different_point - newest;
    let edge_to_c = vertex_c.minkowski_different_point - newest;
    let edge_to_d = vertex_d.minkowski_different_point - newest;
    let to_origin = -newest;

    let normal_abc = glm::cross(&edge_to_b, &edge_to_c);
    let normal_acd = glm::cross(&edge_to_c, &edge_to_d);
    let normal_adb = glm::cross(&edge_to_d, &edge_to_b);

    if is_acute(&normal_abc, &to_origin) {
        simplex.set(&[newest_vertex, vertex_b, vertex_c]);
        return triangle_case(simplex, search_direction);
    }
    if is_acute(&normal_acd, &to_origin) {
        simplex.set(&[newest_vertex, vertex_c, vertex_d]);
        return triangle_case(simplex, search_direction);
    }
    if is_acute(&normal_adb, &to_origin) {
        simplex.set(&[newest_vertex, vertex_d, vertex_b]);
        return triangle_case(simplex, search_direction);
    }

    true
}

fn next_simplex(simplex: &mut Simplex, search_direction: &mut glm::Vec3) -> bool {
    match simplex.size {
        2 => line_case(simplex, search_direction),
        3 => triangle_case(simplex, search_direction),
        4 => tetrahedron_case(simplex, search_direction),
        _ => false,
    }
}

fn gjk(
    collider_a: &WorldCollider<'_>,
    collider_b: &WorldCollider<'_>,
    simplex: &mut Simplex,
) -> bool {
    let mut search_direction = collider_b.position - collider_a.position;

    if glm::dot(&search_direction, &search_direction) < 1e-12 {
        search_direction = glm::vec3(1.0, 0.0, 0.0);
    }

    simplex.push(minkowski(collider_a, collider_b, search_direction));

    search_direction = -simplex.vertices[0].minkowski_different_point;

    if glm::dot(&search_direction, &search_direction) < 1e-12 {
        search_direction = glm::vec3(0.0, 1.0, 0.0);
    }

    for _ in 0..64 {
        let support_vertex = minkowski(collider_a, collider_b, search_direction);

        if glm::dot(&support_vertex.minkowski_different_point, &search_direction) <= 0.0 {
            return false;
        }

        simplex.push(support_vertex);

        if next_simplex(simplex, &mut search_direction) {
            return true;
        }
    }

    false
}

struct Face {
    vertex_indices: [usize; 3],
    normal: glm::Vec3,
    distance: f32,
}

fn make_face(vertices: &[Vertex], index_a: usize, index_b: usize, index_c: usize) -> Face {
    let mut face = Face {
        vertex_indices: [index_a, index_b, index_c],
        normal: glm::Vec3::zeros(),
        distance: f32::MAX,
    };

    let (point_a, point_b, point_c) = (
        vertices[index_a].minkowski_different_point,
        vertices[index_b].minkowski_different_point,
        vertices[index_c].minkowski_different_point,
    );

    let raw_normal = glm::cross(&(point_b - point_a), &(point_c - point_a));
    let length = glm::length(&raw_normal);

    if length < 1e-10 {
        return face;
    }

    face.normal = raw_normal / length;
    face.distance = glm::dot(&face.normal, &point_a);

    if face.distance < 0.0 {
        face.normal = -face.normal;
        face.distance = -face.distance;
        face.vertex_indices.swap(1, 2);
    }

    face
}

fn make_contact(vertices: &[Vertex], face: &Face) -> Contact {
    let (vertex0, vertex1, vertex2) = (
        vertices[face.vertex_indices[0]],
        vertices[face.vertex_indices[1]],
        vertices[face.vertex_indices[2]],
    );

    let closest_point = face.normal * face.distance;

    let (edge0, edge1) = (
        vertex1.minkowski_different_point - vertex0.minkowski_different_point,
        vertex2.minkowski_different_point - vertex0.minkowski_different_point,
    );

    let to_closest = closest_point - vertex0.minkowski_different_point;

    let (dot00, dot01, dot11) = (
        glm::dot(&edge0, &edge0),
        glm::dot(&edge0, &edge1),
        glm::dot(&edge1, &edge1),
    );

    let (dot_closest0, dot_closest1) =
        (glm::dot(&to_closest, &edge0), glm::dot(&to_closest, &edge1));

    let denominator = dot00 * dot11 - dot01 * dot01;

    let (mut weight0, mut weight1, mut weight2) = (1.0, 0.0, 0.0);

    if denominator.abs() > 1e-12 {
        weight1 = (dot11 * dot_closest0 - dot01 * dot_closest1) / denominator;
        weight2 = (dot00 * dot_closest1 - dot01 * dot_closest0) / denominator;
        weight0 = 1.0 - weight1 - weight2;
    }

    let (contact_point_a, contact_point_b) = (
        weight0 * vertex0.single_point_shape_a
            + weight1 * vertex1.single_point_shape_a
            + weight2 * vertex2.single_point_shape_a,
        weight0 * vertex0.single_point_shape_b
            + weight1 * vertex1.single_point_shape_b
            + weight2 * vertex2.single_point_shape_b,
    );

    let contact_point = (contact_point_a + contact_point_b) * 0.5;

    Contact {
        point: contact_point,
        normal: -face.normal,
        depth: face.distance,
    }
}

fn find_closest_face(faces: &[Face]) -> Option<&Face> {
    faces
        .iter()
        .min_by(|lhs, rhs| lhs.distance.total_cmp(&rhs.distance))
}

fn epa(
    collider_a: &WorldCollider<'_>,
    collider_b: &WorldCollider<'_>,
    simplex: &mut Simplex,
    out_contact: &mut Contact,
) -> bool {
    if simplex.size != 4 {
        return false;
    }

    let mut vertices = simplex.vertices.to_vec();

    let mut faces = vec![
        make_face(&vertices, 0, 1, 2),
        make_face(&vertices, 0, 3, 1),
        make_face(&vertices, 0, 2, 3),
        make_face(&vertices, 1, 3, 2),
    ];

    for _ in 0..64 {
        let Some(closest_face) = find_closest_face(&faces) else {
            return false;
        };

        let support_vertex = minkowski(collider_a, collider_b, closest_face.normal);

        if glm::dot(
            &closest_face.normal,
            &support_vertex.minkowski_different_point,
        ) - closest_face.distance
            < 1e-3
        {
            *out_contact = make_contact(&vertices, closest_face);
            return true;
        }

        let mut horizon_edges: Vec<(usize, usize)> = Vec::new();

        let mut add_edge = |from: usize, to: usize| {
            let reversed = horizon_edges
                .iter()
                .position(|&(a, b)| a == to && b == from);

            if let Some(index) = reversed {
                horizon_edges.remove(index);
            } else {
                horizon_edges.push((from, to));
            }
        };

        let mut face_index = 0;
        while face_index < faces.len() {
            let face = &faces[face_index];
            let face_point = vertices[face.vertex_indices[0]].minkowski_different_point;

            if glm::dot(
                &face.normal,
                &(support_vertex.minkowski_different_point - face_point),
            ) > 1e-6
            {
                add_edge(face.vertex_indices[0], face.vertex_indices[1]);
                add_edge(face.vertex_indices[1], face.vertex_indices[2]);
                add_edge(face.vertex_indices[2], face.vertex_indices[0]);
                faces.swap_remove(face_index);
            } else {
                face_index += 1;
            }
        }

        if horizon_edges.is_empty() {
            return false;
        }

        let new_vertex_index = vertices.len();
        vertices.push(support_vertex);
        for (from, to) in horizon_edges {
            faces.push(make_face(&vertices, from, to, new_vertex_index));
        }
    }

    let Some(closes_face) = find_closest_face(&faces) else {
        return false;
    };

    *out_contact = make_contact(&vertices, closes_face);
    true
}

pub fn collide(collider_a: &WorldCollider<'_>, collider_b: &WorldCollider<'_>) -> Option<Contact> {
    let mut simplex = Simplex::new();

    if !gjk(collider_a, collider_b, &mut simplex) {
        return None;
    }

    let mut contact = Contact::default();

    if epa(collider_a, collider_b, &mut simplex, &mut contact) {
        Some(contact)
    } else {
        None
    }
}
