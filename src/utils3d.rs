use alloc::vec::Vec;
use core::fmt::Debug;
use num_traits::float::Float;

#[inline]
fn cross<T: Float + Debug>((ax, ay, az): (T, T, T), (bx, by, bz): (T, T, T)) -> (T, T, T) {
    (ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx)
}

fn normal<T: Float + Debug>(vertices: &[T], dim: usize) -> Option<(T, T, T)> {
    let len = vertices.len();
    if len < dim * 3 {
        // At least 3 vertices required
        return None;
    }
    let last_point = (
        vertices[len - dim],
        vertices[len - dim + 1],
        vertices[len - dim + 2],
    );

    let (sum, _) = vertices.chunks_exact(dim).fold(
        ((T::zero(), T::zero(), T::zero()), last_point),
        |(acc, prev), data| {
            let (x, y, z) = (data[0], data[1], data[2]);
            let c = cross(
                (prev.0 - x, prev.1 - y, prev.2 - z),
                (prev.0 + x, prev.1 + y, prev.2 + z),
            );
            ((acc.0 + c.0, acc.1 + c.1, acc.2 + c.2), (x, y, z))
        },
    );
    let d = (sum.0 * sum.0 + sum.1 * sum.1 + sum.2 * sum.2).sqrt();
    if d < T::from(1e-30).unwrap() {
        return None;
    }
    Some((sum.0 / d, sum.1 / d, sum.2 / d))
}

pub fn project3d_to_2d<T: Float + Debug>(
    vertices: &[T],
    num_outer: usize,
    dim: usize,
    buf: &mut Vec<T>,
) -> bool {
    let Some((nx, ny, nz)) = normal(&vertices[0..num_outer * dim], dim) else {
        return false;
    };
    buf.clear();

    let dd = (nx * nx + ny * ny).sqrt();
    if dd < T::from(1e-15).unwrap() {
        if nz > T::zero() {
            // do nothing
            buf.extend(vertices.chunks_exact(dim).flat_map(|d| [d[0], d[1]]))
        } else {
            // flip
            buf.extend(vertices.chunks_exact(dim).flat_map(|d| [d[1], d[0]]))
        }
    } else {
        // rotation
        let ax = -ny / dd;
        let ay = nx / dd;
        let theta = nz.acos();
        let sint = theta.sin();
        let cost = theta.cos();
        let s = ax * ay * (T::one() - cost);
        let t = ay * sint;
        let u = ax * sint;
        let m11 = ax * ax * (T::one() - cost) + cost;
        let m12 = s;
        let m13 = -t;
        let m21 = s;
        let m22 = ay * ay * (T::one() - cost) + cost;
        let m23 = u;
        buf.extend(vertices.chunks_exact(dim).flat_map(|d| {
            let (x, y, z) = (d[0], d[1], d[2]);
            [(x * m11 + y * m12 + z * m13), (x * m21 + y * m22 + z * m23)]
        }))
    }
    true
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_do_nothing() {
        // dim=3
        let mut buf = Vec::new();
        let vertices = &[0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0];
        assert!(project3d_to_2d(vertices, 3, 3, &mut buf));
        assert!(buf == [0., 0., 2., 0., 2., 2.]);

        // dim=5
        let mut buf = Vec::new();
        let vertices = &[
            0.0, 0.0, 0.0, 99.0, 99.0, 2.0, 0.0, 0.0, 99.0, 99.0, 2.0, 2.0, 0.0, 99.0, 99.0,
        ];
        assert!(project3d_to_2d(vertices, 3, 5, &mut buf));
        assert!(buf == [0., 0., 2., 0., 2., 2.]);
    }

    #[test]
    fn test_flip() {
        // dim=3
        let mut buf = Vec::new();
        let vertices = &[0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 2.0, 0.0, 0.0];
        assert!(project3d_to_2d(vertices, 3, 3, &mut buf));
        assert!(buf == [0., 0., 2., 2., 0., 2.]);

        // dim=4
        let mut buf = Vec::new();
        let vertices = &[
            0.0, 0.0, 0.0, 99.0, 2.0, 2.0, 0.0, 99.0, 2.0, 0.0, 0.0, 99.0,
        ];
        assert!(project3d_to_2d(vertices, 3, 4, &mut buf));
        assert!(buf == [0., 0., 2., 2., 0., 2.]);
    }

    #[test]
    fn test_rotate() {
        // dim=3
        let mut buf = Vec::new();
        let vertices = &[0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 2.0];
        assert!(project3d_to_2d(vertices, 3, 3, &mut buf));
        assert!(buf == [0., 0., 2., 0., 2., 2.]);

        // dim=4
        let mut buf = Vec::new();
        let vertices = &[
            0.0, 0.0, 0.0, 99.0, 0.0, 0.0, 2.0, 99.0, 0.0, 2.0, 2.0, 99.0,
        ];
        assert!(project3d_to_2d(vertices, 3, 4, &mut buf));
        assert!(buf == [0., 0., 2., 0., 2., 2.]);
    }

    #[test]
    fn test_invalid_input1() {
        // dim=3
        let mut buf = Vec::new();
        let vertices = &[0., 0., 1., 1.];
        assert!(!project3d_to_2d(vertices, 1, 3, &mut buf));

        // dim=4
        let mut buf = Vec::new();
        let vertices = &[0., 0., 1., 1., 1.];
        assert!(!project3d_to_2d(vertices, 1, 4, &mut buf));
    }

    #[test]
    fn test_invalid_input2() {
        // when normal is zero vector
        let vertices = &[0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.];
        assert!(normal(vertices, 3).is_none());
    }
}
