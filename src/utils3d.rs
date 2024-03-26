use alloc::vec::Vec;
use num_traits::float::Float;

#[inline]
fn cross<T: Float>([ax, ay, az]: [T; 3], [bx, by, bz]: [T; 3]) -> [T; 3] {
    [ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx]
}

fn normal<T: Float>(vertices: &[[T; 3]]) -> Option<[T; 3]> {
    let len = vertices.len();
    if len < 3 {
        // At least 3 vertices required
        return None;
    }
    let last_point = vertices[len - 1];

    let (sum, _) = vertices.iter().fold(
        ([T::zero(), T::zero(), T::zero()], last_point),
        |(acc, prev), data| {
            let (x, y, z) = (data[0], data[1], data[2]);
            let c = cross(
                [prev[0] - x, prev[1] - y, prev[2] - z],
                [prev[0] + x, prev[1] + y, prev[2] + z],
            );
            ([acc[0] + c[0], acc[1] + c[1], acc[2] + c[2]], [x, y, z])
        },
    );
    let d = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
    if d < T::from(1e-30).unwrap() {
        return None;
    }
    Some([sum[0] / d, sum[1] / d, sum[2] / d])
}

pub fn project3d_to_2d<T: Float>(
    vertices: &[[T; 3]],
    num_outer: usize,
    out_buf: &mut Vec<[T; 2]>,
) -> bool {
    let Some([nx, ny, nz]) = normal(&vertices[0..num_outer]) else {
        return false;
    };
    out_buf.clear();

    let dd = (nx * nx + ny * ny).sqrt();
    if dd < T::from(1e-15).unwrap() {
        if nz > T::zero() {
            // do nothing
            out_buf.extend(vertices.iter().map(|d| [d[0], d[1]]))
        } else {
            // flip
            out_buf.extend(vertices.iter().map(|d| [d[1], d[0]]))
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
        out_buf.extend(vertices.iter().map(|d| {
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
        let mut buf = Vec::new();
        let vertices = &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0]];
        assert!(project3d_to_2d(vertices, 3, &mut buf));
        assert!(buf == [[0., 0.], [2., 0.], [2., 2.]]);
    }

    #[test]
    fn test_flip() {
        let mut buf = Vec::new();
        let vertices = &[[0.0, 0.0, 0.0], [2.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        assert!(project3d_to_2d(vertices, 3, &mut buf));
        assert!(buf == [[0., 0.], [2., 2.], [0., 2.]]);
    }

    #[test]
    fn test_rotate() {
        let mut buf = Vec::new();
        let vertices = &[[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 2.0, 2.0]];
        assert!(project3d_to_2d(vertices, 3, &mut buf));
        assert!(buf == [[0., 0.], [2., 0.], [2., 2.]]);
    }

    #[test]
    fn test_invalid_input1() {
        let mut buf = Vec::new();
        let vertices: &[[f64; 3]; 0] = &[];
        assert!(!project3d_to_2d(vertices, 0, &mut buf));
    }

    #[test]
    fn test_invalid_input2() {
        // when normal is zero vector
        let vertices = &[
            [0., 0., 0.],
            [0., 1., 0.],
            [0., 0., 0.],
            [0., 0., 1.],
            [0., 0., 0.],
        ];
        assert!(normal(vertices).is_none());
    }
}
