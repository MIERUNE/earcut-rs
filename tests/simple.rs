use earcut_rs::{deviation, Earcut};

#[test]
fn test_empty() {
    let mut earcut = Earcut::new();
    let data: &[f64] = &[];
    let hole_indices: &[u32] = &[];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, 2, &mut triangles);
    assert_eq!(triangles.len(), 0);
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}

#[test]
fn test_invalid_point() {
    let mut earcut = Earcut::new();
    let data: &[f64] = &[100.0, 200.0];
    let hole_indices: &[u32] = &[];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, dim, &mut triangles);
    assert_eq!(triangles.len(), 0);
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}

#[test]
fn test_invalid_line() {
    let mut earcut = Earcut::new();
    let data: &[f64] = &[0.0, 0.0, 100.0, 200.0];
    let hole_indices: &[u32] = &[];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, dim, &mut triangles);
    assert_eq!(triangles.len(), 0);
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}

#[test]
fn test_invalid_empty_hole() {
    let mut earcut = Earcut::new();
    let data = &[0.0, 0.0, 100.0, 0.0, 100.0, 100.0];
    let hole_indices: &[u32] = &[3];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, dim, &mut triangles);
    assert_eq!(triangles.len(), 3);
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}

#[test]
fn test_steiner_point_hole() {
    let mut earcut = Earcut::new();
    let data = &[0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 50.0, 30.0];
    let hole_indices: &[u32] = &[3];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, dim, &mut triangles);
    assert_eq!(triangles.len(), 3 * 3);
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}

#[test]
fn test_steiner_line_hole() {
    let mut earcut = Earcut::new();
    let data = &[0., 0., 100., 0., 100., 100., 50., 30., 60., 30.];
    let hole_indices: &[u32] = &[3];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, dim, &mut triangles);
    assert_eq!(triangles.len(), 5 * 3);
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}

#[test]
fn test_square() {
    let mut earcut = Earcut::new();
    let data = &[0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0];
    let hole_indices: &[u32] = &[];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, dim, &mut triangles);
    assert_eq!(triangles, vec![2, 3, 0, 0, 1, 2]);
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}

#[test]
fn test_square_u16() {
    let mut earcut = Earcut::new();
    let data = &[0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0];
    let hole_indices: &[u16] = &[];
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, 2, &mut triangles);
    assert_eq!(triangles, vec![2, 3, 0, 0, 1, 2]);
    assert_eq!(deviation(data, hole_indices, 2, &triangles), 0.0);
}

#[test]
fn test_square_usize() {
    let mut earcut = Earcut::new();
    let data = &[0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0];
    let hole_indices: &[usize] = &[];
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, 2, &mut triangles);
    assert_eq!(triangles, vec![2, 3, 0, 0, 1, 2]);
    assert_eq!(deviation(data, hole_indices, 2, &triangles), 0.0);
}

#[test]
fn test_square_with_square_hole() {
    let mut earcut = Earcut::new();
    let data = &[
        0.0, 0.0, 100.0, 0.0, 100.0, 100.0, 0.0, 100.0, 10.0, 10.0, 90.0, 10.0, 90.0, 90.0, 10.0,
        90.0,
    ];
    let hole_indices: &[u32] = &[4];
    let dim = 2;
    let mut triangles = vec![];
    earcut.earcut(data, hole_indices, dim, &mut triangles);
    assert_eq!(
        triangles,
        vec![0, 4, 7, 5, 4, 0, 3, 0, 7, 5, 0, 1, 2, 3, 7, 6, 5, 1, 2, 7, 6, 6, 1, 2]
    );
    assert_eq!(deviation(data, hole_indices, dim, &triangles), 0.0);
}
