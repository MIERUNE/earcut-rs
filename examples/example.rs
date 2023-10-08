use earcut_rs::{deviation, Earcut};
use serde_json;
use std::fs;

fn load_fixture(name: &str, num_triangles: usize, expected_deviation: f64) {
    // load JSON
    type Coords = Vec<Vec<[f64; 2]>>;
    let s = fs::read_to_string("./tests/fixtures/".to_string() + name + ".json").unwrap();
    let expected = serde_json::from_str::<Coords>(&s).unwrap();

    // prepare input
    let num_holes = expected.len();
    let data: Vec<_> = expected.clone().into_iter().flatten().flatten().collect();
    let hole_indices: Vec<_> = expected
        .into_iter()
        .map(|x| x.len() as u32)
        .scan(0, |sum, e| {
            *sum += e;
            Some(*sum)
        })
        .take(num_holes - 1)
        .collect();

    // earcut
    let mut triangles = vec![];
    let mut earcut = Earcut::new();
    for _ in 0..500 {
        earcut.earcut(&data, &hole_indices, 2, &mut triangles);
    }

    // check
    assert!(triangles.len() == num_triangles * 3);
    if triangles.len() > 0 {
        assert!(deviation(&data, &hole_indices, 2, &triangles) <= expected_deviation);
    }
}

fn main() {
    load_fixture("water", 2482, 0.0008);
}
