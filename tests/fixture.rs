use earcut_rs::{deviation, Earcut};
use serde_json;
use std::fs;

fn test_fixture(path: &str, num_triangles: usize, expected_deviation: f64) {
    // load JSON
    let s = fs::read_to_string(path).unwrap();
    let expected = serde_json::from_str::<Vec<Vec<[f64; 2]>>>(&s).unwrap();

    // prepare input
    let data: Vec<_> = expected.clone().into_iter().flatten().flatten().collect();
    let hole_indices: Vec<_> = expected
        .iter()
        .map(|x| x.len() as u32)
        .scan(0, |sum, e| {
            *sum += e;
            Some(*sum)
        })
        .collect();
    let hole_indices = &hole_indices[0..hole_indices.len() - 1];
    let dim = 2;

    // earcut
    let mut triangles = vec![];
    let mut earcut = Earcut::new();
    earcut.earcut(&data, hole_indices, dim, &mut triangles);

    // check
    assert!(triangles.len() == num_triangles * 3);
    if triangles.len() > 0 {
        assert!(deviation(&data, hole_indices, dim, &triangles) <= expected_deviation);
    }
}

#[test]
fn fixture_building() {
    test_fixture("./tests/fixtures/building.json", 13, 0.0);
}

#[test]
fn fixture_dude() {
    test_fixture("./tests/fixtures/dude.json", 106, 2e-15);
}

#[test]
fn fixture_water() {
    test_fixture("./tests/fixtures/water.json", 2482, 0.0008);
}

#[test]
fn fixture_water2() {
    test_fixture("./tests/fixtures/water2.json", 1212, 0.0);
}

#[test]
fn fixture_water3() {
    test_fixture("./tests/fixtures/water3.json", 197, 0.0);
}

#[test]
fn fixture_water3b() {
    test_fixture("./tests/fixtures/water3b.json", 25, 0.0);
}

#[test]
fn fixture_water4() {
    test_fixture("./tests/fixtures/water4.json", 705, 0.0);
}

#[test]
fn fixture_water_huge() {
    test_fixture("./tests/fixtures/water-huge.json", 5177, 0.0011);
}

#[test]
fn fixture_water_huge2() {
    test_fixture("./tests/fixtures/water-huge2.json", 4462, 0.0028);
}

#[test]
fn fixture_degenerate() {
    test_fixture("./tests/fixtures/degenerate.json", 0, 0.0);
}

#[test]
fn fixture_bad_hole() {
    test_fixture("./tests/fixtures/bad-hole.json", 42, 0.019);
}

#[test]
fn fixture_empty_square() {
    test_fixture("./tests/fixtures/empty-square.json", 0, 0.0);
}

#[test]
fn fixture_issue16() {
    test_fixture("./tests/fixtures/issue16.json", 12, 4e-16);
}

#[test]
fn fixture_issue17() {
    test_fixture("./tests/fixtures/issue17.json", 11, 2e-16);
}

#[test]
fn fixture_steiner() {
    test_fixture("./tests/fixtures/steiner.json", 9, 0.0);
}

#[test]
fn fixture_issue29() {
    test_fixture("./tests/fixtures/issue29.json", 40, 2e-15);
}

#[test]
fn fixture_issue34() {
    test_fixture("./tests/fixtures/issue34.json", 139, 0.0);
}

#[test]
fn fixture_issue35() {
    test_fixture("./tests/fixtures/issue35.json", 844, 0.0);
}

#[test]
fn fixture_self_touching() {
    test_fixture("./tests/fixtures/self-touching.json", 124, 2e-13);
}

#[test]
fn fixture_outside_ring() {
    test_fixture("./tests/fixtures/outside-ring.json", 64, 0.0);
}

#[test]
fn fixture_simplified_us_border() {
    test_fixture("./tests/fixtures/simplified-us-border.json", 120, 0.0);
}

#[test]
fn fixture_touching_holes() {
    test_fixture("./tests/fixtures/touching-holes.json", 57, 0.0);
}

#[test]
fn fixture_hole_touching_outer() {
    test_fixture("./tests/fixtures/hole-touching-outer.json", 77, 0.0);
}

#[test]
fn fixture_hilbert() {
    test_fixture("./tests/fixtures/hilbert.json", 1024, 0.0);
}

#[test]
fn fixture_issue45() {
    test_fixture("./tests/fixtures/issue45.json", 10, 0.0);
}

#[test]
fn fixture_eberly_3() {
    test_fixture("./tests/fixtures/eberly-3.json", 73, 0.0);
}

#[test]
fn fixture_eberly_6() {
    test_fixture("./tests/fixtures/eberly-6.json", 1429, 2e-14);
}

#[test]
fn fixture_issue52() {
    test_fixture("./tests/fixtures/issue52.json", 109, 0.0);
}

#[test]
fn fixture_shared_points() {
    test_fixture("./tests/fixtures/shared-points.json", 4, 0.0);
}

#[test]
fn fixture_bad_diagonals() {
    test_fixture("./tests/fixtures/bad-diagonals.json", 7, 0.0);
}

#[test]
fn fixture_issue83() {
    test_fixture("./tests/fixtures/issue83.json", 0, 0.0);
}

#[test]
fn fixture_issue107() {
    test_fixture("./tests/fixtures/issue107.json", 0, 0.0);
}

#[test]
fn fixture_issue111() {
    test_fixture("./tests/fixtures/issue111.json", 19, 0.0);
}

#[test]
fn fixture_collinear_boxy() {
    test_fixture("./tests/fixtures/boxy.json", 57, 0.0);
}

#[test]
fn fixture_collinear_diagonal() {
    test_fixture("./tests/fixtures/collinear-diagonal.json", 14, 0.0);
}

#[test]
fn fixture_issue119() {
    test_fixture("./tests/fixtures/issue119.json", 18, 0.0);
}

#[test]
fn fixture_hourglass() {
    test_fixture("./tests/fixtures/hourglass.json", 2, 0.0);
}

#[test]
fn fixture_touching2() {
    test_fixture("./tests/fixtures/touching2.json", 8, 0.0);
}

#[test]
fn fixture_touching3() {
    test_fixture("./tests/fixtures/touching3.json", 15, 0.0);
}

#[test]
fn fixture_touching4() {
    test_fixture("./tests/fixtures/touching4.json", 20, 0.0);
}

#[test]
fn fixture_rain() {
    test_fixture("./tests/fixtures/rain.json", 2681, 0.0);
}

#[test]
fn fixture_issue131() {
    test_fixture("./tests/fixtures/issue131.json", 12, 0.0);
}

#[test]
fn fixture_infinite_loop_jhl() {
    test_fixture("./tests/fixtures/infinite-loop-jhl.json", 0, 0.0);
}

#[test]
fn fixture_filtered_bridge_jhl() {
    test_fixture("./tests/fixtures/filtered-bridge-jhl.json", 25, 0.0);
}

#[test]
fn fixture_issue149() {
    test_fixture("./tests/fixtures/issue149.json", 2, 0.0);
}

#[test]
fn fixture_issue142() {
    test_fixture("./tests/fixtures/issue142.json", 4, 0.13);
}
