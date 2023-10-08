use std::fs;

use criterion::{criterion_group, criterion_main, Criterion};

use earcut_rs::Earcut;

fn load_fixture(name: &str) -> (Vec<f64>, Vec<usize>) {
    // load JSON
    type Coords = Vec<Vec<[f64; 2]>>;
    let s = fs::read_to_string("./tests/fixtures/".to_string() + name + ".json").unwrap();
    let expected = serde_json::from_str::<Coords>(&s).unwrap();

    // prepare input
    let num_holes = expected.len();
    let data: Vec<_> = expected.clone().into_iter().flatten().flatten().collect();
    let hole_indices: Vec<_> = expected
        .iter()
        .map(|x| x.len() as usize)
        .scan(0, |sum, e| {
            *sum += e;
            Some(*sum)
        })
        .take(num_holes)
        .collect();

    return (data, hole_indices);
}

fn bench(c: &mut Criterion) {
    let mut earcut = Earcut::new();
    let mut triangles = Vec::new();

    c.bench_function("water", |b| {
        let (data, hole_indices) = load_fixture("water");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });

    c.bench_function("building", |b| {
        let (data, hole_indices) = load_fixture("building");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });

    c.bench_function("water2", |b| {
        let (data, hole_indices) = load_fixture("water2");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });

    c.bench_function("water3", |b| {
        let (data, hole_indices) = load_fixture("water3");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });

    c.bench_function("water3b", |b| {
        let (data, hole_indices) = load_fixture("water3b");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });

    c.bench_function("water-huge", |b| {
        let (data, hole_indices) = load_fixture("water-huge");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
            assert_eq!(triangles.len(), 5177 * 3)
        })
    });

    c.bench_function("water-huge2", |b| {
        let (data, hole_indices) = load_fixture("water-huge2");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });

    c.bench_function("rain", |b| {
        let (data, hole_indices) = load_fixture("rain");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });

    c.bench_function("hilbert", |b| {
        let (data, hole_indices) = load_fixture("hilbert");
        b.iter(|| {
            earcut.earcut(&data, &hole_indices, 2, &mut triangles);
        })
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
