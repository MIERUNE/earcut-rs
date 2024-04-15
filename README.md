# earcut-rs

[![Test](https://github.com/MIERUNE/earcut-rs/actions/workflows/Test.yml/badge.svg)](https://github.com/MIERUNE/earcut-rs/actions/workflows/Test.yml)
[![codecov](https://codecov.io/gh/MIERUNE/earcut-rs/graph/badge.svg?token=thKlQiVjLc)](https://codecov.io/gh/MIERUNE/earcut-rs)
[![Crates.io Version](https://img.shields.io/crates/v/earcut)](https://crates.io/crates/earcut)

A Rust port of the [mapbox/earcut](https://github.com/mapbox/earcut) polygon triangulation library, implemented from scratch with some reference to [donbright/earcutr](https://github.com/donbright/earcutr).

- Based on the latest earcut 2.2.4 release.
- Designed to avoid unnecessary memory allocations. You can reuse the internal buffer and the output index vector for multiple triangulations.
- (Experimental) An additional module, `utils3d`, can rotate 3D coplanar polygons into the 2D plane before triangulation.
- License: ISC

<p align="center">
<img src="./docs/image.png" width="300">
</p>


## Benchmarks

on Macbook Pro (M1 Pro)

| Polygon        | earcut.hpp | earcut-rs (0.4.0) | earcutr (0.4.3) |
|----------------|-------------:|--------------:|-------------:|
| bad_hole       |  3.574 µs/i  |   4.630 µs/i  |  4.415 µs/i  |          
| building       |    397 ns/i  |     197 ns/i  |    604 ns/i  |
| degenerate     |    142 ns/i  |      46 ns/i  |    206 ns/i  |
| dude           |  5.061 µs/i  |   6.516 µs/i  |  8.096 µs/i  |
| empty_square   |    195 ns/i  |      88 ns/i  |    331 ns/i  |
| water          |  459.6 µs/i  |   623.5 µs/i  |  801.3 µs/i  |
| water2         |  334.1 µs/i  |   409.5 µs/i  |  450.3 µs/i  |
| water3         |  13.12 µs/i  |   19.12 µs/i  |  23.46 µs/i  |
| water3b        |  1.340 µs/i  |   1.368 µs/i  |  2.165 µs/i  |
| water4         |  81.48 µs/i  |   122.0 µs/i  |  154.1 µs/i  |
| water_huge     |  6.906 ms/i  |   11.01 ms/i  |  10.90 ms/i  |
| water_huge2    |  15.38 ms/i  |   23.40 ms/i  |  22.35 ms/i  |

(Earcutr is not besed on the latest earcut.)
