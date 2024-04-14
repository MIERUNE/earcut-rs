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
