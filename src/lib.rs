//! A Rust port of the [Earcut](https://github.com/mapbox/earcut) polygon triangulation library.

#![no_std]

extern crate alloc;

pub mod utils3d;

use alloc::vec::Vec;
use core::ptr;
use num_traits::float::Float;

/// Index of a vertex
pub trait Index: Copy {
    fn into_usize(self) -> usize;
    fn from_usize(v: usize) -> Self;
}
impl Index for u32 {
    fn into_usize(self) -> usize {
        self as usize
    }
    fn from_usize(v: usize) -> Self {
        v as Self
    }
}
impl Index for u16 {
    fn into_usize(self) -> usize {
        self as usize
    }
    fn from_usize(v: usize) -> Self {
        v as Self
    }
}
impl Index for usize {
    fn into_usize(self) -> usize {
        self
    }
    fn from_usize(v: usize) -> Self {
        v as Self
    }
}

macro_rules! node {
    ($self:ident.$nodes:ident, $index:expr) => {
        unsafe { $self.$nodes.get_unchecked($index as usize) }
    };
    ($nodes:ident, $index:expr) => {
        unsafe { $nodes.get_unchecked($index as usize) }
    };
}

macro_rules! node_mut {
    ($self:ident.$nodes:ident, $index:expr) => {
        unsafe { $self.$nodes.get_unchecked_mut($index as usize) }
    };
    ($nodes:ident, $index:expr) => {
        unsafe { $nodes.get_unchecked_mut($index as usize) }
    };
}

struct Node<T: Float> {
    /// vertex index in coordinates array
    i: u32,
    /// vertex coordinates x
    x: T,
    /// vertex coordinates y
    y: T,
    /// previous vertex nodes in a polygon ring
    prev_i: u32,
    /// next vertex nodes in a polygon ring
    next_i: u32,
    /// z-order curve value
    z: i32,
    /// previous nodes in z-order
    prev_z_i: Option<u32>,
    /// next nodes in z-order
    next_z_i: Option<u32>,
    /// indicates whether this is a steiner point
    steiner: bool,
}

impl<T: Float> Node<T> {
    fn new(i: u32, x: T, y: T) -> Self {
        Self {
            i,
            x,
            y,
            prev_i: 0,
            next_i: 0,
            z: 0,
            prev_z_i: None,
            next_z_i: None,
            steiner: false,
        }
    }
}

/// Instance of the earcut algorithm.
pub struct Earcut<T: Float> {
    data: Vec<T>,
    nodes: Vec<Node<T>>,
    queue: Vec<u32>,
}

impl<T: Float> Default for Earcut<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Float> Earcut<T> {
    /// Creates a new instance of the earcut algorithm.
    ///
    /// You can reuse a single instance for multiple triangulations to reduce memory allocations.
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            nodes: Vec::new(),
            queue: Vec::new(),
        }
    }

    fn reset(&mut self, capacity: usize) {
        self.nodes.clear();
        self.nodes.reserve(capacity);
    }

    /// Performs the earcut triangulation on a polygon.
    ///
    /// The API is similar to the original JavaScript implementation, except you can provide a vector for the output indices.
    pub fn earcut<N: Index>(
        &mut self,
        data: impl IntoIterator<Item = T>,
        hole_indices: &[N],
        triangles_out: &mut Vec<N>,
    ) {
        self.data.clear();
        self.data.extend(data);

        triangles_out.clear();
        self.reset(self.data.len() / 2 * 3 / 2);

        let has_holes = !hole_indices.is_empty();
        let outer_len: usize = if has_holes {
            hole_indices[0].into_usize() * 2
        } else {
            self.data.len()
        };
        let Some(mut outer_node_i) = self.linked_list(0, outer_len as u32, true) else {
            return;
        };

        let outer_node = node!(self.nodes, outer_node_i);
        if outer_node.next_i == outer_node.prev_i {
            return;
        }

        if has_holes {
            outer_node_i = self.eliminate_holes(hole_indices, outer_node_i);
        }

        let mut min_x = T::zero();
        let mut min_y = T::zero();
        let mut inv_size = T::zero();

        // if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
        if self.data.len() > 80 * 2 {
            let max_x = self.data[2..outer_len]
                .iter()
                .step_by(2)
                .fold(self.data[0], |a, b| T::max(a, *b));
            min_x = self.data[2..outer_len]
                .iter()
                .step_by(2)
                .fold(self.data[0], |a, b| T::min(a, *b));
            let max_y = self.data[2 + 1..outer_len]
                .iter()
                .step_by(2)
                .fold(self.data[1], |a, b| T::max(a, *b));
            min_y = self.data[2 + 1..outer_len]
                .iter()
                .step_by(2)
                .fold(self.data[1], |a, b| T::min(a, *b));

            // minX, minY and invSize are later used to transform coords into integers for z-order calculation
            inv_size = (max_x - min_x).max(max_y - min_y);
            if inv_size != T::zero() {
                inv_size = T::from(32767.0).unwrap() / inv_size;
            }
        }

        self.earcut_linked(outer_node_i, triangles_out, min_x, min_y, inv_size, 0);
    }

    /// create a circular doubly linked list from polygon points in the specified winding order
    fn linked_list(&mut self, start: u32, end: u32, clockwise: bool) -> Option<u32> {
        let mut last_i: Option<u32> = None;
        if start >= end {
            return None;
        }

        if clockwise == (signed_area(&self.data, start, end) > T::zero()) {
            for (i, v) in self.data[start as usize..end as usize]
                .chunks_exact(2)
                .enumerate()
            {
                let idx = start + i as u32 * 2;
                last_i = Some(insert_node(&mut self.nodes, idx, v[0], v[1], last_i));
            }
        } else {
            for (i, v) in self.data[start as usize..end as usize]
                .chunks_exact(2)
                .enumerate()
                .rev()
            {
                let idx = start + i as u32 * 2;
                last_i = Some(insert_node(&mut self.nodes, idx, v[0], v[1], last_i));
            }
        };

        if let Some(li) = last_i {
            let last = &node!(self.nodes, li);
            if equals(last, node!(self.nodes, last.next_i)) {
                let (_, next_i) = remove_node(&mut self.nodes, li);
                last_i = Some(next_i);
            }
        }

        last_i
    }

    /// eliminate colinear or duplicate points
    fn filter_points(&mut self, start_i: u32, end_i: Option<u32>) -> u32 {
        let mut end_i = end_i.unwrap_or(start_i);

        let mut p_i = start_i;
        loop {
            let p = node!(self.nodes, p_i);
            let p_next = node!(self.nodes, p.next_i);
            if !p.steiner
                && (equals(p, p_next) || area(node!(self.nodes, p.prev_i), p, p_next).is_zero())
            {
                let (prev_i, next_i) = remove_node(&mut self.nodes, p_i);
                (p_i, end_i) = (prev_i, prev_i);
                if p_i == next_i {
                    return end_i;
                }
            } else {
                p_i = p.next_i;
                if p_i == end_i {
                    return end_i;
                }
            };
        }
    }

    /// main ear slicing loop which triangulates a polygon (given as a linked list)
    #[allow(clippy::too_many_arguments)]
    fn earcut_linked<N: Index>(
        &mut self,
        ear_i: u32,
        triangles: &mut Vec<N>,
        min_x: T,
        min_y: T,
        inv_size: T,
        pass: usize,
    ) {
        let mut ear_i = ear_i;

        // interlink polygon nodes in z-order
        if pass == 0 && inv_size != T::zero() {
            self.index_curve(ear_i, min_x, min_y, inv_size);
        }

        let mut stop_i = ear_i;

        // iterate through ears, slicing them one by one
        while node!(self.nodes, ear_i).prev_i != node!(self.nodes, ear_i).next_i {
            let ear = node!(self.nodes, ear_i);
            let prev_i = ear.prev_i;
            let next_i = ear.next_i;

            let is_ear = if inv_size != T::zero() {
                self.is_ear_hashed(ear_i, min_x, min_y, inv_size)
            } else {
                self.is_ear(ear_i)
            };
            if is_ear {
                // cut off the triangle
                triangles.extend([
                    N::from_usize(node!(self.nodes, prev_i).i as usize / 2),
                    N::from_usize(ear.i as usize / 2),
                    N::from_usize(node!(self.nodes, next_i).i as usize / 2),
                ]);

                remove_node(&mut self.nodes, ear_i);

                // skipping the next vertex leads to less sliver triangles
                let next_next_i = node!(self.nodes, next_i).next_i;
                (ear_i, stop_i) = (next_next_i, next_next_i);

                continue;
            }

            ear_i = next_i;

            // if we looped through the whole remaining polygon and can't find any more ears
            if ear_i == stop_i {
                if pass == 0 {
                    // try filtering points and slicing again
                    ear_i = self.filter_points(ear_i, None);
                    self.earcut_linked(ear_i, triangles, min_x, min_y, inv_size, 1);
                } else if pass == 1 {
                    // if this didn't work, try curing all small self-intersections locally
                    let filtered = self.filter_points(ear_i, None);
                    ear_i = self.cure_local_intersections(filtered, triangles);
                    self.earcut_linked(ear_i, triangles, min_x, min_y, inv_size, 2);
                } else if pass == 2 {
                    // as a last resort, try splitting the remaining polygon into two
                    self.split_earcut(ear_i, triangles, min_x, min_y, inv_size);
                }
                return;
            }
        }
    }

    /// check whether a polygon node forms a valid ear with adjacent nodes
    fn is_ear(&self, ear_i: u32) -> bool {
        let b = node!(self.nodes, ear_i);
        let a = node!(self.nodes, b.prev_i);
        let c = node!(self.nodes, b.next_i);

        if area(a, b, c) >= T::zero() {
            // reflex, can't be an ear
            return false;
        }

        // now make sure we don't have other points inside the potential ear

        // triangle bbox
        let x0 = a.x.min(b.x.min(c.x));
        let y0 = a.y.min(b.y.min(c.y));
        let x1 = a.x.max(b.x.max(c.x));
        let y1 = a.y.max(b.y.max(c.y));

        let mut p = node!(self.nodes, c.next_i);
        let mut p_prev = node!(self.nodes, p.prev_i);
        while !ptr::eq(p, a) {
            let p_next = node!(self.nodes, p.next_i);
            if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1)
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y)
                && area(p_prev, p, p_next) >= T::zero()
            {
                return false;
            }
            (p_prev, p) = (p, p_next);
        }
        true
    }

    fn is_ear_hashed(&self, ear_i: u32, min_x: T, min_y: T, inv_size: T) -> bool {
        let b = node!(self.nodes, ear_i);
        let a = node!(self.nodes, b.prev_i);
        let c = node!(self.nodes, b.next_i);

        if area(a, b, c) >= T::zero() {
            // reflex, can't be an ear
            return false;
        }

        // triangle bbox
        let x0 = a.x.min(b.x.min(c.x));
        let y0 = a.y.min(b.y.min(c.y));
        let x1 = a.x.max(b.x.max(c.x));
        let y1 = a.y.max(b.y.max(c.y));

        // z-order range for the current triangle bbox;
        let min_z = z_order(x0, y0, min_x, min_y, inv_size);
        let max_z = z_order(x1, y1, min_x, min_y, inv_size);

        let ear = node!(self.nodes, ear_i);
        let mut o_p = ear.prev_z_i.map(|i| node!(self.nodes, i as usize));
        let mut o_n = ear.next_z_i.map(|i| node!(self.nodes, i));

        let ear_prev = node!(self.nodes, ear.prev_i);
        let ear_next = node!(self.nodes, ear.next_i);

        // look for points inside the triangle in both directions
        loop {
            let Some(p) = o_p else { break };
            if p.z < min_z {
                break;
            };
            let Some(n) = o_n else { break };
            if n.z > max_z {
                break;
            };

            if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1)
                && (!ptr::eq(p, a) && !ptr::eq(p, c))
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y)
                && area(node!(self.nodes, p.prev_i), p, node!(self.nodes, p.next_i)) >= T::zero()
            {
                return false;
            }
            o_p = p.prev_z_i.map(|i| node!(self.nodes, i));

            if (n.x >= x0 && n.x <= x1 && n.y >= y0 && n.y <= y1)
                && (!ptr::eq(n, a) && !ptr::eq(n, c))
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y)
                && area(node!(self.nodes, n.prev_i), n, node!(self.nodes, n.next_i)) >= T::zero()
            {
                return false;
            }
            o_n = p.next_z_i.map(|i| node!(self.nodes, i));
        }

        // look for remaining points in decreasing z-order
        while let Some(p) = o_p {
            if p.z < min_z {
                break;
            };
            if (!ptr::eq(p, ear_prev) && !ptr::eq(p, ear_next))
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y)
                && area(node!(self.nodes, p.prev_i), p, node!(self.nodes, p.next_i)) >= T::zero()
            {
                return false;
            }
            o_p = p.prev_z_i.map(|i| node!(self.nodes, i));
        }

        // look for remaining points in increasing z-order
        while let Some(n) = o_n {
            if n.z > max_z {
                break;
            };
            if (!ptr::eq(n, ear_prev) && !ptr::eq(n, ear_next))
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y)
                && area(node!(self.nodes, n.prev_i), n, node!(self.nodes, n.next_i)) >= T::zero()
            {
                return false;
            }
            o_n = n.next_z_i.map(|i| node!(self.nodes, i));
        }

        true
    }

    /// go through all polygon nodes and cure small local self-intersections
    fn cure_local_intersections<N: Index>(
        &mut self,
        mut start_i: u32,
        triangles: &mut Vec<N>,
    ) -> u32 {
        let mut p_i = start_i;
        loop {
            let p_prev_i = node!(self.nodes, p_i).prev_i;
            let p_next_i = node!(self.nodes, p_i).next_i;
            let a_i = p_prev_i;
            let p_next = node!(self.nodes, p_next_i);
            let b_i = p_next.next_i;
            let a = node!(self.nodes, a_i);
            let b = node!(self.nodes, b_i);
            let p = node!(self.nodes, p_i);

            if !equals(a, b)
                && self.intersects(a, p, p_next, b)
                && self.locally_inside(a, b)
                && self.locally_inside(b, a)
            {
                triangles.extend([
                    N::from_usize(a.i as usize / 2),
                    N::from_usize(p.i as usize / 2),
                    N::from_usize(b.i as usize / 2),
                ]);

                remove_node(&mut self.nodes, p_i);
                remove_node(&mut self.nodes, p_next_i);

                (p_i, start_i) = (b_i, b_i);
            }

            p_i = node!(self.nodes, p_i).next_i;
            if p_i == start_i {
                return self.filter_points(p_i, None);
            }
        }
    }

    /// try splitting polygon into two and triangulate them independently
    fn split_earcut<N: Index>(
        &mut self,
        start_i: u32,
        triangles: &mut Vec<N>,
        min_x: T,
        min_y: T,
        inv_size: T,
    ) {
        // look for a valid diagonal that divides the polygon into two
        let mut a_i = start_i;
        loop {
            let ai = node!(self.nodes, a_i).i;
            let a_prev_i = node!(self.nodes, a_i).prev_i;
            let a_next_i = node!(self.nodes, a_i).next_i;
            let mut b_i = (node!(self.nodes, a_next_i)).next_i;

            while b_i != a_prev_i {
                let b = node!(self.nodes, b_i);
                if ai != b.i && self.is_valid_diagonal(node!(self.nodes, a_i), b) {
                    // split the polygon in two by the diagonal
                    let mut c_i = self.split_polygon(a_i, b_i);

                    // filter colinear points around the cuts
                    a_i = self.filter_points(a_i, Some(node!(self.nodes, a_i).next_i));
                    c_i = self.filter_points(c_i, Some(node!(self.nodes, c_i).next_i));

                    // run earcut on each half
                    self.earcut_linked(a_i, triangles, min_x, min_y, inv_size, 0);
                    self.earcut_linked(c_i, triangles, min_x, min_y, inv_size, 0);
                    return;
                }
                b_i = b.next_i;
            }

            a_i = node!(self.nodes, a_i).next_i;
            if a_i == start_i {
                return;
            }
        }
    }

    /// link every hole into the outer loop, producing a single-ring polygon without holes
    fn eliminate_holes<N: Index>(&mut self, hole_indices: &[N], mut outer_node_i: u32) -> u32 {
        self.queue.clear();
        let len = hole_indices.len();
        for (i, hi) in hole_indices.iter().enumerate() {
            let start = (*hi).into_usize() * 2;
            let end = if i < len - 1 {
                hole_indices[i + 1].into_usize() * 2
            } else {
                self.data.len()
            };
            if let Some(list_i) = self.linked_list(start as u32, end as u32, false) {
                let list = &mut node_mut!(self.nodes, list_i);
                if list_i == list.next_i {
                    list.steiner = true;
                }
                self.queue.push(self.get_leftmost(list_i))
            }
        }

        self.queue.sort_unstable_by(|a, b| {
            node!(self.nodes, *a)
                .x
                .partial_cmp(&node!(self.nodes, *b).x)
                .unwrap_or(core::cmp::Ordering::Equal)
        });

        // process holes from left to right
        for i in 0..self.queue.len() {
            outer_node_i = self.eliminate_hole(self.queue[i], outer_node_i);
        }

        outer_node_i
    }

    /// find a bridge between vertices that connects hole with an outer ring and and link it
    fn eliminate_hole(&mut self, hole_i: u32, outer_node_i: u32) -> u32 {
        let Some(bridge_i) = self.find_hole_bridge(node!(self.nodes, hole_i), outer_node_i) else {
            return outer_node_i;
        };
        let bridge_reverse_i = self.split_polygon(bridge_i, hole_i);

        // filter collinear points around the cuts
        self.filter_points(
            bridge_reverse_i,
            Some(node!(self.nodes, bridge_reverse_i).next_i),
        );
        self.filter_points(bridge_i, Some(node!(self.nodes, bridge_i).next_i))
    }

    /// dimavid Eberly's algorithm for finding a bridge between hole and outer polygon
    fn find_hole_bridge(&self, hole: &Node<T>, outer_node_i: u32) -> Option<u32> {
        let mut p_i = outer_node_i;
        let mut qx = T::neg_infinity();
        let mut m_i: Option<u32> = None;

        // find a segment intersected by a ray from the hole's leftmost point to the left;
        // segment's endpoint with lesser x will be potential connection point
        loop {
            let p = node!(self.nodes, p_i);
            let p_next_i = p.next_i;
            let p_next = node!(self.nodes, p_next_i);
            if hole.y <= p.y && hole.y >= p_next.y && p_next.y != p.y {
                let x = p.x + (hole.y - p.y) * (p_next.x - p.x) / (p_next.y - p.y);
                if x <= hole.x && x > qx {
                    qx = x;
                    m_i = Some(if p.x < p_next.x { p_i } else { p_next_i });
                    if x == hole.x {
                        // hole touches outer segment; pick leftmost endpoint
                        return m_i;
                    }
                }
            }
            p_i = p_next_i;
            if p_i == outer_node_i {
                break;
            }
        }

        let mut m_i = m_i?;

        // look for points inside the triangle of hole point, segment intersection and endpoint;
        // if there are no points found, we have a valid connection;
        // otherwise choose the point of the minimum angle with the ray as connection point

        let stop_i = m_i;
        let Node { x: mx, y: my, .. } = *node!(self.nodes, m_i); // must copy
        let mut tan_min = T::infinity();

        p_i = m_i;
        let mut p = node!(self.nodes, p_i);
        let mut m = p;

        loop {
            if (hole.x >= p.x && p.x >= mx && hole.x != p.x)
                && point_in_triangle(
                    if hole.y < my { hole.x } else { qx },
                    hole.y,
                    mx,
                    my,
                    if hole.y < my { qx } else { hole.x },
                    hole.y,
                    p.x,
                    p.y,
                )
            {
                let tan = (hole.y - p.y).abs() / (hole.x - p.x);
                if self.locally_inside(p, hole)
                    && (tan < tan_min
                        || (tan == tan_min
                            && (p.x > m.x || p.x == m.x && self.sector_contains_sector(m, p))))
                {
                    (m_i, m) = (p_i, p);
                    tan_min = tan;
                }
            }

            p_i = p.next_i;
            if p_i == stop_i {
                return Some(m_i);
            }
            p = node!(self.nodes, p_i);
        }
    }

    /// whether sector in vertex m contains sector in vertex p in the same coordinates
    fn sector_contains_sector(&self, m: &Node<T>, p: &Node<T>) -> bool {
        area(node!(self.nodes, m.prev_i), m, node!(self.nodes, p.prev_i)) < T::zero()
            && area(node!(self.nodes, p.next_i), m, node!(self.nodes, m.next_i)) < T::zero()
    }

    /// interlink polygon nodes in z-order
    fn index_curve(&mut self, start_i: u32, min_x: T, min_y: T, inv_size: T) {
        let mut p_i = start_i;

        loop {
            let p = node_mut!(self.nodes, p_i);
            if p.z == 0 {
                p.z = z_order(p.x, p.y, min_x, min_y, inv_size);
            }
            p.prev_z_i = Some(p.prev_i);
            p.next_z_i = Some(p.next_i);
            p_i = p.next_i;
            if p_i == start_i {
                break;
            }
        }

        let p_prev_z_i = node!(self.nodes, p_i).prev_z_i.unwrap();
        node_mut!(self.nodes, p_prev_z_i).next_z_i = None;
        node_mut!(self.nodes, p_i).prev_z_i = None;
        self.sort_linked(p_i);
    }

    /// Simon Tatham's linked list merge sort algorithm
    /// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
    fn sort_linked(&mut self, list_i: u32) {
        let mut in_size: usize = 1;
        let mut list_i = Some(list_i);

        loop {
            let mut p_i = list_i;
            list_i = None;
            let mut tail_i: Option<u32> = None;
            let mut num_merges = 0;

            while let Some(pp) = p_i {
                num_merges += 1;
                let mut q_i = node!(self.nodes, pp).next_z_i;
                let mut p_size: u32 = 1;
                for _ in 1..in_size {
                    if let Some(i) = q_i {
                        p_size += 1;
                        q_i = node!(self.nodes, i).next_z_i;
                    } else {
                        q_i = None;
                        break;
                    }
                }
                let mut q_size = in_size;

                while p_size > 0 || (q_size > 0 && q_i.is_some()) {
                    let (e_i, e) = if p_size == 0 {
                        q_size -= 1;
                        let e_i = q_i.unwrap();
                        let e = node_mut!(self.nodes, e_i);
                        q_i = e.next_z_i;
                        (e_i, e)
                    } else if q_size == 0 {
                        p_size -= 1;
                        let e_i = p_i.unwrap();
                        let e = node_mut!(self.nodes, e_i);
                        p_i = e.next_z_i;
                        (e_i, e)
                    } else {
                        let p_i_s = p_i.unwrap();
                        if let Some(q_i_s) = q_i {
                            if node!(self.nodes, p_i_s).z <= node!(self.nodes, q_i_s).z {
                                p_size -= 1;
                                let e_i = p_i_s;
                                let e = node_mut!(self.nodes, p_i_s);
                                p_i = e.next_z_i;
                                (e_i, e)
                            } else {
                                q_size -= 1;
                                let e_i = q_i_s;
                                let e = node_mut!(self.nodes, q_i_s);
                                q_i = e.next_z_i;
                                (e_i, e)
                            }
                        } else {
                            p_size -= 1;
                            let e_i = p_i_s;
                            let e = node_mut!(self.nodes, e_i);
                            p_i = e.next_z_i;
                            (e_i, e)
                        }
                    };

                    e.prev_z_i = tail_i;

                    if let Some(tail_i) = tail_i {
                        node_mut!(self.nodes, tail_i).next_z_i = Some(e_i);
                    } else {
                        list_i = Some(e_i);
                    }
                    tail_i = Some(e_i);
                }

                p_i = q_i;
            }

            node_mut!(self.nodes, tail_i.unwrap()).next_z_i = None;
            if num_merges <= 1 {
                break;
            }
            in_size *= 2;
        }
    }

    /// find the leftmost node of a polygon ring
    fn get_leftmost(&self, start_i: u32) -> u32 {
        let mut p_i = start_i;
        let mut leftmost_i = start_i;
        let mut p = node!(self.nodes, p_i);
        let mut leftmost = p;

        loop {
            if p.x < leftmost.x || (p.x == leftmost.x && p.y < leftmost.y) {
                (leftmost_i, leftmost) = (p_i, p);
            }
            p_i = p.next_i;
            if p_i == start_i {
                return leftmost_i;
            }
            p = node!(self.nodes, p_i);
        }
    }

    /// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
    fn is_valid_diagonal(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let a_next = node!(self.nodes, a.next_i);
        let a_prev = node!(self.nodes, a.prev_i);
        let b_next = node!(self.nodes, b.next_i);
        let b_prev = node!(self.nodes, b.prev_i);
        // dones't intersect other edges
        (a_next.i != b.i && a_prev.i != b.i && !self.intersects_polygon(a, b))
            // locally visible
            && ((self.locally_inside(a, b) && self.locally_inside(b, a) && self.middle_inside(a, b))
                // does not create opposite-facing sectors
                && (area(a_prev, a, b_prev) != T::zero() || area(a, b_prev, b) != T::zero())
                // special zero-length case
                || equals(a, b)
                    && area(a_prev, a, a_next) > T::zero()
                    && area(b_prev, b, b_next) > T::zero())
    }

    /// check if two segments intersect
    fn intersects(&self, p1: &Node<T>, q1: &Node<T>, p2: &Node<T>, q2: &Node<T>) -> bool {
        let o1 = sign(area(p1, q1, p2));
        let o2 = sign(area(p1, q1, q2));
        let o3 = sign(area(p2, q2, p1));
        let o4 = sign(area(p2, q2, q1));
        (o1 != o2 && o3 != o4) // general case
            || (o1 == 0 && on_segment(p1, p2, q1)) // p1, q1 and p2 are collinear and p2 lies on p1q1
            || (o2 == 0 && on_segment(p1, q2, q1)) // p1, q1 and q2 are collinear and q2 lies on p1q1
            || (o3 == 0 && on_segment(p2, p1, q2)) // p2, q2 and p1 are collinear and p1 lies on p2q2
            || (o4 == 0 && on_segment(p2, q1, q2)) // p2, q2 and q1 are collinear and q1 lies on p2q2
    }

    /// check if a polygon diagonal intersects any polygon segments
    fn intersects_polygon(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let mut p = a;
        loop {
            let p_next = node!(self.nodes, p.next_i);
            if (p.i != a.i && p.i != b.i && p_next.i != a.i && p_next.i != b.i)
                && self.intersects(p, p_next, a, b)
            {
                return true;
            }
            p = p_next;
            if ptr::eq(p, a) {
                return false;
            }
        }
    }

    /// check if a polygon diagonal is locally inside the polygon
    fn locally_inside(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let a_prev = node!(self.nodes, a.prev_i);
        let a_next = node!(self.nodes, a.next_i);
        if area(a_prev, a, a_next) < T::zero() {
            area(a, b, a_next) >= T::zero() && area(a, a_prev, b) >= T::zero()
        } else {
            area(a, b, a_prev) < T::zero() || area(a, a_next, b) < T::zero()
        }
    }

    /// check if the middle point of a polygon diagonal is inside the polygon
    fn middle_inside(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let mut p = a;
        let mut inside = false;
        let two = T::one() + T::one();
        let (px, py) = ((a.x + b.x) / two, (a.y + b.y) / two);
        loop {
            let p_next = node!(self.nodes, p.next_i);
            inside ^= (p.y > py) != (p_next.y > py)
                && p_next.y != p.y
                && (px < (p_next.x - p.x) * (py - p.y) / (p_next.y - p.y) + p.x);
            p = p_next;
            if ptr::eq(p, a) {
                return inside;
            }
        }
    }

    /// link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
    /// if one belongs to the outer ring and another to a hole, it merges it into a single ring
    fn split_polygon(&mut self, a_i: u32, b_i: u32) -> u32 {
        let a2_i = self.nodes.len() as u32;
        let b2_i = a2_i + 1;

        let a = node_mut!(self.nodes, a_i);
        let mut a2 = Node::new(a.i, a.x, a.y);
        let an_i = a.next_i;
        a.next_i = b_i;
        a2.prev_i = b2_i;
        a2.next_i = an_i;

        let b = node_mut!(self.nodes, b_i);
        let mut b2 = Node::new(b.i, b.x, b.y);
        let bp_i = b.prev_i;
        b.prev_i = a_i;
        b2.next_i = a2_i;
        b2.prev_i = bp_i;

        node_mut!(self.nodes, an_i).prev_i = a2_i;
        node_mut!(self.nodes, bp_i).next_i = b2_i;

        self.nodes.push(a2);
        self.nodes.push(b2);

        b2_i
    }
}

/// create a node and optionally link it with previous one (in a circular doubly linked list)
fn insert_node<T: Float>(nodes: &mut Vec<Node<T>>, i: u32, x: T, y: T, last: Option<u32>) -> u32 {
    let mut p = Node::new(i, x, y);
    let p_i = nodes.len() as u32;
    match last {
        Some(last_i) => {
            let last_next_i = node!(nodes, last_i).next_i;
            p.next_i = last_next_i;
            p.prev_i = last_i;
            node_mut!(nodes, last_next_i).prev_i = p_i;
            node_mut!(nodes, last_i).next_i = p_i;
        }
        None => {
            (p.prev_i, p.next_i) = (p_i, p_i);
        }
    }
    nodes.push(p);
    p_i
}

fn remove_node<T: Float>(nodes: &mut [Node<T>], p_i: u32) -> (u32, u32) {
    let p = node!(nodes, p_i);
    let p_next_i = p.next_i;
    let p_prev_i = p.prev_i;
    let p_next_z_i = p.next_z_i;
    let p_prev_z_i = p.prev_z_i;
    node_mut!(nodes, p_next_i).prev_i = p_prev_i;
    node_mut!(nodes, p_prev_i).next_i = p_next_i;
    if let Some(prev_z_i) = p_prev_z_i {
        node_mut!(nodes, prev_z_i).next_z_i = p_next_z_i;
    }
    if let Some(next_z_i) = p_next_z_i {
        node_mut!(nodes, next_z_i).prev_z_i = p_prev_z_i;
    }
    (p_prev_i, p_next_i)
}

/// Returns a percentage difference between the polygon area and its triangulation area;
/// used to verify correctness of triangulation
pub fn deviation<T: Float, N: Index>(
    data: impl IntoIterator<Item = T>,
    hole_indices: &[N],
    triangles: &[N],
) -> T {
    let data = data.into_iter().collect::<Vec<T>>();

    let has_holes = !hole_indices.is_empty();
    let outer_len = match has_holes {
        true => hole_indices[0].into_usize() * 2,
        false => data.len(),
    };
    let mut polygon_area = signed_area(&data, 0, outer_len as u32).abs();
    if has_holes {
        for i in 0..hole_indices.len() {
            let start = hole_indices[i].into_usize() * 2;
            let end = if i < hole_indices.len() - 1 {
                hole_indices[i + 1].into_usize() * 2
            } else {
                data.len()
            };
            polygon_area = polygon_area - signed_area(&data, start as u32, end as u32).abs();
        }
    }

    let mut triangles_area = T::zero();
    for [a, b, c] in triangles
        .chunks_exact(3)
        .map(|idxs| [idxs[0], idxs[1], idxs[2]])
    {
        let a = a.into_usize() * 2;
        let b = b.into_usize() * 2;
        let c = c.into_usize() * 2;
        triangles_area = triangles_area
            + ((data[a] - data[c]) * (data[b + 1] - data[a + 1])
                - (data[a] - data[b]) * (data[c + 1] - data[a + 1]))
                .abs();
    }
    if polygon_area == T::zero() && triangles_area == T::zero() {
        T::zero()
    } else {
        ((polygon_area - triangles_area) / polygon_area).abs()
    }
}

/// check if a point lies within a convex triangle
fn signed_area<T: Float>(data: &[T], start: u32, end: u32) -> T {
    if start == end {
        return T::zero();
    }
    let j = if end > 2 { end - 2 } else { 0 };
    let mut bx = data[j as usize];
    let mut by = data[(j + 1) as usize];
    let mut sum = T::zero();
    for a in data[start as usize..end as usize].chunks_exact(2) {
        let (ax, ay) = (a[0], a[1]);
        sum = sum + (bx - ax) * (ay + by);
        (bx, by) = (ax, ay);
    }
    sum
}

/// z-order of a point given coords and inverse of the longer side of data bbox
fn z_order<T: Float>(x: T, y: T, min_x: T, min_y: T, inv_size: T) -> i32 {
    // coords are transformed into non-negative 15-bit integer range
    let mut x = ((x - min_x) * inv_size).to_i32().unwrap();
    let mut y = ((y - min_y) * inv_size).to_i32().unwrap();
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;
    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;
    x | (y << 1)
}

#[allow(clippy::too_many_arguments)]
fn point_in_triangle<T: Float>(ax: T, ay: T, bx: T, by: T, cx: T, cy: T, px: T, py: T) -> bool {
    (cx - px) * (ay - py) >= (ax - px) * (cy - py)
        && (ax - px) * (by - py) >= (bx - px) * (ay - py)
        && (bx - px) * (cy - py) >= (cx - px) * (by - py)
}

/// signed area of a triangle
fn area<T: Float>(p: &Node<T>, q: &Node<T>, r: &Node<T>) -> T {
    (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
}

/// check if two points are equal
fn equals<T: Float>(p1: &Node<T>, p2: &Node<T>) -> bool {
    p1.x == p2.x && p1.y == p2.y
}

/// for collinear points p, q, r, check if point q lies on segment pr
fn on_segment<T: Float>(p: &Node<T>, q: &Node<T>, r: &Node<T>) -> bool {
    q.x <= p.x.max(r.x) && q.x >= p.x.min(r.x) && q.y <= p.y.max(r.y) && q.y >= p.y.max(r.y)
}

fn sign<T: Float>(v: T) -> i8 {
    (v > T::zero()) as i8 - (v < T::zero()) as i8
}
