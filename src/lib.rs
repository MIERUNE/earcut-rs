//! A Rust port of the [Earcut](https://github.com/mapbox/earcut) polygon triangulation library.

#![no_std]

extern crate alloc;

pub mod utils3d;

use alloc::vec::Vec;
use core::num::NonZeroU32;
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
        unsafe {
            debug_assert!($index.get() < $self.$nodes.len() as u32);
            $self.$nodes.get_unchecked($index.get() as usize)
        }
    };
    ($nodes:ident, $index:expr) => {
        unsafe {
            debug_assert!($index.get() < $nodes.len() as u32);
            $nodes.get_unchecked($index.get() as usize)
        }
    };
}

macro_rules! node_mut {
    ($self:ident.$nodes:ident, $index:expr) => {
        unsafe {
            debug_assert!($index.get() < $self.$nodes.len() as u32);
            $self.$nodes.get_unchecked_mut($index.get() as usize)
        }
    };
    ($nodes:ident, $index:expr) => {
        unsafe {
            debug_assert!($index.get() < $nodes.len() as u32);
            $nodes.get_unchecked_mut($index.get() as usize)
        }
    };
}

type NodeIndex = NonZeroU32;

struct Node<T: Float> {
    /// vertex index in coordinates array
    i: u32,
    /// z-order curve value
    z: i32,
    /// vertex coordinates x
    x: T,
    /// vertex coordinates y
    y: T,
    /// previous vertex nodes in a polygon ring
    prev_i: NodeIndex,
    /// next vertex nodes in a polygon ring
    next_i: NodeIndex,
    /// previous nodes in z-order
    prev_z_i: Option<NodeIndex>,
    /// next nodes in z-order
    next_z_i: Option<NodeIndex>,
    /// indicates whether this is a steiner point
    steiner: bool,
}

impl<T: Float> Node<T> {
    fn new(i: u32, x: T, y: T) -> Self {
        Self {
            i,
            x,
            y,
            prev_i: unsafe { NodeIndex::new_unchecked(1) },
            next_i: unsafe { NodeIndex::new_unchecked(1) },
            z: 0,
            prev_z_i: None,
            next_z_i: None,
            steiner: false,
        }
    }
}

/// Instance of the earcut algorithm.
pub struct Earcut<T: Float> {
    data: Vec<[T; 2]>,
    nodes: Vec<Node<T>>,
    queue: Vec<NodeIndex>,
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
        self.nodes.push(Node::new(0, T::infinity(), T::infinity())); // dummy node
    }

    /// Performs the earcut triangulation on a polygon.
    ///
    /// The API is similar to the original JavaScript implementation, except you can provide a vector for the output indices.
    pub fn earcut<N: Index>(
        &mut self,
        data: impl IntoIterator<Item = [T; 2]>,
        hole_indices: &[N],
        triangles_out: &mut Vec<N>,
    ) {
        self.data.clear();
        self.data.extend(data);
        triangles_out.clear();
        if self.data.len() < 3 {
            return;
        }
        self.earcut_impl(hole_indices, triangles_out);
    }

    pub fn earcut_impl<N: Index>(&mut self, hole_indices: &[N], triangles_out: &mut Vec<N>) {
        triangles_out.reserve(self.data.len() + 1);
        self.reset(self.data.len() / 2 * 3);

        let has_holes = !hole_indices.is_empty();
        let outer_len: usize = if has_holes {
            hole_indices[0].into_usize()
        } else {
            self.data.len()
        };

        // create nodes
        let Some(mut outer_node_i) = self.linked_list(0, outer_len, true) else {
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
        if self.data.len() > 80 {
            let [max_x, max_y] =
                self.data[1..outer_len]
                    .iter()
                    .fold(self.data[0], |[ax, ay], b| {
                        let (bx, by) = (b[0], b[1]);
                        [T::max(ax, bx), T::max(ay, by)]
                    });
            [min_x, min_y] = self.data[1..outer_len]
                .iter()
                .fold(self.data[0], |[ax, ay], b| {
                    let (bx, by) = (b[0], b[1]);
                    [T::min(ax, bx), T::min(ay, by)]
                });
            // minX, minY and invSize are later used to transform coords into integers for z-order calculation
            inv_size = (max_x - min_x).max(max_y - min_y);
            if inv_size != T::zero() {
                inv_size = T::from(32767.0).unwrap() / inv_size;
            }
        }

        earcut_linked(
            &mut self.nodes,
            outer_node_i,
            triangles_out,
            min_x,
            min_y,
            inv_size,
            Pass::P0,
        );
    }

    /// create a circular doubly linked list from polygon points in the specified winding order
    fn linked_list(&mut self, start: usize, end: usize, clockwise: bool) -> Option<NodeIndex> {
        let mut last_i: Option<NodeIndex> = None;
        let iter = self.data[start..end].iter().enumerate();

        if clockwise == (signed_area(&self.data, start, end) > T::zero()) {
            for (i, &[x, y]) in iter {
                let idx = start + i;
                last_i = Some(insert_node(&mut self.nodes, idx as u32, x, y, last_i));
            }
        } else {
            for (i, &[x, y]) in iter.rev() {
                let idx = start + i;
                last_i = Some(insert_node(&mut self.nodes, idx as u32, x, y, last_i));
            }
        };

        if let Some(li) = last_i {
            let last = node!(self.nodes, li);
            if equals(last, node!(self.nodes, last.next_i)) {
                let (_, next_i) = remove_node(&mut self.nodes, li);
                last_i = Some(next_i);
            }
        }

        last_i
    }

    /// link every hole into the outer loop, producing a single-ring polygon without holes
    fn eliminate_holes<N: Index>(
        &mut self,
        hole_indices: &[N],
        mut outer_node_i: NodeIndex,
    ) -> NodeIndex {
        self.queue.clear();
        for (i, hi) in hole_indices.iter().enumerate() {
            let start = (*hi).into_usize();
            let end = if i < hole_indices.len() - 1 {
                hole_indices[i + 1].into_usize()
            } else {
                self.data.len()
            };
            if let Some(list_i) = self.linked_list(start, end, false) {
                let list = &mut node_mut!(self.nodes, list_i);
                if list_i == list.next_i {
                    list.steiner = true;
                }
                self.queue.push(get_leftmost(&self.nodes, list_i))
            }
        }

        self.queue.sort_unstable_by(|a, b| {
            node!(self.nodes, *a)
                .x
                .partial_cmp(&node!(self.nodes, *b).x)
                .unwrap_or(core::cmp::Ordering::Equal)
        });

        // process holes from left to right
        for &q in &self.queue {
            outer_node_i = eliminate_hole(&mut self.nodes, q, outer_node_i);
        }

        outer_node_i
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Pass {
    P0 = 0,
    P1 = 1,
    P2 = 2,
}

/// main ear slicing loop which triangulates a polygon (given as a linked list)
#[allow(clippy::too_many_arguments)]
fn earcut_linked<T: Float, N: Index>(
    nodes: &mut Vec<Node<T>>,
    ear_i: NodeIndex,
    triangles: &mut Vec<N>,
    min_x: T,
    min_y: T,
    inv_size: T,
    pass: Pass,
) {
    let mut ear_i = ear_i;

    // interlink polygon nodes in z-order
    if pass == Pass::P0 && inv_size != T::zero() {
        index_curve(nodes, ear_i, min_x, min_y, inv_size);
    }

    let mut stop_i = ear_i;

    // iterate through ears, slicing them one by one
    loop {
        let ear = node!(nodes, ear_i);
        if ear.prev_i == ear.next_i {
            break;
        }
        let pi = ear.prev_i;
        let ni = ear.next_i;

        let is_ear = if inv_size != T::zero() {
            is_ear_hashed(nodes, ear_i, min_x, min_y, inv_size)
        } else {
            is_ear(nodes, ear_i)
        };
        if is_ear {
            let next = node!(nodes, ni);
            let next_i = next.i;
            let next_next_i = next.next_i;

            // cut off the triangle
            triangles.push(N::from_usize(node!(nodes, pi).i as usize));
            triangles.push(N::from_usize(ear.i as usize));
            triangles.push(N::from_usize(next_i as usize));

            remove_node(nodes, ear_i);

            // skipping the next vertex leads to less sliver triangles
            (ear_i, stop_i) = (next_next_i, next_next_i);

            continue;
        }

        ear_i = ni;

        // if we looped through the whole remaining polygon and can't find any more ears
        if ear_i == stop_i {
            if pass == Pass::P0 {
                // try filtering points and slicing again
                ear_i = filter_points(nodes, ear_i, None);
                earcut_linked(nodes, ear_i, triangles, min_x, min_y, inv_size, Pass::P1);
            } else if pass == Pass::P1 {
                // if this didn't work, try curing all small self-intersections locally
                let filtered = filter_points(nodes, ear_i, None);
                ear_i = cure_local_intersections(nodes, filtered, triangles);
                earcut_linked(nodes, ear_i, triangles, min_x, min_y, inv_size, Pass::P2);
            } else if pass == Pass::P2 {
                // as a last resort, try splitting the remaining polygon into two
                split_earcut(nodes, ear_i, triangles, min_x, min_y, inv_size);
            }
            return;
        }
    }
}

/// check whether a polygon node forms a valid ear with adjacent nodes
fn is_ear<T: Float>(nodes: &[Node<T>], ear_i: NodeIndex) -> bool {
    let b = node!(nodes, ear_i);
    let a = node!(nodes, b.prev_i);
    let c = node!(nodes, b.next_i);

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

    let mut p = node!(nodes, c.next_i);
    let mut p_prev = node!(nodes, p.prev_i);
    while !ptr::eq(p, a) {
        let p_next = node!(nodes, p.next_i);
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

fn is_ear_hashed<T: Float>(
    nodes: &[Node<T>],
    ear_i: NodeIndex,
    min_x: T,
    min_y: T,
    inv_size: T,
) -> bool {
    let b = node!(nodes, ear_i);
    let a = node!(nodes, b.prev_i);
    let c = node!(nodes, b.next_i);

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

    let ear = node!(nodes, ear_i);
    let mut o_p = ear.prev_z_i.map(|i| node!(nodes, i));
    let mut o_n = ear.next_z_i.map(|i| node!(nodes, i));

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
            && area(node!(nodes, p.prev_i), p, node!(nodes, p.next_i)) >= T::zero()
        {
            return false;
        }
        o_p = p.prev_z_i.map(|i| node!(nodes, i));

        if (n.x >= x0 && n.x <= x1 && n.y >= y0 && n.y <= y1)
            && (!ptr::eq(n, a) && !ptr::eq(n, c))
            && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y)
            && area(node!(nodes, n.prev_i), n, node!(nodes, n.next_i)) >= T::zero()
        {
            return false;
        }
        o_n = n.next_z_i.map(|i| node!(nodes, i));
    }

    // look for remaining points in decreasing z-order
    while let Some(p) = o_p {
        if p.z < min_z {
            break;
        };
        if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1)
            && (!ptr::eq(p, a) && !ptr::eq(p, c))
            && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y)
            && area(node!(nodes, p.prev_i), p, node!(nodes, p.next_i)) >= T::zero()
        {
            return false;
        }
        o_p = p.prev_z_i.map(|i| node!(nodes, i));
    }

    // look for remaining points in increasing z-order
    while let Some(n) = o_n {
        if n.z > max_z {
            break;
        };
        if (n.x >= x0 && n.x <= x1 && n.y >= y0 && n.y <= y1)
            && (!ptr::eq(n, a) && !ptr::eq(n, c))
            && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y)
            && area(node!(nodes, n.prev_i), n, node!(nodes, n.next_i)) >= T::zero()
        {
            return false;
        }
        o_n = n.next_z_i.map(|i| node!(nodes, i));
    }

    true
}

/// go through all polygon nodes and cure small local self-intersections
fn cure_local_intersections<T: Float, N: Index>(
    nodes: &mut [Node<T>],
    mut start_i: NodeIndex,
    triangles: &mut Vec<N>,
) -> NodeIndex {
    let mut p_i = start_i;
    loop {
        let p = node!(nodes, p_i);
        let p_next_i = p.next_i;
        let p_next = node!(nodes, p_next_i);
        let b_i = p_next.next_i;
        let a = node!(nodes, p.prev_i);
        let b = node!(nodes, b_i);

        if !equals(a, b)
            && intersects(a, p, p_next, b)
            && locally_inside(nodes, a, b)
            && locally_inside(nodes, b, a)
        {
            triangles.extend([
                N::from_usize(a.i as usize),
                N::from_usize(p.i as usize),
                N::from_usize(b.i as usize),
            ]);

            remove_node(nodes, p_i);
            remove_node(nodes, p_next_i);

            (p_i, start_i) = (b_i, b_i);
        }

        p_i = node!(nodes, p_i).next_i;
        if p_i == start_i {
            return filter_points(nodes, p_i, None);
        }
    }
}

/// try splitting polygon into two and triangulate them independently
fn split_earcut<T: Float, N: Index>(
    nodes: &mut Vec<Node<T>>,
    start_i: NodeIndex,
    triangles: &mut Vec<N>,
    min_x: T,
    min_y: T,
    inv_size: T,
) {
    // look for a valid diagonal that divides the polygon into two
    let mut ai = start_i;
    let mut a = node!(nodes, ai);
    loop {
        let a_next = node!(nodes, a.next_i);
        let a_prev = node!(nodes, a.prev_i);
        let mut bi = a_next.next_i;

        while bi != a.prev_i {
            let b = node!(nodes, bi);
            if a.i != b.i && is_valid_diagonal(nodes, a, b, a_next, a_prev) {
                // split the polygon in two by the diagonal
                let mut ci = split_polygon(nodes, ai, bi);

                // filter colinear points around the cuts
                let end_i = Some(node!(nodes, ai).next_i);
                ai = filter_points(nodes, ai, end_i);
                let end_i = Some(node!(nodes, ci).next_i);
                ci = filter_points(nodes, ci, end_i);

                // run earcut on each half
                earcut_linked(nodes, ai, triangles, min_x, min_y, inv_size, Pass::P0);
                earcut_linked(nodes, ci, triangles, min_x, min_y, inv_size, Pass::P0);
                return;
            }
            bi = b.next_i;
        }

        ai = a.next_i;
        if ai == start_i {
            return;
        }
        a = a_next;
    }
}

/// interlink polygon nodes in z-order
fn index_curve<T: Float>(
    nodes: &mut [Node<T>],
    start_i: NodeIndex,
    min_x: T,
    min_y: T,
    inv_size: T,
) {
    let mut p_i = start_i;
    let mut p = node_mut!(nodes, p_i);

    loop {
        if p.z == 0 {
            p.z = z_order(p.x, p.y, min_x, min_y, inv_size);
        }
        p.prev_z_i = Some(p.prev_i);
        p.next_z_i = Some(p.next_i);
        p_i = p.next_i;
        p = node_mut!(nodes, p_i);
        if p_i == start_i {
            break;
        }
    }

    let p_prev_z_i = p.prev_z_i.take().unwrap();
    node_mut!(nodes, p_prev_z_i).next_z_i = None;
    sort_linked(nodes, p_i);
}

/// Simon Tatham's linked list merge sort algorithm
/// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
fn sort_linked<T: Float>(nodes: &mut [Node<T>], list_i: NodeIndex) {
    let mut in_size: usize = 1;
    let mut list_i = Some(list_i);

    loop {
        let mut p_i = list_i;
        list_i = None;
        let mut tail_i: Option<NodeIndex> = None;
        let mut num_merges = 0;

        while let Some(p_i_s) = p_i {
            num_merges += 1;
            let mut q_i = node!(nodes, p_i_s).next_z_i;
            let mut p_size: u32 = 1;
            for _ in 1..in_size {
                if let Some(i) = q_i {
                    p_size += 1;
                    q_i = node!(nodes, i).next_z_i;
                } else {
                    break;
                }
            }
            let mut q_size = in_size;

            loop {
                let e_i = if p_size > 0 {
                    let Some(p_i_s) = p_i else { break };
                    if q_size > 0 {
                        if let Some(q_i_s) = q_i {
                            if node!(nodes, p_i_s).z <= node!(nodes, q_i_s).z {
                                p_size -= 1;
                                let e = node_mut!(nodes, p_i_s);
                                e.prev_z_i = tail_i;
                                p_i = e.next_z_i;
                                p_i_s
                            } else {
                                q_size -= 1;
                                let e = node_mut!(nodes, q_i_s);
                                e.prev_z_i = tail_i;
                                q_i = e.next_z_i;
                                q_i_s
                            }
                        } else {
                            p_size -= 1;
                            let e = node_mut!(nodes, p_i_s);
                            e.prev_z_i = tail_i;
                            p_i = e.next_z_i;
                            p_i_s
                        }
                    } else {
                        p_size -= 1;
                        let e = node_mut!(nodes, p_i_s);
                        e.prev_z_i = tail_i;
                        p_i = e.next_z_i;
                        p_i_s
                    }
                } else if q_size > 0 {
                    if let Some(q_i_s) = q_i {
                        q_size -= 1;
                        let e = node_mut!(nodes, q_i_s);
                        e.prev_z_i = tail_i;
                        q_i = e.next_z_i;
                        q_i_s
                    } else {
                        break;
                    }
                } else {
                    break;
                };

                if let Some(tail_i) = tail_i {
                    node_mut!(nodes, tail_i).next_z_i = Some(e_i);
                } else {
                    list_i = Some(e_i);
                }
                tail_i = Some(e_i);
            }

            p_i = q_i;
        }

        node_mut!(nodes, tail_i.unwrap()).next_z_i = None;
        if num_merges <= 1 {
            break;
        }
        in_size *= 2;
    }
}

/// find the leftmost node of a polygon ring
fn get_leftmost<T: Float>(nodes: &[Node<T>], start_i: NodeIndex) -> NodeIndex {
    let mut p_i = start_i;
    let mut p = node!(nodes, p_i);
    let mut leftmost_i = start_i;
    let mut leftmost = p;

    loop {
        if p.x < leftmost.x || (p.x == leftmost.x && p.y < leftmost.y) {
            (leftmost_i, leftmost) = (p_i, p);
        }
        p_i = p.next_i;
        if p_i == start_i {
            return leftmost_i;
        }
        p = node!(nodes, p_i);
    }
}

/// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
fn is_valid_diagonal<T: Float>(
    nodes: &[Node<T>],
    a: &Node<T>,
    b: &Node<T>,
    a_next: &Node<T>,
    a_prev: &Node<T>,
) -> bool {
    let b_next = node!(nodes, b.next_i);
    let b_prev = node!(nodes, b.prev_i);
    // dones't intersect other edges
    (a_next.i != b.i && a_prev.i != b.i && !intersects_polygon(nodes, a, b))
            // locally visible
            && ((locally_inside(nodes, a, b) && locally_inside(nodes, b, a) && middle_inside(nodes, a, b))
                // does not create opposite-facing sectors
                && (area(a_prev, a, b_prev) != T::zero() || area(a, b_prev, b) != T::zero())
                // special zero-length case
                || equals(a, b)
                    && area(a_prev, a, a_next) > T::zero()
                    && area(b_prev, b, b_next) > T::zero())
}

/// check if two segments intersect
fn intersects<T: Float>(p1: &Node<T>, q1: &Node<T>, p2: &Node<T>, q2: &Node<T>) -> bool {
    let o1 = sign(area(p1, q1, p2));
    let o2 = sign(area(p1, q1, q2));
    let o3 = sign(area(p2, q2, p1));
    let o4 = sign(area(p2, q2, q1));
    (o1 != o2 && o3 != o4) // general case
        || (o3 == 0 && on_segment(p2, p1, q2)) // p2, q2 and p1 are collinear and p1 lies on p2q2
        || (o4 == 0 && on_segment(p2, q1, q2)) // p2, q2 and q1 are collinear and q1 lies on p2q2
        || (o2 == 0 && on_segment(p1, q2, q1)) // p1, q1 and q2 are collinear and q2 lies on p1q1
        || (o1 == 0 && on_segment(p1, p2, q1)) // p1, q1 and p2 are collinear and p2 lies on p1q1
}

/// check if a polygon diagonal intersects any polygon segments
fn intersects_polygon<T: Float>(nodes: &[Node<T>], a: &Node<T>, b: &Node<T>) -> bool {
    let mut p = a;
    loop {
        let p_next = node!(nodes, p.next_i);
        if (p.i != a.i && p.i != b.i && p_next.i != a.i && p_next.i != b.i)
            && intersects(p, p_next, a, b)
        {
            return true;
        }
        p = p_next;
        if ptr::eq(p, a) {
            return false;
        }
    }
}

/// check if the middle point of a polygon diagonal is inside the polygon
fn middle_inside<T: Float>(nodes: &[Node<T>], a: &Node<T>, b: &Node<T>) -> bool {
    let mut p = a;
    let mut inside = false;
    let two = T::one() + T::one();
    let (px, py) = ((a.x + b.x) / two, (a.y + b.y) / two);
    loop {
        let p_next = node!(nodes, p.next_i);
        inside ^= (p.y > py) != (p_next.y > py)
            && p_next.y != p.y
            && (px < (p_next.x - p.x) * (py - p.y) / (p_next.y - p.y) + p.x);
        p = p_next;
        if ptr::eq(p, a) {
            return inside;
        }
    }
}

/// find a bridge between vertices that connects hole with an outer ring and and link it
fn eliminate_hole<T: Float>(
    nodes: &mut Vec<Node<T>>,
    hole_i: NodeIndex,
    outer_node_i: NodeIndex,
) -> NodeIndex {
    let Some(bridge_i) = find_hole_bridge(nodes, node!(nodes, hole_i), outer_node_i) else {
        return outer_node_i;
    };
    let bridge_reverse_i = split_polygon(nodes, bridge_i, hole_i);

    // filter collinear points around the cuts
    let end_i = Some(node!(nodes, bridge_reverse_i).next_i);
    filter_points(nodes, bridge_reverse_i, end_i);
    let end_i = Some(node!(nodes, bridge_i).next_i);
    filter_points(nodes, bridge_i, end_i)
}

/// check if a polygon diagonal is locally inside the polygon
fn locally_inside<T: Float>(nodes: &[Node<T>], a: &Node<T>, b: &Node<T>) -> bool {
    let a_prev = node!(nodes, a.prev_i);
    let a_next = node!(nodes, a.next_i);
    if area(a_prev, a, a_next) < T::zero() {
        area(a, b, a_next) >= T::zero() && area(a, a_prev, b) >= T::zero()
    } else {
        area(a, b, a_prev) < T::zero() || area(a, a_next, b) < T::zero()
    }
}

/// David Eberly's algorithm for finding a bridge between hole and outer polygon
fn find_hole_bridge<T: Float>(
    nodes: &[Node<T>],
    hole: &Node<T>,
    outer_node_i: NodeIndex,
) -> Option<NodeIndex> {
    let mut p_i = outer_node_i;
    let mut qx = T::neg_infinity();
    let mut m_i: Option<NodeIndex> = None;

    // find a segment intersected by a ray from the hole's leftmost point to the left;
    // segment's endpoint with lesser x will be potential connection point
    let mut p = node!(nodes, p_i);
    loop {
        let p_next = node!(nodes, p.next_i);
        if hole.y <= p.y && hole.y >= p_next.y && p_next.y != p.y {
            let x = p.x + (hole.y - p.y) * (p_next.x - p.x) / (p_next.y - p.y);
            if x <= hole.x && x > qx {
                qx = x;
                m_i = Some(if p.x < p_next.x { p_i } else { p.next_i });
                if x == hole.x {
                    // hole touches outer segment; pick leftmost endpoint
                    return m_i;
                }
            }
        }
        p_i = p.next_i;
        if p_i == outer_node_i {
            break;
        }
        p = p_next;
    }

    let mut m_i = m_i?;

    // look for points inside the triangle of hole point, segment intersection and endpoint;
    // if there are no points found, we have a valid connection;
    // otherwise choose the point of the minimum angle with the ray as connection point

    let stop_i = m_i;
    let Node { x: mx, y: my, .. } = *node!(nodes, m_i); // must copy
    let mut tan_min = T::infinity();

    p_i = m_i;
    let mut p = node!(nodes, p_i);
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
            if locally_inside(nodes, p, hole)
                && (tan < tan_min
                    || (tan == tan_min
                        && (p.x > m.x || (p.x == m.x && sector_contains_sector(nodes, m, p)))))
            {
                (m_i, m) = (p_i, p);
                tan_min = tan;
            }
        }

        p_i = p.next_i;
        if p_i == stop_i {
            return Some(m_i);
        }
        p = node!(nodes, p_i);
    }
}

/// whether sector in vertex m contains sector in vertex p in the same coordinates
fn sector_contains_sector<T: Float>(nodes: &[Node<T>], m: &Node<T>, p: &Node<T>) -> bool {
    area(node!(nodes, m.prev_i), m, node!(nodes, p.prev_i)) < T::zero()
        && area(node!(nodes, p.next_i), m, node!(nodes, m.next_i)) < T::zero()
}

/// eliminate colinear or duplicate points
fn filter_points<T: Float>(
    nodes: &mut [Node<T>],
    start_i: NodeIndex,
    end_i: Option<NodeIndex>,
) -> NodeIndex {
    let mut end_i = end_i.unwrap_or(start_i);

    let mut p_i = start_i;
    let mut p = node!(nodes, p_i);
    loop {
        let p_next = node!(nodes, p.next_i);
        if !p.steiner && (equals(p, p_next) || area(node!(nodes, p.prev_i), p, p_next) == T::zero())
        {
            let (prev_i, next_i) = remove_node(nodes, p_i);
            (p_i, end_i) = (prev_i, prev_i);
            if p_i == next_i {
                return end_i;
            }
            p = node!(nodes, p_i);
        } else {
            p_i = p.next_i;
            if p_i == end_i {
                return end_i;
            }
            p = p_next;
        };
    }
}

/// link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
/// if one belongs to the outer ring and another to a hole, it merges it into a single ring
fn split_polygon<T: Float>(nodes: &mut Vec<Node<T>>, a_i: NodeIndex, b_i: NodeIndex) -> NodeIndex {
    debug_assert!(!nodes.is_empty());
    let a2_i = unsafe { NodeIndex::new_unchecked(nodes.len() as u32) };
    let b2_i = unsafe { NodeIndex::new_unchecked(nodes.len() as u32 + 1) };

    let a = node_mut!(nodes, a_i);
    let mut a2 = Node::new(a.i, a.x, a.y);
    let an_i = a.next_i;
    a.next_i = b_i;
    a2.prev_i = b2_i;
    a2.next_i = an_i;
    node_mut!(nodes, an_i).prev_i = a2_i;

    let b = node_mut!(nodes, b_i);
    let mut b2 = Node::new(b.i, b.x, b.y);
    let bp_i = b.prev_i;
    b.prev_i = a_i;
    b2.next_i = a2_i;
    b2.prev_i = bp_i;
    node_mut!(nodes, bp_i).next_i = b2_i;

    nodes.extend([a2, b2]);

    b2_i
}

/// create a node and optionally link it with previous one (in a circular doubly linked list)
fn insert_node<T: Float>(
    nodes: &mut Vec<Node<T>>,
    i: u32,
    x: T,
    y: T,
    last: Option<NodeIndex>,
) -> NodeIndex {
    let mut p = Node::new(i, x, y);
    let p_i = unsafe { NodeIndex::new_unchecked(nodes.len() as u32) };
    match last {
        Some(last_i) => {
            let last = node_mut!(nodes, last_i);
            let last_next_i = last.next_i;
            (p.next_i, last.next_i) = (last_next_i, p_i);
            p.prev_i = last_i;
            node_mut!(nodes, last_next_i).prev_i = p_i;
        }
        None => {
            (p.prev_i, p.next_i) = (p_i, p_i);
        }
    }
    nodes.push(p);
    p_i
}

fn remove_node<T: Float>(nodes: &mut [Node<T>], p_i: NodeIndex) -> (NodeIndex, NodeIndex) {
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
    data: impl IntoIterator<Item = [T; 2]>,
    hole_indices: &[N],
    triangles: &[N],
) -> T {
    let data = data.into_iter().collect::<Vec<[T; 2]>>();
    let has_holes = !hole_indices.is_empty();
    let outer_len = match has_holes {
        true => hole_indices[0].into_usize(),
        false => data.len(),
    };
    let polygon_area = if data.len() < 3 {
        T::zero()
    } else {
        let mut polygon_area = signed_area(&data, 0, outer_len).abs();
        if has_holes {
            for i in 0..hole_indices.len() {
                let start = hole_indices[i].into_usize();
                let end = if i < hole_indices.len() - 1 {
                    hole_indices[i + 1].into_usize()
                } else {
                    data.len()
                };
                if end - start >= 3 {
                    polygon_area = polygon_area - signed_area(&data, start, end).abs();
                }
            }
        }
        polygon_area
    };

    let mut triangles_area = T::zero();
    for [a, b, c] in triangles
        .chunks_exact(3)
        .map(|idxs| [idxs[0], idxs[1], idxs[2]])
    {
        let a = a.into_usize();
        let b = b.into_usize();
        let c = c.into_usize();
        triangles_area = triangles_area
            + ((data[a][0] - data[c][0]) * (data[b][1] - data[a][1])
                - (data[a][0] - data[b][0]) * (data[c][1] - data[a][1]))
                .abs();
    }
    if polygon_area == T::zero() && triangles_area == T::zero() {
        T::zero()
    } else {
        ((polygon_area - triangles_area) / polygon_area).abs()
    }
}

/// check if a point lies within a convex triangle
fn signed_area<T: Float>(data: &[[T; 2]], start: usize, end: usize) -> T {
    let [mut bx, mut by] = data[end - 1];
    let mut sum = T::zero();
    for &[ax, ay] in &data[start..end] {
        sum = sum + (bx - ax) * (ay + by);
        (bx, by) = (ax, ay);
    }
    sum
}

/// z-order of a point given coords and inverse of the longer side of data bbox
fn z_order<T: Float>(x: T, y: T, min_x: T, min_y: T, inv_size: T) -> i32 {
    // coords are transformed into non-negative 15-bit integer range
    let x = ((x - min_x) * inv_size).to_u32().unwrap();
    let y = ((y - min_y) * inv_size).to_u32().unwrap();
    let mut xy = (x as i64) << 32 | y as i64;
    xy = (xy | (xy << 8)) & 0x00FF00FF00FF00FF;
    xy = (xy | (xy << 4)) & 0x0F0F0F0F0F0F0F0F;
    xy = (xy | (xy << 2)) & 0x3333333333333333;
    xy = (xy | (xy << 1)) & 0x5555555555555555;
    (xy >> 32 | xy << 1) as i32
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

fn sign<T: Float>(v: T) -> i32 {
    (v > T::zero()) as i32 - (v < T::zero()) as i32
}
