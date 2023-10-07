use num_traits::float::Float;

pub trait Index: Copy {
    fn into(self) -> usize;
    fn from(v: usize) -> Self;
}
impl Index for u32 {
    fn into(self) -> usize {
        self as usize
    }
    fn from(v: usize) -> Self {
        v as Self
    }
}
impl Index for u16 {
    fn into(self) -> usize {
        self as usize
    }
    fn from(v: usize) -> Self {
        v as Self
    }
}
impl Index for usize {
    fn into(self) -> usize {
        self as usize
    }
    fn from(v: usize) -> Self {
        v as Self
    }
}

macro_rules! node {
    ($self:expr, $index:expr) => {
        $self.nodes[unsafe { $index.unwrap_unchecked() }]
    };
}

#[derive(Debug)]
struct Node<T: Float> {
    i: usize,
    x: T,
    y: T,
    prev_i: Option<usize>,
    next_i: Option<usize>,
    z: Option<u32>,
    prev_z_i: Option<usize>,
    next_z_i: Option<usize>,
    steiner: bool,
}

impl<T: Float> Node<T> {
    fn new(i: usize, x: T, y: T) -> Self {
        Self {
            i,
            x,
            y,
            prev_i: None,
            next_i: None,
            z: None,
            prev_z_i: None,
            next_z_i: None,
            steiner: false,
        }
    }
}
pub struct Earcut<T: Float> {
    nodes: Vec<Node<T>>,
}

fn sign<T: Float>(num: T) -> i32 {
    if num > T::zero() {
        1
    } else if num < T::zero() {
        -1
    } else {
        0
    }
}

impl<T: Float> Earcut<T> {
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    fn reset(&mut self, capacity: usize) {
        self.nodes.clear();
        self.nodes.reserve(capacity);
    }

    pub fn earcut<N: Index>(
        &mut self,
        data: &[T],
        hole_indices: &[N],
        dim: usize,
        triangles_out: &mut Vec<N>,
    ) {
        triangles_out.clear();
        self.reset(data.len() / dim * 3 / 2);

        let has_holes = hole_indices.len() > 0;
        let outer_len: usize = if has_holes {
            hole_indices[0].into() * dim
        } else {
            data.len()
        };
        let Some(mut outer_node_i) = self.linked_list(data, 0, outer_len, dim, true) else {
            return;
        };

        let outer_node = &self.nodes[outer_node_i];
        if outer_node.next_i == outer_node.prev_i {
            return;
        }

        if has_holes {
            outer_node_i = self.eliminate_holes(data, hole_indices, outer_node_i, dim);
        }

        let mut min_x = T::zero();
        let mut min_y = T::zero();
        let mut inv_size = T::zero();
        if data.len() > 80 * dim {
            let max_x = data[0..outer_len]
                .iter()
                .step_by(dim)
                .fold(data[0], |a, b| T::max(a, *b));
            min_x = data[dim..outer_len]
                .iter()
                .step_by(dim)
                .fold(data[0], |a, b| T::min(a, *b));
            let max_y = data[dim + 1..outer_len]
                .iter()
                .step_by(dim)
                .fold(data[1], |a, b| T::max(a, *b));
            min_y = data[dim + 1..outer_len]
                .iter()
                .step_by(dim)
                .fold(data[1], |a, b| T::min(a, *b));

            inv_size = (max_x - min_x).max(max_y - min_y);
            if inv_size != T::zero() {
                inv_size = T::from(32767.0).unwrap() / inv_size;
            }
        }

        self.earcut_linked(
            Some(outer_node_i),
            triangles_out,
            dim,
            min_x,
            min_y,
            inv_size,
            0,
        );
    }

    fn linked_list(
        &mut self,
        data: &[T],
        start: usize,
        end: usize,
        dim: usize,
        clockwise: bool,
    ) -> Option<usize> {
        let mut last_i: Option<usize> = None;

        if clockwise == (signed_area(data, start, end, dim) > T::zero()) {
            for (i, w) in data[start..end].windows(2).enumerate().step_by(dim) {
                if let [x, y] = *w {
                    last_i = Some(self.insert_node((start + i) as usize, x, y, last_i));
                }
            }
        } else {
            for (i, w) in data[start..end].windows(2).enumerate().step_by(dim).rev() {
                if let [x, y] = *w {
                    last_i = Some(self.insert_node((start + i) as usize, x, y, last_i));
                }
            }
        };

        if let Some(li) = last_i {
            let last = &self.nodes[li];
            if equals(last, &node!(self, last.next_i)) {
                self.remove_node(li);
                last_i = self.nodes[li].next_i;
            }
        }

        return last_i;
    }

    fn filter_points(&mut self, start_i: usize, end_i: Option<usize>) -> usize {
        let mut end_i = end_i.unwrap_or_else(|| start_i);

        let mut p_i = start_i;
        loop {
            let p = &self.nodes[p_i];
            let p_next = &node!(self, p.next_i);
            let again = if !p.steiner
                && (equals(p, &p_next) || area(&node!(self, p.prev_i), p, &p_next).is_zero())
            {
                self.remove_node(p_i);
                let p = &self.nodes[p_i];
                end_i = p.prev_i.unwrap();
                p_i = end_i;
                if p_i == p.next_i.unwrap() {
                    break;
                }
                true
            } else {
                p_i = p.next_i.unwrap();
                false
            };

            if !again && p_i == end_i {
                break;
            }
        }

        return end_i;
    }

    fn earcut_linked<N: Index>(
        &mut self,
        ear_i: Option<usize>,
        triangles: &mut Vec<N>,
        dim: usize,
        min_x: T,
        min_y: T,
        inv_size: T,
        pass: usize,
    ) {
        let Some(mut ear_i) = ear_i else { return };

        if pass == 0 && inv_size != T::zero() {
            self.index_curve(ear_i, min_x, min_y, inv_size);
        }

        let mut stop_i = ear_i;

        while self.nodes[ear_i].prev_i != self.nodes[ear_i].next_i {
            let prev_i = self.nodes[ear_i].prev_i.unwrap();
            let next_i = self.nodes[ear_i].next_i.unwrap();
            let is_ear = if inv_size != T::zero() {
                self.is_ear_hashed(ear_i, min_x, min_y, inv_size)
            } else {
                self.is_ear(ear_i)
            };

            if is_ear {
                let ear = &self.nodes[ear_i];
                triangles.push(N::from(self.nodes[prev_i].i / dim));
                triangles.push(N::from(ear.i / dim));
                triangles.push(N::from(self.nodes[next_i].i / dim));

                self.remove_node(ear_i);

                let next_next_i = self.nodes[next_i].next_i.unwrap();
                ear_i = next_next_i;
                stop_i = next_next_i;

                continue;
            }

            ear_i = next_i;

            if ear_i == stop_i {
                if pass == 0 {
                    let eee_i = Some(self.filter_points(ear_i, None));
                    self.earcut_linked(eee_i, triangles, dim, min_x, min_y, inv_size, 1);
                } else if pass == 1 {
                    let filtered = self.filter_points(ear_i, None);
                    let ear_i = self.cure_local_intersections(filtered, triangles, dim);
                    self.earcut_linked(Some(ear_i), triangles, dim, min_x, min_y, inv_size, 2);
                } else if pass == 2 {
                    self.split_earcut(ear_i, triangles, dim, min_x, min_y, inv_size);
                }

                break;
            }
        }
    }

    fn is_ear(&self, ear_i: usize) -> bool {
        let b = &self.nodes[ear_i];
        let a = &node!(self, b.prev_i);
        let c = &node!(self, b.next_i);
        if area(a, b, c) >= T::zero() {
            return false;
        }

        let x0 = a.x.min(b.x.min(c.x));
        let y0 = a.y.min(b.y.min(c.y));
        let x1 = a.x.max(b.x.max(c.x));
        let y1 = a.y.max(b.y.max(c.y));

        let mut p = &node!(self, c.next_i);
        while !std::ptr::eq(p, a) {
            if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1)
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y)
                && area(&node!(self, p.prev_i), p, &node!(self, p.next_i)) >= T::zero()
            {
                return false;
            }
            p = &node!(self, p.next_i);
        }
        return true;
    }

    fn is_ear_hashed(&self, ear_i: usize, min_x: T, min_y: T, inv_size: T) -> bool {
        let b = &self.nodes[ear_i];
        let a = &node!(self, b.prev_i);
        let c = &node!(self, b.next_i);
        if area(a, b, c) >= T::zero() {
            return false;
        }

        let x0 = a.x.min(b.x.min(c.x));
        let y0 = a.y.min(b.y.min(c.y));
        let x1 = a.x.max(b.x.max(c.x));
        let y1 = a.y.max(b.y.max(c.y));

        let min_z = z_order(x0, y0, min_x, min_y, inv_size);
        let max_z = z_order(x1, y1, min_x, min_y, inv_size);

        let ear = &self.nodes[ear_i];
        let mut p_i = ear.prev_z_i;
        let mut n_i = ear.next_z_i;

        let compare_p = |p_i: Option<usize>| {
            if let Some(p) = p_i.and_then(|p_i| Some(&self.nodes[p_i])) {
                return p.z.unwrap() >= min_z;
            }
            false
        };
        let compare_n = |n_i: Option<usize>| {
            if let Some(n) = n_i {
                return self.nodes[n].z.unwrap() <= max_z;
            }
            false
        };

        while compare_p(p_i) && compare_n(n_i) {
            let p = &node!(self, p_i);
            if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1)
                && !std::ptr::eq(p, a)
                && !std::ptr::eq(p, c)
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y)
                && area(&node!(self, p.prev_i), p, &node!(self, p.next_i)) >= T::zero()
            {
                return false;
            }
            p_i = p.prev_z_i;

            let n = &node!(self, n_i);
            if (n.x >= x0 && n.x <= x1 && n.y >= y0 && n.y <= y1)
                && !std::ptr::eq(n, a)
                && !std::ptr::eq(n, c)
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y)
                && area(&node!(self, n.prev_i), n, &node!(self, n.next_i)) >= T::zero()
            {
                return false;
            }
            n_i = n.next_z_i;
        }

        while compare_p(p_i) {
            let p = &node!(self, p_i);
            if p_i != ear.prev_i
                && p_i != ear.next_i
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y)
                && area(&node!(self, p.prev_i), p, &node!(self, p.next_i)) >= T::zero()
            {
                return false;
            }
            p_i = p.prev_z_i;
        }

        while compare_n(n_i) {
            let n = &node!(self, n_i);
            if n_i != ear.prev_i
                && n_i != ear.next_i
                && point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y)
                && area(&node!(self, n.prev_i), n, &node!(self, n.next_i)) >= T::zero()
            {
                return false;
            }
            n_i = n.next_z_i;
        }

        return true;
    }

    fn cure_local_intersections<N: Index>(
        &mut self,
        mut start_i: usize,
        triangles: &mut Vec<N>,
        dim: usize,
    ) -> usize {
        let mut p_i = start_i;
        loop {
            let Node {
                prev_i: p_prev_i,
                next_i: p_next_i,
                ..
            } = self.nodes[p_i];
            let a_i = p_prev_i.unwrap();
            let p_next = &node!(self, p_next_i);
            let b_i = p_next.next_i.unwrap();
            let a = &self.nodes[a_i];
            let b = &self.nodes[b_i];
            let p = &self.nodes[p_i];

            if !equals(a, b)
                && self.intersects(a, p, p_next, b)
                && self.locally_inside(a, b)
                && self.locally_inside(b, a)
            {
                triangles.push(N::from(a.i / dim));
                triangles.push(N::from(p.i / dim));
                triangles.push(N::from(b.i / dim));

                self.remove_node(p_i);
                self.remove_node(p_next_i.unwrap());

                start_i = b_i;
                p_i = start_i;
            }

            p_i = self.nodes[p_i].next_i.unwrap();
            if p_i == start_i {
                break;
            }
        }

        return self.filter_points(p_i, None);
    }

    fn split_earcut<N: Index>(
        &mut self,
        start_i: usize,
        triangles: &mut Vec<N>,
        dim: usize,
        min_x: T,
        min_y: T,
        inv_size: T,
    ) {
        let mut a_i = start_i;
        loop {
            let Node {
                i: ai,
                prev_i: a_prev_i,
                ..
            } = self.nodes[a_i];
            let mut b_i = (&node!(self, self.nodes[a_i].next_i)).next_i.unwrap();

            while b_i != a_prev_i.unwrap() {
                let b = &self.nodes[b_i];
                if ai != b.i && self.is_valid_diagonal(&self.nodes[a_i], b) {
                    let mut c_i = self.split_polygon(a_i, b_i);
                    a_i = self.filter_points(a_i, Some(self.nodes[a_i].next_i.unwrap()));
                    c_i = self.filter_points(c_i, self.nodes[c_i].next_i);

                    self.earcut_linked(Some(a_i), triangles, dim, min_x, min_y, inv_size, 0);
                    self.earcut_linked(Some(c_i), triangles, dim, min_x, min_y, inv_size, 0);
                    return;
                }
                b_i = b.next_i.unwrap();
            }

            a_i = self.nodes[a_i].next_i.unwrap();
            if a_i == start_i {
                break;
            }
        }
    }

    fn eliminate_holes<N: Index>(
        &mut self,
        data: &[T],
        hole_indices: &[N],
        mut outer_node_i: usize,
        dim: usize,
    ) -> usize {
        let mut queue = Vec::new();
        let len = hole_indices.len();

        for (i, hi) in hole_indices.iter().enumerate() {
            let start = (*hi).into() as usize * dim;
            let end = if i < len - 1 {
                hole_indices[i + 1].into() * dim
            } else {
                data.len()
            };
            if let Some(list_i) = self.linked_list(data, start, end, dim, false) {
                if list_i == self.nodes[list_i].next_i.unwrap() {
                    self.nodes[list_i].steiner = true;
                }
                queue.push(self.get_left_most(list_i))
            }
        }

        queue.sort_by(|a, b| self.nodes[*a].x.partial_cmp(&self.nodes[*b].x).unwrap());

        for q_i in queue.into_iter() {
            outer_node_i = self.eliminate_hole(q_i, outer_node_i);
        }

        return outer_node_i;
    }

    fn eliminate_hole(&mut self, hole_i: usize, outer_node_i: usize) -> usize {
        let Some(bridge_i) = self.find_hole_bridge(&self.nodes[hole_i], outer_node_i) else {
            return outer_node_i;
        };
        let bridge_reverse_i = self.split_polygon(bridge_i, hole_i);
        self.filter_points(bridge_reverse_i, self.nodes[bridge_reverse_i].next_i);
        return self.filter_points(bridge_i, self.nodes[bridge_i].next_i);
    }

    fn find_hole_bridge(&self, hole: &Node<T>, outer_node_i: usize) -> Option<usize> {
        let mut p_i = outer_node_i;
        let mut qx = T::neg_infinity();
        let mut m_i: Option<usize> = None;

        loop {
            let p = &self.nodes[p_i];
            let p_next_i = p.next_i.unwrap();
            let p_next = &self.nodes[p_next_i];
            if hole.y <= p.y && hole.y >= p_next.y && p_next.y != p.y {
                let x = p.x + (hole.y - p.y) * (p_next.x - p.x) / (p_next.y - p.y);
                if x <= hole.x && x > qx {
                    qx = x;
                    m_i = Some(if p.x < p_next.x { p_i } else { p_next_i });
                    if x == hole.x {
                        return m_i;
                    }
                }
            }
            p_i = p_next_i;
            if p_i == outer_node_i {
                break;
            }
        }

        let Some(mut m_i) = m_i else {
            return None;
        };
        if hole.x == qx {
            return Some(m_i);
        }
        let stop_i = m_i;
        let Node { x: mx, y: my, .. } = self.nodes[m_i]; // must copy
        let mut tan_min = T::infinity();

        p_i = m_i;

        loop {
            let p = &self.nodes[p_i];

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
                let m = &self.nodes[m_i];
                if self.locally_inside(p, hole)
                    && (tan < tan_min
                        || (tan == tan_min
                            && (p.x > m.x || p.x == m.x && self.sector_contains_sector(m, p))))
                {
                    m_i = p_i;
                    tan_min = tan;
                }
            }

            p_i = p.next_i.unwrap();
            if p_i == stop_i {
                break;
            }
        }

        return Some(m_i);
    }

    fn sector_contains_sector(&self, m: &Node<T>, p: &Node<T>) -> bool {
        return area(&node!(self, m.prev_i), m, &node!(self, p.prev_i)) < T::zero()
            && area(&node!(self, p.next_i), m, &node!(self, m.next_i)) < T::zero();
    }

    fn index_curve(&mut self, start_i: usize, min_x: T, min_y: T, inv_size: T) {
        let mut p_i = start_i;
        loop {
            let p = &mut self.nodes[p_i];
            if p.z.is_none() {
                p.z = Some(z_order(p.x, p.y, min_x, min_y, inv_size));
            }
            p.prev_z_i = p.prev_i;
            p.next_z_i = p.next_i;
            p_i = p.next_i.unwrap();
            if p_i == start_i {
                break;
            }
        }

        let p_prev_z_i = self.nodes[p_i].prev_z_i.unwrap();
        self.nodes[p_prev_z_i].next_z_i = None;
        self.nodes[p_i].prev_z_i = None;

        self.sort_linked(p_i);
    }

    fn sort_linked(&mut self, _list: usize) {
        let mut in_size = 1;
        let mut _list = Option::Some(_list);

        loop {
            let mut p_i = _list;
            _list = None;
            let mut tail_i: Option<usize> = None;
            let mut num_merges = 0;

            while p_i.is_some() {
                num_merges += 1;
                let mut q_i = p_i;
                let mut p_size = 0;
                for _ in 0..in_size {
                    p_size += 1;
                    q_i = node!(self, q_i).next_z_i;
                    if q_i.is_none() {
                        break;
                    }
                }
                let mut q_size = in_size;

                while p_size > 0 || (q_size > 0 && q_i.is_some()) {
                    let e_i = if p_size != 0
                        && (q_size == 0
                            || q_i.is_none()
                            || &node!(self, p_i).z.unwrap() <= &node!(self, q_i).z.unwrap())
                    {
                        let e_i = p_i.unwrap();
                        p_i = node!(self, p_i).next_z_i;
                        p_size -= 1;
                        e_i
                    } else {
                        let e_i = q_i.unwrap();
                        q_i = node!(self, q_i).next_z_i;
                        q_size -= 1;
                        e_i
                    };

                    if tail_i.is_some() {
                        (&mut node!(self, tail_i)).next_z_i = Some(e_i);
                    } else {
                        _list = Some(e_i);
                    }

                    self.nodes[e_i].prev_z_i = tail_i;
                    tail_i = Some(e_i);
                }

                p_i = q_i;
            }

            (&mut node!(self, tail_i)).next_z_i = None;
            in_size *= 2;

            if num_merges <= 1 {
                break;
            }
        }
    }

    fn get_left_most(&self, start: usize) -> usize {
        let mut p_i = start;
        let mut leftmost_i = start;

        loop {
            let p = &self.nodes[p_i];
            let leftmost = &self.nodes[leftmost_i];
            if p.x < leftmost.x || (p.x == leftmost.x && p.y < leftmost.y) {
                leftmost_i = p_i
            }
            p_i = p.next_i.unwrap();
            if p_i == start {
                break;
            }
        }
        leftmost_i
    }

    fn is_valid_diagonal(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let a_next = &node!(self, a.next_i);
        let a_prev = &node!(self, a.prev_i);
        let b_next = &node!(self, b.next_i);
        let b_prev = &node!(self, b.prev_i);
        (a_next.i != b.i && a_prev.i != b.i)
            && !self.intersects_polygon(a, b)
            && (self.locally_inside(a, b)
                && self.locally_inside(b, a)
                && self.middle_inside(a, b)
                && (area(a_prev, a, b_prev) != T::zero() || area(a, b_prev, b) != T::zero())
                || equals(a, b)
                    && area(a_prev, a, a_next) > T::zero()
                    && area(b_prev, b, b_next) > T::zero())
    }

    fn intersects(&self, p1: &Node<T>, q1: &Node<T>, p2: &Node<T>, q2: &Node<T>) -> bool {
        let o1 = sign(area(p1, q1, p2));
        let o2 = sign(area(p1, q1, q2));
        let o3 = sign(area(p2, q2, p1));
        let o4 = sign(area(p2, q2, q1));
        (o1 != o2 && o3 != o4)
            || (o1 == 0 && on_segment(p1, p2, q1))
            || (o2 == 0 && on_segment(p1, q2, q1))
            || (o3 == 0 && on_segment(p2, p1, q2))
            || (o4 == 0 && on_segment(p2, q1, q2))
    }

    fn intersects_polygon(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let mut p = a;
        loop {
            let p_next_i = p.next_i.unwrap();
            let p_next = &self.nodes[p_next_i];
            if (p.i != a.i && p_next.i != a.i && p.i != b.i && p_next.i != b.i)
                && self.intersects(p, p_next, a, b)
            {
                return true;
            }

            p = &self.nodes[p_next_i];
            if std::ptr::eq(p, a) {
                break;
            }
        }
        false
    }

    #[inline]
    fn locally_inside(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let a_prev = &node!(self, a.prev_i);
        let a_next = &node!(self, a.next_i);
        if area(a_prev, a, a_next) < T::zero() {
            area(a, b, a_next) >= T::zero() && area(a, a_prev, b) >= T::zero()
        } else {
            area(a, b, a_prev) < T::zero() || area(a, a_next, b) < T::zero()
        }
    }

    fn middle_inside(&self, a: &Node<T>, b: &Node<T>) -> bool {
        let mut p = a;
        let mut inside = false;
        let two = T::from(2).unwrap();
        let px = (a.x + b.x) / two;
        let py = (a.y + b.y) / two;
        loop {
            let p_next = &node!(self, p.next_i);
            inside ^= (p.y > py) != (p_next.y > py)
                && p_next.y != p.y
                && (px < (p_next.x - p.x) * (py - p.y) / (p_next.y - p.y) + p.x);
            p = p_next;
            if std::ptr::eq(p, a) {
                break;
            }
        }
        return inside;
    }

    fn split_polygon(&mut self, a_i: usize, b_i: usize) -> usize {
        let a = &self.nodes[a_i];
        let a2_i = self.nodes.len();
        self.nodes.push(Node::new(a.i, a.x, a.y));

        let b = &self.nodes[b_i];
        let b2_i = self.nodes.len();
        self.nodes.push(Node::new(b.i, b.x, b.y));

        let an_i = self.nodes[a_i].next_i.unwrap();
        let bp_i = self.nodes[b_i].prev_i.unwrap();

        self.nodes[a_i].next_i = Some(b_i);
        self.nodes[b_i].prev_i = Some(a_i);
        self.nodes[a2_i].next_i = Some(an_i);
        self.nodes[an_i].prev_i = Some(a2_i);
        self.nodes[b2_i].next_i = Some(a2_i);
        self.nodes[a2_i].prev_i = Some(b2_i);
        self.nodes[bp_i].next_i = Some(b2_i);
        self.nodes[b2_i].prev_i = Some(bp_i);

        return b2_i;
    }

    fn insert_node(&mut self, i: usize, x: T, y: T, last: Option<usize>) -> usize {
        let p = Node::new(i, x, y);
        let p_i = self.nodes.len();
        self.nodes.push(p);
        match last {
            Some(last_i) => {
                let last_next_i = self.nodes[last_i].next_i.unwrap();
                self.nodes[p_i].next_i = Some(last_next_i);
                self.nodes[p_i].prev_i = Some(last_i);
                self.nodes[last_next_i].prev_i = Some(p_i);
                self.nodes[last_i].next_i = Some(p_i);
            }
            None => {
                self.nodes[p_i].prev_i = Some(p_i);
                self.nodes[p_i].next_i = Some(p_i);
            }
        }
        return p_i;
    }

    #[inline(always)]
    fn remove_node(&mut self, p_i: usize) {
        let p = &self.nodes[p_i];
        let p_next_i = p.next_i.unwrap();
        let p_prev_i = p.prev_i.unwrap();
        let p_next_z_i = p.next_z_i;
        let p_prev_z_i = p.prev_z_i;

        self.nodes[p_next_i].prev_i = Some(p_prev_i);
        self.nodes[p_prev_i].next_i = Some(p_next_i);
        if let Some(prev_z_i) = p_prev_z_i {
            self.nodes[prev_z_i].next_z_i = p_next_z_i;
        }
        if let Some(next_z_i) = p_next_z_i {
            self.nodes[next_z_i].prev_z_i = p_prev_z_i;
        }
    }
}

pub fn deviation<T: Float, N: Index>(
    data: &[T],
    hole_indices: &[N],
    dim: usize,
    triangles: &[N],
) -> T {
    let has_holes = hole_indices.len() > 0;
    let outer_len = match has_holes {
        true => hole_indices[0].into() * dim,
        false => data.len(),
    };
    let mut polygon_area = signed_area(data, 0, outer_len, dim).abs();
    if has_holes {
        for i in 0..hole_indices.len() {
            let start = hole_indices[i].into() * dim;
            let end = if i < hole_indices.len() - 1 {
                hole_indices[i + 1].into() * dim
            } else {
                data.len()
            };
            polygon_area = polygon_area - signed_area(data, start, end, dim).abs();
        }
    }

    let mut triangles_area = T::zero();
    for i in (0..triangles.len()).step_by(3) {
        let a = triangles[i].into() * dim;
        let b = triangles[i + 1].into() * dim;
        let c = triangles[i + 2].into() * dim;
        triangles_area = triangles_area
            + ((data[a] - data[c]) * (data[b + 1] - data[a + 1])
                - (data[a] - data[b]) * (data[c + 1] - data[a + 1]))
                .abs();
    }

    return if polygon_area == T::zero() && triangles_area == T::zero() {
        T::zero()
    } else {
        ((polygon_area - triangles_area) / polygon_area).abs()
    };
}

fn z_order<T: Float>(x: T, y: T, min_x: T, min_y: T, inv_size: T) -> u32 {
    // coords are transformed into non-negative 15-bit integer range
    let mut x = ((x - min_x) * inv_size).to_u32().unwrap();
    let mut y = ((y - min_y) * inv_size).to_u32().unwrap();
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;
    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;
    return x | (y << 1);
}

fn point_in_triangle<T: Float>(ax: T, ay: T, bx: T, by: T, cx: T, cy: T, px: T, py: T) -> bool {
    return (cx - px) * (ay - py) >= (ax - px) * (cy - py)
        && (ax - px) * (by - py) >= (bx - px) * (ay - py)
        && (bx - px) * (cy - py) >= (cx - px) * (by - py);
}

#[inline(always)]
fn signed_area<T: Float>(data: &[T], start: usize, end: usize, dim: usize) -> T {
    let mut sum = T::zero();
    let mut j = if end > dim { end - dim } else { 0 };
    for i in (start..end).step_by(dim) {
        sum = sum + (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
        j = i;
    }
    return sum;
}

#[inline(always)]
fn area<T: Float>(p: &Node<T>, q: &Node<T>, r: &Node<T>) -> T {
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
}

#[inline(always)]
fn equals<T: Float>(p1: &Node<T>, p2: &Node<T>) -> bool {
    return p1.x == p2.x && p1.y == p2.y;
}

#[inline(always)]
fn on_segment<T: Float>(p: &Node<T>, q: &Node<T>, r: &Node<T>) -> bool {
    return q.x <= p.x.max(r.x)
        && q.x >= p.x.min(r.x)
        && q.y <= p.y.max(r.y)
        && q.y >= p.y.max(r.y);
}
