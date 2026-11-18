//! Trees: neighbour-joining, midpoint rooting, Newick I/O, branch collapsing and
//! Robinson-Foulds distance (6.4, 6.9, 11).
//!
//! All of this is implemented in the project rather than taken from a crate, as
//! section 12 suggests: each piece is small, and owning them keeps the
//! Windows build free of C dependencies.
//!
//! Every routine here is deterministic. Ties are broken by name or by index,
//! never by iteration order of a hash map, because section 0 rule 3 requires
//! byte-identical output for any thread count.

use std::collections::HashMap;
use std::fmt::Write as _;

use serde::{Deserialize, Serialize};

/// A rooted phylogenetic tree.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Tree {
    Leaf { name: String, len: f64 },
    Internal { len: f64, children: Vec<Tree> },
}

impl Default for Tree {
    /// The empty tree. Used where a `Tree` field must exist before any
    /// clustering has happened; it has no leaves and serialises as `();`.
    fn default() -> Self {
        Tree::Internal {
            len: 0.0,
            children: Vec::new(),
        }
    }
}

impl Tree {
    pub fn leaf(name: impl Into<String>) -> Self {
        Tree::Leaf {
            name: name.into(),
            len: 0.0,
        }
    }

    pub fn len(&self) -> f64 {
        match self {
            Tree::Leaf { len, .. } | Tree::Internal { len, .. } => *len,
        }
    }

    pub fn set_len(&mut self, v: f64) {
        match self {
            Tree::Leaf { len, .. } | Tree::Internal { len, .. } => *len = v,
        }
    }

    pub fn is_leaf(&self) -> bool {
        matches!(self, Tree::Leaf { .. })
    }

    pub fn children(&self) -> &[Tree] {
        match self {
            Tree::Leaf { .. } => &[],
            Tree::Internal { children, .. } => children,
        }
    }

    pub fn children_mut(&mut self) -> &mut [Tree] {
        match self {
            Tree::Leaf { .. } => &mut [],
            Tree::Internal { children, .. } => children,
        }
    }

    /// Leaf names in a deterministic (document) order.
    pub fn leaf_names(&self) -> Vec<String> {
        let mut out = Vec::new();
        self.collect_names(&mut out);
        out
    }

    fn collect_names(&self, out: &mut Vec<String>) {
        match self {
            Tree::Leaf { name, .. } => out.push(name.clone()),
            Tree::Internal { children, .. } => {
                for c in children {
                    c.collect_names(out);
                }
            }
        }
    }

    pub fn leaf_count(&self) -> usize {
        match self {
            Tree::Leaf { .. } => 1,
            Tree::Internal { children, .. } => children.iter().map(|c| c.leaf_count()).sum(),
        }
    }

    /// Replace every negative branch length with zero, as 6.9 requires.
    pub fn clamp_negative_lengths(&mut self) {
        if self.len() < 0.0 || !self.len().is_finite() {
            self.set_len(0.0);
        }
        for c in self.children_mut() {
            c.clamp_negative_lengths();
        }
    }

    /// Total branch length.
    pub fn total_length(&self) -> f64 {
        let mut s = self.len();
        for c in self.children() {
            s += c.total_length();
        }
        s
    }
}

// ---------------------------------------------------------------------------
// Newick
// ---------------------------------------------------------------------------

/// Write a tree as a Newick string, terminated with `;`.
///
/// The root's own branch length is omitted when it is zero, which is the
/// convention every tree reader expects; a spurious `:0` on the root makes some
/// readers treat the tree as unrooted-with-a-stem.
pub fn to_newick(t: &Tree) -> String {
    let mut s = String::new();
    write_newick(t, &mut s, true);
    s.push(';');
    s
}

fn write_newick(t: &Tree, out: &mut String, is_root: bool) {
    let skip_len = is_root && t.len() == 0.0;
    match t {
        Tree::Leaf { name, len } => {
            out.push_str(&escape_name(name));
            if !skip_len {
                let _ = write!(out, ":{}", fmt_len(*len));
            }
        }
        Tree::Internal { len, children } => {
            out.push('(');
            for (i, c) in children.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_newick(c, out, false);
            }
            out.push(')');
            if !skip_len {
                let _ = write!(out, ":{}", fmt_len(*len));
            }
        }
    }
}

/// Format a branch length with enough digits to round-trip but without noise.
fn fmt_len(v: f64) -> String {
    if !v.is_finite() {
        return "0".to_string();
    }
    format!("{:.9}", v)
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

/// Quote a name when it contains a character Newick treats as structural.
fn escape_name(name: &str) -> String {
    if name.is_empty() {
        return "''".to_string();
    }
    let needs_quote = name
        .chars()
        .any(|c| matches!(c, '(' | ')' | '[' | ']' | ',' | ':' | ';' | ' ' | '\''));
    if !needs_quote {
        return name.to_string();
    }
    format!("'{}'", name.replace('\'', "''"))
}

/// Read a Newick string. Square-bracket comments are skipped, as is the
/// optional trailing `;`.
pub fn parse_newick(input: &str) -> Result<Tree, String> {
    let cleaned = strip_comments(input);
    let bytes: Vec<char> = cleaned.chars().collect();
    let mut pos = 0usize;
    skip_ws(&bytes, &mut pos);
    if pos >= bytes.len() {
        return Err("empty Newick string".to_string());
    }
    let t = parse_subtree(&bytes, &mut pos)?;
    skip_ws(&bytes, &mut pos);
    if pos < bytes.len() && bytes[pos] == ';' {
        pos += 1;
    }
    skip_ws(&bytes, &mut pos);
    if pos != bytes.len() {
        return Err(format!(
            "trailing characters at offset {pos} in Newick string"
        ));
    }
    Ok(t)
}

fn strip_comments(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth = 0i32;
    for c in s.chars() {
        match c {
            '[' => depth += 1,
            ']' => depth = (depth - 1).max(0),
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out
}

fn skip_ws(bytes: &[char], pos: &mut usize) {
    while *pos < bytes.len() && bytes[*pos].is_whitespace() {
        *pos += 1;
    }
}

fn parse_subtree(bytes: &[char], pos: &mut usize) -> Result<Tree, String> {
    skip_ws(bytes, pos);
    let node = if *pos < bytes.len() && bytes[*pos] == '(' {
        *pos += 1;
        let mut children = Vec::new();
        loop {
            children.push(parse_subtree(bytes, pos)?);
            skip_ws(bytes, pos);
            if *pos < bytes.len() && bytes[*pos] == ',' {
                *pos += 1;
                continue;
            }
            break;
        }
        skip_ws(bytes, pos);
        if *pos >= bytes.len() || bytes[*pos] != ')' {
            return Err("expected ')' in Newick string".to_string());
        }
        *pos += 1;
        // Internal node labels are parsed (so the offset advances past them)
        // but not kept: `Tree` has no field for them, and nothing in the
        // outline uses them. Section 6.9 identifies cluster subtrees by their
        // pivot leaf name, which is a leaf label and is kept.
        parse_name(bytes, pos);
        let len = parse_len(bytes, pos)?;
        Tree::Internal { len, children }
    } else {
        let name = parse_name(bytes, pos);
        let len = parse_len(bytes, pos)?;
        Tree::Leaf { name, len }
    };
    Ok(node)
}

fn parse_name(bytes: &[char], pos: &mut usize) -> String {
    skip_ws(bytes, pos);
    if *pos < bytes.len() && bytes[*pos] == '\'' {
        *pos += 1;
        let mut s = String::new();
        while *pos < bytes.len() {
            let c = bytes[*pos];
            *pos += 1;
            if c == '\'' {
                if *pos < bytes.len() && bytes[*pos] == '\'' {
                    s.push('\'');
                    *pos += 1;
                } else {
                    break;
                }
            } else {
                s.push(c);
            }
        }
        return s;
    }
    let mut s = String::new();
    while *pos < bytes.len() {
        let c = bytes[*pos];
        if matches!(c, '(' | ')' | ',' | ':' | ';' | '[' | ']') || c.is_whitespace() {
            break;
        }
        s.push(c);
        *pos += 1;
    }
    s
}

fn parse_len(bytes: &[char], pos: &mut usize) -> Result<f64, String> {
    skip_ws(bytes, pos);
    if *pos >= bytes.len() || bytes[*pos] != ':' {
        return Ok(0.0);
    }
    *pos += 1;
    skip_ws(bytes, pos);
    let start = *pos;
    while *pos < bytes.len() {
        let c = bytes[*pos];
        if c.is_ascii_digit() || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E' {
            *pos += 1;
        } else {
            break;
        }
    }
    let txt: String = bytes[start..*pos].iter().collect();
    if txt.is_empty() {
        return Ok(0.0);
    }
    txt.parse::<f64>()
        .map_err(|e| format!("bad branch length {txt:?}: {e}"))
}

// ---------------------------------------------------------------------------
// Neighbour joining
// ---------------------------------------------------------------------------

/// Neighbour-joining on a symmetric distance matrix.
///
/// `names` are the leaf labels and `dist[i][j]` the distance between `i` and
/// `j`. Returns an unrooted tree represented as a root with three children
/// (two when there are only two taxa). Negative branch lengths are left as
/// produced; call [`Tree::clamp_negative_lengths`] when the outline says to
/// zero them (6.9).
pub fn neighbor_joining(names: &[String], dist: &[Vec<f64>]) -> Tree {
    let n = names.len();
    match n {
        0 => Tree::Internal {
            len: 0.0,
            children: Vec::new(),
        },
        1 => Tree::Leaf {
            name: names[0].clone(),
            len: 0.0,
        },
        2 => Tree::Internal {
            len: 0.0,
            children: vec![
                Tree::Leaf {
                    name: names[0].clone(),
                    len: 0.0,
                },
                Tree::Leaf {
                    name: names[1].clone(),
                    len: dist[0][1].max(0.0),
                },
            ],
        },
        _ => nj_general(names, dist),
    }
}

fn nj_general(names: &[String], dist: &[Vec<f64>]) -> Tree {
    let n = names.len();
    // Working matrix on a fixed-size index space; removed taxa are marked.
    let mut d: Vec<Vec<f64>> = dist.to_vec();
    let mut active: Vec<usize> = (0..n).collect();
    // node[i] holds the subtree currently attached to index i.
    let mut node: Vec<Option<Tree>> = (0..n)
        .map(|i| {
            Some(Tree::Leaf {
                name: names[i].clone(),
                len: 0.0,
            })
        })
        .collect();
    let mut next_id = n;

    while active.len() > 2 {
        let m = active.len() as f64;
        // Row sums over active taxa.
        let mut sum = vec![0f64; next_id];
        for &i in &active {
            let mut s = 0.0;
            for &j in &active {
                if i != j {
                    s += get(&d, i, j);
                }
            }
            sum[i] = s;
        }
        // Find the pair minimising D(i,j) - u_i - u_j, ties by index.
        let mut best: Option<(usize, usize, f64)> = None;
        for ai in 0..active.len() {
            for aj in (ai + 1)..active.len() {
                let i = active[ai];
                let j = active[aj];
                let q = get(&d, i, j) - sum[i] / (m - 2.0) - sum[j] / (m - 2.0);
                if best.map(|(_, _, bq)| q < bq - 1e-12).unwrap_or(true) {
                    best = Some((i, j, q));
                }
            }
        }
        let (i, j, _) = best.expect("at least three active taxa");
        let dij = get(&d, i, j);
        let ui = sum[i] / (m - 2.0);
        let uj = sum[j] / (m - 2.0);
        let li = dij / 2.0 + (ui - uj) / 2.0;
        let lj = dij - li;

        let mut ci = node[i].take().unwrap_or_else(|| Tree::leaf(""));
        let mut cj = node[j].take().unwrap_or_else(|| Tree::leaf(""));
        ci.set_len(li);
        cj.set_len(lj);

        // Distances from the new node k to every remaining active taxon.
        let k = next_id;
        next_id += 1;
        grow(&mut d, next_id);
        for &l in &active {
            if l == i || l == j {
                continue;
            }
            let v = (get(&d, i, l) + get(&d, j, l) - dij) / 2.0;
            set(&mut d, k, l, v);
        }
        set(&mut d, k, i, 0.0);
        set(&mut d, k, j, 0.0);
        node.push(Some(Tree::Internal {
            len: 0.0,
            children: vec![ci, cj],
        }));
        active.retain(|&x| x != i && x != j);
        active.push(k);
    }

    // Two taxa left: join them, splitting the remaining distance.
    let (i, j) = (active[0], active[1]);
    let dij = get(&d, i, j);
    let mut ci = node[i].take().unwrap_or_else(|| Tree::leaf(""));
    let mut cj = node[j].take().unwrap_or_else(|| Tree::leaf(""));
    // Deterministic split: give the alphabetically-first leaf the shorter side
    // only when the two are equal, otherwise split evenly. NJ leaves the final
    // edge unidentifiable, so an even split is the convention.
    let half = dij / 2.0;
    ci.set_len(half);
    cj.set_len(dij - half);
    Tree::Internal {
        len: 0.0,
        children: vec![ci, cj],
    }
}

fn get(d: &[Vec<f64>], i: usize, j: usize) -> f64 {
    d.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0)
}

fn set(d: &mut [Vec<f64>], i: usize, j: usize, v: f64) {
    if i < d.len() && j < d[i].len() {
        d[i][j] = v;
    }
    if j < d.len() && i < d[j].len() {
        d[j][i] = v;
    }
}

fn grow(d: &mut Vec<Vec<f64>>, size: usize) {
    while d.len() < size {
        d.push(vec![0.0; d.len()]);
    }
    for row in d.iter_mut() {
        while row.len() < size {
            row.push(0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Midpoint rooting
// ---------------------------------------------------------------------------

/// Root a tree at the midpoint of its longest leaf-to-leaf path (6.4, 6.9).
///
/// The tree is first read as an undirected weighted graph, the diameter is
/// found with two depth-first passes, and the edge carrying the midpoint is
/// split in two. A tree that is already rooted stays a valid tree; this only
/// moves where the root sits.
pub fn midpoint_root(t: Tree) -> Tree {
    let leaves = t.leaf_count();
    if leaves <= 1 {
        return t;
    }
    let Graph { names, adj } = to_graph(&t);
    let leaf_idx: Vec<usize> = (0..names.len()).filter(|&i| adj[i].len() == 1).collect();
    if leaf_idx.len() < 2 {
        return t;
    }

    // Two-pass diameter: farthest from any leaf, then farthest from that.
    let (a, _) = farthest(&adj, leaf_idx[0]);
    let (b, total) = farthest(&adj, a);
    if total <= 0.0 {
        return t;
    }
    let path = path_between(&adj, a, b);
    if path.len() < 2 {
        return t;
    }

    // Walk the path to the midpoint.
    let half = total / 2.0;
    let mut acc = 0.0;
    let mut split_edge: Option<(usize, usize, f64, f64)> = None; // (u, v, len, x)
    for w in path.windows(2) {
        let (u, v) = (w[0], w[1]);
        let l = edge_len(&adj, u, v);
        if acc + l >= half - 1e-12 {
            let x = (half - acc).clamp(0.0, l);
            split_edge = Some((u, v, l, x));
            break;
        }
        acc += l;
    }
    let (u, v, l, x) = match split_edge {
        Some(e) => e,
        None => return t,
    };

    // Insert a root node inside the edge (u, v).
    let root_id = names.len();
    let mut adj2 = adj;
    let mut names2 = names;
    names2.push(None);
    // Remove u-v from both adjacency lists.
    remove_edge(&mut adj2, u, v);
    adj2.push(vec![(u, x), (v, l - x)]);
    adj2[u].push((root_id, x));
    adj2[v].push((root_id, l - x));

    from_graph(&names2, &adj2, root_id, usize::MAX, 0.0)
}

fn edge_len(adj: &[Vec<(usize, f64)>], u: usize, v: usize) -> f64 {
    adj[u]
        .iter()
        .find(|&&(t, _)| t == v)
        .map(|&(_, l)| l)
        .unwrap_or(0.0)
}

fn remove_edge(adj: &mut [Vec<(usize, f64)>], u: usize, v: usize) {
    if let Some(row) = adj.get_mut(u) {
        row.retain(|&(t, _)| t != v);
    }
    if let Some(row) = adj.get_mut(v) {
        row.retain(|&(t, _)| t != u);
    }
}

/// A tree flattened into an undirected graph: node names, plus adjacency lists
/// of `(neighbour, edge length)`. Only leaves are named; internal nodes carry
/// `None`.
struct Graph {
    names: Vec<Option<String>>,
    adj: Vec<Vec<(usize, f64)>>,
}

/// Flatten a tree into an undirected graph (see [`Graph`]).
fn to_graph(t: &Tree) -> Graph {
    let mut g = Graph {
        names: Vec::new(),
        adj: Vec::new(),
    };
    build_graph(t, &mut g.names, &mut g.adj);
    g
}

fn build_graph(
    t: &Tree,
    names: &mut Vec<Option<String>>,
    adj: &mut Vec<Vec<(usize, f64)>>,
) -> usize {
    let id = names.len();
    names.push(match t {
        Tree::Leaf { name, .. } => Some(name.clone()),
        Tree::Internal { .. } => None,
    });
    adj.push(Vec::new());
    for c in t.children() {
        let child = build_graph(c, names, adj);
        let l = c.len();
        adj[id].push((child, l));
        adj[child].push((id, l));
    }
    id
}

/// Depth-first search back into a rooted tree.
fn from_graph(
    names: &[Option<String>],
    adj: &[Vec<(usize, f64)>],
    node: usize,
    parent: usize,
    len: f64,
) -> Tree {
    let mut children = Vec::new();
    for &(t, l) in &adj[node] {
        if t == parent {
            continue;
        }
        children.push(from_graph(names, adj, t, node, l));
    }
    if children.is_empty() {
        Tree::Leaf {
            name: names[node].clone().unwrap_or_default(),
            len,
        }
    } else {
        Tree::Internal { len, children }
    }
}

/// Farthest node from `start`, with its distance. Ties go to the lower index.
fn farthest(adj: &[Vec<(usize, f64)>], start: usize) -> (usize, f64) {
    let n = adj.len();
    let mut dist = vec![f64::NAN; n];
    let mut stack = vec![(start, usize::MAX, 0.0)];
    let mut best = (start, 0.0);
    while let Some((u, parent, d)) = stack.pop() {
        if !dist[u].is_nan() {
            continue;
        }
        dist[u] = d;
        if d > best.1 || (d == best.1 && u < best.0) {
            best = (u, d);
        }
        for &(v, l) in &adj[u] {
            if v != parent && dist[v].is_nan() {
                stack.push((v, u, d + l));
            }
        }
    }
    best
}

/// The node sequence from `a` to `b`, inclusive. Empty when unreachable.
fn path_between(adj: &[Vec<(usize, f64)>], a: usize, b: usize) -> Vec<usize> {
    let n = adj.len();
    let mut prev = vec![usize::MAX; n];
    let mut seen = vec![false; n];
    let mut stack = vec![a];
    seen[a] = true;
    while let Some(u) = stack.pop() {
        if u == b {
            break;
        }
        for &(v, _) in &adj[u] {
            if !seen[v] {
                seen[v] = true;
                prev[v] = u;
                stack.push(v);
            }
        }
    }
    if !seen[b] {
        return Vec::new();
    }
    let mut path = vec![b];
    let mut cur = b;
    while cur != a {
        cur = prev[cur];
        if cur == usize::MAX {
            return Vec::new();
        }
        path.push(cur);
    }
    path.reverse();
    path
}

// ---------------------------------------------------------------------------
// Grafting the full tree together (6.9)
// ---------------------------------------------------------------------------

/// Replace the leaf called `name` with `subtree`. Returns false when there is
/// no such leaf. This is how the backbone tree (leaves are clusters) becomes
/// the full tree (6.9): each pivot leaf is swapped for its cluster tree.
pub fn graft(t: &mut Tree, name: &str, subtree: Tree) -> bool {
    match t {
        Tree::Leaf { name: n, .. } => {
            if n == name {
                *t = subtree;
                true
            } else {
                false
            }
        }
        Tree::Internal { children, .. } => {
            for (i, c) in children.iter_mut().enumerate() {
                if let Tree::Leaf { name: n, .. } = c {
                    if n == name {
                        let mut s = subtree;
                        // The backbone edge length stays: the cluster subtree
                        // hangs below it.
                        let l = c.len();
                        s.set_len(s.len() + l);
                        children[i] = s;
                        return true;
                    }
                }
                if graft(c, name, subtree.clone()) {
                    return true;
                }
            }
            false
        }
    }
}

// ---------------------------------------------------------------------------
// Branch collapsing and Robinson-Foulds
// ---------------------------------------------------------------------------

/// Contract every edge shorter than `min_len`, turning the nodes it leads to
/// into polytomies. Section 11 requires this before comparing trees, so that a
/// branch of less than half a SNP cannot count as a topological difference.
pub fn collapse_short(t: Tree, min_len: f64) -> Tree {
    match t {
        leaf @ Tree::Leaf { .. } => leaf,
        Tree::Internal { len, children } => {
            let mut out = Vec::new();
            for c in children {
                push_collapsed(c, min_len, &mut out);
            }
            Tree::Internal { len, children: out }
        }
    }
}

fn push_collapsed(t: Tree, min_len: f64, out: &mut Vec<Tree>) {
    match t {
        leaf @ Tree::Leaf { .. } => out.push(leaf),
        Tree::Internal { len, children } => {
            if len < min_len && len.is_finite() {
                for c in children {
                    push_collapsed(c, min_len, out);
                }
            } else {
                let mut inner = Vec::new();
                for c in children {
                    push_collapsed(c, min_len, &mut inner);
                }
                out.push(Tree::Internal {
                    len,
                    children: inner,
                });
            }
        }
    }
}

/// The set of bipartitions (splits) a tree implies, as canonical sorted leaf
/// index sets. Each split is stored as the side whose smallest member is
/// smaller, so a split and its complement collapse to the same key. Trivial
/// splits (one leaf, or all but one) are excluded, which is what makes the
/// Robinson-Foulds count below meaningful.
pub fn splits(t: &Tree, leaf_index: &HashMap<String, usize>) -> Vec<Vec<usize>> {
    let total = leaf_index.len();
    let mut out = Vec::new();
    let mut scratch = vec![false; total];
    collect_splits(t, leaf_index, total, &mut out, &mut scratch);
    out.sort();
    out.dedup();
    out
}

fn collect_splits(
    t: &Tree,
    leaf_index: &HashMap<String, usize>,
    total: usize,
    out: &mut Vec<Vec<usize>>,
    scratch: &mut [bool],
) -> Vec<usize> {
    match t {
        Tree::Leaf { name, .. } => match leaf_index.get(name) {
            Some(i) => vec![*i],
            None => Vec::new(),
        },
        Tree::Internal { children, .. } => {
            let mut here = Vec::new();
            for c in children {
                let sub = collect_splits(c, leaf_index, total, out, scratch);
                // A non-trivial bipartition: at least two leaves on this side,
                // and at least two on the other.
                if sub.len() > 1 && sub.len() + 1 < total {
                    for &i in &sub {
                        scratch[i] = true;
                    }
                    let comp: Vec<usize> = (0..total).filter(|i| !scratch[*i]).collect();
                    for &i in &sub {
                        scratch[i] = false;
                    }
                    // Canonical side: the one whose smallest member is smaller.
                    let key = if sub[0] < comp[0] { sub.clone() } else { comp };
                    out.push(key);
                }
                here.extend(sub);
            }
            here.sort_unstable();
            here.dedup();
            here
        }
    }
}

/// Robinson-Foulds symmetric difference between two trees over the same leaf
/// set. Returns `(rf, splits_a, splits_b)`; the normalised distance is
/// `rf / (splits_a + splits_b)`.
pub fn rf_distance(a: &Tree, b: &Tree) -> (usize, usize, usize) {
    let mut leaf_index: HashMap<String, usize> = HashMap::new();
    for n in a.leaf_names() {
        if !leaf_index.contains_key(&n) {
            let i = leaf_index.len();
            leaf_index.insert(n, i);
        }
    }
    for n in b.leaf_names() {
        if !leaf_index.contains_key(&n) {
            let i = leaf_index.len();
            leaf_index.insert(n, i);
        }
    }
    let sa = splits(a, &leaf_index);
    let sb = splits(b, &leaf_index);
    let mut only_a = 0;
    let mut only_b = 0;
    let (mut i, mut j) = (0, 0);
    while i < sa.len() && j < sb.len() {
        if sa[i] == sb[j] {
            i += 1;
            j += 1;
        } else if sa[i] < sb[j] {
            only_a += 1;
            i += 1;
        } else {
            only_b += 1;
            j += 1;
        }
    }
    only_a += sa.len() - i;
    only_b += sb.len() - j;
    (only_a + only_b, sa.len(), sb.len())
}

/// True when `set`, or its complement, is one of the splits the tree implies.
///
/// Callers should not have to know which side of a bipartition [`splits`]
/// canonicalises to; this takes a leaf set in any orientation.
pub fn has_split(t: &Tree, leaf_index: &HashMap<String, usize>, set: &[usize]) -> bool {
    let total = leaf_index.len();
    let mut s = set.to_vec();
    s.sort_unstable();
    s.dedup();
    if s.is_empty() || s.len() >= total || total < 4 {
        return false;
    }
    let mut marked = vec![false; total];
    for &i in &s {
        if i < total {
            marked[i] = true;
        }
    }
    let comp: Vec<usize> = (0..total).filter(|i| !marked[*i]).collect();
    let key = if s[0] < comp[0] { s } else { comp };
    splits(t, leaf_index).contains(&key)
}

/// Normalised Robinson-Foulds distance in `[0, 1]`. `None` when neither tree
/// has any split, so a ratio is never invented from nothing (section 0 rule 8).
pub fn rf_normalised(a: &Tree, b: &Tree) -> Option<f64> {
    let (rf, na, nb) = rf_distance(a, b);
    if na + nb == 0 {
        return None;
    }
    Some(rf as f64 / (na + nb) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names() -> Vec<String> {
        ["a", "b", "c", "d", "e"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    }

    #[test]
    fn newick_round_trips() {
        let t = Tree::Internal {
            len: 0.0,
            children: vec![
                Tree::Leaf {
                    name: "a".into(),
                    len: 1.5,
                },
                Tree::Internal {
                    len: 0.25,
                    children: vec![
                        Tree::Leaf {
                            name: "b".into(),
                            len: 2.0,
                        },
                        Tree::Leaf {
                            name: "c".into(),
                            len: 3.0,
                        },
                    ],
                },
            ],
        };
        let s = to_newick(&t);
        let back = parse_newick(&s).unwrap();
        assert_eq!(back, t);
        assert_eq!(s, "(a:1.5,(b:2,c:3):0.25);");
    }

    #[test]
    fn newick_names_with_odd_characters_are_quoted() {
        let t = Tree::Leaf {
            name: "strain 01".into(),
            len: 1.0,
        };
        let s = to_newick(&t);
        assert!(s.contains("'strain 01'"));
        assert_eq!(parse_newick(&s).unwrap(), t);
    }

    #[test]
    fn newick_parser_tolerates_comments_and_no_lengths() {
        let t = parse_newick("(a,b,(c,d))[&note] ;").unwrap();
        assert_eq!(t.leaf_count(), 4);
        assert_eq!(t.len(), 0.0);
    }

    #[test]
    fn nj_recovers_a_known_additive_tree() {
        // Distances from the tree ((a:1,b:1):1,(c:1,d:1):1) plus e hanging off
        // the root: NJ must group a with b and c with d.
        let n = names();
        //        a    b    c    d    e
        let d = vec![
            vec![0.0, 2.0, 4.0, 4.0, 3.0], // a
            vec![2.0, 0.0, 4.0, 4.0, 3.0], // b
            vec![4.0, 4.0, 0.0, 2.0, 3.0], // c
            vec![4.0, 4.0, 2.0, 0.0, 3.0], // d
            vec![3.0, 3.0, 3.0, 3.0, 0.0], // e
        ];
        let t = neighbor_joining(&n, &d);
        let s = to_newick(&t);
        assert!(s.contains("(a:") || s.contains("(b:"), "{s}");
        let idx = index_of(&t);
        assert!(has_split(&t, &idx, &[0, 1]), "a+b clade missing in {s}");
        assert!(has_split(&t, &idx, &[2, 3]), "c+d clade missing in {s}");
        assert!(
            !has_split(&t, &idx, &[0, 2]),
            "a+c must not be a clade in {s}"
        );
    }

    /// Leaf index by alphabetical name, so that assertions about which leaves
    /// form a clade do not depend on the order the tree happens to be written
    /// in. NJ output order is a property of the algorithm, not of the biology.
    fn index_of(t: &Tree) -> HashMap<String, usize> {
        let mut names = t.leaf_names();
        names.sort();
        let mut m = HashMap::new();
        for (i, n) in names.into_iter().enumerate() {
            m.insert(n, i);
        }
        m
    }

    #[test]
    fn nj_is_deterministic() {
        let n = names();
        let d = vec![
            vec![0.0, 5.0, 9.0, 9.0, 8.0],
            vec![5.0, 0.0, 10.0, 10.0, 9.0],
            vec![9.0, 10.0, 0.0, 8.0, 7.0],
            vec![9.0, 10.0, 8.0, 0.0, 3.0],
            vec![8.0, 9.0, 7.0, 3.0, 0.0],
        ];
        let a = to_newick(&neighbor_joining(&n, &d));
        let b = to_newick(&neighbor_joining(&n, &d));
        assert_eq!(a, b);
    }

    #[test]
    fn nj_handles_one_and_two_taxa() {
        let one = neighbor_joining(&["a".to_string()], &[vec![0.0]]);
        assert_eq!(one.leaf_count(), 1);
        let two = neighbor_joining(
            &["a".to_string(), "b".to_string()],
            &[vec![0.0, 4.0], vec![4.0, 0.0]],
        );
        assert_eq!(two.leaf_count(), 2);
        assert!((two.total_length() - 4.0).abs() < 1e-9);
    }

    #[test]
    fn midpoint_root_puts_the_root_in_the_middle_of_the_longest_path() {
        // A comb where a-b is 10 apart and everything else hangs close.
        let t = Tree::Internal {
            len: 0.0,
            children: vec![
                Tree::Leaf {
                    name: "a".into(),
                    len: 9.0,
                },
                Tree::Internal {
                    len: 1.0,
                    children: vec![
                        Tree::Leaf {
                            name: "b".into(),
                            len: 0.0,
                        },
                        Tree::Leaf {
                            name: "c".into(),
                            len: 0.0,
                        },
                    ],
                },
            ],
        };
        let r = midpoint_root(t);
        assert_eq!(r.leaf_count(), 3);
        // The root must have exactly two children after splitting an edge.
        assert_eq!(r.children().len(), 2);
        let names = r.leaf_names();
        assert!(names.contains(&"a".to_string()));
    }

    #[test]
    fn midpoint_root_preserves_pairwise_distances() {
        let n = names();
        let d = vec![
            vec![0.0, 2.0, 4.0, 4.0, 3.0],
            vec![2.0, 0.0, 4.0, 4.0, 3.0],
            vec![4.0, 4.0, 0.0, 2.0, 3.0],
            vec![4.0, 4.0, 2.0, 0.0, 3.0],
            vec![3.0, 3.0, 3.0, 3.0, 0.0],
        ];
        let t = neighbor_joining(&n, &d);
        let before = pairwise(&t, &n);
        let r = midpoint_root(t.clone());
        let after = pairwise(&r, &n);
        for i in 0..n.len() {
            for j in 0..n.len() {
                assert!(
                    (before[i][j] - after[i][j]).abs() < 1e-9,
                    "rooting changed d({},{}) from {} to {}",
                    n[i],
                    n[j],
                    before[i][j],
                    after[i][j]
                );
            }
        }
    }

    /// Leaf-to-leaf distances of a tree, for the tests above.
    fn pairwise(t: &Tree, names: &[String]) -> Vec<Vec<f64>> {
        let Graph { names: ns, adj } = to_graph(t);
        let mut idx = HashMap::new();
        for (i, n) in names.iter().enumerate() {
            // find node with this name
            for (j, nm) in ns.iter().enumerate() {
                if nm.as_deref() == Some(n.as_str()) {
                    idx.insert(i, j);
                }
            }
        }
        let n = names.len();
        let mut out = vec![vec![0f64; n]; n];
        for i in 0..n {
            let start = idx[&i];
            let mut dist = vec![f64::NAN; adj.len()];
            let mut stack = vec![(start, usize::MAX, 0.0)];
            while let Some((u, p, d)) = stack.pop() {
                if !dist[u].is_nan() {
                    continue;
                }
                dist[u] = d;
                for &(v, l) in &adj[u] {
                    if v != p && dist[v].is_nan() {
                        stack.push((v, u, d + l));
                    }
                }
            }
            for j in 0..n {
                out[i][j] = dist[idx[&j]];
            }
        }
        out
    }

    #[test]
    fn rf_is_zero_for_the_same_topology_and_positive_for_a_different_one() {
        let t1 = parse_newick("((a,b),(c,d),e);").unwrap();
        let t2 = parse_newick("((b,a),(d,c),e);").unwrap();
        assert_eq!(rf_distance(&t1, &t2).0, 0);
        let t3 = parse_newick("((a,c),(b,d),e);").unwrap();
        assert!(rf_distance(&t1, &t3).0 > 0);
        assert_eq!(rf_normalised(&t1, &t1), Some(0.0));
    }

    #[test]
    fn rf_of_two_unresolved_stars_is_undefined_not_zero() {
        let star = parse_newick("(a,b,c);").unwrap();
        // No non-trivial splits at all: the normalised distance must be None
        // rather than a fabricated 0.
        assert_eq!(rf_normalised(&star, &star), None);
    }

    #[test]
    fn collapsing_short_branches_removes_them() {
        // The (a,b) node hangs off a 0.2 edge, below the 0.5 threshold, so a
        // and b must become direct children of the root. The (c,d) node hangs
        // off a 4.0 edge, well above it, so it must survive.
        let t = parse_newick("((a:0.1,b:0.1):0.2,(c:5,d:5):4);").unwrap();
        let c = collapse_short(t, 0.5);
        let s = to_newick(&c);
        assert_eq!(s, "(a:0.1,b:0.1,(c:5,d:5):4);");
        assert_eq!(c.leaf_count(), 4);
    }

    #[test]
    fn collapsing_below_zero_changes_nothing() {
        let t = parse_newick("((a:0.1,b:0.1):0.2,c:1);").unwrap();
        assert_eq!(collapse_short(t.clone(), 0.0), t);
    }

    #[test]
    fn graft_replaces_a_leaf_with_a_subtree() {
        let mut backbone = parse_newick("(pivot1:1,pivot2:1);").unwrap();
        let sub = parse_newick("(m1:0.5,m2:0.5):0.5;").unwrap();
        assert!(graft(&mut backbone, "pivot1", sub));
        assert_eq!(backbone.leaf_count(), 3);
        assert!(!graft(&mut backbone, "nope", Tree::leaf("x")));
    }

    #[test]
    fn negative_lengths_are_clamped() {
        let mut t = parse_newick("(a:-1,b:2);").unwrap();
        t.clamp_negative_lengths();
        let s = to_newick(&t);
        assert!(s.contains("a:0"), "{s}");
        assert!(s.contains("b:2"), "{s}");
    }
}
