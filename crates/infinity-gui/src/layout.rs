//! Turning a tree into positions on screen (section 8, "Tree panel").
//!
//! egui has no tree widget, so the rectangular layout, the collapsed-cluster
//! triangles, the scale bar and the 300 ms move animation are all computed
//! here. Nothing in this module touches `egui::Ui` or a `Painter`: it maps data
//! to data, which is what makes the geometry testable without a window.
//!
//! # Units
//!
//! The layout space is dimensionless. Vertically one leaf row is one unit.
//! Every branch length arrives here already in SNPs (see
//! `UiModel::display_tree`): exact inside a Focus 2 cluster, and "≈" elsewhere
//! (the backbone and Focus 0 trees are sketch distances times the median
//! genome length, Focus 1 trees count difference pieces).
//!
//! Horizontally the tree is still normalised *per region*. Between clusters
//! the distances are thousands of SNPs, inside a cluster a handful; one shared
//! axis would squash every cluster into an invisible stub. The backbone
//! therefore spans a fixed layout width, and opened clusters span another. All
//! exact (Focus 2) clusters share *one* SNP scale, so a tight outbreak looks
//! tighter than a loose one; an approximate cluster gets a scale of its own,
//! because its estimates must not stretch or squash the exact ones. Each
//! [`Region`] carries its scale and its "≈" flag, and the scale bar always
//! names what it is drawing.

use std::collections::HashMap;

use infinity_core::tree::Tree;

use crate::fmt;

/// One leaf row, in layout units.
pub const ROW: f64 = 1.0;
/// Horizontal layout extent of the backbone region.
pub const BACKBONE_W: f64 = 1.0;
/// Horizontal layout extent of one expanded cluster subtree.
pub const CLUSTER_W: f64 = 0.7;
/// How far right a collapsed cluster's triangle reaches, in layout units.
pub const COLLAPSED_SPREAD: f64 = 0.05;
/// Section 8: "nodes move smoothly to their new place (about 300 ms)".
pub const ANIM_SECS: f64 = 0.3;

/// What a node of the display tree is.
#[derive(Debug, Clone, PartialEq)]
pub enum NodeKind {
    /// The invisible root of the whole picture.
    Root,
    /// A branching point with no meaning of its own.
    Internal,
    /// One genome.
    Genome { id: u32, name: String },
    /// A cluster drawn as a triangle with its member count (section 8).
    ClusterCollapsed { id: u32, members: usize, focus: u8 },
    /// The root of an opened cluster subtree.
    ClusterRoot { id: u32, members: usize, focus: u8 },
}

impl NodeKind {
    pub fn name(&self) -> Option<&str> {
        match self {
            NodeKind::Genome { name, .. } => Some(name),
            _ => None,
        }
    }

    /// True when this node occupies a leaf row of its own.
    pub fn is_terminal(&self) -> bool {
        matches!(
            self,
            NodeKind::Genome { .. } | NodeKind::ClusterCollapsed { .. }
        )
    }
}

/// A node of the tree the panel draws: an `infinity_core::tree::Tree` enriched
/// with the cluster and focus information the drawing needs.
#[derive(Debug, Clone, PartialEq)]
pub struct DisplayNode {
    pub kind: NodeKind,
    /// Length of the edge from this node's parent, in the region's own unit.
    pub len: f64,
    pub focus: u8,
    /// True when the number behind `len` is an estimate (section 0 rule 8).
    pub approximate: bool,
    pub children: Vec<DisplayNode>,
}

impl DisplayNode {
    pub fn new(kind: NodeKind, len: f64, focus: u8, approximate: bool) -> Self {
        Self {
            kind,
            len,
            focus,
            approximate,
            children: Vec::new(),
        }
    }
}

impl DisplayNode {
    /// Multiply every branch length in this subtree by `factor`, which is how a
    /// tree in substitutions per site is put into (approximate) SNPs.
    pub fn scale_lengths(&mut self, factor: f64) {
        let f = if factor.is_finite() && factor > 0.0 {
            factor
        } else {
            0.0
        };
        let mut stack: Vec<&mut DisplayNode> = vec![self];
        while let Some(n) = stack.pop() {
            n.len *= f;
            for c in &mut n.children {
                stack.push(c);
            }
        }
    }
}

/// Convert an engine tree into display nodes.
///
/// `ids` maps a leaf name to a genome id; a name that is not in it (which
/// happens when `unique_names` had to de-duplicate two identical file names)
/// gets `u32::MAX`, meaning "no genome behind this label".
pub fn from_core(
    t: &Tree,
    focus: u8,
    approximate: bool,
    ids: &HashMap<String, u32>,
) -> DisplayNode {
    match t {
        Tree::Leaf { name, len } => {
            let id = ids.get(name).copied().unwrap_or(u32::MAX);
            DisplayNode::new(
                NodeKind::Genome {
                    id,
                    name: name.clone(),
                },
                *len,
                focus,
                approximate,
            )
        }
        Tree::Internal { len, children } => {
            let mut node = DisplayNode::new(NodeKind::Internal, *len, focus, approximate);
            node.children = children
                .iter()
                .map(|c| from_core(c, focus, approximate, ids))
                .collect();
            node
        }
    }
}

/// A horizontal stretch of the picture that shares one unit and one scale.
#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    /// Layout x of the region's root.
    pub x0: f64,
    /// Layout width the region is normalised into.
    pub width: f64,
    /// Real units spanned by the region's deepest leaf.
    pub raw_max: f64,
    /// Name of the unit `raw_max` is measured in.
    pub unit: &'static str,
    pub cluster: Option<u32>,
    pub focus: u8,
    /// True when the lengths in this region are estimates ("≈").
    pub approximate: bool,
}

impl Region {
    /// Real units per layout unit, for the scale bar and the tooltips.
    pub fn units_per_layout(&self) -> f64 {
        if self.width > 0.0 {
            self.raw_max / self.width
        } else {
            0.0
        }
    }

    pub fn contains_x(&self, x: f64) -> bool {
        x >= self.x0 - 1e-9 && x <= self.x0 + self.width + 1e-9
    }
}

/// One laid-out node.
#[derive(Debug, Clone, PartialEq)]
pub struct LayoutNode {
    pub idx: usize,
    /// `Some(p)` with `p < idx`: parents always come first, which the animation
    /// relies on.
    pub parent: Option<usize>,
    pub kind: NodeKind,
    pub x: f64,
    pub y: f64,
    /// Row span of the subtree, for drawing brackets and triangles.
    pub row_min: f64,
    pub row_max: f64,
    pub leaf_count: usize,
    /// Length of the edge from the parent, in the parent region's unit.
    pub len: f64,
    pub focus: u8,
    pub approximate: bool,
    /// Index into [`TreeLayout::regions`].
    pub region: usize,
    /// Identity for the move animation: stable while the leaf set below the
    /// node is stable, and different as soon as it changes.
    pub key: u64,
    /// The opened cluster this node belongs to: the cluster root itself and
    /// everything below it. `None` on the backbone.
    pub cluster: Option<u32>,
    /// Edges between this node and the root of its region (the cluster root
    /// inside a cluster, the tree root on the backbone). The unfold animation
    /// staggers branches by it.
    pub depth: u32,
    /// Focus and "≈" of the edge *into* this node. The same as the node's own
    /// everywhere except at a cluster root, whose edge is a backbone branch: a
    /// sketch estimate even when the cluster below it is exact, and drawn so.
    pub edge_focus: u8,
    pub edge_approximate: bool,
}

/// A whole laid-out tree.
#[derive(Debug, Clone, Default)]
pub struct TreeLayout {
    pub nodes: Vec<LayoutNode>,
    pub root: usize,
    pub regions: Vec<Region>,
    /// Number of leaf rows.
    pub rows: f64,
    pub x_max: f64,
    /// Indices of the terminal nodes, in row order.
    pub leaves: Vec<usize>,
    /// Leaf name to node index, for the search box.
    pub by_name: HashMap<String, usize>,
}

impl TreeLayout {
    /// The region a layout x falls in, preferring a cluster over the backbone
    /// (a cluster region starts at a backbone x, so the two overlap there).
    pub fn region_at(&self, x: f64) -> Option<&Region> {
        self.regions
            .iter()
            .find(|r| r.cluster.is_some() && r.contains_x(x))
            .or_else(|| self.regions.first())
    }

    /// The node a leaf name resolves to, for "search box to find a genome".
    pub fn find(&self, name: &str) -> Option<usize> {
        self.by_name.get(name).copied()
    }
}

/// Lay a display tree out into [`TreeLayout`] coordinates.
pub fn layout(root: &DisplayNode) -> TreeLayout {
    let backbone_raw = region_raw(root);
    let bb_scale = scale_for(backbone_raw, BACKBONE_W);
    let mut out = TreeLayout {
        regions: vec![Region {
            x0: 0.0,
            width: BACKBONE_W,
            raw_max: backbone_raw,
            unit: fmt::unit_label(root.focus),
            cluster: None,
            focus: root.focus,
            approximate: root.approximate,
        }],
        ..Default::default()
    };
    let root_idx;
    let rows;
    {
        let mut b = Builder {
            rows: 0.0,
            exact_raw: shared_exact_raw(root),
            out: &mut out,
        };
        root_idx = b.build(root, None, 0.0, bb_scale, 0, None, 0);
        rows = b.rows;
    }
    out.root = root_idx;
    out.rows = rows.max(1.0);
    // A backbone with no length (one cluster drawn on its own) reserves no
    // width: the cluster then fills the panel like any standard phylogram.
    let skip_backbone = out.regions.len() > 1 && out.regions[0].raw_max <= 0.0;
    let x_max = out
        .nodes
        .iter()
        .map(|n| n.x)
        .chain(
            out.regions
                .iter()
                .enumerate()
                .filter(|(i, _)| !(skip_backbone && *i == 0))
                .map(|(_, r)| r.x0 + r.width),
        )
        .fold(0.0f64, f64::max);
    out.x_max = x_max;
    out
}

/// Layout scale for a region: how many layout units one real unit is worth.
fn scale_for(raw_max: f64, width: f64) -> f64 {
    if raw_max > 0.0 && raw_max.is_finite() {
        width / raw_max
    } else {
        0.0
    }
}

/// Deepest accumulated edge length inside one region. A `ClusterRoot` child
/// ends the region: its own edge belongs to the parent region, its contents to
/// a new one.
fn region_raw(node: &DisplayNode) -> f64 {
    let mut best = 0.0f64;
    for c in &node.children {
        let l = finite(c.len);
        let d = if matches!(c.kind, NodeKind::ClusterRoot { .. }) {
            l
        } else {
            l + region_raw(c)
        };
        best = best.max(d);
    }
    best
}

/// The widest exact cluster: every Focus 2 cluster is drawn on this one SNP
/// scale. `None` when no exact cluster is open.
fn shared_exact_raw(root: &DisplayNode) -> Option<f64> {
    let mut best: Option<f64> = None;
    let mut stack: Vec<&DisplayNode> = vec![root];
    while let Some(n) = stack.pop() {
        if let NodeKind::ClusterRoot { focus, .. } = n.kind {
            if focus == fmt::FOCUS_EXACT && !n.approximate {
                let raw = cluster_raw(n);
                if raw > 0.0 {
                    best = Some(best.map_or(raw, |b: f64| b.max(raw)));
                }
            }
            continue;
        }
        stack.extend(n.children.iter());
    }
    best
}

/// Deepest accumulated edge length below a cluster root, all of which lies in
/// the cluster's own region.
fn cluster_raw(node: &DisplayNode) -> f64 {
    let mut best = 0.0f64;
    for c in &node.children {
        best = best.max(finite(c.len) + cluster_raw(c));
    }
    best
}

/// Branch lengths are floats that come out of neighbour joining; a NaN must not
/// poison the whole layout.
fn finite(v: f64) -> f64 {
    if v.is_finite() {
        v.max(0.0)
    } else {
        0.0
    }
}

struct Builder<'o> {
    rows: f64,
    /// See [`shared_exact_raw`].
    exact_raw: Option<f64>,
    out: &'o mut TreeLayout,
}

impl<'o> Builder<'o> {
    // The recursion carries the region and cluster context down the tree; a
    // struct for it would only move the same eight values elsewhere.
    #[allow(clippy::too_many_arguments)]
    fn build(
        &mut self,
        node: &DisplayNode,
        parent: Option<usize>,
        x_parent: f64,
        scale: f64,
        region: usize,
        cluster: Option<u32>,
        depth: u32,
    ) -> usize {
        let len = finite(node.len);
        let x = x_parent + len * scale;

        // An opened cluster starts a new region with its own scale: the shared
        // exact scale when it is exact, one of its own otherwise.
        let (child_scale, child_region, cluster, depth) =
            if let NodeKind::ClusterRoot { id, focus, .. } = node.kind {
                let own = cluster_raw(node);
                let exact = focus == fmt::FOCUS_EXACT && !node.approximate;
                let raw = match self.exact_raw {
                    Some(shared) if exact => shared,
                    _ => own,
                };
                self.out.regions.push(Region {
                    x0: x,
                    width: CLUSTER_W,
                    raw_max: raw,
                    unit: fmt::unit_label(focus),
                    cluster: Some(id),
                    focus,
                    approximate: node.approximate,
                });
                (
                    scale_for(raw, CLUSTER_W),
                    self.out.regions.len() - 1,
                    Some(id),
                    0,
                )
            } else {
                (scale, region, cluster, depth)
            };

        let (edge_focus, edge_approximate) = match (&node.kind, parent) {
            (NodeKind::ClusterRoot { .. }, Some(p)) => {
                (self.out.nodes[p].focus, self.out.nodes[p].approximate)
            }
            _ => (node.focus, node.approximate),
        };
        let idx = self.out.nodes.len();
        self.out.nodes.push(LayoutNode {
            idx,
            parent,
            kind: node.kind.clone(),
            x,
            y: 0.0,
            row_min: 0.0,
            row_max: 0.0,
            leaf_count: 0,
            len,
            focus: node.focus,
            approximate: node.approximate,
            region: child_region,
            key: 0,
            cluster,
            depth,
            edge_focus,
            edge_approximate,
        });

        if node.children.is_empty() || node.kind.is_terminal() {
            let y = self.rows;
            self.rows += ROW;
            let key = node_key(&node.kind, &[]);
            let n = &mut self.out.nodes[idx];
            n.y = y;
            n.row_min = y;
            n.row_max = y;
            n.leaf_count = 1;
            n.key = key;
            if let Some(name) = node.kind.name() {
                self.out.by_name.entry(name.to_string()).or_insert(idx);
            }
            self.out.leaves.push(idx);
        } else {
            let mut kids: Vec<usize> = Vec::with_capacity(node.children.len());
            for c in &node.children {
                kids.push(self.build(
                    c,
                    Some(idx),
                    x,
                    child_scale,
                    child_region,
                    cluster,
                    depth + 1,
                ));
            }
            let row_min = self.out.nodes[kids[0]].row_min;
            let row_max = self.out.nodes[kids[kids.len() - 1]].row_max;
            let leaf_count: usize = kids.iter().map(|&k| self.out.nodes[k].leaf_count).sum();
            let keys: Vec<u64> = kids.iter().map(|&k| self.out.nodes[k].key).collect();
            let kind = self.out.nodes[idx].kind.clone();
            let key = node_key(&kind, &keys);
            let n = &mut self.out.nodes[idx];
            n.y = (row_min + row_max) / 2.0;
            n.row_min = row_min;
            n.row_max = row_max;
            n.leaf_count = leaf_count;
            n.key = key;
        }

        idx
    }
}

/// A 64-bit mixer, so the animation identity needs no hashing crate.
fn mix(h: u64, v: u64) -> u64 {
    let mut x = h ^ v.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^= x >> 33;
    x
}

fn hash_str(s: &str) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for b in s.as_bytes() {
        h = mix(h, u64::from(*b));
    }
    h
}

/// Animation identity: a leaf is identified by its label, an internal node by
/// the ordered identities of its children. A node therefore keeps its identity
/// while the leaves below it are unchanged, which is exactly when a smooth move
/// rather than a redraw is the honest response to a new tree.
pub fn node_key(kind: &NodeKind, child_keys: &[u64]) -> u64 {
    match kind {
        NodeKind::Genome { name, .. } => mix(1, hash_str(name)),
        NodeKind::ClusterCollapsed { id, .. } => mix(2, u64::from(*id)),
        NodeKind::ClusterRoot { id, .. } => mix(3, u64::from(*id)),
        NodeKind::Root => 4,
        NodeKind::Internal => {
            let mut h = 5u64;
            for k in child_keys {
                h = mix(h, *k);
            }
            h
        }
    }
}

/// Half-height, in rows, of the triangle a collapsed cluster is drawn as. It
/// grows with the member count so that a cluster of 500 is visibly heavier than
/// a pair, but stays under half a row so that neighbours never overlap.
pub fn collapsed_half_rows(members: usize) -> f64 {
    let m = members.max(1) as f64;
    (0.25 + m.ln() * 0.05).clamp(0.28, 0.46)
}

/// Corners of a collapsed-cluster triangle in layout space: the apex at the
/// node, the base to the right.
pub fn triangle(x: f64, y: f64, half_rows: f64) -> [(f64, f64); 3] {
    [
        (x, y),
        (x + COLLAPSED_SPREAD, y - half_rows),
        (x + COLLAPSED_SPREAD, y + half_rows),
    ]
}

// ---------------------------------------------------------------------------
// Animation
// ---------------------------------------------------------------------------

/// Node positions keyed by animation identity.
pub type Positions = HashMap<u64, (f64, f64)>;

/// Key a list of positions (parallel to `l.nodes`) by animation identity. This
/// is how an interrupted animation hands its current positions to the next one.
pub fn positions_map(l: &TreeLayout, positions: &[(f64, f64)]) -> Positions {
    l.nodes
        .iter()
        .zip(positions.iter())
        .map(|(n, p)| (n.key, *p))
        .collect()
}

/// Remember where every node of a layout is, to animate away from it.
pub fn snapshot(l: &TreeLayout) -> Positions {
    let target: Vec<(f64, f64)> = l.nodes.iter().map(|n| (n.x, n.y)).collect();
    positions_map(l, &target)
}

/// The ease used for the move: slow at both ends, which reads as a physical
/// object settling rather than as a linear slide.
pub fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Positions of every node of `target` part-way between `from` and `target`.
///
/// `t` is the raw fraction of the animation that has elapsed; it is eased here.
/// A node that has no counterpart in `from` (a genome that just appeared)
/// starts at its parent's *starting* position, so a whole new subtree grows out
/// of a fixed place instead of chasing its parent across the screen. Every node
/// therefore interpolates between two fixed points and never overshoots.
pub fn animated(target: &TreeLayout, from: &Positions, t: f64) -> Vec<(f64, f64)> {
    // Parents come before children in `target.nodes`, so one pass resolves the
    // start of a node whose parent was itself unmatched.
    let mut starts: Vec<(f64, f64)> = Vec::with_capacity(target.nodes.len());
    for n in &target.nodes {
        let s = match from.get(&n.key) {
            Some(p) => *p,
            None => match n.parent {
                Some(p) if p < starts.len() => starts[p],
                _ => (n.x, n.y),
            },
        };
        starts.push(s);
    }
    let t = smoothstep(t);
    starts
        .iter()
        .zip(target.nodes.iter())
        .map(|(s, n)| (s.0 + (n.x - s.0) * t, s.1 + (n.y - s.1) * t))
        .collect()
}

// ---------------------------------------------------------------------------
// Scale bar
// ---------------------------------------------------------------------------

/// A scale bar: a layout width and the real distance it stands for.
#[derive(Debug, Clone, PartialEq)]
pub struct ScaleBar {
    pub layout_width: f64,
    pub value: f64,
    pub unit: &'static str,
    pub label: String,
}

/// Pick a scale bar for a region that is about `target` layout units wide.
///
/// The value is rounded down to 1, 2 or 5 times a power of ten, so the label is
/// a number a reader can divide by in their head.
pub fn scale_bar(region: &Region, target: f64) -> Option<ScaleBar> {
    let per = region.units_per_layout();
    let usable = per.is_finite() && per > 0.0 && target.is_finite() && target > 0.0;
    if !usable {
        return None;
    }
    let value = nice_below(per * target)?;
    let width = value / per;
    Some(ScaleBar {
        layout_width: width,
        value,
        unit: region.unit,
        label: format!(
            "{}{} {}",
            fmt::approx_mark(region.approximate || fmt::is_approximate_focus(region.focus)),
            fmt::snps(value),
            region.unit
        ),
    })
}

/// The largest 1/2/5 x 10^n that is not greater than `v`.
pub fn nice_below(v: f64) -> Option<f64> {
    if !v.is_finite() || v <= 0.0 {
        return None;
    }
    let e = v.log10().floor();
    let base = 10f64.powf(e);
    let m = v / base;
    let step = if m >= 5.0 {
        5.0
    } else if m >= 2.0 {
        2.0
    } else {
        1.0
    };
    Some(step * base)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn leaf(name: &str, len: f64) -> DisplayNode {
        DisplayNode::new(
            NodeKind::Genome {
                id: 0,
                name: name.to_string(),
            },
            len,
            2,
            false,
        )
    }

    /// Backbone of two clusters: `a` collapsed with 3 members, `b` opened into
    /// two genomes.
    fn sample() -> DisplayNode {
        let mut b = DisplayNode::new(
            NodeKind::ClusterRoot {
                id: 1,
                members: 2,
                focus: 2,
            },
            0.4,
            0,
            true,
        );
        b.children = vec![leaf("b1", 0.1), leaf("b2", 0.3)];
        let mut root = DisplayNode::new(NodeKind::Root, 0.0, 0, true);
        root.children = vec![
            DisplayNode::new(
                NodeKind::ClusterCollapsed {
                    id: 0,
                    members: 3,
                    focus: 0,
                },
                0.2,
                0,
                true,
            ),
            b,
        ];
        root
    }

    #[test]
    fn every_leaf_gets_its_own_row_in_traversal_order() {
        let l = layout(&sample());
        assert_eq!(l.rows, 3.0);
        assert_eq!(l.leaves.len(), 3);
        let rows: Vec<f64> = l.leaves.iter().map(|&i| l.nodes[i].y).collect();
        assert_eq!(rows, vec![0.0, 1.0, 2.0]);
    }

    #[test]
    fn an_internal_node_sits_at_the_midpoint_of_its_subtree() {
        let l = layout(&sample());
        let b = l
            .nodes
            .iter()
            .find(|n| matches!(n.kind, NodeKind::ClusterRoot { id: 1, .. }))
            .expect("cluster root");
        assert_eq!(b.row_min, 1.0);
        assert_eq!(b.row_max, 2.0);
        assert_eq!(b.y, 1.5);
        assert_eq!(b.leaf_count, 2);
    }

    #[test]
    fn parents_come_before_children_so_the_animation_can_chain() {
        let l = layout(&sample());
        for n in &l.nodes {
            if let Some(p) = n.parent {
                assert!(p < n.idx, "node {} has a later parent {}", n.idx, p);
            }
        }
        assert_eq!(l.root, 0);
        // Every node got a key, including the terminals.
        assert!(l.nodes.iter().all(|n| n.key != 0));
    }

    #[test]
    fn an_opened_cluster_gets_its_own_normalised_region() {
        let l = layout(&sample());
        assert_eq!(l.regions.len(), 2, "backbone plus one opened cluster");
        let c = &l.regions[1];
        assert_eq!(c.cluster, Some(1));
        assert_eq!(c.unit, "SNPs");
        // The cluster's deepest leaf is 0.3 units from its root, and the region
        // is CLUSTER_W wide, so 0.3 maps onto the full width.
        assert!((c.raw_max - 0.3).abs() < 1e-9);
        let b2 = l
            .nodes
            .iter()
            .find(|n| n.kind.name() == Some("b2"))
            .unwrap();
        assert!((b2.x - (c.x0 + CLUSTER_W)).abs() < 1e-9);
    }

    #[test]
    fn the_backbone_region_spans_the_backbone_only() {
        let l = layout(&sample());
        let bb = &l.regions[0];
        // Deepest backbone leaf is the cluster root at 0.4, not b2 at 0.7.
        assert!((bb.raw_max - 0.4).abs() < 1e-9);
        assert_eq!(bb.unit, "SNPs");
        assert!(bb.approximate, "the backbone is a sketch estimate");
        assert_eq!(bb.cluster, None);
        assert!(l.x_max >= c_x_max(&l));
    }

    fn c_x_max(l: &TreeLayout) -> f64 {
        l.regions.iter().map(|r| r.x0 + r.width).fold(0.0, f64::max)
    }

    #[test]
    fn a_zero_length_tree_does_not_divide_by_zero() {
        let mut root = DisplayNode::new(NodeKind::Root, 0.0, 0, true);
        root.children = vec![leaf("x", 0.0), leaf("y", 0.0)];
        let l = layout(&root);
        assert_eq!(l.regions[0].units_per_layout(), 0.0);
        assert_eq!(l.nodes.iter().filter(|n| n.x != 0.0).count(), 0);
        assert!(scale_bar(&l.regions[0], 0.2).is_none());
    }

    #[test]
    fn negative_and_infinite_lengths_are_clamped_not_propagated() {
        let mut root = DisplayNode::new(NodeKind::Root, 0.0, 0, true);
        root.children = vec![
            leaf("x", -5.0),
            leaf("y", f64::NAN),
            leaf("z", f64::INFINITY),
        ];
        let l = layout(&root);
        assert!(l.nodes.iter().all(|n| n.x.is_finite() && n.x >= 0.0));
    }

    #[test]
    fn a_childless_internal_node_still_occupies_a_row() {
        // `Tree::default()` is an empty internal node; it must not disappear.
        let root = DisplayNode::new(NodeKind::Internal, 0.0, 0, true);
        let l = layout(&root);
        assert_eq!(l.nodes.len(), 1);
        assert_eq!(l.rows, 1.0);
        assert_eq!(l.leaves, vec![0]);
    }

    #[test]
    fn region_lookup_prefers_the_cluster_over_the_backbone() {
        let l = layout(&sample());
        let c = l.regions[1].clone();
        assert_eq!(
            l.region_at(c.x0 + CLUSTER_W * 0.5).unwrap().cluster,
            Some(1)
        );
        // On the backbone, left of every cluster.
        assert_eq!(l.region_at(-0.1).unwrap().cluster, None);
    }

    #[test]
    fn collapsed_triangles_grow_with_the_member_count_but_stay_apart() {
        assert!(collapsed_half_rows(2) < collapsed_half_rows(500));
        for m in [1usize, 2, 10, 41, 88, 1000, 100_000] {
            assert!(collapsed_half_rows(m) < 0.5, "{m} members overflow the row");
            assert!(collapsed_half_rows(m) > 0.2);
        }
        let t = triangle(1.0, 2.0, 0.4);
        assert_eq!(t[0], (1.0, 2.0));
        assert!((t[1].0 - t[2].0).abs() < 1e-12);
        assert!(t[1].1 < t[2].1);
        assert!((t[1].0 - 1.0 - COLLAPSED_SPREAD).abs() < 1e-12);
    }

    #[test]
    fn the_scale_bar_picks_a_round_number_and_matches_the_region() {
        let r = Region {
            x0: 0.0,
            width: 1.0,
            raw_max: 4700.0,
            unit: "SNPs",
            cluster: None,
            focus: 0,
            approximate: true,
        };
        let b = scale_bar(&r, 0.25).expect("a positive region always gives a bar");
        assert!((b.value - 1000.0).abs() < 1e-9, "got {}", b.value);
        assert!(b.layout_width <= 0.25 + 1e-12);
        assert_eq!(b.label, "≈ 1\u{2009}000 SNPs");
        // The bar's width, converted back, is exactly its value.
        assert!((b.layout_width * r.units_per_layout() - b.value).abs() < 1e-12);
    }

    #[test]
    fn a_focus_1_scale_bar_says_it_is_approximate() {
        let r = Region {
            x0: 0.0,
            width: 1.0,
            raw_max: 12.0,
            unit: fmt::unit_label(1),
            cluster: Some(3),
            focus: 1,
            approximate: true,
        };
        let b = scale_bar(&r, 0.3).unwrap();
        assert!(b.label.starts_with('≈'), "got {:?}", b.label);
        assert_eq!(b.unit, "SNPs");
        // An exact cluster's bar is a plain count of SNPs.
        let exact = Region {
            focus: 2,
            approximate: false,
            ..r
        };
        assert_eq!(scale_bar(&exact, 0.3).unwrap().label, "2 SNPs");
    }

    /// Two opened exact clusters, one ten times wider than the other.
    fn two_exact_clusters(approx_second: bool) -> DisplayNode {
        let mut a = DisplayNode::new(
            NodeKind::ClusterRoot {
                id: 0,
                members: 2,
                focus: 2,
            },
            100.0,
            2,
            false,
        );
        a.children = vec![leaf("a1", 1.0), leaf("a2", 2.0)];
        let mut b = DisplayNode::new(
            NodeKind::ClusterRoot {
                id: 1,
                members: 2,
                focus: if approx_second { 1 } else { 2 },
            },
            100.0,
            2,
            approx_second,
        );
        b.children = vec![leaf("b1", 10.0), leaf("b2", 20.0)];
        let mut root = DisplayNode::new(NodeKind::Root, 0.0, 0, true);
        root.children = vec![a, b];
        root
    }

    #[test]
    fn exact_clusters_share_one_snp_scale() {
        let l = layout(&two_exact_clusters(false));
        let ra = l.regions.iter().find(|r| r.cluster == Some(0)).unwrap();
        let rb = l.regions.iter().find(|r| r.cluster == Some(1)).unwrap();
        assert_eq!(ra.units_per_layout(), rb.units_per_layout());
        // So the tight cluster is drawn ten times narrower than the loose one.
        let x = |name: &str| {
            let n = l
                .nodes
                .iter()
                .find(|n| n.kind.name() == Some(name))
                .unwrap();
            n.x - l.regions[n.region].x0
        };
        assert!((x("b2") / x("a2") - 10.0).abs() < 1e-9);
        // Nodes know their cluster and their depth below its root.
        let a1 = l
            .nodes
            .iter()
            .find(|n| n.kind.name() == Some("a1"))
            .unwrap();
        assert_eq!(a1.cluster, Some(0));
        assert_eq!(a1.depth, 1);
        assert_eq!(l.nodes[l.root].cluster, None);
        // The edge into an exact cluster is still a backbone estimate.
        let root_a = l
            .nodes
            .iter()
            .find(|n| matches!(n.kind, NodeKind::ClusterRoot { id: 0, .. }))
            .unwrap();
        assert_eq!(root_a.focus, 2);
        assert_eq!(root_a.edge_focus, 0);
        assert!(root_a.edge_approximate && !root_a.approximate);
        assert_eq!(a1.edge_focus, 2);
    }

    #[test]
    fn an_approximate_cluster_keeps_a_scale_of_its_own() {
        let l = layout(&two_exact_clusters(true));
        let ra = l.regions.iter().find(|r| r.cluster == Some(0)).unwrap();
        let rb = l.regions.iter().find(|r| r.cluster == Some(1)).unwrap();
        assert!(rb.approximate && !ra.approximate);
        // Each fills the cluster width with its own deepest leaf.
        assert!((ra.raw_max - 2.0).abs() < 1e-9);
        assert!((rb.raw_max - 20.0).abs() < 1e-9);
    }

    #[test]
    fn scaling_lengths_reaches_every_branch() {
        let mut t = sample();
        t.scale_lengths(1000.0);
        assert!((t.children[0].len - 200.0).abs() < 1e-9);
        assert!((t.children[1].children[1].len - 300.0).abs() < 1e-9);
        // A nonsense factor zeroes rather than poisoning the layout.
        t.scale_lengths(f64::NAN);
        assert_eq!(t.children[1].children[1].len, 0.0);
    }

    #[test]
    fn nice_below_lands_on_one_two_or_five() {
        assert_eq!(nice_below(9.9), Some(5.0));
        assert_eq!(nice_below(5.0), Some(5.0));
        assert_eq!(nice_below(4.9), Some(2.0));
        assert_eq!(nice_below(1.0), Some(1.0));
        assert!((nice_below(0.0037).unwrap() - 0.002).abs() < 1e-12);
        assert_eq!(nice_below(0.0), None);
        assert_eq!(nice_below(-1.0), None);
    }

    #[test]
    fn animation_starts_at_the_old_place_and_ends_at_the_new_one() {
        let l = layout(&sample());
        // Swap the two backbone children: the same leaves, a different shape.
        let mut moved = sample();
        moved.children.swap(0, 1);
        let l2 = layout(&moved);
        let from = snapshot(&l);
        let start = animated(&l2, &from, 0.0);
        let end = animated(&l2, &from, 1.0);
        assert_eq!(start.len(), l2.nodes.len());
        for (i, n) in l2.nodes.iter().enumerate() {
            if let Some(&(fx, fy)) = from.get(&n.key) {
                assert!((start[i].0 - fx).abs() < 1e-9);
                assert!((start[i].1 - fy).abs() < 1e-9);
            }
            assert!((end[i].0 - n.x).abs() < 1e-9);
            assert!((end[i].1 - n.y).abs() < 1e-9);
        }
        // Halfway is not the end yet.
        let mid = animated(&l2, &from, 0.5);
        assert!(mid
            .iter()
            .zip(end.iter())
            .any(|(m, e)| (m.1 - e.1).abs() > 1e-9 || (m.0 - e.0).abs() > 1e-9));
    }

    #[test]
    fn a_new_node_grows_out_of_its_parent_instead_of_popping_in() {
        let mut a = DisplayNode::new(NodeKind::Root, 0.0, 0, true);
        a.children = vec![leaf("x", 1.0)];
        let la = layout(&a);
        let from = snapshot(&la);

        let mut b = DisplayNode::new(NodeKind::Root, 0.0, 0, true);
        let mut inner = DisplayNode::new(NodeKind::Internal, 0.5, 2, false);
        inner.children = vec![leaf("x", 0.5), leaf("y", 0.5)];
        b.children = vec![inner];
        let lb = layout(&b);

        let start = animated(&lb, &from, 0.0);
        let y = lb
            .nodes
            .iter()
            .position(|n| n.kind.name() == Some("y"))
            .unwrap();
        let parent = lb.nodes[y].parent.unwrap();
        assert!((start[y].0 - start[parent].0).abs() < 1e-9);
        assert!((start[y].1 - start[parent].1).abs() < 1e-9);
        // "x" kept its identity, so it starts where it was.
        let x = lb
            .nodes
            .iter()
            .position(|n| n.kind.name() == Some("x"))
            .unwrap();
        assert!((start[x].0 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn smoothstep_is_clamped_and_symmetric() {
        assert_eq!(smoothstep(-1.0), 0.0);
        assert_eq!(smoothstep(2.0), 1.0);
        assert_eq!(smoothstep(0.0), 0.0);
        assert_eq!(smoothstep(1.0), 1.0);
        assert!((smoothstep(0.5) - 0.5).abs() < 1e-12);
        for i in 0..10 {
            let t = i as f64 / 10.0;
            assert!((smoothstep(t) + smoothstep(1.0 - t) - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn keys_are_stable_per_leaf_and_change_with_the_shape() {
        let a = node_key(&NodeKind::Internal, &[10, 20]);
        let b = node_key(&NodeKind::Internal, &[20, 10]);
        assert_ne!(a, b, "children are visited in row order");
        assert_eq!(a, node_key(&NodeKind::Internal, &[10, 20]));
        assert_ne!(
            node_key(
                &NodeKind::Genome {
                    id: 1,
                    name: "x".into()
                },
                &[]
            ),
            node_key(
                &NodeKind::Genome {
                    id: 2,
                    name: "y".into()
                },
                &[]
            )
        );
        // The genome id must not leak into the key: a renumbered genome is the
        // same leaf and should animate, not jump.
        assert_eq!(
            node_key(
                &NodeKind::Genome {
                    id: 1,
                    name: "x".into()
                },
                &[]
            ),
            node_key(
                &NodeKind::Genome {
                    id: 9,
                    name: "x".into()
                },
                &[]
            )
        );
    }

    #[test]
    fn from_core_maps_names_to_ids_and_keeps_the_focus() {
        let t = Tree::Internal {
            len: 0.0,
            children: vec![
                Tree::Leaf {
                    name: "a".into(),
                    len: 0.1,
                },
                Tree::Leaf {
                    name: "b".into(),
                    len: 0.2,
                },
            ],
        };
        let mut ids = HashMap::new();
        ids.insert("a".to_string(), 7u32);
        let d = from_core(&t, 2, false, &ids);
        assert_eq!(d.children.len(), 2);
        assert_eq!(
            d.children[0].kind,
            NodeKind::Genome {
                id: 7,
                name: "a".into()
            }
        );
        // "b" is not in the map: it gets the sentinel rather than a wrong id.
        assert_eq!(
            d.children[1].kind,
            NodeKind::Genome {
                id: u32::MAX,
                name: "b".into()
            }
        );
        assert_eq!(d.children[0].focus, 2);
        assert!(!d.children[0].approximate);
    }

    #[test]
    fn searching_finds_a_leaf_by_name() {
        let l = layout(&sample());
        assert_eq!(l.find("b2"), l.by_name.get("b2").copied());
        assert!(l.find("nope").is_none());
    }

    #[test]
    fn a_thousand_leaves_lay_out_without_recursion_problems() {
        let mut root = DisplayNode::new(NodeKind::Root, 0.0, 0, true);
        for i in 0..1000 {
            root.children.push(leaf(&format!("g{i}"), (i % 7) as f64));
        }
        let l = layout(&root);
        assert_eq!(l.rows, 1000.0);
        assert_eq!(l.leaves.len(), 1000);
        assert_eq!(l.by_name.len(), 1000);
    }
}
