//! Random binary trees with Poisson branch lengths (section 10.1).
//!
//! A cluster's pedigree is built by repeatedly picking two of the still-free
//! nodes and joining them under a new parent, which is the standard coalescent
//! construction of a random binary tree. Branch lengths are then drawn
//! independently from a Poisson distribution, and they *are* the SNP counts: the
//! same number decides how many substitutions the branch applies and what the
//! truth Newick says, so the two cannot drift apart.

use infinity_core::tree::{to_newick, Tree};

use crate::rng::Rng;

/// One node of a pedigree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PedigreeNode {
    /// Children, empty for a leaf. Always zero or two: the tree is binary.
    pub children: Vec<usize>,
    /// Which leaf this is, by index, when it is a leaf.
    pub leaf: Option<usize>,
    /// SNPs on the branch from the parent. Zero for the root.
    pub snps: u64,
}

/// A rooted binary tree over `n` leaves.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Pedigree {
    pub nodes: Vec<PedigreeNode>,
    pub root: usize,
    pub leaves: usize,
}

impl Pedigree {
    /// Draw a random binary topology with `leaves` leaves, then a Poisson branch
    /// length for every edge below the root.
    ///
    /// `leaves == 1` gives a single node that is both root and leaf; its branch
    /// length is zero because there is no branch.
    pub fn random(rng: &mut Rng, leaves: usize, mean_snps: f64) -> Self {
        assert!(leaves >= 1, "a pedigree needs at least one leaf");
        let mut nodes: Vec<PedigreeNode> = (0..leaves)
            .map(|i| PedigreeNode {
                children: Vec::new(),
                leaf: Some(i),
                snps: 0,
            })
            .collect();
        if leaves == 1 {
            return Self {
                nodes,
                root: 0,
                leaves,
            };
        }
        let mut active: Vec<usize> = (0..leaves).collect();
        while active.len() > 1 {
            let a = active.swap_remove(rng.below_usize(active.len()));
            let b = active.swap_remove(rng.below_usize(active.len()));
            let id = nodes.len();
            // Children in index order: the topology is already fixed by which
            // pairs joined, and a canonical left/right keeps the Newick stable.
            nodes.push(PedigreeNode {
                children: vec![a.min(b), a.max(b)],
                leaf: None,
                snps: 0,
            });
            active.push(id);
        }
        let root = active[0];
        // Branch lengths in node order, root excepted. The draw order is part of
        // the deterministic recipe: changing it changes every dataset.
        for (id, node) in nodes.iter_mut().enumerate() {
            if id != root {
                node.snps = rng.poisson(mean_snps) as u64;
            }
        }
        Self {
            nodes,
            root,
            leaves,
        }
    }

    /// Node id of each leaf, indexed by leaf number.
    pub fn leaf_nodes(&self) -> Vec<usize> {
        let mut out = vec![usize::MAX; self.leaves];
        for (id, node) in self.nodes.iter().enumerate() {
            if let Some(i) = node.leaf {
                out[i] = id;
            }
        }
        out
    }

    /// Parent of every node, `usize::MAX` for the root.
    pub fn parents(&self) -> Vec<usize> {
        let mut parent = vec![usize::MAX; self.nodes.len()];
        for (id, node) in self.nodes.iter().enumerate() {
            for &c in &node.children {
                parent[c] = id;
            }
        }
        parent
    }

    /// Total SNPs from the root to every node, and the number of edges.
    ///
    /// Walked from the root down: a parent always has a higher index than its
    /// children here, because internal nodes are created after the nodes they
    /// join.
    pub fn depths(&self) -> (Vec<u64>, Vec<u32>) {
        let parent = self.parents();
        let mut snps = vec![0u64; self.nodes.len()];
        let mut edges = vec![0u32; self.nodes.len()];
        for id in (0..self.nodes.len()).rev() {
            let p = parent[id];
            if p != usize::MAX {
                snps[id] = snps[p] + self.nodes[id].snps;
                edges[id] = edges[p] + 1;
            }
        }
        (snps, edges)
    }

    /// Prepare an ancestor walk once, so that many pairwise distances can be
    /// taken without recomputing it. Building the matrix is quadratic in the
    /// leaf count; this keeps the per-pair cost at the depth of the tree.
    pub fn distance_index(&self) -> DistanceIndex {
        let parent = self.parents();
        let (snps, edges) = self.depths();
        DistanceIndex {
            parent,
            snps,
            edges,
        }
    }

    /// Convert to the engine's tree type, labelling leaves from `names`.
    pub fn to_core_tree(&self, names: &[String]) -> Tree {
        self.build(self.root, names)
    }

    fn build(&self, node: usize, names: &[String]) -> Tree {
        let n = &self.nodes[node];
        let len = n.snps as f64;
        match n.leaf {
            Some(i) => Tree::Leaf {
                name: names[i].clone(),
                len,
            },
            None => Tree::Internal {
                len,
                children: n.children.iter().map(|&c| self.build(c, names)).collect(),
            },
        }
    }

    /// The truth Newick for this pedigree, branch lengths in SNPs.
    pub fn to_newick(&self, names: &[String]) -> String {
        to_newick(&self.to_core_tree(names))
    }
}

/// A pedigree flattened for fast lowest-common-ancestor walks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DistanceIndex {
    parent: Vec<usize>,
    snps: Vec<u64>,
    edges: Vec<u32>,
}

impl DistanceIndex {
    /// Sum of branch lengths on the path between two nodes.
    pub fn distance(&self, a: usize, b: usize) -> u64 {
        let (mut x, mut y) = (a, b);
        while self.edges[x] > self.edges[y] {
            x = self.parent[x];
        }
        while self.edges[y] > self.edges[x] {
            y = self.parent[y];
        }
        while x != y {
            x = self.parent[x];
            y = self.parent[y];
        }
        self.snps[a] + self.snps[b] - 2 * self.snps[x]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::tree::parse_newick;

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("leaf{i:02}")).collect()
    }

    #[test]
    fn single_leaf_pedigree() {
        let mut rng = Rng::from_seed(1);
        let p = Pedigree::random(&mut rng, 1, 3.0);
        assert_eq!(p.nodes.len(), 1);
        assert_eq!(p.root, 0);
        assert_eq!(p.nodes[0].snps, 0);
        assert_eq!(p.to_newick(&names(1)), "leaf00;");
        assert_eq!(p.distance_index().distance(0, 0), 0);
    }

    #[test]
    fn two_leaf_pedigree_distance_is_the_sum() {
        let mut rng = Rng::from_seed(2);
        let p = Pedigree::random(&mut rng, 2, 5.0);
        assert_eq!(p.nodes.len(), 3);
        let leaves: Vec<usize> = (0..3).filter(|&i| p.nodes[i].leaf.is_some()).collect();
        assert_eq!(leaves.len(), 2);
        let want = p.nodes[leaves[0]].snps + p.nodes[leaves[1]].snps;
        assert_eq!(p.distance_index().distance(leaves[0], leaves[1]), want);
        assert_eq!(p.nodes[p.root].snps, 0);
    }

    #[test]
    fn topology_is_binary_and_complete() {
        let ns = names(17);
        for seed in 0..20u64 {
            let mut rng = Rng::from_seed(seed);
            let p = Pedigree::random(&mut rng, 17, 3.0);
            assert_eq!(p.nodes.len(), 2 * 17 - 1);
            assert_eq!(p.nodes[p.root].snps, 0);
            let parent = p.parents();
            for (id, n) in p.nodes.iter().enumerate() {
                assert!(n.children.is_empty() || n.children.len() == 2, "node {id}");
                if id != p.root {
                    assert_ne!(parent[id], usize::MAX);
                }
            }
            // Every leaf index appears exactly once.
            let mut seen: Vec<usize> = p.nodes.iter().filter_map(|n| n.leaf).collect();
            seen.sort_unstable();
            assert_eq!(seen, (0..17).collect::<Vec<_>>());
            // The Newick round-trips through the engine's parser and keeps the
            // integer branch lengths exactly.
            let nwk = p.to_newick(&ns);
            let back = parse_newick(&nwk).unwrap();
            assert_eq!(to_newick(&back), nwk);
            assert_eq!(back.leaf_count(), 17);
            assert_eq!(
                back.total_length(),
                p.nodes.iter().map(|n| n.snps as f64).sum::<f64>()
            );
        }
    }

    #[test]
    fn deterministic_for_a_seed() {
        let mut a = Rng::from_seed(42);
        let mut b = Rng::from_seed(42);
        let ta = Pedigree::random(&mut a, 25, 4.0);
        let tb = Pedigree::random(&mut b, 25, 4.0);
        assert_eq!(ta, tb);
        let mut c = Rng::from_seed(43);
        assert_ne!(ta, Pedigree::random(&mut c, 25, 4.0));
    }

    #[test]
    fn distance_is_symmetric_and_metric() {
        let mut rng = Rng::from_seed(7);
        let p = Pedigree::random(&mut rng, 12, 3.0);
        let d = p.distance_index();
        let leaves = p.leaf_nodes();
        assert_eq!(leaves.len(), 12);
        for &a in &leaves {
            assert_eq!(d.distance(a, a), 0);
            for &b in &leaves {
                assert_eq!(d.distance(a, b), d.distance(b, a));
                for &c in &leaves {
                    assert!(d.distance(a, b) + d.distance(b, c) >= d.distance(a, c));
                }
            }
        }
    }

    #[test]
    fn mean_branch_length_tracks_the_requested_mean() {
        let mut rng = Rng::from_seed(99);
        let p = Pedigree::random(&mut rng, 200, 3.0);
        let branches: Vec<u64> = p
            .nodes
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != p.root)
            .map(|(_, n)| n.snps)
            .collect();
        let mean = branches.iter().sum::<u64>() as f64 / branches.len() as f64;
        assert!((mean - 3.0).abs() < 0.5, "mean = {mean}");
        assert_eq!(branches.len(), 2 * 200 - 2);
    }
}
