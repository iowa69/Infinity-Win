//! Annotated tree exports for other tree viewers: NEXUS with per-node
//! annotations for FigTree, and colour-strip / label-colour datasets for iTOL.
//!
//! Both describe the tree the window draws: every cluster opened, branch
//! lengths in SNPs (exact inside Focus 2 clusters, estimates elsewhere; see
//! `UiModel::display_tree`), and each genome in its cluster colour from
//! [`crate::colors::cluster_color`] on the light theme, which is what both
//! viewers draw on.
//!
//! Pure string building: no egui, no files. The export module writes the files.

use std::collections::HashSet;
use std::fmt::Write as _;

use crate::colors;
use crate::layout::{DisplayNode, NodeKind};
use crate::state::UiModel;

/// The whole tree with every cluster opened, in SNPs.
pub fn full_display_tree(model: &UiModel) -> Option<DisplayNode> {
    let all: HashSet<u32> = model.clusters.iter().map(|c| c.id).collect();
    model.display_tree(&all)
}

/// A NEXUS file with one rooted tree, readable by FigTree.
///
/// Leaves carry `[&cluster="c4",focus=2,snps_to_pivot=17,pivot="name"]`
/// (`snps_to_pivot` only once the cluster is exact), cluster roots carry
/// `[&cluster="c4",focus=2,members=30]`, and every branch that is an estimate
/// carries `approximate=true`, so a reader of the file can tell measured SNPs
/// from estimated ones exactly as the window does.
pub fn nexus(model: &UiModel) -> Result<String, String> {
    let tree = full_display_tree(model)
        .ok_or_else(|| "no tree yet: the backbone has not been built".to_string())?;
    let mut taxa: Vec<String> = Vec::new();
    collect_taxa(&tree, &mut taxa);
    let mut out = String::new();
    out.push_str("#NEXUS\n");
    out.push_str(
        "[Written by Infinity-Win. Branch lengths are in SNPs; \
                  approximate=true marks estimated lengths.]\n\n",
    );
    out.push_str("BEGIN TAXA;\n");
    let _ = writeln!(out, "  DIMENSIONS NTAX={};", taxa.len());
    out.push_str("  TAXLABELS\n");
    for t in &taxa {
        let _ = writeln!(out, "    {}", quote(t));
    }
    out.push_str("  ;\nEND;\n\n");
    out.push_str("BEGIN TREES;\n");
    let mut body = String::new();
    write_node(&tree, model, None, true, true, &mut body);
    let _ = writeln!(out, "  TREE infinity = [&R] {body};");
    out.push_str("END;\n");
    Ok(out)
}

fn collect_taxa(n: &DisplayNode, out: &mut Vec<String>) {
    match &n.kind {
        NodeKind::Genome { name, .. } => out.push(name.clone()),
        NodeKind::ClusterCollapsed { id, members, .. } => {
            out.push(collapsed_label(*id, *members));
        }
        _ => {
            for c in &n.children {
                collect_taxa(c, out);
            }
        }
    }
}

fn collapsed_label(id: u32, members: usize) -> String {
    format!("c{id} ({members} genomes)")
}

/// Write one node; `cluster` is the opened cluster the node lies in.
fn write_node(
    n: &DisplayNode,
    model: &UiModel,
    cluster: Option<u32>,
    parent_approx: bool,
    is_root: bool,
    out: &mut String,
) {
    let cluster = match n.kind {
        NodeKind::ClusterRoot { id, .. } => Some(id),
        _ => cluster,
    };
    let mut ann: Vec<String> = Vec::new();
    match &n.kind {
        NodeKind::Genome { id, name } => {
            out.push_str(&quote(name));
            let g = model.genome(*id);
            let cid = g.and_then(|g| g.cluster).or(cluster);
            if let Some(c) = cid {
                ann.push(format!("cluster=\"c{c}\""));
            }
            ann.push(format!("focus={}", n.focus));
            if let Some(c) = cid.and_then(|c| model.cluster(c)) {
                if c.focus >= crate::fmt::FOCUS_EXACT {
                    let snps = if *id == c.pivot {
                        Some(0)
                    } else {
                        c.numbers.get(id).map(|m| m.snps)
                    };
                    if let Some(s) = snps {
                        ann.push(format!("snps_to_pivot={s}"));
                    }
                }
                ann.push(format!("pivot=\"{}\"", escape_ann(&c.pivot_name)));
            }
        }
        NodeKind::ClusterCollapsed { id, members, focus } => {
            out.push_str(&quote(&collapsed_label(*id, *members)));
            ann.push(format!("cluster=\"c{id}\""));
            ann.push(format!("focus={focus}"));
            ann.push(format!("members={members}"));
        }
        NodeKind::ClusterRoot { id, members, focus } => {
            write_children(n, model, cluster, out);
            ann.push(format!("cluster=\"c{id}\""));
            ann.push(format!("focus={focus}"));
            ann.push(format!("members={members}"));
        }
        NodeKind::Internal | NodeKind::Root => {
            write_children(n, model, cluster, out);
            if let Some(c) = cluster {
                ann.push(format!("cluster=\"c{c}\""));
            }
        }
    }
    // The edge into a cluster root is a backbone branch: as approximate as
    // the backbone, whatever the cluster below it is.
    let edge_approx = match n.kind {
        NodeKind::ClusterRoot { .. } => parent_approx,
        _ => n.approximate,
    };
    if edge_approx && !is_root {
        ann.push("approximate=true".to_string());
    }
    if !ann.is_empty() {
        let _ = write!(out, "[&{}]", ann.join(","));
    }
    if !is_root {
        let _ = write!(out, ":{}", len(n.len));
    }
}

fn write_children(n: &DisplayNode, model: &UiModel, cluster: Option<u32>, out: &mut String) {
    out.push('(');
    for (i, c) in n.children.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        write_node(c, model, cluster, n.approximate, false, out);
    }
    out.push(')');
}

/// A branch length with enough digits to round-trip a half SNP, and no noise.
fn len(v: f64) -> String {
    if !v.is_finite() || v <= 0.0 {
        return "0".to_string();
    }
    let s = format!("{v:.6}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// A NEXUS/Newick name, single-quoted when it holds anything structural.
pub fn quote(name: &str) -> String {
    let plain = !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.' | '|' | '#'));
    if plain {
        name.to_string()
    } else {
        format!("'{}'", name.replace('\'', "''"))
    }
}

/// A string value inside a `[&...]` annotation.
fn escape_ann(s: &str) -> String {
    s.replace('"', "'").replace([']', '['], "")
}

/// `#rrggbb` for a colour.
pub fn hex(c: egui::Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b())
}

/// Genome name and cluster id of every clustered genome, by cluster then name.
fn clustered(model: &UiModel) -> Vec<(String, u32)> {
    let mut v: Vec<(String, u32)> = Vec::new();
    for c in &model.clusters {
        for name in &c.member_names {
            v.push((name.clone(), c.id));
        }
    }
    v.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
    v
}

/// iTOL `DATASET_COLORSTRIP`: one coloured cell per genome, in its cluster
/// colour, with a legend. Drag the file onto a tree in iTOL.
pub fn itol_color_strip(model: &UiModel) -> String {
    let mut out = String::new();
    out.push_str("DATASET_COLORSTRIP\nSEPARATOR TAB\n");
    out.push_str("DATASET_LABEL\tInfinity-Win clusters\n");
    out.push_str("COLOR\t#3fb8bf\n");
    out.push_str("STRIP_WIDTH\t24\n");
    out.push_str("LEGEND_TITLE\tCluster\n");
    let ids: Vec<u32> = {
        let mut v: Vec<u32> = model.clusters.iter().map(|c| c.id).collect();
        v.sort_unstable();
        v
    };
    let shapes: Vec<&str> = ids.iter().map(|_| "1").collect();
    let cols: Vec<String> = ids
        .iter()
        .map(|&i| hex(colors::cluster_color(i, false)))
        .collect();
    let labels: Vec<String> = ids.iter().map(|i| format!("c{i}")).collect();
    let _ = writeln!(out, "LEGEND_SHAPES\t{}", shapes.join("\t"));
    let _ = writeln!(out, "LEGEND_COLORS\t{}", cols.join("\t"));
    let _ = writeln!(out, "LEGEND_LABELS\t{}", labels.join("\t"));
    out.push_str("DATA\n");
    for (name, c) in clustered(model) {
        let _ = writeln!(
            out,
            "{}\t{}\tc{c}",
            itol_id(&name),
            hex(colors::cluster_color(c, false))
        );
    }
    out
}

/// iTOL `TREE_COLORS`: every tip label in its cluster colour.
pub fn itol_label_colors(model: &UiModel) -> String {
    let mut out = String::new();
    out.push_str("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
    for (name, c) in clustered(model) {
        let _ = writeln!(
            out,
            "{}\tlabel\t{}\tbold",
            itol_id(&name),
            hex(colors::cluster_color(c, false))
        );
    }
    out
}

/// iTOL matches node ids against the tree's labels as written; a tab or a
/// newline would break the file, so those become spaces.
fn itol_id(name: &str) -> String {
    name.replace(['\t', '\n', '\r'], " ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::events::{ClusterFocusInfo, Event, MemberNumbers};

    fn model() -> UiModel {
        let mut m = UiModel::default();
        for i in 0..4u32 {
            m.apply(0.0, crate::state::scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2]), (1, 3, &[3])]),
        );
        m.apply(
            0.0,
            Event::BackboneTree {
                newick: "(g0:0.002,g3:0.004);".into(),
                pivots: vec!["g0".into(), "g3".into()],
            },
        );
        m.apply(
            0.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(g0:1,(g1:2.5,g2:4):1);".into(),
                members: (0..3)
                    .map(|i| MemberNumbers {
                        id: i,
                        name: format!("g{i}"),
                        diff_pieces: 0,
                        snps: u64::from(i) * 7,
                        indels: 0,
                        uncovered_bp: 0,
                        novel_bp: 0,
                    })
                    .collect(),
            }),
        );
        m
    }

    #[test]
    fn the_nexus_tree_has_every_genome_and_figtree_annotations() {
        let m = model();
        let nx = nexus(&m).expect("a tree exists");
        assert!(nx.starts_with("#NEXUS\n"));
        assert!(nx.contains("DIMENSIONS NTAX=4;"), "{nx}");
        assert!(nx.contains("TREE infinity = [&R] ("));
        // A member of an exact cluster: its SNPs to the pivot, and the pivot.
        assert!(
            nx.contains("g1[&cluster=\"c0\",focus=2,snps_to_pivot=7,pivot=\"g0\"]:2.5"),
            "{nx}"
        );
        // The pivot is 0 SNPs from itself.
        assert!(nx.contains("g0[&cluster=\"c0\",focus=2,snps_to_pivot=0,pivot=\"g0\"]:1"));
        // The cluster root is annotated and its backbone edge is ≈ SNPs.
        assert!(
            nx.contains("[&cluster=\"c0\",focus=2,members=3,approximate=true]:2"),
            "{nx}"
        );
        // The tree line parses back with the engine's own reader, which skips
        // the bracketed annotations, and has all four genomes.
        let line = nx.lines().find(|l| l.contains("TREE infinity")).unwrap();
        let nwk = line.split("[&R]").nth(1).unwrap().trim();
        let t = infinity_core::tree::parse_newick(nwk).expect("parses");
        let mut names = t.leaf_names();
        names.sort();
        assert_eq!(names, vec!["g0", "g1", "g2", "g3"]);
    }

    #[test]
    fn a_model_without_a_backbone_has_no_nexus() {
        assert!(nexus(&UiModel::default()).is_err());
    }

    #[test]
    fn names_are_quoted_only_when_they_must_be() {
        assert_eq!(quote("sim_c00_s001"), "sim_c00_s001");
        assert_eq!(quote("E. coli K-12"), "'E. coli K-12'");
        assert_eq!(quote("it's"), "'it''s'");
        assert_eq!(quote(""), "''");
        assert_eq!(len(2.5), "2.5");
        assert_eq!(len(3.0), "3");
        assert_eq!(len(-1.0), "0");
    }

    #[test]
    fn itol_files_colour_every_genome_like_the_window() {
        let m = model();
        let strip = itol_color_strip(&m);
        assert!(strip.starts_with("DATASET_COLORSTRIP\nSEPARATOR TAB\n"));
        let c0 = hex(colors::cluster_color(0, false));
        let c1 = hex(colors::cluster_color(1, false));
        assert!(
            strip.contains(&format!("LEGEND_COLORS\t{c0}\t{c1}\n")),
            "{strip}"
        );
        assert!(strip.contains(&format!("g1\t{c0}\tc0\n")));
        assert!(strip.contains(&format!("g3\t{c1}\tc1\n")));
        let data_rows = strip.split("DATA\n").nth(1).unwrap().lines().count();
        assert_eq!(data_rows, 4);
        let labels = itol_label_colors(&m);
        assert!(labels.contains(&format!("g2\tlabel\t{c0}\tbold\n")));
        assert_eq!(hex(egui::Color32::from_rgb(255, 0, 16)), "#ff0010");
    }
}
