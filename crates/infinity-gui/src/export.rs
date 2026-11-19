//! The export menu of section 8, in three kinds:
//!
//! * **Written by the window** from its own model: the tree as Newick, the
//!   annotated tree as NEXUS (FigTree) and as iTOL datasets, the clusters, the
//!   genomes, the distance matrices, the variants, and the tree picture as SVG
//!   (vector) and PNG (a screenshot of the tree panel).
//! * **Copied from the run's output folder**: the files only the engine writes
//!   (the Excel workbook, the square SNP matrices, VCF, the core-SNP
//!   alignment). The window does not reimplement them; it copies them from the
//!   folder it ran into or opened, and says plainly when an older run does not
//!   have them.
//!
//! The tab-separated files use the same columns and the same 1-based positions as
//! the engine's own output folder (section 9), so a file exported from the
//! window and the same file written by `infinity run` can be compared line by line.
//! The one thing the interface cannot reproduce is a member's contig names: they
//! are not in the event stream, so the member contig column carries the index.
//! When the run finished in this window, `RunResults` has them and is used
//! instead.

use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use infinity_core::engine::RunResults;
use infinity_core::tree::{graft, to_newick, Tree};

use crate::fmt;
use crate::pairwise::{self, MemberProfile};
use crate::state::UiModel;
use crate::svg::SvgOptions;

/// Above this many members per cluster, pairwise distances are taken from the
/// engine's matrix only. Deriving them here is quadratic in the members, and an
/// export must not freeze the window (section 8).
pub const DERIVED_PAIR_CAP: usize = 300;

/// What the export menu offers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Item {
    Tree,
    Nexus,
    Itol,
    Clusters,
    Genomes,
    Distances,
    Variants,
    SnpGroups,
    Svg,
    Png,
    Excel,
    SnpMatrices,
    Vcf,
    SnpAlignment,
}

impl Item {
    pub fn label(&self) -> &'static str {
        match self {
            Item::Tree => "Tree (Newick)",
            Item::Nexus => "Annotated tree (NEXUS, for FigTree)",
            Item::Itol => "Cluster colours (iTOL datasets)",
            Item::Clusters => "Clusters (TSV)",
            Item::Genomes => "Genomes (TSV)",
            Item::Distances => "Distance matrices (TSV)",
            Item::Variants => "Variants (TSV)",
            Item::SnpGroups => "SNP groups at the slider's threshold (CSV)",
            Item::Svg => "Tree picture (SVG, vector, on white)",
            Item::Png => "Tree picture (PNG, as on screen)",
            Item::Excel => "Excel workbook (from the run folder)",
            Item::SnpMatrices => "Square SNP matrices (from the run folder)",
            Item::Vcf => "VCF per cluster (from the run folder)",
            Item::SnpAlignment => "Core-SNP alignment (from the run folder)",
        }
    }

    /// A button-sized name.
    pub fn short_label(&self) -> &'static str {
        match self {
            Item::Tree => "tree.nwk",
            Item::Nexus => "NEXUS",
            Item::Itol => "iTOL",
            Item::Clusters => "clusters",
            Item::Genomes => "genomes",
            Item::Distances => "distances",
            Item::Variants => "variants",
            Item::SnpGroups => "SNP groups",
            Item::Svg => "SVG",
            Item::Png => "PNG",
            Item::Excel => "Excel",
            Item::SnpMatrices => "SNP matrices",
            Item::Vcf => "VCF",
            Item::SnpAlignment => "SNP alignment",
        }
    }

    /// Everything except the PNG, which needs a screenshot callback.
    pub fn files() -> Vec<Item> {
        Self::menu()
            .into_iter()
            .filter(|i| *i != Item::Png)
            .collect()
    }

    /// Every item, in the order the Export menu lists them.
    pub fn menu() -> Vec<Item> {
        Self::groups()
            .iter()
            .flat_map(|(_, items)| items.iter().copied())
            .collect()
    }

    /// The items grouped the way the Run panel shows them.
    pub fn groups() -> &'static [(&'static str, &'static [Item])] {
        &[
            ("trees", &[Item::Tree, Item::Nexus, Item::Itol]),
            (
                "tables",
                &[
                    Item::Clusters,
                    Item::Genomes,
                    Item::Distances,
                    Item::Variants,
                    Item::SnpGroups,
                ],
            ),
            (
                "from the run folder",
                &[
                    Item::Excel,
                    Item::SnpMatrices,
                    Item::Vcf,
                    Item::SnpAlignment,
                ],
            ),
            ("pictures", &[Item::Svg, Item::Png]),
        ]
    }

    /// The engine-written files this item copies, as `*` patterns. Empty for
    /// the items the window writes itself.
    pub fn run_files(&self) -> &'static [&'static str] {
        match self {
            Item::Excel => &["results.xlsx"],
            Item::SnpMatrices => &["cluster_*.snp_matrix.tsv", "snp_clusters.tsv"],
            Item::Vcf => &["cluster_*.vcf"],
            Item::SnpAlignment => &["cluster_*.core_snps.fasta"],
            _ => &[],
        }
    }
}

/// Everything an export may need besides the target folder.
#[derive(Clone, Copy)]
pub struct ExportCtx<'a> {
    pub model: &'a UiModel,
    pub results: Option<&'a RunResults>,
    /// The run's own output folder: where the engine-written files are copied
    /// from. `None` when no run or opened folder is behind the model.
    pub run_dir: Option<&'a Path>,
    /// The tree exactly as the panel shows it (opened clusters, selection,
    /// SNP groups), for the SVG. Without one the SVG draws the whole tree
    /// with every cluster opened.
    pub figure: Option<&'a crate::pubtree::Figure>,
    pub svg: SvgOptions,
    /// The SNP groups at the slider's threshold, when they are switched on.
    pub groups: Option<&'a crate::groups::GroupSet>,
}

impl<'a> ExportCtx<'a> {
    /// A context with only the model: what the tests and a bare export need.
    pub fn new(model: &'a UiModel, results: Option<&'a RunResults>) -> Self {
        Self {
            model,
            results,
            run_dir: None,
            figure: None,
            svg: SvgOptions::default(),
            groups: None,
        }
    }
}

/// What an export did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Report {
    pub written: Vec<PathBuf>,
    pub notes: Vec<String>,
    pub errors: Vec<String>,
}

impl Report {
    pub fn ok(&self) -> bool {
        self.errors.is_empty()
    }

    /// One line for the status bar.
    pub fn summary(&self) -> String {
        if self.errors.is_empty() {
            format!("exported {} file(s)", self.written.len())
        } else {
            format!(
                "exported {} file(s), {} problem(s): {}",
                self.written.len(),
                self.errors.len(),
                self.errors[0]
            )
        }
    }
}

/// Write every exportable file into `dir`, with only the model to go on.
#[cfg(test)]
pub fn write_all(dir: &Path, model: &UiModel, results: Option<&RunResults>) -> Report {
    write_all_with(dir, &ExportCtx::new(model, results))
}

/// Write every exportable file into `dir` (everything but the PNG). Files the
/// run folder does not have are a note, not an error: "everything" means
/// everything this run has.
pub fn write_all_with(dir: &Path, ctx: &ExportCtx<'_>) -> Report {
    let mut report = Report::default();
    if ctx.model.genome_count() == 0 {
        report.errors.push("nothing to export yet".to_string());
        return report;
    }
    if let Err(e) = std::fs::create_dir_all(dir) {
        report
            .errors
            .push(format!("cannot create {}: {e}", dir.display()));
        return report;
    }
    for item in Item::files() {
        match write_item(item, dir, ctx, &mut report.notes) {
            Ok(paths) => report.written.extend(paths),
            Err(e) if !item.run_files().is_empty() || item == Item::SnpGroups => {
                report.notes.push(e)
            }
            Err(e) => report.errors.push(format!("{}: {e}", item.label())),
        }
    }
    report
}

/// Write one item, which may be several files (one per cluster).
pub fn write_item(
    item: Item,
    dir: &Path,
    ctx: &ExportCtx<'_>,
    notes: &mut Vec<String>,
) -> Result<Vec<PathBuf>, String> {
    let model = ctx.model;
    let results = ctx.results;
    let mut out = Vec::new();
    match item {
        Item::Tree => {
            let p = dir.join("tree.nwk");
            write_file(&p, |w| tree_newick(w, model, results))?;
            out.push(p);
        }
        Item::Nexus => {
            let text = crate::annotate::nexus(model)?;
            let p = dir.join("tree.nexus");
            write_text(&p, &text)?;
            out.push(p);
        }
        Item::Itol => {
            let strip = dir.join("itol_cluster_strip.txt");
            write_text(&strip, &crate::annotate::itol_color_strip(model))?;
            let labels = dir.join("itol_label_colors.txt");
            write_text(&labels, &crate::annotate::itol_label_colors(model))?;
            out.push(strip);
            out.push(labels);
        }
        Item::Svg => {
            let owned;
            let f = match ctx.figure {
                Some(f) => f,
                None => {
                    let t = crate::annotate::full_display_tree(model).ok_or_else(|| {
                        "no tree yet: the backbone has not been built".to_string()
                    })?;
                    owned = crate::pubtree::Figure {
                        layout: crate::layout::layout(&t),
                        marks: crate::pubtree::marks_from_model(model),
                    };
                    &owned
                }
            };
            let p = dir.join("tree.svg");
            write_text(&p, &crate::svg::render(f, &ctx.svg))?;
            out.push(p);
        }
        Item::SnpGroups => {
            let set = ctx
                .groups
                .filter(|s| !s.covered.is_empty())
                .ok_or_else(|| {
                    "SNP groups: none yet (they need an exact cluster, and SNP groups switched on \
                 in the tree header)"
                        .to_string()
                })?;
            let p = dir.join(format!("snp_groups_le{}.csv", set.threshold));
            let text = crate::groups::csv(set, |g| {
                model
                    .genome(g)
                    .map(|r| r.name.clone())
                    .unwrap_or_else(|| format!("genome {g}"))
            });
            write_text(&p, &text)?;
            out.push(p);
        }
        Item::Clusters => {
            let p = dir.join("clusters.tsv");
            write_file(&p, |w| clusters_tsv(w, model))?;
            out.push(p);
        }
        Item::Genomes => {
            let p = dir.join("genomes.tsv");
            write_file(&p, |w| genomes_tsv(w, model))?;
            out.push(p);
        }
        Item::Variants => {
            for id in cluster_ids(model) {
                let p = dir.join(format!("cluster_{id}.variants.tsv"));
                match write_file(&p, |w| variants_tsv(w, model, results, id)) {
                    Ok(0) => {
                        std::fs::remove_file(&p).ok();
                    }
                    Ok(_) => out.push(p),
                    Err(e) => {
                        std::fs::remove_file(&p).ok();
                        return Err(e);
                    }
                }
            }
            if out.is_empty() {
                notes.push("no variants to export: no cluster has reached Focus 2".to_string());
            }
        }
        Item::Distances => {
            for id in cluster_ids(model) {
                let p = dir.join(format!("cluster_{id}.dist.tsv"));
                let derived = derived_needed(model, results, id);
                match write_file(&p, |w| distances_tsv(w, model, results, id)) {
                    Ok(0) => {
                        std::fs::remove_file(&p).ok();
                    }
                    Ok(_) => {
                        if derived {
                            notes.push(format!(
                                "cluster {id}: distances were recomputed from the variant lists, \
                                 not taken from the engine's matrix"
                            ));
                        }
                        out.push(p);
                    }
                    Err(e) => {
                        // A header-only file would look like a result.
                        std::fs::remove_file(&p).ok();
                        return Err(e);
                    }
                }
            }
            if out.is_empty() {
                notes.push("no distances to export yet".to_string());
            }
        }
        Item::Excel | Item::SnpMatrices | Item::Vcf | Item::SnpAlignment => {
            out = copy_from_run(item, dir, ctx.run_dir)?;
        }
        Item::Png => {
            // Needs a screenshot, which only the interface can take.
            return Err("the picture is exported from the tree panel".to_string());
        }
    }
    Ok(out)
}

/// Copy the engine-written files of `item` from the run folder into `dir`.
fn copy_from_run(item: Item, dir: &Path, run_dir: Option<&Path>) -> Result<Vec<PathBuf>, String> {
    let patterns = item.run_files();
    let what = patterns.join(", ");
    let Some(src) = run_dir else {
        return Err(format!(
            "{what}: there is no run folder to copy from; run the genomes or open a \
             finished output folder first"
        ));
    };
    let entries = std::fs::read_dir(src).map_err(|e| format!("reading {}: {e}", src.display()))?;
    let mut names: Vec<String> = entries
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().map(|t| t.is_file()).unwrap_or(false))
        .filter_map(|e| e.file_name().into_string().ok())
        .filter(|n| patterns.iter().any(|p| glob_match(p, n)))
        .collect();
    names.sort();
    if names.is_empty() {
        return Err(format!(
            "{what} not found in {}: this run was made by an engine that did not write \
             it yet; run again to get it",
            src.display()
        ));
    }
    let same_dir = match (src.canonicalize(), dir.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    };
    let mut out = Vec::with_capacity(names.len());
    for n in names {
        let to = dir.join(&n);
        if !same_dir {
            std::fs::copy(src.join(&n), &to).map_err(|e| format!("copying {n}: {e}"))?;
        }
        out.push(to);
    }
    Ok(out)
}

/// A file-name match with at most one `*`, which is all the export patterns
/// need (`cluster_*.vcf`).
pub fn glob_match(pattern: &str, name: &str) -> bool {
    match pattern.split_once('*') {
        None => pattern == name,
        Some((pre, suf)) => {
            name.len() >= pre.len() + suf.len() && name.starts_with(pre) && name.ends_with(suf)
        }
    }
}

fn write_text(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("writing {}: {e}", path.display()))
}

fn cluster_ids(model: &UiModel) -> Vec<u32> {
    let mut ids: Vec<u32> = model.clusters.iter().map(|c| c.id).collect();
    ids.sort_unstable();
    ids
}

fn derived_needed(model: &UiModel, results: Option<&RunResults>, cluster: u32) -> bool {
    results.and_then(|r| r.matrices.get(&cluster)).is_none()
        && !model.mapped_members(cluster).is_empty()
}

/// Run a writer into a fresh file, returning how many data rows it wrote.
fn write_file<F>(path: &Path, body: F) -> Result<usize, String>
where
    F: FnOnce(&mut RowCount) -> io::Result<()>,
{
    let file =
        std::fs::File::create(path).map_err(|e| format!("creating {}: {e}", path.display()))?;
    let mut w = BufWriter::new(file);
    let mut rc = RowCount::new(&mut w);
    body(&mut rc).map_err(|e| format!("writing {}: {e}", path.display()))?;
    rc.finish()
        .map_err(|e| format!("writing {}: {e}", path.display()))
}

/// A writer that counts the data rows it produced, so an empty export can be
/// removed instead of leaving a header-only file behind.
pub struct RowCount<'w> {
    w: &'w mut dyn Write,
    rows: usize,
}

impl<'w> RowCount<'w> {
    pub fn new(w: &'w mut dyn Write) -> Self {
        Self { w, rows: 0 }
    }

    pub fn line(&mut self, args: std::fmt::Arguments<'_>) -> io::Result<()> {
        self.rows += 1;
        writeln!(self.w, "{args}")
    }

    pub fn raw(&mut self, s: &str) -> io::Result<()> {
        self.w.write_all(s.as_bytes())
    }

    fn finish(&mut self) -> io::Result<usize> {
        self.w.flush()?;
        Ok(self.rows)
    }
}

// ---------------------------------------------------------------------------
// Content
// ---------------------------------------------------------------------------

/// The full tree: the backbone with every cluster subtree grafted in at its
/// pivot, exactly as `engine::build_full_tree` does.
pub fn full_tree(model: &UiModel, results: Option<&RunResults>) -> Option<Tree> {
    if let Some(r) = results {
        if r.full_tree.leaf_count() > 0 {
            return Some(r.full_tree.clone());
        }
    }
    let mut t = model.backbone.clone()?;
    for id in cluster_ids(model) {
        let c = model.cluster(id)?;
        if let Some(sub) = &c.tree {
            graft(&mut t, &c.pivot_name, sub.clone());
        }
    }
    Some(t)
}

pub fn tree_newick(
    w: &mut RowCount<'_>,
    model: &UiModel,
    results: Option<&RunResults>,
) -> io::Result<()> {
    match full_tree(model, results) {
        Some(t) => w.line(format_args!("{}", to_newick(&t))),
        None => Err(io::Error::new(
            io::ErrorKind::NotFound,
            "no tree yet: the backbone has not been built",
        )),
    }
}

pub fn clusters_tsv(w: &mut RowCount<'_>, model: &UiModel) -> io::Result<()> {
    w.raw("cluster_id\tpivot\tpivot_name\tsize\tfocus\tapproximate\n")?;
    for id in cluster_ids(model) {
        let Some(c) = model.cluster(id) else { continue };
        w.line(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            c.id,
            c.pivot,
            cell(&c.pivot_name),
            c.size(),
            c.focus,
            c.approximate
        ))?;
    }
    Ok(())
}

pub fn genomes_tsv(w: &mut RowCount<'_>, model: &UiModel) -> io::Result<()> {
    w.raw(
        "id\tname\tlength\tcontigs\tn50\tgc\tcluster\tis_pivot\tfocus\t\
         diff_pieces\tsnps\tindels\tuncovered_bp\tnovel_bp\n",
    )?;
    for g in &model.genomes {
        let n = g.numbers.as_ref();
        let num = |f: fn(&infinity_core::events::MemberNumbers) -> u64| -> String {
            n.map(|x| f(x).to_string())
                .unwrap_or_else(|| "NA".to_string())
        };
        w.line(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.id,
            cell(&g.name),
            g.length,
            g.contigs,
            g.n50,
            g.gc,
            g.cluster
                .map(|c| c.to_string())
                .unwrap_or_else(|| "NA".into()),
            g.is_pivot,
            g.focus,
            num(|m| m.diff_pieces),
            num(|m| m.snps),
            num(|m| m.indels),
            num(|m| m.uncovered_bp),
            num(|m| m.novel_bp),
        ))?;
    }
    Ok(())
}

pub fn variants_tsv(
    w: &mut RowCount<'_>,
    model: &UiModel,
    results: Option<&RunResults>,
    cluster: u32,
) -> io::Result<()> {
    w.raw(
        "genome\tgenome_name\tpivot_contig\tpivot_pos\tref\talt\ttype\t\
         genome_contig\tgenome_pos\tstrand\tflags\n",
    )?;
    let empty: Vec<String> = Vec::new();
    for gid in model.mapped_members(cluster) {
        let Some(m) = model.map(cluster, gid) else {
            continue;
        };
        let name = model.genome(gid).map(|g| g.name.as_str()).unwrap_or("");
        let member_contigs = results
            .and_then(|r| r.contig_names.get(&gid))
            .unwrap_or(&empty);
        for v in pairwise::sorted_variants(m) {
            w.line(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                gid,
                cell(name),
                name_or_index(&m.pivot_contigs, v.pivot_contig),
                fmt::one_based(v.pivot_start, v.kind),
                or_dash(&v.ref_allele),
                or_dash(&v.alt_allele),
                v.kind.as_str(),
                name_or_index(member_contigs, v.genome_contig),
                v.genome_start + 1,
                v.strand,
                v.flags_string(),
            ))?;
        }
    }
    Ok(())
}

pub fn distances_tsv(
    w: &mut RowCount<'_>,
    model: &UiModel,
    results: Option<&RunResults>,
    cluster: u32,
) -> io::Result<()> {
    w.raw("genome_a\tname_a\tgenome_b\tname_b\tsnps\tindels\tcomparable_bp\tsubs_per_site\n")?;
    if let Some(r) = results {
        if let Some(m) = r.matrices.get(&cluster) {
            let n = m.ids.len();
            for i in 0..n {
                for j in (i + 1)..n {
                    w.line(format_args!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        m.ids[i],
                        cell(&name_of(model, results, m.ids[i])),
                        m.ids[j],
                        cell(&name_of(model, results, m.ids[j])),
                        m.snps[i][j],
                        m.indels[i][j],
                        m.comparable[i][j],
                        dist_cell(m.comparable[i][j], m.dist[i][j]),
                    ))?;
                }
            }
            return Ok(());
        }
    }

    // No matrix: recompute 6.9 from the variant lists of the map events.
    let members = model.mapped_members(cluster);
    if members.len() > DERIVED_PAIR_CAP {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            format!(
                "cluster {cluster} has {} members, more than the {DERIVED_PAIR_CAP} this window \
                 will recompute; run with an output folder to get the engine's matrix",
                members.len()
            ),
        ));
    }
    let focus = model.cluster(cluster).map(|c| c.focus).unwrap_or(0);
    if focus < fmt::FOCUS_EXACT {
        return Ok(());
    }
    let mut profiles: Vec<MemberProfile> = Vec::with_capacity(members.len());
    for &gid in &members {
        if let Some(m) = model.map(cluster, gid) {
            profiles.push(pairwise::profile(m, gid));
        }
    }
    for i in 0..profiles.len() {
        for j in (i + 1)..profiles.len() {
            let p = pairwise::pair_from_profiles(&profiles[i], &profiles[j], focus);
            w.line(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                p.a,
                cell(&name_of(model, results, p.a)),
                p.b,
                cell(&name_of(model, results, p.b)),
                p.snps,
                p.indels,
                p.comparable,
                dist_cell(p.comparable, p.dist.unwrap_or(0.0)),
            ))?;
        }
    }
    Ok(())
}

fn dist_cell(comparable: u64, d: f64) -> String {
    if comparable == 0 {
        "NA".to_string()
    } else {
        format!("{d:.12}")
    }
}

fn name_of(model: &UiModel, results: Option<&RunResults>, id: u32) -> String {
    if let Some(g) = model.genome(id) {
        return g.name.clone();
    }
    results
        .and_then(|r| r.names.get(id as usize).cloned())
        .unwrap_or_else(|| format!("genome {id}"))
}

/// A cell that cannot break the tab-separated layout.
pub fn cell(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c == '\t' || c == '\n' || c == '\r' {
                ' '
            } else {
                c
            }
        })
        .collect()
}

fn or_dash(s: &str) -> &str {
    if s.is_empty() {
        "-"
    } else {
        s
    }
}

fn name_or_index(names: &[String], index: u32) -> String {
    names
        .get(index as usize)
        .cloned()
        .unwrap_or_else(|| index.to_string())
}

// ---------------------------------------------------------------------------
// The picture
// ---------------------------------------------------------------------------

/// Write a screenshot as a PNG, cropped to `crop` when one is given.
///
/// `pixels_per_point` converts the logical rectangle the tree panel occupies
/// into texels of the captured image.
pub fn write_png(
    path: &Path,
    image: &egui::ColorImage,
    crop: Option<(egui::Rect, f32)>,
) -> io::Result<()> {
    let img = match crop {
        Some((rect, ppp)) => image.region(&rect, Some(ppp)),
        None => image.clone(),
    };
    let [w, h] = img.size;
    if w == 0 || h == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the tree panel has no size to capture",
        ));
    }
    let mut rgba = Vec::with_capacity(w * h * 4);
    for p in &img.pixels {
        let [r, g, b, a] = p.to_array();
        rgba.extend_from_slice(&[r, g, b, a]);
    }
    let buf = image::RgbaImage::from_raw(w as u32, h as u32, rgba).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "the captured size does not match the pixels",
        )
    })?;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    buf.save(path)
        .map_err(|e| io::Error::other(format!("writing {}: {e}", path.display())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::events::{
        ClusterFocusInfo, ClusterInfo, Event, GenomeScannedInfo, MapUpdateInfo, PieceView,
    };
    use infinity_core::variant::{Variant, VariantKind};

    fn model_with_maps(members: u32) -> UiModel {
        let mut m = UiModel::default();
        for i in 0..members {
            m.apply(
                0.0,
                Event::GenomeScanned(GenomeScannedInfo {
                    id: i,
                    name: format!("g{i}"),
                    length: 1000,
                    contigs: 1,
                    n50: 1000,
                    gc: 0.5,
                    sketch_size: 10,
                }),
            );
        }
        let ids: Vec<u32> = (0..members).collect();
        m.apply(
            0.0,
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![ClusterInfo {
                    id: 0,
                    pivot: 0,
                    pivot_name: "g0".into(),
                    members: ids.clone(),
                    member_names: ids.iter().map(|i| format!("g{i}")).collect(),
                    focus: 0,
                    cluster_dist: 0.005,
                }],
            },
        );
        m.apply(
            0.0,
            Event::BackboneTree {
                newick: "(g0:0.001);".into(),
                pivots: vec!["g0".into()],
            },
        );
        m.apply(
            0.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(g0:1,g1:2);".into(),
                members: (0..members)
                    .map(|i| infinity_core::events::MemberNumbers {
                        id: i,
                        name: format!("g{i}"),
                        diff_pieces: u64::from(i),
                        snps: u64::from(i) * 2,
                        indels: 0,
                        uncovered_bp: 50,
                        novel_bp: 0,
                    })
                    .collect(),
            }),
        );
        for i in 1..members {
            m.apply(
                0.0,
                Event::MapUpdate(Box::new(MapUpdateInfo {
                    cluster: 0,
                    genome: i,
                    genome_name: format!("g{i}"),
                    pivot_name: "g0".into(),
                    pieces: vec![PieceView {
                        pivot_contig: 0,
                        pivot_contig_name: "c0".into(),
                        start: 100,
                        end: 140,
                        resolved: true,
                    }],
                    variants: vec![Variant {
                        kind: VariantKind::Snp,
                        pivot_contig: 0,
                        pivot_start: 10 + u64::from(i),
                        pivot_end: 11 + u64::from(i),
                        ref_allele: "A".into(),
                        alt_allele: "G".into(),
                        genome_contig: 0,
                        genome_start: 10,
                        genome_end: 11,
                        strand: '+',
                        flags: 0,
                    }],
                    uncovered: vec![(900, 950)],
                    pivot_contigs: vec!["c0".into()],
                    pivot_contig_lens: vec![1000],
                })),
            );
        }
        m
    }

    /// Render a writer into a string.
    fn render<F>(f: F) -> String
    where
        F: FnOnce(&mut RowCount<'_>) -> io::Result<()>,
    {
        let mut buf: Vec<u8> = Vec::new();
        let mut w = RowCount::new(&mut buf);
        f(&mut w).expect("writes");
        w.finish().expect("flushes");
        String::from_utf8(buf).expect("utf-8")
    }

    #[test]
    fn the_genomes_file_has_one_row_per_genome_and_section_nine_columns() {
        let m = model_with_maps(3);
        let s = render(|w| genomes_tsv(w, &m));
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 4, "header plus three genomes");
        assert!(lines[0].starts_with("id\tname\tlength\tcontigs\tn50\tgc\tcluster"));
        // The pivot: exact, no differences against itself.
        assert_eq!(
            lines[1],
            "0\tg0\t1000\t1\t1000\t0.500000\t0\ttrue\t2\t0\t0\t0\t50\t0"
        );
        // A member: its numbers came with the ClusterFocus event.
        assert_eq!(
            lines[2],
            "1\tg1\t1000\t1\t1000\t0.500000\t0\tfalse\t2\t1\t2\t0\t50\t0"
        );
    }

    #[test]
    fn a_genome_the_engine_never_reported_numbers_for_stays_na() {
        let mut m = model_with_maps(2);
        m.genomes[1].numbers = None;
        let s = render(|w| genomes_tsv(w, &m));
        let lines: Vec<&str> = s.lines().collect();
        assert!(
            lines[2].ends_with("\tNA\tNA\tNA\tNA\tNA"),
            "an unknown number must not be written as zero: {}",
            lines[2]
        );
    }

    #[test]
    fn the_clusters_file_reports_focus_and_whether_it_is_approximate() {
        let m = model_with_maps(3);
        let s = render(|w| clusters_tsv(w, &m));
        assert_eq!(s.lines().count(), 2);
        assert!(s.contains("cluster_id\tpivot\tpivot_name\tsize\tfocus\tapproximate"));
        assert!(s.contains("0\t0\tg0\t3\t2\tfalse"));
    }

    #[test]
    fn the_variants_file_converts_to_one_based_positions() {
        let m = model_with_maps(3);
        let s = render(|w| variants_tsv(w, &m, None, 0));
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 3, "header plus one variant per member");
        // genome 1 has a SNP at 0-based 11, which is 1-based 12.
        assert!(lines[1].contains("\tc0\t12\tA\tG\tsnp\t"), "{}", lines[1]);
        assert!(lines[1].starts_with("1\tg1\t"));
    }

    #[test]
    fn distances_are_derived_from_the_variant_lists_when_there_is_no_matrix() {
        let m = model_with_maps(4);
        let s = render(|w| distances_tsv(w, &m, None, 0));
        let lines: Vec<&str> = s.lines().collect();
        // 3 mapped members (1, 2, 3) -> 3 pairs.
        assert_eq!(lines.len(), 4, "{}", s);
        assert!(lines[0].starts_with("genome_a\tname_a\tgenome_b\tname_b\tsnps"));
        // Each member has a private SNP, so every pair is 2 SNPs apart.
        for l in &lines[1..] {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f[4], "2", "{l}");
            assert_eq!(f[6], "950", "comparable = 1000 - the shared 50 bp gap");
        }
    }

    #[test]
    fn a_cluster_below_focus_two_has_no_distances_to_export() {
        let mut m = model_with_maps(3);
        m.clusters[0].focus = 1;
        let s = render(|w| distances_tsv(w, &m, None, 0));
        assert_eq!(s.lines().count(), 1, "header only");
    }

    #[test]
    fn a_cluster_too_big_to_recompute_says_so_instead_of_hanging() {
        let m = model_with_maps(DERIVED_PAIR_CAP as u32 + 2);
        let mut buf: Vec<u8> = Vec::new();
        let mut w = RowCount::new(&mut buf);
        let err = distances_tsv(&mut w, &m, None, 0).expect_err("too many members");
        assert!(err.to_string().contains("more than"), "{err}");
    }

    #[test]
    fn the_tree_is_the_backbone_with_the_clusters_grafted_in() {
        let m = model_with_maps(3);
        let t = full_tree(&m, None).expect("a backbone exists");
        let nwk = to_newick(&t);
        // The backbone leaf "g0" became the cluster's own tree.
        assert!(nwk.contains("g1"), "{nwk}");
        assert_eq!(t.leaf_count(), 2);
    }

    #[test]
    fn the_tree_file_fails_honestly_when_there_is_no_backbone() {
        let m = UiModel::default();
        let mut buf: Vec<u8> = Vec::new();
        let mut w = RowCount::new(&mut buf);
        assert!(tree_newick(&mut w, &m, None).is_err());
    }

    #[test]
    fn tabs_and_newlines_in_a_name_cannot_break_the_layout() {
        assert_eq!(cell("a\tb\nc\r"), "a b c ");
        assert_eq!(cell("plain"), "plain");
    }

    #[test]
    fn writing_to_a_real_folder_produces_the_expected_files() {
        let dir = std::env::temp_dir().join(format!("infinity-gui-export-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        let m = model_with_maps(3);
        let report = write_all(&dir, &m, None);
        assert!(report.ok(), "{:?}", report.errors);
        assert!(dir.join("tree.nwk").is_file());
        assert!(dir.join("clusters.tsv").is_file());
        assert!(dir.join("genomes.tsv").is_file());
        assert!(dir.join("cluster_0.variants.tsv").is_file());
        assert!(dir.join("cluster_0.dist.tsv").is_file());
        // A cluster with no map data leaves no header-only file behind.
        assert!(!dir.join("cluster_1.variants.tsv").exists());
        assert!(!report.written.is_empty());
        // The window's own tree files come with "everything".
        assert!(dir.join("tree.nexus").is_file());
        assert!(dir.join("tree.svg").is_file());
        assert!(dir.join("itol_cluster_strip.txt").is_file());
        assert!(dir.join("itol_label_colors.txt").is_file());
        // With no run folder behind the model, the engine-written files are a
        // note, not an error.
        assert!(
            report.notes.iter().any(|n| n.contains("results.xlsx")),
            "{:?}",
            report.notes
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn engine_files_are_copied_from_the_run_folder_or_reported_missing() {
        let base = std::env::temp_dir().join(format!("infinity-gui-copy-{}", std::process::id()));
        std::fs::remove_dir_all(&base).ok();
        let run = base.join("run");
        let out = base.join("out");
        std::fs::create_dir_all(&run).unwrap();
        std::fs::create_dir_all(&out).unwrap();
        std::fs::write(run.join("results.xlsx"), b"PK fake").unwrap();
        std::fs::write(run.join("cluster_0.snp_matrix.tsv"), "\tg0\ng0\t0\n").unwrap();
        std::fs::write(run.join("cluster_3.snp_matrix.tsv"), "x").unwrap();
        std::fs::write(run.join("snp_clusters.tsv"), "x").unwrap();
        std::fs::write(run.join("cluster_0.dist.tsv"), "not copied").unwrap();
        let m = model_with_maps(2);
        let ctx = ExportCtx {
            run_dir: Some(&run),
            ..ExportCtx::new(&m, None)
        };
        let mut notes = Vec::new();
        let got = write_item(Item::Excel, &out, &ctx, &mut notes).expect("copied");
        assert_eq!(got, vec![out.join("results.xlsx")]);
        assert_eq!(std::fs::read(out.join("results.xlsx")).unwrap(), b"PK fake");
        let got = write_item(Item::SnpMatrices, &out, &ctx, &mut notes).expect("copied");
        assert_eq!(got.len(), 3, "two matrices and the SNP clusters: {got:?}");
        assert!(!out.join("cluster_0.dist.tsv").exists());
        // An older run has no VCF: a clear message naming what is missing.
        let err = write_item(Item::Vcf, &out, &ctx, &mut notes).expect_err("no vcf");
        assert!(
            err.contains("cluster_*.vcf") && err.contains("not found"),
            "{err}"
        );
        // No run folder at all is said differently.
        let bare = ExportCtx::new(&m, None);
        let err = write_item(Item::Excel, &out, &bare, &mut notes).expect_err("no folder");
        assert!(err.contains("no run folder"), "{err}");
        // Exporting into the run folder itself does not copy a file onto itself.
        let ctx_same = ExportCtx {
            run_dir: Some(&run),
            ..ExportCtx::new(&m, None)
        };
        let got = write_item(Item::Excel, &run, &ctx_same, &mut notes).expect("same folder");
        assert_eq!(got.len(), 1);
        assert_eq!(std::fs::read(run.join("results.xlsx")).unwrap(), b"PK fake");
        std::fs::remove_dir_all(&base).ok();
    }

    #[test]
    fn globs_match_one_star() {
        assert!(glob_match("cluster_*.vcf", "cluster_12.vcf"));
        assert!(!glob_match("cluster_*.vcf", "cluster_12.vcf.gz"));
        assert!(!glob_match("cluster_*.vcf", "cluster.vcf"));
        assert!(glob_match("results.xlsx", "results.xlsx"));
        assert!(!glob_match("results.xlsx", "old_results.xlsx"));
    }

    #[test]
    fn every_menu_item_has_words_and_only_engine_items_copy() {
        for i in Item::menu() {
            assert!(!i.label().is_empty() && !i.short_label().is_empty());
            let copies = !i.run_files().is_empty();
            assert_eq!(
                copies,
                matches!(
                    i,
                    Item::Excel | Item::SnpMatrices | Item::Vcf | Item::SnpAlignment
                )
            );
        }
        assert!(!Item::files().contains(&Item::Png));
        assert!(Item::menu().contains(&Item::Png));
    }

    #[test]
    fn exporting_nothing_reports_nothing_rather_than_an_empty_folder() {
        let dir =
            std::env::temp_dir().join(format!("infinity-gui-export-empty-{}", std::process::id()));
        let m = UiModel::default();
        let report = write_all(&dir, &m, None);
        assert!(!report.ok());
        assert!(report.summary().contains("nothing to export"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_png_round_trips_through_the_encoder() {
        let dir = std::env::temp_dir().join(format!("infinity-gui-png-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("tree.png");
        let img = egui::ColorImage::filled([4, 3], egui::Color32::from_rgb(10, 20, 30));
        write_png(&path, &img, None).expect("writes a png");
        // The PNG signature, then a real decode.
        let raw = std::fs::read(&path).unwrap();
        assert_eq!(&raw[..8], b"\x89PNG\r\n\x1a\n");
        let back = image::open(&path).expect("reads it back").to_rgba8();
        assert_eq!(back.dimensions(), (4, 3));
        assert_eq!(back.get_pixel(0, 0).0, [10, 20, 30, 255]);
        // A crop of the same picture keeps only the asked-for texels.
        let cropped = egui::Rect::from_min_size(egui::pos2(0.0, 0.0), egui::vec2(2.0, 1.0));
        write_png(&path, &img, Some((cropped, 1.0))).expect("writes a cropped png");
        let back = image::open(&path).expect("reads it back").to_rgba8();
        assert_eq!(back.dimensions(), (2, 1));
        // An empty capture is refused rather than written as a broken file.
        let empty = egui::ColorImage::filled([0, 0], egui::Color32::BLACK);
        assert!(write_png(&path, &empty, None).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
