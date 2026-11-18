//! Remembered folders and parameters, kept in a small JSON file under the user's
//! configuration directory.
//!
//! The program deliberately has no native file dialog: the GTK-based crates
//! that provide one need system `-dev` packages at build time, which section 12
//! rules out for the Windows build. Paths are typed or dropped onto the window
//! instead, and then remembered here, so the second run of the day is one click
//! long.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use infinity_core::params::Parameters;

use crate::settings::ResolutionSettings;

/// Folder name under the platform's configuration directory.
pub const APP_DIR: &str = "infinity";
/// File name inside it.
pub const FILE: &str = "gui_prefs.json";
/// How many recent paths to remember per list.
pub const MAX_RECENT: usize = 12;

/// Everything the interface remembers between sessions.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct Prefs {
    /// Genome files or folders, most recent first.
    pub recent_inputs: Vec<String>,
    /// Output folders, most recent first.
    pub recent_outputs: Vec<String>,
    /// Last folder the export menu wrote to.
    pub export_dir: String,
    /// Last parameter set, so a re-run does not mean re-typing.
    pub params: Parameters,
    /// The Resolution and SNP-rule settings (see `settings`).
    pub resolution: ResolutionSettings,
    /// "dark", "light" or "system". Empty means dark, the default look.
    pub theme: String,
}

/// The platform's configuration directory, resolved without a `dirs` crate so
/// that the build stays free of system dependencies.
pub fn config_dir() -> Option<PathBuf> {
    #[cfg(target_os = "windows")]
    {
        std::env::var_os("APPDATA").map(PathBuf::from)
    }
    #[cfg(not(target_os = "windows"))]
    {
        if let Some(xdg) = std::env::var_os("XDG_CONFIG_HOME") {
            let p = PathBuf::from(xdg);
            if p.is_absolute() {
                return Some(p);
            }
        }
        std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".config"))
    }
}

/// Where the preferences file lives, or `None` when there is no home directory
/// to put it in (which must not stop the program from running).
pub fn prefs_path() -> Option<PathBuf> {
    config_dir().map(|d| d.join(APP_DIR).join(FILE))
}

/// Read the preferences, falling back to the defaults. A file written by a
/// future version, or one the user edited badly, must not stop the program.
pub fn load() -> Prefs {
    let Some(path) = prefs_path() else {
        return Prefs::default();
    };
    load_from(&path)
}

/// Read preferences from an explicit path. Split out from [`load`] so the
/// fallback behaviour can be tested.
pub fn load_from(path: &Path) -> Prefs {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Prefs::default();
    };
    if let Ok(p) = serde_json::from_str::<Prefs>(&text) {
        return p;
    }
    // A file written by a version whose `Parameters` had different fields cannot
    // be read whole. The remembered paths are still valid and worth keeping, so
    // drop just the parameter set and fall back to the defaults for it, rather
    // than throwing the whole file away.
    if let Ok(mut v) = serde_json::from_str::<serde_json::Value>(&text) {
        if let Some(obj) = v.as_object_mut() {
            obj.remove("params");
            if let Ok(p) = serde_json::from_value::<Prefs>(v.clone()) {
                return p;
            }
            // The same goes for the resolution settings, which mirror engine
            // parameters that are still changing shape.
            if let Some(obj) = v.as_object_mut() {
                obj.remove("resolution");
            }
            if let Ok(p) = serde_json::from_value::<Prefs>(v) {
                return p;
            }
        }
    }
    Prefs::default()
}

/// Write the preferences. Failure is not fatal and is reported by the caller.
pub fn save(prefs: &Prefs) -> std::io::Result<()> {
    let Some(path) = prefs_path() else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "no configuration directory",
        ));
    };
    save_to(prefs, &path)
}

/// Write preferences to an explicit path, creating the folder if needed.
pub fn save_to(prefs: &Prefs, path: &Path) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let json = serde_json::to_string_pretty(prefs)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    std::fs::write(path, json)
}

/// Put `item` at the front of a recent list, dropping an older copy and the
/// tail beyond [`MAX_RECENT`].
pub fn push_recent(list: &mut Vec<String>, item: &str, max: usize) {
    let item = item.trim();
    if item.is_empty() {
        return;
    }
    list.retain(|p| p != item);
    list.insert(0, item.to_string());
    list.truncate(max);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!("infinity-gui-prefs-{}-{name}", std::process::id()))
    }

    #[test]
    fn a_recent_list_moves_a_repeated_entry_to_the_front() {
        let mut v = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        push_recent(&mut v, "c", MAX_RECENT);
        assert_eq!(v, vec!["c", "a", "b"]);
        push_recent(&mut v, "d", MAX_RECENT);
        assert_eq!(v, vec!["d", "c", "a", "b"]);
    }

    #[test]
    fn a_recent_list_is_capped_and_ignores_blanks() {
        let mut v: Vec<String> = Vec::new();
        for i in 0..(MAX_RECENT + 5) {
            push_recent(&mut v, &format!("p{i}"), MAX_RECENT);
        }
        assert_eq!(v.len(), MAX_RECENT);
        assert_eq!(v[0], format!("p{}", MAX_RECENT + 4));
        push_recent(&mut v, "   ", MAX_RECENT);
        push_recent(&mut v, "", MAX_RECENT);
        assert_eq!(v.len(), MAX_RECENT);
    }

    #[test]
    fn preferences_round_trip_through_disk() {
        let path = tmp("roundtrip");
        let mut p = Prefs::default();
        p.recent_inputs.push("/data/genomes".to_string());
        p.recent_outputs.push("/data/out".to_string());
        p.export_dir = "/data/export".to_string();
        p.params.k = 21;
        p.params.cluster_dist = 0.01;
        p.resolution.snp_thresholds = vec![3, 12];
        p.resolution.snp_filter.enabled = false;
        p.theme = "light".to_string();
        save_to(&p, &path).expect("writes");
        let back = load_from(&path);
        assert_eq!(back, p);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn a_missing_or_broken_file_gives_the_defaults() {
        let missing = tmp("does-not-exist");
        assert_eq!(load_from(&missing), Prefs::default());

        let broken = tmp("broken");
        std::fs::write(&broken, "{ not json").unwrap();
        assert_eq!(load_from(&broken), Prefs::default());

        // A file from a version that knew fewer fields still loads.
        let partial = tmp("partial");
        std::fs::write(&partial, r#"{"recent_inputs":["/x"]}"#).unwrap();
        let p = load_from(&partial);
        assert_eq!(p.recent_inputs, vec!["/x".to_string()]);
        assert_eq!(p.params, Parameters::default());

        // A parameter set from a different schema cannot be read; the remembered
        // paths survive and only the parameters fall back.
        let foreign = tmp("foreign-params");
        std::fs::write(
            &foreign,
            r#"{"recent_inputs":["/data/x"],"params":{"k":30,"nope":1}}"#,
        )
        .unwrap();
        let p = load_from(&foreign);
        assert_eq!(p.recent_inputs, vec!["/data/x".to_string()]);
        assert_eq!(p.params, Parameters::default());
        p.params.validate().expect("the fallback must be valid");

        // Resolution settings of an unknown shape fall back the same way.
        let foreign_res = tmp("foreign-resolution");
        std::fs::write(
            &foreign_res,
            r#"{"recent_inputs":["/data/y"],"resolution":{"resolution":"warp-speed"}}"#,
        )
        .unwrap();
        let p = load_from(&foreign_res);
        assert_eq!(p.recent_inputs, vec!["/data/y".to_string()]);
        assert_eq!(p.resolution, ResolutionSettings::default());

        for f in [broken, partial, foreign, foreign_res] {
            std::fs::remove_file(f).ok();
        }
    }

    #[test]
    fn the_config_dir_is_absolute_or_absent() {
        if let Some(d) = config_dir() {
            assert!(d.is_absolute(), "{d:?} is not absolute");
            assert_eq!(prefs_path().unwrap().file_name().unwrap(), FILE);
        }
    }
}
