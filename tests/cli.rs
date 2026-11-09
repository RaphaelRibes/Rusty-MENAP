//! End-to-end checks of the binary: arguments, outputs, exit codes.

use std::path::{Path, PathBuf};
use std::process::Command;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("rusty-menap-cli-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn table(dir: &Path, body: &str) -> PathBuf {
    let path = dir.join("table.txt");
    std::fs::write(&path, body).unwrap();
    path
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rusty-menap"))
        .args(args)
        .output()
        .expect("spawn binary")
}

fn arg(path: &Path) -> &str {
    path.to_str().unwrap()
}

#[test]
fn fixed_threshold_writes_a_consistent_edge_list() {
    let dir = scratch("fixed");
    let input = table(
        &dir,
        "ID\tS1\tS2\tS3\tS4\n\
         OTU1\t1\t2\t3\t4\n\
         OTU2\t2\t4\t6\t8\n\
         OTU3\t4\t3\t2\t1\n\
         OTU4\t1\t5\t2\t8\n",
    );
    let out = dir.join("out");
    let res = run(&["--st", "0.9", "-i", arg(&input), "-o", arg(&out)]);
    assert!(
        res.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&res.stderr)
    );

    let summary = std::fs::read_to_string(out.join("summary.json")).unwrap();
    assert!(summary.contains("\"mode\": \"fixed\""));
    let edges: usize = std::fs::read_to_string(out.join("network.tsv"))
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();
    assert!(edges > 0);
    assert!(summary.contains(&format!("\"edges\": {edges}")));
}

#[test]
fn ragged_input_fails_with_the_offending_line() {
    let dir = scratch("ragged");
    let input = table(&dir, "ID\tA\tB\nOTU1\t1\t2\nOTU2\t3\t4\t5\n");
    let res = run(&["-i", arg(&input), "-o", arg(&dir.join("out"))]);
    assert_eq!(res.status.code(), Some(1));
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("line 3"), "stderr: {err}");
}

#[test]
fn unknown_option_exits_with_code_two() {
    let res = run(&["--nope"]);
    assert_eq!(res.status.code(), Some(2));
}

#[test]
fn empty_network_at_a_fixed_threshold_is_an_error() {
    let dir = scratch("empty");
    // Constant rows carry no correlation, so no pair can reach the threshold.
    let input = table(&dir, "ID\tA\tB\nOTU1\t5\t5\nOTU2\t7\t7\n");
    let res = run(&[
        "--st",
        "0.9",
        "-i",
        arg(&input),
        "-o",
        arg(&dir.join("out")),
    ]);
    assert!(!res.status.success());
    assert!(
        String::from_utf8_lossy(&res.stderr).contains("network is empty"),
        "stderr: {}",
        String::from_utf8_lossy(&res.stderr)
    );
}
