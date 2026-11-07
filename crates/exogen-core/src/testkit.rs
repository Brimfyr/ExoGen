// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! What the tests share: a whole-body run built once and read by every check that reads it,
//! each check under its own name, so that a failing check is known by it and every failing
//! check is reported, not only the first.

/// Runs each check, labelled, and fails with every failing check's label and message.
pub(crate) fn checks(list: &[(&str, &dyn Fn())]) {
    let failed: Vec<String> = list
        .iter()
        .filter_map(|&(label, check)| {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(check))
                .err()
                .map(|e| {
                    let message = e
                        .downcast_ref::<String>()
                        .cloned()
                        .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                        .unwrap_or_default();
                    format!("{label}: {message}")
                })
        })
        .collect();
    assert!(
        failed.is_empty(),
        "{} of {} checks failed:\n{}",
        failed.len(),
        list.len(),
        failed.join("\n")
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every failing check is reported by its label; the passing ones are not.
    #[test]
    fn every_failing_check_is_reported_by_its_name() {
        let result = std::panic::catch_unwind(|| {
            checks(&[
                ("passes", &|| {}),
                ("fails first", &|| panic!("one")),
                ("fails second", &|| panic!("two {}", 2)),
            ])
        });
        let message = result
            .err()
            .and_then(|e| e.downcast_ref::<String>().cloned())
            .unwrap_or_default();
        assert!(message.contains("2 of 3 checks failed"), "{message}");
        assert!(message.contains("fails first: one"), "{message}");
        assert!(message.contains("fails second: two 2"), "{message}");
        assert!(!message.contains("passes:"), "{message}");
    }
}
