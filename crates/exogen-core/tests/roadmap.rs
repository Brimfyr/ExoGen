// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! The README's Status stands where PLAN.md's roadmap says the work stands (the user,
//! 2026-11-19): every phase and sub-phase Done, Under way or Planned alike in both, and its
//! "Now" line naming the sub-phases under way. A gate marked passed in PLAN.md fails this check
//! until the README follows it, so the README's Status changes at each gate and at no other time.
//! The roadmap writes each phase (`### Phase 5: ...`) and sub-phase (`#### 5A: ...`) the same
//! way, its state on its own `**State:**` line under its heading.

use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum State {
    Done,
    UnderWay,
    Planned,
}

/// The roadmap as read: each phase's and sub-phase's state, and each sub-phase's phase.
struct Roadmap {
    states: BTreeMap<String, State>,
    parents: BTreeMap<String, String>,
}

fn read(name: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// A state by its words: Done, Under way or Planned, and nothing else.
fn state_of(text: &str) -> Option<State> {
    let text = text.trim_start();
    if text.starts_with("Done") {
        Some(State::Done)
    } else if text.starts_with("Under way") {
        Some(State::UnderWay)
    } else if text.starts_with("Planned") {
        Some(State::Planned)
    } else {
        None
    }
}

/// PLAN.md's roadmap: each `### Phase` and `####` heading's id and the state its `**State:**`
/// line gives, the first after it; a `####` heading a sub-phase of the phase above it. A heading
/// with no State line before the next heading is an error.
fn roadmap(plan: &str) -> Roadmap {
    let mut states = BTreeMap::new();
    let mut parents = BTreeMap::new();
    let mut phase: Option<String> = None;
    let mut open: Option<String> = None;
    for line in plan.lines() {
        let (id, sub) = if let Some(rest) = line.strip_prefix("### Phase ") {
            (rest.split(':').next(), false)
        } else if let Some(rest) = line.strip_prefix("#### ") {
            (rest.split(':').next(), true)
        } else {
            (None, false)
        };
        if line.starts_with("## ") || line.starts_with("### ") || line.starts_with("#### ") {
            assert!(open.is_none(), "{open:?} has no State line");
            if !line.starts_with("#### ") {
                phase = None;
            }
        }
        if let Some(id) = id.map(str::trim).filter(|i| !i.is_empty()) {
            if sub {
                let p = phase
                    .clone()
                    .unwrap_or_else(|| panic!("{id} under no phase"));
                parents.insert(id.to_string(), p);
            } else {
                phase = Some(id.to_string());
            }
            open = Some(id.to_string());
        } else if let Some(rest) = line.strip_prefix("**State:**")
            && let Some(id) = open.take()
        {
            let state = state_of(rest).unwrap_or_else(|| panic!("{id}'s State: {rest:?}"));
            states.insert(id, state);
        }
    }
    assert!(open.is_none(), "{open:?} has no State line");
    Roadmap { states, parents }
}

/// The README's Status: its table's phases and states, and its "Now" line.
fn readme_states(readme: &str) -> (BTreeMap<String, State>, String) {
    let status = readme
        .split("\n## Status")
        .nth(1)
        .expect("the README's Status section");
    let status = status.split("\n## ").next().unwrap_or("");
    let mut states = BTreeMap::new();
    let mut now = String::new();
    for line in status.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("**Now:**") {
            now = rest.to_string();
        }
        if !line.starts_with('|') {
            continue;
        }
        let cells: Vec<&str> = line.trim_matches('|').split('|').map(str::trim).collect();
        let Some(state) = cells.last().and_then(|c| state_of(c)) else {
            continue;
        };
        states.insert(cells[0].to_string(), state);
    }
    (states, now)
}

#[test]
fn the_readme_s_status_stands_where_the_plan_does() {
    let plan = roadmap(&read("PLAN.md"));
    let (readme, now) = readme_states(&read("README.md"));
    assert!(
        plan.states.len() >= 10,
        "PLAN.md's roadmap read as {:?}",
        plan.states
    );
    // A phase with sub-phases stands as they do: done when all are, planned when none has
    // begun, under way else.
    let mut subs: BTreeMap<&String, Vec<State>> = BTreeMap::new();
    for (sub, p) in &plan.parents {
        subs.entry(p).or_default().push(plan.states[sub]);
    }
    for (p, s) in &subs {
        let state = if s.iter().all(|&x| x == State::Done) {
            State::Done
        } else if s.iter().all(|&x| x == State::Planned) {
            State::Planned
        } else {
            State::UnderWay
        };
        assert_eq!(
            plan.states[*p], state,
            "phase {p} against its sub-phases {s:?}"
        );
    }
    assert_eq!(
        readme, plan.states,
        "the README's Status table (left) must state each phase as PLAN.md's roadmap does (right): update the README at the gate"
    );
    // The "Now" line names every sub-phase, or phase without sub-phases, that is under way,
    // and nothing that is not.
    let named: Vec<&str> = now
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|w| plan.states.contains_key(*w))
        .collect();
    for (id, state) in &plan.states {
        let leaf = !subs.contains_key(id);
        let listed = named.contains(&id.as_str());
        assert_eq!(
            listed,
            leaf && *state == State::UnderWay,
            "the README's Now line ({now:?}) and {id}, {state:?}"
        );
    }
}

/// A roadmap the check reads: each heading's State line, 1F a phase of its own, 2A to 2C
/// phase 2's sub-phases.
#[test]
fn the_roadmap_s_states_are_read_from_its_state_lines() {
    let plan = "## 10. Roadmap\n\n### Phase 1: one (L)\n\n**State:** Done.\n\n**Builds:**\n- a\n\n\
                ### Phase 1F: one again (L)\n\n**State:** Done 2026-01-02, its gate passed.\n\n\
                ### Phase 2: two (L)\n\n**State:** Under way since 2026-01-01.\n\n\
                #### 2A: first\n\n**State:** Done 2026-01-01, its gate passed.\n\n\
                #### 2B: second\n\n**State:** Under way since 2026-01-02.\n\n\
                #### 2C: third\n\n**State:** Planned.\n\n### Risks and mitigations\n\n| a | b |\n\n\
                ## 11. Decisions\n";
    let r = roadmap(plan);
    let expected: BTreeMap<String, State> = [
        ("1", State::Done),
        ("1F", State::Done),
        ("2", State::UnderWay),
        ("2A", State::Done),
        ("2B", State::UnderWay),
        ("2C", State::Planned),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_string(), v))
    .collect();
    assert_eq!(r.states, expected);
    let parents: Vec<(&str, &str)> = r
        .parents
        .iter()
        .map(|(a, b)| (a.as_str(), b.as_str()))
        .collect();
    assert_eq!(parents, [("2A", "2"), ("2B", "2"), ("2C", "2")]);
    let readme = "# X\n\n## Status\n\n**Now:** 2B, second.\n\n| Phase | What | State |\n|---|---|---|\n\
                  | 1 | one | Done |\n| 2B | two | Under way |\n\n## Outputs\n\n| 9 | x | Done |\n";
    let (states, now) = readme_states(readme);
    assert_eq!(states.len(), 2, "{states:?}");
    assert_eq!(states["2B"], State::UnderWay);
    assert!(now.contains("2B"));
}

/// A heading without its State line is an error, not a phase taken as planned.
#[test]
#[should_panic(expected = "has no State line")]
fn a_roadmap_entry_without_its_state_line_is_refused() {
    roadmap("### Phase 1: one (L)\n\n- an item\n\n### Phase 2: two (L)\n\n**State:** Planned.\n");
}
