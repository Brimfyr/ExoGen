// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Brimfyr

//! Background worker: owns the stage runner and its caches, runs pipelines, previews and
//! exports off the UI thread, and reports through a channel.

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};

use exogen_core::export::{ExportSettings, export};
use exogen_core::params::ParamGraph;
use exogen_core::stages::{
    FinishOutput, LoopOutput, NoProgress, PipelineResult, Progress, ProgressSink, Runner, finish,
};

pub enum Job {
    Run(ParamGraph),
    Preview {
        graph: ParamGraph,
        loop_output: Arc<LoopOutput>,
        index: usize,
        face: u32,
    },
    Export {
        graph: ParamGraph,
        finish: Arc<FinishOutput>,
        dir: PathBuf,
        width: usize,
        name: String,
    },
}

pub enum Msg {
    Progress(Progress),
    RunDone(Arc<PipelineResult>),
    PreviewDone(usize, Arc<FinishOutput>),
    ExportDone(String),
    Error(String),
}

struct ChannelSink {
    tx: Sender<Msg>,
    cancel: Arc<AtomicBool>,
}

impl ProgressSink for ChannelSink {
    fn report(&self, p: Progress) {
        let _ = self.tx.send(Msg::Progress(p));
    }
    fn cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }
}

pub struct Worker {
    jobs: Sender<Job>,
    msgs: Receiver<Msg>,
    cancel: Arc<AtomicBool>,
    pub busy: bool,
}

impl Worker {
    pub fn start(repaint: impl Fn() + Send + 'static) -> Self {
        let (job_tx, job_rx) = channel::<Job>();
        let (msg_tx, msg_rx) = channel::<Msg>();
        let cancel = Arc::new(AtomicBool::new(false));
        let cancel_thread = cancel.clone();
        std::thread::Builder::new()
            .name("exogen-worker".into())
            .spawn(move || {
                let mut runner = Runner::new();
                for job in job_rx {
                    cancel_thread.store(false, Ordering::Relaxed);
                    let sink = ChannelSink {
                        tx: msg_tx.clone(),
                        cancel: cancel_thread.clone(),
                    };
                    let result = match job {
                        Job::Run(graph) => match runner.run(&graph, &sink) {
                            Ok(r) => Msg::RunDone(Arc::new(r)),
                            Err(e) => Msg::Error(e.to_string()),
                        },
                        Job::Preview {
                            graph,
                            loop_output,
                            index,
                            face,
                        } => match finish(&graph, &loop_output, face, index, &NoProgress) {
                            Ok(f) => Msg::PreviewDone(index, Arc::new(f)),
                            Err(e) => Msg::Error(e.to_string()),
                        },
                        Job::Export {
                            graph,
                            finish,
                            dir,
                            width,
                            name,
                        } => {
                            let settings = ExportSettings {
                                width,
                                dir: dir.clone(),
                                project_name: name,
                            };
                            match export(&graph, &finish, &settings) {
                                Ok(s) => Msg::ExportDone(format!(
                                    "wrote {} maps to {}",
                                    s.maps.len(),
                                    dir.display()
                                )),
                                Err(e) => Msg::Error(e.to_string()),
                            }
                        }
                    };
                    let _ = msg_tx.send(result);
                    repaint();
                }
            })
            .expect("worker thread");
        Self {
            jobs: job_tx,
            msgs: msg_rx,
            cancel,
            busy: false,
        }
    }

    pub fn submit(&mut self, job: Job) {
        self.busy = true;
        let _ = self.jobs.send(job);
    }

    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    pub fn poll(&mut self) -> Vec<Msg> {
        let mut out = Vec::new();
        while let Ok(m) = self.msgs.try_recv() {
            if !matches!(m, Msg::Progress(_)) {
                self.busy = false;
            }
            out.push(m);
        }
        out
    }
}
