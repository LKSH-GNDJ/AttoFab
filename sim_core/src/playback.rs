//! Animation transport state machine (PRD §3).
//!
//! Playback position is measured in *segments* (one per stage); each segment
//! has `FRAMES_PER_STAGE + 1` keyframes. The transport never computes
//! geometry: it only says *which frame* to show, and the caller replays that
//! frame deterministically (see `session::FrameCache`). Pausing therefore
//! cannot desynchronise the animation from the simulation state.

use crate::pipeline::FRAMES_PER_STAGE;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlayState { Stopped, Playing, Paused }

#[derive(Debug, Clone)]
pub struct Playback {
    pub state: PlayState,
    /// 0 ..= segments.
    pub position: f64,
    pub segments: usize,
    pub seconds_per_segment: f64,
    pub speed: f64,
}

impl Default for Playback {
    fn default() -> Self { Self { state: PlayState::Stopped, position: 0.0, segments: 0, seconds_per_segment: 4.0, speed: 1.0 } }
}

impl Playback {
    pub fn new(segments: usize) -> Self { Self { segments, ..Default::default() } }

    /// Replace the timeline length; the position is clamped. Called when the recipe, preview or run changes.
    pub fn set_segments(&mut self, n: usize) {
        self.segments = n;
        self.position = self.position.clamp(0.0, n as f64);
        if n == 0 { self.state = PlayState::Stopped; }
    }

    pub fn is_playing(&self) -> bool { self.state == PlayState::Playing }
    pub fn at_end(&self) -> bool { self.position >= self.segments as f64 - 1e-12 }

    pub fn play(&mut self) {
        if self.segments == 0 { return; }
        if self.at_end() { self.position = 0.0; }
        self.state = PlayState::Playing;
    }
    pub fn pause(&mut self) { if self.state == PlayState::Playing { self.state = PlayState::Paused; } }
    pub fn resume(&mut self) { if self.state == PlayState::Paused { self.state = PlayState::Playing; } }
    pub fn toggle(&mut self) { if self.is_playing() { self.pause() } else if self.state == PlayState::Paused { self.resume() } else { self.play() } }
    pub fn reset(&mut self) { self.position = 0.0; self.state = PlayState::Stopped; }

    fn frame_step() -> f64 { 1.0 / FRAMES_PER_STAGE as f64 }

    /// Step one keyframe forward; pauses playback.
    pub fn step_forward(&mut self) {
        let (seg, k) = self.at_keyframe();
        let _ = seg;
        self.position = (((self.position / Self::frame_step()).round() + if k >= FRAMES_PER_STAGE && self.position < self.segments as f64 { 0.0 } else { 1.0 }) * Self::frame_step()).min(self.segments as f64);
        self.state = if self.segments > 0 { PlayState::Paused } else { PlayState::Stopped };
    }
    pub fn step_back(&mut self) {
        self.position = (((self.position / Self::frame_step()).round() - 1.0) * Self::frame_step()).max(0.0);
        self.state = if self.segments > 0 { PlayState::Paused } else { PlayState::Stopped };
    }
    pub fn seek(&mut self, position: f64) {
        if position.is_finite() { self.position = position.clamp(0.0, self.segments as f64); }
    }

    /// Advance by `dt` wall-clock seconds. Only moves while `Playing`.
    pub fn advance(&mut self, dt: f64) -> bool {
        if self.state != PlayState::Playing || !dt.is_finite() || dt <= 0.0 { return false; }
        self.position += dt * self.speed / self.seconds_per_segment.max(1e-3);
        if self.position >= self.segments as f64 { self.position = self.segments as f64; self.state = PlayState::Stopped; }
        true
    }

    /// (segment, keyframe) to display. At the very end this is the last frame of the last segment.
    pub fn at_keyframe(&self) -> (usize, usize) {
        if self.segments == 0 { return (0, 0); }
        let seg = (self.position.floor() as usize).min(self.segments - 1);
        let local = (self.position - seg as f64).clamp(0.0, 1.0);
        (seg, (local * FRAMES_PER_STAGE as f64).round() as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn anim01_pause_freezes_the_frame_and_resume_continues_from_it() {
        let mut p = Playback::new(3);
        p.play();
        assert!(p.advance(2.0));
        p.pause();
        let (pos, frame) = (p.position, p.at_keyframe());
        assert!(!p.advance(10.0), "paused playback must not advance");
        assert_eq!((p.position, p.at_keyframe()), (pos, frame));
        p.resume();
        assert!(p.advance(1.0));
        assert!(p.position > pos);
        assert!(p.at_keyframe() >= frame);
    }

    #[test]
    fn stepping_pauses_and_walks_one_keyframe_at_a_time() {
        let mut p = Playback::new(2);
        p.step_forward();
        assert_eq!(p.state, PlayState::Paused);
        assert_eq!(p.at_keyframe(), (0, 1));
        p.step_forward(); p.step_forward();
        assert_eq!(p.at_keyframe(), (0, 3));
        p.step_back();
        assert_eq!(p.at_keyframe(), (0, 2));
        for _ in 0..FRAMES_PER_STAGE + 3 { p.step_forward(); }
        assert!(p.at_keyframe().0 >= 1, "crosses into the next segment");
        p.reset();
        assert_eq!((p.position, p.state), (0.0, PlayState::Stopped));
    }

    #[test]
    fn playback_stops_at_the_end_and_replays_from_the_start() {
        let mut p = Playback::new(1);
        p.play();
        p.advance(100.0);
        assert!(p.at_end() && p.state == PlayState::Stopped);
        assert_eq!(p.at_keyframe(), (0, FRAMES_PER_STAGE), "last frame of the last segment");
        p.play();
        assert_eq!(p.position, 0.0);
    }

    #[test]
    fn empty_timeline_and_bad_inputs_are_safe() {
        let mut p = Playback::default();
        p.play(); assert_eq!(p.state, PlayState::Stopped);
        p.set_segments(2); p.seek(f64::NAN); assert_eq!(p.position, 0.0);
        p.seek(9.0); assert_eq!(p.position, 2.0);
        p.set_segments(1); assert_eq!(p.position, 1.0);
        assert!(!p.advance(f64::NAN));
    }
}
