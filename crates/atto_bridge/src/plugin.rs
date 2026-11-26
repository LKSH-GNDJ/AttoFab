// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A C ABI for an external solver, with a reference implementation (W-10).
//!
//! The point of the boundary is that a third party can replace one model without rebuilding the
//! simulator. That needs three things, and only two of them belong in this crate:
//!
//! * a **versioned, `repr(C)` request and response**, defined here;
//! * a **reference implementation** of the entry point, also here, so the ABI is exercised rather
//!   than merely declared;
//! * **loading** the shared object, which is file I/O and belongs in an application shell — purity
//!   (A-02) forbids it in a domain crate, and `dlopen` is not available in every target.
//!
//! The version check is the part that matters in practice: a plugin built against a different
//! struct layout must be refused, not misinterpreted.

/// The ABI version. Bump this whenever either struct changes.
pub const ABI_VERSION: u32 = 1;

/// Status codes a plugin may return.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
/// The outcome of a plugin call, as an integer the C side can read.
pub enum PluginStatus {
    /// The request was served.
    Ok = 0,
    /// The plugin was built against a different ABI.
    VersionMismatch = 1,
    /// The request was malformed.
    BadRequest = 2,
    /// The plugin could not complete the solve.
    Failed = 3,
}

impl PluginStatus {
    /// The status a raw integer names, or `None` if it is not one of ours.
    pub fn from_code(code: i32) -> Option<Self> {
        match code {
            0 => Some(PluginStatus::Ok),
            1 => Some(PluginStatus::VersionMismatch),
            2 => Some(PluginStatus::BadRequest),
            3 => Some(PluginStatus::Failed),
            _ => None,
        }
    }
}

/// What the simulator asks a plugin to compute.
///
/// `repr(C)` with fixed-width fields only, so the layout is the same in Rust and C.
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(C)]
/// What the simulator asks a plugin to compute.
pub struct PluginRequest {
    /// Must equal [`ABI_VERSION`].
    pub version: u32,
    /// Number of points in the input array.
    pub n: u32,
    /// Input values.
    pub input: *const f64,
}

/// What a plugin returns.
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(C)]
/// What a plugin returns.
pub struct PluginResponse {
    /// One of [`PluginStatus`] as an integer.
    pub status: i32,
    /// Number of values written to `output`.
    pub n_out: u32,
    /// Output values.
    pub output: *mut f64,
}

/// The entry point every plugin must export.
///
/// # Safety
///
/// The caller must ensure both pointers are valid for the duration of the call and that
/// `response.output` has room for `response.n_out` values.
pub type PluginEntry = unsafe extern "C" fn(*const PluginRequest, *mut PluginResponse) -> i32;

/// Check a request's version before anything else touches it.
///
/// # Safety
///
/// `req` must be non-null and point to a readable `PluginRequest`.
pub unsafe fn check_version(req: *const PluginRequest) -> PluginStatus {
    if req.is_null() {
        return PluginStatus::BadRequest;
    }
    if (*req).version != ABI_VERSION {
        return PluginStatus::VersionMismatch;
    }
    PluginStatus::Ok
}

/// A reference plugin: doubles every input value.
///
/// This exists so the ABI is exercised, not just declared. It is a plain `extern "C"` function, so
/// the same code could be built into a shared object and loaded by a shell — but the test here
/// calls it directly, which proves the signature and the struct layouts without `dlopen`.
///
/// # Safety
///
/// Same contract as [`PluginEntry`].
pub unsafe extern "C" fn reference_double(req: *const PluginRequest, resp: *mut PluginResponse) -> i32 {
    let status = check_version(req);
    if status != PluginStatus::Ok {
        return status as i32;
    }
    if resp.is_null() || (*req).input.is_null() || (*resp).output.is_null() {
        return PluginStatus::BadRequest as i32;
    }
    let n = (*req).n as usize;
    let room = (*resp).n_out as usize;
    if n > room {
        // the caller promised room for n_out values; refusing is better than overrunning
        return PluginStatus::BadRequest as i32;
    }
    let input = std::slice::from_raw_parts((*req).input, n);
    let output = std::slice::from_raw_parts_mut((*resp).output, n);
    for (o, i) in output.iter_mut().zip(input.iter()) {
        *o = i * 2.0;
    }
    (*resp).status = PluginStatus::Ok as i32;
    (*resp).n_out = n as u32;
    PluginStatus::Ok as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reference_plugin_runs_and_returns_a_structure() {
        // W-10's acceptance, without dlopen: the ABI is exercised through the real signature
        let input = [1.0f64, 2.5, -3.0];
        let mut output = [0.0f64; 3];
        let req = PluginRequest { version: ABI_VERSION, n: 3, input: input.as_ptr() };
        let mut resp = PluginResponse { status: -1, n_out: 3, output: output.as_mut_ptr() };
        let code = unsafe { reference_double(&req, &mut resp) };
        assert_eq!(code, PluginStatus::Ok as i32, "the call succeeded");
        assert_eq!(resp.status, PluginStatus::Ok as i32, "and said so in the response");
        assert_eq!(resp.n_out, 3, "three values written");
        assert_eq!(output, [2.0, 5.0, -6.0], "and they are the doubled inputs: {:?}", output);
    }

    #[test]
    fn a_plugin_built_against_another_abi_is_refused() {
        // the check that matters: a mismatched layout must not be misinterpreted
        let input = [1.0f64];
        let mut output = [0.0f64; 1];
        let req = PluginRequest { version: ABI_VERSION + 1, n: 1, input: input.as_ptr() };
        let mut resp = PluginResponse { status: -1, n_out: 1, output: output.as_mut_ptr() };
        let code = unsafe { reference_double(&req, &mut resp) };
        assert_eq!(code, PluginStatus::VersionMismatch as i32, "refused, not run");
        assert_eq!(output, [0.0], "and nothing was written");
    }

    #[test]
    fn a_null_pointer_is_refused_rather_than_dereferenced() {
        let code = unsafe { reference_double(std::ptr::null(), std::ptr::null_mut()) };
        assert_eq!(code, PluginStatus::BadRequest as i32);
    }

    #[test]
    fn a_response_without_room_is_refused() {
        let input = [1.0f64, 2.0, 3.0];
        let mut output = [0.0f64; 1];
        let req = PluginRequest { version: ABI_VERSION, n: 3, input: input.as_ptr() };
        let mut resp = PluginResponse { status: -1, n_out: 1, output: output.as_mut_ptr() };
        let code = unsafe { reference_double(&req, &mut resp) };
        assert_eq!(code, PluginStatus::BadRequest as i32, "three values into room for one");
    }

    #[test]
    fn the_structs_have_a_fixed_layout() {
        // repr(C) with fixed-width fields: a pointer and two u32s is 16 bytes on 64-bit
        assert_eq!(std::mem::size_of::<u32>(), 4);
        assert_eq!(std::mem::size_of::<PluginRequest>(), 16, "two u32 and a pointer");
        assert_eq!(std::mem::size_of::<PluginResponse>(), 16, "an i32, a u32 and a pointer");
        assert_eq!(std::mem::align_of::<PluginRequest>(), 8);
    }

    #[test]
    fn every_status_code_round_trips() {
        for s in [PluginStatus::Ok, PluginStatus::VersionMismatch, PluginStatus::BadRequest, PluginStatus::Failed] {
            assert_eq!(PluginStatus::from_code(s as i32), Some(s));
        }
        assert_eq!(PluginStatus::from_code(99), None, "an unknown code is not invented");
    }
}
