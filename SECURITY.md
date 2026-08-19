# Security Policy

The security and integrity of **WinCalc-rs** are fundamental to the project's mission of delivering a memory-safe, reliable calculation engine.

---

### Core Security Guarantees

* **Strict Memory Safety:** The codebase enforces `#![forbid(unsafe_code)]` at the root and across all submodules. Raw pointers, manual pointer arithmetic, and unverified memory operations are strictly prohibited.
* **Arithmetic Correctness:** Arbitrary-precision calculations, base conversions, and rational operations must handle potential overflows, underflows, and edge cases safely without causing panics or undefined behavior.
* **Deterministic Execution:** Memory allocation and state mutations managed via `Rc` and `RefCell` must preserve invariants and prevent borrow panics in single-threaded environments.

---

### Scope of Security Vulnerabilities

We consider the following to be actionable security issues:

* Any mechanism, compiler exploit, or dependency artifact that bypasses `#![forbid(unsafe_code)]` or triggers undefined behavior.
* Uncontrolled panics or resource exhaustion denial-of-service (DoS) triggered by malformed mathematical expressions or edge-case numerical inputs.
* Memory leaks, cyclic references, or state corruption in the calculation engine lifecycle.
* Flaws in foreign function interface (FFI) boundaries that compromise host runtime safety.

---

### Reporting a Vulnerability

If you discover a security vulnerability or memory safety flaw, please adhere to responsible disclosure practices:

1. **Private Disclosure:** Do not open a public GitHub issue, pull request, or discussion thread.
2. **Submission Method:** Submit your findings privately through the **GitHub Security Advisory** tab on this repository.
3. **Report Contents:**
* A clear and concise description of the issue.
* A minimal, reproducible example (PoC) or input payload triggering the behavior.
* Details on the affected component (e.g., Ratpack, CEngine, FFI bindings).
* Any potential impact analysis or proposed mitigations.



---

### Handling and Response Process

* **Acknowledgment:** Maintainers will acknowledge receipt of the vulnerability report within **48 hours**.
* **Assessment & Fix:** The issue will be triaged, validated, and resolved in a private advisory workspace.
* **Disclosure:** Once a patch is finalized and merged, a formal security advisory will be published, giving full credit to the researcher (unless anonymity is requested).
