# Contributing to WinCalc-rs

First off, thanks for taking the time to contribute! 🎉

This project is a 100% Safe Rust translation of the Microsoft Windows Calculator core engine, built with strictly enforced memory safety.

---

### Ground Rules

* **Zero Unsafe Code:** `#![forbid(unsafe_code)]` is enforced globally. No `unsafe` blocks or raw pointer arithmetic are allowed under any circumstances.
* **Idiomatic Rust:** Prefer clean, modern Rust patterns (enums, pattern matching, safe abstractions, `Rc`/`RefCell` for state management) over literal line-by-line C++ translation artifacts.
* **Logic Parity:** The mathematical behavior and calculation accuracy must match the original Windows Calculator behavior.

---

### How to Contribute

1. **Fork** the repository.
2. **Clone** your fork locally:
```bash
git clone https://github.com/YOUR_USERNAME/WinCalc-rs.git

```


3. **Create a branch** for your feature or bugfix:
```bash
git checkout -b feature/amazing-feature

```


4. **Commit your changes:**
```bash
git commit -m "feat: add amazing feature"

```


5. **Push to the branch:**
```bash
git commit -m "feat: add amazing feature"

```


6. Open a **Pull Request**.

---

### Code Standards

* Run `cargo check` and `cargo test` locally before submitting your code.
* Ensure code formatting complies with standard `rustfmt` rules:
```bash
cargo fmt

```


* Run linter checks to catch potential issues:
```bash
cargo clippy

```



---

### License

By contributing, you agree that your contributions will be licensed under the MIT License.
