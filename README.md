# WinCalc-rs

<p align="center">
  <img src="https://img.shields.io/badge/Language-Safe%20Rust-orange.svg?style=for-the-badge&logo=rust" alt="Rust" />
  <img src="https://img.shields.io/badge/Unsafe-Forbidden-success.svg?style=for-the-badge" alt="No Unsafe" />
  <img src="https://img.shields.io/badge/Fork%20Of-Windows%20Calculator-0078D4.svg?style=for-the-badge&logo=windows" alt="Windows Calculator" />
  <img src="https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge" alt="License" />
</p>

A pure Safe Rust translation of the Microsoft Windows Calculator calculation engine.

---

### Features

* **Memory Safe:** Built with strictly enforced `#![forbid(unsafe_code)]`.
* **Complete Logic Parity:** Fully preserving the mathematical behaviors and features of the original Windows Calculator.
* **Drop-in Core:** Dedicated computational backend ready for desktop integrations.

---

### Getting Started

#### Prerequisites
* [Rust](https://www.rust-lang.org/) (latest stable toolchain)
* Cargo

#### Build
```bash
cargo build --release
```
#### Run Tests
```Bash
cargo test
```

### License
This project is licensed under the MIT License.
