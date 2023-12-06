# la-por implementation in Rust

This is an implementation of the Init and Audit methods of https://github.com/dsroche/la-por in Rust. We are omitting the construction and usage of a Merkle tree, since we do not implement the Update method.

This project is aimed to be used as a library (WIP). Therefore, `main.rs` consists of the implementation of the algorithm and an example program how to use it to initialize the states for the server/client and to run a verification.

# Prerequisites

You need Rust in order to run the project. You can install it via [rustup](https://rustup.rs/).

# Installation

```sh
git clone https://github.com/luchev/rust-lapor
cd rust-lapor
cargo install --path .
```

# Running tests

```sh
cargo test
```
