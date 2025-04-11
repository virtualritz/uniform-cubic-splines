fn main() {
    // FIXME: use autocfg to emit it
    println!("cargo:rustc-cfg=has_f16");
    println!("cargo:rustc-check-cfg=cfg(has_f16)");
    println!("cargo:rustc-cfg=has_f128");
    println!("cargo:rustc-check-cfg=cfg(has_f128)");
}
