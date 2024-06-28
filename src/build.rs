fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            println!("cargo:rustc-cfg=feature=\"avx2\"");
        }
    }
}
