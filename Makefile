.PHONY: tests precommit setup

# This will check for dead dependencies
check_dead_dep:
	cargo +nightly udeps --all-targets
