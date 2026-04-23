CARGO ?= cargo
PREFIX ?= /usr/local
BINARY ?= fmsi-rust
TARGET_DIR ?= target

RUSTFLAGS_SUPER := -C target-cpu=native -C lto=fat -C codegen-units=1 -C panic=abort

.PHONY: all release superrelease install uninstall test clean fmt clippy

all: release

release:
	$(CARGO) build --release

superrelease:
	RUSTFLAGS="$(RUSTFLAGS_SUPER)" $(CARGO) build --release

install: release
	install -d $(PREFIX)/bin
	install -m 0755 $(TARGET_DIR)/release/$(BINARY) $(PREFIX)/bin/$(BINARY)

uninstall:
	rm -f $(PREFIX)/bin/$(BINARY)

test:
	$(CARGO) test --release

fmt:
	$(CARGO) fmt

clippy:
	$(CARGO) clippy --release -- -D warnings

clean:
	$(CARGO) clean
