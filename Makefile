CARGO ?= cargo

.PHONY: all build test clean

all: build

build:
	$(CARGO) build --release

test:
	$(CARGO) test --release

clean:
	$(CARGO) clean
