#MAIN=main.for
MAIN=maim1.for

SRC=discret.for dop_prots.for les1.for $(MAIN) pechat.for proizv_les.for

.PHONY: all test build

all: test

test: build
	./les

build: les

les: $(SRC)
	pwd
	gfortran -ffixed-form -std=legacy $(SRC) -o les
