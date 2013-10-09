#MAIN=main.for
MAIN=maim1.for

SRC=discret.for dop_prots.for les1.for $(MAIN) pechat.for proizv_les.for

.PHONY: all test build wrap

all: test

test: build
	./les

build: les

les: $(SRC)
	gfortran -ffixed-form -std=legacy -g $(SRC) -o les

wrap: $(SRC)
	f2py -c --f77flags="-ffixed-form -std=legacy" $(SRC)
