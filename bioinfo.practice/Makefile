PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd  check clean

rd:
	Rscript -e 'library(methods); devtools::document()'


build:
	Rscript -e 'devtools::build()'

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check:
	Rscript -e 'devtools::check()'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/
