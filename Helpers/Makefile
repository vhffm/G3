all: ph
ph:
	f2py -c -m potential_helpers potential_helpers.f90
	f95 -o potential_helpers potential_helpers.f90
clean:
	rm potential_helpers.so potential_helpers
