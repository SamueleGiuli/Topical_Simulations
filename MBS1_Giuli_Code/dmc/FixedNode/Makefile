GFORT := gfortran

Task: ./dmc_fn.x ./dmc_fn_read.x

dmc_fn.x: ./DMC_FN.f90
	$(GFORT) DMC_FN.f90 -o dmc_fn.x


dmc_fn_read.x: ./DMC_FN_read.f90
	$(GFORT) DMC_FN_read.f90 -o dmc_fn_read.x

clean:
	@rm *.x
	@rm modello.mod
