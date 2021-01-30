GFORT := ifort

Task: ./dmc.x ./dmc_read.x

dmc.x: ./DMC.f90
	$(GFORT) DMC.f90 -o dmc.x


dmc_read.x: ./DMC_read.f90
	$(GFORT) DMC_read.f90 -o dmc_read.x

clean:
	@rm *.x
	@rm modello.mod
