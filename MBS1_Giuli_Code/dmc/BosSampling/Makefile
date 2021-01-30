GFORT := ifort

Task: ./dmc_bs.x ./dmc_bs_read.x

dmc_bs.x: ./DMC_BS.f90
	$(GFORT) DMC_BS.f90 -o dmc_bs.x


dmc_bs_read.x: ./DMC_BS_read.f90
	$(GFORT) DMC_BS_read.f90 -o dmc_bs_read.x

clean:
	@rm *.x
	@rm modello.mod
	@rm Process.dat
