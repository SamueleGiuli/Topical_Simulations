GFORT := gfortran

Task: ./vmc.x ./vmc_sd_read.x

vmc.x: ./VMC.f90
	$(GFORT) VMC.f90 -o vmc.x


vmc_sd_read.x: ./VMC_SD_read.f90
	$(GFORT) VMC_SD_read.f90 -o vmc_sd_read.x

clean:
	@rm *.x
	@rm modello.mod
