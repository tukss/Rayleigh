Module Input
	Use Grid, Only : grid_namelist
	Implicit None
Contains

	Subroutine Main_Input()
		Implicit None
		Open(unit=20, file="main_input", status="old", position="rewind")
		Read(unit=20, nml=grid_namelist)
		Close(20)
	End Subroutine Main_Input


End Module Input
