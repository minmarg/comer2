# install perl scripts
install(FILES 
	makepro_and_run_comer.pl 
	DESTINATION ${LOCALBINDIR} 
	PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_READ WORLD_READ
)

