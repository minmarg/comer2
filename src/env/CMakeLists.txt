## install all configuration files
install(FILES
	comer.par
	comer.par.org
	cvs2s
	hdp1.par
	hdps
	hdps33
	hdps36
	hdps66
	ssss
	DESTINATION ${LOCALSTATEDIR}
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
)
## install options file
install(FILES
	options.txt
	DESTINATION ${LOCALSTATEDIR}
	PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
)

