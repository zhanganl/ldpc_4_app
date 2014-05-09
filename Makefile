MAKE = make -s

.PHONY: 
all: ldpc


.PHONY: ldpc
ldpc:
	@cd src; ${MAKE} cleanall; ${MAKE}
	@cd demos; ${MAKE} cleanall; ${MAKE}

.PHONY: clean
clean:
	@cd src; ${MAKE} clean
	@cd demos; ${MAKE} clean


.PHONY: cleanall
cleanall:
	@cd src; ${MAKE} cleanall
	@cd demos; ${MAKE} cleanall

