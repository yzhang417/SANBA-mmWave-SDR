## Comments
This folder [Experimental_Validation/src/array_control] contains all the 
related Matlab functions and the interface functions (in C language) to 
control the phased array by UART. 

The control of the phased arrays is mainly used for calibration and data
collection.

The controlling interface via UART is available by a library in C language, 
which is used in the function set PART.c. 

The Matlab functions will call the C functions written in PART.c to 
indirectly control the phased arrays.