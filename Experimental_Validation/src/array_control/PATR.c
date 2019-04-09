#include "PATR.h"
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>    /* File Control Definitions           */
#include <termios.h>  /* POSIX Terminal Control Definitions */
#include <unistd.h>   /* UNIX Standard Definitions          */ 
#include <errno.h>    /* ERROR Number Definitions           */
#include <unistd.h>
#include <string.h>
#include <time.h>

//Define a global variable for file descriptor
int PATX = -1; 
int PARX = -1;
int PortOpenTx = -1;
int PortOpenRx = -1;
int operation_hold_time = 800;
int long_operation_hold_time = 1000000;


/*---------- Open the serial port using termios structure --------- */
int OpenSerialPortRx(const char *portName)
{
	/*Check whether the port is open already*/
	if(PARX == -1)
    {
        /* "PARX" - file descriptor pointing to the opened serial port */
        PARX = open(portName, O_RDWR | O_NOCTTY | O_NDELAY); 
        PortOpenRx = 1;
        /* O_RDWR Read/Write access to serial port           */
        /* O_NOCTTY - No terminal will control the process   */
        /* O_NDELAY -Non Blocking Mode,Does not care about-  */
        /* -the status of DCD line, Open() returns immediatly */				
    }    
    
    /* Error Checking */
    if(PARX == -1)
    {
    	mexPrintf("\nError In Opening Rx Phased Array!\n");
    	return -1;
    }
    else
    {
        mexPrintf("\nRx Phased Array Is Successfully Opened!\n");
    }
    	
	/*---------- Setting the Attributes of the serial port using termios structure --------- */
	struct termios SerialPortSettings;	/* Create the structure                          */
	tcgetattr(PARX, &SerialPortSettings);	/* Get the current attributes of the Serial port */
	cfsetispeed(&SerialPortSettings,B115200); /* Set Read  Speed as 115200                       */
	cfsetospeed(&SerialPortSettings,B115200); /* Set Write Speed as 115200                       */
	SerialPortSettings.c_cflag &= ~PARENB;   /* Disables the Parity Enable bit(PARENB),So No Parity   */
	SerialPortSettings.c_cflag &= ~CSTOPB;   /* CSTOPB = 2 Stop bits,here it is cleared so 1 Stop bit */
	SerialPortSettings.c_cflag &= ~CSIZE;	 /* Clears the mask for setting the data size             */
	SerialPortSettings.c_cflag |=  CS8;      /* Set the data bits = 8                                 */
	SerialPortSettings.c_cflag &= ~CRTSCTS;       /* No Hardware flow Control                         */
	SerialPortSettings.c_cflag |= CREAD | CLOCAL; /* Enable receiver,Ignore Modem Control lines       */ 
	SerialPortSettings.c_iflag &= ~(IXON | IXOFF | IXANY);          /* Disable XON/XOFF flow control both i/p and o/p */
	SerialPortSettings.c_iflag &= ~(ICANON | ECHO | ECHOE | ISIG);  /* Non Cannonical mode                            */
	SerialPortSettings.c_oflag &= ~OPOST;/*No Output Processing*/
	//SerialPortSettings.c_cflag &= ~(PARENB | PARODD);    /* Disables the Parity Enable bit(PARENB),So No Parity   */
    //SerialPortSettings.c_iflag &= ~IGNBRK;          // disable break processing
    SerialPortSettings.c_lflag = 0;                 
    SerialPortSettings.c_oflag = 0;                // no remapping, no delays
	//SerialPortSettings.c_cc[VMIN] = 10; /* Read at least 10 characters */
	//SerialPortSettings.c_cc[VTIME] = 0; /* Wait indefinetly   */
    
	if((tcsetattr(PARX,TCSANOW,&SerialPortSettings)) != 0) /* Set the attributes to the termios structure*/
	{
		mexPrintf("\Error In Setting Attributes For Rx Phased Array.\n");
		return -1;
	}    
    return 0;
}	  


/*---------- Close the serial port using termios structure --------- */
int CloseSerialPortRx(void)
{
    if(PARX == -1)
    {
        mexPrintf("\nError: No Opening Serial Port!\n");
    }
    else
    {
        mexPrintf("\nRx Phased Array Is Successfully Closed!\n");
        close(PARX);
        PARX = -1;
    }
    return 0;
}


/*------------------------------- Write data to serial port -----------------------------*/
int WriteCommandRx(const char *input_command, const char *portName)
{   
	int len_Std_Command = strlen(input_command)+2;
    char *Std_Command = (char*)malloc(len_Std_Command*sizeof(unsigned char));
    Std_Command[0] = '\n';
    Std_Command[len_Std_Command-1] = '\r';
        
    /*Copy command from input_command*/
    int i;
    for (i=1; i<=len_Std_Command-2; i++)
    {
        Std_Command[i] = input_command[i-1];
    }
    
    /*Check whether the serial port is opened*/
    int reopen;
	int bytes_written  = 0;  
    if (PortOpenRx != 1)
    {
        reopen = OpenSerialPortRx(portName);
    }
    
    /*Clean the read buffer before writting*/    
    tcflush(PARX, TCIFLUSH);                  
   
    /*Write Command*/
    bytes_written = write(PARX, Std_Command, len_Std_Command);
    //mexPrintf("Bytes written: %d", bytes_written);
    /*char input_command_4letters[4];
    for (i=0; i<=3; i++)
    {
        input_command_4letters[i] = input_command[i];
    }
    if (strstr(input_command_4letters, "rxon") != NULL)
    {
        //mexPrintf("\nSleeping for 1 milliseconds for turining on and off antennas: %s\n", input_command);
        usleep(operation_hold_time);
    }
    else
    {
        //mexPrintf("\nSleeping for 1 milliseconds for changing the phase shifter.\n");
        usleep(operation_hold_time);
    }
    usleep(operation_hold_time*2);*/
	return 0;	
}


/*------------------------------- Read data from serial port -----------------------------*/
int ReadCommandRx(void)
{   
	char read_buffer[300];
	int bytes_read = 0;              
	bytes_read = read(PARX, &read_buffer, 300);
    return bytes_read;
}


/*------------------------------- Read data from serial port -----------------------------*/
int ReadCommandRxPrint(void)
{   
	char read_buffer[300];
	int bytes_read = 0;              
	bytes_read = read(PARX, &read_buffer, 300);
    /*Print the number of read bytes and contents*/	
    /*mexPrintf("\nBytes Read: %d\n",bytes_read);*/
	int i = 0;
	for( i = 0; i < bytes_read; i++)	 		
    {
    	mexPrintf("%c",read_buffer[i]);
    } 	
    mexPrintf("\n");
    return bytes_read;
}


/*---------- Open the serial port using termios structure --------- */
int OpenSerialPortTx(const char *portName)
{
	/*Check whether the port is open already*/
	if(PATX == -1)
    {
        /* "PATX" - file descriptor pointing to the opened serial port */
        PATX = open(portName, O_RDWR | O_NOCTTY | O_NDELAY); 
        PortOpenTx = 1;
        /* O_RDWR Read/Write access to serial port           */
        /* O_NOCTTY - No terminal will control the process   */
        /* O_NDELAY -Non Blocking Mode,Does not care about-  */
        /* -the status of DCD line, Open() returns immediatly */				
    }    
    
    /* Error Checking */
    if(PATX == -1)
    {
    	mexPrintf("\nError In Opening Tx Phased Array.\n");
    	return -1;
    }
    else
    {
        mexPrintf("\nTx Phased Array Is Successfully Opened!\n");
    }
    	
	/*---------- Setting the Attributes of the serial port using termios structure --------- */
	struct termios SerialPortSettings;	/* Create the structure                          */
	tcgetattr(PATX, &SerialPortSettings);	/* Get the current attributes of the Serial port */
	cfsetispeed(&SerialPortSettings,B115200); /* Set Read  Speed as 115200                       */
	cfsetospeed(&SerialPortSettings,B115200); /* Set Write Speed as 115200                       */
	SerialPortSettings.c_cflag &= ~PARENB;   /* Disables the Parity Enable bit(PARENB),So No Parity   */
	SerialPortSettings.c_cflag &= ~CSTOPB;   /* CSTOPB = 2 Stop bits,here it is cleared so 1 Stop bit */
	SerialPortSettings.c_cflag &= ~CSIZE;	 /* Clears the mask for setting the data size             */
	SerialPortSettings.c_cflag |=  CS8;      /* Set the data bits = 8                                 */
	SerialPortSettings.c_cflag &= ~CRTSCTS;       /* No Hardware flow Control                         */
	SerialPortSettings.c_cflag |= CREAD | CLOCAL; /* Enable receiver,Ignore Modem Control lines       */ 
	SerialPortSettings.c_iflag &= ~(IXON | IXOFF | IXANY);          /* Disable XON/XOFF flow control both i/p and o/p */
	SerialPortSettings.c_iflag &= ~(ICANON | ECHO | ECHOE | ISIG);  /* Non Cannonical mode                            */
	SerialPortSettings.c_oflag &= ~OPOST;/*No Output Processing*/
	//SerialPortSettings.c_cflag &= ~(PARENB | PARODD);    /* Disables the Parity Enable bit(PARENB),So No Parity   */
    //SerialPortSettings.c_iflag &= ~IGNBRK;          // disable break processing
    SerialPortSettings.c_lflag = 0;                 
    SerialPortSettings.c_oflag = 0;                // no remapping, no delays
	//SerialPortSettings.c_cc[VMIN] = 10; /* Read at least 10 characters */
	//SerialPortSettings.c_cc[VTIME] = 0; /* Wait indefinetly   */
    
	if((tcsetattr(PATX,TCSANOW,&SerialPortSettings)) != 0) /* Set the attributes to the termios structure*/
	{
        mexPrintf("\nError In Setting Attributes For Tx Phased Array.\n");
		return -1;
	}
    return 0;
}	  


/*---------- Close the serial port using termios structure --------- */
int CloseSerialPortTx(void)
{
    if(PATX == -1)
    {
        mexPrintf("\nErrpr: No Opening Serial Port!\n");
    }
    else
    {
        mexPrintf("\nTx Phased Array Is Successfully Closed!\n");
        close(PATX);
        PATX = -1;
    }
    return 0;
}


/*------------------------------- Write data to serial port -----------------------------*/
int WriteCommandTx(const char *input_command, const char *portName)
{   
	int len_Std_Command = strlen(input_command)+2;
    char *Std_Command = (char*)malloc(len_Std_Command*sizeof(unsigned char));
    Std_Command[0] = '\n';
    Std_Command[len_Std_Command-1] = '\r';
        
    /*Copy command from input_command*/
    int i;
    for (i=1; i<=len_Std_Command-2; i++)
    {
        Std_Command[i] = input_command[i-1];
    }
    
    /*Check whether the serial port is opened*/
    int reopen;
	int bytes_written  = 0;  
    if (PortOpenTx != 1)
    {
        reopen = OpenSerialPortTx(portName);
    }
    
    /*Clean the read buffer before writting*/    
    tcflush(PATX, TCIFLUSH);                  
   
    /*Write Command*/
    bytes_written = write(PATX, Std_Command, len_Std_Command);
    //mexPrintf("Bytes written: %d", bytes_written);
    /*char input_command_4letters[4];
    for (i=0; i<=3; i++)
    {
        input_command_4letters[i] = input_command[i];
    }
    if (strstr(input_command_4letters, "txon") != NULL)
    {
        //mexPrintf("\nSleeping for 1 milliseconds for turining on and off antennas: %s\n", input_command);
        usleep(operation_hold_time);
    }
    else
    {
        //mexPrintf("\nSleeping for 1 milliseconds for changing the phase shifter.\n");
        usleep(operation_hold_time);
    }
    usleep(operation_hold_time*2);*/
	return 0;	
}


/*------------------------------- Read data from serial port -----------------------------*/
int ReadCommandTx(void)
{   
	char read_buffer[300];
	int bytes_read = 0;              
	bytes_read = read(PATX, &read_buffer, 300);
    return bytes_read;
}


/*------------------------------- Read data from serial port -----------------------------*/
int ReadCommandTxPrint(void)
{   
	char read_buffer[300];
	int bytes_read = 0;              
	bytes_read = read(PATX, &read_buffer, 300);
    /*Print the number of read bytes and contents*/	
    /*mexPrintf("\nBytes Read: %d\n",bytes_read);*/
	int i = 0;
	for( i = 0; i < bytes_read; i++)	 		
    {
    	mexPrintf("%c",read_buffer[i]);
    }
    mexPrintf("\n");
    return bytes_read;
}

int Operation_Interval_Time(int nMicroSecond)
{
    usleep(nMicroSecond);
    return 0;
}
